import sys
import os
from dataclasses import dataclass, field

from .util import check_order, check_4site_water

import mdtraj as md

from ..system import *

@dataclass
class Mapbuilder:
   itp: list = field(default_factory=lambda: ['topol.itp'])
   gro:  str = "confout.gro" 
   xyz:  str = field(default_factory=lambda: ["solvent.xyz","solute.xyz"])
   xtc:  str = "traj.xtc"
   top:  str = "topology.top"
   method: str = "B3LYP"
   basis:  str = "6-311++G(d,p)"
   vib:  str = "normal"
   nframes: int = 1
   software: str = "Gaussian"
   dir: str = os.getcwd()

   def createJobs(self):
     self.solute, self.solvent = self.extract_solute_solvent() 

   def extract_solute_solvent(self):
     # create a system object
     s = System(self.itp,self.top,self.gro,self.xyz)
     self.atoms, self.molecules, self.atoms_in_mol, self.chem_labels, self.molecule_list = s.read(read_xyz=True)

     # analyze information about the system
     # assuming we have one organic molecule in many solvent molecules
     print(f" >>>>> Looking for a molecule for which the map will be built.")
     res_num = [ x[1] for x in self.molecules ]
     if res_num.count("1") == 1:
        res_index = res_num.index("1")
        s_resid = self.molecules[res_index][0]
        print(f"      The following molecule has been found: {s_resid}")
     else:
        sys.exit(f" In the current implementation MapBuilder can only do 1 molecule in a solvent but you have: {self.molecules} ")

     # find which molecules from xyz files matches 
     # note that this can also be 4-site water so, later, we would
     # need to incorporate this functionality.
     chem_labels_selected_mol=[]
     res_atoms = [x[2] for x in self.molecule_list[res_index]]
     for n in range(len(self.chem_labels)):
        if len(self.chem_labels[n]) == self.atoms_in_mol[res_index]:
           if check_order(res_atoms, self.chem_labels[n]):
              chem_labels_selected_mol.append(s_resid)
              chem_labels_selected_mol.append(self.chem_labels[n])
              chem_labels_selected_mol.append([])
              print(f"      This molecule is matched with the following atoms from {self.xyz[n]} file: {chem_labels_selected_mol} ") 
           else:
              sys.exit(f" Check order of atoms in xyz file for '{self.molecules[res_index][0]}' should be compatible with configuration file: {res_atoms}")
     if not chem_labels_selected_mol:
        sys.exit(f" Cannot find the corresponding xyz file. ")

     
     print(f" >>>>> Looking for a solvent/environment.") 
     # check if we have more than just solute
     assert len(self.molecules) > 1, f" Did not find any molecules besides {self.molecules} "

     # loop over all molecules looking for a solvent
     solv_list = [val[0] for x, val in enumerate(self.molecules) if val[0] not in s_resid ]
     print(f"      The following molecules were found: {solv_list}") 

     # find q chem representation of all solvent molecules
     assert len(solv_list) == 1, f"  So far only one molecule type can be solvent (more complex envirnments will be implemened in the future"

     # find out what solvent it is
     res_names = [ x[0] for x in self.molecules ]
     solvent_atoms=[]
     for solvent in solv_list:
        solv_index = res_names.index(solvent)
        solv_atoms = [ x[2] for x in self.molecule_list[solv_index] ] 
        Msite_ind = check_4site_water(solv_atoms, s.msite_list)
        if Msite_ind > 0:
           print (f"      {solvent} appears to be 4-site water: {solv_atoms}.")
           for n in range(len(self.chem_labels)):
              if check_order(solv_atoms[:Msite_ind]+solv_atoms[Msite_ind+1:], self.chem_labels[n]):
                 solvent_atoms.append(solvent)
                 solvent_atoms.append(self.chem_labels[n])
                 solvent_atoms.append(Msite_ind)
                 print(f"      {solvent} is matched with the following atoms from {self.xyz[n]} file: {self.chem_labels[n]}")
        else:
           # just add like a normal solvent
           for n in range(len(self.chem_labels)):
              if check_order(solv_atoms, self.chem_labels[n]):
                 solvent_atoms.append(solvent)
                 solvent_atoms.append(self.chem_labels[n])
                 solvent_atoms.append([])  

     return chem_labels_selected_mol, solvent_atoms
     

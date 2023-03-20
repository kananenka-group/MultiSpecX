import sys
from dataclasses import dataclass, field

from .util import check_order

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
   soft: str = "Gaussian"

   def createJobs(self):
     # create a system object
     s = System(self.itp,self.top,self.gro,self.xyz)
     self.atoms, self.molecules, self.atoms_in_mol, self.chem_labels, self.molecule_list = s.read(read_xyz=True)

     # analyze information about the system
     # assuming we have one organic molecule in many solvent molecules
     print(f">>>>> Looking for a molecule for which the map will be built.")
     res_num = [ x[1] for x in self.molecules ]
     if res_num.count("1") == 1:
        res_index = res_num.index("1")
        s_resid = self.molecules[res_index][0]
        print(f"      The following molecule has been found: {s_resid}")
     else:
        sys.exit(f" In the current implementation MapBuilder can only do 1 molecule in a solvent but you have: {self.molecules} ")

     # find which molecules from xyz files matches 
     chem_labels_selected_mol=[]
     res_atoms = [x[2] for x in self.molecule_list[res_index]]
     for n in range(len(self.chem_labels)):
        if len(self.chem_labels[n]) == self.atoms_in_mol[res_index]:
           if check_order(res_atoms, self.chem_labels[n]):
              chem_labels_selected_mol = self.chem_labels[n]
           else:
              sys.exit(f" Check order of atoms in xyz file for '{self.molecules[res_index][0]}' should be compatible with configuration file: {res_atoms}")
     if not chem_labels_selected_mol:
        sys.exit(f" Cannot find the corresponding xyz file. ")

     # print:
     print(f"      This molecule is matched with the following atoms from xyz files: {chem_labels_selected_mol} ") 
     print(f">>>>> Looking for a solvent/environment.") 
     
     # check if we have more than just solute
     assert len(self.molecules) > 1, f" Did not find any molecules besides {self.molecules} "

     # loop over all molecules looking for a solvent
     solv_list = [val[0] for x, val in enumerate(self.molecules) if val[0] not in s_resid ]
     print(f"      The following molecules were found: {solv_list}") 

 

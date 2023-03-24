import sys
import os
from pathlib import Path
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
   ncores: int = 1
   wdir: str = os.getcwd()
   # cut-offs for explicit solvent and point charges
   # default values are taken for water from Skinner papers
   cut_off1: float = 4.0
   cut_off2: float = 7.831
   # transform system:
   transform: list = field(default_factory=[])
   

   def createJobs(self):
     """
        Read snapshot from MD simulation, turn MD configuration into
        input file for Qchem calculation
     """
     self.solute, self.solvent, self.solu_ind, self.solv_ind = self.extract_solute_solvent() 

     t = md.load(self.xtc, top=self.gro)
     for frame in range(self.nframes):
        xyz = 10.0*t.xyz[frame,:,:]
        solu_xyz_raw = xyz[self.solu_ind,:]
        solv_xyz_raw = xyz[self.solv_ind,:]

        # transform here
        solu_xyz, solv_xyz = self.transformXYZ(solu_xyz_raw, solv_xyz_raw)

        # calculate COM for solvent...

        # input file
        self.QC_input_file(solu_xyz, solv_xyz, frame)

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
                 ignore=[]
                 ignore.append(Msite_ind)
                 solvent_atoms.append(ignore)
                 print(f"      {solvent} is matched with the following atoms from {self.xyz[n]} file: {self.chem_labels[n]}")
        else:
           # just add like a normal solvent
           for n in range(len(self.chem_labels)):
              if check_order(solv_atoms, self.chem_labels[n]):
                 solvent_atoms.append(solvent)
                 solvent_atoms.append(self.chem_labels[n])
                 solvent_atoms.append([])  

     # find out indices of solute and solvent in the configuration 
     solu_idx = [ int(atom[0])-1 for atom in self.atoms if atom[7] == s_resid]

     # this line below works for a single solvent, need to extend it in the future
     # to support more than one type of molecule in the environment
     solv_idx = [ int(atom[0])-1 for atom in self.atoms if atom[7] in solv_list]
     
     return chem_labels_selected_mol, solvent_atoms, solu_idx, solv_idx

     
   def QC_input_file(self, solute, solvent, frame):
      """
         Here we will write an input file for electronic structure calculation
      """
      path = Path(os.path.join(self.wdir,f"{frame}"))
      path.mkdir(parents=True, exist_ok=True)

      if self.software.lower() == "gaussian":
         self.write_Gaussian_input(solute, solvent, path)
      else:
         sys.exit(" Unrecognized option: software. At this time only Gaussian is supported")



   def write_Gaussian_input(self, su_xyz, sv_xyz, path):
      """
         Gaussian DFT job:

         Normal modes:
         -------------
         Freeze all solvent atoms, perform geometry optimization and frequency calculation
        
         iop(7/33=1) calculates transition dipole moments w.r.t. normal modes
                     output units are (km/mol)^1/2 read here
                     https://mattermodeling.stackexchange.com/questions/5021/what-units-are-used-in-gaussian-16-for-dipole-derivatives-output
                     to get it in [D A^(-1) u^(-1/2)] units divide by sqrt(42.2561)

      """      
      input_file = path/"input.com"
      solute_atoms_to_ignore = self.solute[2]
      solute_atoms_list = iter(self.solute[1])

      solvent_atoms_ignore = self.solvent[2]
      n_solvent_atoms = len(self.solvent[1])+len(self.solvent[2])
      n_solvent_mols = sv_xyz.shape[0] // n_solvent_atoms

      with open(input_file,"w") as f:
         f.write(f"%nprocshared={self.ncores}\n")
         f.write(f"# Opt Freq {self.method}/{self.basis} iop(7/33=1) NoSymm Int=Ultrafine SCF=tight Test\n")
         f.write(" \n")
         f.write(f"{self.solute[0]} in solvent \n")
         f.write(" \n")
         f.write("0 1\n")
         for n in range(su_xyz.shape[0]):
            if n not in solute_atoms_to_ignore:
               f.write(f"  {next(solute_atoms_list)}   {su_xyz[n,0]:.4f}   {su_xyz[n,1]:.4f}   {su_xyz[n,2]:.4f}\n")

         for n in range(n_solvent_mols):
            solvent_atoms_list = iter(self.solvent[1])
            for m in range(n_solvent_atoms):
               atom_index=n_solvent_atoms*n+m
               if m not in solvent_atoms_ignore:
                  f.write(f"  {next(solvent_atoms_list)}   {sv_xyz[atom_index,0]:.4f}   {sv_xyz[atom_index,1]:.4f}   {sv_xyz[atom_index,2]:.4f}\n")
 
         f.write(" \n")

   def transformXYZ(self, solu_xyz, solv_xyz):
      """
          Here input coordinates for the solute and solvent will be
          transformed
      """
      solu_xyz_t = np.copy(solu_xyz)
      solv_xyz_t = np.copy(solv_xyz)

      for item in self.transform:
         # center frame on a given atom
         if item[0] == "center":
            atom_center=item[1]-1
            xyz_shift = solu_xyz[atom_center,:]
            solu_xyz_t -= xyz_shift
            solv_xyz_t -= xyz_shift
      return solu_xyz_t, solv_xyz_t

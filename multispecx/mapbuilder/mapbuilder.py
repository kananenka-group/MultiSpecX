import sys
import os
from pathlib import Path
from dataclasses import dataclass, field

from .util import check_order, check_4site_water, rotation_matrix, minImage, PBC, distance, inBox

import mdtraj as md

from ..system import *

@dataclass
class Mapbuilder:
   solv_ref_atom: str
   solu_ref_atom: str

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
   opt_cycles: int = 200
   wdir: str = os.getcwd()
   # cut-offs for explicit solvent and point charges
   # default values are taken for water from Skinner papers
   cut_off1: float = 4.0
   cut_off2: float = 7.831
   transform: list = field(default_factory=[])
   

   def createJobs(self):
     """
        Read snapshot from MD simulation, turn MD configuration into
        input file for Qchem calculation
     """
     self.solute, self.solvent, self.solu_ind, self.solv_ind, self.solv_charge, self.solu_charge = self.extract_solute_solvent() 
     self.print_transform_info()

     t = md.load(self.xtc, top=self.gro)
     ws_count=0
     for frame in range(self.nframes):
        xyz = 10.0*t.xyz[frame,:,:]
        box = 10.0*t.unitcell_lengths[frame,:]

        solu_xyz_raw = xyz[self.solu_ind,:]
        solv_xyz_raw = xyz[self.solv_ind,:]

        # center box on selected atom
        solu_xyz, solv_xyz = self.centerBox(solu_xyz_raw, solv_xyz_raw, box)

        # explicit or emplicit solvent
        solv_e_xyz, solv_p_xyz = self.SplitSolvent(solu_xyz, solv_xyz, box)
        ws_count+=solv_e_xyz.shape[0]

        # transform here
        solu_xyzT, solv_e_xyzT, solv_p_xyzT = self.transformXYZ(solu_xyz, solv_e_xyz, solv_p_xyz, box)

        # input file
        self.QC_input_file(solu_xyzT, solv_e_xyzT, solv_p_xyzT, frame)

     print(f"      Average number of solvent molecules within {self.cut_off1} A cutoff is {ws_count/self.nframes}")

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
        print(f"       The following molecule has been found: {s_resid}")
     else:
        sys.exit(f" In the current implementation MapBuilder can only do 1 molecule in a solvent but you have: {self.molecules} ")

     # find which molecules from xyz files matches 
     chem_labels_selected_mol=[]
     res_atoms  = [x[3] for x in self.molecule_list[res_index]]
     atom_types = [x[0] for x in self.molecule_list[res_index]]
     for n in range(len(self.chem_labels)):
        if len(self.chem_labels[n]) == self.atoms_in_mol[res_index]:
           if check_order(res_atoms, self.chem_labels[n]):
              chem_labels_selected_mol.append(s_resid)
              chem_labels_selected_mol.append(res_atoms)
              chem_labels_selected_mol.append(atom_types)
              chem_labels_selected_mol.append(self.chem_labels[n])
              #chem_labels_selected_mol.append([])
              print(f"       This molecule is matched with the following atoms from {self.xyz[n]} file: {chem_labels_selected_mol[2]} ") 
           else:
              sys.exit(f" Check order of atoms in xyz file for '{self.molecules[res_index][0]}' should be compatible with configuration file: {res_atoms}")
     if not chem_labels_selected_mol:
        sys.exit(f" Cannot find the corresponding xyz file. ")
     
     print(f" >>>>> Looking for a solvent/environment.") 
     # check if we have more than just solute
     assert len(self.molecules) > 1, f" Did not find any molecules besides {self.molecules} "

     # loop over all molecules looking for a solvent
     solv_list = [val[0] for x, val in enumerate(self.molecules) if val[0] not in s_resid ]
     print(f"       The following molecules were found: {solv_list}") 

     # find q chem representation of all solvent molecules
     assert len(solv_list) == 1, f"  So far only one molecule type can be solvent (more complex envirnments will be implemened in the future"

     # find out what solvent it is
     res_names = [ x[0] for x in self.molecules ]
     solvent_atoms=[]
     for solvent in solv_list:
        solv_index = res_names.index(solvent)
        solv_atoms = [ x[3] for x in self.molecule_list[solv_index] ] 
        solv_atom_types = [ x[0] for x in self.molecule_list[solv_index] ] 
        Msite_ind = check_4site_water(solv_atoms, s.msite_list)
        if Msite_ind > 0:
           print (f"      {solvent} appears to be 4-site water: {solv_atoms}.")
           sys.exit(" Error! Since so far only ONIOM with amber FF has been implemented it is not clear how to map dummy atoms to real ones.\n Use TIP3P or SPC water.")
           for n in range(len(self.chem_labels)):
              if check_order(solv_atoms[:Msite_ind]+solv_atoms[Msite_ind+1:], self.chem_labels[n]):
                 solvent_atoms.append(solvent)
                 solvent_atoms.append(solv_atoms)
                 solvent_atoms.append(self.chem_labels[n])
                 ignore=[]
                 ignore.append(Msite_ind)
                 solvent_atoms.append(ignore)
                 # add solvent charges here...
                 print(f"      {solvent} is matched with the following atoms from {self.xyz[n]} file: {self.chem_labels[n]}")
        else:
           # just add like a normal solvent
           for n in range(len(self.chem_labels)):
              if check_order(solv_atoms, self.chem_labels[n]):
                 solvent_atoms.append(solvent)
                 solvent_atoms.append(solv_atoms)
                 solvent_atoms.append(solv_atom_types)
                 solvent_atoms.append(self.chem_labels[n])
                 #solvent_atoms.append([])  

     # find out indices of solute and solvent in the configuration 
     solu_idx = [ int(atom[0])-1 for atom in self.atoms if atom[8] == s_resid]

     # this line below works for a single solvent, need to extend it in the future
     # to support more than one type of molecule in the environment
     solv_idx = [ int(atom[0])-1 for atom in self.atoms if atom[8] in solv_list]

     # get solvent charges
     solv_charge = [ float(atom[6]) for atom in self.atoms if atom[8] in solv_list ]

     # get solute charges
     solu_charge = [ float(atom[6]) for atom in self.atoms if atom[8] == s_resid ]

     # number of atoms in solvent molecule:
     return chem_labels_selected_mol, solvent_atoms, solu_idx, solv_idx, solv_charge, solu_charge

     
   def QC_input_file(self, solute, solvent1, solvent2, frame):
      """
         Here we will write an input file for electronic structure calculation
      """
      path = Path(os.path.join(self.wdir,f"{frame}"))
      path.mkdir(parents=True, exist_ok=True)

      if self.software.lower() == "gaussian":
         self.write_Gaussian_input(solute, solvent1, solvent2, path)
      else:
         sys.exit(" Unrecognized option: software. At this time only Gaussian is supported")



   def write_Gaussian_input(self, su_xyz, sv1_xyz, sv2_xyz, path):
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
      freeze = -1
      nofreeze=0

      input_file = path/"input.com"

      # calculation type 1 normal modes
      #
      if self.vib.lower() == "normal":
         with open(input_file,"w") as f:
            # Step 1. Geometry optimization without point charges
            f.write(f"%nprocshared={self.ncores}\n")
            f.write("%chk=job.chk\n")
            f.write("%mem=20Gb\n")
            #f.write(f"#p Opt(MaxCycles={self.opt_cycles},CalcFC) oniom({self.method}/{self.basis}:amber) NoSymm \n")
            f.write(f"#p Opt(MaxCycles={self.opt_cycles},CalcFC) {self.method}/{self.basis} \n")
            f.write(" \n")
            f.write(f"{self.solute[0]} in solvent \n")
            f.write(" \n")
            #f.write("0 1 0 1 0 1\n")
            f.write("0 1\n")
            # high-level
            #for n in range(su_xyz.shape[0]):
            #   f.write(f"  {self.solute[3][n].upper()}-{self.solute[2][n].upper()}-{self.solu_charge[n]}   {nofreeze}   {su_xyz[n,0]:.4f}   {su_xyz[n,1]:.4f}   {su_xyz[n,2]:.4f}  H \n")
            for n in range(su_xyz.shape[0]):
               f.write(f"  {self.solute[3][n].upper()}   {nofreeze}   {su_xyz[n,0]:.4f}   {su_xyz[n,1]:.4f}   {su_xyz[n,2]:.4f}  \n")

            for n in range(sv1_xyz.shape[0]):
               for m in range(sv1_xyz.shape[1]):
                  #f.write(f"  {self.solvent[3][m].upper()}-{self.solvent[2][m].upper()}-{self.solv_charge[m]}   {freeze}   {sv1_xyz[n,m,0]:.4f}   {sv1_xyz[n,m,1]:.4f}   {sv1_xyz[n,m,2]:.4f} H \n")
                  f.write(f"  {self.solvent[3][m].upper()}   {freeze}   {sv1_xyz[n,m,0]:.4f}   {sv1_xyz[n,m,1]:.4f}   {sv1_xyz[n,m,2]:.4f} \n")
            
            #f.write(" \n")
            # low-level
            #for n in range(sv2_xyz.shape[0]):
            #   for m in range(sv2_xyz.shape[1]):
            #      f.write(f"  {self.solvent[3][m].upper()}-{self.solvent[2][m].upper()}-{self.solv_charge[m]}  {freeze} {sv2_xyz[n,m,0]:.4f}   {sv2_xyz[n,m,1]:.4f}   {sv2_xyz[n,m,2]:.4f} L \n")

            f.write(" \n")
            f.write("--Link1--\n")
            f.write(f"%nprocshared={self.ncores}\n")
            f.write("%chk=job.chk\n")
            f.write(f"#p Freq {self.method} ChkBasis iop(7/33=1) Geom=AllCheck Guess=Read NoSymm \n")
            f.write(" \n") 
      elif self.vib.lower() == "local":
         sys.exit(" Local mode calculations have not been implemented.")
      else:
         sys.exit(f" Unrecognized vib option: {self.vib}.")    

   def transformXYZ(self, solu_xyz, solv_e_xyz, solv_p_xyz, box):
      """
          Here input coordinates for the solute and solvent will be
          transformed
      """
      thresh=1.0e-4
      solu_xyz_t   = np.copy(solu_xyz)
      solv_e_xyz_t = np.copy(solv_e_xyz)
      solv_p_xyz_t = np.copy(solv_p_xyz)

      for item in self.transform:
         # center frame on a given atom
         if item[0].lower() == "center":
            atom_center=item[1]-1
            xyz_shift = np.copy(solu_xyz_t[atom_center,:])
            if np.linalg.norm(xyz_shift) > 1.0e-4:
               solu_xyz_t[:,:]     = np.subtract(solu_xyz_t[:,:],xyz_shift)
               solv_e_xyz_t[:,:,:] = np.subtract(solv_e_xyz_t[:,:,:],xyz_shift)
               solv_p_xyz_t[:,:,:] = np.subtract(solv_p_xyz_t[:,:,:],xyz_shift)
         # align w.r.t particular axis
         elif item[0].lower() == "rotate":
            atomn = item[1]-1
            atomp = item[2]-1
            va = np.subtract(solu_xyz_t[atomp,:],solu_xyz_t[atomn,:])
            if item[3].lower() == "z":
               vb = np.array([0.0,0.0,1.0])
            elif item[3].lower() == "y":
               vb = np.array([0.0,1.0,0.0])
            elif item[3].lower() == "x":
               vb = np.array([1.0,0.0,0.0])
            else:
               sys.exit(f" Cannot recognize the rotation axis {item[3]} can only be 'x', 'y', or 'z'") 
            Rot = rotation_matrix(va,vb) 

            # rotate all atoms here ...
            if np.linalg.norm(Rot-np.eye(3)) > 1.0e-4:
               for n in range(solu_xyz_t.shape[0]):
                  solu_xyz_t[n,:] = Rot.dot(solu_xyz_t[n,:])
               for n in range(solv_e_xyz_t.shape[0]):
                  for m in range(solv_e_xyz_t.shape[1]):
                     solv_e_xyz_t[n,m,:] = Rot.dot(solv_e_xyz_t[n,m,:])
               for n in range(solv_p_xyz_t.shape[0]):
                  for m in range(solv_p_xyz_t.shape[1]):
                     solv_p_xyz_t[n,m,:] = Rot.dot(solv_p_xyz_t[n,m,:])
         else:
            sys.exit(" Cannot recognize this transformation operation {item}")

      # check the distance w.r.t ref atom before and after transformaton:
      su_er = np.linalg.norm(solu_xyz[atom_center,:]-solu_xyz[:,:]) - np.linalg.norm(solu_xyz_t[atom_center,:]-solu_xyz_t[:,:])
      if np.max(np.abs(su_er)) > thresh:
         sys.exit(f" Error with solvent coordinate transformation {su_er}")

      sv_er1 = np.linalg.norm(solu_xyz[atom_center,:]-solv_e_xyz[:,:,:]) - np.linalg.norm(solu_xyz_t[atom_center,:]-solv_e_xyz_t[:,:,:])
      if np.max(np.abs(sv_er1)) > thresh:
         sys.exit(f" Error with solvent coordinate transformation {sv_er1}")

      sv_er2 = np.linalg.norm(solu_xyz[atom_center,:]-solv_p_xyz[:,:,:]) - np.linalg.norm(solu_xyz_t[atom_center,:]-solv_p_xyz_t[:,:,:])
      if np.max(np.abs(sv_er2)) > thresh:
         sys.exit(f" Error with solvent coordinate transformation {sv_er2}")

      return solu_xyz_t, solv_e_xyz_t, solv_p_xyz_t

   def print_transform_info(self):
      """
         Just print transformation operations to be applied
      """ 
      if self.transform:
         print(f">>>>> Coordinate transformation:")

      for item in self.transform:
         if item[0].lower() == "center":
            print(f"      The frames will be translated to put atom {self.solute[3][item[1]-1]}({item[1]}) at the center of the box.")
         if item[0].lower() == "rotate":
            print(f"      The frames will be rotated such that vector {self.solute[3][item[1]-1]}({item[1]}) -> {self.solute[3][item[2]-1]}({item[2]}) will be aligned with the positive {item[3]} axis")

 
   def SplitSolvent(self, solu_xyz, solv_xyz, box):
      """
         Split solvent molecules into 2 groups:
         - explicit solvent 
         - solvent to be used as point charges
      """
      # find reference atoms
      sura = self.solute[1].index(self.solu_ref_atom)
      svra = self.solvent[1].index(self.solv_ref_atom)
      n_solv_atoms = len(self.solvent[1])
      n_solv_mols  = solv_xyz.shape[0] // n_solv_atoms

      sv_qm = []  
      sv_pc = []
      for n in range(n_solv_mols):
         vnr = solv_xyz[n_solv_atoms*n+svra] - solu_xyz[sura,:]
         vnn = minImage(vnr,box)      
         dnn = np.sqrt(vnn.dot(vnn))
         if dnn < self.cut_off2:
            if dnn < self.cut_off1:
               sv_qm.append(solv_xyz[n_solv_atoms*n:n_solv_atoms*n+n_solv_atoms,:])
            else:
               sv_pc.append(solv_xyz[n_solv_atoms*n:n_solv_atoms*n+n_solv_atoms,:])

      return np.asarray(sv_qm,dtype=np.float32), np.asarray(sv_pc,dtype=np.float32)
             
   def centerBox(self, solu_xyz, solv_xyz, box):
      """
         Center box on the reference atom
      """
      su_xyz = np.copy(solu_xyz)
      sv_xyz = np.copy(solv_xyz)

      sura = self.solute[1].index(self.solu_ref_atom)
      cbox = box/2

      shift = solu_xyz[sura,:] - cbox
      for n in range(solu_xyz.shape[0]):
         su_xyz[n,:] = inBox(su_xyz[n,:]-shift,box)
      for n in range(solv_xyz.shape[0]):
         sv_xyz[n,:] = inBox(sv_xyz[n,:]-shift,box)

      return su_xyz, sv_xyz         


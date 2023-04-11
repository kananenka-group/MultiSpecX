import mdtraj as md
import time
import os
from dataclasses import dataclass, field
from pathlib import Path

from ..system import *
from ..amideI import *

@dataclass
class Ester:
   nframes: int = None
   type: list = field(default_factory=lambda: ['dipole'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   ester_unit: list = field(default_factory=lambda: ['C','O','OE','CA'])
   gro:  str = "confout.gro"
   xtc:  str = "traj.xtc" 
   top:  str = "topology.top"
   outDir: str = os.getcwd()
   isotope_labels: list = field(default_factory=lambda: [])
   transform: list = field(default_factory=lambda: [])
   start: int = 1
   elFmap: str  = "Edington2016"
   freq_shift: float = 0.0

   def checkTrajectory(self) -> int:
      """
         Perform some checks here.

         Input:
         ------
         system information (gro, top, xtc files)

         Output:
         -------
         number of frames in MD trajectory
      """
      start_time = time.time()
      printDT("starts")
      print(f" >>>>> Checking system and MD xtc trajectory.")

      s = System(self.itp,self.top,self.gro)
      self.atoms, self.molecules, self.atoms_in_mol, _, _ = s.read()

      t = md.load(self.xtc, top=self.gro)
      nframes = t.xyz.shape[0]
      print(f" >>>>> Reading frames from {self.xtc} file")
      print(f"       Total number of frames to read: {nframes}")

      end_time = time.time()
      print(f" >>>>> Execution time: {(end_time - start_time)/60:.1f} minutes")
      printDT("ends")
      return nframes


   def generateHamiltonian(self) -> None: 
      start_time = time.time()
      printDT("starts")

      emap_cut = 21.0  # units=Angstrom

      # create a system object
      s = System(self.itp,self.top,self.gro)
      self.atoms, self.molecules, self.atoms_in_mol, _, _ = s.read()

      # find indices of ester groups in each molecule
      chrom_idx, ester_list_idx, n_ester_mol = chromList(self.isotope_labels, self.ester_unit, self.atoms, self.atoms_in_mol)
      if not chrom_idx:
        print(f" Did not find any ester groups like this: {self.ester_unit}")
        sys.exit(" exiting...")
      else:
        print(f"       Found {len(chrom_idx)} ester chromophores: ")
        [print (f"       {id} in {molid[0]} ") for (id,molid) in zip(n_ester_mol,self.molecules) if id>0]

      # determine where charge groups start
      cgS = chargeGroupSt(self.atoms)

      # build coordiate transformation here
      self.transform_internal = getInternalTransformXYZ(self.transform, self.ester_unit, chrom_idx)

      # figure out how to use map:
      map_w0, Elst_map, ef_atoms, ef_proj = self.getMap()

      # prepare some variables 
      charges = np.array([ x[6] for x in self.atoms ],dtype=np.float32)
      masses  = np.array([ x[7] for x in self.atoms ],dtype=np.float32)

      # create directory
      path = Path(self.outDir)
      path.mkdir(parents=True, exist_ok=True)
      print(f" >>>>> Output files will be written into {path.absolute()} directory.")

      # loop over all frames
      t = md.load(self.xtc, top=self.gro)
      nframes = t.n_frames

      if not self.nframes:
        self.nframes = nframes

      if t.n_atoms != len(self.atoms):
         raise ValueError(f" Number of atoms in xtc file {t.n_atoms} does not match number of atoms in config. files {len(self.atoms)}")

      Energy = np.zeros((self.nframes,len(chrom_idx),len(chrom_idx)),dtype=np.float32)  
      Dipole = np.zeros((self.nframes,len(chrom_idx)*3),dtype=np.float32)

      print(f" >>>>> Reading frames from {self.xtc} file") 
      print(f"       Total number of frames to read: {self.nframes}")
      print(f"       Timestep: {t.timestep:.3f} [ps]")

      for frame in range(self.nframes):

         xyz_raw = NMTOA*t.xyz[frame,:,:]                #units=A
         box     = NMTOA*t.unitcell_lengths[frame,:]     #units=A

         tdv_f  = np.zeros((len(chrom_idx),3),dtype=np.float32)
         tdp_f  = np.zeros((len(chrom_idx),3),dtype=np.float32)

         # loop over all chromphores
         for chind, chrom in enumerate(chrom_idx):

            xyz_chrom_raw = xyz_raw[chrom,:]

            # calculate transition dipole moments before we do anything with the box
            tdv_f[chind,:], tdp_f[chind,:] = self.ester_TDC_Wang20(xyz_chrom_raw, box)
            Dipole[frame,3*chind:3*chind+3] = np.copy(tdv_f[chind,:])

            # re-center box at the COM of selected atoms
            com_raw = getCOM(xyz_chrom_raw, masses[chrom])
            xyz = centerBox(xyz_raw, com_raw, box)
            xyz_chrom = xyz[chrom,:]
            com = getCOM(xyz_chrom, masses[chrom])

            # determine COM for all charge groups
            comCg = getCOMChg(xyz, cgS, masses)

            # save atoms that contribute to efield of this chromophore 
            cg_atoms_include = include_CG_atoms(comCg, com, emap_cut, cgS)
            # remove atoms that belong to the chromophore
            atoms_include = exclude_chrom_atoms(cg_atoms_include, chrom)
            # remove other atoms here?
           
            # coordinate transformation for all included_atoms and ester unit
            ester_t, envr_t = transformXYZ(self.transform_internal[chind], xyz_chrom, xyz[atoms_include,:])

            # calculate electric field components at selected atoms
            efc, efp = calcEf(ef_atoms, ef_proj, envr_t, ester_t, charges[atoms_include])
           
            # map goes here
            Energy[frame,chind,chind] = map_w0 + self.freq_shift + np.dot(Elst_map, efc)
          
         # calculate TDC:
         for i1 in range(len(chrom_idx)):
            for i2 in range(i1):
               Energy[frame,i1,i2] = Energy[frame,i2,i1] = TDC(tdv_f[i1], tdv_f[i2], tdp_f[i1], tdp_f[i2], box)
    
      # print Hamiltonian and TD into file
      printEnergy(path,Energy) 
      printDipole(path,Dipole)
     
      # finish here
      end_time = time.time()
      print(f" >>>>> Execution time: {(end_time - start_time)/60:.1f} minutes")
      printDT("ends")

   def ester_TDC_Wang20(self, xyz, box):
      """
         Here we will use TD for the ester developed by Lu Wang in 
         J. Chem. Phys. 153, 035101 (2020)
         
      """
      tdAngle: float = 15.1*np.pi/180.0
      tdMag: float   = 2.10   # units are D*A^{-1}*u^{-1/2}, u is a.m.u.

      vcod = minImage(np.subtract(xyz[1,:],xyz[0,:]),box)
      vcoe = minImage(np.subtract(xyz[2,:],xyz[0,:]),box)

      vC = np.linalg.norm(vcod)

      vcod /= np.linalg.norm(vcod)
      vcoe /= np.linalg.norm(vcoe)

      n1  = np.cross(vcod,vcoe)
      n2  = np.cross(n1,vcod)
      n2 /= np.linalg.norm(n2) 

      tdv  = np.sin(tdAngle)*n2 - np.cos(tdAngle)*vcod
      tdv /= np.linalg.norm(tdv)
      tdv *= tdMag

      # td position is midway CO bond, see Lu Wang paper
      tdp = xyz[0,:] + 0.5*vcod*vC   
 
      return tdv, tdp

   def getMap(self):
      """
         Return map 
      """
      map_w0: float
      ef_atoms = []
      ef_projection = []   
      Elst_map = []

      if self.elFmap == "Edington2016":
         print(f' >>>>> Ester frequency map from Edington et al., J. Phys. Chem. A 2016, 10, 3888-3896 will be used.')
         map_w0 = 1745.0
         Elst_map = [1967.6, -640.4, -835.4, 1154.6, -1964.2, 0.0, 0.0, -2776.0, 0.0]
         ef_atoms = [0, 1, 2]
         ef_projection = []
      else:
         sys.exit(f" Only Baiz 2016 map is currently implemented.")

      return map_w0, np.array(Elst_map,dtype=np.float32), ef_atoms, ef_projection

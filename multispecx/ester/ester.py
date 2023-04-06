import mdtraj as md
import time
from ..system import *
from dataclasses import dataclass, field

from ..amideI import *

@dataclass
class Ester:
   nframes: int = None
   type: list = field(default_factory=lambda: ['dipole'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   ester_unit: list = field(default_factory=lambda: ['C','O','C','O'])
   gro:  str = "confout.gro"
   xtc:  str = "traj.xtc" 
   top:  str = "topology.top"
   isotope_labels: list = field(default_factory=lambda: [])
   transform: list = field(default_factory=lambda: [])
   start: int = 1
   elFmap: str  = "Baiz2016"
   
   freq_shift: float = 0.0

   def generateHamiltonian(self): 
     start_time = time.time()
     printDT("starts")

     emap_cut = 2.1*NMTOAU

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
     cgO = chargeGroupSt(self.atoms)
     cgS = np.append(cgO,[len(self.atoms)])

     # build coordiate transformation here
     self.transform_internal = getInternalTransformXYZ(self.transform, self.ester_unit, chrom_idx)

     # figure out how to use map:
     #self.getMap()

     # prepare some variables for fast processing
     charges = np.asarray([ x[6] for x in self.atoms ],dtype=np.float32)
     masses  = np.asarray([ x[7] for x in self.atoms ],dtype=np.float32)

     # loop over all frames
     t = md.load(self.xtc, top=self.gro)
     nframes = t.xyz.shape[0]

     if not self.nframes:
        self.nframes = nframes

     Energy = np.zeros((self.nframes,len(chrom_idx),len(chrom_idx)),dtype=np.float32)  
     Dipole = np.zeros((self.nframes,len(chrom_idx)*3),dtype=np.float32)

     print(f" >>>>> Reading frames from {self.xtc} file") 
     print(f"       Total number of frames to read: {self.nframes}")

     for frame in range(self.nframes):
        xyz_raw = NMTOAU*t.xyz[frame,:,:]
        box     = NMTOAU*t.unitcell_lengths[frame,:]

        tdv_f  = np.zeros((len(chrom_idx),3),dtype=np.float32)
        tdp_f  = np.zeros((len(chrom_idx),3),dtype=np.float32)
        # loop over all chromphores
        for chind, chrom in enumerate(chrom_idx):

           xyz_chrom_raw = xyz_raw[chrom,:]

           # calculate transition dipole moments before we do anything with the box
           tdv_f[chind,:], tdp_f[chind,:], tdMag = self.ester_TDC_Wang20(xyz_chrom_raw, box)
           Dipole[frame,3*chind:3*chind+3] = np.copy(tdv_f[chind,:])

           # re-center box at the COM of selected atoms
           com_raw = getCOM(xyz_chrom_raw, masses[chrom])
           xyz = centerBox(xyz_raw, com_raw, box)
           xyz_chrom = xyz[chrom,:]
           com = getCOM(xyz_chrom, masses[chrom])

           # determine COM for all charge groups
           comCg = getCOMChg(xyz, cgS, masses)

           # save atoms that contribute to efield of this chromophore 
           atoms_include = AinF(comCg, chrom, com, emap_cut, cgS)
           
           # coordinate transformation for all included_atoms and ester unit
           ester_t, envr_t = transformXYZ(self.transform_internal[chind], xyz_chrom, xyz[atoms_include,:])

           # calculate electric field components at selected atoms
           ef_atoms_num = [0, 1, 2] 
           efc, efp = calcEf(ef_atoms_num, envr_t, ester_t, charges[atoms_include])
           
           # map goes here
           map_w0: float = 1745.0
           Elst_map = np.array([1967.6, -640.4, -835.4, 1154.6, -1964.2, 0.0, 0.0, -2776.0, 0.0])
           w = map_w0 + self.freq_shift + np.dot(Elst_map, efc)
          
           Energy[frame,chind,chind] = w

        # calculate TDC:
        for i1 in range(len(chrom_idx)):
           for i2 in range(i1):
              Energy[frame,i1,i2] = Energy[frame,i2,i1] = TDC(tdv_f[i1], tdv_f[i2], tdp_f[i1], tdp_f[i2], tdMag)
    
     # print Hamiltonian into file
     printEnergy(Energy) 
     printDipole(Dipole)
     
     # finish here
     end_time = time.time()
     print(f" >>>>> The execution time: {(end_time - start_time)/60:.1f} minutes")
     printDT("ends")

   def ester_TDC_Wang20(self, xyz, box):
      """
         Here we will use TD for the ester developed by Lu Wang in J. Chem. Phys. 153, 035101 (2020)
         see also Wei Zhuang paper
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
 
      return tdv, tdp, tdMag


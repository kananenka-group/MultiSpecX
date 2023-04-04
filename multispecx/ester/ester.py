import mdtraj as md
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
     emap_cut = 20.0

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
     Dipole = np.zeros((self.nframes,len(chrom_idx),3),dtype=np.float32)

     print(f" >>>>> Reading frames from {self.xtc} file") 
     print(f"       Total number of frames to read: {self.nframes}")

     w_avg: float = 0.0
     for frame in range(self.nframes):
        xyz_raw = 10.0*t.xyz[frame,:,:]
        box     = 10.0*t.unitcell_lengths[frame,:]

        TDC_f  = np.zeros((len(chrom_idx),3),dtype=np.float32)
        # loop over all chromphores
        for chind, chrom in enumerate(chrom_idx):

           # re-center box at the COM of selected atoms
           xyz_chrom_raw = xyz_raw[chrom,:]
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
           Baiz_map = np.array([1967.6, -640.4, -835.4, 1154.6, -1964.2, 0.0, 0.0, -2776.0, 0.0])
           w = 1745.0 + 1967.6*efc[0,0] - 640.4*efc[0,1] - 835.4*efc[0,2] + 1154.6*efc[1,0] - 1964.2*efc[1,1] -2776.0*efc[2,1]
           #w = 1745.0 + self.freq_shift + np.multiply(Baiz_map, efc)
           w_avg += w
          
           Energy[frame,chind,chind] = w

           # calc. transition dipole vector for this chrom.
           TDC_f[chind,:] = self.ester_TDC(ester_t)
     
     print (f" Average frequency {w_avg/(len(chrom_idx)*(frame+1))}")

   def ester_TDC(self,xyz):
      tdcv = np.zeros((3))
      return tdcv

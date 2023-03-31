import mdtraj as md
from ..system import *
from dataclasses import dataclass, field

from ..amideI import *

@dataclass
class Ester:
   type: list = field(default_factory=lambda: ['dipole'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   ester_unit: list = field(default_factory=lambda: ['C','O','C','O'])
   gro:  str = "confout.gro"
   xtc:  str = "traj.xtc" 
   top:  str = "topology.top"
   nframes: int = 1
   isotope_labels: list = field(default_factory=lambda: [])
   transform: list = field(default_factory=lambda: [])
   start: int = 1
   elFmap: str  = "BaizXXX"
   freq_shift: float = 0.0

   def generateHamiltonian(self): 
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

     # prepare some variables for fast processing
     charges = np.asarray([ x[6] for x in self.atoms ],dtype=np.float32)
     masses  = np.asarray([ x[7] for x in self.atoms ],dtype=np.float32)

     # loop over all frames
     t = md.load(self.xtc, top=self.gro)
     nframes = t.xyz.shape[0]

     print(f" >>>>> Reading frames from {self.xtc} file") 
     print(f"       Total number of frames to read: {nframes}")
     for frame in range(nframes):
        xyz_raw = 10.0*t.xyz[frame,:,:]
        box     = 10.0*t.unitcell_lengths[frame,:]

        # loop over all chromphores
        for chrom in chrom_idx:

           # re-center box at the COM of selected atoms
           xyz_chrom_raw = xyz_raw[chrom,:]
           com = getCOM(xyz_chrom_raw, masses[chrom])
           xyz = centerBox(xyz_raw, com, box)
           xyz_chrom = xyz[chrom,:]

           print (xyz_chrom)
           dddd

           # determine COM for all charge groups

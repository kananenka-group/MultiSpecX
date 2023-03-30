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
   start: int = 1
   elFmap: str  = "BaizXXX"
   freq_shift: float = 0.0

   def generateHamiltonian(self): 
     # create a system object
     s = System(self.itp,self.top,self.gro)
     self.atoms, self.molecules, self.atoms_in_mol, _, _ = s.read()

     # find indices of ester groups in each molecule
     chrom_start_idx, ester_list_idx, n_ester_mol = chromList(self.isotope_labels, self.ester_unit, self.atoms, self.atoms_in_mol)
     if not chrom_start_idx:
        print(f" Did not find any ester groups like this: {self.ester_unit}")
        sys.exit(" exiting...")
     else:
        print(f"       Found {len(chrom_start_idx)} ester chromophores: ")
        [print (f"       {id} in {molid[0]} ") for (id,molid) in zip(n_ester_mol,self.molecules) if id>0]


     cgS = chargeGroupSt(self.atoms)

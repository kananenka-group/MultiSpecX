import mdtraj as md
from ..system import *
from dataclasses import dataclass, field

@dataclass
class Ester:
   type: list = field(default_factory=lambda: ['dipole'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   ester_unit: list = field(default_factory=lambda: ['C','O','C','O'])
   gro:  str = "confout.gro"
   xtc:  str = "traj.xtc" 
   top:  str = "topology.top"
   nframes: int = 1
   start: int = 1
   elFmap: str  = "BaizXXX"
   freq_shift: float = 0.0

   def generateHamiltonian(self): 
     # create a system object
     s = System(self.itp)
     self.atoms, self.molecules, self.atoms_in_mol, _, _ = s.readITP()

     # find indices of ester units in the system

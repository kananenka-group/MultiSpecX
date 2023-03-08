import mdtraj as md
from ..system import *
from dataclasses import dataclass, field

@dataclass
class Ester:
   type: list = field(default_factory=lambda: ['dipole'])
   isotope_labels: list = field(default_factory=lambda: ['none'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   xtc:  str = "traj.xtc" 
   top:  str = "topology.top"
   nframes: int = 1
   start: int = 1
   elFmap: str  = "BaizXXX"
   freq_shift: float = 0.0

   def generateHamiltonian(self): 
     # create a system object
     s = System(self.itp)
     s.readITP()

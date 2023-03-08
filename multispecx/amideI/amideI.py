import mdtraj as md
from ..system import *
from dataclasses import dataclass, field

@dataclass
class AmideI:
   type: list = field(default_factory=lambda: ['dipole'])
   isotope_labels: list = field(default_factory=lambda: ['none'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   xtc:  str = "traj.xtc" 
   top:  str = "topology.top"
   nframes: int = 1
   start: int = 1
   nnFmap: str = "Jansen2006"
   nnCmap: str = "Jansen2006"
   Efmap: str  = "Wang2011"
   isotope_shift: float = -66.0

   def generateHamiltonian(self): 
     # create a system object
     s = System(self.itp)
     s.readITP()


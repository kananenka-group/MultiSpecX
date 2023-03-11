import mdtraj as md
from ..system import *
from dataclasses import dataclass, field

@dataclass
class Water:
   type: list = field(default_factory=lambda: ['dipole'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   gro:  str = "confout.gro" 
   xtc:  str = "traj.xtc" 
   top:  str = "topology.top"
   nframes: int = 1
   start: int = 1

   def generateHamiltonian(self): 
     # create a system object
     s = System(self.itp,self.top,self.gro)
     mol = s.read()
     print (mol,len(mol))


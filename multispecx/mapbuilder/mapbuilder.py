import sys
from dataclasses import dataclass, field

import mdtraj as md

from ..system import *

@dataclass
class Mapbuilder:
   itp: list = field(default_factory=lambda: ['topol.itp'])
   gro:  str = "confout.gro" 
   xyz:  str = "mol.xyz"
   xtc:  str = "traj.xtc"
   top:  str = "topology.top"
   method: str = "B3LYP"
   basis:  str = "6-311++G(d,p)"
   vib:  str = "normal"
   nframes: int = 1
   soft: str = "Gaussian"

   def createJobs(self):
     # create a system object
     s = System(self.itp,self.top,self.gro,self.xyz)
     self.atoms, self.molecules, self.atoms_in_mol, self.chem_labels = s.read(read_xyz=True)

     print (self.chem_labels)

     
     
     

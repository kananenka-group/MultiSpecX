import mdtraj as md
from ..system import *
from dataclasses import dataclass, field

@dataclass
class AmideI:
   type: list = field(default_factory=lambda: ['dipole'])
   isotope_labels: list = field(default_factory=lambda: ['none'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   amide_unit: list = field(default_factory=lambda: ['C','O','N','H'])
   gro:  str = "confout.gro" 
   xtc:  str = "traj.xtc" 
   top:  str = "topology.top"
   nframes: int = 1
   start: int = 1
   nnFmap: str = "Jansen2006"
   nnCmap: str = "Jansen2006"
   elFmap: str  = "Wang2011"
   freq_shift: float = -66.0

   def generateHamiltonian(self): 
     # create a system object
     s = System(self.itp,self.top,self.gro)
     self.atoms, self.molecules, self.atoms_in_mol = s.read()

     # find indices of peptide groups in each molecule.
     self.searchAmide()



   def searchAmide(self):
      print(f" >>>>> Searching amide I units defined as {self.amide_unit}")
      print (self.atoms_in_mol)

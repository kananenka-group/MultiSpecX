import mdtraj as md
from ..system import *
from dataclasses import dataclass, field
from .util import getIndex

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
     chrom_start_idx, n_amideI_mol = self.chromList()
     print(f"    Found {len(chrom_start_idx)} amide I chromophores: ")
     for nm, moln in zip(n_amideI_mol,self.molecules):
        if nm>0:
           print(f"  {nm} in {moln[0]} ") 

   def chromList(self):
     print(f" >>>>> Searching amide I units defined as {self.amide_unit} for residues {self.isotope_labels}")
     res_num = [ x[1] for x in self.atoms ]
     res_nam = [ x[2] for x in self.atoms ]
     atm_nam = [ x[3] for x in self.atoms ]
     ist=0
     cind=[]
     nas=[]
     for numa in self.atoms_in_mol:
        cout = getIndex(self.amide_unit, self.isotope_labels, res_num[ist:ist+numa], res_nam[ist:ist+numa], atm_nam[ist:ist+numa])
        nas.append(len(cout))
        if cout:
           cind.append(cout)
        ist+=numa
     return cind,nas
         
      
         

import sys
from dataclasses import dataclass, field

import mdtraj as md

from ..system import *
from .util import getIndex

@dataclass
class AmideI:
   type: list = field(default_factory=lambda: ['dipole'])
   isotope_labels: list = field(default_factory=lambda: ['none'])
   itp: list = field(default_factory=lambda: ['topol.itp'])
   amideI_unit: list = field(default_factory=lambda: ['C','O','N','H'])
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
     if not chrom_start_idx:
        print(f" Did not find any amide I groups like this: {self.amideI_unit}")
        sys.exit(" exiting...")
     else:
        print(f"       Found {len(chrom_start_idx)} amide I {self.isotope_labels} chromophores: ")
        [print (f"       {id} in {molid[0]} ") for (id,molid) in zip(n_amideI_mol,self.molecules) if id>0]

     # run loop over MD snapshots to get the Hamiltonian
     # 1. convert to numpy arrays for efficiency

   def chromList(self):
     print(f" >>>>> Searching amide I units defined as {self.amideI_unit} for residues {self.isotope_labels}")
     res_num = [ x[1] for x in self.atoms ]
     res_nam = [ x[2] for x in self.atoms ]
     atm_nam = [ x[3] for x in self.atoms ]
     ist=0
     cind=[]
     nas=[]
     for numa in self.atoms_in_mol:
        cout = getIndex(self.amideI_unit, self.isotope_labels, res_num[ist:ist+numa], res_nam[ist:ist+numa], atm_nam[ist:ist+numa])
        nas.append(len(cout))
        if cout:
           cind.append(cout)
        ist+=numa
     return cind,nas
         
      
         

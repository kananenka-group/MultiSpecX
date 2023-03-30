import sys
from dataclasses import dataclass, field

import mdtraj as md

from ..system import *
from .util import getIndex, chromList

@dataclass
class AmideI:
   type: list = field(default_factory=lambda: ['dipole'])
   isotope_labels: list = field(default_factory=lambda: [])
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
   freq_shift: float = 0.0

   def generateHamiltonian(self): 
     # create a system object
     s = System(self.itp,self.top,self.gro)
     self.atoms, self.molecules, self.atoms_in_mol, _, _ = s.read()

     # find indices of peptide groups in each molecule.
     chrom_idx, amideI_list_idx, n_amideI_mol = chromList(self.isotope_labels, self.amideI_unit, self.atoms, self.atoms_in_mol)
     if not chrom_idx:
        print(f" Did not find any amide I groups like this: {self.amideI_unit}")
        sys.exit(" exiting...")
     else:
        print(f"       Found {len(chrom_idx)} amide I chromophores: ")
        [print (f"       {id} in {molid[0]} ") for (id,molid) in zip(n_amideI_mol,self.molecules) if id>0]

     # run loop over MD snapshots to get the Hamiltonian
     # 1. convert to numpy arrays for efficiency


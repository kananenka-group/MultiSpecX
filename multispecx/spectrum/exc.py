from dataclasses import dataclass

@dataclass
class Sim():
   type: str  = "FTIR"
   hfile: str = "Energy.txt"
   dfile: str = "Dipole.txt"
   pfile: str = "Polarizability.txt"
   dt: float  = 0.02 # in ps

     

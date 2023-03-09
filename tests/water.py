import sys
sys.path.append("..")
from multispecx import water

s = water.Water(itp=["tip4p2005.itp"],top="water.top",gro="water.gro")
s.generateHamiltonian()



import sys
import os
sys.path.append("../..")
from multispecx import mapbuilder

s = mapbuilder.Mapbuilder(itp=["etoac.itp","tip4p.itp"],top="etoac.top",gro="prod.gro",method="B3LYP",basis="6-311++G(d,p)",vib="normal",nframes=2,xyz=["etoac.xyz","water.xyz"],software="Gaussian",dir=os.getcwd())
s.createJobs()


import sys
import os
sys.path.append("../..")
from multispecx import mapbuilder

work_dir = "/Users/akanane/Documents/Research/ester"
transform = [['center',3],['rotate',3,4,'x'],['rotate',3,1,'y']]

s = mapbuilder.Mapbuilder(itp=["etoac.itp","tip4p.itp"],top="etoac.top",gro="prod.gro",method="B3LYP",basis="6-311++G(d,p)",vib="normal",nframes=2,xyz=["etoac.xyz","water.xyz"],software="Gaussian",wdir=work_dir,xtc="prod.xtc",transform=transform,solv_ref_atom="OW",solu_ref_atom="C4")
s.createJobs()


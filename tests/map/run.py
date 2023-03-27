import sys
import os
sys.path.append("../..")
from multispecx import mapbuilder

work_dir = "/Users/akanane/Documents/Research/ester"
transform = [['center',3],['rotate',3,4,'x'],['rotate',3,1,'y']]
basis="Def2TZVP"
nframes=2
cut_off1=5.0

s = mapbuilder.Mapbuilder(itp=["etoac.itp","tip4p.itp"],top="etoac.top",gro="prod.gro",method="B3LYP",basis=basis,vib="normal",nframes=nframes,xyz=["etoac.xyz","water.xyz"],software="Gaussian",wdir=work_dir,xtc="prod.xtc",transform=transform,solv_ref_atom="OW",solu_ref_atom="C4",ncores=40,cut_off1=cut_off1)
s.createJobs()


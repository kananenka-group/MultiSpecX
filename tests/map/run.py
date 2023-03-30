import sys
import os
sys.path.append("../..")
from multispecx import mapbuilder

work_dir = "/Users/akanane/Documents/Research/ester"
transform = [['center',2],['rotate',2,6,'x'],['rotate',2,1,'y']]
basis="6-311++G(d,p)"
nframes=2
cut_off1=5.0
cut_off2=8.5

s = mapbuilder.Mapbuilder(itp=["etoac_tip3p.top"],top="etoac_tip3p.top",gro="etoac_tip3p.gro",method="B3LYP",basis=basis,vib="normal",nframes=nframes,xyz=["etoac.xyz","water.xyz"],software="Gaussian",wdir=work_dir,xtc="etoac_tip3p.xtc",transform=transform,solv_ref_atom="O",solu_ref_atom="O2",ncores=40,cut_off1=cut_off1,cut_off2=cut_off2)
s.createJobs()


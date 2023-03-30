import sys
sys.path.append("/Users/akanane/Documents/Software/MultiSpecX/")
from multispecx import spectrum
from multispecx import ester

itp_files = ["topol_Protein.itp","topol_Protein2.itp","topol_Protein3.itp","topol_Protein4.itp","topol_Other5.itp","topol_Other6.itp","topol_Other7.itp","tip3p.itp"]
ester_unit=['C','O','OE','CA']

s = ester.Ester(itp=itp_files,top="topol.top",gro="conf.gro",ester_unit=ester_unit)
s.generateHamiltonian()



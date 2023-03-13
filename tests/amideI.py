import sys
sys.path.append("..")
from multispecx import spectrum
from multispecx import amideI

itp_files = ["topol_Protein1.itp","topol_Protein2.itp","topol_Protein3.itp","topol_Protein4.itp","topol_DPPC.itp","topol_POTion.itp","topol_bindion.itp","topol_CLion.itp","spce.itp"]
amideI_unit=['C','O','N','H']

s = amideI.AmideI(itp=itp_files,top="topol.top",gro="s0s2s3.gro",isotope_labels=['VAL76'],amideI_unit=amideI_unit)
s.generateHamiltonian()



import sys
sys.path.append("..")
from multispecx import spectrum
from multispecx import amideI

s = amideI.AmideI(itp=["topol_Protein1.itp","topol_Protein2.itp","topol_Protein2.itp","topol_Protein4.itp"])
s.generateHamiltonian()



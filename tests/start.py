import sys
sys.path.append("..")
from multispecx import spectrum
from multispecx import amideI

s = amideI.AmideI(itp=["topol_Protein1.itp","topol_Protein2.itp","topol_Protein3.itp","topol_Protein4.itp","topol_DPPC.itp","topol_POTion.itp","topol_bindion.itp","spce.itp","topol_CLion.itp"],top="topol.top")
s.generateHamiltonian()



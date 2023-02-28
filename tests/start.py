import sys
sys.path.append("..")
from multispecx import spectrum

s = spectrum.Sim(type="2DIR",dt=0.02)

print (s.dt)
print (s)

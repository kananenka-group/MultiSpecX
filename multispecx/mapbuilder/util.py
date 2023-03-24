import numpy as np

def rotation_matrix(va, uvb):
   uva = va/np.linalg.norm(va)
   v = np.cross(uva, uvb)
   s = v/np.linalg.norm(v)
   c = np.dot(uva,uvb)
   Vx = np.array([[    0, -v[2],   v[1]],
                  [ v[2],     0,  -v[0]],
                  [-v[1],  v[0],      0]])
   Vx2 = Vx @ Vx
   fc = 1.0/(1.0 + c)
   Rot = np.eye(3) + Vx + fc*Vx2
   return Rot
 
def check_order(md, tr):
   for n in range(len(md)):
     if tr[n] not in md[n]:
        return False

   return True           

def check_4site_water(solv_atoms, Mlist):
   for ind, atom in enumerate(solv_atoms):
     if atom in Mlist:
        return ind
   return 0
   

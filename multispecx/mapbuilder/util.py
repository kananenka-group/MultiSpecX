import numpy as np

def rotation_matrix(va, vb):
   """
      Rotation matrix:
      you want to find a rotation matrix R that rotates unit vector va onto unit vector vb
      
      The rotation matrix is given by:
      R = I + [v]x + [v]x^2*(1-c)/s^2

      where:
      ------
      I is 3x3 identity matrix
      [v]x is the skew-symmetric cross-product matrix of v which is a cross product of va and vb
      c is a dot product of a and b (cosine of angle)
      s is the magnitude of a cross product between a and b
      [v]x is shown below (Vx) variable
      1-c/s^2 can also be written as 1-c/s^2 = 1-c/(1-c^2) = 1/(1+c)

   """
   uva = va/np.linalg.norm(va)
   uvb = vb/np.linalg.norm(vb)

   v = np.cross(uva, uvb)
   s = v/np.linalg.norm(v)
   c = np.dot(uva,uvb)

   if c > -1.0:
      Vx = np.array([[    0, -v[2],   v[1]],
                     [ v[2],     0,  -v[0]],
                     [-v[1],  v[0],      0]])
      Vx2 = Vx @ Vx
      fc = 1.0/(1.0 + c)
      Rot = np.eye(3) + Vx + fc*Vx2
   else:
      indl = np.where(uvb == 1.0)[0]
      Rot = np.eye(3)
      Rot[:,indl] *=-1
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
   

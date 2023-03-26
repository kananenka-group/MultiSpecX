import numpy as np

def PBC(a_xyz,ref_xyz,box):
   vv = ref_xyz - a_xyz
   sv = vv - minImage(vv, box) 
   return a_xyz + sv

def minImage(v, box):
   return v - np.multiply(box,np.rint(np.divide(v,box)))

def rotation_matrix(va, vb):
   """
      Generate rotation matrix:
      find a rotation matrix R that rotates vector va onto vector vb
      
      The rotation matrix is given by:
      R = I + [v]x + [v]x^2*(1-c)/s^2

      where:
      ------
      I       -  3x3 identity matrix
      [v]x    - the skew-symmetric cross-product matrix of v which is a cross product of va and vb
      c       - a dot product of a and b (cosine of angle)
      s       -  the magnitude of a cross product between a and b

      [v]x    is (Vx) variable below
      1-c/s^2 can also be written as 1-c/s^2 = 1-c/(1-c^2) = 1/(1+c)

      This approach works unless c=-1 which corresponds to
      uva and uvb pointing in exactly opposite direction. This special
      case is addressed below too
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
   

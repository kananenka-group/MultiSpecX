import numpy as np
import sys
from typing import List


def calcEf(atoms: List[int], xyz, xyz_ref, charges):
   """
      Caclulate electric fields on atoms and 
      optionally projections on selected vectors
   """
   assert xyz.shape[0] == len(charges), f" something is wrong in calcEf {xyz_ref.shape[0]} vs {len(charges)}"

   eF = np.zeros((len(atoms),3),dtype=np.float32)

   for atom in atoms:
      eFa = np.zeros((3))
      for ai in range(xyz.shape[0]):
         vij = xyz[ai,:] - xyz_ref[atom,:]
         vij *= 1.88973   # convert to a.u.
         dij = np.linalg.norm(vij)
         dij3 = dij**3
         eFa += vij*charges[ai]/dij3
      eF[atom,:] = np.copy(eFa) 
   # make it 1D for the electric field
   return np.reshape(eF,(len(atoms*3))), eF

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

   fc: float 
   indl: int 

   if c > -1.0:
      Vx = np.array([[    0, -v[2],   v[1]],
                     [ v[2],     0,  -v[0]],
                     [-v[1],  v[0],     0]],dtype=np.float32)
      Vx2 = Vx @ Vx
      fc = 1.0/(1.0 + c)
      Rot = np.eye(3,dtype=np.float32) + Vx + fc*Vx2
   else:
      indl = np.where(uvb == 1.0)[0]
      Rot = np.eye(3)
      Rot[:,indl] *=-1

   return Rot
 
def transformXYZ(transform, solu_xyz, solv_xyz):
   """
      Here input coordinates for the solute and solvent will be
      transformed

      throw in a warning message when coordiates are distorted by 
      1% (distance change w.r.t. reference atom)
   """
   thresh: float = 0.01
   solu_xyz_t = np.copy(solu_xyz)
   solv_xyz_t = np.copy(solv_xyz)

   for item in transform:
      # center frame on a given atom
      if item[0].lower() == "center":
         atom_center=item[1]-1
         xyz_shift = np.copy(solu_xyz_t[atom_center,:])
         for n in range(solu_xyz_t.shape[0]):
            solu_xyz_t[n,:] = np.subtract(solu_xyz_t[n,:],xyz_shift)
         for n in range(solv_xyz_t.shape[0]):
            solv_xyz_t[n,:] = np.subtract(solv_xyz_t[n,:],xyz_shift)
      # align w.r.t particular axis
      elif item[0].lower() == "rotate":
         atomn = item[1]-1
         atomp = item[2]-1
         va = np.subtract(solu_xyz_t[atomp,:],solu_xyz_t[atomn,:])
         if item[3].lower() == "z":
            vb = np.array([0.0,0.0,1.0],dtype=np.float32)
         elif item[3].lower() == "y":
            vb = np.array([0.0,1.0,0.0],dtype=np.float32)
         elif item[3].lower() == "x":
            vb = np.array([1.0,0.0,0.0],dtype=np.float32)
         else:
            sys.exit(f" Cannot recognize the rotation axis {item[3]} can only be 'x', 'y', or 'z'") 
         Rot = rotation_matrix(va,vb) 
         # rotate all atoms here ...
         for n in range(solu_xyz_t.shape[0]):
            solu_xyz_t[n,:] = Rot.dot(solu_xyz_t[n,:])
         for n in range(solv_xyz_t.shape[0]):
            solv_xyz_t[n,:] = Rot.dot(solv_xyz_t[n,:])
      else:
         sys.exit(" Cannot recognize this transformation operation {item}")

   # check the distance w.r.t ref atom before and after transformaton:
   for n in range(solu_xyz.shape[0]):
      bf = np.linalg.norm(np.subtract(solu_xyz[atom_center,:],solu_xyz[n,:]))
      af = np.linalg.norm(np.subtract(solu_xyz_t[atom_center,:],solu_xyz_t[n,:]))
      su_er = (bf-af)/bf if bf > 0.0 else 0.0
      if np.max(np.abs(su_er)) > thresh:
         print(f" Warning. Potential problem with coordinate transformation: ")
         print(" Before = ",solu_xyz[n,:])
         print(" After = ",solu_xyz_t[n,:])
         print(f" Error (w.r.t. solute reference atom {atom_center}) = {100*sv_er} %.")

   for n in range(solv_xyz.shape[0]):
      bf = np.linalg.norm(np.subtract(solu_xyz[atom_center,:],solv_xyz[n,:]))
      af = np.linalg.norm(np.subtract(solu_xyz_t[atom_center,:],solv_xyz_t[n,:]))
      sv_er = (bf-af)/bf if bf > 0.0 else 0.0
      if np.max(np.abs(sv_er)) > thresh:
         print(f" Warning. Potential problem with coordinate transformation: ")
         print(" Before = ",solv_xyz[n,:])
         print(" After = ",solv_xyz_t[n,:])
         print(f" Error (w.r.t. solute reference atom {atom_center}) = {100*sv_er} %.")

   return solu_xyz_t, solv_xyz_t

   def print_transform_info(self):
      """
         Just print transformation operations to be applied
      """ 
      if self.transform:
         print(f">>>>> Coordinate transformation:")

      for item in self.transform:
         if item[0].lower() == "center":
            print(f"      The frames will be translated to put atom {self.solute[3][item[1]-1]}({item[1]}) at the center of the box.")
         if item[0].lower() == "rotate":
            print(f"      The frames will be rotated such that vector {self.solute[3][item[1]-1]}({item[1]}) -> {self.solute[3][item[2]-1]}({item[2]}) will be aligned with the positive {item[3]} axis")

def AinF(xyz, atoms_exclude, xyz_ref, cut, cgS):

   atoms_include:list(int) = []
   for n in range(xyz.shape[0]):
      if np.linalg.norm(xyz[n,:]-xyz_ref) < cut:
         atoms = list(range(cgS[n],cgS[n+1]))
         atoms_include.extend(atoms)
   [ atoms_include.remove(atr) for atr in atoms_exclude ]
   
   return np.asarray(atoms_include,dtype=np.int32)

def getCOMChg(xyz, cgS, masses):
   """
      Returns center of masses of all charge groups
   """
   nchg = len(cgS)-1
   natoms = len(masses)

   com = np.zeros((nchg,3),dtype=np.float32)
   for n in range(nchg):
     gs = cgS[n]
     ge = cgS[n+1]
     mg = masses[gs:ge]
     cg = xyz[gs:ge,:]
     com[n,:] = getCOM(cg, mg)
   return com

def inBox(v, box):
   return v - np.multiply(box,np.floor(np.divide(v,box)))

def centerBox(xyz, ref, box):
   """
      Center simulation box at a given loc.
   """
   xyzc = np.copy(xyz)

   shift = np.subtract(ref,box/2.0)
   for n in range(xyzc.shape[0]):
      xyzc[n,:] = inBox(np.subtract(xyzc[n,:],shift),box)

   return xyzc

def getCOM(xyz, mass):
   """
      Calculate center of mass
   """
   com = np.zeros((3),dtype=np.float32)
   tmas = np.sum(mass)
   for n in range(3):
      com[n] = np.sum(np.multiply(xyz[:,n],mass))/tmas
   return com

def getInternalTransformXYZ(transform_in, atom_names, chrom_idx):
    print(f" >>>>> Building coordinate transformation matrix.")
    transform_out = []
    for chrom_idx_in in chrom_idx:
      this_chrom=[]
      for item in transform_in:
         if item[0] == 'center': 
            tloc = []
            tloc.append('center')
            loc = atom_names.index(item[1])
            tloc.append(int(loc+1)) #chrom_idx_in[loc])
         elif item[0] == 'rotate':
            tloc = []
            tloc.append('rotate')
            loc1 = atom_names.index(item[1])
            loc2 = atom_names.index(item[2])
            tloc.append(int(loc1)+1) #chrom_idx_in[loc1])
            tloc.append(int(loc2)+1) #chrom_idx_in[loc2]) 
            tloc.append(item[3])
         this_chrom.append(tloc)
      transform_out.append(this_chrom)

    return transform_out 

def chargeGroupSt(atoms):
   """
      This function determines were each charge group starts.
   """
   chg = [0]
   allCg = [ x[5] for x in atoms ]
   allMo = [ x[10] for x in atoms ] 
   
   for mol_idx in range(1,len(allMo)):
      if allMo[mol_idx] == allMo[mol_idx-1]:
         if allCg[mol_idx] == allCg[mol_idx-1]:
            continue
         else:
            chg.append(mol_idx)
      else:
         chg.append(mol_idx)
   return np.asarray(chg,dtype=int)

def chromList(isotope_labels, search_unit, atoms, atoms_in_mol):

     if isotope_labels:
        print(f" >>>>> Searching {search_unit} in {isotope_labels} residues")
     else:
        print(f" >>>>> Searching {search_unit} in all residues")

     res_num = [ x[2] for x in atoms ]
     res_nam = [ x[3] for x in atoms ]
     atm_nam = [ x[4] for x in atoms ]
     ist=0
     cind=[]
     aIind=[]
     nas=[]
  
     for numa in atoms_in_mol:
        cout, amd = getIndex(search_unit, isotope_labels, res_num[ist:ist+numa], res_nam[ist:ist+numa], atm_nam[ist:ist+numa])
        cout2 = [ [ y+ist for y in x ] for x in cout]
        amd2  = [ [ y+ist for y in x ] for x in amd]
        aIind.append(amd2)
        nas.append(len(cout))
        if cout:
           cind.append(cout2[0])
        ist+=numa
     return cind, aIind, nas

def getIndex(a_unit, isotope_labels, res_num, res_nam, atm_nam):
   """
      This function will return starting index of all chromophores.
   """
   natoms = len(atm_nam)
   nfind  = len(a_unit)

   chrom_list=[]
   amideI_list=[]
 
   for n in range(natoms-nfind):
      if a_unit == atm_nam[n:n+nfind]:
         resid = res_nam[n]+res_num[n]
         amideI_list.append(list(range(n,n+nfind)))
         if not isotope_labels:
            chrom_list.append(list(range(n,n+nfind)))
         else:
            if resid in isotope_labels:
                chrom_list.append(list(range(n,n+nfind)))
             
   return chrom_list, amideI_list
       

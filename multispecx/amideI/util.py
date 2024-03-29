import numpy as np
import warnings
from os import sys
from itertools import chain

#from colorama import Fore
# pip install colorama
from datetime import datetime
from typing import List

from .constants import *

def printChromList(chrom_idx, chrom_unit, molecules, n_chrom_mol, chrom_name):
   """
      Print chromophores found in the system
   """
   if not chrom_idx:
      print(f" Did not find any {chrom_name} groups like this: {chrom_unit}")
      sys.exit(" exiting...")
   else:
      print(f"       Found {len(chrom_idx)} ester chromophores: ")
      [print (f"       {id} in {molid[0]} ") for (id,molid) in zip(n_chrom_mol,molecules) if id>0]

def printDipole(path,Dipole) -> None:
   """
      Print transition dipoles into Dipole.txt file
   """
   outFile = path/"Dipole.txt"
   row_n = np.arange(0,Dipole.shape[0],dtype=np.int32)
   fmt = ["%.7f"]*Dipole.shape[1]
   fmt.insert(0,"%i")
   np.savetxt(outFile,np.column_stack((row_n,Dipole)),fmt=fmt)
   print(f" >>>>> Transition dipole moments have been saved to Dipole.txt file.")

def printEnergy(path,Energy) -> None:
   """
      Print upper triangle Hamiltonian into Energy.txt file
   """ 
   outFile = path/"Energy.txt"
   n: int = Energy.shape[1]*(Energy.shape[1]+1)//2
   Energy_out = np.zeros((Energy.shape[0],n),dtype=np.float32)

   for frame in range(Energy.shape[0]):
      Energy_out[frame,:] = Energy[frame,:,:][np.triu_indices(Energy.shape[1])]
   
   row_n = np.arange(0,Energy.shape[0],dtype=np.int32)
   fmt = ["%.7f"]*n 
   fmt.insert(0,"%i")
   np.savetxt(outFile,np.column_stack((row_n,Energy_out)),fmt=fmt)
   print(f" >>>>> Hamiltonian has been saved to Energy.txt file.")

def minImage(v, box):
   return v - np.multiply(box,np.round(np.divide(v,box)))

def printDT(pt):
   print(f" >>>>> Simulation {pt}: {datetime.today():%B %d, %Y %H:%M:%S}") 

def TDC(tdv_i, tdv_j, tdp_i, tdp_j, box) -> float:
   """
      Calculate transition dipole coupling in cm-1

      Input:
      ------
      tdv_i and tdv_j - transition dipole moments in D*A^{-1}*a.m.u.^{-1/2}
      tdp_i and tdp_j - location of dipoles in A    
      box             - box dimensions in A
  
      Use scaling factor 84861.9/1650.0 which converts to cm-1
         from: J. Phys. Chem. B 2006, 110, 3362-3374 
          and  J. Phys. Chem. B 2011, 115, 3713–3724
 
   """
   scale: float = 84861.9/1650.0  
   rij = minImage(np.subtract(tdp_j,tdp_i),box)
   dij = 1.0/np.linalg.norm(rij)
   dij3 = dij**3
   dij5 = dij3*dij*dij
   omega = dij3*np.dot(tdv_i,tdv_j) - 3.0*dij5*np.dot(tdv_i,rij)*np.dot(tdv_j,rij)
   return omega*scale

def calcEf(atoms: List[int], projections, xyz, xyz_ref, charges):
   """
      Caclulate electric fields on atoms and 
      projections on selected bonds.

      Input
      -----
      atoms       - List of atoms for which electric field is to be calculated
      projections - List of pairs of atoms connected by a bond to project electric field on
      xyz         - Cartesian coordinates of all atoms in the environemnt, in A
      xyz_ref     - Cartesian coordinates of the reference group of atoms, in A
      charges     - Atomic charges in a.u.
   """
   if xyz.shape[0] != len(charges):
      raise ValueError(f""" Wrong dimensions in calcEf xyz shape {xyz_ref.shape[0]} vs 
                            charges shape {len(charges)}""")

   eF = np.zeros((len(atoms),3),dtype=np.float32)

   #for atom in atoms:
   #   eFa = np.zeros((3))
   #   for ai in range(xyz.shape[0]):
   #      vij = ATOAU*(xyz[ai,:] - xyz_ref[atom,:])  #units=a.u. because map requires this
   #      dij = np.linalg.norm(vij)
   #      dij3 = dij**3
   #      eFa += vij*charges[ai]/dij3
   #      #print (" atom ",vij*charges[ai]/dij3)
   #   eF[atom,:] = np.copy(eFa) 
   
   for atom in atoms:
      #units=a.u. because map requires this
      vij = np.array([ ATOAU*(xyz[x,:] - xyz_ref[atom,:]) for x in range(xyz.shape[0]) ])
      qvij= np.array([ charges[x]*vij[x,:] for x in range(vij.shape[0]) ])  
      dij = np.array([ np.linalg.norm(vij[x,:]) for x in range(vij.shape[0]) ])
      dij3 = dij**3
      eF[atom,:] = np.sum(np.array([ np.divide(qvij[x,:],dij3[x]) for x in range(qvij.shape[0])]),axis=0)

   # projections:
   eFp = np.zeros((len(projections)),dtype=np.float32)
   for proja in projections:
      ...

   return np.reshape(eF,(len(atoms*3))), eFp

def rotFlip(uvb):
   """
      In some cases rotation needed is 180o flip
      handle it here
   """
   indl = np.where(uvb == 1.0)[0]
   Rot = np.eye(3,dtype=np.float32)
   Rot[:,indl] *=-1.0
   return Rot

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
   fc: float 

   uva = va/np.linalg.norm(va)
   uvb = vb/np.linalg.norm(vb)

   v = np.cross(uva, uvb)

   norm_v = np.linalg.norm(v)
   if np.abs(norm_v) < 1.0e-7:
      if norm_v > 0.0:
         # two vectors are parallel, no rotation is needed
         warnings.warn(f" Norm of uva x uvb is small {norm_v}. Vectors are parallel? {uva} and {uvb}. No rotation needed.")
         Rot = np.eye(3,dtype=np.float32)
         return Rot
      else:
         # two vectors are antiparallel flip would be needed
         warnings.warn(f" Norm of uva x uvb is small {norm_v}. Vectors are antiparallel? {uva} and {uvb}. Flipping...")
         return rotFlip(uvb)
   else:
      s = v/norm_v
      c = np.dot(uva,uvb)
   
      # work on this below, make return statement for c=-1
      if c > -1.0:
         Vx = np.array([[    0, -v[2],   v[1]],
                        [ v[2],     0,  -v[0]],
                        [-v[1],  v[0],     0]],dtype=np.float32)
         Vx2 = Vx @ Vx
         fc = 1.0/(1.0 + c)
         Rot = np.eye(3,dtype=np.float32) + Vx + fc*Vx2
      else:
         Rot = rotFlip(uvb)

   return Rot

def CartDir(inp: str):
   """
      Return unit vector in specified direction.

      Input:
      ------
      str - "x", "y", or "z"
   
   """
   match inp:
      case 'x':
         return np.array([1.0,0.0,0.0],dtype=np.float32)
      case 'y':
         return np.array([0.0,1.0,0.0],dtype=np.float32)
      case 'z':
         return np.array([0.0,0.0,1.0],dtype=np.float32)
      case _:
         raise ValueError (f" Cannot recognize the rotation axis {inp} can only be 'x', 'y', or 'z'") 

def transformXYZ(transform, solu_xyz, solv_xyz):
   """
      Transform solute and solvent coordinates according to
      given transformations

      Input:
      ------
      transform - nested list of all required transformations (shifts and rotations)
      solu_xyz  - Cartesian coordinates of the solute molecule, in A
      solv_xyz  - Cartesian coordinates of the solvent molecules, in A 

      Notes:
      ------
      Check the distance between each solute/solvent atom and the reference solute
      atom and throw in a warning when coordiates are changes by more than thresh
      (currently set to 1.0%)
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
         vb = CartDir(item[3].lower())
         Rot = rotation_matrix(va,vb) 
         # rotate all atoms here ...
         for n in range(solu_xyz_t.shape[0]):
            solu_xyz_t[n,:] = Rot.dot(solu_xyz_t[n,:])
         for n in range(solv_xyz_t.shape[0]):
            solv_xyz_t[n,:] = Rot.dot(solv_xyz_t[n,:])
      else:
         raise ValueError(f" Cannot recognize this transformation operation {item}")

   # check the distance w.r.t ref atom before and after the transformaton
   # and throw in a warning if distance changes by more than the defined threshold
   for n in range(solu_xyz.shape[0]):
      bf = np.linalg.norm(np.subtract(solu_xyz[atom_center,:],solu_xyz[n,:]))
      af = np.linalg.norm(np.subtract(solu_xyz_t[atom_center,:],solu_xyz_t[n,:]))
      su_er = (bf-af)/bf if bf > 0.0 else 0.0
      if np.max(np.abs(su_er)) > thresh:
         print(" Before = ",solu_xyz[n,:])
         print(" After = ",solu_xyz_t[n,:])
         su_er*=100.0
         print(f" {su_er:.2f} % ")
         warnings.warn(f" Coordinate transformation error (w.r.t. solute reference atom {atom_center}) = {su_er:.2f} %.")

   for n in range(solv_xyz.shape[0]):
      bf = np.linalg.norm(np.subtract(solu_xyz[atom_center,:],solv_xyz[n,:]))
      af = np.linalg.norm(np.subtract(solu_xyz_t[atom_center,:],solv_xyz_t[n,:]))
      sv_er = (bf-af)/bf if bf > 0.0 else 0.0
      if np.max(np.abs(sv_er)) > thresh:
         print(" Before = ",solv_xyz[n,:])
         print(" After = ",solv_xyz_t[n,:])
         sv_er*=100.0
         print(f" {sv_er:.2f} % ")
         warnings.warn(f" Coordinate transformation error (w.r.t. solute reference atom {atom_center}) = {sv_er:.2f} %.")

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

def include_CG_atoms(xyz, xyz_ref, cut, cgS):
   """
      Make a list of atoms that are close enough to be
      included into electric field calculation
   """
   atoms_include:list(int) = []

   for n in range(xyz.shape[0]):
      if np.linalg.norm(xyz[n,:]-xyz_ref) < cut:
         # range below is iterator use it like the one to get atoms...
         atoms = list(range(cgS[n],cgS[n+1]))
         atoms_include.extend(atoms)
   
   return list(atoms_include)

def exclude_atoms_from_list(atoms, atoms_exclude) -> List[int]:
   """
      Exclude atoms of the chromophore from atoms list
   """
   return list(filter(lambda x: x not in chain(*atoms_exclude), atoms))
   
def getCOMChg(xyz, cgS, masses):
   """
      Returns center of masses of all charge groups
   """
   nchg = len(cgS)-1

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
      Center simulation box at a given location
   """
   xyzc = np.copy(xyz)

   shift = np.subtract(ref,box/2.0)
   for n in range(xyzc.shape[0]):
      xyzc[n,:] = inBox(np.subtract(xyzc[n,:],shift),box)

   #xyzc = np.array([ inBox((xyzc[n,:]-shift),box) for n in range(xyzc.shape[0]) ])

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
   """
      Translate input transformation instructions into
      internal transformation procedure
   """
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
            tloc.append(int(loc1)+1) 
            tloc.append(int(loc2)+1) 
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
   # add the last atom to terminate the count
   chg.append(len(atoms))

   return np.array(chg,dtype=np.int32)

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
      This function will return the starting index of all chromophores.
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
       

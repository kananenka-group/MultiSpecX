import numpy as np

def getCOMChg(xyz, cgS, masses):
   """
      Returns center of masses of all charge groups
   """
   nchg = len(cgS)-1
   natoms = len(masses)

   com = np.zeros((nchg,3))
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

   shift = np.subtract(ref,box/2)
   for n in range(xyzc.shape[0]):
      xyzc[n,:] = inBox(np.subtract(xyzc[n,:],shift),box)

   return xyzc

def getCOM(xyz, mass):
   """
      Calculate center of mass
   """
   com = np.zeros((3))
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
            tloc.append(chrom_idx_in[loc])
         elif item[0] == 'rotate':
            tloc = []
            tloc.append('rotate')
            loc1 = atom_names.index(item[1])
            loc2 = atom_names.index(item[2])
            tloc.append(chrom_idx_in[loc1])
            tloc.append(chrom_idx_in[loc2]) 
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
       

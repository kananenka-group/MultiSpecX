import numpy as np

def chargeGroupSt(atoms):
   """
      This function determines were each charge group starts.
   """
   chg = [0]
   allCg = [ x[4] for x in atoms ]
   allMo = [ x[9] for x in atoms ] 
   
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

     res_num = [ x[1] for x in atoms ]
     res_nam = [ x[2] for x in atoms ]
     atm_nam = [ x[3] for x in atoms ]
     ist=0
     cind=[]
     aIind=[]
     nas=[]
     for numa in atoms_in_mol:
        cout, amd = getIndex(search_unit, isotope_labels, res_num[ist:ist+numa], res_nam[ist:ist+numa], atm_nam[ist:ist+numa])
        aIind.append(amd)
        nas.append(len(cout))
        if cout:
           cind.append(cout)
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
         amideI_list.append(n)
         if not isotope_labels:
            chrom_list.append(n)
         else:
            if resid in isotope_labels:
                chrom_list.append(n)
             
   return chrom_list, amideI_list
       

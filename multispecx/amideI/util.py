
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
         if resid in isotope_labels:
            chrom_list.append(n)
             

   return chrom_list, amideI_list
       

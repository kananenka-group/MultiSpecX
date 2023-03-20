
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
   

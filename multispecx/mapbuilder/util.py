
def check_order(md, tr):
    
    for n in range(len(md)):
      if tr[n] not in md[n]:
         return False

    return True           

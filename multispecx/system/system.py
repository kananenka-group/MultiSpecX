from dataclasses import dataclass, field

@dataclass
class System:
   itp_files: list = field(default_factory=lambda: ['topol.itp'])

   # read itp files:
   def readITP(self):
      print(" reading ITP files ")
      for itp_file in self.itp_files:
         with open(itp_file,"r") as f:
            for line in f:
               if not line.lstrip().startswith(';'):
                   print (line)

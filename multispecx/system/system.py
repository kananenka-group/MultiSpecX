from dataclasses import dataclass, field
import numpy as np
import sys

@dataclass
class System:
   itp_files: list = field(default_factory=lambda: ['topol.itp'])
   top_file:  str = "topol.top"
   gro_file:  str = "confout.gro"

   # reading gromacs files to learn about the system
   def read(self):
      # reading topology file
      self.molecules, self.molnum = self.readTOP()
      print (f" Found following molecules in {self.top_file} file")
      [print(f" {line[0]}  {line[1]}") for line in self.molnum]

      # reading itp files
      self.molecule_list, self.molecule_find = self.readITP()
      if np.any(self.molecule_find==False):
         print(" The following systems were not found in *.itp files:")
         x = np.transpose(np.argwhere(self.molecule_find==False))[0]
         [print (f" {self.molecules[i]}") for i in x]
         sys.exit(" exiting...")

      # reading GRO file
      self.natoms, self.system = self.readGRO()
      print (f" Total number of atoms: {self.natoms}")

      # do the match
      out = self.match()
      
   
   def match(self):
      """
         match all the information and output for each atom in the
         configuration its charge, charge group, mass, and a list to which
         molecule this atom belongs
      """
      return 0

   def readGRO(self):
      """
         Reading GROMACS configuration files to be able to match atoms
         in the system and topology.
      """
      system=[]
      with open(self.gro_file,"r") as f:
         print (f" >>>>> Reading configuration file: {self.gro_file}")
         next(f)
         line=f.readline()
         natoms=int(line.split()[0])
         for n in range(natoms):
             line=f.readline()
             system.append([line[0:5].strip(),line[5:10].strip(),line[10:15].strip()])
      return natoms, system

   def readTOP(self):
      """
         Reading GROMACS *.top file to learn what to read from
         *.itp files. 
      """
      print (f" >>>>> Reading topology file: {self.top_file}")
      record=False
      molecules=[]
      mol_numbers=[]
      with open(self.top_file,"r") as f:
         for line in f: 
            lines=line.strip()
            if not lines.startswith(';'): #or not lines.startswith('#'):
               if record:
                  if not lines:
                     break
                  else:
                     addM = []
                     addM.append(line.split()[0])
                     addM.append(line.split()[1])
                     mol_numbers.append(addM)
                     molecules.append(line.split()[0])
               if lines=="[ molecules ]":
                  record=True
      return molecules, mol_numbers

   def readITP(self):
      """
         Read *.itp files
      """
      molecule_list = [[] for i in range(len(self.molecules))]
      molecule_find = np.full(len(self.molecules), False)

      for itp_file in self.itp_files:
         with open(itp_file,"r") as f:
            print(f" >>>>> Reading itp file: {itp_file} ")
            record=False
            record2=False
            for line in f:
               lines=line.strip()
               if not lines.startswith(';'): #or not lines.startswith('#'):
                  if record2:
                     if not lines:
                        record=False
                        record2=False
                     else:
                        atom=[]
                        if len(lines.split())==7:
                          atom.append(lines.split()[2:7])
                        elif len(lines.split())>7:
                          atom.append(lines.split()[2:8])
                        else:
                          sys.exit(f" Cannot read this line{lines} in {itp_file} wrong number of entries.")
                        molecule_list[index].append(atom[0])
                  if len(lines.split())>1 and lines.split()[0] in self.molecules:
                     record=True
                     index=self.molecules.index(lines.split()[0])
                     molecule_find[index]=True
                  if record and lines=="[ atoms ]":
                     record2=True 
                      
      return molecule_list, molecule_find  

from dataclasses import dataclass, field
import numpy as np
import sys

@dataclass
class System:
   # defaults:
   itp_files: list = field(default_factory=lambda: ['topol.itp'])
   top_file:  str = "topol.top"
   gro_file:  str = "confout.gro"
   xyz_files:  str = field(default_factory=lambda: ['solute.xyz','solvent.xyz'])

   hydrogen_list = ['HW1','HW2','HW','HWT4']
   oxygen_list   = ['OW','OWT4']
   msite_list    = ['MW','MWT4']

   # reading gromacs files to learn about the system
   def read(self,read_xyz=False):
      # reading topology file
      self.molecules, self.molnum = self.readTOP()
      print (f"       Found following molecules in {self.top_file} file")
      [print(f"       {line[0]}  {line[1]}") for line in self.molnum]

      # reading itp files
      self.molecule_list, self.molecule_find = self.readITP()
      if np.any(self.molecule_find==False):
         print(" The following systems were not found in *.itp files:")
         x = np.transpose(np.argwhere(self.molecule_find==False))[0]
         [print (f" {self.molecules[i]}") for i in x]
         sys.exit(" exiting...")

      # reading GRO file
      self.natoms, self.system = self.readGRO()
      print (f"       Total number of atoms: {self.natoms}")

      # do the match
      atoms_out, atoms_in_mol = self.match()
      assert self.natoms == len(atoms_out), f" Assignment problem: total atoms {self.natoms}, assigned = {len(mol_out)}"
 
      print (" >>>>> Done reading GROMACS files.")
 
      atom_labels=[]
      if read_xyz:
         atom_labels = self.readXYZ()

      return atoms_out, self.molnum, atoms_in_mol, atom_labels, self.molecule_list
   
   def match(self):
      """
         Match all the information and output for each atom in the
         configuration its charge, charge group, mass, and a list to which
         molecule this atom belongs.
         
         This match function is based on GROMACS rule that states:
         'The order of the blocks of molecule types and the numbers of 
          such molecules must match the coordinate file that accompanies 
          the topology when supplied to grompp.' 
          [from GROMACS manual]:
          https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html

          This function outputs a list containing the following information 
          about each atom:

          0. atom number
          1. atom type
          2. residue number
          3. molecule name
          4. atom name
          5. charge group
          6. atomic charge
          7. mass
          8. molecule name
          9. id of a molecule this atoms belongs to (w.r.t the total number of molecules)

      """
      #hydrogen_list = ['HW1','HW2','HW','HWT4']
      #oxygen_list   = ['OW','OWT4']
      #msite_list    = ['MW','MWT4']
      h_mass        = 1.008
      o_mass        = 15.999
      m_mass        = 0.0

      add_o = add_h = add_m = False
      h_labels=[]
      o_labels=[]
      m_labels=[]
      atoms_in_mol=[]
     
      print (f" >>>>> Generating atom info")
      start=0 

      atoms_info=[]
      this_mol=0
      for mol_id, (mol, atoms) in enumerate(zip(self.molnum, self.molecule_list)):
         mol_name = mol[0]
         mol_nums = int(mol[1])
         for s in range(mol_nums):
            for ind, atom in enumerate(atoms):
               if atom[2:4] == self.system[start+ind][1:3]:
                  ats = atom[1:7].copy()
                  ats.insert(0,self.system[start+ind][0])
                  ats.insert(0,self.system[start+ind][3])
                  #if len(ats)==5:
                  if len(ats)==7:
                     val=ats[-3]
                     if val in self.hydrogen_list:
                        ats.append(h_mass)
                        add_h=True
                        h_labels.append(val)
                     elif val in self.oxygen_list: 
                        ats.append(o_mass)
                        add_o=True
                        o_labels.append(val)
                     elif val in self.msite_list:
                        ats.append(m_mass)
                        add_m=True
                        m_labels.append(val)
                  ats.append(mol[0])
                  ats.append(mol_id)
                  ats.append(this_mol)
                  atoms_info.append(ats)
            this_mol+=1
            atoms_in_mol.append(len(atoms_info)-start)
            start=len(atoms_info)

      if add_h:
          hf = list(set(h_labels))
          print(f"       {len(h_labels)} water hydrogen atoms were found, types={hf} and assigned mass {h_mass}")
      if add_o:
          of = list(set(o_labels))
          print(f"       {len(o_labels)} water oxygen atoms were found, types={of} and assigned mass {o_mass}")
      if add_m:
          mf = list(set(m_labels))
          print(f"       {len(m_labels)} water m-site atoms were found, types={mf} and assigned mass {m_mass}") 

      return atoms_info, atoms_in_mol
             

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
             system.append([line[0:5].strip(),line[5:10].strip(),line[10:15].strip(),line[15:20].strip()])
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
                          atom.append(lines.split()[1:7])
                        elif len(lines.split())>7:
                          atom.append(lines.split()[1:8])
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

   def readXYZ(self):
      """
         Reading xyz files
      """
      molecules=[]
      for xyz_file in self.xyz_files:
         print(f" >>>>> Reading xyz file: {xyz_file} ")
         with open(xyz_file,"r") as f:
            atoms=[]
            natoms = int(f.readline().split()[0])
            xyz = np.zeros((natoms,3))
            f.readline()
            for na in range(natoms):
               atom_line = f.readline().split()
               atoms.append(atom_line[0])
            molecules.append(atoms)
      return molecules


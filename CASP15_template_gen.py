#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:24:43 2022

@author: aaron
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 18:33:53 2022
@author: Aaron Sweeney
File writes out a template for an sdf file given a smiles string
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import sys
def get_zeroed_mol(smi = None,inchi = None, removeHs = False):
  ''' Get a mol with all atoms in position 0.0, 0.0, 0.0'''
  if not inchi is None:
    mol = Chem.MolFromInchi(inchi)
  elif not smi is None:
    mol = Chem.MolFromSmiles(smi)
  elif smi is None and inchi is None:
    print('SMILES and Inchi input are both None')
    sys.exit()
  mol = Chem.AddHs(mol)
  AllChem.EmbedMolecule(mol)
  if removeHs:
    mol = Chem.RemoveHs(mol)
  #return mol
  conf = mol.GetConformer()
  for i in range(mol.GetNumAtoms()):
    conf.SetAtomPosition(i,Point3D(0.0, 0.0, 0.0))
  return mol
def convert(option ):
  if option.upper() == 'TRUE':
    return True
  elif option.upper() == 'FALSE':
    return False
  else:
    print('Unknown option input { }'.format(option))
    sys.exit()
if __name__ == '__main__':
  if len(sys.argv) <= 1:
    print('\nUsage template_gen.py <Options>')
    print('Options: \n\t-n : Number of copies to write to a file, default = 1')
    print('\t-s : Use smiles string, use -s <smiles string>')
    print('\t-i : Use Inchi string, use -s <Inchi string>')
    print('\t-h : remove Hs upon writing template, default = False')
    print('\t-f : file to write to, default = ./TEMPLATE_<smiles string>.sdf\n')
    sys.exit()
  else:
    #smi = Chem.CanonSmiles(sys.argv[1])
    if '-s' in sys.argv and '-i' in sys.argv:
      print('Both SMILES string and Inchi used as input pick one!!')
      sys.exit()
    if '-s' in sys.argv:
      # n = number of copies in file
      smi = Chem.CanonSmiles(sys.argv[sys.argv.index('-s') + 1])
      use_smi = True
    else:
      use_smi = False
    if '-n' in sys.argv:
      # n = number of copies in file
      n = int(sys.argv[sys.argv.index('-n') + 1])
    else:
      n = 1
    if '-h' in sys.argv:
      removeHs = convert(sys.argv[sys.argv.index('-h') + 1].upper())
    else:
      #Are we giving them H's in the template file?
      removeHs = False
    if '-f' in sys.argv:
      fn = sys.argv[sys.argv.index('-f') + 1]
    else:
      #Are we giving them H's in the template file?
      fn = 'TEMPLATE_{}.sdf'.format(smi)
    if '-i' in sys.argv:
      inchi = sys.argv[sys.argv.index('-i') + 1]
      use_inchi = True
    else:
      #Are we giving them H's in the template file?
      use_inchi = False
  if use_inchi:
    mol = get_zeroed_mol(inchi = inchi, removeHs = removeHs)
  elif use_smi:
    mol = get_zeroed_mol(smi, removeHs = removeHs)
  else:
    print('No SMILES or Inchi string detected!!')
    sys.exit()
  #mol = get_zeroed_mol(smi, removeHs = removeHs)
  writer = Chem.SDWriter(fn)
  for _ in range(n):
    writer.write(mol)

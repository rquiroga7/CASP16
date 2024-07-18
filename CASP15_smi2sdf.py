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
import sys

def MolFromSmiles(smi, removeHs=False):
    ''' Get a mol from SMILES '''
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print('Invalid SMILES string')
        sys.exit()
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    if removeHs:
        mol = Chem.RemoveHs(mol)
    return mol

def convert(option):
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
        if '-s' in sys.argv and '-i' in sys.argv:
            print('Both SMILES string and Inchi used as input pick one!!')
            sys.exit()
        if '-s' in sys.argv:
            smi = Chem.CanonSmiles(sys.argv[sys.argv.index('-s') + 1])
            use_smi = True
        else:
            use_smi = False
        if '-n' in sys.argv:
            n = int(sys.argv[sys.argv.index('-n') + 1])
        else:
            n = 1
        if '-h' in sys.argv:
            removeHs = convert(sys.argv[sys.argv.index('-h') + 1].upper())
        else:
            removeHs = False
        if '-f' in sys.argv:
            fn = sys.argv[sys.argv.index('-f') + 1]
        else:
            fn = 'TEMPLATE_{}.sdf'.format(smi)
        if '-i' in sys.argv:
            inchi = sys.argv[sys.argv.index('-i') + 1]
            use_inchi = True
        else:
            use_inchi = False

    if use_inchi:
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            print('Invalid InChI string')
            sys.exit()
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        if removeHs:
            mol = Chem.RemoveHs(mol)
    elif use_smi:
        mol = MolFromSmiles(smi, removeHs=removeHs)
    else:
        print('No SMILES or Inchi string detected!!')
        sys.exit()

    writer = Chem.SDWriter(fn)
    for _ in range(n):
        writer.write(mol)


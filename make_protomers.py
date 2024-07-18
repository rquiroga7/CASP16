import sys
from dimorphite_dl import DimorphiteDL
from rdkit import Chem
from rdkit.Chem import Draw
import os

# Initialize DimorphiteDL with the desired parameters
dimorphite_dl = DimorphiteDL(
    min_ph=7.2,
    max_ph=7.6,
    max_variants=1,
    label_states=True,
    pka_precision=1.0
)

# Default values
output_file_base = 'protomers'
sdf_output_file_base = 'protomers'

# Check command line arguments
if '-s' not in sys.argv:
    print('Usage: python script.py -s <SMILES> [-f <output_file_base>]')
    sys.exit(1)

# Get the SMILES string from command line argument
smiles_index = sys.argv.index('-s') + 1
smiles_string = sys.argv[smiles_index]

# Check if an output file base name is specified
if '-f' in sys.argv:
    file_index = sys.argv.index('-f') + 1
    output_file_base = sys.argv[file_index]
    sdf_output_file_base = output_file_base

# Protonate the molecule using DimorphiteDL
prot_lst = dimorphite_dl.protonate(smiles_string)

# List to store RDKit molecules
mol_lst = []

# Process each protonated SMILES
for idx, smiles in enumerate(prot_lst):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f'Error processing SMILES: {smiles}')
        continue
    
    # Calculate and print formal charge
    charge = Chem.rdmolops.GetFormalCharge(mol)
    print(f"SMILES: {smiles}, Charge: {charge}")

    # Append molecule to list
    mol_lst.append(mol)

    # Write molecule to SDF file with index
    sdf_output_file = f'{sdf_output_file_base}_{idx}.sdf'
    w = Chem.SDWriter(sdf_output_file)
    w.write(mol)
    w.close()
    print(f'Molecule {idx} saved to SDF file: {sdf_output_file}')

# Draw molecules
img = Draw.MolsToGridImage(mol_lst, molsPerRow=2, subImgSize=(500, 500))

# Save to PNG
output_file = f'{output_file_base}.png'
img.save(output_file)
print(f'Image saved to: {output_file}')


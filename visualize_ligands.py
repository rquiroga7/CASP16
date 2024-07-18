from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Read SMILES from file
with open('tmp.smi', 'r') as f:
    smiles = f.readlines()

# Convert SMILES to RDKit Molecule objects
mols = [Chem.MolFromSmiles(smile) for smile in smiles]

# Generate a 2D coordinates for visualization
for mol in mols:
    Chem.rdDepictor.Compute2DCoords(mol)

# Draw molecules
img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(400, 400))

# Save to PNG
img.save('molecules.png')

# Display the image
#img = mpimg.imread('molecules.png')
#imgplot = plt.imshow(img)
#plt.show()

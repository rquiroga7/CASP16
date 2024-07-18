from rdkit import Chem
from rdkit.Chem import Draw
import glob

# Function to read all SDF files in a folder
def read_sdf_files(folder):
    sdf_files = glob.glob(folder + '/*.sdf')
    return sdf_files

# Function to generate 2D representations of molecules and save as PNG
def generate_png_from_sdf(sdf_files, output_file):
    mols = []
    for sdf_file in sdf_files:
        suppl = Chem.SDMolSupplier(sdf_file)
        for mol in suppl:
            if mol:
                mols.append(mol)

    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200))
    img.save(output_file)
    print(f'Image saved to: {output_file}')

# Main function
if __name__ == '__main__':
    # Replace with the path to your folder containing .sdf files
    folder_path = '.'

    # Read all .sdf files in the folder
    sdf_files = read_sdf_files(folder_path)

    # Specify the output file for the PNG image
    output_png = 'output.png'

    # Generate the PNG image from SDF files
    generate_png_from_sdf(sdf_files, output_png)


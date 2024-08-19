# Append the PocketOptimizer Code
import sys
sys.path.append('/home/rquiroga/github/pocketoptimizer')

# Import the pocketoptimizer module
import pocketoptimizer as po

project_dir="/home/rquiroga/github/CASP16/L1000_exper/PDB"
# Initialize a new design pipeline
design = po.DesignPipeline(work_dir=project_dir,         # Path to working directory containing scaffold and ligand subdirectory
                           ph=7,                         # pH used for protein and ligand protonation
                           forcefield='amber_ff14SB',    # forcefield used for all energy computations (Use Amber as it is better tested!)
                           ncpus=8)                      # Number of CPUs for multiprocessing


design.parameterize_ligand(
input_ligand='L1012_ligand_1.mol2', # Input ligand structure file could be .mol2/.sdf
addHs=True                              # Whether to add hydrogen atoms to the input structure
)

design.prepare_protein(
    protein_structure='4k2y_protein.pdb',  # Input PDB
    keep_chains=['A'],  # Specific protein chain to keep
    backbone_restraint=True, # Restrains the backbone during the minimization
    cuda=False,              # Performs minimization on CPU instead of GPU
    discard_mols=[]          # Special molecules to exclude. Per default everything, but peptides have to be defined manually
    )
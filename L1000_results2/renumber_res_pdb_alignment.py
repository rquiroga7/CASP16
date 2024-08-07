from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO
from Bio import Align
from Bio.Align import MultipleSeqAlignment
import os
import sys
import warnings
warnings.filterwarnings("ignore")

# Read in the PDB file
pdb_file = sys.argv[1]
#pdb_file= 'data/4k69_protein.pdb'
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure("protein", pdb_file)

# Extract the one letter aminoacid PDB sequence from the PDB file
pdb_sequence = SeqIO.read(pdb_file, "pdb-atom").seq
#Remove X characters
pdb_sequence = str(pdb_sequence).replace('X', '')

# Read in the protein sequence
#seq_file = sys.argv[2]
seq_file = 'L1000_chymase.fasta'
with open(seq_file) as f:
    protein_sequence = SeqIO.read(f, "fasta").seq

#Align the pdb_sequence and protein_sequence sequences using Bio.Align.PairwiseAligner
aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -4
aligner.extend_gap_score = -1
alignments = aligner.align(protein_sequence,pdb_sequence)
alignment = alignments[0]
print(alignment)
aligned_res_number=alignment.inverse_indices[1]

#Write code to loop over models, chains and residues in the PDB file
io = PDBIO()
index=0
for model in structure: #WARNING! Works for single model pdbs
    for chain in model:
        for residue in chain:
            # Create a new tuple with the modified residue number
            res_id = list(residue.id)
            res_id[1] = aligned_res_number[index] + 1
            res_id[2] = ' '
            residue.id = tuple(res_id)
            #print(residue.get_full_id())
            index += 1
    io.set_structure(model)
    io.save(os.path.join("renumber",pdb_file))
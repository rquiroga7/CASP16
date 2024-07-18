#PRINT HEADER
"PFRMAT     LG
TARGET     L1000
AUTHOR     3419-7309-1619
METHOD     Beta version of 2Vinardo, uses empirical docking scoring function F36 to predict pose and ligand binging affinity
"
#Gets protein pdb file coordinates from: C:\Users\Usuario\Documents\GitHub\CASP16\L1000_pdbs\PDBBIND\2hsv_protein.pdb 

#Gets ligand ID and name from first two columns of C:\Users\Usuario\Documents\GitHub\CASP16\L1000_ligands\L1001.tsv , for example "LIGAND 0	201"

# Given energy value in kcal/mol
delta_G <- -11
# Calculate Kd in M
Kd_nM <- exp(delta_G * 1000 / (1.987 * 298))* 1e9
# Print the result
Kd_nM


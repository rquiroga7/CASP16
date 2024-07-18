#for each of L1001 to L1016
for ligand in paste0("L",{1001..1016})

#PRINT HEADER
"PFRMAT     LG
TARGET     L1000
AUTHOR     3419-7309-1619
METHOD     Beta version of 2Vinardo, uses empirical docking scoring function F36 to predict pose and ligand binging affinity
"

#for MODEL x in 1 to 5 print
"MODEL x
PARENT 2hvx"

#Gets protein pdb file coordinates from: C:\Users\Usuario\Documents\GitHub\CASP16\L1000_results\Prot\2hvx_protein.pdb prints them here
#Gets ligand ID and name from first two columns of C:\Users\Usuario\Documents\GitHub\CASP16\L1000_ligands\L1001.tsv , for example "LIGAND 0	201" and prints them here

#Get delta_G for each model from C:\Users\Usuario\Documents\GitHub\CASP16\L1000_results\DOCK.F36.P2hvx\L1001_ligand.pdbt using grep to get lines with "REMARK 921   NORMALIZED FREE ENERGY" and printing the 7th column in each of those lines. First match is model 1, second is model 2 and so forth
# For each model, calculate Kd_nM using the delta_G from above
delta_G <- -11
# Calculate Kd in M
Kd_nM <- exp(delta_G * 1000 / (1.987 * 298))* 1e9
# Print the result as "AFFNTY x aa" where x is the Kd in nM, expressed in non scientific notation     



#print all to file L1001LG363_{MODEL}.txt


#Create tar
tar -czf L1000LG363.tgz ./L3000_models
#HEADER
"PFRMAT     LG
TARGET     L1000
AUTHOR     3419-7309-1619
METHOD     Beta version of 2Vinardo empirical docking scoring function used to predict pose and ligand bdinging affinity
"


# Given energy value in kcal/mol
delta_G <- -11
# Calculate Kd in M
Kd_nM <- exp(delta_G * 1000 / (1.987 * 298))* 1e9
# Print the result
Kd_nM
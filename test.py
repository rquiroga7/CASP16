import os
import pandas as pd
import numpy as np

path='C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/'
# Define the gas constant and temperature
R = float(1.987)
T = float(298)

header = "\n".join([
    "PFRMAT     LG",
    "TARGET     L1000",
    "AUTHOR     3419-7309-1619",
    "METHOD     Beta version of 2Vinardo, uses empirical docking scoring function F36 to predict pose and ligand binding affinity"
])





#Define list of ligands to process (L1001 to L1016)
lig_list = ['L1001', 'L1002', 'L1003', 'L1004', 'L1005', 'L1006', 'L1007', 'L1008', 'L1009', 'L1010', 'L1011', 'L1012', 'L1013', 'L1014', 'L1015', 'L1016']
#lig_list = ['L1001']
#Read the three files for each ligand, separate into the 5 models
for lig in lig_list:
    #Define the path for the ligand
    lig_path1 = path + 'DOCK.F36.P2hvx/' + lig + '_ligand.sdf'
    lig_path2 = path + 'DOCK.F36.P1t31/' + lig + '_ligand.sdf'
    lig_path3 = path + 'DOCK.F36.P3n7o/' + lig + '_ligand.sdf'
    #Open the file
    lig_file = open(lig_path1, 'r')
    #Read the file
    lig_data = lig_file.read()
    #Close the file
    lig_file = open(lig_path2, 'r')
    lig_data += lig_file.read()
    lig_file = open(lig_path3, 'r')
    lig_data += lig_file.read()
    lig_file.close()
    #Split the file into the 5 models
    models = lig_data.split('$$$$')
    #For models 0,5 and 10, grep the value of the NORMALIZED FREE ENERGY
    all_energies = []

    # Extract energy values from all models
    for model in models:
        lines = model.split('\n')
        for line in lines:
            if 'NORMALIZED FREE ENERGY' in line:
                all_energies.append(float(line.split()[5]))
                break
    # Compare energies with index 0, 5, and 10
    selected_indices = [0, 5, 10]
    selected_energies = [all_energies[i] for i in selected_indices if i < len(all_energies)]
    
    lowest_energy = min(selected_energies)
    lowest_energy_index = selected_energies.index(lowest_energy)
    #if lowest_energy_index == 0: then prot = '2hvx', if 1 then prot = '1t31', if 2 then prot = '3n7o'
    if lowest_energy_index == 0:
        prot = '2hvx'
        #save models 0 to 4 in models_low
        models_low = models[0:5]
        energies_low = all_energies[0:5]
    elif lowest_energy_index == 1:
        prot = '1t31'
        models_low = models[5:10]
        energies_low = all_energies[5:10]
    elif lowest_energy_index == 2:
        prot = '3n7o'
        models_low = models[10:15]
        energies_low = all_energies[10:15]
    # Define the paths
    protein_pdb_path = os.path.join("C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/Prot", f"{prot}_protein.pdb")
    ligand_info_path = os.path.join("C:/Users/Usuario/Documents/GitHub/CASP16/L1000_ligands", f"{lig}.tsv")
    ligand_pdbt_path = os.path.join(f"C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/DOCK.F36.P{prot}", f"{lig}_ligand.pdbt")

# Get protein pdb file coordinates
    with open(protein_pdb_path, 'r') as file:
        protein_pdb = file.readlines()

    # Get ligand ID and name
    ligand_info = pd.read_csv(ligand_info_path, sep='\t', header=None)
    ligand_id = ligand_info.iloc[1, 0]
    ligand_name = ligand_info.iloc[1, 1]

    # Get delta_G for each model in nM Kd equivalent
    Kd_nM = np.round(np.exp(np.array(energies_low) * 1000 / (R * T)) * 1e9,3)
    
    #for models 1 to 5
    for index in range(5):
        model=index+1
        # Create the model information string
        model_info1 = "\n".join([
            f"MODEL {model}",
            f"PARENT {prot}",
        ])
        model_info2 = "\n".join([
            f"LIGAND {ligand_id} {ligand_name}",
            f"{ligand_name}"
        ])
        end = "\n".join([
            "M  END"
            f"AFFNTY {Kd_nM[index]} aa",
            "END"
        ])

        # Create the output file name
        output_file = os.path.join("C:/Users/Usuario/Documents/GitHub/CASP16/L1000_models/", f"{lig}LG363_{model}")
        # Write the output file
        with open(output_file, 'w') as file:
            file.write("\n".join([header, model_info1,"".join(protein_pdb), model_info2, models_low[index], end]))
            file.close()

#Process C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/DOCK.F36.P2hvx/L1001_ligand.sdf, C:/Users/Usuario/Documents/GitHub/CASP16/L1000_results/DOCK.F36.P/L1001_ligand.sdf 
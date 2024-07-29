import os
import pandas as pd
import numpy as np

path='.'
# Define the gas constant and temperature
R = float(1.987)
T = float(298)

results_folder='./L1000_results2'
#protein list
prot_list = ['2hvx', '3s0n', '1t31', '4k2y', '4k5z', '4k60', '4k69', '4kp0']


#Define list of ligands to process (L1001 to L1016)
lig_list = ['L1001', 'L1002', 'L1003', 'L1004', 'L1005', 'L1006', 'L1007', 'L1008', 'L1009', 'L1010', 'L1011', 'L1012', 'L1013', 'L1014', 'L1015', 'L1016','L1017']
#lig_list = ['L1001']
#Read the three files for each ligand, separate into the 5 models
for lig in lig_list:
    header = "\n".join([
    "PFRMAT LG",
    f"TARGET {lig}",
    "AUTHOR 3419-7309-1619",
    "METHOD Beta version of 2Vinardo, uses empirical docking scoring function F40 to predict ligand pose and binding affinity"
    ])
    #Define the path for the ligand
    lig_path = {}
    lig_data = ""
    for prot in prot_list:
    #create a dictionary to store the ligand path for each protein use results_folder
        lig_path[prot] = os.path.join(results_folder, f"DOCK.F40.P{prot}", f"{lig}_ligand.sdf")
        #lig_path[prot] = path + '/L1000_results/DOCK.F40.P' + prot + '/' + lig + '_ligand.sdf'

    #lig_path1 = path + '/L1000_results/DOCK.F40.P' + f{prot} + '/' + lig + '_ligand.sdf'
    #lig_path2 = path + '/L1000_results/DOCK.F40.P' + f'1t31' + '/' + lig + '_ligand.sdf'
    #lig_path3 = path + '/L1000_results/DOCK.F40.P' + f'3n7o' + '/' + lig + '_ligand.sdf'
    #Open the file
        lig_file = open(lig_path[prot], 'r')
    #Read the file
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
    #for each model remove all lines after "M END"
    for i in range(len(models)):
        models[i] = models[i].split("M  END")[0] + "M  END"
        lines = models[i].split("\n")
        if len(lines) > 2:
            models[i] = "\n".join(lines[2:])
        else:
            models[i] = ""
    #remove empty lines in models
    models = [model for model in models if model]
    # Compare energies with index 0, 5, and 10 and so forth, use length of prot_list * 5
    selected_indices = [i for i in range(len(prot_list)*5) if i % 5 == 0]
    #selected_indices = [0, 5, 10]
    selected_energies = [all_energies[i] for i in selected_indices if i < len(all_energies)]
    
    lowest_energy = min(selected_energies)
    lowest_energy_index = selected_energies.index(lowest_energy)
    #Take lowest energy index and assign the protein
    prot = prot_list[lowest_energy_index]
    #Take the models and energies for the lowest energy
    models_low = models[lowest_energy_index*5:lowest_energy_index*5+5]
    energies_low = all_energies[lowest_energy_index*5:lowest_energy_index*5+5]
    #if lowest_energy_index == 0:
    #    prot = '2hvx'
    #    #save models 0 to 4 in models_low
    #    models_low = models[0:5]
    #    energies_low = all_energies[0:5]
    #elif lowest_energy_index == 1:
    #    prot = '1t31'
    #    models_low = models[5:10]
    #    energies_low = all_energies[5:10]
    #elif lowest_energy_index == 2:
    #    prot = '3n7o'
    #    models_low = models[10:15]
    #    energies_low = all_energies[10:15]
    # Define the paths
    #protein_pdb_path = os.path.join("./L1000_results2/data/", f"{prot}_protein.pdb")
    #Now using results_folder
    protein_pdb_path = os.path.join(results_folder, f"data/{prot}_protein.pdb")
    ligand_info_path = os.path.join("./L1000_ligands", f"{lig}.tsv")
    #ligand_pdbt_path = os.path.join(f"./L1000_results2/DOCK.F40.P{prot}", f"{lig}_ligand.pdbt")
    ligand_pdbt_path = os.path.join(results_folder, f"DOCK.F40.P{prot}/{lig}_ligand.pdbt")

# Get protein pdb file coordinates
    with open(protein_pdb_path, 'r') as file:
        protein_pdb = file.readlines()

    # Get ligand ID and name
    ligand_info = pd.read_csv(ligand_info_path, sep='\t', header=None)
    ligand_id = ligand_info.iloc[1, 0]
    ligand_name = ligand_info.iloc[1, 1]

    # Get delta_G for each model in nM Kd equivalent
    Kd_nM = np.round((np.exp(np.array(energies_low) * 1000 / (R * T)) * 1e9),3)
    
    #for models 1 to 5
    for index in range(5):
        model=index+1
        # Create the model information string
        model_info1 = "\n".join([
            "MODEL 1",
            f"PARENT {prot}",
        ])
        model_info2 = "\n".join([
            f"LIGAND {ligand_id} {ligand_name}",
            f"{ligand_name}"
        ])
        end = "\n".join([
            f"AFFNTY {Kd_nM[index]} aa",
            "END"
        ])

        # Create the output file name
        output_file = os.path.join("./L1000_models_corrected/", f"{lig}LG363_{model}")
        # Write the output file
        with open(output_file, 'w') as file:
            protpdb="".join(protein_pdb)
            #remove last item in prot
            protpdb = protpdb[:-1]
            file.write("\n".join([header, model_info1,protpdb, model_info2, models_low[index], end]))
            file.close()


import numpy as np

R = float(1.987)   # Universal gas constant in J/(mol*K)
T = 298  # Temperature in Kelvin

with open('L1000_stage2/ana.F14.dG.dat', 'r') as file:
    f_14 = file.readlines()

# Process each line
f_14_processed = []
for line in f_14:
    columns = line.split()
    first_col = columns[0]
    second_col = float(columns[1])
    
    # Perform the calculation
    calculated_value = np.round((np.exp(second_col * 1000 / (R * T)) * 1e9), 3)
    
    # Replace the second column with the calculated value
    new_line = f"{first_col} {calculated_value} aa\n"
    f_14_processed.append(new_line)

# Save the processed lines back to the file
with open('LG363.affinities', 'w') as file:
    file.writelines(f_14_processed)
#!/bin/bash

i_range=(L1001 L1002 L1003 L1004 L1005 L1006 L1007 L1008 L1009 L1010 L1011 L1012 L1013 L1014 L1015 L1016 L1017)
k_range=(F10 F36 F40)
j_range=(1t31 2hvx 3s0n 4k2y 4k5z 4k60 4k69 4kp0)

# Create a file to store the results
echo "" > energies.txt

# Collect energy values
for i in "${i_range[@]}"; do
    for k in "${k_range[@]}"; do
        for j in "${j_range[@]}"; do
            ener=$(grep -m 1 "NORMALIZED FREE ENERGY" DOCK.${k}.P${j}/${i}_ligand.pdbt | awk '{print $7}')
            
            # Check if the energy value is a valid number
            if [[ -n "$ener" && "$ener" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
                # Output each $i $j $k combination with its energy value
                echo "$i $j $k $ener" >> energies.txt
            fi
        done
    done
done

# Process the file to find the most negative energy for each combination of $i and $k
awk '
{
    key = $1" "$3  # Combine $i and $k as key
    if (!(key in min_energies) || $4 < min_energies[key]) {
        min_energies[key] = $4
        best_j[key] = $2  # Store the corresponding $j value for the minimum energy
    }
}
END {
    for (key in min_energies) {
        # Print the combination of $i, $k, best $j, and the most negative energy
        print key, best_j[key], min_energies[key]
    }
}
' energies.txt

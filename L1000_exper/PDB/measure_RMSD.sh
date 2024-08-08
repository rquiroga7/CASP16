#!/bin/bash

#!/bin/bash
rm rmsd.txt
# Write a bash script using obrms to calculate the RMSD between L1001_ligand.pdb and L1001_ligand_*.pdbt, then do the same for L1002:L1017
for i in L1001 L1002 L1003 L1004 L1005 L1006 L1007 L1008 L1009 L1010 L1011 L1012 L1013 L1014 L1015 L1016 L1017; do
    for j in 1 2 3 4 5; do
        echo -n "$i $j " >> rmsd.txt
        obrms ${i}_ligand.pdb ${i}_ligand_${j}.pdbt 2> /dev/null >> rmsd.txt
    done
done

rm results.txt
unset ID_array
unset min_values
unset min_lines
declare -A min_values
declare -A min_lines
declare -a ID_array

while read line ; do
    ID=$(echo $line | awk '{print $1}')
    num=$(echo $line | awk '{print $2}')
    val=$(echo $line | awk '{print $5}')
    #echo $ID $num $val
    #If the RMSD is less than 2, and ID is not in the ID_array
    if [[ $(echo "$val < 2" | bc -l) -eq 1 ]]; then
        #echo "RMSD less than 2";
        if [[ ! " ${ID_array[@]} " =~ " ${ID} " ]]; then
            #echo "ID not in ID_array";
            echo $line >> results.txt
            ID_array+=($ID)
        fi
    else
        if [[ ! ${min_values[$ID]} || $(echo "$val < ${min_values[$ID]}" | bc -l) -eq 1 ]]; then
            min_values[$ID]=$val
            min_lines[$ID]=$line
        fi
    fi
done < rmsd.txt

for ID in "${!min_values[@]}"; do
    if [[ ! " ${ID_array[@]} " =~ " ${ID} " ]]; then
        echo ${min_lines[$ID]} >> results.txt
    fi
done


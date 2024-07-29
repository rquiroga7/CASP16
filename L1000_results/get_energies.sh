for i in L1001 L1002 L1003 L1004 L1005 L1006 L1007 L1008 L1009 L1010 L1011 L1012 L1013 L1014 L1015 L1016 L1017; do
    echo "" >> energies.txt
    echo -n $i >> energies.txt
    for k in F10 F22 F36; do
        for j in 1t31 2hvx 3n7o; do
            ener=$(grep -m 1 "NORMALIZED FREE ENERGY" DOCK.${k}.P${j}/${i}_ligand.pdbt | awk '{print $7}')
            echo -n " $ener" >> energies.txt
        done
    done
done
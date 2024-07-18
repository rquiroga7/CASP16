for i in L1000;do
	awk '{if(NR > 1) print $0}' all_${i}_ligands.tsv > tmp_all.smiles
	while IFS= read -r line; do
	    # Extract the first and fourth columns
	    h=$(echo "$line" | awk '{print $1}')
	    j=$(echo "$line" | awk '{print $4}')
	    python CASP15_smi2sdf.py -s ${j} -f ./${i}_ligands/$h.sdf
	done < tmp_all.smiles
done


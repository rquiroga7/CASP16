for i in L1000 L2000 L3000 L4000; do
	awk '{print $4}' all_${i}_ligands.tsv > tmp.smi
	tail -n +2 tmp.smi > tmp2.smi
	mv tmp2.smi tmp.smi
	python visualize_ligands.py
	mv molecules.png ${i}.png
done

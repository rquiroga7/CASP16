#/bin/bash

for i in L1000_ligands L2000_ligands L3000_ligands L4000_ligands; do
	cd $i
	rm ../all_${i}.tmp
	for j in L*.tsv; do
		awk '{if ($1 == "ID"){print ("TARGET ",$0)} else {print substr(FILENAME,1,5), $0}}' ${j} >> ../all_${i}.tmp
	done
	cat ../all_${i}.tmp | awk '!x[$0]++' > ../all_${i}.tsv
	rm ../all_${i}.tmp
	cd ..
done

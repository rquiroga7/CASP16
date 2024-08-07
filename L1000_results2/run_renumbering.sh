for i in 1t31  3s0n 4k5z  4k69 2hvx  4k2y  4k60  4kp0; do
	cp data/${i}_protein.pdb .
	python renumber_res_pdb_alignment.py ${i}_protein.pdb
	rm ${i}_protein.pdb
done


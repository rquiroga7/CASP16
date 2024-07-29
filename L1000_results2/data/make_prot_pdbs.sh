for i in 2hvx 3s0n 1t31 4k2y 4k5z 4k60 4k69 4kp0; do
	obabel -ipdbt ${i}_protein.pdbt -opdb -O test.pdb ; egrep "ATOM" test.pdb > ${i}_protein.pdb
	echo "TER" >> ${i}_protein.pdb
done

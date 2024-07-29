#!/bin/bash


while read -r code ; do

     echo $code 
     obabel  -ipdb  ${code}_ligand.pdb   -opdbt -O ${code}_ligand.pdbt          >> OUT.$code.lig.txt
     obabel  -ipdb  ${code}_protein.pdb  -opdbt -O ${code}_protein.pdbt -xr -xp >> OUT.$code.pro.txt

done < lista.complex.txt


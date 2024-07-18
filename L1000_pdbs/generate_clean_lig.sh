grep -v "#" INDEX_general_PL_data.2020 | awk '{print $1,$8}' | sed 's/[()]//g' > INDEX_clean_lig.txt
grep -Fwf acc.list INDEX_clean_lig.txt > lig_names.all
awk '{print $2}' lig_names.all > ligs.txt

#for file in *.pdb; do
#	filename=$(basename "$file" .pdb)
#	grep -Fwf ligs.txt $file > ${filename}_lig.pdb
#done

direc1='/media/rquiroga/dae292c1-10dd-4169-83ec-c1360228407b/home/rodrigo/Datasets/PDBBIND_2020/PDBbind_v2020_refined/refined-set'
direc2='/media/rquiroga/dae292c1-10dd-4169-83ec-c1360228407b/home/rodrigo/Datasets/PDBBIND_2020/PDBbind_v2020_other_PL'


for file in *.pdb; do
    filename=$(basename "$file" .pdb)
    
    if [ -d "${direc1}/${filename}" ]; then
        cp "${direc1}/${filename}/${filename}_ligand.sdf" ./PDBBIND/
        cp "${direc1}/${filename}/${filename}_protein.pdb" ./PDBBIND/
    elif [ -d "${direc2}/${filename}" ]; then
        cp "${direc2}/${filename}/${filename}_ligand.sdf" ./PDBBIND/
        cp "${direc2}/${filename}/${filename}_protein.pdb" ./PDBBIND/
    fi
done

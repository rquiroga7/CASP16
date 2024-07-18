for i in DOCK.* ; do
    cd $i
    rm *.sdf
    rm *_protein.*
    for j in *.pdbt ; do
        name=$(basename $j .pdbt)
        obabel -ipdbt $name.pdbt -osdf -O $name.sdf
    done
    cd ..
done
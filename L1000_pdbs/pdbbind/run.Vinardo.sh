#!/bin/bash

function set_params {

    VINARDO_EXE=2vinardo-at21-HB_1_-1_ct6_G2_TXT_LR

    RMSDoSUCOS=RMSD
    cutoff_rmsd_ana=2.0

    DIR=./

    RUN=36

     DIR_FUNC=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/CASP16/Train/Rcte-Orde-0/HB_1.0_1.0/Ref-Xray/RepuCte_RMSD-2.0-2.0_VS-A-1.0_SR-LS-3x_noforza_ct6_G2/Err_rmsd_TxT_LR/Ronda.$RUN/Fit
     #DIR_FUNC=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/All/Rcte-Orde-0/HB_1.0_1.0/Ref-Xray/RepuCte_RMSD-2.0-2.0_VS-A-1.0_SR-LS-3x_noforza_ct6_G2/Err_rmsd/Ronda.$RUN/Fit

    DIR_PDBT_REF=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/CASP16/Exercise/pdbbind/data
}


function make_config_1 {

   BOX=${3}

   minX=`grep -E "ATOM|HETATM"  ${2}/${1}_ligand.pdbt | cut -c 31-55 | awk '{print $1}' | sort -nr | tail -n 1`
   maxX=`grep -E "ATOM|HETATM"  ${2}/${1}_ligand.pdbt | cut -c 31-55 | awk '{print $1}' | sort -n  | tail -n 1`

   minY=`grep -E "ATOM|HETATM"  ${2}/${1}_ligand.pdbt | cut -c 31-55 | awk '{print $2}' | sort -nr | tail -n 1`
   maxY=`grep -E "ATOM|HETATM"  ${2}/${1}_ligand.pdbt | cut -c 31-55 | awk '{print $2}' | sort -n  | tail -n 1`

   minZ=`grep -E "ATOM|HETATM"  ${2}/${1}_ligand.pdbt | cut -c 31-55 | awk '{print $3}' | sort -nr | tail -n 1`
   maxZ=`grep -E "ATOM|HETATM"  ${2}/${1}_ligand.pdbt | cut -c 31-55 | awk '{print $3}' | sort -n  | tail -n 1`

   MinX=`echo "$minX - $BOX " | bc `
   MaxX=`echo "$maxX + $BOX " | bc `

   MinY=`echo "$minY - $BOX " | bc `
   MaxY=`echo "$maxY + $BOX " | bc `

   MinZ=`echo "$minZ - $BOX " | bc `
   MaxZ=`echo "$maxZ + $BOX " | bc `

   cX=`echo "scale=2 ; ($MinX + $MaxX)/2.0" | bc `
   cY=`echo "scale=2 ; ($MinY + $MaxY)/2.0" | bc `
   cZ=`echo "scale=2 ; ($MinZ + $MaxZ)/2.0" | bc `

   lX=`echo " $MinX - $MaxX " | bc  | tr -d - `
   lY=`echo " $MinY - $MaxY " | bc  | tr -d - `
   lZ=`echo " $MinZ - $MaxZ " | bc  | tr -d - `

   echo "
   receptor = ${2}/${1}_protein.pdbt
   ligand   = ${2}/${1}_ligand.pdbt

   center_x=$cX
   center_y=$cY
   center_z=$cZ

   size_x=$lX
   size_y=$lY
   size_z=$lZ
   " > config
}

function make_config_2 {

  read   wRP sRP sSR sLR sHB gHB scale_intra scale_intra_14  scale_tors  offset  wRPHB < ${1}
  local table=${2}
  local tableTxT=${3}

  echo "

  out = ./DOCK

  granularity = 0.20
  conformations = 5

  wRP =  $wRP
  wSR = -1.0
  wLR = -1.0
  wHB = -1.0

  wRPHB = $wRPHB

  sRP =  $sRP
  sSR =  $sSR
  sLR =  $sLR
  sHB =  $sHB
  gHB =  $gHB

  scale_intra    = $scale_intra
  scale_intra_14 = $scale_intra_14
  scale_tors     = $scale_tors

  offset   = $offset

  table    = $table
  tableTxT = $tableTxT

  " >> config
}

function make_config {
	local code=$1
	local dir=$2
	local box=$3
	local scfun=$4
	local table=$5
    local tableTxT=$6
    make_config_1 $code $dir $box
	make_config_2 $scfun $table $tableTxT
}


function score_min_dock  {

    local nmol=$(wc -l $DIR/lista.txt | awk '{print $1}')
    local cont=1

    echo -e "\nScoring, Minimization, Docking"

    while read code ; do

	echo -ne "Complejo $code ( $(( cont++ )) / $nmol) \r"

        #### Score only
        #echo "   Score only"
        make_config $code $DIR_PDBT_REF 2.0  $DIR_FUNC/scfun.dat $DIR_FUNC/param.dat  $DIR_FUNC/param.TxT.dat
	$VINARDO_EXE  --threads 24 --score_only   --config config --useINTRA --out SCORE  > $code.score.out #2> /dev/null

	#### Minimzacion
	#echo "   Minimización"
        make_config $code $DIR_PDBT_REF 4.0  $DIR_FUNC/scfun.dat $DIR_FUNC/param.dat  $DIR_FUNC/param.TxT.dat
	$VINARDO_EXE  --threads 24   --minimizeBFGS  --config config  --useINTRA --out MINIM  > $code.minim.out #2> /dev/null

	#### Docking
	#echo "   Docking"
        make_config $code $DIR_PDBT_REF 6.0  $DIR_FUNC/scfun.dat $DIR_FUNC/param.dat  $DIR_FUNC/param.TxT.dat
	$VINARDO_EXE  --threads 24   --task 2000  --config config  --useINTRA --out DOCK  > $code.dock.out #2> /dev/null

	#cp config config.$code.txt
	if [ -e config ] ; then rm config ; fi


    done < $DIR/lista.txt


    echo -e "\nFin del Docking"

}


function aux_rmsd_sucos {

    local score_o_minim=${1} #Esto vale SCORE o MINIM
    local code=${2}
    local pose=${3}

	#rmsd=`DockRMSD SCORE/${code}.ligand.score.mol2  DOCK/${code}.ligand.dock$i.mol2 -s`
	rmsd=$(obrms $score_o_minim/$code.ligand.$score_o_minim.mol2  DOCK/$code.ligand.DOCK$pose.mol2 | awk '{print $3}')
	#sucos=`python3 /usr/local/SuCOS/calc_SuCOS_mol2_Ro.py --lig1 SCORE/${code}.ligand.score.mol2 --lig2  DOCK/${code}.ligand.dock$i.mol2`
	sucos=1.0
	inter_E=$(grep "NORMALIZED FREE ENERGY"  DOCK/${code}_ligand.pdbt  |  awk '{print $7}' | awk -v i=$pose '{if(NR==i) print}')
    printf "%-20s %4i %7.2f  %6.3f %6.3f\n" $code  $pose  $inter_E  $rmsd  $sucos > DvS.$score_o_minim.$code.$pose.tmp
}


function aux_mol2 {
    obabel -ipdbt SCORE/${1}_ligand.pdbt -omol2 -O SCORE/${1}.ligand.SCORE.mol2    2> /dev/null
	obabel -ipdbt MINIM/${1}_ligand.pdbt -omol2 -O MINIM/${1}.ligand.MINIM.mol2    2> /dev/null
    obabel -ipdbt  DOCK/${1}_ligand.pdbt  -omol2 -O DOCK/${1}.ligand.DOCK.mol2  -m 2> /dev/null
}


function calc_rmsd_sucos {

    echo -e "\nPasando todos los resultados a mol2"

    while read code ; do
	    aux_mol2 $code &
    done < $DIR/lista.txt
    wait

    echo "Calculando el RMSD y SuCOS de docking_vs_xray y docking_vs_minim"
    while read code ; do

	    Ns=$(ls DOCK/${code}.ligand.DOCK*.mol2 | wc -l)

        for pose in `seq 1 1 $Ns` ; do
	        aux_rmsd_sucos SCORE $code $pose &
	        aux_rmsd_sucos MINIM $code $pose &
        done

    done < $DIR/lista.txt
    wait

    cat DvS.SCORE.*.tmp > rmsd.docking_vs_xray.txt
    cat DvS.MINIM.*.tmp > rmsd.docking_vs_minim.txt

    local inf_xray=$( grep inf rmsd.docking_vs_xray.txt  | wc -l)
    local inf_minim=$(grep inf rmsd.docking_vs_minim.txt | wc -l)

    echo "RMSD infinito dock_vs_xray = $inf_xray    rmsd infinito dock_vs_minim = $inf_minim"

    rm DvS.*.tmp

}


function aux_procesa_dG {

    local code=$1

    dG_score=$(tail -n 1 $code.score.out | awk '{print $4}')
    dG_minim=$(tail -n 1 $code.minim.out | awk '{print $4}')
    dG_dock=$(tail -n 1 $code.dock.out  | awk '{print $4}')

    #RMSD de la minimizacion
    RMSminim=$(obrms  SCORE/$code.ligand.SCORE.mol2  MINIM/$code.ligand.MINIM.mol2 | awk '{print $3}')
    RMSdock=$(obrms   SCORE/$code.ligand.SCORE.mol2   DOCK/$code.ligand.DOCK1.mol2 | awk '{print $3}')
    #SUCOS=`python3 /usr/local/SuCOS/calc_SuCOS_mol2_Ro.py --lig1  SCORE/${code}.ligand.score.mol2  --lig2 MINIM/${code}.ligand.minim.mol2`
    SUCOS=1.0

    dG_ref=$(grep $code Eref.dat | awk '{print $2}')

    printf "%-10s  %7.2f  %7.2f %7.2f  %7.2f %7.2f %7.2f %7.2f \n" $code $dG_ref $dG_score $dG_minim $dG_dock $RMSminim $RMSdock  $SUCOS | grep -v "valid" > dG.$code.tmp
}



function procesa_dG {

    echo -e "\nProceso los dG "

    while read code ; do
	    aux_procesa_dG $code &
    done < $DIR/lista.txt
    wait

    cat dG.*.tmp > Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt
    
    awk '{print $2, $3, "#" $1}' Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt > dG_Ref.vs.dG_score.dat
    awk '{print $2, $4, "#" $1}' Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt > dG_Ref.vs.dG_minim.dat
    awk '{print $2, $5, "#" $1}' Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt > dG_Ref.vs.dG_dock.dat
    awk '{print $4, $5, "#" $1}' Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt > dG_minim.vs.dG_dock.dat
    awk '{print $6    , "#" $1}' Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt > Rmsd.minim.dat
    awk '{print $7    , "#" $1}' Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt > Rmsd.dock.top1.dat
    awk '{print $8    , "#" $1}' Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt > Sucos.minim.dat
    awk '{print $6, $8, "#" $1}' Tabla.code.dGref.dGscore.dGminim.dGdock.RMSDminim.RMSDdock.sucos.txt > Rmsd.vs.Sucos.minim.dat

    rm  dG.*.tmp *.score.out *.minim.out *.dock.out
}


function ana_dock {

    echo -e "\nAnalizando el docking. Uso las estructuras de X-ray como referencia "
    
    local rmsd_file=${1}
    grep -v inf $rmsd_file > TMP

    Npdb=`awk '{print $1}' TMP | sort -u | wc -l`

    TOP1=` awk -v c=$cutoff_rmsd_ana '{if ($4 <=c && $2<=1 && $2!=0  ) print }'  TMP | uniq -w 4 | wc -l | awk -v N=$Npdb '{printf("%6.1f\n", 100*($1/N))}'`
    TOP3=` awk -v c=$cutoff_rmsd_ana '{if ($4 <=c && $2<=3 && $2!=0  ) print }'  TMP | uniq -w 4 | wc -l | awk -v N=$Npdb '{printf("%6.1f\n", 100*($1/N))}'`
	TOP5=` awk -v c=$cutoff_rmsd_ana '{if ($4 <=c && $2<=5 && $2!=0  ) print }'  TMP | uniq -w 4 | wc -l | awk -v N=$Npdb '{printf("%6.1f\n", 100*($1/N))}'`
    S=`echo "$TOP1*5 + $TOP3*3 + $TOP5" | bc`
	printf " %-7s %-7s %-7s %-7s %4i \n" "$TOP1" "$TOP3" "$TOP5" "$S" "$Npdb" >  Dock.rmsd.out

    TOP1=` awk -v c=$cutoff_sucos_ana '{if ($5 >=c && $2<=1 && $2!=0  ) print }'  TMP | uniq -w 4 | wc -l | awk -v N=$Npdb '{printf("%6.1f\n", 100*($1/N))}'`
    TOP3=` awk -v c=$cutoff_sucos_ana '{if ($5 >=c && $2<=3 && $2!=0  ) print }'  TMP | uniq -w 4 | wc -l | awk -v N=$Npdb '{printf("%6.1f\n", 100*($1/N))}'`
	TOP5=` awk -v c=$cutoff_sucos_ana '{if ($5 >=c && $2<=5 && $2!=0  ) print }'  TMP | uniq -w 4 | wc -l | awk -v N=$Npdb '{printf("%6.1f\n", 100*($1/N))}'`
    S=`echo "$TOP1*5 + $TOP3*3 + $TOP5" | bc`
    printf " %-7s %-7s %-7s %-7s %4i \n" "$TOP1" "$TOP3" "$TOP5" "$S" "$Npdb" >  Dock.sucos.out

    rm TMP
}



set_params $1 $2 $3

score_min_dock

calc_rmsd_sucos

#procesa_dG

ana_dock  rmsd.docking_vs_xray.txt

mv Dock.rmsd.out Dock.rmsd.$RUN.out
mv rmsd.docking_vs_xray.txt  rmsd.docking_vs_xray.$RUN.txt

rm *.dock.out *.minim.out *.score.out Dock.sucos.out  rmsd.docking_vs_minim.txt
rm -r MINIM SCORE
mv DOCK  DOCK.$RUN



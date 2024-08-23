#!/bin/bash

function set_params {

    VINARDO_EXE=2vinardo-at21-HB_1_1_ct6_G2_rmsdOK
    #VINARDO_EXE=2vinardo-at21-HB_1_-1_ct6_G2_TXT_LR_rmsdOK
    DOCK_BOX=esta_fija_en_el_codigo

    RMSDoSUCOS=RMSD
    cutoff_rmsd_ana=2.0

    local run=$1

    DIR_LIG_REF=./data
    #DIR_FUNC=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/CASP16/Train/Rcte-Orde-0/HB_1.0_1.0/Ref-Xray/RepuCte_RMSD-2.0-2.0_VS-A-1.0_SR-LS-3x_noforza_ct6_G2/Err_rmsd_TxT_LR/Ronda.$RUN/Fit
     DIR_FUNC=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/All/Rcte-Orde-0/HB_1.0_1.0/Ref-Xray/RepuCte_RMSD-2.0-2.0_VS-A-1.0_SR-LS-3x_noforza_ct6_G2/Err_rmsd/Ronda.$run/Fit
    #DIR_FUNC=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/All/Rcte-Orde-0/HB_1.0_1.0/Ref-Xray/RepuCte_RMSD-2.0-2.0_VS-A-1.0_SR-LS-3x_noforza_ct6_G2/Err_rmsd_TxT_LR/Ronda.$run/Fit
    #DIR_FUNC=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/CASP16/Train.small2/Rcte-Orde-0/HB_1.0_1.0/Ref-Xray/RepuCte_RMSD-2.0-2.0_VS-A-1.0_SR-LS-3x_noforza_ct6_G2/Err_rmsd_TxT_LR/Ronda.$run/Fit
    DIR_DATA=./data

}


function make_config_1 {

   local LIG=${1}
   local PRO=${1}
   local DAT=${2}
   local BOX=${3}

   minX=`grep -E "ATOM|HETATM"  ${DAT}/${LIG}_ligand.pdbt | cut -c 31-55 | awk '{print $1}' | sort -nr | tail -n 1`
   maxX=`grep -E "ATOM|HETATM"  ${DAT}/${LIG}_ligand.pdbt | cut -c 31-55 | awk '{print $1}' | sort -n  | tail -n 1`

   minY=`grep -E "ATOM|HETATM"  ${DAT}/${LIG}_ligand.pdbt | cut -c 31-55 | awk '{print $2}' | sort -nr | tail -n 1`
   maxY=`grep -E "ATOM|HETATM"  ${DAT}/${LIG}_ligand.pdbt | cut -c 31-55 | awk '{print $2}' | sort -n  | tail -n 1`

   minZ=`grep -E "ATOM|HETATM"  ${DAT}/${LIG}_ligand.pdbt | cut -c 31-55 | awk '{print $3}' | sort -nr | tail -n 1`
   maxZ=`grep -E "ATOM|HETATM"  ${DAT}/${LIG}_ligand.pdbt | cut -c 31-55 | awk '{print $3}' | sort -n  | tail -n 1`

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
   receptor = ${DAT}/${PRO}_protein.pdbt
   ligand   = ${DAT}/${LIG}_ligand.pdbt

   center_x=-24.0
   center_y=20.0
   center_z=6.0

   size_x=22
   size_y=22
   size_z=20

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
    local data=$2
    local box=$3
    local scfun=$4
    local table=$5
    local tableTxT=$6
    make_config_1 $code $data $box
    make_config_2 $scfun $table $tableTxT
}


function dock  {

    local run=$1
    local nmol=$(wc -l lista.txt | awk '{print $1}')
    local cont=1

    echo -e "\nScoring, Minimization, Docking"


    while read code ; do

        echo -ne "Complejo $code ( $(( cont++ )) / $nmol) \r"

	make_config $code $DIR_DATA  $DOCK_BOX  $DIR_FUNC/scfun.dat $DIR_FUNC/param.dat  $DIR_FUNC/param.TxT.dat
        $VINARDO_EXE  --threads 24   --task 1500  --config config  --useINTRA --out DOCK  > $code.dock.out #2> /dev/null
	clina $run

	if [ -e config ] ; then rm config ; fi

    done < lista.txt

    echo -e "\nFin del Docking"

}


function calc_sucos {

    local run=$1

    while read code ; do
	    obabel -ipdbt  DOCK.F$run/${code}_ligand.pdbt  -omol2 -O DOCK.F$run/${code}.ligand.DOCK.mol2 -p 7 -m 2> /dev/null
	    obabel -ipdbt  DOCK.F$run/${code}_ligand.pdbt  -osdf  -O DOCK.F$run/${code}.ligand.DOCK.sdf  -p 7 -m 2> /dev/null
    done < lista.txt

    if [ -e DOCK.F$run/sucos.txt ] ; then rm DOCK.F$run/sucos.txt ; fi
    if [ -e DOCK.F$run/rmsd.txt  ] ; then rm DOCK.F$run/rmsd.txt  ; fi

    while read code ; do
	for i in 1 2 3 4 5 ; do
	    sucos=$( python3 /usr/local/SuCOS/calc_SuCOS_normalized.py --lig1 $DIR_LIG_REF/${code}_ligand.sdf  --lig2 DOCK.F$run/$code.ligand.DOCK$i.sdf | awk '{print $3}')
	    echo $code $run $i $sucos
	done >> DOCK.F$run/sucos.txt
    done < lista.txt

    while read code ; do
	for i in 1 2 3 4 5 ; do
	    rmsd=$( obrms  $DIR_LIG_REF/${code}_ligand.mol2  DOCK.F$run/$code.ligand.DOCK$i.mol2 | awk '{print $3}' )
	    echo $code $run $i $rmsd
	done >> DOCK.F$run/rmsd.txt
    done < lista.txt



}

function ana_dG_rmsd_sucos {

    local run=$1
    while read code ; do
        G=$(grep -A1 "MODEL        1" DOCK.F$run/${code}_ligand.pdbt | awk '{print $7}')
	echo $code $G
    done <lista.txt > ana.F$run.dG.dat

    awk '{if($3==1) print }' DOCK.F$run/rmsd.txt  > ana.F$run.rmsd.dat
    awk '{if($3==1) print }' DOCK.F$run/sucos.txt > ana.F$run.sucos.dat

}


function clina {

    local run=$1

    if [ ! -d DOCK.F$run ] ; then mkdir DOCK.F$run ; fi
    cp DOCK/*  DOCK.F$run
    rm *.dock.out
    rm -r DOCK/

}


for RUN in 10 ; do
	set_params $RUN
        dock       $RUN
	calc_sucos $RUN
	ana_dG_rmsd_sucos $RUN
done


#!/bin/bash

function set_params {

    VINARDO_EXE=2vinardo-at21-HB_1_-1_ct6_G2_TXT_LR
    DOCK_BOX=5.0

    RMSDoSUCOS=RMSD
    cutoff_rmsd_ana=2.0

    local pro=$1
    local run=$2

    DIR_LIG_REF=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/CASP16/Exercise/pdbbind/data
    DIR_FUNC=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/CASP16/Train/Rcte-Orde-0/HB_1.0_1.0/Ref-Xray/RepuCte_RMSD-2.0-2.0_VS-A-1.0_SR-LS-3x_noforza_ct6_G2/Err_rmsd_TxT_LR/Ronda.$RUN/Fit
    #DIR_FUNC=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/All/Rcte-Orde-0/HB_1.0_1.0/Ref-Xray/RepuCte_RMSD-2.0-2.0_VS-A-1.0_SR-LS-3x_noforza_ct6_G2/Err_rmsd/Ronda.$run/Fit
    DIR_DATA=/home/marcos/Calculos/2Vinardo-fit/Refined-2020/Long.tiposOK.rings/CASP16/Exercise/L1000/data

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
   receptor = ${2}/${4}_protein.pdbt
   ligand   = ${2}/${1}_ligand.pdbt

   center_x=0.0
   center_y=38.0
   center_z=22.0

   size_x=24
   size_y=24
   size_z=24
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
    local pro=$7
    make_config_1 $code $data $box $pro
    make_config_2 $scfun $table $tableTxT
}


function dock  {

    local pro=$1
    local nmol=$(wc -l lista.txt | awk '{print $1}')
    local cont=1

    echo -e "\nScoring, Minimization, Docking"

    while read code ; do

	echo -ne "Complejo $code ( $(( cont++ )) / $nmol) \r"

        make_config $code $DIR_DATA  $DOCK_BOX  $DIR_FUNC/scfun.dat $DIR_FUNC/param.dat  $DIR_FUNC/param.TxT.dat $pro
	$VINARDO_EXE  --threads 24   --task 2000  --config config  --useINTRA --out DOCK  > $code.dock.out #2> /dev/null

	if [ -e config ] ; then rm config ; fi

    done < lista.txt

    echo -e "\nFin del Docking"

}


function calc_sucos {

    local pro=$1
    local run=$2

    while read code ; do
	    obabel -ipdbt  DOCK.F$run.P$pro/${code}_ligand.pdbt  -omol2 -O DOCK.F$run.P$pro/${code}.ligand.DOCK.mol2  -m 2> /dev/null
    done < lista.txt

    if [ -e DOCK.F$run.P$pro/sucos.txt ] ; then rm DOCK.F$run.P$pro/sucos.txt ; fi

    while read code ; do
	while read ref ; do
	    sucos=$( python3 /usr/local/SuCOS/calc_SuCOS_mol2_Ro.py --lig1 $DIR_LIG_REF/${ref}_ligand.mol2  --lig2 DOCK.F$run.P$pro/$code.ligand.DOCK1.mol2 )
	    echo $pro $run $code $ref $sucos
        done < lista.lig.pdbbind.txt | sort -nk5 | tail -n -1 >> DOCK.F$run.P$pro/sucos.txt
    done < lista.txt

}



function aver_sucos {

    local pro=$1
    local run=$2

    if [ -e DOCK.F$run.P$pro/sucos.txt ] ; then
	S=$(awk '{ a+=$5 } END {print a/NR }' DOCK.F$run.P$pro/sucos.txt  )
	echo $pro $run $S
    fi

}


function clina {

    local pro=$1
    local run=$2

    rm *.dock.out
    mv DOCK  DOCK.F$RUN.P$pro

}


for RUN in 41 ; do
    for PRO in 1t31 2hvx 3n7o 3s0n 4k2y 4k5z 4k60 4k69 4kp0 5yjm 5yjp ; do
	set_params $PRO $RUN
        dock $PRO
	clina      $PRO $RUN
	calc_sucos $PRO $RUN
	aver_sucos $PRO $RUN
    done
done


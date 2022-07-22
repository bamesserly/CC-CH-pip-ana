#!/bin/bash

#Usage: runTransWarp.sh runEventLoopMC.root warped.root

declare -a VARIABLE=("adphi" "cosadtheta" "ehad" "enu" "pimuAngle" "pmu" "ptmu" "pzmu" "q2" "thetamu_deg" "thetapi_deg" "tpi" "wexp")
declare -a warps=("WARP1" "WARP2" "WARP3" "WARP4")
#declare -a VARIABLE=("tpi")
#declare -a warps=("WARP1")
Ximaxaxis=(1500 10000 1000 100 650 75 4000 100 100 50 550 275 800 10000 10000 200 100 2200 75 4000 100 50 50 850 150 150 1100 10000 1200 100 450 75 4500 100 75 50 200 50 550 1000 10000 1400 100 450 75 5000 100 250 50 500 100 5000)
#Ximaxaxis=(300)

#MIGRATION_FILE=$1
#TRUE_HIST=effnum_${VARIABLE}
#WARPED_FILE=$2
#RECO_HIST=selection_mc_${VARIABLE}

OUTFILE_NAME="/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/WarpingStudies/"
#OUTFILE_NAME=$(basename $2)
counter=0

for TAG in "${warps[@]}"; do
  for v in "${VARIABLE[@]}"; do
    MIGRATION_FILE="MCXSecInputs_20220407_NOMINAL.root"
    WARPED_FILE="MCXSecInputs_20220407_${TAG}.root"

    TransWarpExtraction --output_file ${OUTFILE_NAME}Warping_${TAG}_${v}.root --data effnum_${v} --data_file $WARPED_FILE --data_truth effnum_${v}_true --data_truth_file $WARPED_FILE --migration migration_${v} --migration_file $MIGRATION_FILE --reco effnum_${v} --reco_file $MIGRATION_FILE --truth effnum_${v}_true --truth_file $MIGRATION_FILE --num_uni 500 --max_chi2 ${Ximaxaxis[${counter}]} --step_chi2 0.5 --num_iter 0,1,2,3,4,5,10,20,30,50,100,200 --log_scale -C 0.5
#echo "Variable ${v} Warp ${TAG} Xi Y Axis ${Ximaxaxis[${counter}]}"
    cd ../MAT/macros/
    python PrintWarpingStudy.py -i ${OUTFILE_NAME}Warping_${TAG}_${v}.root -o ${OUTFILE_NAME}${v}_${TAG} -L
    cd -
let counter=counter+1
  done # vars
done #warps



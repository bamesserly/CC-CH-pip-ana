#!/bin/bash

#Usage: runTransWarp.sh runEventLoopMC.root warped.root

declare -a VARIABLE=("tpi" "thetapi_deg" "pmu" "thetamu_deg" "enu" "q2" "wexp" "ptmu" "pzmu" "ehad" "cosadtheta" "adphi" "pimuAngle")
declare -a warps=("WARP1" "WARP2" "WARP3" "WARP4")
#MIGRATION_FILE=$1
#TRUE_HIST=effnum_${VARIABLE}
#WARPED_FILE=$2
#RECO_HIST=selection_mc_${VARIABLE}

OUTFILE_NAME="/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/WarpingStudies/"
#OUTFILE_NAME=$(basename $2)

for TAG in "${warps[@]}"; do
  for v in "${VARIABLE[@]}"; do
    MIGRATION_FILE="MCXSecInputs_20220407_${TAG}.root"
    WARPED_FILE="MCXSecInputs_20220407_${TAG}.root"

    TransWarpExtraction --output_file ${OUTFILE_NAME}Warping_${TAG}_${v}.root --data effnum_${v} --data_file $WARPED_FILE --data_truth effnum_${v}_true --data_truth_file $WARPED_FILE --migration migration_${v} --migration_file $MIGRATION_FILE --reco effnum_${v} --reco_file $MIGRATION_FILE --truth effnum_${v}_true --truth_file $MIGRATION_FILE --num_uni 500 --max_chi2 25 --step_chi2 0.5 --num_iter 0,1,2,3,4,5,10,20,30,50,100,200 --log_scale -C 25

    cd ../MAT/macros/
    python PrintWarpingStudy.py -i ${OUTFILE_NAME}Warping_${TAG}_${v}.root -o ${OUTFILE_NAME}${v}_${TAG} -L
    cd -
  done # vars
done #warps



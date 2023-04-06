#!/bin/bash

#Usage: runTransWarp.sh runEventLoopMC.root warped.root

#declare -a VARIABLE=("adphi" "cosadtheta" "ehad" "enu" "pimuAngle" "pmu" "ptmu" "pzmu" "q2" "thetamu_deg" "thetapi_deg" "tpi" "wexp")
declare -a VARIABLE=("pzmu%s_vs_ptmu%s" "pmu%s_vs_thetamu_deg%s" "tpi%s_vs_pmu%s" "tpi%s_vs_thetapi_deg%s" "ptmu%s_vs_tpi%s")
#declare -a VARIABLE=("pzmu%s_vs_ptmu%s" "pmu%s_vs_thetamu_deg%s")
#declare -a warps=("NOMINAL")
#declare -a VARIABLE=("tpi")
declare -a warps=("WARP1")
#Ximaxaxis=(1500 10000 1000 100 650 75 4000 100 100 50 550 275 800 10000 10000 200 100 2200 75 4000 100 50 50 850 150 150 1100 10000 1200 100 450 75 4500 100 75 50 200 50 550 1000 10000 1400 100 450 75 5000 100 250 50 500 100 5000)
Ximaxaxis=(5000 5000 5000 5000)
true_tag="_true"
reco_tag=""
#MIGRATION_FILE=$1
#TRUE_HIST=effnum_${VARIABLE}
#WARPED_FILE=$2
#RECO_HIST=selection_mc_${VARIABLE}

OUTFILE_NAME="/minerva/data/users/granados/WarpingStudies/Warping/2DWarping/"
#OUTFILE_NAME=$(basename $2)
counter=0

for TAG in "${warps[@]}"; do
  for v in "${VARIABLE[@]}"; do
    MIGRATION_FILE="MCXSecInputs_20230328_NOMINAL.root"
    WARPED_FILE="MCXSecInputs_20230328_${TAG}.root"
    reco_var=$(printf "$v" $reco_tag $reco_tag)
    true_var=$(printf "$v" $true_tag $true_tag)
    TransWarpExtraction --output_file ${OUTFILE_NAME}Warping_2D${TAG}_${reco_var}.root --data effnum_${reco_var} --data_file $WARPED_FILE --data_truth effnum_${true_var} --data_truth_file $WARPED_FILE --migration Migration2d_${reco_var}_migration --migration_file $MIGRATION_FILE --reco effnum_${reco_var} --reco_file $MIGRATION_FILE --truth effnum_${true_var} --truth_file $MIGRATION_FILE --num_uni 500 --max_chi2 ${Ximaxaxis[${counter}]} --step_chi2 0.5 --num_iter 0,1,2,3,4,5,10,20,30,50,100 --log_scale -C 0.5 --num_dim 2
    echo "Variable ${reco_var} ${true_var} Warp ${TAG} Xi Y Axis ${Ximaxaxis[${counter}]}"
    cd ../MAT/macros/
    python PrintWarpingStudy.py -i ${OUTFILE_NAME}Warping_2D${TAG}_${reco_var}.root -o ${OUTFILE_NAME}${reco_var}_${TAG} -L
    cd -
    let counter=counter+1
  done # vars
done #warps

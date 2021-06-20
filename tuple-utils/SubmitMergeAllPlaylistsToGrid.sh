#!/bin/bash

INDATE="20210307"
OUTDATE="20210307"
RELEASE_DIR="v22r1p1"

indir_format_str="pions/%s/%s/%s/"
#outdir_format_str="pions/%s/merged/%s/%s"
outdir_format_str="pions/%s/merged/%s/%s"

#declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D" "ME1E" "ME1F" "ME1G" "ME1L" "ME1M" "ME1N" "ME1O" "ME1P")
declare -a playlists=("ME1B" "ME1C" "ME1D")
#declare -a playlists=("ME1A")
declare -a datamc=("mc")

for playlist in "${playlists[@]}"; do
  for dm in "${datamc[@]}"; do
    indir=$(printf "$indir_format_str" $INDATE $dm $playlist)
    outdir=$(printf "$outdir_format_str" $OUTDATE $dm $playlist)
    #if [ ! -d $outdir ]; then
    #  mkdir -p $outdir;
    #fi
    #echo ${indir}
    #echo ${outdir}
    python $PLOTUTILSROOT/scripts/submitMergeToGrid.py --tool CCNuPionInc --release ${RELEASE_DIR} --inputdir ${indir} --outputdir $outdir --${dm} --release_dir ${RELEASE_DIR} --tracker_only
  done # data-mc
done # playlists

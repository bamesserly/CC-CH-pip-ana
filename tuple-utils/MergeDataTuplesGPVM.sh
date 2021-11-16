#!/bin/bash
#
# This script calls the C script that merges data tuples. I (Ben) use this
# method to merge data tuples. It loops playlists and calls mergeDataFiles, my
# special script for merging all data from a playlist into a single file. The
# MAT/PU submitMergeToGrid.py script can't do this, and it's bad; event loops
# are WAY faster per-playlist than per run.
#
# It now works for MAD tuples.
#
# I used to use this script for MC, but now I use
# SubmitMergeAllPlaylistsToGrid.sh which calls submitMergeToGrid. Can't merge
# mc by-playlist.

IN_DATE="20211012"
OUT_DATE="20211115"
INV="v22r1p1"
MAD_VERSION="v49p1"
IN_TOPDIR="/pnfs/minerva/persistent/users/granados/"
OUT_TOPDIR="/pnfs/minerva/persistent/users/bmesserl/pions/"

data_indir_format_str="${IN_TOPDIR}/%s_MAD%s_${MAD_VERSION}/grid/minerva/ana/numibeam/${INV}/" # data-mc, plist
outdir_format_str="${OUT_TOPDIR}/${OUT_DATE}/merged/%s/%s" # data-mc, plist
log_format_str="%s_merge_%s_%s.txt"


### declare playlists and datamc arrays to loop over.
declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D" "ME1E" "ME1F" "ME1G" "ME1L" "ME1M" "ME1N" "ME1O" "ME1P")
declare -a datamc=("data")


# Loop playlists
for playlist in "${playlists[@]}"; do
  # Loop data-mc
  for dm in "${datamc[@]}"; do
    outdir=$(printf "$outdir_format_str" $dm $playlist)
    log=$(printf "$log_format_str" $IN_DATE $dm $playlist)
    # create outdir if it doesn't exist
    if [ ! -d $outdir ]; then
      mkdir -p $outdir;
    fi
    # data/mc
    if [ "$dm" = "data" ]; then
      indir=$(printf "$data_indir_format_str" $dm $playlist)
      echo `nohup root -b -q loadLibs.C+ "tuple-utils/mergeDataFiles.C+(\"${indir}\", \"${outdir}\")" > ${log}&`
    #elif [ "$dm" = "mc" ]; then
    #  indir=$(printf "$mc_indir_format_str" $IN_DATE $playlist)
    #  nohup `for ((i=first_run;i<=last_run;i++)); do root -b -q loadLibs.C+ "file-merging/mergeMCRun2.C+( \"${indir}\", \"${outdir}\", $i )"; done >> ${log}`&
    else
      echo "This script is for merging data only."
      echo "If you want to merge MC, you should use PU/MAT's submitMergeToGrid.py."
    fi
  done # end datamc loop
done # end playlist loop

echo "Done submitting"

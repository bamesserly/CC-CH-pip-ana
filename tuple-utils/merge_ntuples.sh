#!/bin/bash

DATE="20210307"
OUTDATE="20210307"
INV="v22r1p1"
PNFS_TOPDIR="/pnfs/minerva/persistent/users/bmesserl/pions/"
BLUEARC_TOPDIR="/minerva/data/users/bmesserl/pionproduction/"

data_indir_format_str="${PNFS_TOPDIR}/%s/data/%s/grid/minerva/ana/numibeam/${INV}/"
mc_indir_format_str="${PNFS_TOPDIR}/%s/mc/%s/grid/central_value/minerva/ana/${INV}/"
outdir_format_str="${BLUEARC_TOPDIR}/%s/merged/%s/%s"
log_format_str="%s_merge_%s_%s.txt"


### declare playlists and datamc arrays to loop over.
#declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D" "ME1E" "ME1F" "ME1G" "ME1L" "ME1M" "ME1N" "ME1O" "ME1P")
declare -a playlists=("ME1A" "ME1D" )
declare -a datamc=("data")


# Loop playlists
for playlist in "${playlists[@]}"; do
  # set playlist run numbers (for mc)
  first_run=0
  last_run=0
  if [ "$playlist" = "ME1A" ]; then 
    first_run=110000
    last_run=110020
  elif [ "$playlist" = "ME1B" ]; then  
    first_run=111000
    last_run=111000
    #last_run=111005
  elif [ "$playlist" = "ME1C" ]; then  
    first_run=111030
    last_run=111040
  elif [ "$playlist" = "ME1D" ]; then  
    first_run=111100
    last_run=111130
  elif [ "$playlist" = "ME1E" ]; then
    first_run=111325
    last_run=111350
  elif [ "$playlist" = "ME1F" ]; then
    first_run=111490
    last_run=111525
  elif [ "$playlist" = "ME1G" ]; then
    first_run=110150
    last_run=110180
  else
    echo "Invalid playlist"
  fi

  # Loop data-mc
  for dm in "${datamc[@]}"; do
    outdir=$(printf "$outdir_format_str" $OUTDATE $dm $playlist)
    log=$(printf "$log_format_str" $DATE $dm $playlist)
    # create outdir if it doesn't exist
    if [ ! -d $outdir ]; then
      mkdir -p $outdir;
    fi
    # data/mc
    if [ "$dm" = "data" ]; then
      indir=$(printf "$data_indir_format_str" $DATE $playlist)
      echo `nohup root -b -q loadLibs.C+ "file-merging/mergeDataFiles.C+(\"${indir}\", \"${outdir}\")" > ${log}&`
    elif [ "$dm" = "mc" ]; then
      indir=$(printf "$mc_indir_format_str" $DATE $playlist)
      nohup `for ((i=first_run;i<=last_run;i++)); do root -b -q loadLibs.C+ "file-merging/mergeMCRun2.C+( \"${indir}\", \"${outdir}\", $i )"; done >> ${log}`&
    else
      echo "Invalid data-mc"
    fi
  done # end datamc loop
done # end playlist loop

echo "Done submitting"

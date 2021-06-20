#!/bin/bash

DATE="20210307"

merged_format_str="/pnfs/minerva/persistent/users/bmesserl/pions/%s/merged/%s/%s/audit"
tag_format_str="%s_%s"

#declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D" "ME1E" "ME1F" "ME1G" "ME1L" "ME1M" "ME1N" "ME1O" "ME1P")
declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D")
declare -a datamc=("mc")

for playlist in "${playlists[@]}"; do
  for dm in "${datamc[@]}"; do
    merged=$(printf "$merged_format_str" $DATE $dm $playlist)
    tag=$(printf "$tag_format_str" $dm $playlist)
    python $PLOTUTILSROOT/scripts/SummarizeAudit.py $merged $tag
  done # data-mc
done # playlists

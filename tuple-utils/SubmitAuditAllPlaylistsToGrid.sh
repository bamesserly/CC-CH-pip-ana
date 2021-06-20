#!/bin/bash

UNMERGEDDATE="20210307"
MERGEDDATE="20210307"
RELEASE_DIR="v22r1p1"

unmerged_format_str="pions/%s/%s/%s/"
merged_format_str="pions/%s/merged/%s/%s/"

#declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D" "ME1E" "ME1F" "ME1G" "ME1L" "ME1M" "ME1N" "ME1O" "ME1P")
declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D")
declare -a datamc=("mc")

for playlist in "${playlists[@]}"; do
  for dm in "${datamc[@]}"; do
    unmerged=$(printf "$unmerged_format_str" $UNMERGEDDATE $dm $playlist)
    merged=$(printf "$merged_format_str" $MERGEDDATE $dm $playlist)
    python $PLOTUTILSROOT/scripts/submitAuditToGrid.py --tool CCNuPionInc --release ${RELEASE_DIR} --unmerged_dir ${unmerged} --merged_dir $merged --${dm} --release_dir ${RELEASE_DIR}
  done # data-mc
done # playlists

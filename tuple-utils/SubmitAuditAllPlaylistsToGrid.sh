#!/bin/bash

UNMERGEDDATE="20211012"
MERGEDDATE="20211115"
RELEASE_DIR="v22r1p1_UpMAD"

unmerged_format_str="/%s_MAD%s_v49p1/"
merged_format_str="MADtuplas/merged/%s/%s/%s/"

declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D" "ME1E" "ME1F" "ME1G" "ME1L" "ME1M" "ME1N" "ME1O" "ME1P")
#declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D")
declare -a datamc=("mc")

for playlist in "${playlists[@]}"; do
  for dm in "${datamc[@]}"; do
    unmerged=$(printf "$unmerged_format_str" $dm $playlist)
    merged=$(printf "$merged_format_str" $MERGEDDATE $dm $playlist)
    python $PLOTUTILSSCRIPTS/submitAuditToGrid.py --tool MasterAnaDev --release ${RELEASE_DIR} --unmerged_dir ${unmerged} --merged_dir $merged --${dm} --release_dir ${RELEASE_DIR} --memory 500 --disk 15 
  done # data-mc
done # playlists

#!/bin/bash

# From merged anatuples, make text playlists

DATE="20200713"
PNFS_TOPDIR="/pnfs/minerva/persistent/users/bmesserl/pions/"
BLUE_TOPDIR="/minerva/data/users/bmesserl/"
INV="v21r1p1"

indir_format_str="${PNFS_TOPDIR}/%s/merged/%s/%s"
outfile_format_str="${BLUE_TOPDIR}/MECCCHpip_ana_plists/%s/%s_%s_plist.txt"
outfile_xrootd_format_str="${BLUE_TOPDIR}/MECCCHpip_ana_plists/%s/%s_%s_xrootd_plist.txt"

#declare -a playlists=("ME1A" "ME1B" "ME1C" "ME1D" "ME1E" "ME1F" "ME1G" "ME1L" "ME1M" "ME1N" "ME1O" "ME1P")
declare -a playlists=("ME1A")
declare -a datamc=("data" "mc")

for playlist in "${playlists[@]}"; do
  for dm in "${datamc[@]}"; do
      indir=$(printf "$indir_format_str" $DATE $dm $playlist)
      outfile=$(printf "$outfile_format_str" $DATE $dm $playlist)
      outfile_xrootd=$(printf "$outfile_xrootd_format_str" $DATE $dm $playlist)
      if [ ! -f $outfile ]; then
        find ${indir} -name "*.root" -printf "%p\n" | sort >  ${outfile}
      fi
      if [ ! -f $outfile_xrootd ]; then
        sed 's|^/pnfs/minerva/\(.*\)\.root|root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr/minerva/\1\.root|g' $outfile > $outfile_xrootd
      fi
  done # datamc loop
done # playlist loop

# And for data, make an "ALL" playlist
indir=$(printf "$indir_format_str" $DATE "data")
outfile=$(printf "$outfile_format_str" $DATE "data" "ALL")
outfile_xrootd=$(printf "$outfile_xrootd_format_str" $DATE "data" "ALL")
if [ ! -f $outfile ]; then
  find ${indir} -name "*.root" -printf "%p\n" | sort >  ${outfile}
fi
if [ ! -f $outfile_xrootd ]; then
  sed 's|^/pnfs/minerva/\(.*\)\.root|root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr/minerva/\1\.root|g' $outfile > $outfile_xrootd
fi

#!/bin/bash

#===============================================================================
# Setup env variables, cmt config, source setup.sh packages
#===============================================================================
# The -n option [to export]causes the export property to be removed from each
# name.
#
# CONDOR_DIR_INPUT is the area where `-f` jobsub arguments (i.e. our tarball)
# are dropped.
#
echo "======== Set HOME = TOPDIR = CONDOR_DIR_INPUT ========"
export -n HOME 
export -n TOPDIR
export HOME=${CONDOR_DIR_INPUT}
export TOPDIR=${CONDOR_DIR_INPUT}

echo
echo "======== cd to HOME/TOPDIR/CONDOR_DIR_INPUT ========"
cd $HOME

echo
echo "======== pwd, ls -a ========"
pwd
ls -a

echo
echo "======== Untarring... ========"
tar xvzf ${TARFILE} -C ./ > /dev/null

echo
echo "======== ls -a ========"
ls -a

echo
echo "======== Setting Up Ana Packages ========"
cd ${CONDOR_DIR_INPUT}/Ana/CCPionInc/cmt/; cmt config; source setup.sh
cd ${CONDOR_DIR_INPUT}/Ana/PlotUtils/cmt; cmt config; source setup.sh
cd ${CONDOR_DIR_INPUT}/Ana/UnfoldUtils/cmt; cmt config; source setup.sh

echo
echo "======== cd to MAD_ME_1pi ========"
cd ${HOME}/Ana/CCPionInc/ana/MAD_ME_1pi
pwd

echo
echo "======== rm any pre-existing root files that got passed to the grid ========"
rm *.root

# had some trouble when *.so files already existed
echo
echo "======== clear out pre-existing .so, .d, .o files ========"
find . -type f -name '*.o' -delete
find . -type f -name '*.so' -delete
find . -type f -name '*.d' -delete

#===============================================================================
# Tell root (via the .rootrc) that every time root is open/run, it should first
# run rootlogon_grid.C, which is located in the cd.
# 
# rootlogon_grid contains the PU setup code.
#===============================================================================
echo
echo "======== make ./.rootrc ========"
echo Rint.Logon: ./rootlogon_grid.C > ./.rootrc

echo
echo "======== ls -a ========"
ls -a

echo
echo "======== cat .rootrc ========"
cat .rootrc

echo
echo "======== cat .rootrc ========"
cat rootlogon_grid.C

# To make sure we're using the right rootlogon.C 
echo
echo "======== gEnv->Print() ========"
echo 'gEnv->Print(); gSystem->Exit(0);' | root -b -l | grep Logon

#===============================================================================
# Ship it
# The "++" loadLibs might be redundant with removing the .so's, but just to be safe
#===============================================================================
echo
echo "======== MACRO ========"
echo $MACRO

echo 
echo "======== root.exe -b -q loadLibs.C++ MACRO ========"
root.exe -b -q -l loadLibs.C++ ${MACRO}

echo
echo "======== ls *.root ========"
ls *.root

echo
echo "Move *.root to CONDOR_DIR_OUT"
mv *.root $CONDOR_DIR_OUT

echo
echo "======== ls CONDOR_DIR_OUT ========"
ls $CONDOR_DIR_OUT

echo
echo "done"

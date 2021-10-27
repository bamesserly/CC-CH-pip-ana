//******************************************************************************
//*
//* Merge all data tuples in a given topdir into single file, regardless of run
//* number or playlist.
//* Original author - Laura Fields, taken from Phil Rodrigues.
//*
//* I use this script for my official analysis because it is much faster than
//* the PlotUtils script, which merges data tuples run-by-run. I use this
//* script to merge data tuples playlist-by-playlist (single file for each
//* playlist).
//*
//* Now streaming the files with xrootd.
//* 
//******************************************************************************

#ifndef __CINT__
#include "glob.h"
#endif

#include "merge_common.h"

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TList.h"

#include <iostream>
#include <string>

using namespace std;

//======================================================================
// Like mergeTrees, but keeps all the branches in CCQEAntiNuTool, and
// only merges one run.
void mergeDataFiles(const char* inDirBase, const char* outDir, const char* tag="MasterAnaDev") {
  //******************************************************************
  //* Set Location of output files
  //******************************************************************
  TString output=TString::Format("%s/merged_MAD.root", outDir);

  //******************************************************************
  //* Load Input Ntuples
  //******************************************************************
  TChain inChain("MasterAnaDev");

  TString inGlob(TString::Format("%s/*/*/*/*/MV_*%s*.root", inDirBase, tag));
  //TString inGlob(TString::Format("%s/*/*/*/*/MV_00016523_Subruns_000*%s*.root", inDirBase, tag));

  cout << "Filename glob is " << inGlob << endl;
  cout << "Output filename is " << output << endl;
  const std::string XROOTD_PREFIX = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/";
  glob_t g;
  glob(inGlob.Data(), 0, 0, &g);
  for(int i=0; i<(int)g.gl_pathc; ++i){
    std::string tup(g.gl_pathv[i]);
    tup = XROOTD_PREFIX + tup.substr(5,tup.length()); // replace "/pnfs/" with xrootd prefix
    const char* t = tup.c_str();
    inChain.Add(t);
  }

  // For summing up the POT totals from the Meta tree
  double sumPOTUsed=getTChainPOTXROOTD(inChain);

  int nFiles=g.gl_pathc;
  cout << "Added " << nFiles << " files of data" << endl;
  cout << "POT total used=" << sumPOTUsed << endl;
  globfree(&g);

  if(nFiles==0){
    cout << "No files added, nothing to do..." << endl;
    return;
  }
  
  TStopwatch ts;
  TFile* fout=new TFile(output, "RECREATE");
  cout << "Merging ana tree" << endl;
  fout->cd(); // Just in case the surrounding lines get separated
  inChain.Merge(fout, 32000, "keep SortBasketsByBranch");

  fout->cd();
  TTree* newMetaTree=new TTree("Meta", "Titles are stupid");
  newMetaTree->Branch("POT_Used", &sumPOTUsed);
  newMetaTree->Fill();
  newMetaTree->Write();
  ts.Stop();
  cout << "Merging time:" << endl;
  ts.Print();
}

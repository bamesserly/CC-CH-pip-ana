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
void mergeDataFiles(const char* inDirBase, const char* outDir, const char* tag="CCNuPionInc") {
  //******************************************************************
  //* Set Location of output files
  //******************************************************************
  TString output=TString::Format("%s/merged_CCPionAnaTool.root", outDir);

  //******************************************************************
  //* Load Input Ntuples
  //******************************************************************
  TChain inChain("CCNuPionInc");

  //TString inGlob(TString::Format("%s/??/??/??/??/MV_*%s*.root", inDirBase, tag));
  TString inGlob(TString::Format("%s/*/*/*/*/MV_*%s*.root", inDirBase, tag));

  cout << "Filename glob is " << inGlob << endl;
  cout << "Output filename is " << output << endl;

  glob_t g;
  glob(inGlob.Data(), 0, 0, &g);
  for(int i=0; i<(int)g.gl_pathc; ++i){
    inChain.Add(g.gl_pathv[i]);
  }

  // For summing up the POT totals from the Meta tree
  double sumPOTUsed=getTChainPOT(inChain);

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

//******************************************************************************
//*
//* Merges data tuples, run-by-run
//* Original author - Laura Fields, taken from Phil Rodrigues.
//*
//* This script is not maintained, nor used. Instead I merge entire data
//* playlists together. It's much faster, and can be looped on a gpvm under
//* normal circumstances in ~15-30 minutes.
//*
//******************************************************************************

#include "glob.h"

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
void mergeDataPlaylist(const char* inDirBase, const char* outDir, const char* playlist, const char* tag="MECAnaTool",
                       const char* treeName="MECAna") {
  //******************************************************************
  //* Set Location of output files
  //******************************************************************
  TString output=TString::Format("%s/%s_merged_%s.root", outDir, playlist, treeName);
  cout << "Output filename is " << output << endl;
  //******************************************************************
  //* Load Input Ntuples
  //******************************************************************
  TChain inChain(treeName);

  map<string, pair<int,int> > playlistRunRanges;

  // ----------------------------
  // Low energy playlists
  // ----------------------------
  playlistRunRanges["minerva1"]=make_pair(2000, 2397  );
  playlistRunRanges["minerva5"]=make_pair(2650, 2856  );
  playlistRunRanges["minerva7"]=make_pair(3131, 3160  );
  playlistRunRanges["minerva7a"]=make_pair(3086, 3098  );
  // 7 -> 7B happens in the middle of a run, so just put the runs in minerva7 to avoid double-counting
  // playlistRunRanges["minerva7B"]=make_pair(3152, 3160  );
  playlistRunRanges["minerva9"]=make_pair(3345, 3357  );
  playlistRunRanges["minerva13a"]=make_pair(3526, 3542  );
  playlistRunRanges["minerva13b"]=make_pair(3543, 3591  );
  playlistRunRanges["minerva13c"]=make_pair(3596,3816  );
  // Combine 13D and E for the same reason
  playlistRunRanges["minerva13d"]=make_pair(3816, 3901 );

  // ----------------------------
  // Medium energy playlists
  // ----------------------------
  playlistRunRanges["minervame1a"]=make_pair(6038,  10066);
  playlistRunRanges["minervame1b"]=make_pair(10068, 10128);
  playlistRunRanges["minervame1c"]=make_pair(10129, 10254);
  playlistRunRanges["minervame1d"]=make_pair(10255, 10662);
  playlistRunRanges["minervame1e"]=make_pair(16003, 16522);
  playlistRunRanges["minervame1f"]=make_pair(16523, 16843);
  
  int minrun=playlistRunRanges[playlist].first;
  int maxrun=playlistRunRanges[playlist].second;

  int nFiles=0;
  for(int run=minrun; run<=maxrun; ++run){
    string runStr(TString::Format("%08d", run));
    string runStrParts[4];
    for(int i=0; i<4; ++i) runStrParts[i]=runStr.substr(i*2, 2);
    TString inGlob(TString::Format("%s/%s/%s/%s/%s/MV_*%s_*_%s*.root",
                                   inDirBase,
                                   runStrParts[0].c_str(),
                                   runStrParts[1].c_str(),
                                   runStrParts[2].c_str(),
                                   runStrParts[3].c_str(),
                                   runStr.c_str(),
                                   tag));
    
    if(run==minrun){
      cout << "First filename glob is " << inGlob << endl;
    }

    glob_t g;
    glob(inGlob.Data(), 0, 0, &g);
    nFiles+=g.gl_pathc;
    for(int i=0; i<(int)g.gl_pathc; ++i){
      inChain.Add(g.gl_pathv[i]);
    }
    globfree(&g);
  }

  // For summing up the POT totals from the Meta tree
  double sumPOTUsed=getTChainPOT(inChain);

  cout << "Added " << nFiles << " files of data" << endl;
  cout << "POT total used=" << sumPOTUsed << endl;


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


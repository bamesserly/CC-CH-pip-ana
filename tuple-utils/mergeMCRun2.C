//******************************************************************************
//*
//* Merge mc tuples run-by-run
//* Original author - Laura Fields, taken from Phil Rodrigues.
//*
//* I used to use this because it made a fiducial cut on the truth tree,
//* removing events that cannot possibly be used, even after considering
//* systematics and efficiency correction. This significantly reduced the size
//* of the tuples.
//*
//* TTree* outTreeTruth=inChainTruth.CopyTree("mc_vtx[2]>5891 && mc_vtx[2]<8439");
//*
//* I stopped using this script because this fiducial truth cut was worked into
//* the official PlotUtils merging script.
//*
//******************************************************************************

#ifndef __CINT__
#include "glob.h"
#endif

#include "merge_common.h"

#include <iostream>
#include <string>

using namespace std;

//======================================================================
// Like mergeTrees, but keeps all the branches in CCQEAntiNuTool, and
// only merges one run.
void mergeMCRun2(const char* inDirBase, const char* outDir, int run, const char* tag="CCNuPionInc", const char* treeName="CCNuPionInc") {
  //******************************************************************
  //* Set Location of output files
  //******************************************************************
  TString output=TString::Format("%s/merged_%s_run%08d.root", outDir, tag, run);

  //******************************************************************
  //* Load Input Ntuples
  //******************************************************************
  TChain inChain(treeName);
  TChain inChainTruth("Truth");

  string runStr(TString::Format("%08d", run));
  string runStrParts[4];
  for(int i=0; i<4; ++i) runStrParts[i]=runStr.substr(i*2, 2);
  TString inGlob(TString::Format("%s/%s/%s/%s/%s/SIM_*%s_*_%s*.root",
                                 inDirBase,
                                 runStrParts[0].c_str(),
                                 runStrParts[1].c_str(),
                                 runStrParts[2].c_str(),
                                 runStrParts[3].c_str(),
                                 runStr.c_str(),
                                 tag));
  
  cout << "Filename glob is " << inGlob << endl;
  cout << "Output filename is " << output << endl;

  glob_t g;
  glob(inGlob.Data(), 0, 0, &g);

  cout << "Total files " << g.gl_pathc << ". Adding good files" << endl;

  int nFiles=0;
  for(int i=0; i<(int)g.gl_pathc; ++i){
    if(i%100==0) cout << i << " " << flush;
    const char* filename=g.gl_pathv[i];
    if(isGoodFile(filename)){
      inChain.Add(filename);
      inChainTruth.Add(filename);
      ++nFiles;
    }
    else{
      cout << "Skipping " << filename << endl;
    }
  }
  cout << endl;

  // For summing up the POT totals from the Meta tree
  double sumPOTUsed=getTChainPOT(inChainTruth, "POT_Used");
  double sumPOTTotal=getTChainPOT(inChainTruth, "POT_Total");
  int nFilesTotal=g.gl_pathc;
  cout << "Added " << nFiles << " files out of " << nFilesTotal << " in run " << run << endl;
  cout << "POT totals: Total=" << sumPOTTotal << " Used=" << sumPOTUsed << endl;
  globfree(&g);

  if(nFiles==0){
    cout << "No files added, nothing to do..." << endl;
    return;
  }
  
  TStopwatch ts;
  TFile* fout=new TFile(output, "RECREATE");
  cout << "Merging ana tree" << endl;
  setBranchStatuses(inChain);
  fout->cd(); // Just in case the surrounding lines get separated
  inChain.Merge(fout, 32000, "keep SortBasketsByBranch");

  cout << "Merging truth tree" << endl;
  setBranchStatuses(inChainTruth);
  fout->cd();
  TTree* outTreeTruth=inChainTruth.CopyTree("mc_vtx[2]>5891 && mc_vtx[2]<8439");
  outTreeTruth->Write();
  
  // inChainTruth.Merge(fout, 32000, "keep SortBasketsByBranch");

  fout->cd();
  TTree* newMetaTree=new TTree("Meta", "Titles are stupid");
  newMetaTree->Branch("POT_Used", &sumPOTUsed);
  newMetaTree->Branch("POT_Total", &sumPOTTotal);
  //if(!mc) newMetaTree->Branch("POT_Unanalyzable", &sumPOTUnanalyzable);
  newMetaTree->Fill();
  newMetaTree->Write();
  ts.Stop();
  cout << "Merging time:" << endl;
  ts.Print();
}


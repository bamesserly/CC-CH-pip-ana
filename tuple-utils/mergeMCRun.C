//******************************************************************************
//*
//* Don't know exactly what this script does. I've never used it.
//* Original author - Laura Fields, taken from Phil Rodrigues.
//*
//******************************************************************************

#include "glob.h"

#include "merge_common.h"

#include <iostream>
#include <string>

using namespace std;

//======================================================================
// Like mergeTrees, but keeps all the branches in CCQEAntiNuTool, and
// only merges one run.
void mergeMCRun(const char* inDirBase, const char* outDir, int run, const char* tag="MECAnaTool") {
  //******************************************************************
  //* Set Location of output files
  //******************************************************************
  TString output=TString::Format("%s/merged_MECAnaTool_run%08d.root", outDir, run);

  //******************************************************************
  //* Load Input Ntuples
  //******************************************************************
  TChain inChain("MECAna");
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
  // For summing up the POT totals from the Meta tree
  double sumPOTUsed=0;
  double sumPOTTotal=0;

  glob_t g;
  glob(inGlob.Data(), 0, 0, &g);
  for(int i=0; i<(int)g.gl_pathc; ++i){
    inChain.Add(g.gl_pathv[i]);
    inChainTruth.Add(g.gl_pathv[i]);

    double thisPOT=getResurrectionMCFilePOT(g.gl_pathv[i]);

    sumPOTUsed+=thisPOT;
    sumPOTTotal+=thisPOT;
  }
  int nFiles=g.gl_pathc;
  cout << "Added " << nFiles << " files of MC in run " << run << endl;
  cout << "POT totals: Total=" << sumPOTTotal << " Used=" << sumPOTUsed << endl;
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

  cout << "Merging truth tree" << endl;
  fout->cd();
  inChainTruth.Merge(fout, 32000, "keep SortBasketsByBranch");

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


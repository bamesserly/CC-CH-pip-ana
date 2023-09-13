#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>     
#include <map>
using namespace std;

//This does an event by event comparison of the branches in "names" and prints out the event IDs where they don't match.
void RecoilComparison(){
  TString oldFileName = "/pnfs/minerva/persistent/users/bmesserl/pions/20210307/merged/mc/ME1A/CCNuPionInc_mc_AnaTuple_run00110030_Playlist.root";
  TString newFileName = "/pnfs/minerva/persistent/users/granados/MADtuplas/merged/20230329/mc/ME1A/MasterAnaDev_mc_AnaTuple_run00110030_Playlist.root";

  cout << "opening: " << oldFileName << endl;
  TFile* oldFile = new TFile(oldFileName,"READ");
  cout << "opening: " << newFileName << endl;
  TFile* newFile = new TFile(newFileName,"READ");

  //TTree* oldTree = (TTree*)oldFile->Get("CCQENu");
  TTree* oldTree = (TTree*)oldFile->Get("CCNuPionInc");
  TTree* newTree = (TTree*)newFile->Get("MasterAnaDev");

  double newID;
  newTree->SetBranchAddress("eventID",&newID);
  int newEntries = newTree->GetEntries();

  map<double,int> indexNew;

  cout << newEntries << endl;
  for (int i=0; i < newEntries; ++i){
    if (i%(newEntries/100)==0) cout << i <<"/"<< newEntries << endl;
    newTree->GetEntry(i);
    if (newID < 0) continue;
    //cout << "NEW Non-negative ID" << endl;
    indexNew[newID]=i;
  }

  cout << "Setting NHit" << endl;
  double passive_em_ID_Old;
  oldTree->SetBranchAddress("part_response_recoil_passive_em_id",&passive_em_ID_Old);
  TH1D* h_passive_em_ID = new TH1D("passive_em_ID_Old","passive_em_ID Comparison;Diff.;Events",80,-20,20);

  double prong_Score_Old;
  oldTree->SetBranchAddress("prong_part_score",&prong_Score_Old);
  TH1D* h_prong_part_score = new TH1D("prong_part_score","prong_part_score;Diff.;Events",80,-20,20);

  int n_prongs_Old;
  oldTree->SetBranchAddress("n_prongs",&n_prongs_Old);
  TH1D* h_n_prongs = new TH1D("n_prongs","n_prongs MAD - CCNuPionInc;Diff.;Events",20,-10,10);

  int n_long_tracks_Old;
  oldTree->SetBranchAddress("n_long_tracks",&n_long_tracks_Old);
  TH1D* h_n_long_tracks = new TH1D("n_long_tracks","n_long_tracks MAD - CCNuPionInc;Diff.;Events",20,-10,10);

  int n_short_tracks_Old;
  oldTree->SetBranchAddress("n_short_tracks",&n_short_tracks_Old);
  TH1D* h_n_short_tracks = new TH1D("n_short_tracks","n_short_tracks MAD - CCNuPionInc;Diff.;Events",20,-10,10);


/*
  cout << "Stteing NVtx" << endl;
  int nVtxOld;
  oldTree->SetBranchAddress("event_vertex_time_diff_sz",&nVtxOld);

  cout << "Setting Vtx Diff" << endl;
  double vtxDiffOld[200];
  oldTree->SetBranchAddress("event_vertex_time_diff",&vtxDiffOld);
  TH1D* hCompareVtxDiff = new TH1D("hCompareVtxDiff","Vertex 0 Time Diff. Comparison;Diff.;Events",50,-25,25);

  cout << "Setting Recoil." << endl;
  double oldRecoil[200];
  oldTree->SetBranchAddress("recoil_summed_energy",&oldRecoil);
  TH1D* hCompareRecoil = new TH1D("hCompareRecoil","Recoil Summed Comparison;Percent Diff.;Events",201,-100.5,100.5);
  TH2D* h2DRecoil = new TH2D("h2DRecoil","Recoil Summed 2D Comparison;CCQENu;Percent Diff.",200,0,10000,200,-100,100);

  cout << "Setting Recoil Non" << endl;
  double oldRecNon;
  oldTree->SetBranchAddress("recoil_energy_nonmuon_nonvtx100mm",&oldRecNon);
  TH1D* hCompareRecNon = new TH1D("hCompareRecNon","Recoil NonVtx Comparison;Percent Diff.;Events",201,-100.5,100.5);
  TH2D* h2DRecNon = new TH2D("h2DRecNon","Recoil NonVtx 2D Comparison;CCQENu;Percent Diff.",200,0,10000,200,-100,100);

  cout << "Setting Recoil Vtx" << endl;
  double oldRecVtx;
  oldTree->SetBranchAddress("recoil_energy_nonmuon_vtx100mm",&oldRecVtx);
  TH1D* hCompareRecVtx = new TH1D("hCompareRecVtx","Recoil Vtx Comparison;Percent Diff.;Events",201,-100.5,100.5);
  TH2D* h2DRecVtx = new TH2D("h2DRecVtx","Recoil Vtx 2D Comparison;CCQENu;Percent Diff.",200,0,10000,200,-100,100);

  cout << "Added Comparison" << endl;
  TH1D* hCompareRecAdded = new TH1D("hCompareRecAdded","Recoil Vtx+Non-Vtx Comparison;Percent Diff.;Events",201,-100.5,100.5);
*/
  cout << "Setting New NHit" << endl;
  double passive_em_ID_New;
  newTree->SetBranchAddress("part_response_recoil_passive_em_id",&passive_em_ID_New);

  double prong_part_score_New;
  newTree->SetBranchAddress("prong_part_score",&prong_part_score_New);

  int n_prongs_New;
  newTree->SetBranchAddress("n_prongs",&n_prongs_New);

  int n_long_tracks_New;
  newTree->SetBranchAddress("n_long_tracks",&n_long_tracks_New);

  int n_short_tracks_New;
  newTree->SetBranchAddress("n_short_tracks",&n_short_tracks_New);


/*
  cout << "Setting New Vtx Diff" << endl;
  double vtxDiffNew[200];
  newTree->SetBranchAddress("event_vertex_time_diff",&vtxDiffNew);

  cout << "Setting New Recoil" << endl;
  double newRecoil[200];
  newTree->SetBranchAddress("recoil_summed_energy",&newRecoil);

  cout << "Setting New Recoil Non" << endl;
  double newRecNon;
  newTree->SetBranchAddress("recoil_energy_nonmuon_nonvtx100mm",&newRecNon);

  cout << "Setting New Recoil Vtx" << endl;
  double newRecVtx;
  newTree->SetBranchAddress("recoil_energy_nonmuon_vtx100mm",&newRecVtx);
*/
  double oldID;
  oldTree->SetBranchAddress("eventID",&oldID);
  int oldEntries = oldTree->GetEntries();
  TCanvas* c1 = new TCanvas();

  for (int i=0; i < oldEntries; ++i){
    if (i%(oldEntries/100)==0) cout << i << "/" << oldEntries << endl;
    oldTree->GetEntry(i);
    if (oldID < 0) continue;
//    if (oldRecoil[0]<=0) continue;
    double ID = oldID;

    if (indexNew[ID] ){
      newTree->GetEntry(indexNew[ID]);
      if (n_prongs_New - n_prongs_Old == 0.) h_passive_em_ID->Fill(passive_em_ID_New - passive_em_ID_Old);
      h_n_prongs->Fill(n_prongs_New - n_prongs_Old);
      h_n_long_tracks->Fill(n_long_tracks_New - n_long_tracks_Old);
      h_n_short_tracks->Fill(n_short_tracks_New - n_short_tracks_Old);
      //hCompareNHits->Fill(nHitsNew-nHitsOld);
      //for (int j=0; j<nVtxOld; ++j)hCompareVtxDiff->Fill(vtxDiffNew[j]-vtxDiffOld[j]);
/*      hCompareRecoil->Fill(100.0*(newRecoil[0]-oldRecoil[0])/oldRecoil[0]);
      double oldRecAdded = oldRecNon + oldRecVtx;
      double newRecAdded = newRecNon + newRecVtx;
      if (oldRecAdded > 0) hCompareRecAdded->Fill(100.0*(newRecAdded-oldRecAdded)/oldRecAdded);
      if (oldRecNon >0) hCompareRecNon->Fill(100.0*(newRecNon-oldRecNon)/oldRecNon);
      if (oldRecVtx >0) hCompareRecVtx->Fill(100.0*(newRecVtx-oldRecVtx)/oldRecVtx);
      h2DRecoil->Fill(oldRecoil[0],100.0*(newRecoil[0]-oldRecoil[0])/oldRecoil[0]);
      if (oldRecNon > 0) h2DRecNon->Fill(oldRecNon,100.0*(newRecNon-oldRecNon)/oldRecNon);
      if (oldRecVtx > 0) h2DRecVtx->Fill(oldRecVtx,100.0*(newRecVtx-oldRecVtx)/oldRecVtx);*/
    } 
    else continue;
  }
  h_n_prongs->Draw("HIST");
  c1->SetLogy();
  c1->Print("n_prongs_compare.png");

  TCanvas* c4 = new TCanvas();
  h_passive_em_ID->Draw("HIST");
  c4->SetLogy();
  c4->Print("passive_em_ID.png");

  TCanvas* c2 = new TCanvas();
  h_n_long_tracks->Draw("HIST");
  c2->SetLogy();
  c2->Print("n_long_tracks.png");

  TCanvas* c3 = new TCanvas();
  h_n_short_tracks->Draw("HIST");
  c3->SetLogy();
  c3->Print("n_short_tracks.png");



  cout << "DONE" << endl;
}

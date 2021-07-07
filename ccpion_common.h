#ifndef CCNUPION_COMMON_H
#define CCNUPION_COMMON_H

#include "util/util.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TSystem.h"
#include "TLorentzVector.h"
#include "TH2.h"
#include "TFile.h"
#include "TSpline.h"
#include "TChain.h"
#include "TRandom.h"

#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/MnvFluxConstraint.h"
#ifndef __CINT__
#include "PlotUtils/ChainWrapper.h"
#endif
#include "PlotUtils/MnvH1D.h"

// Should we enable the TTreeCache in GetChainWrapper?
bool gUseTTreeCache=true;

enum EProcessing
{
  kResurrection, kEroica, kNX
};

//==============================================================================

enum EPlaylist
{
  kAllPlaylists = 0xFF, // Just the FHC playlists. I'll explicitly ask for minerva5 when I want it
  kMinerva1     = 0x1,
  kMinerva7     = 0x2, // Includes minerva7B, because changeover was in middle of a run
  kMinerva7a    = 0x4,
  kMinerva9     = 0x8,
  kMinerva13a   = 0x10,
  kMinerva13b   = 0x20,
  kMinerva13c   = 0x40,
  kMinerva13d   = 0x80, // Includes minerva13E, for the same reason

  kMinerva5     = 0x100,

  kMinervaME1a  = 0x120,

  kNPlaylists   = 10 // Can I do this? It's the same value as an existing enum value. Let's guess yes
};

//==============================================================================

namespace ccpion_common {
  enum EDataMC
  {
    kData, kMC, kFakeData //, kMECMC
  };
}

//==============================================================================

const std::map<int, std::string>& GetPlaylistStrings()
{
  static std::map<int, std::string> playlistStrings;

  if (playlistStrings.empty()) {
    playlistStrings[kMinerva1]    = "minerva1";
    playlistStrings[kMinerva5]    = "minerva5";
    playlistStrings[kMinerva7]    = "minerva7";
    playlistStrings[kMinerva7a]   = "minerva7a";
    playlistStrings[kMinerva9]    = "minerva9";
    playlistStrings[kMinerva13a]  = "minerva13a";
    playlistStrings[kMinerva13b]  = "minerva13b";
    playlistStrings[kMinerva13c]  = "minerva13c";
    playlistStrings[kMinerva13d]  = "minerva13d";
    playlistStrings[kMinervaME1a] = "minervame1a";
  }

  return playlistStrings;
}

//==============================================================================

#ifndef __CINT__
std::string GetPlaylistFile(std::string plist, bool is_mc, bool use_xrootd = true) {
  const std::string processing_date = "20210520";//20210306
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  std::string topdir = "/minerva/data/users/granados/MAD_ana_plists/" + processing_date;
  std::string playlist_file = use_xrootd ?
      Form("%s/%s_%s_xrootd_plist.txt", topdir.c_str(), is_mc_str.c_str(), plist.c_str()) :
      Form("%s/%s_%s_plist.txt",        topdir.c_str(), is_mc_str.c_str(), plist.c_str());
  return playlist_file;
}
#endif // __CINT__

#ifndef __CINT__
std::string GetPlaylistFileCCPi(std::string plist, bool is_mc, bool use_xrootd = true) {
  const std::string processing_date = "20200713";
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  std::string topdir = "/minerva/data/users/bmesserl/MECCCHpip_ana_plists/" + processing_date;
  std::string playlist_file = use_xrootd ?
      Form("%s/%s_%s_xrootd_plist.txt", topdir.c_str(), is_mc_str.c_str(), plist.c_str()) :
      Form("%s/%s_%s_plist.txt",        topdir.c_str(), is_mc_str.c_str(), plist.c_str());
  return playlist_file;
}
#endif // __CINT__

//==============================================================================


#ifndef __CINT__
// playlist-by-playlist + grid + no-grid version
PlotUtils::ChainWrapper* 
    GetChainWrapperPtr(ccpion_common::EDataMC datamc, 
                       std::string plist       = "ME1A",
                       const bool do_grid      = false,
                       const char* tree_name   = "MasterAnaDev") {
  // Use a pointer because I can't be bothered working out whether my
  // copy constructor does what I want. This'll leak but irrelevant
  // since its lifetime is the whole process anyway
  PlotUtils::ChainWrapper* chw=new PlotUtils::ChainWrapper(tree_name);
  if (gUseTTreeCache) {
    chw->GetChain()->SetCacheSize(100000000);
    chw->GetChain()->SetCacheLearnEntries(1000);
  }

  // upper case
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);

  //const TString latest_processing("/pnfs/minerva/persistent/users/bmesserl/pions/20190824");// thesis
  const TString latest_processing("/pnfs/minerva/persistent/users/granados/MADtuplas/20210520");//MAD tuplas

  TString baseDir, fileGlob;
  if (do_grid) {
    baseDir = TString::Format("%s", gSystem->Getenv("CONDOR_DIR_INPUT"));
    baseDir+="/";
    fileGlob = TString::Format("%s/*CC*.root", baseDir.Data());
    std::cout << "Adding " << fileGlob << std::endl;
    chw->Add(fileGlob.Data());
  } // on grid
  else {
    baseDir = TString(latest_processing);
    baseDir+="/";
    baseDir+="/merged";
    if (datamc==ccpion_common::kData)
      baseDir+="/data";
    else if (datamc==ccpion_common::kMC || datamc==ccpion_common::kFakeData)
      baseDir+="/mc";

    std::vector<std::string> plists = {"ME1A", "ME1B", "ME1C", "ME1D", "ME1E",
                                       "ME1F", "ME1G",
                                       "ME1L", "ME1M", "ME1N",
                                       "ME1O", "ME1P"
                                      };

    for (auto p : plists) {
      if (plist == p || (datamc==ccpion_common::kData && plist == "ALL")) {
        fileGlob = TString((TString::Format("%s/%s/*CC*.root", baseDir.Data(), p.c_str())));
        std::cout << "Adding " << fileGlob << std::endl;
        chw->Add(fileGlob.Data());
      }
    }
  } // not on grid

  return chw;
}

#endif

//==============================================================================
#ifndef __CINT__
// xrootd, individual file, should work on and off grid
PlotUtils::ChainWrapper* 
    GetChainWrapperPtr(std::string xrootd_full_path,
                       const char* tree_name   = "MasterAnaDev") {
  // Use a pointer because I can't be bothered working out whether my
  // copy constructor does what I want. This'll leak but irrelevant
  // since its lifetime is the whole process anyway
  PlotUtils::ChainWrapper* chw=new PlotUtils::ChainWrapper(tree_name);

  if (gUseTTreeCache) {
    chw->GetChain()->SetCacheSize(100000000);
    chw->GetChain()->SetCacheLearnEntries(1000);
  }

  std::cout << "Adding " << xrootd_full_path << std::endl;
  chw->Add(xrootd_full_path);

  return chw;
}

#endif
//==============================================================================

#ifndef __CINT__
double GetChainPOTFromFileList(const std::vector<std::string>& file_list) {
  double pot = 0.;
  for(uint i=0; i<file_list.size(); ++i) {
    TFile f(file_list[i].c_str());
    assert(!f.IsZombie());
    TTree* t=(TTree*)f.Get("Meta");
    assert(t);
    assert(t->GetEntries()==1);
    t->GetEntry(0);

    TLeaf* lUsed=t->GetLeaf("POT_Used");

    pot+=lUsed->GetValue();
  }
  return pot;
}


double GetChainPOT(PlotUtils::ChainWrapper& chw) {
  return GetChainPOTFromFileList(chw.GetListOfFiles());
}


double GetChainPOT(PlotUtils::ChainWrapper* chw) {
  return GetChainPOTFromFileList(chw->GetListOfFiles());
}
#endif

#endif

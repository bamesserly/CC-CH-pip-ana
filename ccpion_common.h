#ifndef ccpion_common_h
#define ccpion_common_h

#include <algorithm>  // transform
#include <string>

#include "TString.h"  // Form

std::string GetPlaylistFile(std::string plist, bool is_mc,
                            bool use_xrootd = true) {

  // small sample with michel branches
  // /pnfs/minerva/persistent/users/granados/MADtuplas/merged/mc/Michel/ME1A/MasterAnaDev_mc_AnaTuple_run00110000_Playlist.root 
  // /pnfs/minerva/persistent/users/granados/MADtuplas/merged/data/Michel/

  // const std::string processing_date = "20200713"; // new short tracking branches
  const std::string processing_date = "20211115";  // new recoil energy branches
  //  const std::string processing_date = "test"; // For test with small MAD tuplas
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  std::string topdir =
      is_mc ? "/minerva/data/users/granados/MAD_ana_plists/"
            : "/minerva/data/users/granados/MAD_ana_plists/";  // correct merging method
  topdir += processing_date;
  std::string playlist_file =
      use_xrootd ? Form("%s/%s_%s_xrootd_plist.txt", topdir.c_str(),
                        is_mc_str.c_str(), plist.c_str())
                 : Form("%s/%s_%s_plist.txt", topdir.c_str(), is_mc_str.c_str(),
                        plist.c_str());
  return playlist_file;
}

std::string GetPlaylistFileCCPi(std::string plist, bool is_mc,
                                bool use_xrootd = true) {
  const std::string processing_date = "20210307";
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  std::string topdir =
      "/minerva/data/users/bmesserl/MECCCHpip_ana_plists/" + processing_date;
  std::string playlist_file =
      use_xrootd ? Form("%s/%s_%s_xrootd_plist.txt", topdir.c_str(),
                        is_mc_str.c_str(), plist.c_str())
                 : Form("%s/%s_%s_plist.txt", topdir.c_str(), is_mc_str.c_str(),
                        plist.c_str());
  return playlist_file;
}

#endif  // ccpion_common_h

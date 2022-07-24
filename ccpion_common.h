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
  // /pnfs/minerva/persistent/users/granados/MADtuplas/merged/20220101/mc/ME1A/
  // /pnfs/minerva/persistent/users/granados/MADtuplas/merged/data/Michel/ME1A/ 
  // /minerva/data/users/bmesserl/MECCCHpip_ana_plists/20220101/

  // const std::string processing_date = "20200713"; // new short tracking branches
  // const std::string processing_date = "20211012"; // new recoil energy branches
  const std::string processing_date = "20220101"; // includes mehreen's michel branches
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  std::string topdir = "/minerva/data/users/bmesserl/MECCCHpip_ana_plists/";
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

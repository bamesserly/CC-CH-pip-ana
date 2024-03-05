#ifndef ccpion_common_h
#define ccpion_common_h

#include <algorithm>  // transform
#include <string>

#include "TString.h"  // Form

std::string GetPlaylistFile(std::string plist, bool is_mc,
                            bool use_xrootd = true) {
  // const std::string processing_date = "20200713"; // new short tracking
  // branches
  const std::string processing_date =
      "InPrepTop4";  // new recoil energy branches
  std::string old = Form("/minerva/data/users/granados/MAD_ana_plists/%s",
                         processing_date.c_str());

  std::transform(plist.begin(), plist.end(), plist.begin(),
                 ::tolower);  // Convert entire plisting to lowercase
  plist.back() =
      std::toupper(plist.back());  // Convert last character to uppercase
  const std::string mc_plist = Form(
      "/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/"
      "minerva%s_MC.txt",
      plist.c_str());
  const std::string data_plist = Form(
      "/pnfs/minerva/persistent/DataPreservation/p4/Data/minerva%s_Data.txt",
      plist.c_str());
  return is_mc ? mc_plist : data_plist;
}

std::string GetTestPlaylist(bool is_mc) {
  return is_mc ? "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/cache/"
                 "mc_ME1B_p4.txt"
               : "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/cache/"
                 "data_ME1B_p4.txt";
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

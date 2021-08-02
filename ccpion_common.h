#ifndef ccpion_common_h
#define ccpion_common_h

#include <string>
#include <algorithm> // transform
#include "TString.h" // Form

std::string GetPlaylistFile(std::string plist, bool is_mc, bool use_xrootd = true) {
#ifndef __CINT__
  //const std::string processing_date = "20200713"; // new short tracking branches
  const std::string processing_date = "20210730"; // new recoil energy branches
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  std::string topdir = "/minerva/data/users/granados/MAD_ana_plists/" + processing_date;
  std::string playlist_file = use_xrootd ?
      Form("%s/%s_%s_xrootd_plist.txt", topdir.c_str(), is_mc_str.c_str(), plist.c_str()) :
      Form("%s/%s_%s_plist.txt",        topdir.c_str(), is_mc_str.c_str(), plist.c_str());
  return playlist_file;
#endif // __CINT__
}

#endif // ccpion_common_h

#ifndef CCPiMacroUtil_h
#define CCPiMacroUtil_h

//==============================================================================
// Extension of PlotUtils::MacroUtil
// * Add do_[data,mc,truth,systematics] bools
// * Add SignalDefinition
// * Add a member CVUniverse and MCReco and Truth universes (should be in PU)
// * Extend PrintMacroConfiguration to print all the above
// Helper functions:
// SetupLoop
//==============================================================================
#include <cassert>

#include "Constants.h"  // EDataMC for the SetupLoop function
#include "PlotUtils/MacroUtil.h"
#include "SignalDefinition.h"

namespace CCPi {

std::string GetPlaylistFile(std::string plist, const bool is_mc,
                            const bool do_test_playlist,
                            const bool use_xrootd) {
  if (do_test_playlist) {
    return is_mc ? "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/cache/"
                   "mc_ME1A_p3.txt"
                 : "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/cache/"
                   "data_ME1A_p3.txt";
  }
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  const std::string processing_date = "production_p3";
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::string topdir = is_mc ? "/minerva/data/users/granados/MAD_ana_plists/"
                             : "/minerva/data/users/granados/MAD_ana_plists/";
  topdir += processing_date;
  std::string playlist_file =
      use_xrootd ? Form("%s/%s_%s_xrootd_plist.txt", topdir.c_str(),
                        is_mc_str.c_str(), plist.c_str())
                 : Form("%s/%s_%s_plist.txt", topdir.c_str(), is_mc_str.c_str(),
                        plist.c_str());
  return playlist_file;
}

class MacroUtil : public PlotUtils::MacroUtil {
 public:
  // Data
  MacroUtil(const int signal_definition, const std::string& data_file_list,
            const std::string& plist, const bool is_grid);

  // MC (and Truth)
  MacroUtil(const int signal_definition, const std::string& mc_file_list,
            const std::string& plist, const bool do_truth, const bool is_grid,
            const bool do_systematics);

  // Data, MC (and Truth)
  MacroUtil(const int signal_definition, const std::string& mc_file_list,
            const std::string& data_file_list, const std::string& plist,
            const bool do_truth, const bool is_grid, const bool do_systematics);

  bool m_do_data;
  bool m_do_mc;
  bool m_do_truth;
  bool m_do_systematics;
  bool m_is_grid;
  SignalDefinition m_signal_definition;
  CVUniverse* m_data_universe;
  UniverseMap m_error_bands;
  UniverseMap m_error_bands_truth;
  double m_pot_scale;  // For now, only used in xsecDataFromFile
  std::string m_name;
  void PrintMacroConfiguration(std::string macro_name = "") override;

 private:
  void Init();
  void InitSystematics();
};  // class MacroUtil
}  // namespace CCPi

void SetupLoop(const EDataMCTruth& type, const CCPi::MacroUtil& util,
               bool& is_mc, bool& is_truth, Long64_t& n_entries);

#endif  // CCPiMacroUtil_h

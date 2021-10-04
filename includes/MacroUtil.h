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

#include "PlotUtils/MacroUtil.h"
#include "SignalDefinition.h"
#include "Constants.h" // EDataMC for the SetupLoop function

namespace CCPi {
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
  SignalDefinition m_signal_definition;
  CVUniverse* m_data_universe;
  UniverseMap m_error_bands;
  UniverseMap m_error_bands_truth;
  double m_pot_scale; // For now, only used in xsecDataFromFile
//#ifndef __CINT__
  void PrintMacroConfiguration(std::string macro_name = "") override;
//#endif

 private:
  void Init(const int signal_definition);
  void InitSystematics();
};
}  // namespace CCPi

void SetupLoop(const EDataMCTruth& type, const CCPi::MacroUtil& util,
               bool& is_mc, bool& is_truth, Long64_t& n_entries);

#endif // CCPiMacroUtil_h

#ifndef CCPiMacroUtil_cxx
#define CCPiMacroUtil_cxx

#include "MacroUtil.h"
#include "Systematics.h" // GetSystematicUniversesMap
#include "myPlotStyle.h" // Load my plot style in Init

// CTOR Data
CCPi::MacroUtil::MacroUtil(const int signal_definition,
                           const std::string& data_file_list,
                           const std::string& plist, const bool is_grid)
    : PlotUtils::MacroUtil("MasterAnaDev", data_file_list, plist),
      m_do_data(true),
      m_do_mc(false),
      m_do_truth(false),
      m_do_systematics(false),
      m_is_grid(is_grid) {
  Init(signal_definition);
}

// CTOR MC (and Truth)
CCPi::MacroUtil::MacroUtil(const int signal_definition,
                           const std::string& mc_file_list,
                           const std::string& plist, const bool do_truth,
                           const bool is_grid, const bool do_systematics)
    : PlotUtils::MacroUtil("MasterAnaDev", mc_file_list, plist, do_truth),
      m_do_data(false),
      m_do_mc(true),
      m_do_truth(do_truth),
      m_do_systematics(do_systematics),
      m_is_grid(is_grid) {
  Init(signal_definition);
}

// CTOR Data, MC (and Truth)
CCPi::MacroUtil::MacroUtil(const int signal_definition,
                           const std::string& mc_file_list,
                           const std::string& data_file_list,
                           const std::string& plist, const bool do_truth,
                           const bool is_grid, const bool do_systematics)
    : PlotUtils::MacroUtil("MasterAnaDev", mc_file_list, data_file_list, plist,
                           do_truth),
      m_do_data(false),
      m_do_mc(true),
      m_do_truth(do_truth),
      m_do_systematics(do_systematics),
      m_is_grid(is_grid) {
  Init(signal_definition);
}

// Extend
void CCPi::MacroUtil::PrintMacroConfiguration(std::string macro_name) {
  PlotUtils::MacroUtil::PrintMacroConfiguration(macro_name);
  std::cout << "\n** Fill truth = " << m_do_truth
            << "\n** is grid = " << m_is_grid
            << "\n** Do full systematics = " << m_do_systematics
            << "\n** UseNuEConstraint = "
            << MinervaUniverse::UseNuEConstraint()
            << "\n** AnalysisNuPDG = " << MinervaUniverse::GetAnalysisNuPDG()
            << "\n** Use NonResPi Weight as CV = "
            << MinervaUniverse::UseNonResPiReweight()
            << "\n** NFluxUniverses = "
            << MinervaUniverse::GetNFluxUniverses() << "\n\n";
}

// Private/internal
void CCPi::MacroUtil::Init(const int signal_definition) {
  m_signal_definition = static_cast<SignalDefinition>(signal_definition);
  myPlotStyle();
  TH1::SetDefaultSumw2();
  std::cout << "\n== Data POT " << m_data_pot << ", MC POT " << m_mc_pot
            << "\n";
  InitSystematics();
}

// Call GetSystematicsUniversesMap from Systematics.h
// Initialize all the booleans that DCVU asks for
void CCPi::MacroUtil::InitSystematics() {
  MinervaUniverse::SetZExpansionFaReweight(false); // before genie

  m_data_universe = new CVUniverse(m_data);

  m_error_bands =
      systematics::GetSystematicUniversesMap(m_mc, false, m_do_systematics);

  m_error_bands_truth =
      systematics::GetSystematicUniversesMap(m_truth, true, m_do_systematics);

  // Set Various Constant Parameters
  using namespace PlotUtils;
  MinervaUniverse::SetNuEConstraint(CCNuPionIncConsts::kUseNueConstraint);
  MinervaUniverse::SetAnalysisNuPDG(CCNuPionIncConsts::kAnaNuPDG);
  MinervaUniverse::SetNonResPiReweight(CCNuPionIncConsts::kUseNonResPiWgt);
  MinervaUniverse::SetNFluxUniverses(CCNuPionIncConsts::kNFluxUniverses);
  MinervaUniverse::SetDeuteriumGeniePiTune(false);

  // Set playlist -- for systematics, flux, and other stuff(?)
  // If we're only doing data, we don't care what playlist FRW wants to use
  // (Indeed, this further helps us because we want to loop over ALL data in
  // one loop)
  std::string plist = "minerva" + m_plist_string;
  std::transform(plist.begin(), plist.end(), plist.begin(), ::tolower);
  if (plist == "minervaall"){
    std::cout << "WARNING: Systematics not set up correctly with \'all\' ";
    std::cout << "playlist.\n";

    std::cout << "    Defaulting to me1a.\n";

    std::cout << "    If this is mc or truth processing, this is definitely ";
    std::cout << "bad, and you'll get the wrong answers.\n";

    std::cout << "    If this is data or xsec calculation, this is probably ";
    std::cout << "bad, and should be fixed.\n";

    std::cout << "    TODO: figure out which systs, universe calcs, flux ";
    std::cout << "things, or data/xsec calculations depend on playlist.\n";

    plist = "minervame1a";
  }
  MinervaUniverse::SetPlaylist(plist);
}

// Helper -- maybe this belongs somewhere else
void SetupLoop(const EDataMCTruth& type, const CCPi::MacroUtil& util,
               bool& is_mc, bool& is_truth, Long64_t& n_entries) {
  switch (type) {
    case kData:
      is_mc = false;
      is_truth = false;
      n_entries = util.GetDataEntries();
      std::cout << "*** Looping Data ***\n";
      break;
    case kMC:
      is_mc = true;
      is_truth = false;
      n_entries = util.GetMCEntries();
      std::cout << " *** Looping MC ***\n";
      break;
    case kTruth:
      is_mc = true;
      is_truth = true;
      n_entries = util.GetTruthEntries();
      std::cout << "*** Looping Truth ***\n";
      break;
    default:
      std::cerr << "MacroUtils::SetupLoop failure\n";
  }
}

#endif // CCPiMacroUtil_cxx

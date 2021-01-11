#ifndef CCPiMacroUtil_cxx
#define CCPiMacroUtil_cxx

#include "MacroUtil.h"
#include "Systematics.h" // GetSystematicUniversesMap
#include "../util/plot/myPlotStyle.h" // Load my plot style in Init

// CTOR Data
CCPi::MacroUtil::MacroUtil(const int signal_definition,
                           const std::string& data_file_list,
                           const std::string& plist, const bool is_grid)
    : PlotUtils::MacroUtil("CCNuPionInc", data_file_list, plist, is_grid),
      m_do_data(true),
      m_do_mc(false),
      m_do_truth(false),
      m_do_systematics(false) {
  Init(signal_definition);
}

// CTOR MC (and Truth)
CCPi::MacroUtil::MacroUtil(const int signal_definition,
                           const std::string& mc_file_list,
                           const std::string& plist, const bool do_truth,
                           const bool is_grid, const bool do_systematics)
    : PlotUtils::MacroUtil("CCNuPionInc", mc_file_list, plist, do_truth,
                           is_grid),
      m_do_data(false),
      m_do_mc(true),
      m_do_truth(do_truth),
      m_do_systematics(do_systematics) {
  Init(signal_definition);
}

// CTOR Data, MC (and Truth)
CCPi::MacroUtil::MacroUtil(const int signal_definition,
                           const std::string& mc_file_list,
                           const std::string& data_file_list,
                           const std::string& plist, const bool do_truth,
                           const bool is_grid, const bool do_systematics)
    : PlotUtils::MacroUtil("CCNuPionInc", mc_file_list, data_file_list, plist,
                           do_truth, is_grid),
      m_do_data(false),
      m_do_mc(true),
      m_do_truth(do_truth),
      m_do_systematics(do_systematics) {
  Init(signal_definition);
}

// Extend
void CCPi::MacroUtil::PrintMacroConfiguration(std::string macro_name) {
  PlotUtils::MacroUtil::PrintMacroConfiguration(macro_name);
  std::cout << "\n** Fill truth = " << m_do_truth
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
  // If we're only doing data, we don't care what playlist FRW wants to use
  // (Indeed, this further helps us because we want to loop over ALL data in
  // one loop)
  if (m_do_mc || m_do_truth)
    MinervaUniverse::SetPlaylist("minerva" + m_plist_string);
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
  }
}

#endif // CCPiMacroUtil_cxx
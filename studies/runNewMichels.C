//==============================================================================
// Template for a generic loop data/mc and plot script.
// Assume that NO systematics are being analyzed.
// Good for stacked histograms, branch validation, residuals, etc.
//==============================================================================
#include <iostream>
#include <vector>

#include "ccpion_common.h" // GetPlaylistFile
#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/MacroUtil.h"
#include "includes/HadronVariable.h"
#include "includes/Variable.h"
#include "plotting_functions.h"

#include "includes/common_functions.h"

// Forward declare my variables because we're hiding the header.
class Variable;
class HadronVariable;

namespace run_new_michels{
//==============================================================================
// Do some event processing (e.g. make cuts, get best pion) and fill hists
//==============================================================================
void FillVars(CCPiEvent& event, const std::vector<Variable*>& variables) {
  const CVUniverse* universe = event.m_universe;
  const double wgt = event.m_weight;
  const bool is_mc = event.m_is_mc;
  const SignalDefinition sd = event.m_signal_definition;

  if (universe->ShortName() != "cv") return;

  event.m_passes_cuts             = PassesCuts(event, event.m_is_w_sideband);
  event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
  if(event.m_passes_cuts)
    ccpi_event::FillStackedHists(event, variables);
}

//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  typedef Variable Var;
  typedef HadronVariable HVar;
  HVar* thetapi_deg = new HVar("thetapi_deg", "#theta_{#pi}", "deg", CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);
  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"), &CVUniverse::GetPmu);
  std::vector<Var*> variables = {thetapi_deg, pmu};
  return variables;
}
} // namespace run_new_michels

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFill(const CCPi::MacroUtil& util, CVUniverse* universe,
                        const EDataMCTruth& type,
                        std::vector<Variable*>& variables) {

  std::cout << "Loop and Fill CutVars\n";
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  //for(Long64_t i_event=0; i_event < n_entries; ++i_event){
  for(Long64_t i_event=0; i_event < 5000; ++i_event){
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
    universe->SetEntry(i_event);

    // For mc, get weight, check signal, and sideband
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    std::tie(event.m_passes_cuts, event.m_is_w_sideband, event.m_reco_pion_candidate_idxs) = PassesCuts(event);

    // WRITE THE FILL FUNCTION
    //run_new_michels::FillVars(event, variables);
  } // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runNewMichels(std::string plist = "ME1A") {
  //=========================================
  // Input tuples
  //=========================================
    bool is_mc = true;
    std::string mc_file_list, data_file_list;
    bool use_xrootd = true;
    mc_file_list = GetPlaylistFile(plist, is_mc, use_xrootd);
    is_mc = false;
    data_file_list = GetPlaylistFile(plist, is_mc, use_xrootd);
    
  //=========================================
  // Init macro utility
  //=========================================
    const int signal_definition_int = 0;
    const std::string macro("runNewMichels");
    const bool is_grid = false;
    const bool do_truth = false;
    const bool do_systematics = false;

    CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list, plist, do_truth, is_grid, do_systematics);
    util.PrintMacroConfiguration(macro);

  //=========================================
  // Get variables and initialize their hists
  //=========================================
  std::vector<Variable*> variables = run_new_michels::GetVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  //=========================================
  // Loop and Fill
  //=========================================
  LoopAndFill(util, util.m_error_bands.at("cv").at(0), kMC,   variables);
/*
  LoopAndFill(util, util.m_data_universe,              kData, variables);

  for (auto v : variables) {
    std::string tag = v->Name();
    double ymax = -1;
    bool do_bwn = true;
    std::cout << "Plotting" << std::endl;
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kS), util.m_data_pot, util.m_mc_pot, util.m_signal_definition, v->Name(),"SSB", ymax, do_bwn);
  }
*/

  std::cout << "Success" << std::endl;
}

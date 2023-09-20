//==============================================================================
// Template for a generic loop data/mc and plot script.
// Assume that NO systematics are being analyzed.
// Good for stacked histograms, branch validation, residuals, etc.
//==============================================================================
#include <iostream>
#include <vector>

#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/MacroUtil.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
//#include "plotting_functions.h"

/*
namespace run_study_template {
//==============================================================================
// Do some event processing (e.g. make cuts, get best pion) and fill hists
//==============================================================================
void FillVars(CCPiEvent& event, const std::vector<Variable*>& variables) {
  CVUniverse* universe = event.m_universe;
  const double wgt = event.m_weight;
  const bool is_mc = event.m_is_mc;
  const SignalDefinition sd = event.m_signal_definition;

  if (universe->ShortName() != "cv") return;

  // Process Event

  // Check cuts, check is_sideband, and get pion candidates
  PassesCutsInfo cv_cuts_info = PassesCuts(event);

  // Set all of this info to the event and/or universe
  std::tie(event.m_passes_cuts,
           event.m_is_w_sideband,
           event.m_passes_all_cuts_except_w,
           event.m_reco_pion_candidate_idxs) = cv_cuts_info.GetAll();
  event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
  universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);

  // Need to re-call this because the node cut efficiency systematic
  // needs a pion candidate to calculate its weight.
  event.m_weight = universe->GetWeight();

  if (event.m_passes_cuts) ccpi_event::FillStackedHists(event, variables);
}

//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  typedef Variable Var;
  typedef HadronVariable HVar;
  HVar* thetapi_deg =
      new HVar("thetapi_deg", "#theta_{#pi}", "deg",
               CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);
  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);
  std::vector<Var*> variables = {thetapi_deg, pmu};
  return variables;
}

}  // namespace run_study_template

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFill(const CCPi::MacroUtil& util, CVUniverse* universe,
                 const EDataMCTruth& type, std::vector<Variable*>& variables) {
  std::cout << "Loop and Fill CutVars\n";
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);

    // For mc, get weight, check signal, and sideband
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    // WRITE THE FILL FUNCTION
    run_study_template::FillVars(event, variables);
  }  // events
  std::cout << "*** Done ***\n\n";
}
*/

void LoopAndFill(CVUniverse* universe,
                 const Long64_t n_entries, const EDataMCTruth& type,
                 const SignalDefinition& signal_definition,
                 std::vector<VariableBase*>& variables) {
  bool is_mc = type == kMC || type == kTruth;
  bool is_truth = type == kTruth;

  std::cout << "MC? Truth?\n";
  std::cout << is_mc << "  " << is_truth << "\n";

  universe->SetTruth(is_truth);

  for (Long64_t i_event = 0; i_event < 100; ++i_event) {
    if (i_event % (n_entries / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;

    // Save vertical-only universe info across universes for optimization --
    // there's no need to recheck vert univ cuts.
    VertUniverseInfo vertical_universe_info;

    universe->SetEntry(i_event);

    std::cout << universe->ShortName() << "\n";

    // CCPiEvent keeps track of lots of event properties
    CCPiEvent event(is_mc, is_truth, signal_definition, universe);

    // Macro-level event reco (computationally intensive).
    // 
    // Construct michels/pions and check cuts.
    //
    // As we fail cuts (and we're not sideband either), don't waste time
    // continuing to process the event.
    LoopStatusCode status;
    std::tie(vertical_universe_info, status) = event.Process(vertical_universe_info);

    if (status == LoopStatusCode::SKIP)
      continue;

    for (auto v : variables) {
      std::cout << v->Name() << "  " << v->GetValue(event) << "\n";
    }

    //// Fill reco -- enforce cuts, internally update the histograms owned by variables
    //ccpi_event::FillRecoEvent(event, variables);
  } // entries 
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runRefactorTest(std::string plist = "ME1L") {
  //=========================================
  // Input tuples
  //=========================================
  bool is_mc = true;
  std::string mc_file_list, data_file_list;
  // mc_file_list = GetPlaylistFile(plist, is_mc);
  // is_mc = false;
  // data_file_list = GetPlaylistFile(plist, is_mc);

  mc_file_list = GetTestPlaylist(is_mc);
  is_mc = false;
  data_file_list = GetTestPlaylist(is_mc);

  //=========================================
  // Init macro utility
  //=========================================
  const int signal_definition_int = SignalDefinition::OnePiTracked().m_id;
  const bool is_grid = false;
  const bool do_truth = false;
  const bool do_systematics = false;

  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.m_name = "runRefactorTest";
  util.PrintMacroConfiguration();

  EventVariable* thetapi_deg = new EventVariable("thetapi_deg", "#theta_{#pi}", "deg", CCPi::GetBinning("thetapi_deg"), &CCPiEvent::GetDummyVar);
  CVUVariable* pmu = new CVUVariable("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"), &CVUniverse::GetPmu);

  std::vector<VariableBase*> variables = {thetapi_deg, pmu};
  
  //std::vector<Variable*> variables = run_study_template::GetVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  LoopAndFill(util.m_data_universe, util.GetDataEntries(), kData, util.m_signal_definition, variables);
  LoopAndFill(util.m_error_bands.at("cv").at(0), util.GetMCEntries(),  kMC, util.m_signal_definition, variables);

  /*
  for (auto v : variables) {
    std::string tag = v->Name();
    double ymax = -1;
    bool do_bwn = true;
    std::cout << "Plotting" << std::endl;
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kS),
               util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
               v->Name(), "SSB", ymax, do_bwn);
  }
  */

  std::cout << "Success" << std::endl;
}

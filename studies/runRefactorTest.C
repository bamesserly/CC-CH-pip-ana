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
#include "includes/common_var_functions.h"
#include "studies/plotting_functions.h"

void FillStackedHists(const CCPiEvent& event, Variable* v, double fill_val) {
  if (!event.m_is_mc) return;

  v->GetStackComponentHist(GetFSParticleType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetChannelType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetHadronType(*event.m_universe, 0))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetNPionsType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetNPi0Type(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetNPipType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetSignalBackgroundType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetWSidebandType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetMesonBackgroundType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetWBackgroundType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetTruthWType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetCoherentType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);
}

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

void LoopAndFill(CVUniverse* universe, const Long64_t n_entries,
                 const EDataMCTruth& type,
                 const SignalDefinition& signal_definition,
                 std::vector<Variable*>& variables) {
  bool is_mc = type == kMC || type == kTruth;
  bool is_truth = type == kTruth;

  std::cout << "MC? Truth?\n";
  std::cout << is_mc << "  " << is_truth << "\n";

  universe->SetTruth(is_truth);

  int n_pass = 0;

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % (n_entries / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;

    // Save vertical-only universe info across universes for optimization --
    // there's no need to recheck vert univ cuts.
    VerticalUniverseInfo vertical_universe_info;

    universe->SetEntry(i_event);

    // CCPiEvent keeps track of lots of event properties
    CCPiEvent event(is_mc, is_truth, signal_definition, universe);

    // Get untracked michels reco, is the SD calls for it.
    //
    // Do this only once for all universes -- it's not impacted by different
    // universes.
    if (signal_definition.m_do_untracked_michel_reco)
      std::tie(event.m_untracked_michels, event.m_untracked_michels_pass) =
          GetUntrackedMichels(universe);

    // Macro-level event reco (computationally intensive).
    //
    // Construct michels/pions and check cuts.
    //
    // As we fail cuts (and we're not sideband either), don't waste time
    // continuing to process the event.
    LoopStatusCode status;
    std::tie(vertical_universe_info, status) =
        event.Process(vertical_universe_info);

    if (status == LoopStatusCode::SKIP) continue;

    if (event.m_passes_cuts) n_pass++;

    // Current study: events that have made it this far pass basic cuts and
    // have at least one good michel of either reco method. No W cut, no track
    // cuts being applied.

    if (event.m_tracked_michels.empty()) continue;

    Variable* var = GetVar(variables, "tracked_tpi");
    assert(!event.m_tracked_michels.empty());
    assert(!event.m_tracked_michels.begin()->second.had_idx != -107);
    RecoPionIdx best_pi_idx = event.m_tracked_michels.begin()->second.had_idx;
    double fill_val = universe->GetTpi(best_pi_idx);

    {  // Fill selected
      // Sanity Checks
      assert(!(var->m_is_true && !event.m_is_mc));  // truth, but not MC?

      // total = signal & background, together
      if (event.m_is_mc) {
        var->m_hists.m_selection_mc.FillUniverse(*event.m_universe, fill_val,
                                                 event.m_weight);
      } else {
        var->m_hists.m_selection_data->Fill(fill_val);
      }

      // MC stuff
      // done with data. Mc-specific stuff
      if (event.m_is_mc) {
        // signal and background individually
        if (event.m_is_signal) {
          var->m_hists.m_effnum.FillUniverse(*event.m_universe, fill_val,
                                             event.m_weight);
        } else {
          var->m_hists.m_bg.FillUniverse(*event.m_universe, fill_val,
                                         event.m_weight);

          // Fill bg by W sideband category
          switch (event.m_w_type) {
            case kWSideband_Signal:
              break;
            case kWSideband_Low:
              var->m_hists.m_bg_loW.FillUniverse(*event.m_universe, fill_val,
                                                 event.m_weight);
              break;
            case kWSideband_Mid:
              var->m_hists.m_bg_midW.FillUniverse(*event.m_universe, fill_val,
                                                  event.m_weight);
              break;
            case kWSideband_High:
              var->m_hists.m_bg_hiW.FillUniverse(*event.m_universe, fill_val,
                                                 event.m_weight);
              break;
            default:
              std::cerr << "FillBackgrounds: no such W category\n";
              std::exit(2);
          }
        }
        FillStackedHists(event, var, fill_val);
      }
    }

    // for (auto v : variables) {
    //  std::cout << v->Name() << "  " << v->GetValue(event) << "\n";
    //}

    //// Fill reco -- enforce cuts, internally update the histograms owned by
    /// variables
    // ccpi_event::FillRecoEvent(event, variables);
  }  // entries
  std::cout << n_pass << " pass\n";
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

  mc_file_list = CCPi::GetPlaylistFile(plist, is_mc, true, true);
  is_mc = false;
  data_file_list = CCPi::GetPlaylistFile(plist, is_mc, true, true);

  //=========================================
  // Init macro utility
  //=========================================
  // const int signal_definition_int = SignalDefinition::OnePiTracked().m_id;
  const int signal_definition_int = SignalDefinition::OnePi().m_id;
  const bool is_grid = false;
  const bool do_truth = false;
  const bool do_systematics = false;

  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.m_name = "runRefactorTest";
  util.PrintMacroConfiguration();

  // VARIABLES
  CVUVariable* tracked_tpi =
      new CVUVariable("tracked_tpi", "T_{#pi}", "MeV", CCPi::GetBinning("tpi"));
  CVUVariable* untracked_tpi = new CVUVariable("untracked_tpi", "T_{#pi}",
                                               "MeV", CCPi::GetBinning("tpi"));
  CVUVariable* all_tpi =
      new CVUVariable("all_tpi", "T_{#pi}", "MeV", CCPi::GetBinning("tpi"));

  EventVariable* thetapi_deg = new EventVariable(
      "thetapi_deg", "#theta_{#pi}", "deg", CCPi::GetBinning("thetapi_deg"),
      &CCPiEvent::GetDummyVar);

  CVUVariable* pmu = new CVUVariable(
      "pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"), &CVUniverse::GetPmu);

  std::vector<Variable*> variables = {pmu, tracked_tpi, untracked_tpi, all_tpi};

  // std::vector<Variable*> variables = run_study_template::GetVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  // DATA
  LoopAndFill(util.m_data_universe, util.GetDataEntries(), kData,
              util.m_signal_definition, variables);

  std::cout << n_more_untracked << "  " << n_more_tracked << "  " << n_same
            << "\n";
  std::cout << n_untracked_multipi << "  " << n_tracked_multipi << "  "
            << n_multipi_agree << "\n";

  n_more_untracked = n_more_tracked = n_same = 0;
  n_untracked_multipi = n_tracked_multipi = n_multipi_agree = 0;

  // MC
  LoopAndFill(util.m_error_bands.at("cv").at(0), util.GetMCEntries(), kMC,
              util.m_signal_definition, variables);

  std::cout << n_more_untracked << "  " << n_more_tracked << "  " << n_same
            << "\n";
  std::cout << n_untracked_multipi << "  " << n_tracked_multipi << "  "
            << n_multipi_agree << "\n";

  // PLOTTING
  Variable* v = GetVar(variables, "tracked_tpi");

  // Stacks
  int ymax = -1;
  bool do_bwn = true;
  std::string tag("no w cut, no track cuts");
  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kS),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "SSB", ymax, do_bwn);

  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePion),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "NPions", ymax, do_bwn);

  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kOtherInt),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "FSPart", ymax, do_bwn);

  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kRES),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "Channel", ymax, do_bwn);

  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kPim),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "Hadron", ymax, do_bwn);

  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePi0),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "NPi0", ymax, do_bwn);

  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePip),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "NPip", ymax, do_bwn);

  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kLowW),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "Wtrue", ymax, do_bwn);

  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kCOHERENT_S),
             util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag,
             "Coherent", ymax, do_bwn);

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

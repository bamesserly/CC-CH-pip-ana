//==============================================================================
// This script visualizes the sideband region and the quality of the sideband
// fit. It creates stacked plots of analysis variables from sideband region
// events, both before and after the fit.
// Plots are broken down into the truth categories that are also used to
// perform the fit (regions of Wexptrue).
//==============================================================================
#ifndef runSidebands_C
#define runSidebands_C

#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/HadronVariable.h"
#include "includes/MacroUtil.h"
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/Variable.h"
#include "includes/common_functions.h"      // GetVar
#include "xsec/crossSectionDataFromFile.C"  // DoWSidebandTune
#include "xsec/plotting_functions.h"

class Variable;
class HadronVariable;

//==============================================================================
// Get Variables
//==============================================================================
namespace run_sidebands {
typedef Variable Var;
typedef HadronVariable HVar;

std::vector<Variable*> GetOnePiVariables(bool include_truth_vars = false) {
  HVar* tpi = new HVar("tpi", "T_{#pi}", "MeV", CCPi::GetBinning("tpi"),
                       &CVUniverse::GetTpi);

  HVar* tpi_mbr = new HVar("tpi_mbr", "T_{#pi} (MBR)", "MeV",
                           CCPi::GetBinning("tpi"), &CVUniverse::GetTpiMBR);

  HVar* thetapi_deg =
      new HVar("thetapi_deg", "#theta_{#pi}", "deg",
               CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);

  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);

  Var* thetamu_deg =
      new Var("thetamu_deg", "#theta_{#mu}", "deg",
              CCPi::GetBinning("thetamu_deg"), &CVUniverse::GetThetamuDeg);

  Var* enu = new Var("enu", "E_{#nu}", "MeV", CCPi::GetBinning("enu"),
                     &CVUniverse::GetEnu);

  Var* q2 = new Var("q2", "Q^{2}", "#frac{MeV^{2}}{c^{2}})",
                    CCPi::GetBinning("q2"), &CVUniverse::GetQ2);

  Var* wexp = new Var("wexp", "W_{exp}", "MeV", CCPi::GetBinning("wexp"),
                      &CVUniverse::GetWexp);

  Var* wexp_fit = new Var(sidebands::kFitVarString, wexp->m_hists.m_xlabel,
                          wexp->m_units, 32, 0.e3, 3.2e3, &CVUniverse::GetWexp);

  std::vector<Var*> variables = {tpi,
                                 // tpi_mbr,
                                 thetapi_deg, pmu, thetamu_deg, enu, q2, wexp,
                                 wexp_fit};

  return variables;
}

std::map<std::string, Variable*> GetOnePiVariables_Map(
    bool include_truth_vars = true) {
  std::map<std::string, Var*> var_map;
  std::vector<Var*> var_vec = GetOnePiVariables(include_truth_vars);
  for (auto v : var_vec) var_map[v->m_label] = v;
  return var_map;
}
}  // namespace run_sidebands

std::vector<Variable*> GetSidebandVariables(SignalDefinition signal_definition,
                                            bool include_truth_vars = false) {
  std::vector<Variable*> variables;
  switch (signal_definition) {
    case kOnePi:
      variables = run_sidebands::GetOnePiVariables(include_truth_vars);
      break;
    default:
      std::cerr << "Variables for other SDs not yet implemented.\n";
      std::exit(1);
  }
  return variables;
}

//==============================================================================
// Loop
//==============================================================================
void FillWSideband(const CCPi::MacroUtil& util, CVUniverse* universe,
                   const EDataMCTruth& type,
                   std::vector<Variable*>& variables) {
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    // First -- Fill histos
    // With events in the sideband region (passes all cuts except, W > 1.5),
    // fill histograms for all variables including the tune variable (Wexp
    // reco), broken down by the (3) truth categories (
    //   (A)       Wexptrue < 1.4,
    //   (B) 1.4 < Wexptrue < 1.8,
    //   (C) 1.8 < Wexptrue
    // )
    // that we'll use to perform the fit.
    // For variables other than the fit var, if the data-mc agreement is bad,
    // it undermines the sideband tune.
    //
    // Second -- visualize/study the sideband region
    // With events that pass all cuts except for a W cut, fill W so you can
    // visualize the whole spectrum.
    //
    // PassesCuts returns is_w_sideband in the process of checking all cuts.
    bool passes_all_cuts =
        PassesCuts(*universe, event.m_reco_pion_candidate_idxs, is_mc,
                   util.m_signal_definition, event.m_is_w_sideband);

    if (event.m_is_w_sideband) ccpi_event::FillWSideband(event, variables);

    ccpi_event::FillWSideband_Study(event, variables);
  }  // end event loop

  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Sync
//==============================================================================
namespace run_sidebands {
void SyncAllHists(Variable& var) {
  var.m_hists.m_wsidebandfit_sig.SyncCVHistos();
  var.m_hists.m_wsidebandfit_loW.SyncCVHistos();
  var.m_hists.m_wsidebandfit_midW.SyncCVHistos();
  var.m_hists.m_wsidebandfit_hiW.SyncCVHistos();
}
}  // namespace run_sidebands

//==============================================================================
// Main
//==============================================================================
void runSidebands(int signal_definition_int = 0, const char* plist = "ME1A",
                  int do_systematics = 0) {
  // INIT MACRO UTILITY OBJECT
  std::string mc_file_list = GetPlaylistFile(plist, true /*is mc*/);
  std::string data_file_list = GetPlaylistFile(plist, false);

  const std::string macro("runCutVariables");
  bool do_truth = false, is_grid = false;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // INIT VARS, HISTOS, AND EVENT COUNTERS
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetSidebandVariables(util.m_signal_definition, do_truth_vars);

  for (auto var : variables) {
    var->InitializeSidebandHists(util.m_error_bands);
    var->InitializeStackedHists();
    var->InitializeDataHists();
  }

  FillWSideband(util, util.m_data_universe, kData, variables);
  FillWSideband(util, util.m_error_bands.at("cv").at(0), kMC, variables);

  for (auto var : variables) run_sidebands::SyncAllHists(*var);

  //============================================================================
  // End event selection
  // Begin sideband tune & BG Sub
  //============================================================================
  PlotUtils::HistWrapper<CVUniverse> hw_loW_fit_wgt =
      PlotUtils::HistWrapper<CVUniverse>("h_loW_fit_wgt",
                                         "W Sideband Fit Weight -- low W", 1,
                                         0., 15., util.m_error_bands);

  PlotUtils::HistWrapper<CVUniverse> hw_midW_fit_wgt =
      PlotUtils::HistWrapper<CVUniverse>("h_midW_fit_wgt",
                                         "W Sideband Fit Weight -- mid W", 1,
                                         0., 15., util.m_error_bands);

  PlotUtils::HistWrapper<CVUniverse> hw_hiW_fit_wgt =
      PlotUtils::HistWrapper<CVUniverse>("h_hiW_fit_wgt",
                                         "W Sideband Fit Weight -- high W", 1,
                                         0., 15., util.m_error_bands);

  // Sideband tune
  // Fill the fit parameter hists (by reference)
  DoWSidebandTune(util, GetVar(variables, sidebands::kFitVarString),
                  hw_loW_fit_wgt, hw_midW_fit_wgt, hw_hiW_fit_wgt);

  //============================================================================
  // Plot
  //============================================================================
  std::string tag;
  double ymax = -1;

  // Plot W before fit, with no W cut
  if (1) {
    Variable* var = GetVar(variables, sidebands::kFitVarString);
    PlotWSidebandStacked(var, var->m_hists.m_wsideband_data,
                         var->GetStackArray(static_cast<WSidebandType>(0)),
                         util.m_data_pot, util.m_mc_pot,
                         util.m_signal_definition, tag, ymax);
  }

  // Plot all vars W before and after fit
  if (1) {
    for (auto var : variables) {
      tag = "SidebandRegion";
      bool do_prefit = true;
      bool do_bin_width_norm = true;
      CVUniverse* universe = util.m_error_bands.at("cv").at(0);
      PlotFittedW(var, *universe, hw_loW_fit_wgt, hw_midW_fit_wgt,
                  hw_hiW_fit_wgt, util.m_data_pot, util.m_mc_pot,
                  util.m_signal_definition, do_prefit, tag, ymax,
                  do_bin_width_norm);
      do_prefit = false;
      PlotFittedW(var, *universe, hw_loW_fit_wgt, hw_midW_fit_wgt,
                  hw_hiW_fit_wgt, util.m_data_pot, util.m_mc_pot,
                  util.m_signal_definition, do_prefit, tag, ymax,
                  do_bin_width_norm);
    }
  }
}

#endif  // runSidebands_C

//==============================================================================
// This script visualizes the sideband region and the quality of the sideband
// fit.
// It creates stacked plots of analysis variables from sideband region events,
// both before and after the fit.
//
// Plots are broken down into the truth categories that are also used to
// perform the fit (regions of Wexptrue).
//
// Does not perform the fit in systematic universes.
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
#include "xsec/makeCrossSectionMCInputs.C" // GetAnalysisVariables
#include "PlotUtils/LowRecoilPionReco.h"
#include "PlotUtils/LowRecoilPionCuts.h"

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

  Var* wexp_fit =
      new Var(sidebands::kFitVarString, wexp->m_hists.m_xlabel, wexp->m_units,
              CCPi::GetBinning("wexp_fit"), &CVUniverse::GetWexp);

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
void FillWSideband(const CCPi::MacroUtil& util, const EDataMCTruth& type,
                   std::vector<Variable*>& variables) {
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  // typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;
  UniverseMap error_bands;
  switch (type) {
    case kData:
      error_bands.insert(
          std::make_pair(util.m_data_universe->ShortName(),
                         std::vector<CVUniverse*>{util.m_data_universe}));
      break;
    case kMC:
      error_bands = util.m_error_bands;
      break;
    case kTruth:
      error_bands = util.m_error_bands_truth;
      break;
    default:
      std::cerr << "invalid tuple type\n";
  }
  // util.m_error_bands.at("cv").at(0)

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    if (i_event == 10000) break;
    for (auto error_band : error_bands) {
      std::vector<CVUniverse*> universes = error_band.second;
      for (auto universe : universes) {
        universe->SetEntry(i_event);
        CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
        universe->SetTruth(is_mc);
        LowRecoilPion::Cluster d;
        LowRecoilPion::Cluster c(*universe,0);
        LowRecoilPion::Michel<CVUniverse> m(*universe,0);
        LowRecoilPion::MichelEvent<CVUniverse> trackless_michels;
        bool good_trackless_michels = LowRecoilPion::hasMichel<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::hasMichelCut(*universe, trackless_michels);
        // good_trackless_michels = BestMichelDistance2DCut(*universe, trackless_michels);
        good_trackless_michels = good_trackless_michels && LowRecoilPion::BestMichelDistance2D<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::BestMichelDistance2DCut(*universe, trackless_michels);
        // good_trackless_michels = MichelRangeCut(*universe, trackless_michels);
        good_trackless_michels = good_trackless_michels && LowRecoilPion::GetClosestMichel<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::GetClosestMichelCut(*universe, trackless_michels); 
        // First -- Fill histos
        // With events in the sideband region (passes all cuts except, W > 1.5),
        // fill histograms for all variables including the tune variable (Wexp
        // reco), broken down by the (3) truth categories (
        //   (A)       Wexptrue < 1.4,
        //   (B) 1.4 < Wexptrue < 1.8,
        //   (C) 1.8 < Wexptrue
        // )
        // that we'll use to perform the fit.
        // For variables other than the fit var, if the data-mc agreement is
        // bad, it undermines the sideband tune.
        //
        // Second -- visualize/study the sideband region
        // With events that pass all cuts except for a W cut, fill W so you can
        // visualize the whole spectrum.
        //
        
        universe->SetVtxMichels(trackless_michels);
        bool pass = true; 
        pass = pass && universe->GetNMichels() == 1;
        pass = pass && universe->GetTpiTrackless() < 350.;
        pass = pass && universe->GetPmu() > 1500.;
        pass = pass && universe->GetPmu() < 20000.;
        pass = pass && universe->GetNIsoProngs() < 2; 
        pass = pass && universe->IsInHexagon(universe->GetVecElem("vtx", 0), universe->GetVecElem("vtx", 1), 850.);
        pass = pass && universe->GetVecElem("vtx", 2) > 5990.;
        pass = pass && universe->GetVecElem("vtx", 2) < 8340.;
        pass = pass && universe->GetBool("isMinosMatchTrack");  
        pass = pass && universe->GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
        pass = pass && universe->GetThetamuDeg() < 20;
        // PassesCuts returns is_w_sideband in the process of checking all cuts.
        PassesCutsInfo cuts_info = PassesCuts(event);
        std::tie(event.m_passes_cuts, event.m_is_w_sideband,
                 event.m_passes_all_cuts_except_w,
                 event.m_reco_pion_candidate_idxs) = cuts_info.GetAll();
        event.m_highest_energy_pion_idx =
            GetHighestEnergyPionCandidateIndex(event);
        universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);
        universe->SetVtxMichels(trackless_michels);
        if (is_mc) event.m_weight = universe->GetWeight();
        event.m_passes_trackless_cuts_except_w = pass;
        event.m_passes_trackless_sideband = false;
        if (pass && universe->GetTracklessWexp() > 1400.){
          if (universe->GetTracklessWexp() > 1500.){
            event.m_passes_trackless_sideband = true;}
          pass = false;
        }
        event.m_passes_trackless_cuts = good_trackless_michels && pass;
        event.m_passes_trackless_sideband = event.m_passes_trackless_sideband && good_trackless_michels;
        event.m_passes_trackless_cuts_except_w = event.m_passes_trackless_cuts_except_w && good_trackless_michels;
        universe->SetPassesTrakedTracklessCuts(event.m_passes_cuts,
                   event.m_passes_trackless_cuts, event.m_is_w_sideband,
                   event.m_passes_trackless_sideband, event.m_passes_all_cuts_except_w,
                   event.m_passes_trackless_cuts_except_w);


        // Fill histograms of all variables with events in the sideband region.
        // Each variable has 4 such histograms for signal, low-w sb, med-w sb,
        // and high-w sb
        //
        // Fit will be performed in the wexp_fit variable. Other variables
        // will be filled in order to visualize the impact of the fit via the
        // PlotFittedW function.
        //
        // If do_systematics, then fill sideband events for all variables in
        // all universes. Filling the for variables other than wexp_fit in any
        // universe other than the CV is not used.
        //
        // Technically, we're filling m_wsidebandfit_[sig, lo, mid, hi, data]
        // for all variables.
        if (event.m_is_w_sideband || event.m_passes_trackless_sideband)
          ccpi_event::FillWSideband(event, variables);

        // Fill events in and out of the sideband for the wexp_fit variable.
        //
        // To be visualized by PlotWSidebandStacked.
        //
        // Filling wexp_fit->GetStackComponentHist(event.m_w_type) and
        // wexp_fit->m_hists.m_wsideband_data
        if ((event.m_passes_all_cuts_except_w ||
            event.m_passes_trackless_cuts_except_w) &&
            event.m_universe->ShortName() == "cv")
          ccpi_event::FillWSideband_Study(event, variables);
      }  // universes
    }    // error bands
  }      // end event loop

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
void runSidebands(int signal_definition_int = 0, const char* plist = "ME1L",
                  int do_systematics = 0) {
  // INIT MACRO UTILITY OBJECT
  bool use_xrootd = false;
  std::string mc_file_list = GetPlaylistFile(plist, true /*is mc*/, use_xrootd);
  std::string data_file_list = GetPlaylistFile(plist, false, use_xrootd);
  //std::string data_file_list = GetTestPlaylist(false);
  //std::string mc_file_list = GetTestPlaylist(true);

  const std::string macro("runSidebands");
  bool do_truth = false, is_grid = false;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.m_pot_scale = util.m_data_pot / util.m_mc_pot;
  // util.m_pot_scale = 1.;
  util.PrintMacroConfiguration(macro);

  // INIT VARS, HISTOS, AND EVENT COUNTERS
  const bool do_truth_vars = false;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars); 
//      GetSidebandVariables(util.m_signal_definition, do_truth_vars);

  for (auto var : variables) {
    var->InitializeSidebandHists(util.m_error_bands);
    var->InitializeStackedHists();
    var->InitializeDataHists();
  }

  FillWSideband(util, kData, variables);
  FillWSideband(util, kMC, variables);
  // util.m_error_bands.at("cv").at(0)

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
  // This plot is filled by FillWSideband_study
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

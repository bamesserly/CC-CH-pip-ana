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
#include "includes/common_functions.h"
#include "PlotUtils/LowRecoilPionReco.h"
#include "PlotUtils/LowRecoilPionCuts.h"
#include "plotting_functions.h"

// Forward declare my variables because we're hiding the header.
class Variable;
class HadronVariable;

namespace run_study_template {
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
  if(event.m_passes_trackless_cuts)// && !event.m_is_signal)
    ccpi_event::FillStackedHists(event, variables);
}

//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  typedef Variable Var;
  typedef HadronVariable HVar;
//HVar* thetapi_deg = new HVar("thetapi_deg", "#theta_{#pi}", "deg", CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);
  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"), &CVUniverse::GetPmu);
  Var* tpi_trackless = new Var("mtpi", "Mehreen T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
                      &CVUniverse::GetTpiTrackless);
  std::vector<Var*> variables = {tpi_trackless, pmu};
  return variables;
}
} // namespace run_study_template

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFill(const CCPi::MacroUtil& util, CVUniverse* universe,
                        const EDataMCTruth& type,
                        std::vector<Variable*>& variables,
                        double& signal, double& bg) {

  std::cout << "Loop and Fill CutVars\n";
  bool is_mc, is_truth = false;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  for(Long64_t i_event=0; i_event < n_entries; ++i_event){
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    universe->SetTruth(is_truth);
//    if (i_event == 1000) break;   

    LowRecoilPion::Cluster d;
    LowRecoilPion::Cluster c(*universe,0);
    LowRecoilPion::Michel<CVUniverse> m(*universe,0);
    LowRecoilPion::MichelEvent<CVUniverse> trackless_michels;

    // MEHREEN CUTS -- all of these functions fill/modify the trackless_michels.

    // basically we have this:
    // MichelEvent GetQualityTracklessMichels(*universe) {
    //   MichelEvent trackless_michels;
    //   bool pass = HasMichelCut(*universe, trackless_michels);
    //   pass = BestMichelDistance2DCut(*universe, trackless_michels);
    //   pass = MichelRangeCut(*universe, trackless_michels);
    //     return {trackless_michels, pass}
    //   }

    // code location: MAT-MINERvA/utilities/LowRecoilPionReco.h
    //                MAT-MINERvA/calculators/LowRecoilPionFunctions.h
    //                MAT-MINERvA/calculators/LowRecoilPionCuts.h

    //bool pass = HasMichelCut(*universe, trackless_michels);
    bool good_trackless_michels = LowRecoilPion::hasMichel<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::HasMichelCut(*universe, trackless_michels);

    // good_trackless_michels = BestMichelDistance2DCut(*universe, trackless_michels);
    good_trackless_michels = good_trackless_michels && LowRecoilPion::BestMichelDistance2D<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::BestMichelDistance2DCut(*universe, trackless_michels);

    // good_trackless_michels = MichelRangeCut(*universe, trackless_michels);
    good_trackless_michels = good_trackless_michels && LowRecoilPion::GetClosestMichel<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::MichelRangeCut(*universe, trackless_michels);

    universe->SetVtxMichels(trackless_michels);
 
/*    if (good_trackless_michels){
      std::cout << "trackless_michels.m_bestdist = " << trackless_michels.m_bestdist << "\n";
      std::cout << "universe->m_vtx_michels.m_bestdist = " << universe->GetBestDistance() << "\n";
      std::cout << "universe->GetTpiUntracked(trackless_michels.m_bestdist) = " << universe->GetTpiUntracked(trackless_michels.m_bestdist) << "\n";
      std::cout << "universe->GetTpiTrackless() = " << universe->GetTpiTrackless() << "\n";
    }*/
    // For mc, get weight, check signal, and sideband

    // This cuts are used to have a similar topology to the 1 pi analysis
    bool pass = true;
    pass = pass && universe->GetNMichels() == 1;
    pass = pass && universe->GetTpiTrackless() < 350.;
    pass = pass && universe->GetWexp() < 1400.;
    pass = pass && universe->GetPmu() > 1500.;
    pass = pass && universe->GetPmu() < 20000.;
    pass = pass && universe->GetNIsoProngs() < 2; 
    pass = pass && universe->IsInHexagon(universe->GetVecElem("vtx", 0), universe->GetVecElem("vtx", 1), 850.);
    pass = pass && universe->GetVecElem("vtx", 2) > 5990.;
    pass = pass && universe->GetVecElem("vtx", 2) < 8340.;
    pass = pass && universe->GetBool("isMinosMatchTrack");  
    pass = pass && universe->GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
    pass = pass && universe->GetThetamuDeg() < 20;

    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);


    // if event passes ccpi cuts
    //   if good_trackless_michels // "mix"
    //     plot trackless_michels.m_best_dist
    event.m_passes_trackless_cuts = good_trackless_michels && pass;
    if (event.m_passes_trackless_cuts && is_mc){
//      std::cout << "Pass Cuts \n";
      if (event.m_is_signal){
        signal = signal + 1.;
//        std::cout << "Is signal \n";
      }  
      else{
        bg = bg + 1.; 
//        std::cout << "Is bg \n";
      }
    }
    // WRITE THE FILL FUNCTION
    run_study_template::FillVars(event, variables);
  } // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runMATMichels(std::string plist = "ME1A") {
  //=========================================
  // Input tuples
  //=========================================
    bool is_mc = true;
    std::string mc_file_list, data_file_list;
    mc_file_list = GetPlaylistFile(plist, is_mc);
    is_mc = false;
    data_file_list = GetPlaylistFile(plist, is_mc);
    
  //=========================================
  // Init macro utility
  //=========================================
    const int signal_definition_int = 0;
    const std::string macro("runStudyTemplate");
    const bool is_grid = false;
    const bool do_truth = false;
    const bool do_systematics = false;

    CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list, plist, do_truth, is_grid, do_systematics);
    util.PrintMacroConfiguration(macro);

  //=========================================
  // Get variables and initialize their hists
  //=========================================
  std::vector<Variable*> variables = run_study_template::GetVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  //=========================================
  // Loop and Fill
  //=========================================
  double signal = 1., bg = 1.;
  LoopAndFill(util, util.m_data_universe,              kData, variables, signal, bg);
  signal = 0.;
  bg = 0.;
  LoopAndFill(util, util.m_error_bands.at("cv").at(0), kMC,   variables, signal, bg);

  std::cout << "Signal = " << signal << "\n"
            << "Background = " << bg << "\n"
            << "Purity = " << signal/(signal+bg) << "\n";

  for (auto v : variables) {
    std::string tag = v->Name();
    double ymax = -1;
    bool do_bwn = true;
    bool draw_arrow = v->Name() == "Wexp" ? true : false;
    std::cout << "Plotting" << std::endl;
    std::string study = "Sel_1piCuts";
    double data_pot = util.m_data_pot; 
    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kOtherInt),
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "FSP", ymax, draw_arrow, study);

    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kCCQE),    
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Int",  ymax, draw_arrow, study);

    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kPim),     
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Hadrons", ymax, draw_arrow, study);

    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePion), 
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npi",  ymax, draw_arrow, study);

    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePi0),  
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npi0", ymax, draw_arrow, study);

    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePip),  
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npip", ymax, draw_arrow, study);

    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kWSideband_Low),  
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "WSB", ymax, draw_arrow, study);

    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kB_Meson), 
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Msn",  ymax, draw_arrow, study);

    PlotBreakdown(v, v->m_hists.m_selection_data, v->GetStackArray(kB_HighW), 
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "WBG", ymax, draw_arrow, study);
  }

  std::cout << "Success" << std::endl;
}

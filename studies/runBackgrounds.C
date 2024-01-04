//==============================================================================
// This script visualizes diagnostic variables and the variables that we cut on
// before and after the cut is performed.
// Plots are broken down into their truth categories (includes/TruthCategories).
//==============================================================================
#ifndef runBackgrounds_C
#define runBackgrounds_C

#include "includes/MacroUtil.h"
#include "includes/CCPiEvent.h"
#include "includes/Variable.h"
#include "includes/HadronVariable.h"
#include "includes/TruthMatching.h" //GetTruthCategory functions
#include "plotting_functions.h"
#include "xsec/makeCrossSectionMCInputs.C" // GetAnalysisVariables
#include "PlotUtils/LowRecoilPionReco.h"
#include "PlotUtils/LowRecoilPionCuts.h"

class Variable;
class HadronVariable;

//==============================================================================
// Loop
//==============================================================================
void LoopAndFillBackgrounds(const CCPi::MacroUtil& util, CVUniverse* universe,
                     std::vector<Variable*>& variables) {
  bool is_mc = true;
  bool is_truth = false;

  std::cout << " *** Looping MC to Fill Backgrounds ***\n";
  for(Long64_t i_event=0; i_event < util.GetMCEntries(); ++i_event) {
  //for(Long64_t i_event=0; i_event < 5000; ++i_event) {
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
//    if (i_event == 10000)break;
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    bool is_w_sideband = false;
    LowRecoilPion::Cluster d;
    LowRecoilPion::Cluster c(*universe,0);
    LowRecoilPion::Michel<CVUniverse> m(*universe,0);
    LowRecoilPion::MichelEvent<CVUniverse> trackless_michels;
    bool good_trackless_michels = LowRecoilPion::hasMichel<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::hasMichelCut(*universe, trackless_michels);
    good_trackless_michels = good_trackless_michels && LowRecoilPion::BestMichelDistance2D<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::BestMichelDistance2DCut(*universe, trackless_michels);
    good_trackless_michels = good_trackless_michels && LowRecoilPion::GetClosestMichel<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::GetClosestMichelCut(*universe, trackless_michels);
    // Get Quality Michels

    universe->SetVtxMichels(trackless_michels);
    bool pass = true; 
    pass = pass && universe->GetNMichels() == 1;
    pass = pass && universe->GetTpiTrackless() < 350.;
    pass = pass && universe->GetTpiTrackless() > 0.;
    pass = pass && universe->GetPmu() > 1500.;
    pass = pass && universe->GetPmu() < 20000.;
    pass = pass && universe->GetNIsoProngs() < 2; 
    pass = pass && universe->IsInHexagon(universe->GetVecElem("vtx", 0), universe->GetVecElem("vtx", 1), 850.);
    pass = pass && universe->GetVecElem("vtx", 2) > 5990.;
    pass = pass && universe->GetVecElem("vtx", 2) < 8340.;
    pass = pass && universe->GetInt("isMinosMatchTrack")==1;  
    pass = pass && universe->GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
    pass = pass && universe->GetThetamu() < CCNuPionIncConsts::kThetamuMaxCutVal;
    pass = pass && universe->GetTracklessWexp() > 0.;

    PassesCutsInfo cuts_info = PassesCuts(event);
    std::tie(event.m_passes_cuts, event.m_is_w_sideband, event.m_passes_all_cuts_except_w, event.m_reco_pion_candidate_idxs) = cuts_info.GetAll();
    event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);

    universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);
    universe->SetVtxMichels(trackless_michels);
    event.m_passes_trackless_cuts_except_w = pass;
    event.m_passes_trackless_sideband = false;
    event.m_weight = universe->GetWeight();
    if (pass && universe->GetTracklessWexp() > 1400){
      if (universe->GetTracklessWexp() >= sidebands::kSidebandCutVal) event.m_passes_trackless_sideband = true;
      pass = false;
    }
    event.m_passes_trackless_cuts = good_trackless_michels && pass;
    universe->SetPassesTrakedTracklessCuts(event.m_passes_cuts, event.m_passes_trackless_cuts, event.m_is_w_sideband, event.m_passes_trackless_sideband, event.m_passes_all_cuts_except_w, event.m_passes_trackless_cuts_except_w);
//    std::cout << "Event = " << i_event << "\n"; 
//    std::cout << "Pass Tracked cuts" << event.m_passes_cuts << "\n";
//    std::cout << "Pass Trackless cuts" << event.m_passes_trackless_cuts << "\n";
    if ((event.m_passes_cuts || event.m_passes_trackless_cuts) && !event.m_is_signal) {
      ccpi_event::FillStackedHists(event, variables);
    }
  } // events
  std::cout << "*** Done ***\n\n";
}


// Plot
void PlotAllBackgrounds(Variable* v, const CCPi::MacroUtil& util) {
  std::string tag;
  double ymax = -1;
  bool draw_arrow = v->Name() == "Wexp" ? true : false;

  // Not plotting data, so do this to make the scaling factor = 1.
  // There are other ways to solve this problem, but this is the easiest.
  double data_pot = util.m_mc_pot;

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kOtherInt),
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "FSP", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kCCQE),    
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Int",  ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kPim),     
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Hadrons", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePion), 
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npi",  ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePi0),  
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npi0", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePip),  
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npip", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kWSideband_Low),  
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "WSB", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kB_Meson), 
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Msn",  ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kB_HighW), 
                 data_pot, util.m_mc_pot, util.m_signal_definition,
                 "WBG", ymax, draw_arrow);
}


//==============================================================================
// Main
//==============================================================================
void runBackgrounds(int signal_definition_int = 0, 
                     const char* plist = "ALL", bool is_grid = false,
 		      std::string input_file = "", int run = 0) 
 {

  // INPUT TUPLES
//  std::string input_file = "";
//  bool is_grid = false;
  const bool is_mc = true;
  std::string mc_file_list;
  assert(!(is_grid && input_file.empty()) &&
         "On the grid, infile must be specified.");
  // const bool use_xrootd = false;
  mc_file_list = input_file.empty()
                     ? GetPlaylistFile(plist, is_mc /*, use_xrootd*/)
                     : input_file;
  TFile fout(Form("Background_Breakdown_%s_%d.root", plist, run), "RECREATE"); 

  // Init macro utility object
  const std::string macro("runBackgrounds");
  bool do_data = false;
  bool do_mc = true;
  bool do_truth = false;
  bool do_systematics = false;
  bool do_grid = false;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, plist, do_truth,
                       is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // Init vars, histos, and event counters
  const bool do_truth_vars = false;

  std::vector<Variable*> variables = 
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  for (auto var : variables) {
    var->InitializeStackedHists(); 
    var->InitializeDataHists(); 
  }


  // Fill
  CVUniverse* cvu = util.m_error_bands.at("cv").at(0);
  LoopAndFillBackgrounds(util, cvu, variables);


  // Plot
  fout.cd();
  for (auto v : variables){
    v->GetStackArray(kOtherInt).Write(Form("%s_FSP", v->Name().c_str()));
    v->GetStackArray(kCCQE).Write(Form("%s_Int", v->Name().c_str()));
    v->GetStackArray(kPim).Write(Form("%s_Hadrons", v->Name().c_str()));
    v->GetStackArray(kOnePion).Write(Form("%s_Npi", v->Name().c_str()));
    v->GetStackArray(kOnePi0).Write(Form("%s_Npi0", v->Name().c_str()));
    v->GetStackArray(kOnePip).Write(Form("%s_Npip", v->Name().c_str()));
    v->GetStackArray(kWSideband_Low).Write(Form("%s_WSB", v->Name().c_str()));
    v->GetStackArray(kB_Meson).Write(Form("%s_Msn", v->Name().c_str()));
    v->GetStackArray(kB_HighW).Write(Form("%s_WBG", v->Name().c_str()));
    if (!is_grid) PlotAllBackgrounds(v, util);
  }

}

#endif // runBackgrounds_C

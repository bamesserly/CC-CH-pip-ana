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
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    bool is_w_sideband = false;
    event.m_passes_cuts = PassesCuts(event, is_w_sideband);
    event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
    if (event.m_passes_cuts && !event.m_is_signal) {
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
                     const char* plist = "ME1L") {

  // INPUT TUPLES
  std::string input_file = "";
  bool is_grid = false;
  const bool is_mc = true;
  std::string mc_file_list;
  assert(!(is_grid && input_file.empty()) &&
         "On the grid, infile must be specified.");
  // const bool use_xrootd = false;
  mc_file_list = input_file.empty()
                     ? GetPlaylistFile(plist, is_mc /*, use_xrootd*/)
                     : input_file;


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
  for (auto v : variables)
    PlotAllBackgrounds(v, util);
}

#endif // runBackgrounds_C

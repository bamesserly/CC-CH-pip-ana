//==============================================================================
// This script visualizes diagnostic variables and the variables that we cut on
// before and after the cut is performed.
// Plots are broken down into their truth categories (includes/TruthCategories).
//==============================================================================
#ifndef runBackgrounds_C
#define runBackgrounds_C

#include "includes/CCPiMacroUtil.h"
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
void LoopAndFillBackgrounds(const CCPiMacroUtil& util, CVUniverse* universe,
                     std::vector<Variable*>& variables) {
  bool is_mc = true;
  bool is_truth = false;

  std::cout << " *** Looping MC to Fill Backgrounds ***\n";
  for(Long64_t i_event=0; i_event < util.GetMCEntries(); ++i_event) {
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    event.m_passes_cuts = PassesCuts(event);
    if (event.m_passes_cuts && !event.m_is_signal)
      ccpi_event::FillStackedHists(event, variables);
  } // events
  std::cout << "*** Done ***\n\n";
}


// Plot
void PlotAllBackgrounds(Variable* v, const CCPiMacroUtil& util) {
  std::string tag;
  double ymax = -1;
  bool draw_arrow = v->Name() == "Wexp" ? true : false;

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kOtherInt),
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "FSP", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kCCQE),    
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Int",  ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kPim),     
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Hadrons", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePion), 
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npi",  ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePi0),  
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npi0", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kOnePip),  
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Npip", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kWSideband_Low),  
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "WSB", ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kB_Meson), 
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "Msn",  ymax, draw_arrow);

  PlotBackground(v, v->m_hists.m_selection_data, v->GetStackArray(kB_HighW), 
                 util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                 "WBG", ymax, draw_arrow);
}


//==============================================================================
// Main
//==============================================================================
void runBackgrounds(int signal_definition_int = 0, 
                     const char* plist = "ME1B") {
  // Init macro utility object
  const std::string macro("runBackgrounds");
  bool do_data = false, do_mc = true, do_truth = false;
  bool do_systematics = false, do_grid = false;
  CCPiMacroUtil util(signal_definition_int, plist, do_data, do_mc, do_truth,
                     do_systematics, do_grid);
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

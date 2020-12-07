#ifndef runCutVariables_C
#define runCutVariables_C

#include "includes/MacroUtil.h"
#include "includes/CCPiEvent.h"
#include "includes/Variable.h"
#include "includes/HadronVariable.h"
#include "includes/common_functions.h" // GetVar
#include "includes/Constants.h" // EDataMCTruth
#include "plotting_functions.h"
#include "xsec/makeCrossSectionMCInputs.C" // GetAnalysisVariables

//==============================================================================
// This script visualizes diagnostic variables and the variables that we cut on
// before and after the cut is performed.
// Plots are broken down into their truth categories (includes/TruthCategories).
//==============================================================================

class Variable;
class HadronVariable;

//==============================================================================
// Get Variables
//==============================================================================
namespace run_cut_variables {
  typedef Variable Var;
  typedef HadronVariable HVar;

  std::vector<Variable*> GetOnePiVariables(bool include_truth_vars = false) {
    Var* n_had_tracks = new Var("n_had_tracks", "n_had_tracks", "",    5,  0., 5.);
    Var* wexp_cut     = new Var("wexp_cut",     "W_{exp}",  "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* michel_count = new Var("michel_count", "michel_count", "",    5,  0., 5.);

    HVar* llr     = new HVar("llr",     "LLR",     "", 40,  -60., 60., &CVUniverse::GetLLRScore);
    HVar* enode01 = new HVar("enode01", "Enode01", "", 70,  0.,   70., &CVUniverse::GetEnode01);
    HVar* enode2  = new HVar("enode2",  "Enode2",  "", 40, 0.,   40., &CVUniverse::GetEnode2);
    HVar* enode3  = new HVar("enode3",  "Enode3",  "", 40, 0.,   40., &CVUniverse::GetEnode3);
    HVar* enode4  = new HVar("enode4",  "Enode4",  "", 40, 0.,   40., &CVUniverse::GetEnode4);
    HVar* enode5  = new HVar("enode5",  "Enode5",  "", 40, 0.,   40., &CVUniverse::GetEnode5);

    Var* n_iso_prongs = new Var("n_iso_prongs", "n_iso_prongs", "",  5, 0., 5.,    &CVUniverse::GetNIsoProngs);
    Var* n_pions      = new Var("n_pions",      "NPions",       "",  5, 0,  5.);

    Var* n_short_tracks = new Var("n_short_tracks", "n_short_tracks", "", 5, 0., 5., &CVUniverse::GetNAnchoredShortTracks);
    Var* n_long_tracks = new Var("n_long_tracks", "n_long_tracks", "", 5, 0., 5., &CVUniverse::GetNAnchoredLongTracks);

    Var* fit_vtx_x = new Var("fit_vtx_x", "fit_vtx_x", "cm", 50, -100, 100, &CVUniverse::GetFitVtxX);
    Var* fit_vtx_y = new Var("fit_vtx_y", "fit_vtx_y", "cm", 50, -100, 100, &CVUniverse::GetFitVtxY);
    Var* fit_vtx_z = new Var("fit_vtx_z", "fit_vtx_z", "cm", 50, -1000, 1000, &CVUniverse::GetFitVtxZ);
    
    HVar* track_reco_meth = new HVar("track_reco_meth", "trac_reco_meth", "", 8, 0., 8., &CVUniverse::GetTrackReconstructionMethod);
    HVar* n_nodes         = new HVar("n_nodes", "n_nodes", "", 50, 0., 50., &CVUniverse::GetNNodes);

    Var* wexp0        = new Var("wexp0", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp1        = new Var("wexp1", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp2        = new Var("wexp2", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp3        = new Var("wexp3", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp4        = new Var("wexp4", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp5        = new Var("wexp5", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp6        = new Var("wexp6", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp7        = new Var("wexp7", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp8        = new Var("wexp8", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);
    Var* wexp9        = new Var("wexp9", "W_{exp}", "MeV", 20, 0., 2.0e3, &CVUniverse::GetWexp);

    std::vector<Var*> variables = { 
      n_had_tracks, wexp_cut, michel_count, llr,
      enode01, enode2, enode3, enode4, enode5,
      n_iso_prongs, n_pions,
      //wexp0, wexp1, wexp2, wexp3, wexp4, wexp5, wexp6, wexp7, wexp8, wexp9,
      //fit_vtx_x, fit_vtx_y, fit_vtx_z, 
      n_short_tracks, n_long_tracks,track_reco_meth, n_nodes,
    };

    return variables;
  }
} // namespace run_cut_variables 


std::vector<Variable*> GetCutVariables(
    SignalDefinition signal_definition, bool include_truth_vars = false) {
  std::vector<Variable*> variables;
  switch (signal_definition) {
    case kOnePi:
      variables = run_cut_variables::GetOnePiVariables(include_truth_vars);
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
void LoopAndFillCutVars(const CCPi::MacroUtil& util, CVUniverse* universe,
                        const EDataMCTruth& type,
                        std::vector<Variable*>& variables) {

  std::cout << "Loop and Fill CutVars\n";
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  for(Long64_t i_event=0; i_event < n_entries; ++i_event){
  //for(Long64_t i_event=0; i_event < 10000; ++i_event){
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    ccpi_event::FillCutVars(event, variables); // this function does a lot of work
  } // events
  std::cout << "*** Done ***\n\n";
}


//==============================================================================
// Main
//==============================================================================
void runCutVariables(int signal_definition_int = 0, 
                     const char* plist = "ME1A") {

  // Macro Utility
  std::string data_file_list = GetPlaylistFile(plist, false);
  std::string mc_file_list = GetPlaylistFile(plist, true); 
  const std::string macro("runCutVariables");
  bool do_truth = false, is_grid = false, do_systematics = false;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // INIT VARS, HISTOS, AND EVENT COUNTERS
    const bool do_truth_vars = false;
    std::vector<Variable*> cut_variables = GetCutVariables(util.m_signal_definition,
                                                           do_truth_vars);
    //std::vector<Variable*> cut_variables;

    std::vector<Variable*> ana_variables = GetAnalysisVariables(util.m_signal_definition,
                                                                do_truth_vars);

    std::vector<Variable*> variables;

    variables.reserve( cut_variables.size() + ana_variables.size() ); // preallocate memory
    variables.insert( variables.end(), cut_variables.begin(), cut_variables.end() );
    variables.insert( variables.end(), ana_variables.begin(), ana_variables.end() );

    for (auto var : variables) {
      var->InitializeStackedHists(); 
      var->InitializeDataHists();
    }

  LoopAndFillCutVars(util, util.m_data_universe,              kData, variables);
  LoopAndFillCutVars(util, util.m_error_bands.at("cv").at(0), kMC,   variables);

  // PLOT STUFF
  for (auto v : variables){
    std::string tag;
    if(v->Name() == "wexp0")
      tag = "After: no cuts. Before: Precuts.";
    else if(v->Name() == "wexp1")
      tag = "After: Precuts. Before: MINOS.";
    else if(v->Name() == "wexp2")
      tag = "After: MINOS. Before: N Had Tracks.";
    else if(v->Name() == "wexp3")
      tag = "After: N Had Tracks. Before: Wexp.";
    else if(v->Name() == "wexp4")
      tag = "After: Wexp. Before: Michel.";
    else if(v->Name() == "wexp5")
      tag = "After: Michel. Before: LLR.";
    else if(v->Name() == "wexp6")
      tag = "After: Michel. Before: Node.";
    else if(v->Name() == "wexp7")
      tag = "After: Node. Before: N Iso Prongs.";
    else if(v->Name() == "wexp8")
      tag = "After: N Iso Prongs. Before: Pion Multiplicity.";
    else if(v->Name() == "wexp9")
      tag = "After: All Cuts.";
    double ymax = -1;
    bool do_bwn = true;
    std::cout << "Plotting" << std::endl;
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
  }
}


#endif // runCutVariables.C

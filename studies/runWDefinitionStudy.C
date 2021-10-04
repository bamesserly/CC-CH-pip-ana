//==============================================================================
// What is mc_w? is it Wexp?
// Assume that NO systematics are being analyzed.
//==============================================================================
#include <iostream>
#include <vector>

#include "ccpion_common.h" // GetPlaylistFile
#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/MacroUtil.h"
//#ifndef __CINT__
#include "includes/HadronVariable.h"
#include "includes/Variable.h"
//#endif
//#include "plotting_functions.h"
#include "../xsec/plotting_functions.h"

#include "includes/common_functions.h"
///////
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector4D.h"

// Forward declare my variables because we're hiding the header.
class Variable;
class HadronVariable;

namespace run_w_definition_study {
//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  typedef Variable Var;
  typedef HadronVariable HVar;
  Var* pmu       = new Var("pmu",       "p_{#mu}",      "MeV",         CCPi::GetBinning("pmu"),    &CVUniverse::GetPmu);
  Var* wexp      = new Var("wexp",      "W_{exp}",      "MeV",         CCPi::GetBinning("wexp"),   &CVUniverse::GetWexp);

  bool is_true_var = true;
  Var* mc_w      = new Var("mc_w",      "mc_w",         wexp->m_units, 100, 0., 20.e2, &CVUniverse::GetWgenie,   is_true_var);
  Var* wexp_true = new Var("wexp_true", "W_{exp} True", wexp->m_units, mc_w->m_hists.m_bins_array, &CVUniverse::GetWexpTrue, is_true_var);

  std::vector<Var*> variables = {wexp_true, mc_w};
  return variables;
}
} // namespace run_w_definition_study

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
  std::cout << "nentries " << n_entries << "\n";
  //11531117

  for(Long64_t i_event=0; i_event < 100000; ++i_event){
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
    universe->SetEntry(i_event);

    // For mc, get weight, check signal, and sideband
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    //if (type == kTruth) {
    ccpi_event::FillEfficiencyDenominator(event, variables);
    //  continue;
    //}

    //universe->GetVecElem("mc_initNucVec", 3);
    

  } // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runWDefinitionStudy(std::string plist = "ME1A") {

  /*
  //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>
  using namespace ROOT::Math;
  PxPyPzE4D<double> a(1, 2, 3, 4);
  PxPyPzE4D<double> b(4, 5, 6, 7);
  std::cout << a.P() << "\n";
  std::cout << a.P2() << "\n";
  std::cout << a.Perp2() << "\n";
  //std::cout << a.Dot(b)<< "\n";
  */

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
    const std::string macro("runWDefinitionStudy");
    const bool is_grid = false;
    const bool do_truth = true;
    const bool do_systematics = false;

    CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list, plist, do_truth, is_grid, do_systematics);
    util.PrintMacroConfiguration(macro);

  //=========================================
  // Get variables and initialize their hists
  //=========================================
  std::vector<Variable*> variables = run_w_definition_study::GetVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);



  //=========================================
  // Loop and Fill
  //=========================================
  LoopAndFill(util, util.m_error_bands_truth.at("cv").at(0), kTruth, variables);

  //int PlotTogether(TH1* h1, std::string label1, TH1* h2, std::string label2,
  //             std::string tag, double ymax = -1, bool do_log_scale = false,
  //             bool do_fit = false) {
  PlotTogether(GetVar(variables,  "wexp_true")->m_hists.m_effden.hist, "wexp_true", GetVar(variables, "mc_w")->m_hists.m_effden.hist, "mc_w", "w_comparison_2", 1800.);

  for (auto v : variables) {
    std::string tag = v->Name();
    double ymax = -1;
    bool do_bwn = true;
    std::cout << "Plotting" << std::endl;
    PlotTH1_1(v->m_hists.m_effden.hist, tag, 2000);
    
    //void PlotTH1_1(TH1* h1, std::string tag, double ymax = -1, bool do_log_scale = false, bool do_fit = false) {
    //PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kS), util.m_data_pot, util.m_mc_pot, util.m_signal_definition, v->Name(),"SSB", ymax, do_bwn);
  }

  std::cout << "Success" << std::endl;
}

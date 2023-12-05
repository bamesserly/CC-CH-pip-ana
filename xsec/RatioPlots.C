
#ifndef RatioPlots_C
#define RatioPlots_C

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMatrixD.h"
#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif  // MNVROOT6

#include "includes/MacroUtil.h"
#include "includes/Plotter.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "plotting_functions.h"
/*
void SetPOT(TFile& fin, CCPi::MacroUtil& util) {
  util.m_mc_pot = -1;
  util.m_data_pot = -1;

  // Data
  // TODO make this safer
  PlotUtils::MnvH1D* h_data_pot = (PlotUtils::MnvH1D*)fin.Get("data_pot");
  float data_pot = h_data_pot->GetBinContent(1);
  util.m_data_pot = data_pot;

  // MC
  PlotUtils::MnvH1D* h_mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float mc_pot = h_mc_pot->GetBinContent(1);
  util.m_mc_pot = mc_pot;

  // Ratio
  util.m_pot_scale = util.m_data_pot / util.m_mc_pot;
}
*/
//==============================================================================
// Main
//==============================================================================
void RatioPlots(int signal_definition_int = 0,
                              int plot_errors = 1) {
  // Infiles
  TFile fin("MCXSecInputs_20231203_ME1A-D_mixed_noSys_noIsoProngCut_p4.root", "READ");
  TFile fin1("MCXSecInputs_20231203_ME1A-D_mixed_noSys_noIsoProngCut_p3.root", "READ");
//TFile fin("DataXSecInputs_20231204_ME1A-D_mixed_noSys_p4.root", "READ");
//TFile fin1("DataXSecInputs_20231204_ME1A-D_mixed_noSys_p3.root", "READ");
  cout << "Reading input from " << fin.GetName() << endl;

  // Set up macro utility object...which gets the list of systematics for us...
  // which we need in order to read in HistWrappers...which we don't need at
  // this point...indeed we only need MnvH1Ds...so that's a TODO: write a
  // function that only loads in MnvH1D's from a file, not HWs.
  // Most of these options aren't used in this script. TODO make a CTOR that
  // doesn't require them.
  // INPUT TUPLES
  // Don't actually use the MC chain, only load it to indirectly access it's
  // systematics
  const std::string plist = "ME1A";
  std::string data_file_list = GetPlaylistFile(plist, false);
  std::string mc_file_list = GetPlaylistFile(plist, true);
  //std::string data_file_list = GetTestPlaylist(false);
  //std::string mc_file_list = GetTestPlaylist(true);
 
  // Macro Utility
  const std::string macro("RatiosPlots");
  bool do_truth = false, is_grid = false, do_systematics = true;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  // Get POT from file, not from any chain
//  SetPOT(fin, util);
  util.PrintMacroConfiguration(macro);
  std::cout << "Que pex 0\n";
  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);
/*
  for (auto var : variables) {
    if (var->Name() == sidebands::kFitVarString) {
      var->m_hists.m_stacked_wsideband = StackedHistogram<WSidebandType>(
          var->m_hists.m_label, var->m_hists.m_xlabel,
          var->m_hists.m_bins_array, kNWSidebandTypes,
          sidebands::kWSideband_ColorScheme);
    }
    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    var->LoadDataHistsFromFile(fin);

    // if(var->Name() == "ptmu")
    //  PrintUniverseContent(var->m_hists.m_cross_section);
  }*/

  {  // remove unwanted variables
    ContainerEraser::erase_if(variables, [](Variable* v) {
      return v->Name() == "tpi_mbr";  // || v->Name() == "wexp_fit";
    });

    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "tpi" || v->Name() == "enu"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "thetapi_deg" || v->Name() == "thetamu_deg"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "q2" || v->Name() == "wexp"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "ptmu" || v->Name() == "pzmu"; });
  }

  // Validation 
  if (true){
//    PlotUtils::MnvH1D* data_pot = (PlotUtils::MnvH1D*)fin.Get("data_pot");
//    PlotUtils::MnvH1D* data_pot_1 = (PlotUtils::MnvH1D*)fin1.Get("data_pot");
    PlotUtils::MnvH1D* mc_pot_1 = (PlotUtils::MnvH1D*)fin1.Get("mc_pot");
    PlotUtils::MnvH1D* mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
    double mc_norm = mc_pot->GetBinContent(1)/mc_pot_1->GetBinContent(1);
//    double data_norm = data_pot->GetBinContent(1)/data_pot_1->GetBinContent(1);
    std::cout << "mc Norm = " << mc_norm << "\n";
//    std::cout << "data Norm = " << data_norm << "\n";
    for (auto var : variables) {
      std::string name = var->Name();
      if (var->m_is_true) continue;
      if (name ==  "wexp_fit") continue;
/*      PlotUtils::MnvH1D* num_data_sel = (PlotUtils::MnvH1D*)fin.Get(Form("selection_data_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_data_sel = (PlotUtils::MnvH1D*)fin1.Get(Form("selection_data_%s", name.c_str()));
*/      PlotUtils::MnvH1D* num_mc_sel = (PlotUtils::MnvH1D*)fin.Get(Form("selection_mc_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_mc_sel = (PlotUtils::MnvH1D*)fin1.Get(Form("selection_mc_%s", name.c_str()));
      PlotUtils::MnvH1D* num_mc_noWcut = (PlotUtils::MnvH1D*)fin.Get(Form("noWcut_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_mc_noWcut = (PlotUtils::MnvH1D*)fin1.Get(Form("noWcut_%s", name.c_str()));
/*      PlotUtils::MnvH1D* num_mc_Effnum = (PlotUtils::MnvH1D*)fin.Get(Form("effnum_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_mc_Effnum = (PlotUtils::MnvH1D*)fin1.Get(Form("effnum_%s", name.c_str()));
      PlotUtils::MnvH1D* num_data_BGsub = (PlotUtils::MnvH1D*)fin.Get(Form("bg_subbed_data_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_data_BGsub = (PlotUtils::MnvH1D*)fin1.Get(Form("bg_subbed_data_%s", name.c_str()));
      PlotUtils::MnvH1D* num_data_unfolded = (PlotUtils::MnvH1D*)fin.Get(Form("unfolded_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_data_unfolded = (PlotUtils::MnvH1D*)fin1.Get(Form("unfolded_%s", name.c_str()));
 

     
      PlotUtils::MnvH1D* num_mc_Eff = (PlotUtils::MnvH1D*)fin.Get(Form("efficiency_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_mc_Eff = (PlotUtils::MnvH1D*)fin1.Get(Form("efficiency_%s", name.c_str()));
      PlotUtils::MnvH1D* num_data_xsec = (PlotUtils::MnvH1D*)fin.Get(Form("cross_section_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_data_xsec = (PlotUtils::MnvH1D*)fin1.Get(Form("cross_section_%s", name.c_str()));
      PlotUtils::MnvH1D* num_mc_xsec = (PlotUtils::MnvH1D*)fin.Get(Form("mc_cross_section_%s", name.c_str()));
      PlotUtils::MnvH1D* denom_mc_xsec = (PlotUtils::MnvH1D*)fin1.Get(Form("mc_cross_section_%s", name.c_str()));
*///      PlotRatio(num_data_sel, denom_data_sel, name, data_norm, "Data_Sel", false);
      denom_mc_sel->Scale(mc_norm);
      PlotRatio(num_mc_sel, denom_mc_sel, name, 1., "MC_Sel", false);
      PlotRatio(num_mc_noWcut, denom_mc_noWcut, name, mc_norm, "MC_noWcut", false);
/*      PlotRatio(num_mc_Eff, denom_mc_Eff, name, 1., "Efficiency", false);
      PlotRatio(num_mc_Effnum, denom_mc_Effnum, name, mc_norm, "Effnum", false);
      PlotRatio(num_data_BGsub, denom_data_BGsub, name, data_norm, "BGSub", false);
      PlotRatio(num_data_unfolded, denom_data_unfolded, name, data_norm, "Unfolding", false);

      PlotRatio(num_data_xsec, denom_data_xsec, name, 1., "Data_XSec", false);
      PlotRatio(num_mc_xsec, denom_mc_xsec, name, 1., "MC_XSec", false);
*/    }
  }
}
#endif  // plotCrossSectionFromFile_C

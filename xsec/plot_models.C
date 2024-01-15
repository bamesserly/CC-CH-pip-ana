#ifndef plotCrossSectionFromFile_C
#define plotCrossSectionFromFile_C

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
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "xsec/makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "xsec/plotting_functions.h"

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

//==============================================================================
// Main
//==============================================================================
void plot_models(int signal_definition_int = 1,
                              int plot_errors = 0) {

  // Infiles
  TFile fin("DataXSecInputs_20231230_ALL_mixed_Sys_p4.root", "READ");
  cout << "Reading input from " << fin.GetName() << endl;

  TFile fin2("nuisance_ME_FHC_tracker_v3_0_6_G18_02b_02_11a_CH_50M.root", "READ");

  // Set up macro utility object...which gets the list of systematics for us...
  // which we need in order to read in HistWrappers...which we don't need at
  // this point...indeed we only need MnvH1Ds...so that's a TODO: write a
  // function that only loads in MnvH1D's from a file, not HWs.
  // Most of these options aren't used in this script. TODO make a CTOR that
  // doesn't require them.
  // INPUT TUPLES
  // Don't actually use the MC chain, only load it to indirectly access it's
  // systematics
  const std::string plist = "ME1L";
  std::string data_file_list = GetPlaylistFile(plist, false);
  std::string mc_file_list = GetPlaylistFile(plist, true);
  // std::string data_file_list = GetTestPlaylist(false);
  // std::string mc_file_list = GetTestPlaylist(true);

  // Macro Utility
  const std::string macro("PlotCrossSectionFromFile");
  bool do_truth = false, is_grid = false, do_systematics = true;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  // Get POT from file, not from any chain
  SetPOT(fin, util);
  util.PrintMacroConfiguration(macro);

  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

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
  }

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

  // PLOT cross section
  if (true) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    const double ymax = -1.;
    const bool do_log_scale = false;
    const bool do_bin_width_norm = true;
    for (auto reco_var : variables) {
      if (reco_var->Name() != "mixtpi") continue;
      if (reco_var->Name() == "wexp_fit") continue;
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
                        do_cov_area_norm, include_stat,
                        util.m_signal_definition);

      TH1D* m_mc_cross_section = (TH1D*)fin2.Get("tpi");

      // std::vector<std::string> x_bands =
      // reco_var->m_hists.m_cross_section->GetVertErrorBandNames();
      // if(reco_var->Name() == "ptmu")
      //  for (auto s : x_bands) std::cout << s << "\n";

      // std::cout << reco_var->Name() << "\n";

      Plot_Models(plot_info, reco_var->m_hists.m_cross_section, m_mc_cross_section);

      if (plot_errors)
        PlotCrossSection_ErrorSummary(
            plot_info);  // Adds chi2 label and prints out assumed binning.

      // Various error matrices
      // TMatrixD
      // scaled_mtx(reco_var->m_hists.m_cross_section->GetTotalErrorMatrix());
      // TMatrixD scaled_mtx(m_mc_cross_section->GetStatErrorMatrix());

      // TMatrixD
      // data_stat_err_mtx(reco_var->m_hists.m_cross_section->GetStatErrorMatrix());
      // TMatrixD
      // data_tot_err_mtx(reco_var->m_hists.m_cross_section->GetTotalErrorMatrix());
      // if (plot_errors) PlotMatrix(data_stat_err_mtx, reco_var->Name(),
      // "Data_Stat_Err_Mtx"); if (plot_errors) PlotMatrix(data_tot_err_mtx,
      // reco_var->Name(), "Data_Tot_Err_Mtx");

      // TMatrixD mc_stat_err_mtx(m_mc_cross_section->GetStatErrorMatrix());
      // TMatrixD mc_tot_err_mtx(m_mc_cross_section->GetTotalErrorMatrix());
      // if (plot_errors) PlotMatrix(mc_stat_err_mtx,   reco_var->Name(),
      // "MC_Stat_Err_Mtx"); if (plot_errors) PlotMatrix(mc_tot_err_mtx,
      // reco_var->Name(), "MC_Tot_Err_Mtx");

      // TMatrixD
      // data_sys_err_mtx(reco_var->m_hists.m_cross_section->GetSysErrorMatrix());
      // TMatrixD
      // data_sys_corr_mtx(reco_var->m_hists.m_cross_section->GetSysCorrelationMatrix
      // ()); if (plot_errors) PlotMatrix(data_sys_err_mtx,  reco_var->Name(),
      // "Data_Sys_Err_Mtx"); if (plot_errors) PlotMatrix(data_sys_cor_mtx,
      // reco_var->Name(), "Data_Sys_Cor_Mtx");

      // TMatrixD mc_sys_err_mtx (m_mc_cross_section->GetSysErrorMatrix());
      // TMatrixD mc_sys_cor_mtx(m_mc_cross_section->GetSysCorrelationMatrix());
      // if (plot_errors) PlotMatrix(mc_sys_err_mtx,    reco_var->Name(),
      // "MC_Sys_Err_Mtx"); if (plot_errors) PlotMatrix(mc_sys_cor_mtx,
      // reco_var->Name(), "MC_Sys_Cor_Mtx");

      // if (plot_errors) PlotCovarianceMatrix(scaled_mtx, reco_var->Name(),
      // "Data_StatErrMtx"); for(int c = 0; c < scaled_mtx.GetNcols(); ++c ) {
      //  for(int r = 0; r < scaled_mtx.GetNrows(); ++r ) {
      //    std::cout << scaled_mtx(r, c) << "  ";
      //  }
      //  std::cout << "\n";
      //}

      // for(int c = 0; c < scaled_mtx2.GetNcols(); ++c ) {
      //  for(int r = 0; r < scaled_mtx2.GetNrows(); ++r ) {
      //    std::cout << scaled_mtx2(r, c) << "  ";
      //  }
      //  std::cout << "\n";
      //}

      // std::cout << "data (" << scaled_mtx.GetNrows() << "," <<
      // scaled_mtx.GetNcols() << ")\n"; std::cout << "mc (" <<
      // scaled_mtx2.GetNrows() << "," << scaled_mtx2.GetNcols() << ")\n";

      // PrintChi2Info(plot_info, reco_var->m_hists.m_cross_section,
      // m_mc_cross_section); // this one works
    }
  }

  //============================================================================
}

#endif  // plotCrossSectionFromFile_C

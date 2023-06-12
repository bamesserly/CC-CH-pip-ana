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
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "plotting_functions.h"

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
void plotCrossSectionFromFile(int signal_definition_int = 0,
                              int plot_errors = 0) {
  // Infiles
  TFile fin("DataXSecInputs_Mpions_20230606_ME1A_Multimode_tracked.root", "READ");
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
    if (var->Name() == "bkdtrackedtpi" && var->Name() == "bkdtracklesstpi" &&
        var->Name() == "bkdmixtpi" && var->Name() == "bkdtrackedtpi_true" &&
        var->Name() == "bkdtracklesstpi_true" && var->Name() == "bkdmixtpi_true")
      continue;
    std::cout << "NE pe ne\n";
    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    std::cout << "NE pe ne x 2\n";
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

  PlotUtils::MnvH1D* h_bkdtrackedtpi = (PlotUtils::MnvH1D*)fin.Get("selection_mc_bkdtrackedtpi");
  PlotUtils::MnvH1D* h_bkdtracklesstpi = (PlotUtils::MnvH1D*)fin.Get("selection_mc_bkdtracklesstpi");
  PlotUtils::MnvH1D* h_bkdmixtpi = (PlotUtils::MnvH1D*)fin.Get("selection_mc_bkdmixtpi");

  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  TCanvas cE ("c1","c1");
  TObjArray* stack = new TObjArray();
  h_bkdtrackedtpi->SetTitle("Tracked");
  h_bkdtracklesstpi->SetTitle("Tracless");
  h_bkdmixtpi->SetTitle("Mix");
  h_bkdtrackedtpi->GetYaxis()->SetTitle("Events/MeV");
  h_bkdtrackedtpi->Scale(util.m_data_pot/util.m_mc_pot, "width");
  h_bkdtracklesstpi->Scale(util.m_data_pot/util.m_mc_pot, "width");
  h_bkdmixtpi->Scale(util.m_data_pot/util.m_mc_pot, "width");
  cE.SetLogx();


  stack->Add(h_bkdtrackedtpi);
  stack->Add(h_bkdmixtpi);
  stack->Add(h_bkdtracklesstpi);
  mnvPlotter.DrawStackedMC(stack, 1.0, "TR", 2, 1, 3001, "Tpi (MeV)");
  mnvPlotter.AddHistoTitle("Tpi Breakdown", 0.05);

  std::string plotname = "Stacked_Tpi";
  mnvPlotter.MultiPrint(&cE, plotname , "png");


  // PLOT Event Selection, BGs (error)
  if (true) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    bool do_cov_area_norm = false;
    for (auto var : variables) {
      if (var->Name() ==  "wexp_fit") continue;
      std::cout << var->Name() << "\n";
      do_cov_area_norm = false;
      Plotter plot_info(var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      bool do_bin_width_norm = true, do_log_scale = false, do_bg = true;
      bool do_tuned_bg = false;
//      if (var->Name() == "mtpi") do_log_scale = true;
      // selection and BG error before tuning
      if (!var->m_is_true) {
        PlotVar_Selection(plot_info, -1., do_log_scale, do_bg, do_tuned_bg,
                          do_bin_width_norm);
        if (plot_errors) PlotVar_ErrorSummary(plot_info);
        if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

        // selection and BG error after tuning
        do_tuned_bg = true;
        PlotVar_Selection(plot_info, -1., do_log_scale, do_bg, do_tuned_bg,
                          do_bin_width_norm);
        if (plot_errors) PlotVar_ErrorSummary(plot_info);
        if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);
      }
    }
  }

  // PLOT Efficiency & Migration
  if (true) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;

    for (auto var : variables) {
      if (var->Name() ==  "wexp_fit") continue;
      // var->LoadMCHistsFromFile(fin, util.m_error_bands);

      const Plotter plot_info(
          var, util.m_mc_pot, util.m_data_pot, do_frac_unc, do_cov_area_norm,
          include_stat, util.m_signal_definition);

      // Efficiency
      if (var->m_is_true) {
        PlotUtils::MnvH1D* eff =
            (PlotUtils::MnvH1D*)var->m_hists.m_effnum.hist->Clone(uniq());
        PlotUtils::MnvH1D* effnum =
            (PlotUtils::MnvH1D*)var->m_hists.m_effnum.hist->Clone(uniq());
        PlotUtils::MnvH1D* effden =
            (PlotUtils::MnvH1D*)var->m_hists.m_effden.hist->Clone(uniq());
        eff->Divide(effnum, effden);
        var->m_hists.m_efficiency = (PlotUtils::MnvH1D*)eff->Clone(uniq());

        double ymax = -1.;
        bool do_log_scale = false;
        bool do_bg = true;
        bool do_tuned_bg = true;
        PlotMC(eff, plot_info, Form("Efficiency_%s", var->Name().c_str()),
               0.075, "Efficiency");
        if (plot_errors) PlotEfficiency_ErrorSummary(plot_info);
      }

      // Migration
      if (!var->m_is_true) {
        PlotUtils::MnvH2D* mig =
            (PlotUtils::MnvH2D*)var->m_hists.m_migration.hist->Clone(uniq());
        PlotMigration_AbsoluteBins(mig, var->Name());
        PlotMigration_VariableBins(mig, var->Name());
      }
    }
  }

  // PLOT Background Subtraction
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    for (auto var : variables) {
      if (var->Name() ==  "wexp_fit") continue;
      if (var->m_is_true) continue;

      Plotter plot_info(var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      double ymax = -1.;
      bool do_log_scale = false;
      bool do_bg = true;
      bool do_tuned_bg = true;
      bool do_bin_width_norm = true;
      if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

      do_tuned_bg = false;
      if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

      Plot_BGSub(plot_info, ".", ymax, do_log_scale, do_bin_width_norm);

      if (plot_errors) PlotBGSub_ErrorSummary(plot_info);
    }
  }

  // PLOT W Sideband Fit
  if (false) {
    const bool do_frac_unc = true;
    const bool do_cov_area_norm = false;
    const bool include_stat = true;

    Plotter plot_info(util.m_mc_pot, util.m_data_pot,
                                     do_frac_unc, do_cov_area_norm,
                                     include_stat, util.m_signal_definition);

    PlotUtils::MnvH1D* loW_fit_wgt = (PlotUtils::MnvH1D*)fin.Get("loW_fit_wgt");
    PlotUtils::MnvH1D* midW_fit_wgt =
        (PlotUtils::MnvH1D*)fin.Get("midW_fit_wgt");
    PlotUtils::MnvH1D* hiW_fit_wgt = (PlotUtils::MnvH1D*)fin.Get("hiW_fit_wgt");

    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, loW_fit_wgt, "WSidebandFit_loW");

    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, midW_fit_wgt,
                                    "WSidebandFit_midW");
    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, hiW_fit_wgt, "WSidebandFit_hiW");

    // Wexp distribution, stacked, all cuts except for W, pre-fit
    std::string tag;
    double ymax = -1;
    Variable* var = GetVar(variables, sidebands::kFitVarString);
    PlotWSidebandStacked(var, var->m_hists.m_wsideband_data,
                         var->GetStackArray(static_cast<WSidebandType>(0)),
                         util.m_data_pot, util.m_mc_pot,
                         util.m_signal_definition, tag, ymax);

    // TODO plot pre/postfit
    // for (auto var : variables) {
    //  tag = "SidebandRegion";
    //  bool do_prefit = true;
    //  bool do_bin_width_norm = true;
    //  CVUniverse* universe = util.m_error_bands.at("cv").at(0);
    //  PlotFittedW(var, *universe, hw_loW_fit_wgt, hw_midW_fit_wgt,
    //              hw_hiW_fit_wgt, util.m_data_pot, util.m_mc_pot,
    //              util.m_signal_definition, do_prefit, tag, ymax,
    //              do_bin_width_norm);
    //  do_prefit = false;
    //  PlotFittedW(var, *universe, hw_loW_fit_wgt, hw_midW_fit_wgt,
    //              hw_hiW_fit_wgt, util.m_data_pot, util.m_mc_pot,
    //              util.m_signal_definition, do_prefit, tag, ymax,
    //              do_bin_width_norm);
    //}
  }

  // PLOT unfolded
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    const double ymax = -1.;
    const bool do_log_scale = false;
    const bool do_bin_width_norm = true;
    for (auto reco_var : variables) {
      if (reco_var->Name() ==  "wexp_fit") continue;
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      Plot_Unfolded(plot_info, reco_var->m_hists.m_unfolded,
                    true_var->m_hists.m_effnum.hist);
      if (plot_errors) PlotUnfolded_ErrorSummary(plot_info);
    }
  }

  // PLOT cross section
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    const double ymax = -1.;
    const bool do_log_scale = false;
    const bool do_bin_width_norm = true;
    for (auto reco_var : variables) {
      if (reco_var->Name() ==  "wexp_fit") continue;
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      PlotUtils::MnvH1D* m_mc_cross_section = (PlotUtils::MnvH1D*)fin.Get(
          Form("mc_cross_section_%s", reco_var->Name().c_str()));

      // std::vector<std::string> x_bands =
      // reco_var->m_hists.m_cross_section->GetVertErrorBandNames();
      // if(reco_var->Name() == "ptmu")
      //  for (auto s : x_bands) std::cout << s << "\n";

      // std::cout << reco_var->Name() << "\n";

      Plot_CrossSection(plot_info, reco_var->m_hists.m_cross_section,
                        m_mc_cross_section);
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

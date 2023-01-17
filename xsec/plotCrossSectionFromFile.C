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
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "plotting_functions.h"
#include "includes/Binning.h"

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
                              int plot_errors = 1) {
  // Infiles
  TFile fin("DataXSecInputs_0010_ME1A_0_2023-01-10_BenSigFixedWeight.root", "READ");
  cout << "Reading input from " << fin.GetName() << endl;

  TFile finCCPi("DataXSecInputs_0010_ME1A_0_2023-01-10_AaronSigfixedWeight.root", "READ");
  //    TFile
  //    finCCPi("/minerva/app/users/granados/cmtuser/Minerva_v22r1p1_CCPionInc/Ana/CCPionInc/ana/ME_CCNuPionInc_Ana/DataXSec_20211010_NewTupla.root",
  //    "READ");

  //    TFile finCCPi("../ME_CCNuPionInc_Ana/DataXSec_20210901_CCPi.root",
  //    "READ");

  // Set up macro utility object...which gets the list of systematics for us...
  // which we need in order to read in HistWrappers...which we don't need at
  // this point...indeed we only need MnvH1Ds...so that's a TODO: write a
  // function that only loads in MnvH1D's from a file, not HWs.
  // Most of these options aren't used in this script. TODO make a CTOR that
  // doesn't require them.
  // INPUT TUPLES
  // Don't actually use the MC chain, only load it to indirectly access it's
  // systematics
  const std::string plist = "ME1B";
  std::string data_file_list = GetPlaylistFile(plist, false);
  std::string mc_file_list = GetPlaylistFile(plist, true);

  // Macro Utility
  const std::string macro("PlotCrossSectionFromFile");
  bool do_truth = false, is_grid = false, do_systematics = true;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // Get POT from file, not from any chain
  SetPOT(fin, util);

  std::string data_file_list_CCPi = GetPlaylistFileCCPi(plist, false);
  std::string mc_file_list_CCPi = GetPlaylistFileCCPi(plist, true);

  CCPi::MacroUtil utilCCPi(signal_definition_int, mc_file_list_CCPi,
                           data_file_list_CCPi, plist, do_truth, is_grid,
                           do_systematics);
  utilCCPi.PrintMacroConfiguration(macro);
  SetPOT(finCCPi, utilCCPi);

  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  for (auto var : variables) {
    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    var->LoadDataHistsFromFile(fin);

    // if(var->Name() == "ptmu")
    //  PrintUniverseContent(var->m_hists.m_cross_section);
  }

  {  // remove unwanted variables
    ContainerEraser::erase_if(variables, [](Variable* v) {
      return v->Name() == "tpi_mbr" || v->Name() == "wexp_fit";
    });
    /*
      ContainerEraser::erase_if(variables, [](Variable* v) {
          return v->Name() == "tpi" || v->Name() == "enu"; });
      ContainerEraser::erase_if(variables, [](Variable* v) {
          return v->Name() == "thetapi_deg" || v->Name() == "thetamu_deg"; });
      ContainerEraser::erase_if(variables, [](Variable* v) {
          return v->Name() == "q2" || v->Name() == "wexp"; });
      ContainerEraser::erase_if(variables, [](Variable* v) {
          return v->Name() == "ptmu" || v->Name() == "pzmu"; });
    */
  }

  // Ratios MAD and CCPionInc
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = false;
    bool do_cov_area_norm = false;
    bool fixRange = true;
    PlotUtils::MnvH1D* MADpotMC = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
    PlotUtils::MnvH1D* MADpotdata = (PlotUtils::MnvH1D*)fin.Get("data_pot");
    PlotUtils::MnvH1D* CCPipotMC = (PlotUtils::MnvH1D*)finCCPi.Get("mc_pot");
    PlotUtils::MnvH1D* CCPipotdata =
        (PlotUtils::MnvH1D*)finCCPi.Get("data_pot");
    double MADnorm = MADpotdata->GetBinContent(1) / MADpotMC->GetBinContent(1);
    double CCPinorm =
        CCPipotdata->GetBinContent(1) / CCPipotMC->GetBinContent(1);
    double MC_MADCCPinorm =
        MADpotMC->GetBinContent(1) / CCPipotMC->GetBinContent(1);
    double data_MADCCPinorm =
        MADpotdata->GetBinContent(1) / CCPipotdata->GetBinContent(1);

    std::cout << "MADnorm = " << MADnorm << "\n";
    std::cout << "CCPinorm = " << CCPinorm << "\n";
    std::cout << "mc_MADCCPinorm = " << MC_MADCCPinorm << "\n";
    std::cout << "data_MADCCPinorm = " << data_MADCCPinorm << "\n";

    std::cout << "MADPOT data from Util = " << util.m_data_pot << "\n";
    std::cout << "MADPOT mc from Util = " << util.m_mc_pot << "\n";
    std::cout << "CCPiPOT data from UtilCCPi = " << utilCCPi.m_data_pot << "\n";
    std::cout << "CCPiPOT mc from UtilCCPi = " << utilCCPi.m_mc_pot << "\n";

    std::cout << "mc Norm = " << utilCCPi.m_mc_pot / util.m_mc_pot << "\n";
    std::cout << "data Norm = " << utilCCPi.m_data_pot / util.m_data_pot
              << "\n";

    std::cout << "POT Scale = " << util.m_pot_scale << "\n";
    util.m_pot_scale = util.m_pot_scale * MC_MADCCPinorm;
    util.m_mc_pot = util.m_mc_pot * MC_MADCCPinorm;
    util.m_data_pot = util.m_data_pot * data_MADCCPinorm;

    std::vector<std::string> sec;
    sec.push_back("selection_data");
    sec.push_back("selection_mc");
    sec.push_back("BGSub_data");
    sec.push_back("BGSub_MC");
    sec.push_back("Unfolded_Data");
    sec.push_back("Unfolded_MC");
    sec.push_back("efficiency");
    sec.push_back("cross_section");
    sec.push_back("mc_cross_section");

    for (auto var : variables) {
      auto reco_var = var;
      if (var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      PlotUtils::MnvH1D* Num;
      PlotUtils::MnvH1D* Denom;

      for (auto s : sec) {
        std::string curr = s + "_" + var->Name();
        std::cout << curr << "\t" << var->Name() << "\n";
        double Norm = 1.;
        if (s == "selection_data" || s == "BGSub_data" ||
            s == "cross_section" || s == "Unfolded_Data") {
          Norm = data_MADCCPinorm;
        }
        if (s == "selection_mc" || s == "BGSub_MC" || s == "mc_cross_section" ||
            s == "Unfolded_MC" || s == "efficiency") {
          Norm = MC_MADCCPinorm;
        }
        if (s == "BGSub_data") {
          curr = "bg_subbed_data_" + var->Name();
          Num = (PlotUtils::MnvH1D*)fin.Get(Form("%s", curr.c_str()));
          Denom = (PlotUtils::MnvH1D*)finCCPi.Get(Form("%s", curr.c_str()));
        } else if (s == "BGSub_MC") {
          curr = "effnum_" + var->Name();
          Num = (PlotUtils::MnvH1D*)fin.Get(Form("%s", curr.c_str()));
          Denom = (PlotUtils::MnvH1D*)finCCPi.Get(Form("%s", curr.c_str()));
        }
        if (s == "Unfolded_Data") {
          curr = "unfolded_" + var->Name();
          Num = (PlotUtils::MnvH1D*)fin.Get(Form("%s", curr.c_str()));
          Denom = (PlotUtils::MnvH1D*)finCCPi.Get(Form("%s", curr.c_str()));
        } else if (s == "Unfolded_MC") {
          curr = "effnum_" + true_var->Name();
          Num = (PlotUtils::MnvH1D*)fin.Get(Form("%s", curr.c_str()));
          Denom = (PlotUtils::MnvH1D*)finCCPi.Get(Form("%s", curr.c_str()));
        } else {
          Num = (PlotUtils::MnvH1D*)fin.Get(Form("%s", curr.c_str()));
          Denom = (PlotUtils::MnvH1D*)finCCPi.Get(Form("%s", curr.c_str()));
        }
        std::cout << "Norm = " << Norm << "\n";
        //	Norm = 1;
        //        std::cout << "Norm = " << Norm << "\n";
        PlotRatio(Num, Denom, var->Name(), Norm, s, fixRange);
      }
    }
  }

  // Comparing with Aaron's Cross Section
  if (true){
    int nbins = CCPi::GetBinning("q2_GeV").GetSize()-1;  
    const Double_t* bins = CCPi::GetBinning("q2_GeV").GetArray();
    std::cout << "Numero de bins " << nbins << "\n";
//    for (int i = 0; i < nbins; ++i)
//      std::cout << "Min " << bins[i] << " max " << bins[i+1] << "\n";  
    PlotUtils::MnvH1D *q2_Aaron_tesis = new PlotUtils::MnvH1D("q2_Aaron_tesis","q2_Aaron_tesis",nbins, CCPi::GetBinning("q2_GeV").GetArray());
    
    std::vector<double> xsec;
    xsec.push_back(90.26e-42); 
    xsec.push_back(104.81e-42);
    xsec.push_back(219.20e-42);
    xsec.push_back(402.76e-42);
    xsec.push_back(344.86e-42);
    xsec.push_back(272.35e-42);
    xsec.push_back(208.85e-42);
    xsec.push_back(279.50e-42);
    xsec.push_back(212.70e-42);
    xsec.push_back(98.25-42);
    xsec.push_back(62.76e-42);
    xsec.push_back(20.81e-42);
//    for (int i = 0; i < nbins; ++i)
//      std::cout << "Min " << bins[i] << " max " << bins[i+1] << "\n"; 
    for (int i = 1; i <= nbins; ++i)
      q2_Aaron_tesis->SetBinContent(i, xsec[i-1]);

//    for (int i = 0; i < nbins; ++i)
//      std::cout << "Min " << q2_Aaron_tesis->GetBinLowEdge(i+1) << " max " << q2_Aaron_tesis->GetXaxis()->GetBinLowEdge(i+1) + q2_Aaron_tesis->GetXaxis()->GetBinWidth(i+1)  << " BinContent " << q2_Aaron_tesis->GetBinContent(i+1) << "\n";

    PlotUtils::MnvH1D* BenXSecMC = (PlotUtils::MnvH1D*)fin.Get("mc_cross_section_q2_GeV");
    PlotUtils::MnvH1D* AaronXSecMC = (PlotUtils::MnvH1D*)finCCPi.Get("mc_cross_section_q2_GeV");
    PlotUtils::MnvH1D* BenXSecdata = (PlotUtils::MnvH1D*)fin.Get("cross_section_q2_GeV");
    PlotUtils::MnvH1D* AaronXSecdata = (PlotUtils::MnvH1D*)finCCPi.Get("cross_section_q2_GeV");    
    
    PlotRatio2(BenXSecMC, AaronXSecMC, q2_Aaron_tesis, "Q^2", 1., "mc", false, true);
    PlotRatio2(BenXSecdata, AaronXSecdata, q2_Aaron_tesis, "Q^2", 1., "data", false, true);  
    PlotRatio2(BenXSecMC, AaronXSecMC, q2_Aaron_tesis, "Q^2", 1., "mc", false, true);    
//    for (int i = 0; i < nbins; ++i)
//      std::cout << q2_Aaron_tesis->GetBinContent(i+1) << " " << BenXSecMC->GetBinContent(i+1) << " " << AaronXSecMC->GetBinContent(i+1) << " " << BenXSecdata->GetBinContent(i+1) << " " << AaronXSecdata->GetBinContent(i+1) << "\n";
  }

  // PLOT Event Selection, BGs (error)
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    bool do_cov_area_norm = false;
    for (auto var : variables) {
      std::cout << var->Name() << "\n";
      do_cov_area_norm = false;
      EventSelectionPlotInfo plot_info(var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      bool do_bin_width_norm = true, do_log_scale = false, do_bg = true;
      bool do_tuned_bg = false;

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
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;

    for (auto var : variables) {
      // var->LoadMCHistsFromFile(fin, util.m_error_bands);

      const EventSelectionPlotInfo plot_info(
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
      if (var->m_is_true) continue;

      EventSelectionPlotInfo plot_info(var, util.m_mc_pot, util.m_data_pot,
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

    EventSelectionPlotInfo plot_info(util.m_mc_pot, util.m_data_pot,
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
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      EventSelectionPlotInfo plot_info(reco_var, util.m_mc_pot, util.m_data_pot,
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
    bool do_log_scale = false;
    const bool do_bin_width_norm = true;
    for (auto reco_var : variables) {
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));
      if (reco_var->Name() == "q2_GeV") do_log_scale = true;
      else do_log_scale = false;
      EventSelectionPlotInfo plot_info(reco_var, util.m_mc_pot, util.m_data_pot,
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
                        m_mc_cross_section, ".", -1, do_log_scale);
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

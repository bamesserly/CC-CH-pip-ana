#ifndef plotNuisance_C
#define plotNuisance_C

#include <cassert>  // !! must compile in debug mode root script.C++g
#include <iostream>
#include <sstream>

#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif

#include "Constants.h"
#include "PlotUtils/MnvColors.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvVertErrorBand.h"
#include "Plotter.h"
#include "SignalDefinition.h"
#include "Systematics.h"
#include "TArrayD.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TCollection.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TList.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TText.h"
#include "Variable.h"
#include "includes/MacroUtil.h"
#include "includes/Plotter.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "includes/myPlotStyle.h"
#include "xsec/makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "xsec/plotting_functions.h"        // RebinQ2Plot

void plot_all_models(Plotter p, MnvH1D* data,
                     std::map<std::string, PlotUtils::MnvH1D*> mc,
                     std::string outdir = ".", double ymax = -1,
                     bool do_log_scale = false, bool do_bin_width_norm = true) {
  std::cout << "Plotting CrossSection " << p.m_variable->Name() << "\n";

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  for (const auto& i : mc) assert(i.second);

  // Clone the input hists. Don't modify originals.
  PlotUtils::MnvH1D* data_xsec = nullptr;
  std::map<std::string, PlotUtils::MnvH1D*> mc_xsec;

  // rebin q2
  if (p.m_variable->Name() == "q2") {
    data_xsec = RebinQ2Plot(*data);
    for (auto i : mc) {
      mc_xsec[i.first] = RebinQ2Plot(*i.second);
    }
  } else {
    data_xsec = (PlotUtils::MnvH1D*)data->Clone("data");
    for (const auto& i : mc)
      mc_xsec[i.first] = (PlotUtils::MnvH1D*)i.second->Clone(i.first.c_str());
  }

  data_xsec->SetTitle("Data");

  TCanvas canvas("c1", "c1");

  // Get Hists
  TH1D* data_xsec_w_tot_error = new TH1D(data_xsec->GetCVHistoWithError());
  TH1D* data_xsec_w_stat_error = new TH1D(data_xsec->GetCVHistoWithStatError());

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
  if (p.m_variable->Name() == "q2") {
    canvas.SetLogx();
  }

  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // Y-label
  p.m_mnv_plotter.axis_title_offset_y = 1.5;

  // X label
  p.SetXLabel(data_xsec_w_tot_error);
  p.SetXLabel(data_xsec_w_stat_error);
  p.SetXLabel(data_xsec);
  for (auto i : mc_xsec) p.SetXLabel(i.second);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = 1.;
  }

  // Bin Width Normalization, Y-axis label, and 10^-42 shift
  if (do_bin_width_norm) {
    // data_xsec_w_tot_error ->Scale(1.e38, "width");
    // data_xsec_w_stat_error->Scale(1.e38, "width");
    // mc_xsec_w_stat_error  ->Scale(1.e38, "width");

    data_xsec_w_tot_error->Scale(1.e42, "width");
    data_xsec_w_stat_error->Scale(1.e42, "width");
    data_xsec->Scale(1.e42, "width");
    for (auto i : mc_xsec) i.second->Scale(1.e42, "width");

    // Y label
    // std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel + "
    // (10^{-38} cm^{2}/" + p.m_variable->m_units + "/nucleon)";
    std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel +
                        " (10^{-42} cm^{2}/" + p.m_variable->m_units +
                        "/nucleon)";
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    for (auto i : mc_xsec) i.second->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Create a TObjArray from the input vector of model xsecs
  TObjArray* mc_xsec_arr = new TObjArray;
  for (const auto& i : mc_xsec) {
    mc_xsec_arr->Add(i.second);
  }

  // Plot it
  const bool use_hist_titles = true;
  p.m_mnv_plotter.DrawDataMCVariations(data_xsec, mc_xsec_arr, pot_scale, "TR",
                                       use_hist_titles);

  /*
  // Add chi2 label
  {
    const bool use_data_error_mtx = true;
    const bool use_only_shape_errors = false;
    const bool use_model_stat =
        false;  // model statistical errors -- these are small if you use a
                // priori effden, but atm I'm not. So keep this off.
    // p.m_mnv_plotter.AddChi2Label(data_xsec, mc_xsec, pot_scale, "TR", 0.03,
    // -0.175, use_data_error_mtx, use_only_shape_errors); // this auto turns on
    // model stat errors, I've manually turned it off within the function.

    // std::cout << "use_overflow is " << p.m_mnv_plotter.chi2_use_overflow_err
    // << "\n";
    // Check chi2
    int ndf = -1;
    // double pot_scale = 1.;
    double chi2 = p.m_mnv_plotter.Chi2DataMC(
        data_xsec, mc_xsec, ndf, pot_scale, use_data_error_mtx,
        use_only_shape_errors, use_model_stat);

    // std::cout << p.m_variable->Name() << "\n";
    // std::cout << "   chi2 = "         << chi2     << "\n";
    // std::cout << "   ndf = "          << ndf      << "\n";
    // std::cout << "   chi2/ndf = "     << chi2/ndf << "\n";

    // add label manually
    if (p.m_variable->Name() == "tpi") ndf = 6;
    if (p.m_variable->Name() == "enu") ndf = 6;
    if (p.m_variable->Name() == "pzmu") ndf = 9;
    if (p.m_variable->Name() == "pmu") ndf = 8;
    if (p.m_variable->Name() == "wexp") ndf = 4;
    char* words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chi2, ndf,
                       chi2 / (Double_t)ndf);
    int align = 33;
    p.m_mnv_plotter.AddPlotLabel(words, 0.8, 0.745, 0.03, 1, 62, align);
  }
  */

  // POT info
  // -1 --> don't do mc POT
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, -1, 0.3, 0.88);
  else {
    // p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, -1, 0.3, 0.88);
    p.m_mnv_plotter.WriteNorm("POT-Normalized", 0.3, 0.88, 0.03);
    p.m_mnv_plotter.WriteNorm(Form("Data POT: %.2E", p.m_data_pot), 0.3,
                              0.88 - 0.03, 0.03);
  }

  // Change max number of y-axis digits
  // std::cout << "Old max digits = " << TGaxis::GetMaxDigits() << "\n";
  // TGaxis::SetMaxDigits(4);
  // std::cout << "New max digits = " << TGaxis::GetMaxDigits() << "\n";

  // Plot Title
  // minerva plots don't use titles
  // p.m_mnv_plotter.title_size = 0.05;
  // p.SetTitle("Cross Section " + GetSignalName(p.m_signal_definition));

  // Print .png
  std::string logy_str = do_log_scale ? "_logscale" : "";
  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("%s/CrossSection_%s_%s_%s%s%s", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
  delete data_xsec;
  delete data_xsec_w_tot_error;
  delete data_xsec_w_stat_error;
}

void plot_one_model(Plotter p, MnvH1D* data, TH1D* mc, std::string outdir = ".",
                    double ymax = -1, bool do_log_scale = false,
                    bool do_bin_width_norm = true) {
  std::cout << "Plotting CrossSection " << p.m_variable->Name() << "\n";

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

  PlotUtils::MnvH1D* data_xsec = nullptr;
  TH1D* mc_xsec = nullptr;

  if (p.m_variable->Name() == "q2") {
    data_xsec = RebinQ2Plot(*data);
    mc_xsec = RebinQ2Plot(*mc);
  } else {
    data_xsec = (PlotUtils::MnvH1D*)data->Clone("data");
    mc_xsec = dynamic_cast<TH1D*>(mc->Clone("mc"));
  }

  data_xsec->SetTitle("Data");

  TCanvas canvas("c1", "c1");

  // Get Hists
  TH1D* data_xsec_w_tot_error = new TH1D(data_xsec->GetCVHistoWithError());
  TH1D* data_xsec_w_stat_error = new TH1D(data_xsec->GetCVHistoWithStatError());

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
  if (p.m_variable->Name() == "q2") {
    canvas.SetLogx();
  }

  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // Y-label
  p.m_mnv_plotter.axis_title_offset_y = 1.5;

  // X label
  p.SetXLabel(data_xsec_w_tot_error);
  p.SetXLabel(data_xsec_w_stat_error);
  p.SetXLabel(mc_xsec);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = 1.;
  }

  // Bin Width Normalization, Y-axis label, and 10^-42 shift
  if (do_bin_width_norm) {
    // data_xsec_w_tot_error ->Scale(1.e38, "width");
    // data_xsec_w_stat_error->Scale(1.e38, "width");
    // mc_xsec_w_stat_error  ->Scale(1.e38, "width");

    data_xsec_w_tot_error->Scale(1.e42, "width");
    data_xsec_w_stat_error->Scale(1.e42, "width");
    mc_xsec->Scale(1.e42, "width");

    // Y label
    // std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel + "
    // (10^{-38} cm^{2}/" + p.m_variable->m_units + "/nucleon)";
    std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel +
                        " (10^{-42} cm^{2}/" + p.m_variable->m_units +
                        "/nucleon)";
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    mc_xsec->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Draw
  const bool use_hist_titles = true;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(data_xsec_w_tot_error, mc_xsec,
                                          pot_scale, "TR", use_hist_titles);

  /*
  // Add chi2 label
  {
    const bool use_data_error_mtx = true;
    const bool use_only_shape_errors = false;
    const bool use_model_stat =
        false;  // model statistical errors -- these are small if you use a
                // priori effden, but atm I'm not. So keep this off.
    // p.m_mnv_plotter.AddChi2Label(data_xsec, mc_xsec, pot_scale, "TR", 0.03,
    // -0.175, use_data_error_mtx, use_only_shape_errors); // this auto turns on
    // model stat errors, I've manually turned it off within the function.

    // std::cout << "use_overflow is " << p.m_mnv_plotter.chi2_use_overflow_err
    // << "\n";
    // Check chi2
    int ndf = -1;
    // double pot_scale = 1.;
    double chi2 = p.m_mnv_plotter.Chi2DataMC(
        data_xsec, mc_xsec, ndf, pot_scale, use_data_error_mtx,
        use_only_shape_errors, use_model_stat);

    // std::cout << p.m_variable->Name() << "\n";
    // std::cout << "   chi2 = "         << chi2     << "\n";
    // std::cout << "   ndf = "          << ndf      << "\n";
    // std::cout << "   chi2/ndf = "     << chi2/ndf << "\n";

    // add label manually
    if (p.m_variable->Name() == "tpi") ndf = 6;
    if (p.m_variable->Name() == "enu") ndf = 6;
    if (p.m_variable->Name() == "pzmu") ndf = 9;
    if (p.m_variable->Name() == "pmu") ndf = 8;
    if (p.m_variable->Name() == "wexp") ndf = 4;
    char* words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chi2, ndf,
                       chi2 / (Double_t)ndf);
    int align = 33;
    p.m_mnv_plotter.AddPlotLabel(words, 0.8, 0.745, 0.03, 1, 62, align);
  }
  */

  // POT info
  // -1 --> don't do mc POT
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, -1, 0.3, 0.88);
  else {
    // p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, -1, 0.3, 0.88);
    p.m_mnv_plotter.WriteNorm("POT-Normalized", 0.3, 0.88, 0.03);
    p.m_mnv_plotter.WriteNorm(Form("Data POT: %.2E", p.m_data_pot), 0.3,
                              0.88 - 0.03, 0.03);
  }

  //// Label re: error bars
  // p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.6, 0.772, 0.03);

  // Change max number of y-axis digits
  // std::cout << "Old max digits = " << TGaxis::GetMaxDigits() << "\n";
  // TGaxis::SetMaxDigits(4);
  // std::cout << "New max digits = " << TGaxis::GetMaxDigits() << "\n";

  // Plot Title
  // minerva plots don't use titles
  // p.m_mnv_plotter.title_size = 0.05;
  // p.SetTitle("Cross Section " + GetSignalName(p.m_signal_definition));

  // Print .png
  std::string logy_str = do_log_scale ? "_logscale" : "";
  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("%s/CrossSection_%s_%s_%s%s%s", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
  delete data_xsec;
  delete mc_xsec;
  delete data_xsec_w_tot_error;
  delete data_xsec_w_stat_error;
}

void set_POT(TFile& fin, CCPi::MacroUtil& util) {
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
void plotNuisance(int signal_definition_int = 1, int plot_errors = 0) {
  // Data xsec file input
  TFile fin1("DataXSecInputs_20231230_ALL_mixed_Sys_p4.root", "READ");
  std::cout << "Reading data input from " << fin1.GetName() << "\n";

  std::map<std::string, std::string> models;
  models["GENIEv3_G18_02a"] =
      "nuisance/"
      "nuisance_ME_FHC_tracker_GENIE_v3_0_6_G18_02a_02_11a_CH_50M.root";
  models["GENIEv3_G18_02b"] =
      "nuisance/"
      "nuisance_ME_FHC_tracker_GENIE_v3_0_6_G18_02b_02_11a_CH_50M.root";

  const std::string plist = "ME1L";
  std::string data_file_list = GetPlaylistFile(plist, false);
  std::string mc_file_list = GetPlaylistFile(plist, true);

  // Macro Utility
  const std::string macro("Plot Nuisance");
  bool do_truth = false, is_grid = false, do_systematics = true;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);

  // Get POT from file, not from any chain
  set_POT(fin1, util);
  util.PrintMacroConfiguration(macro);

  // Initialize variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  // Load histogram contents into variables from input file.
  for (auto var : variables) {
    if (var->Name() == sidebands::kFitVarString) {
      var->m_hists.m_stacked_wsideband = StackedHistogram<WSidebandType>(
          var->m_hists.m_label, var->m_hists.m_xlabel,
          var->m_hists.m_bins_array, kNWSidebandTypes,
          sidebands::kWSideband_ColorScheme);
    }
    var->LoadMCHistsFromFile(fin1, util.m_error_bands);
    var->LoadDataHistsFromFile(fin1);
  }

  // remove unwanted variables
  ContainerEraser::erase_if(variables, [](Variable* v) {
    return v->Name() == "tpi_mbr" || v->Name() == "wexp_fit" ||
           v->Name() == "thetamu_deg" || v->Name() == "ehad";
  });

  const bool do_frac_unc = true;
  const bool include_stat = true;
  const bool do_cov_area_norm = false;
  const double ymax = -1.;
  const bool do_log_scale = false;
  const bool do_bin_width_norm = true;
  for (auto reco_var : variables) {
    if (reco_var->m_is_true) continue;

    // Set up plotting properties
    Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
                      do_cov_area_norm, include_stat, util.m_signal_definition);

    // Get nuisance prediction from input file
    std::string var_name;
    if (reco_var->Name() == "mixtpi")
      var_name = "tpi";
    else if (reco_var->Name() == "mixthetapi_deg")
      var_name = "thpi";
    else
      var_name = reco_var->Name();

    // Plot data vs a single model
    if (false) {
      auto model = models.begin();
      TFile fin(model->second.c_str(), "READ");
      TH1D* mc_cross_section = (TH1D*)fin.Get(var_name.c_str());
      assert(mc_cross_section);
      mc_cross_section->SetTitle(model->first.c_str());
      plot_one_model(plot_info, reco_var->m_hists.m_cross_section,
                     mc_cross_section);
      delete mc_cross_section;
      fin.Close()
    }

    // Plot data vs all models on a single plot
    if (true) {
      // First get the model xsec plots from their files
      std::map<std::string, PlotUtils::MnvH1D*> mc_cross_sections;
      for (const auto& model : models) {
        std::string description = model.first;
        TFile fin(model.second.c_str(), "READ");
        TH1D* mc_cross_section = (TH1D*)fin.Get(var_name.c_str());
        assert(mc_cross_section);
        mc_cross_section->SetTitle(description.c_str());
        PlotUtils::MnvH1D* mc_cross_section_mnvh1d =
            new PlotUtils::MnvH1D(*mc_cross_section);
        mc_cross_sections[model.first] = mc_cross_section_mnvh1d;
        fin.Close();
      }
      // Plot
      plot_all_models(plot_info, reco_var->m_hists.m_cross_section,
                      mc_cross_sections);
    }
  }
  fin1.Close();

  //============================================================================
}

#endif  // plotNuisance_C

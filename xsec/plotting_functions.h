#ifndef plotting_functions_h
#define plotting_functions_h

#include <iostream>
#include <sstream>

#include "../includes/myPlotStyle.h"
#include "PlotUtils/MnvColors.h"
#include "PlotUtils/MnvH1D.h"
#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif
#include "PlotUtils/MnvVertErrorBand.h"
#include "SignalDefinition.h"
#include "Systematics.h"  // namespace systematics
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TList.h"
#include "TPad.h"
#include "TPaveStats.h"
//#include "TStyle.h"
#include "TText.h"
#include "Variable.h"

class Variable;

//==============================================================================
// Container class
//==============================================================================
class EventSelectionPlotInfo {
 public:
  // Constructor with Var
  EventSelectionPlotInfo(Variable* variable, float mc_pot, float data_pot,
                         bool do_frac_unc, bool do_cov_area_norm,
                         bool include_stat, SignalDefinition signal_definition)
      : m_mnv_plotter(kCCNuPionIncStyle),
        m_variable(variable),
        m_mc_pot(mc_pot),
        m_data_pot(data_pot),
        m_do_frac_unc(do_frac_unc),
        m_do_cov_area_norm(do_cov_area_norm),
        m_include_stat(include_stat),
        m_signal_definition(signal_definition) {
    m_do_frac_unc_str = m_do_frac_unc ? "Frac" : "Abs";
    m_do_cov_area_norm_str = m_do_cov_area_norm ? "CovAreaNorm" : "";
  }

  // Constructor without var
  EventSelectionPlotInfo(float mc_pot, float data_pot, bool do_frac_unc,
                         bool do_cov_area_norm, bool include_stat,
                         SignalDefinition signal_definition)
      : m_mnv_plotter(kCCNuPionIncStyle),
        m_variable(nullptr),
        m_mc_pot(mc_pot),
        m_data_pot(data_pot),
        m_do_frac_unc(do_frac_unc),
        m_do_cov_area_norm(do_cov_area_norm),
        m_include_stat(include_stat),
        m_signal_definition(signal_definition) {
    m_do_frac_unc_str = m_do_frac_unc ? "Frac" : "Abs";
    m_do_cov_area_norm_str = m_do_cov_area_norm ? "CovAreaNorm" : "";
  }

  // Members
  MnvPlotter m_mnv_plotter;
  Variable* m_variable;
  float m_mc_pot;
  float m_data_pot;
  bool m_do_frac_unc;
  bool m_do_cov_area_norm;
  bool m_include_stat;
  SignalDefinition m_signal_definition;
  std::string m_do_frac_unc_str;
  std::string m_do_cov_area_norm_str;

  // Add X label
  void SetXLabel(PlotUtils::MnvH1D* hist) {
    std::string label =
        m_variable->m_hists.m_xlabel + " (" + m_variable->m_units + ")";
    if (hist) hist->GetXaxis()->SetTitle(label.c_str());
  }

  // Add X label
  void SetXLabel(TH1* hist) {
    std::string label =
        m_variable->m_hists.m_xlabel + " (" + m_variable->m_units + ")";
    if (hist) hist->GetXaxis()->SetTitle(label.c_str());
  }

  // Add title to the MnvPlotter
  void SetTitle() {
    if (!gPad)
      throw std::runtime_error("Need a TCanvas. Please make one first.");
    std::string title = GetSignalName(m_signal_definition);
    m_mnv_plotter.AddHistoTitle(title.c_str());
  }

  void SetTitle(std::string title) {
    if (!gPad)
      throw std::runtime_error("Need a TCanvas. Please make one first.");
    m_mnv_plotter.AddHistoTitle(title.c_str());
  }
};

//==============================================================================
// Some Systematics General Functions
//==============================================================================
void SetErrorGroups(MnvPlotter& mnv_plotter) {
  mnv_plotter.error_summary_group_map.clear();
  mnv_plotter.error_summary_group_map["Flux"].push_back("Flux");
  mnv_plotter.error_summary_group_map["NonResPi"].push_back("GENIE_Rvn1pi");
  mnv_plotter.error_summary_group_map["NonResPi"].push_back("GENIE_Rvp1pi");
  mnv_plotter.error_summary_group_map["NonResPi"].push_back("GENIE_Rvn2pi");
  mnv_plotter.error_summary_group_map["NonResPi"].push_back("GENIE_Rvp2pi");
  mnv_plotter.error_summary_group_map["2p2h"].push_back("Low_Recoil_2p2h_Tune");
  mnv_plotter.error_summary_group_map["LowQ2Pi"].push_back("LowQ2Pi");
  mnv_plotter.error_summary_group_map["Muon"].push_back("Muon_Energy_MINOS");
  mnv_plotter.error_summary_group_map["Muon"].push_back("Muon_Energy_MINERvA");
  mnv_plotter.error_summary_group_map["Muon"].push_back("MINOS_Reconstruction_Efficiency");
  mnv_plotter.error_summary_group_map["Muon"].push_back("MuonAngleXResolution");
  mnv_plotter.error_summary_group_map["Muon"].push_back("MuonAngleYResolution");
  mnv_plotter.error_summary_group_map["Muon"].push_back("MuonResolution");
  mnv_plotter.error_summary_group_map["Michel"].push_back("MichelEfficiency");
  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_D2_MaRES");
  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_EP_MvRES");
  mnv_plotter.error_summary_group_map["Target"].push_back("Target_Mass_CH"); 
  mnv_plotter.error_summary_group_map["Response"].push_back("response_em");
  mnv_plotter.error_summary_group_map["Response"].push_back("response_meson");
  mnv_plotter.error_summary_group_map["Response"].push_back("response_other");
  mnv_plotter.error_summary_group_map["Response"].push_back("response_proton");
  mnv_plotter.error_summary_group_map["Diffractive"].push_back("DiffractiveModelUnc");
  // for(auto g : systematics::kGenieSystematics_FSI)
  //  mnv_plotter.error_summary_group_map["Genie_FSI"].push_back(g);

  for (auto g : systematics::kGenieSystematics_FSI_nucleons)
    mnv_plotter.error_summary_group_map["Genie_FSI_nucleons"].push_back(g);

  for (auto g : systematics::kGenieSystematics_FSI_pions)
    mnv_plotter.error_summary_group_map["Genie_FSI_pions"].push_back(g);

  for (auto g : systematics::kGenieSystematics_InteractionModel)
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(g);

  mnv_plotter.error_summary_group_map["Detector"].push_back("EmuRangeCurve");
  mnv_plotter.error_summary_group_map["Detector"].push_back("Birks");
  mnv_plotter.error_summary_group_map["Detector"].push_back("BetheBloch");
  mnv_plotter.error_summary_group_map["Detector"].push_back("Mass");
  mnv_plotter.error_summary_group_map["Detector"].push_back("PartResp");
  mnv_plotter.error_summary_group_map["Detector"].push_back("TrackAngle");
  mnv_plotter.error_summary_group_map["Detector"].push_back("BeamAngle");
  mnv_plotter.error_summary_group_map["Detector"].push_back("NodeCutEff");
  mnv_plotter.error_summary_group_map["Detector"].push_back("BeamAngleX");
  mnv_plotter.error_summary_group_map["Detector"].push_back("BeamAngleY");

  mnv_plotter.error_summary_group_map["RPA"].push_back("RPA_LowQ2");
  mnv_plotter.error_summary_group_map["RPA"].push_back("RPA_HighQ2");

  //-- define colors of the standard errors
  mnv_plotter.error_color_map.clear();

  /*
    //Systematic color scheme
    mnv_plotter.error_color_map["Flux"]                   = kYellow-3;
    mnv_plotter.error_color_map["Genie_FSI"]              = kGreen+2;
    mnv_plotter.error_color_map["Genie_InteractionModel"] = kPink+2;
    mnv_plotter.error_color_map["Detector"]               = kCyan+2;
    mnv_plotter.error_color_map["RPA"]                    = kOrange+2;
    mnv_plotter.error_color_map["NonResPi"]               = kRed+2;
    mnv_plotter.error_color_map["Low_Recoil_2p2h_Tune"]   = kViolet+2;
  */
}

void Plot_ErrorGroup(EventSelectionPlotInfo p, PlotUtils::MnvH1D* h,
                     std::string error_group_name, std::string tag,
                     double ignore_threshold = 0., double ymax = -1.) {
  TCanvas canvas("c1", "c1");

  p.m_mnv_plotter.good_colors = MnvColors::GetColors(MnvColors::k36Palette);

  // Clone hist
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");

  // X label
  p.SetXLabel(hist);

  // For groups MnvPlotter sets the ymax to headroom*axis_maximum_group
  // Or for summary, ymax = axis_maximum
  if (ymax != -1.) {
    p.m_mnv_plotter.axis_maximum = ymax;
    p.m_mnv_plotter.axis_maximum_group = ymax;
    p.m_mnv_plotter.headroom = 1.;
  }

  const char* legend_position = error_group_name == "" ? "N" : "TR";

  if (error_group_name == "LEGENDONLY") {
    p.m_mnv_plotter.axis_maximum = 1000;
    p.m_mnv_plotter.axis_maximum_group = 1000;
    p.m_mnv_plotter.headroom = 1.;

    p.m_mnv_plotter.DrawErrorSummary(hist, legend_position, p.m_include_stat,
                                     true, ignore_threshold,
                                     p.m_do_cov_area_norm, "", p.m_do_frac_unc);
  } else {
    // XXX WARNING: problems when do_cov_area_norm = true !!
    p.m_mnv_plotter.DrawErrorSummary(
        hist, legend_position, p.m_include_stat, true, ignore_threshold,
        p.m_do_cov_area_norm, error_group_name, p.m_do_frac_unc);
  }

  // Plot Title (For systematics, this has to be after drawing)
  p.SetTitle(tag);

  std::string outfile_name =
      Form("ErrorSummary_%s_%s_%s_%s_%s_%s", tag.c_str(),
           p.m_variable->Name().c_str(), p.m_do_frac_unc_str.c_str(),
           p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(),
           error_group_name.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

// Deprecated?
void Plot_ErrorSummary(EventSelectionPlotInfo p, PlotUtils::MnvH1D* hist,
                       std::string tag) {
  SetErrorGroups(p.m_mnv_plotter);

  Plot_ErrorGroup(p, hist, "", tag.c_str(), 0.0, 0.7);
  Plot_ErrorGroup(p, hist, "Flux", tag.c_str(), 0.0, 0.3);
  Plot_ErrorGroup(p, hist, "Detector", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "Genie_FSI", tag.c_str(), 0.04, 0.15);
  Plot_ErrorGroup(p, hist, "Genie_InteractionModel", tag.c_str(), 0.04, 0.4);
  Plot_ErrorGroup(p, hist, "NonResPi", tag.c_str(), 0.0, 0.1);
  Plot_ErrorGroup(p, hist, "2p2h", tag.c_str(), 0.0, 0.1);
  Plot_ErrorGroup(p, hist, "RPA", tag.c_str(), 0.0, 0.1);
  Plot_ErrorGroup(p, hist, "Michel", tag.c_str(), 0.0, 0.3);
  Plot_ErrorGroup(p, hist, "GENIE", tag.c_str(), 0.0, 0.3);
  Plot_ErrorGroup(p, hist, "Target", tag.c_str(), 0.0, 0.3);
  Plot_ErrorGroup(p, hist, "Response", tag.c_str(), 0.0, 0.3);
}

//==============================================================================
// Event Selection Plots
//==============================================================================
void PlotVar_Selection(EventSelectionPlotInfo p, double ymax = -1.,
                       bool do_log_scale = false, bool do_bg = true,
                       bool do_tuned_bg = false,
                       bool do_bin_width_norm = true) {
  std::cout << "Plotting Selection " << p.m_variable->Name() << std::endl;
  TCanvas canvas("c1", "c1");

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_selection_data);
  assert(p.m_variable->m_hists.m_selection_mc.hist);

  // Get Hists
  PlotUtils::MnvH1D* data = nullptr;
  if (!p.m_variable->m_is_true)
    data = (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_selection_data->Clone(
        "data");
  PlotUtils::MnvH1D* mc =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_selection_mc.hist->Clone(
          "mc");

  // Background
  PlotUtils::MnvH1D* tmp_bg = nullptr;
  if (do_bg) {
    if (do_tuned_bg)
      tmp_bg = (PlotUtils::MnvH1D*)(p.m_variable->m_hists.m_tuned_bg)
                   ->Clone("bg_tmp");
    else
      tmp_bg = (PlotUtils::MnvH1D*)(p.m_variable->m_hists.m_bg.hist)
                   ->Clone("bg_tmp");
  } else
    tmp_bg = NULL;

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }

  // Y-axis limit
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  // p.SetXLabel(p.m_variable->m_hists.m_selection_mc.hist);
  p.SetXLabel(mc);
  // p.SetXLabel(data);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm)
    pot_scale = 1.;
  else
    pot_scale = p.m_data_pot / p.m_mc_pot;

  // Bin Width Normalization
  if (do_bin_width_norm) {
    if (tmp_bg) tmp_bg->Scale(1., "width");
    if (data) data->Scale(1., "width");
    mc->Scale(1., "width");
    // Y label
    std::string yaxis = "N Events / " + p.m_variable->m_units;
    mc->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Draw
  const bool use_hist_titles = false;
  // note: this function applies a POT scale to the bg hist
  p.m_mnv_plotter.DrawDataMCWithErrorBand(
      data, mc, pot_scale, "TR", use_hist_titles, tmp_bg, NULL,
      p.m_do_cov_area_norm, p.m_include_stat);

  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);
  else
    p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);

  // Plot Title
  p.m_mnv_plotter.title_size = 0.05;
  p.SetTitle("Selection " + GetSignalName(p.m_signal_definition));

  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bg_str = do_tuned_bg ? "_tunedBG" : "_untunedBG";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("Selection_%s_%s_%s%s%s%s", p.m_variable->Name().c_str(),
           p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bg_str.c_str(), bwn_str.c_str());
  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotVar_ErrorSummary(EventSelectionPlotInfo p) {
  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_selection_mc.hist);

  SetErrorGroups(p.m_mnv_plotter);

  PlotUtils::MnvH1D* sel =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_selection_mc.hist->Clone(
          uniq());
  Plot_ErrorGroup(p, sel, "", "Sel", 0.0, 0.25);
  Plot_ErrorGroup(p, sel, "LEGENDONLY", "Sel", 0.0, 0.2);
  Plot_ErrorGroup(p, sel, "2p2h", "Sel", 0.0, 0.01);
  Plot_ErrorGroup(p, sel, "Detector", "Sel", 0.0, 0.15);
  Plot_ErrorGroup(p, sel, "Flux", "Sel", 0.0, 0.15);
  Plot_ErrorGroup(p, sel, "Genie_FSI_nucleons", "Sel", 0.004, 0.03);
  Plot_ErrorGroup(p, sel, "Genie_FSI_pions", "Sel", 0.01, 0.1);
  Plot_ErrorGroup(p, sel, "Genie_InteractionModel", "Sel", 0.01, 0.2);
  Plot_ErrorGroup(p, sel, "Muon", "Sel", 0.0, 0.06);
  Plot_ErrorGroup(p, sel, "NonResPi", "Sel", 0.0, 0.06);
  Plot_ErrorGroup(p, sel, "RPA", "Sel", 0.0, 0.01);
  Plot_ErrorGroup(p, sel, "Michel", "Sel", 0.0, 0.15);
  Plot_ErrorGroup(p, sel, "GENIE", "Sel", 0.0, 0.15);
  Plot_ErrorGroup(p, sel, "Target", "Sel", 0.0, 0.15);
  Plot_ErrorGroup(p, sel, "Response", "Sel", 0.0, 0.15);
  Plot_ErrorGroup(p, sel, "Diffractive", "Sel", 0.0, 0.15);
}

//==============================================================================
// Background-Subtracted
//==============================================================================
void Plot_BGSub(EventSelectionPlotInfo p, std::string outdir = ".",
                double ymax = -1, bool do_log_scale = false,
                bool do_bin_width_norm = true) {
  std::cout << "Plotting BG-subtracted Data " << p.m_variable->Name()
            << std::endl;

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_bg_subbed_data);
  assert(p.m_variable->m_hists.m_bg_subbed_data);
  assert(p.m_variable->m_hists.m_effnum.hist);

  TCanvas canvas("c1", "c1");

  // Get Hists
  TH1D* bg_sub_data_w_tot_error =
      new TH1D(p.m_variable->m_hists.m_bg_subbed_data->GetCVHistoWithError());
  TH1D* bg_sub_data_w_stat_error = new TH1D(
      p.m_variable->m_hists.m_bg_subbed_data->GetCVHistoWithStatError());
  TH1D* effnum_w_stat_error =
      new TH1D(p.m_variable->m_hists.m_effnum.hist->GetCVHistoWithStatError());

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }

  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  p.SetXLabel(bg_sub_data_w_tot_error);
  p.SetXLabel(bg_sub_data_w_stat_error);
  p.SetXLabel(effnum_w_stat_error);

  // Overall Normalization
  double pot_scale = -99.;
  // PlotUtils::MnvH1D* tmp_bg = nullptr;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  if (do_bin_width_norm) {
    bg_sub_data_w_tot_error->Scale(1., "width");
    bg_sub_data_w_stat_error->Scale(1., "width");
    effnum_w_stat_error->Scale(1., "width");

    // Y label
    std::string yaxis = "N Events / " + p.m_variable->m_units;
    effnum_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Draw
  const bool use_hist_titles = false;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_tot_error,
                                          effnum_w_stat_error, pot_scale, "TR",
                                          use_hist_titles);

  // p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_stat_error,
  // effnum_w_stat_error, pot_scale, "TR");
  // p.m_mnv_plotter.DrawDataMC(bg_sub_data_w_tot_error, effnum_w_stat_error,
  // pot_scale, "TR");
  // p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_error, effnum,
  // pot_scale, "TR",
  //                                        use_hist_titles, NULL, NULL,
  //                                        p.m_do_cov_area_norm,
  //                                        p.m_include_stat);

  // POT info
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);
  else
    p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);

  // Label re: error bars
  p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.7, 0.75, 0.03);

  // Plot Title
  p.m_mnv_plotter.title_size = 0.04;
  p.SetTitle("Background Subtracted " + GetSignalName(p.m_signal_definition));

  // Print .png
  std::string logy_str = do_log_scale ? "_logscale" : "";
  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("%s/BGSubData_%s_%s_%s%s%s", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());
  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotBGSub_ErrorGroup(EventSelectionPlotInfo p,
                          std::string error_group_name,
                          double ignore_threshold = 0., double ymax = -1.) {
  TCanvas canvas("c1", "c1");

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_bg_subbed_data);

  // Get the histo
  PlotUtils::MnvH1D* bg_sub_data =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg_subbed_data->Clone(
          "bg_sub_data");

  // X label
  p.SetXLabel(bg_sub_data);

  // Y-axis limit
  // For groups MnvPlotter sets the ymax to headroom*axis_maximum_group
  // Or for summary, ymax = axis_maximum
  if (ymax != -1.) {
    p.m_mnv_plotter.axis_maximum = ymax;
    p.m_mnv_plotter.axis_maximum_group = ymax;
    p.m_mnv_plotter.headroom = 1.;
  }

  // Draw
  // XXX WARNING: potential problems when do_cov_area_norm = true
  p.m_mnv_plotter.DrawErrorSummary(bg_sub_data, "TR", p.m_include_stat, true,
                                   ignore_threshold, p.m_do_cov_area_norm,
                                   error_group_name, p.m_do_frac_unc);

  // Plot Title (For systematics, this has to be after drawing)
  p.SetTitle("Background-Subtracted Data");

  // Save to file
  std::string outfile_name = Form(
      "ErrorSummary_BGSubData_%s_%s_%s_%s_%s", p.m_variable->Name().c_str(),
      p.m_do_frac_unc_str.c_str(), p.m_do_cov_area_norm_str.c_str(),
      GetSignalFileTag(p.m_signal_definition).c_str(),
      error_group_name.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotBGSub_ErrorSummary(EventSelectionPlotInfo p) {
  SetErrorGroups(p.m_mnv_plotter);

  PlotUtils::MnvH1D* bg_sub_data =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg_subbed_data->Clone(
          "bg_sub_data");

  // name, ignore threshold, ymax
  // PlotBGSub_ErrorGroup(p, "",                       0.0,  -1.);// // plot all
  // groups together PlotBGSub_ErrorGroup(p, "Flux",                   0.0,
  // 0.1);// PlotBGSub_ErrorGroup(p, "Detector",               0.0,   0.3);//
  // PlotBGSub_ErrorGroup(p, "Genie_FSI",              0.01,  0.1);//
  // PlotBGSub_ErrorGroup(p, "Genie_InteractionModel", 0.02,   0.2);//
  // PlotBGSub_ErrorGroup(p, "NonResPi",               0.0,   0.1);//
  // PlotBGSub_ErrorGroup(p, "2p2h",                   0.0,   0.1);//
  // PlotBGSub_ErrorGroup(p, "RPA",                    0.0,   0.1);//
  Plot_ErrorGroup(p, bg_sub_data, "LEGENDONLY", "BGSub", 0.0);
  Plot_ErrorGroup(p, bg_sub_data, "", "BGSub", 0.0,
                  -1.);  // // plot all groups together
  Plot_ErrorGroup(p, bg_sub_data, "Flux", "BGSub", 0.0, 0.1);                 //
  Plot_ErrorGroup(p, bg_sub_data, "Detector", "BGSub", 0.0, 0.3);             //
  Plot_ErrorGroup(p, bg_sub_data, "Genie_FSI_nucleons", "BGSub", 0.01, 0.1);  //
  Plot_ErrorGroup(p, bg_sub_data, "Genie_FSI_pions", "BGSub", 0.01, 0.1);     //
  Plot_ErrorGroup(p, bg_sub_data, "Genie_InteractionModel", "BGSub", 0.02,
                  0.2);                                            //
  Plot_ErrorGroup(p, bg_sub_data, "NonResPi", "BGSub", 0.0, 0.1);  //
  Plot_ErrorGroup(p, bg_sub_data, "2p2h", "BGSub", 0.0, 0.1);      //
  Plot_ErrorGroup(p, bg_sub_data, "RPA", "BGSub", 0.0, 0.1);       //
  Plot_ErrorGroup(p, bg_sub_data, "Michel", "BGSub", 0.0, 0.3);
  Plot_ErrorGroup(p, bg_sub_data, "GENIE", "BGSub", 0.0, 0.3);
  Plot_ErrorGroup(p, bg_sub_data, "Target", "BGSub", 0.0, 0.3);
  Plot_ErrorGroup(p, bg_sub_data, "Response", "BGSub", 0.0, 0.3);
}

//==============================================================================
// Unfolded
//==============================================================================
void Plot_Unfolded(EventSelectionPlotInfo p, MnvH1D* data, MnvH1D* mc,
                   std::string outdir = ".", double ymax = -1,
                   bool do_log_scale = false, bool do_bin_width_norm = true) {
  std::cout << "Plotting Unfolded " << p.m_variable->Name() << std::endl;

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

  TCanvas canvas("c1", "c1");

  // Get Hists
  PlotUtils::MnvH1D* unfolded = (PlotUtils::MnvH1D*)data->Clone("unfolded");
  PlotUtils::MnvH1D* effnum_true = (PlotUtils::MnvH1D*)mc->Clone("effnum_true");

  TH1D* unfolded_w_tot_error = new TH1D(unfolded->GetCVHistoWithError());
  TH1D* unfolded_w_stat_error = new TH1D(unfolded->GetCVHistoWithStatError());
  TH1D* effnum_true_w_stat_error =
      new TH1D(effnum_true->GetCVHistoWithStatError());

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }

  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  p.SetXLabel(unfolded_w_tot_error);
  p.SetXLabel(unfolded_w_stat_error);
  p.SetXLabel(effnum_true_w_stat_error);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  if (do_bin_width_norm) {
    unfolded_w_tot_error->Scale(1., "width");
    unfolded_w_stat_error->Scale(1., "width");
    effnum_true_w_stat_error->Scale(1., "width");

    // Y label
    std::string yaxis = "N Events / " + p.m_variable->m_units;
    effnum_true_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Draw
  const bool use_hist_titles = false;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(unfolded_w_tot_error,
                                          effnum_true_w_stat_error, pot_scale,
                                          "TR", use_hist_titles);

  // POT info
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);
  else
    p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);

  // Label re: error bars
  p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.7, 0.75, 0.03);

  // Plot Title
  p.m_mnv_plotter.title_size = 0.05;
  p.SetTitle("Unfolded " + GetSignalName(p.m_signal_definition));

  // Print .png
  std::string logy_str = do_log_scale ? "_logscale" : "";
  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("%s/Unfolded_%s_%s_%s%s%s", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotUnfolded_ErrorSummary(EventSelectionPlotInfo p) {
  SetErrorGroups(p.m_mnv_plotter);
  PlotUtils::MnvH1D* unf =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_unfolded->Clone(uniq());
  Plot_ErrorGroup(p, unf, "LEGENDONLY", "Unfolded", 0.0);
  Plot_ErrorGroup(p, unf, "", "Unfolded", 0.0);
  Plot_ErrorGroup(p, unf, "2p2h", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "Detector", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "Flux", "Unfolded", 0.0, 0.12);
  Plot_ErrorGroup(p, unf, "Genie_FSI_nucleons", "Unfolded", 0.01, 0.1);
  Plot_ErrorGroup(p, unf, "Genie_FSI_pions", "Unfolded", 0.01, 0.1);
  Plot_ErrorGroup(p, unf, "Genie_InteractionModel", "Unfolded", 0.01, 0.1);
  Plot_ErrorGroup(p, unf, "NonResPi", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "RPA", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "Michel", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "GENIE", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "Target", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "Response", "Unfolded", 0.0, 0.1);
}

//==============================================================================
// Cross Section
//==============================================================================
void Plot_CrossSection(EventSelectionPlotInfo p, MnvH1D* data, MnvH1D* mc,
                       std::string outdir = ".", double ymax = -1,
                       bool do_log_scale = false,
                       bool do_bin_width_norm = true) {
  std::cout << "Plotting CrossSection " << p.m_variable->Name() << std::endl;

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

  PlotUtils::MnvH1D* data_xsec = (PlotUtils::MnvH1D*)data->Clone("data");
  PlotUtils::MnvH1D* mc_xsec = (PlotUtils::MnvH1D*)mc->Clone("mc");

  TCanvas canvas("c1", "c1");

  // Get Hists
  TH1D* data_xsec_w_tot_error = new TH1D(data_xsec->GetCVHistoWithError());
  TH1D* data_xsec_w_stat_error = new TH1D(data_xsec->GetCVHistoWithStatError());
  TH1D* mc_xsec_w_stat_error = new TH1D(mc_xsec->GetCVHistoWithStatError());

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }

  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // Y-label
  p.m_mnv_plotter.axis_title_offset_y = 1.5;

  // X label
  p.SetXLabel(data_xsec_w_tot_error);
  p.SetXLabel(data_xsec_w_stat_error);
  p.SetXLabel(mc_xsec_w_stat_error);

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
    mc_xsec_w_stat_error->Scale(1.e42, "width");

    // Y label
    // std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel + "
    // (10^{-38} cm^{2}/" + p.m_variable->m_units + "/nucleon)";
    std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel +
                        " (10^{-42} cm^{2}/" + p.m_variable->m_units +
                        "/nucleon)";
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    mc_xsec_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }

  /*
  // Print xsec and error for each bin (AFTER BWN)
  int low_edge = -99;
  int up_edge = -99;
  double val = -99.;
  double err = -99.;
  double frac_err = -99.;
  for (int i = 0; i <= data_xsec_w_tot_error->GetNbinsX(); ++i ){
    low_edge = data_xsec_w_tot_error->GetXaxis()->GetBinLowEdge(i);
    up_edge  = data_xsec_w_tot_error->GetXaxis()->GetBinUpEdge(i);
    val      = data_xsec_w_tot_error->GetBinContent(i);
    err      = data_xsec_w_tot_error->GetBinError(i);
    frac_err = err/val;

    std::cout << i << "  " << low_edge << "  " << up_edge << "  " << val << "  "
  << err << "  " << frac_err << "\n";
  }
  */

  // Draw
  const bool use_hist_titles = false;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(data_xsec_w_tot_error,
                                          mc_xsec_w_stat_error, pot_scale, "TR",
                                          use_hist_titles);

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

  // POT info
  // -1 --> don't do mc POT
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, -1, 0.3, 0.88);
  else
    p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, -1, 0.3, 0.88);

  // Label re: error bars
  p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.7, 0.772, 0.03);

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
}

void PlotCrossSection_ErrorSummary(EventSelectionPlotInfo p) {
  SetErrorGroups(p.m_mnv_plotter);
  PlotUtils::MnvH1D* xsec =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_cross_section->Clone(uniq());

  // for (auto b : xsec->GetErrorBandNames()) std::cout << b << "\n";

  double detector_threshold = 0.0, detector_ymax = .15;
  double FSI_threshold = 0.0, FSI_ymax = 0.1;
  double Int_threshold = 0.015, Int_ymax = 0.15;
  std::string name = p.m_variable->Name();
  if (name == "enu") {
    detector_ymax = 0.2;
    FSI_threshold = 0.0085;
    FSI_ymax = 0.045;
    Int_ymax = 0.08;
  } else if (name == "pmu") {
    detector_ymax = 0.2;
    FSI_threshold = 0.011;
    FSI_ymax = 0.045;
    Int_ymax = 0.08;
  } else if (name == "ptmu") {
    detector_ymax = 0.2;
    FSI_threshold = 0.015;
    FSI_ymax = 0.045;
    Int_threshold = 0.025;
    Int_ymax = 0.2;
  } else if (name == "pzmu") {
    detector_ymax = 0.2;
    FSI_threshold = 0.015;
    FSI_ymax = 0.045;
    Int_ymax = 0.09;
  } else if (name == "q2") {
    detector_ymax = 0.25;
    FSI_threshold = 0.015;
    FSI_ymax = 0.065;
    Int_threshold = 0.025;
    Int_ymax = 0.3;
  } else if (name == "thetamu_deg") {
    detector_ymax = 0.25;
    FSI_threshold = 0.015;
    FSI_ymax = 0.065;
    Int_threshold = 0.025;
    Int_ymax = 0.15;
  } else if (name == "thetapi_deg") {
    FSI_threshold = 0.015;
    FSI_ymax = 0.065;
    Int_threshold = 0.02;
    Int_ymax = 0.12;
  } else if (name == "tpi") {
    detector_ymax = 0.3;
    FSI_threshold = 0.01;
    FSI_ymax = 0.085;
    Int_ymax = 0.12;
  } else if (name == "wexp") {
    detector_ymax = 0.6;
    FSI_threshold = 0.008;
    FSI_ymax = 0.06;
    Int_ymax = 0.1;
  }

  const double ymax = 0.5;
  // name, ignore threshold, ymax
  Plot_ErrorGroup(p, xsec, "LEGENDONLY", "CrossSection", 0.0, 0.1);
  Plot_ErrorGroup(p, xsec, "", "CrossSection", 0.0, 0.3);
  Plot_ErrorGroup(p, xsec, "2p2h", "CrossSection", 0.0, 0.01);
  Plot_ErrorGroup(p, xsec, "Detector", "CrossSection", 0.0, 0.15);
  Plot_ErrorGroup(p, xsec, "Flux", "CrossSection", 0.0, 0.2);
  Plot_ErrorGroup(p, xsec, "Genie_FSI_nucleons", "CrossSection", 4e-3, 0.03);
  Plot_ErrorGroup(p, xsec, "Genie_FSI_pions", "CrossSection", 4e-3, 0.05);
  Plot_ErrorGroup(p, xsec, "Genie_InteractionModel", "CrossSection", 4e-3,
                  0.06);
  Plot_ErrorGroup(p, xsec, "Muon", "CrossSection", 0.0, 0.02);
  Plot_ErrorGroup(p, xsec, "NonResPi", "CrossSection", 0.0, 0.04);
  Plot_ErrorGroup(p, xsec, "RPA", "CrossSection", 0.0, 0.015);
  Plot_ErrorGroup(p, xsec, "Michel", "CrossSection", 0.0, 0.015);
  Plot_ErrorGroup(p, xsec, "GENIE", "CrossSection", 0.0, 0.015);
  Plot_ErrorGroup(p, xsec, "Target", "CrossSection", 0.0, 0.015);
  Plot_ErrorGroup(p, xsec, "Response", "CrossSection", 0.0, 0.015);
}

void PlotMatrix(TMatrixD mtx, std::string name, std::string tag) {
  myPlotStyle();
  setCorrelationPalette();
  TCanvas* c = new TCanvas(uniq(), "scaled matrix");
  mtx.Draw("colz");
  c->Print(Form("%s_%s.png", name.c_str(), tag.c_str()));
  // c->Print(Form("CovMatrix_%s_%s.png", name.c_str(), tag.c_str()));
}

void PrintChi2Info(EventSelectionPlotInfo p, MnvH1D* data, MnvH1D* mc) {
  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

  PlotUtils::MnvH1D* data_xsec = (PlotUtils::MnvH1D*)data->Clone("data");
  PlotUtils::MnvH1D* mc_xsec = (PlotUtils::MnvH1D*)mc->Clone("mc");

  TCanvas canvas("c1", "c1");

  myPlotStyle();
  setCorrelationPalette();

  const bool use_data_error_mtx = true;
  const bool use_only_shape_errors = false;
  bool use_model_stat = false;  // model statistical errors, should be small

  int ndf = -1;
  double pot_scale = 1.;
  double chi2 = p.m_mnv_plotter.Chi2DataMC(
      data_xsec, mc_xsec, ndf, pot_scale, use_data_error_mtx,
      use_only_shape_errors, use_model_stat);
  use_model_stat = true;
  double chi2_mcstat = p.m_mnv_plotter.Chi2DataMC(
      data_xsec, mc_xsec, ndf, pot_scale, use_data_error_mtx,
      use_only_shape_errors, use_model_stat);
  // std::cout << p.m_variable->Name() << "\n";
  // std::cout << "   chi2 = "         << chi2     << "\n";
  // std::cout << "   ndf = "          << ndf      << "\n";
  //"   chi2/ndf = "     <<
  std::cout << p.m_variable->Name() << "  " << chi2 / ndf << "  "
            << chi2_mcstat / ndf << "\n";
}

//==============================================================================
// W Sideband Fit
//==============================================================================
void PlotWSidebandFit_ErrorGroup(EventSelectionPlotInfo p,
                                 std::string error_group_name,
                                 PlotUtils::MnvH1D* h, std::string tag) {
  TCanvas canvas("c1", "c1");

  // Make sure we remembered to load the source histos from the input file.
  assert(h);

  // Get the hist
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");

  // if (error_group_name != "")
  p.m_mnv_plotter.axis_maximum = 0.6;

  // X label
  // hist->GetXaxis()->SetTitle();

  double ignore_threshold;
  if (error_group_name == "")
    ignore_threshold = 0.;
  else
    ignore_threshold = 0.01;

  // XXX WARNING: potential problems when do_cov_area_norm = true
  p.m_mnv_plotter.DrawErrorSummary(hist, "TR", p.m_include_stat, true,
                                   ignore_threshold, p.m_do_cov_area_norm,
                                   error_group_name, p.m_do_frac_unc);

  // Plot Title (For systematics, this has to be after drawing)
  p.SetTitle("W Sideband Fit");

  std::string outfile_name =
      Form("ErrorSummary_%s_%s_%s_%s_%s", tag.c_str(),
           p.m_do_frac_unc_str.c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(),
           error_group_name.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");

  p.m_mnv_plotter.axis_maximum = 0.6;
}

void PlotWSidebandFit_ErrorSummary(EventSelectionPlotInfo p,
                                   PlotUtils::MnvH1D* hist, std::string tag) {
  SetErrorGroups(p.m_mnv_plotter);

  PlotWSidebandFit_ErrorGroup(p, "", hist, tag);  // plot all groups together
  PlotWSidebandFit_ErrorGroup(p, "Flux", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Detector", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Genie_FSI_nucleons", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Genie_FSI_pions", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Genie_InteractionModel", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "NonResPi", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "2p2h", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "RPA", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Michel", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "GENIE", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Target", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Response", hist, tag);
}

void PlotWSidebandStacked(const Variable* variable,
                          const PlotUtils::MnvH1D* h_data,
                          const TObjArray& array_mc, float data_pot,
                          float mc_pot, SignalDefinition signal_definition,
                          std::string tag = "", double ymax = -1,
                          bool do_bin_width_norm = true) {
  // Never don't clone when plotting
  PlotUtils::MnvH1D* data = (PlotUtils::MnvH1D*)h_data->Clone("data");
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  std::cout << "Plotting " << variable->Name() << std::endl;
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  if (ymax > 0) mnvPlotter.axis_maximum = ymax;

  double pot_scale = data_pot / mc_pot;
  std::string label =
      Form("Breakdown_WSideband_%s_%s_PN_%s", variable->m_label.c_str(),
           GetSignalFileTag(signal_definition).c_str(), tag.c_str());
  TCanvas cE("c1", "c1");

  std::string y_label = "Events";
  // bin width norm
  if (do_bin_width_norm) {
    data->Scale(1., "width");
    for (auto h : array)
      dynamic_cast<PlotUtils::MnvH1D*>(h)->Scale(1., "width");
    y_label = "Events / MeV";
  }

  mnvPlotter.DrawDataStackedMC(data, &array, pot_scale, "TR", "Data", -1, -1,
                               1001, variable->m_hists.m_xlabel.c_str(),
                               y_label.c_str());

  double arrow_height = data->GetBinContent(data->GetMaximumBin()) *
                        data->GetNormBinWidth() /
                        data->GetBinWidth(data->GetMaximumBin());
  double arrow_location = (signal_definition == kOnePi) ? 1500 : 1800;
  mnvPlotter.AddCutArrow(arrow_location, 0.0, arrow_height, 200., "R");
  mnvPlotter.WritePreliminary("TL");
  mnvPlotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);
  mnvPlotter.AddHistoTitle(tag.c_str());
  mnvPlotter.MultiPrint(&cE, label, "png");
}

void PlotFittedW(const Variable* variable, const CVUniverse& universe,
                 const PlotUtils::HistWrapper<CVUniverse>& loW_fit,
                 const PlotUtils::HistWrapper<CVUniverse>& midW_fit,
                 const PlotUtils::HistWrapper<CVUniverse>& hiW_fit,
                 float data_pot, float mc_pot,
                 SignalDefinition signal_definition, bool do_prefit = false,
                 std::string tag = "", double ymax = -1,
                 bool do_bin_width_norm = true) {
  // Setup
  std::cout << "Plotting " << variable->Name() << std::endl;
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);

  double pot_scale = data_pot / mc_pot;
  std::string fit_str = do_prefit ? "Pre" : "Post";
  std::string label =
      Form("%sWFit_%s_%s_PN_%s", fit_str.c_str(), variable->m_label.c_str(),
           GetSignalFileTag(signal_definition).c_str(), tag.c_str());
  TCanvas cE("c1", "c1");

  // Never don't clone when plotting
  PlotUtils::MnvH1D* h_data =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_data->Clone("data");
  PlotUtils::MnvH1D* h_sig =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_sig
          .univHist(&universe)
          ->Clone("sig");
  PlotUtils::MnvH1D* h_loW =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_loW
          .univHist(&universe)
          ->Clone("loW");
  PlotUtils::MnvH1D* h_midW =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_midW
          .univHist(&universe)
          ->Clone("midW");
  PlotUtils::MnvH1D* h_hiW =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_hiW
          .univHist(&universe)
          ->Clone("hiW");

  // Apply fit
  if (do_prefit) {
    ;
  } else {
    h_loW->Scale(loW_fit.univHist(&universe)->GetBinContent(1));
    h_midW->Scale(midW_fit.univHist(&universe)->GetBinContent(1));
    h_hiW->Scale(hiW_fit.univHist(&universe)->GetBinContent(1));
  }

  std::string y_label = "Events";

  // bin width norm
  if (do_bin_width_norm) {
    h_data->Scale(1., "width");
    h_sig->Scale(1., "width");
    h_loW->Scale(1., "width");
    h_midW->Scale(1., "width");
    h_hiW->Scale(1., "width");
    y_label = "Events / MeV";
  }

  if (ymax < 0) ymax = h_data->GetBinContent(h_data->GetMaximumBin()) * 1.6;
  if (ymax > 0) mnvPlotter.axis_maximum = ymax;

  // Prepare stack
  std::string legend_name =
      GetTruthClassification_LegendLabel(kWSideband_Signal);
  h_sig->SetTitle(legend_name.c_str());

  legend_name = GetTruthClassification_LegendLabel(kWSideband_Low);
  h_loW->SetTitle(legend_name.c_str());

  legend_name = GetTruthClassification_LegendLabel(kWSideband_Mid);
  h_midW->SetTitle(legend_name.c_str());

  legend_name = GetTruthClassification_LegendLabel(kWSideband_High);
  h_hiW->SetTitle(legend_name.c_str());

  SetHistColorScheme(h_sig, int(kWSideband_Signal),
                     sidebands::kWSideband_ColorScheme);
  SetHistColorScheme(h_loW, int(kWSideband_Low),
                     sidebands::kWSideband_ColorScheme);
  SetHistColorScheme(h_midW, int(kWSideband_Mid),
                     sidebands::kWSideband_ColorScheme);
  SetHistColorScheme(h_hiW, int(kWSideband_High),
                     sidebands::kWSideband_ColorScheme);

  TObjArray* array = new TObjArray();
  array->Add(h_sig);
  array->Add(h_loW);
  array->Add(h_midW);
  array->Add(h_hiW);

  // Draw
  mnvPlotter.DrawDataStackedMC(h_data, array, pot_scale, "TR", "Data", -1, -1,
                               1001, variable->m_hists.m_xlabel.c_str(),
                               y_label.c_str());

  mnvPlotter.WritePreliminary("TL");
  mnvPlotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);

  std::ostringstream oss;
  oss << fit_str << "fit " << tag;
  std::string title = oss.str();

  mnvPlotter.AddHistoTitle(title.c_str());
  mnvPlotter.MultiPrint(&cE, label, "png");

  delete h_data;
  delete array;
}

//==============================================================================
// Backgrounds
//==============================================================================
/*
void PlotBG_ErrorGroup(EventSelectionPlotInfo p, std::string error_group_name,
                       bool do_tuned = false, double ignore_threshold = 0.,
                       double ymax = -1.) {
  TCanvas canvas ("c1","c1");

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_tuned_bg);
  assert(p.m_variable->m_hists.m_bg.hist);

  PlotUtils::MnvH1D* tuned_bg =
(PlotUtils::MnvH1D*)p.m_variable->m_hists.m_tuned_bg->Clone("tuned_bg");
  PlotUtils::MnvH1D* bg       =
(PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg.hist->Clone("bg");

  // X label
  p.SetXLabel(tuned_bg);
  p.SetXLabel(bg);

  // For groups MnvPlotter sets the ymax to headroom*axis_maximum_group
  // Or for summary, ymax = axis_maximum
  if (ymax != -1.) {
    p.m_mnv_plotter.axis_maximum = ymax;
    p.m_mnv_plotter.axis_maximum_group = ymax;
    p.m_mnv_plotter.headroom = 1.;
  }

  // Draw
  // XXX WARNING: potential problems when do_cov_area_norm = true
  if (do_tuned) {
    p.m_mnv_plotter.DrawErrorSummary(tuned_bg, "TR", p.m_include_stat, true,
                                     ignore_threshold, p.m_do_cov_area_norm,
                                     error_group_name, p.m_do_frac_unc);
  }
  else {
    p.m_mnv_plotter.DrawErrorSummary(bg, "TR", p.m_include_stat, true,
                                     ignore_threshold, p.m_do_cov_area_norm,
                                     error_group_name, p.m_do_frac_unc);
  }

  // Plot Title (For systematics, this has to be after drawing)
  if (do_tuned)
    p.SetTitle("Tuned Background");
  else
    p.SetTitle("Untuned Background");


  std::string tuned_str = do_tuned ? "BGTuned" : "BGUntuned";

  std::string outfile_name = Form("ErrorSummary_%s_%s_%s_%s_%s_%s",
                                  tuned_str.c_str(),
                                  p.m_variable->Name().c_str(),
                                  p.m_do_frac_unc_str.c_str(),
                                  p.m_do_cov_area_norm_str.c_str(),
                                  GetSignalFileTag(p.m_signal_definition).c_str(),
                                  error_group_name.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}
*/

void PlotBG_ErrorSummary(EventSelectionPlotInfo p, bool do_tuned = false) {
  SetErrorGroups(p.m_mnv_plotter);
  PlotUtils::MnvH1D* bg;

  if (do_tuned)
    bg =
        (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_tuned_bg->Clone("tuned_bg");
  else
    bg = (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg.hist->Clone("bg");

  double detector_threshold = 0.075, detector_ymax = 0.2;
  double FSI_threshold = 0.01, FSI_ymax = 0.06;
  double Int_threshold = 0.015, Int_ymax = 0.15;
  std::string name = p.m_variable->Name();
  if (name == "enu") {
    Int_threshold = 0.01;
    Int_ymax = 0.12;
  } else if (name == "pmu") {
    Int_threshold = 0.01;
    Int_ymax = 0.12;
  } else if (name == "ptmu") {
    Int_threshold = 0.01;
    Int_ymax = 0.3;
  } else if (name == "pzmu") {
    Int_threshold = 0.01;
    Int_ymax = 0.1;
  } else if (name == "q2") {
    Int_threshold = 0.01;
    Int_ymax = 0.2;
  } else if (name == "thetamu_deg") {
    Int_threshold = 0.01;
    Int_ymax = 0.15;
  } else if (name == "thetapi_deg") {
    Int_threshold = 0.01;
    Int_ymax = 0.12;
  } else if (name == "tpi") {
    Int_threshold = 0.01;
    Int_ymax = 0.12;
  } else if (name == "wexp") {
    Int_threshold = 0.01;
    Int_ymax = 0.2;
  }

  std::string tuned_str = do_tuned ? "BGTuned" : "BGUntuned";

  // name, ignore threshold, ymax
  Plot_ErrorGroup(p, bg, "LEGENDONLY", tuned_str, 0.0);
  Plot_ErrorGroup(p, bg, "", tuned_str, 0.0);
  Plot_ErrorGroup(p, bg, "2p2h", tuned_str, 0.0, 0.05);
  Plot_ErrorGroup(p, bg, "Detector", tuned_str, detector_threshold,
                  detector_ymax);
  Plot_ErrorGroup(p, bg, "Flux", tuned_str, 0.0, 0.15);
  Plot_ErrorGroup(p, bg, "Genie_FSI_pions", tuned_str, FSI_threshold, FSI_ymax);
  Plot_ErrorGroup(p, bg, "Genie_FSI_nucleons", tuned_str, FSI_threshold,
                  FSI_ymax);
  Plot_ErrorGroup(p, bg, "Genie_InteractionModel", tuned_str, Int_threshold,
                  Int_ymax);
  Plot_ErrorGroup(p, bg, "NonResPi", tuned_str, 0.0, 0.05);
  Plot_ErrorGroup(p, bg, "RPA", tuned_str, 0.0, 0.05);
  Plot_ErrorGroup(p, bg, "Michel", tuned_str, 0.0, 0.05);
  Plot_ErrorGroup(p, bg, "GENIE", tuned_str, 0.0, 0.05);
  Plot_ErrorGroup(p, bg, "Target", tuned_str, 0.0, 0.05);
  Plot_ErrorGroup(p, bg, "Response", tuned_str, 0.0, 0.05);
}

/*
void PlotBG_ErrorSummary(EventSelectionPlotInfo p, bool do_tuned = false) {
  SetErrorGroups(p.m_mnv_plotter);

  //name, ignore threshold, ymax
  PlotBG_ErrorGroup(p, "",                       do_tuned, 0.0,  0.6); // plot
all groups together PlotBG_ErrorGroup(p, "Flux",                   do_tuned,
0.0,  0.2); PlotBG_ErrorGroup(p, "Detector",               do_tuned, 0.0,  0.6);
  PlotBG_ErrorGroup(p, "Genie_FSI",              do_tuned, 0.02, 0.6);
  PlotBG_ErrorGroup(p, "Genie_InteractionModel", do_tuned, 0.02, 0.5);
  PlotBG_ErrorGroup(p, "NonResPi",               do_tuned, 0.0,  0.1);
  PlotBG_ErrorGroup(p, "2p2h",                   do_tuned, 0.0,  0.1);
  PlotBG_ErrorGroup(p, "RPA",                    do_tuned, 0.0,  0.1);
  //PlotBG_ErrorGroup(p, "LowQ2Pi");
}
*/

//==============================================================================
// Hack-y functions
//==============================================================================
void PlotTotalError(PlotUtils::MnvH1D* hist, std::string method_str,
                    EventSelectionPlotInfo p) {
  TCanvas cF("c4", "c4");
  TH1D* hTotalErr = (TH1D*)hist
                        ->GetTotalError(p.m_include_stat, p.m_do_frac_unc,
                                        p.m_do_cov_area_norm)
                        .Clone(Form("h_total_err_errSum_%d", __LINE__));
  hTotalErr->SetTitle(Form("Total Uncertainty (%s)", method_str.c_str()));
  if (p.m_do_frac_unc) hTotalErr->GetYaxis()->SetRangeUser(0, 1);
  // hTotalErr->Scale(0.408);
  hTotalErr->Draw();
  cF.Print(Form("TotalUncertainty_%s_%s_%s.png", p.m_do_frac_unc_str.c_str(),
                p.m_do_cov_area_norm_str.c_str(), method_str.c_str()));
}

void PlotStatError(PlotUtils::MnvH1D* hist, EventSelectionPlotInfo p,
                   std::string tag, double ymax = -1.,
                   std::string ylabel = "") {
  TCanvas canvas("c1", "c1");
  double pot_scale = p.m_data_pot / p.m_mc_pot;
  p.SetXLabel(hist);
  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;
  // Y-axis label
  if (ylabel != "") hist->GetYaxis()->SetTitle(ylabel.c_str());

  TH1D* stat_error = (TH1D*)hist->GetStatError(p.m_do_frac_unc).Clone(uniq());

  stat_error->Draw();

  p.m_mnv_plotter.MultiPrint(&canvas, tag.c_str(), "png");
}

void PlotVertBand(std::string band, std::string method_str,
                  PlotUtils::MnvH1D* hist, EventSelectionPlotInfo p) {
  TCanvas cF("c4", "c4");
  TH1* h1 = (TH1*)hist->GetVertErrorBand(band.c_str())
                ->GetErrorBand(p.m_do_frac_unc, p.m_do_cov_area_norm)
                .Clone(Form("Pmu_%s_%s", band.c_str(), method_str.c_str()));
  h1->SetTitle(Form("%s Uncertainty (%s)", band.c_str(), method_str.c_str()));
  p.SetTitle();
  p.SetXLabel(hist);
  // h1->Scale(0.408);
  h1->Draw("h");
  cF.Print(Form("%s_band_%s.png", band.c_str(), method_str.c_str()));
}

void PlotVertBandAllUniverses(std::string band, std::string method_str,
                              PlotUtils::MnvH1D* hist,
                              EventSelectionPlotInfo p) {
  TCanvas cF("c4", "c4");
  p.SetTitle();
  p.SetXLabel(hist);
  hist->GetVertErrorBand(band.c_str())->DrawAll("", true);
  cF.Print(
      Form("%s_band_all_universes_%s.png", band.c_str(), method_str.c_str()));
}

void PlotVertUniverse(std::string band, unsigned int universe,
                      std::string method_str, PlotUtils::MnvH1D* hist) {
  TCanvas cF("c1", "c1");
  TH1D* h1 = hist->GetVertErrorBand(band.c_str())->GetHist(universe);

  h1->SetLineColor(kBlack);
  h1->SetLineStyle(1);
  h1->Draw("hist");
  cF.Print(Form("Enu_%s_band_universe%i_%s.png", band.c_str(), universe + 1,
                method_str.c_str()));
}

void PlotDataMC(PlotUtils::MnvH1D* mc, PlotUtils::MnvH1D* data,
                EventSelectionPlotInfo p, std::string tag) {
  TCanvas canvas("c1", "c1");
  double pot_scale = p.m_data_pot / p.m_mc_pot;
  p.m_mnv_plotter.DrawDataMC(data, mc, pot_scale, "TR");
  p.m_mnv_plotter.MultiPrint(&canvas, tag.c_str(), "png");
}

void PlotDataMCWithError(PlotUtils::MnvH1D* mc, PlotUtils::MnvH1D* data,
                         EventSelectionPlotInfo p, std::string tag) {
  std::cout << "Plotting\n";
  TCanvas canvas("c1", "c1");
  double pot_scale = p.m_data_pot / p.m_mc_pot;
  const bool use_hist_titles = true;
  // PlotUtils::MnvH1D* tmp_bg = nullptr;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(data, mc, pot_scale, "TR");
  //, use_hist_titles,
  // p.m_variable->m_hists.m_bg.hist,
  // tmp_bg, NULL, p.m_do_cov_area_norm, p.m_include_stat);
  p.m_mnv_plotter.MultiPrint(&canvas, tag.c_str(), "png");
}

void PlotTH1_1(TH1* h1, std::string tag, double ymax = -1,
               bool do_log_scale = false, bool do_fit = false) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  // TH1::SetDefaultSumw2();

  TCanvas cF("c4", "c4");

  h1->SetTitle(tag.c_str());
  h1->Draw("HIST");

  if (do_log_scale) cF.SetLogy();

  cF.Update();

  if (ymax > 0) h1->GetYaxis()->SetRangeUser(0, ymax);

  if (do_fit) {
    // gaussian fit
    h1->Fit("gaus", "0");
    h1->Draw("");
    cF.Update();
    TF1* fit1 = (TF1*)h1->GetFunction("gaus");
    fit1->SetLineColor(kRed);
    fit1->Draw("SAME");
  }

  /*
  // draw a line
  TLine *line = new TLine(0,0,0,cF.GetUymax());
  line->SetLineColor(kBlue);
  line->Draw();
  */

  cF.Print(Form("%s.png", tag.c_str()));
  // cF.Print(Form("%s.eps", tag.c_str()));
}

int PlotTogether(TH1* h1, std::string label1, TH1* h2, std::string label2,
                 std::string tag, double ymax = -1, bool do_log_scale = false,
                 bool do_fit = false) {
  std::cout << "PlotTogether" << std::endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  // TH1::SetDefaultSumw2();

  TCanvas cF("c4", "c4");
  cF.Update();

  TLegend* leg = new TLegend(0.75, 0.85, 0.95, 0.95, NULL, "brNDC");

  if (h1->GetMaximum() > h2->GetMaximum()) {
    std::cout << "h1 > h2  " << h1->GetMaximum() << "  " << h2->GetMaximum()
              << "\n";
    h1->SetLineWidth(3);
    h1->SetLineColor(kBlack);
    h1->Draw("HIST");

    cF.Update();

    if (ymax > 0) h1->GetYaxis()->SetRangeUser(0, ymax);

    h2->SetLineColor(kRed);
    h2->SetLineWidth(3);
    h2->Draw("HISTSAME");
  } else {
    std::cout << "h1 < h2  " << h1->GetMaximum() << "  " << h2->GetMaximum()
              << "\n";
    h2->SetLineWidth(3);
    h2->SetLineColor(kRed);
    h2->Draw("HIST");

    cF.Update();

    if (ymax > 0) h2->GetYaxis()->SetRangeUser(0, ymax);

    h1->SetLineColor(kBlack);
    h1->SetLineWidth(3);
    h1->Draw("HISTSAME");
  }

  if (do_log_scale) cF.SetLogy();

  cF.Update();

  //// draw a line
  // TLine *line = new TLine(0,0,0,cF.GetUymax());
  // line->SetLineColor(kBlue);
  // line->Draw();

  // legend
  TLegendEntry* entry = leg->AddEntry("NULL", "", "h");
  leg->SetTextSize(0.035);
  {
    entry = leg->AddEntry("entry", label1.c_str(), "l");
    entry->SetLineWidth(5);
    entry->SetLineColor(kBlack);
  }
  {
    entry = leg->AddEntry("entry", label2.c_str(), "l");
    entry->SetLineWidth(5);
    entry->SetLineColor(kRed);
  }
  leg->Draw();

  cF.Print(Form("%s.png", tag.c_str()));
  // cF.Print(Form("%s.eps", tag.c_str()));

  delete leg;

  return 0;
}

void PlotMC(PlotUtils::MnvH1D* hist, EventSelectionPlotInfo p, std::string tag,
            double ymax = -1., std::string ylabel = "") {
  TCanvas canvas("c1", "c1");
  double pot_scale = p.m_data_pot / p.m_mc_pot;
  p.SetXLabel(hist);
  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;
  // Y-axis label
  if (ylabel != "") hist->GetYaxis()->SetTitle(ylabel.c_str());
  // PlotUtils::MnvH1D* tmp_bg = nullptr;
  p.m_mnv_plotter.DrawMCWithErrorBand(
      hist);  // I think that this call only shows stat errors.
  p.m_mnv_plotter.MultiPrint(&canvas, tag.c_str(), "png");
}

  void PlotRatio(PlotUtils::MnvH1D* num, PlotUtils::MnvH1D* denom, std::string v, double norm, std::string l, bool fixRange) {
    char* vchar = &v[0];
    std::string label(Form("Ratio_%s",vchar));
    char* labchar = &label[0];
    const bool drawSysLines = false;
    const bool drawOneLine  = true;
    double Min = -1., Max = -1.;
    if (fixRange){
      Min = 0.0;
      Max = 1.4;
    }
    const double plotMin = Min;
    const double plotMax = Max;
    const bool covAreaNormalize = false;
    double titleSize = 0.05;
    char* Title = &label[0] ;


    cout << "Plotting ratio " << label << endl;

      TCanvas *c2 = new TCanvas();
      const char* yaxisLabel = "MAD/CCPionInc";
      PlotUtils::MnvPlotter* ratio = new PlotUtils::MnvPlotter();
      ratio->PlotUtils::MnvPlotter::DrawDataMCRatio(num, denom, norm, drawSysLines, drawOneLine, plotMin, plotMax, yaxisLabel, covAreaNormalize);
      ratio->AddHistoTitle(Form("%s %s",Title, l.c_str()), titleSize);
      c2->Print(Form("%s_%s.png",labchar,l.c_str()));

  }


//==============================================================================
// Migration & Efficiency
//==============================================================================
// Helpers
TH2D* GetHistWithUnderOverFlow(TH2D* h) {
  // Create new binning with under/overflow
  UInt_t nx = h->GetXaxis()->GetNbins() + 2;
  Double_t* xbins = new Double_t[nx + 1];
  for (UInt_t i = 0; i < nx; i++)
    xbins[i] = h->GetXaxis()->GetBinLowEdge(i + 0);
  xbins[nx] = xbins[nx - 1] + h->GetXaxis()->GetBinWidth(nx);

  UInt_t ny = h->GetYaxis()->GetNbins() + 2;
  Double_t* ybins = new Double_t[ny + 1];
  for (UInt_t i = 0; i < ny; i++)
    ybins[i] = h->GetYaxis()->GetBinLowEdge(i + 0);
  ybins[ny] = ybins[ny - 1] + h->GetYaxis()->GetBinWidth(ny);

  // Create new histogram with under/overflow
  TH2D* htmp = new TH2D(h->GetName(), h->GetTitle(), nx, xbins, ny, ybins);
  htmp->Sumw2();

  // Fill the new histogram including the overflows
  for (UInt_t i = 1; i <= nx; i++) {
    double xbincenter = htmp->GetXaxis()->GetBinCenter(i);
    for (UInt_t j = 1; j <= ny; j++) {
      double ybincenter = htmp->GetYaxis()->GetBinCenter(j);

      int binnumber = htmp->FindBin(xbincenter, ybincenter);
      int oldbinnumber = h->FindBin(xbincenter, ybincenter);
      // std::cout << binnumber << ":(" << xbincenter << "," << ybincenter << ")
      // ";
      double bincontent = h->GetBinContent(oldbinnumber);
      double binerror = h->GetBinError(oldbinnumber);
      htmp->SetBinContent(binnumber, bincontent);
      htmp->SetBinError(binnumber, binerror);
    }
    // std::cout << "\n";
  }

  // Fill underflow specially
  double xval_firstbin =
      h->GetXaxis()->GetBinLowEdge(1);  // i.e. high edge of underflow
  double yval_firstbin = h->GetYaxis()->GetBinLowEdge(1);

  // Pretty weird way of finding new underflow bin number IMO.
  // I found it on the internet. Might depend on binning or on SetCanExtend?
  double xval_underflow = xval_firstbin - 1.;
  double yval_underflow = yval_firstbin - 1.;

  // New hist's bin number, corresponding to old hist's underflow bin
  int underflowbin = htmp->FindBin(xval_underflow, yval_underflow);

  // Get underflow content and set
  double underflowbincontent = h->GetBinContent(0);
  double underflowbinerror = h->GetBinError(0);
  htmp->SetBinContent(underflowbin, underflowbincontent);
  htmp->SetBinError(underflowbin, underflowbinerror);

  // Restore the number of entries (for when using weights!=1)
  htmp->SetEntries(h->GetEffectiveEntries());
  return htmp;
}

TH2D* RowNormalize(TH2D* h) {
  int first_bin = 0;
  int last_bin_X = h->GetXaxis()->GetNbins() + 1;
  int last_bin_Y = h->GetYaxis()->GetNbins() + 1;

  TH2D* tmp = (TH2D*)h->Clone();
  tmp->Reset();

  for (int y = first_bin; y <= last_bin_Y; ++y) {
    Double_t norm = 0.;
    for (int x = first_bin; x <= last_bin_X; ++x)
      norm += h->GetBinContent(x, y);
    if (fabs(norm) > 1E-8) {
      for (int x = first_bin; x <= last_bin_X; ++x) {
        double percentage = 100 * h->GetBinContent(x, y) / norm;
        tmp->SetBinContent(x, y, percentage);
      }
    }
  }
  return tmp;
}

void PlotMigration_AbsoluteBins(PlotUtils::MnvH2D* hist, std::string name) {
  TCanvas c("c1", "c1");
  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  mnv_plotter.SetRedHeatPalette();
  bool draw_as_matrix = true;
  gStyle->SetHistMinimumZero(kFALSE);
  mnv_plotter.DrawNormalizedMigrationHistogram(hist, draw_as_matrix);
  c.Update();
  c.Print(Form("Migration_AbsBins_%s.png", name.c_str()));
}

void PlotMigration_VariableBins(PlotUtils::MnvH2D* hist, std::string name) {
  TGaxis::SetExponentOffset(-0.035, -0.048, "x");
  TH2D* htmp = GetHistWithUnderOverFlow(hist);
  TH2D* htmp2 = RowNormalize(htmp);
  htmp2->GetXaxis()->SetTitle(Form("%s %s", name.c_str(), "Reco"));
  htmp2->GetYaxis()->SetTitle(Form("%s %s", name.c_str(), "True"));
  TCanvas c("c1", "c1");
  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  mnv_plotter.SetRedHeatPalette();
  bool draw_as_matrix = false;
  gStyle->SetHistMinimumZero(kFALSE);
  mnv_plotter.DrawNormalizedMigrationHistogram(htmp2, draw_as_matrix);
  // gStyle->SetPaintTextFormat("2.0f");
  // htmp2->SetMarkerSize(2);
  // htmp2->Draw("colz text");
  c.Update();
  c.Print(Form("Migration_VarBins_%s.png", name.c_str()));
  // c.SetLogz();
  // c.Update();
  // c.Print("WMigrationMatrix_Wbins_logz.png");
  TGaxis::SetExponentOffset(0, 0, "x");
}

void PlotEfficiency_ErrorSummary(EventSelectionPlotInfo p) {
  SetErrorGroups(p.m_mnv_plotter);
  PlotUtils::MnvH1D* eff =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_efficiency->Clone(uniq());
  Plot_ErrorGroup(p, eff, "", "Eff", 0.0, 0.3);
  Plot_ErrorGroup(p, eff, "Flux", "Eff", 0.0, 0.1);
  Plot_ErrorGroup(p, eff, "Detector", "Eff", 0.0, 0.15);
  Plot_ErrorGroup(p, eff, "Genie_FSI_pions", "Eff", 0.01, 0.2);
  Plot_ErrorGroup(p, eff, "Genie_FSI_nucleons", "Eff", 0.01, 0.2);
  Plot_ErrorGroup(p, eff, "Genie_InteractionModel", "Eff", 0.01, 0.2);
  Plot_ErrorGroup(p, eff, "NonResPi", "Eff", 0.0, 0.1);
  Plot_ErrorGroup(p, eff, "2p2h", "Eff", 0.0, 0.1);
  Plot_ErrorGroup(p, eff, "RPA", "Eff", 0.0, 0.1);
  Plot_ErrorGroup(p, eff, "Michel", "Eff", 0.0, 0.15);
  Plot_ErrorGroup(p, eff, "GENIE", "Eff", 0.0, 0.15);
  Plot_ErrorGroup(p, eff, "Target", "Eff", 0.0, 0.15);
  Plot_ErrorGroup(p, eff, "Response", "Eff", 0.0, 0.15);
}

#endif  // plotting_functions_h

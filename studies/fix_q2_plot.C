// This script was to build/test a function that re-binned an MnvH1D.
// Specifically the q2 plot, to look like Aarons with good xaxis.
// This study was finished and the resulting function is:
// xsec/plotting_functions.h: MnvH1D* RebinQ2Plot(const MnvH1D&).
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "PlotUtils/MnvH1D.h"
#include "TArrayD.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TList.h"
#include "TPad.h"
#include "TStyle.h"
#include "ccpion_common.h"        // GetTestPlaylist
#include "includes/Constants.h"   // CVHW, UniverseMap
#include "includes/Histograms.h"  // LoadHWFromFile
#include "includes/MacroUtil.h"
#include "includes/Plotter.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"

void PlotTH1_1(TH1* h, std::string tag, double ymax = -1,
               bool do_log_scale = false, bool do_fit = false) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TH1* h1 = (TH1*)h->Clone(h->GetName());

  // TH1::SetDefaultSumw2();

  TCanvas cF("c4", "c4");

  h1->SetTitle(tag.c_str());
  h1->Draw("HIST");

  if (do_log_scale) cF.SetLogx();

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

  delete h1;
}

void Plot_CrossSection(Plotter p, MnvH1D* data, MnvH1D* mc,
                       std::string var_name, std::string xlabel,
                       std::string units, std::string outdir = ".",
                       double ymax = -1, bool do_log_scale = false,
                       bool do_bin_width_norm = true) {
  std::cout << "Plotting CrossSection " << var_name << std::endl;

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
  // if (do_log_scale) {
  //  canvas.SetLogy();
  //  p.m_mnv_plotter.axis_minimum = 1;
  //}
  if (do_log_scale) {
    canvas.SetLogx();
  }

  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // Y-label
  p.m_mnv_plotter.axis_title_offset_y = 1.5;

  //// X label
  // p.SetXLabel(data_xsec_w_tot_error);
  // p.SetXLabel(data_xsec_w_stat_error);
  // p.SetXLabel(mc_xsec_w_stat_error);

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
    std::string yaxis =
        "d#sigma/d" + xlabel + " (10^{-42} cm^{2}/" + units + "/nucleon)";
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    mc_xsec_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // // Print xsec and error for each bin (AFTER BWN)
  // int low_edge = -99;
  // int up_edge = -99;
  // double val = -99.;
  // double err = -99.;
  // double frac_err = -99.;
  // for (int i = 0; i <= data_xsec_w_tot_error->GetNbinsX(); ++i ){
  //   low_edge = data_xsec_w_tot_error->GetXaxis()->GetBinLowEdge(i);
  //   up_edge  = data_xsec_w_tot_error->GetXaxis()->GetBinUpEdge(i);
  //   val      = data_xsec_w_tot_error->GetBinContent(i);
  //   err      = data_xsec_w_tot_error->GetBinError(i);
  //   frac_err = err/val;

  //   std::cout << i << "  " << low_edge << "  " << up_edge << "  " << val << "
  //   "
  // << err << "  " << frac_err << "\n";
  // }

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

    //// add label manually
    // if (p.m_variable->Name() == "tpi") ndf = 6;
    // if (p.m_variable->Name() == "enu") ndf = 6;
    // if (p.m_variable->Name() == "pzmu") ndf = 9;
    // if (p.m_variable->Name() == "pmu") ndf = 8;
    // if (p.m_variable->Name() == "wexp") ndf = 4;

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
      Form("%s/CrossSection_%s_%s_%s%s%s", outdir.c_str(), var_name.c_str(),
           p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");

  delete data_xsec;
  delete mc_xsec;
  delete data_xsec_w_tot_error;
  delete data_xsec_w_stat_error;
  delete mc_xsec_w_stat_error;
}

void TH1_rebin_test() {
  TFile fin(
      "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/"
      "DataXSecInputs_2023-02-21.root",
      "READ");

  PlotUtils::MnvH1D* old_hist =
      (PlotUtils::MnvH1D*)fin.Get("selection_data_q2");

  // The current binning scheme in MeV^2 and GeV^2
  //{0, 0.025e6, 0.05e6, 0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6,
  // 0.7e6, 1.0e6, 1.3e6, 2.0e6, 3.0e6}; {0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4,
  // 0.5, 0.7, 1.0, 1.3, 2.0, 3.0};

  // GET CURRENT BIN EDGES
  const TArrayD old_bins_array = *(old_hist->GetXaxis()->GetXbins());
  const int n_old_bins = old_bins_array.GetSize();
  std::cout << "original bins\n";
  for (int i = 0; i < n_old_bins; i++) {
    std::cout << i << "  " << old_bins_array[i] << "  "
              << old_hist->GetBinContent(i) << "\n";
  }
  std::cout << "\n";

  // PLOT THE CURRENT (BAD) SITUATION
  std::string tag = "Broken";
  double ymax = -1;
  bool do_log_scale = true;
  bool do_fit = false;
  PlotTH1_1(old_hist, tag, ymax, do_log_scale, do_fit);

  // MAKE A NEW BIN ARRAY AND CONVERT MEV^2 to GEV^2
  TArrayD new_bins_array = old_bins_array;
  new_bins_array.Set(n_old_bins + 1);  // increase number of bins by 1
  int n_new_bins = new_bins_array.GetSize();
  for (int i = n_new_bins - 1; i >= 2; i--) {
    new_bins_array[i] = new_bins_array[i - 1] * 1.e-6;
  }
  new_bins_array[1] = 6.e-3;  // <-- THIS DETERMINES HOW BIG FIRST BIN APPEARS
  new_bins_array[0] = 0.;

  // MAKE NEW HIST WITH NEW BINS
  std::string new_name = Form("%s_%s", old_hist->GetName(), "1");
  PlotUtils::MnvH1D* new_hist = new PlotUtils::MnvH1D(
      new_name.c_str(), old_hist->GetTitle(), new_bins_array.GetSize() - 1,
      new_bins_array.GetArray());

  // COPY OLD BIN CONTENT TO NEW HIST
  // WARNING: THIS IS A PLOTTING HACK. THIS HIST'S 0TH AND 1ST BINS ARE NO
  // LONGER ACCURATE.
  for (int i = 0; i < n_old_bins; i++) {
    int new_bin_idx = i + 1;
    new_hist->SetBinContent(new_bin_idx, old_hist->GetBinContent(i));
  }
  std::cout << "final bins and bin content\n";
  for (int i = 0; i < n_new_bins; i++) {
    std::cout << i << "  " << new_bins_array[i] << "  "
              << new_hist->GetBinContent(i) << "\n";
  }
  std::cout << "\n";

  // PLOT NEW SITUATION
  tag = "Fixed";
  PlotTH1_1(new_hist, tag, ymax, do_log_scale, do_fit);

  /////////////////////////////////////////////////////////

  bool do_frac_unc = true;
  bool do_cov_area_norm = true;
  bool include_stat = true;
  PlotUtils::MnvH1D* h_data_pot = (PlotUtils::MnvH1D*)fin.Get("data_pot");
  PlotUtils::MnvH1D* h_mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float data_pot = h_data_pot->GetBinContent(1);
  float mc_pot = h_mc_pot->GetBinContent(1);
  Plotter plot_info(mc_pot, data_pot, do_frac_unc, do_cov_area_norm,
                    include_stat, SignalDefinition::OnePi());

  PlotUtils::MnvH1D* h_data = (PlotUtils::MnvH1D*)fin.Get("cross_section_q2");
  PlotUtils::MnvH1D* h_mc = (PlotUtils::MnvH1D*)fin.Get("mc_cross_section_q2");

  // Plot_CrossSection(plot_info, h_data, h_mc, "q2", "q2", "mev2");
  //                     ,std::string outdir = ".", double ymax = -1,
  //                     bool do_log_scale = false,
  //                     bool do_bin_width_norm = true) {

  std::cout << "DONE\n";
}

PlotUtils::MnvH1D* RebinQ2Plot(const PlotUtils::MnvH1D& old_hist) {
  std::cout << "RebinQ2Plot\n";
  // make a new bin array and convert mev^2 TO gev^2
  // old bins
  const TArrayD old_bins_array = *(old_hist.GetXaxis()->GetXbins());
  const int n_old_bins = old_bins_array.GetSize();

  std::cout << "size of old hist bin array " << n_old_bins << "\n";

  TArrayD new_bins_array = old_bins_array;

  // increase number of bins by 1.
  // This "Set" adds a new 0 at the beginning of the array.
  new_bins_array.Set(n_old_bins + 1);
  const int n_new_bins = new_bins_array.GetSize();
  std::cout << "size of new hist bin array " << n_new_bins << "\n";

  // Scale to GeV^2
  for (int i = n_new_bins - 1; i >= 2; i--) {
    new_bins_array[i] = new_bins_array[i - 1] * 1.e-6;
  }

  // Setting the first bin edge to 0.006 makes the q2 plot look good in GeV^2
  // and log scale. Also, it's what Aaron uses.
  new_bins_array[1] = 6.e-3;  // <-- THIS DETERMINES HOW BIG FIRST BIN APPEARS
  new_bins_array[0] = 0.;

  // manually make the new MH1D
  std::string new_name = Form("%s_%s", old_hist.GetName(), "_rebin");
  PlotUtils::MnvH1D* new_hist = new PlotUtils::MnvH1D(
      new_name.c_str(), old_hist.GetTitle(), new_bins_array.GetSize() - 1,
      new_bins_array.GetArray());

  new_hist->SetLineColor(old_hist.GetLineColor());
  new_hist->SetLineStyle(old_hist.GetLineStyle());
  new_hist->SetLineWidth(old_hist.GetLineWidth());

  new_hist->SetMarkerColor(old_hist.GetMarkerColor());
  new_hist->SetMarkerStyle(old_hist.GetMarkerStyle());
  new_hist->SetMarkerSize(old_hist.GetMarkerSize());

  new_hist->SetTitle(old_hist.GetTitle());
  new_hist->GetXaxis()->SetTitle(old_hist.GetXaxis()->GetTitle());
  new_hist->GetYaxis()->SetTitle(old_hist.GetYaxis()->GetTitle());

  // finally, move contents, bin-by-bin, universe-by-universe from old to new
  // WARNING: THIS IS A PLOTTING HACK. THIS HIST'S 0TH AND 1ST BINS ARE NO
  // TECHNICALLY CORRECY
  // CV
  for (int i = 0; i < n_old_bins; i++) {
    int new_bin_idx = i + 1;
    new_hist->SetBinContent(new_bin_idx, old_hist.GetBinContent(i));
    new_hist->SetBinError(new_bin_idx, old_hist.GetBinError(i));
  }
  new_hist->SetBinContent(0, 0.);
  new_hist->SetBinError(0, 0.);

  // ASSERT CV
  for (int i = 0; i < n_old_bins; i++) {
    int new_bin_idx = i + 1;
    assert(new_hist->GetBinContent(new_bin_idx) == old_hist.GetBinContent(i));
    assert(new_hist->GetBinError(new_bin_idx) == old_hist.GetBinError(i));
  }

  // Universes
  for (auto error_name : old_hist.GetVertErrorBandNames()) {
    int n_univs = old_hist.GetVertErrorBand(error_name)->GetNHists();
    new_hist->AddVertErrorBand(error_name, n_univs);

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      TH1* univ_i_hist_new =
          new_hist->GetVertErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetVertErrorBand(error_name)->GetHist(univ_i));

      for (int i = 0; i < n_old_bins; i++) {
        int new_bin_idx = i + 1;
        univ_i_hist_new->SetBinContent(new_bin_idx,
                                       univ_i_hist_old->GetBinContent(i));
        univ_i_hist_new->SetBinError(new_bin_idx,
                                     univ_i_hist_old->GetBinError(i));
      }

      univ_i_hist_new->SetBinContent(0, 0.);
      univ_i_hist_new->SetBinError(0, 0.);
      delete univ_i_hist_old;
    }
  }

  // ASSERT UNIVERSES
  for (auto error_name : old_hist.GetVertErrorBandNames()) {
    int n_univs = old_hist.GetVertErrorBand(error_name)->GetNHists();
    // std::cout << error_name << "\n";

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      // std::cout << "  " << univ_i << "\n";

      TH1* univ_i_hist_new =
          new_hist->GetVertErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetVertErrorBand(error_name)->GetHist(univ_i));

      for (int i = 0; i < n_old_bins; i++) {
        int new_bin_idx = i + 1;
        // std::cout << "    " << univ_i_hist_new->GetBinContent(new_bin_idx) <<
        // " = " << univ_i_hist_old->GetBinContent(i) <<  " | "
        //          << univ_i_hist_new->GetBinError(new_bin_idx) << " = " <<
        //          univ_i_hist_old->GetBinError(i) << "\n";
        assert(univ_i_hist_new->GetBinContent(new_bin_idx) ==
               univ_i_hist_old->GetBinContent(i));
        assert(univ_i_hist_new->GetBinError(new_bin_idx) ==
               univ_i_hist_old->GetBinError(i));
      }
      delete univ_i_hist_old;
    }
  }

  return new_hist;
}

/*
void MH1_rebin_test() {
  // I/O
  TFile
fin("/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/DataXSecInputs_2023-02-21.root",
"READ"); std::cout << "Reading input from " << fin.GetName() << endl;

  // INPUT TUPLES
  // Don't actually use the MC chain, only load it to indirectly access its
  // systematics
  //std::string data_file_list = GetPlaylistFile(plist, false);
  //std::string mc_file_list = GetPlaylistFile("ME1A", true);
  std::string data_file_list = GetTestPlaylist(false);
  std::string mc_file_list = GetTestPlaylist(true);

  // Macro Utility
  const std::string macro("fix_q2_plot");
  bool do_truth = false, is_grid = false, do_systematics = true;
  int signal_definition_int = 0;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       "ME1A", do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);


  // load all histograms and HW and MH1s from input file into the q2 variable
  Variable* q2 = new Variable("q2", "Q^{2}", "MeV^{2}", CCPi::GetBinning("q2"),
&CVUniverse::GetQ2);
  //std::vector<Variable*> variables = {q2};
  //for (auto v : variables) {
  std::cout << "Loading hists for variable " << q2->Name() << "\n";
  q2->LoadMCHistsFromFile(fin, util.m_error_bands);
  q2->InitializeDataHists();
  //}

  // var->m_hists.m_cross_section
  // var->m_hists.m_mc_cross_section


  const UniverseMap error_bands = util.m_error_bands;
  // universe loop
  for (auto error_band : error_bands) {
    std::vector<CVUniverse*> universes = error_band.second;
    for (auto universe : universes) {
      TH1* hist_u = q2->m_hists.m_mc_cross_section.univHist(universe);
      std::cout << universe->ShortName() << "  " << hist_u->GetName() << "  " <<
hist_u->GetTitle() << "\n";
    }
  }
}
*/

void fix_q2_plot() {
  TFile fin(
      "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/"
      "DataXSecInputs_2023-02-21.root",
      "READ");
  std::cout << "Reading input from " << fin.GetName() << endl;
  PlotUtils::MnvH1D* old_hist_data =
      (PlotUtils::MnvH1D*)fin.Get("cross_section_q2");
  PlotUtils::MnvH1D* old_hist_mc =
      (PlotUtils::MnvH1D*)fin.Get("mc_cross_section_q2");
  assert(old_hist_data);
  assert(old_hist_mc);
  PlotUtils::MnvH1D* new_hist_data = RebinQ2Plot(*old_hist_data);
  PlotUtils::MnvH1D* new_hist_mc = RebinQ2Plot(*old_hist_mc);

  bool do_frac_unc = true;
  bool do_cov_area_norm = true;
  bool include_stat = true;
  PlotUtils::MnvH1D* h_data_pot = (PlotUtils::MnvH1D*)fin.Get("data_pot");
  PlotUtils::MnvH1D* h_mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float data_pot = h_data_pot->GetBinContent(1);
  float mc_pot = h_mc_pot->GetBinContent(1);
  Plotter plot_info(mc_pot, data_pot, do_frac_unc, do_cov_area_norm,
                    include_stat, SignalDefinition::OnePi());

  PlotUtils::MnvH1D* h_data = (PlotUtils::MnvH1D*)fin.Get("cross_section_q2");
  PlotUtils::MnvH1D* h_mc = (PlotUtils::MnvH1D*)fin.Get("mc_cross_section_q2");

  Plot_CrossSection(plot_info, new_hist_data, new_hist_mc, "q2", "q2", "mev2",
                    ".", -1, true);
  //                     ,std::string outdir = ".", double ymax = -1,
  //                     bool do_log_scale = false,
  //                     bool do_bin_width_norm = true) {

  delete old_hist_data;
  delete old_hist_mc;
  delete new_hist_data;
  delete new_hist_mc;
}

#ifndef studies_plotting_functions_h
#define studies_plotting_functions_h

#include <iostream>
#include <sstream>

#include "../includes/myPlotStyle.h"
#include "PlotUtils/MnvColors.h"
#include "PlotUtils/MnvH1D.h"
#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif

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
//#include "myPlotStyle.h"

#include "Constants.h" // enum SignalDefinition
#include "SignalDefinition.h" // GetSignalFileTag
#include "Variable.h"
#include "Plotter.h"

void PlotCutVar(const Variable* variable, const PlotUtils::MnvH1D* h_data,
                const TObjArray& array_mc, float data_pot, float mc_pot,
                SignalDefinition signal_definition, std::string title = "",
                std::string truth_cat = "", double ymax = -1,
                bool do_bin_width_norm = true) {
  // Never don't clone when plotting
  PlotUtils::MnvH1D* data = (PlotUtils::MnvH1D*)h_data->Clone("data");
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  if (ymax > 0) mnv_plotter.axis_maximum = ymax;
  if (truth_cat == "Hadrons") mnv_plotter.legend_offset_x = .11;
  if (truth_cat == "FSPart") mnv_plotter.legend_offset_x = .15;
  double pot_scale = data_pot / mc_pot;
  std::string label =
      Form("Breakdown_CutVar_%s_%s_PN_%s", variable->m_label.c_str(),
           GetSignalFileTag(signal_definition).c_str(), truth_cat.c_str());
  TCanvas cE("c1", "c1");

  std::string y_label = "Events";
  // bin width norm
  if (do_bin_width_norm) {
    data->Scale(1., "width");
    for (auto h : array)
      dynamic_cast<PlotUtils::MnvH1D*>(h)->Scale(1., "width");
    y_label = "Events / MeV";
  }

  mnv_plotter.DrawDataStackedMC(data, &array, pot_scale, "TR", "Data", -1, -1,
                                1001, variable->m_hists.m_xlabel.c_str(),
                                y_label.c_str());
  mnv_plotter.WritePreliminary("TL");
  mnv_plotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);
  mnv_plotter.AddHistoTitle(title.c_str());
  mnv_plotter.MultiPrint(&cE, label, "png");
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
  double arrow_location = signal_definition.m_w_max + 100.;
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

void PlotBackground(Variable* variable, const PlotUtils::MnvH1D* h_data,
                    const TObjArray& array_mc, float data_pot, float mc_pot,
                    SignalDefinition signal_definition, std::string tag = "",
                    double ymax = -1, bool draw_arrow = false,
                    bool do_bin_width_norm = true) {
  // Never don't clone when plotting
  PlotUtils::MnvH1D* data = (PlotUtils::MnvH1D*)h_data->Clone("data");
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  if (ymax > 0) mnvPlotter.axis_maximum = ymax;
  if (tag == "FSP") mnvPlotter.legend_offset_x = .15;
  if (tag == "Hadrons") mnvPlotter.legend_offset_x = .06;
  if (tag == "Wexp") mnvPlotter.legend_offset_x = .06;
  mnvPlotter.legend_text_size = 0.0405;

  double pot_scale = data_pot / mc_pot;
  std::string label =
      Form("Breakdown_BG_%s_%s_%s", GetSignalFileTag(signal_definition).c_str(),
           variable->m_label.c_str(), tag.c_str());

  std::string y_label = "Events";

  // bin width norm
  if (do_bin_width_norm) {
    data->Scale(1., "width");
    for (auto h : array)
      dynamic_cast<PlotUtils::MnvH1D*>(h)->Scale(1., "width");
    y_label = "Events / MeV";
  }

  TCanvas cE("c1", "c1");
  mnvPlotter.DrawDataStackedMC(data, &array, pot_scale, "TR", "Data", -1, -1,
                               1001, variable->m_hists.m_xlabel.c_str(),
                               y_label.c_str());

  if (draw_arrow) {
    // double arrow_height = data->GetBinContent(data->GetMaximumBin()) *
    //                      data->GetNormBinWidth()/data->GetBinWidth(data->GetMaximumBin());
    double arrow_height = 250;
    double arrow_location = signal_definition.m_w_max;
    mnvPlotter.AddCutArrow(arrow_location, 0.0, arrow_height, 200., "L");
  }

  mnvPlotter.WritePreliminary("TL");
  mnvPlotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);
  mnvPlotter.AddHistoTitle("Background");
  mnvPlotter.MultiPrint(&cE, label, "png");
  // mnvPlotter.MultiPrint(&cE, label, "eps");
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

void PlotMigration_AbsoluteBins(PlotUtils::MnvH2D* hist, std::string name,
                                double zmax = -1) {
  TCanvas c("c1", "c1");
  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  mnv_plotter.SetRedHeatPalette();
  bool draw_as_matrix = true;
  gStyle->SetHistMinimumZero(kFALSE);
  PlotUtils::MnvH2D* h = (PlotUtils::MnvH2D*)hist->Clone("h");
  if (zmax > 0) h->SetMaximum(zmax);
  mnv_plotter.DrawNormalizedMigrationHistogram(h, draw_as_matrix, false, true,
                                               true);
  c.Update();
  c.Print(Form("Migration_AbsBins_%s.png", name.c_str()));
}

void PlotMigration_VariableBins(PlotUtils::MnvH2D* hist, std::string name,
                                double zmax = -1) {
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
  if (zmax > 0) htmp2->SetMaximum(zmax);
  mnv_plotter.DrawNormalizedMigrationHistogram(htmp2, draw_as_matrix, false,
                                               true, true);
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

PlotUtils::MnvH1D* RebinQ2Plot(const PlotUtils::MnvH1D& old_hist) {
  // make a new bin array and convert mev^2 TO gev^2
  // old bins
  const TArrayD old_bins_array = *(old_hist.GetXaxis()->GetXbins());
  const int n_old_bins = old_bins_array.GetSize();

  // std::cout << "size of old hist bin array " << n_old_bins << "\n";

  TArrayD new_bins_array = old_bins_array;

  // increase number of bins by 1.
  // This "Set" adds a new 0 at the beginning of the array.
  new_bins_array.Set(n_old_bins + 1);
  const int n_new_bins = new_bins_array.GetSize();
  // std::cout << "size of new hist bin array " << n_new_bins << "\n";

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

void PlotVar_Selection(Plotter p, double ymax = -1., bool do_log_scale = false,
                       bool do_bg = true, bool do_tuned_bg = false,
                       bool do_bin_width_norm = true) {
  std::cout << "Plotting Selection " << p.m_variable->Name() << std::endl;
  TCanvas canvas("c1", "c1");

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_selection_data);
  assert(p.m_variable->m_hists.m_selection_mc.hist);

  // Get Hists
  // Selection
  PlotUtils::MnvH1D* mc =
      p.m_variable->Name() == "q2"
          ? RebinQ2Plot(*p.m_variable->m_hists.m_selection_mc.hist)
          : (PlotUtils::MnvH1D*)
                p.m_variable->m_hists.m_selection_mc.hist->Clone("mc");

  PlotUtils::MnvH1D* data = nullptr;
  if (!p.m_variable->m_is_true) {
    data =
        p.m_variable->Name() == "q2"
            ? RebinQ2Plot(*p.m_variable->m_hists.m_selection_data)
            : (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_selection_data->Clone(
                  "data");
  }

  // Background
  PlotUtils::MnvH1D* tmp_bg = nullptr;
  if (do_bg) {
    if (do_tuned_bg) {
      tmp_bg =
          p.m_variable->Name() == "q2"
              ? RebinQ2Plot(*p.m_variable->m_hists.m_tuned_bg)
              : (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_tuned_bg->Clone(
                    "bg_tmp");
    } else {
      tmp_bg = p.m_variable->Name() == "q2"
                   ? RebinQ2Plot(*p.m_variable->m_hists.m_bg.hist)
                   : (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg.hist->Clone(
                         "bg_tmp");
    }
  }

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
  if (p.m_variable->Name() == "q2") {
    canvas.SetLogx();
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

#endif  // studies_plotting_functions_h

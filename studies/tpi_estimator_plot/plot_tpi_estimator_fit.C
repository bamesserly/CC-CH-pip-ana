#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "../includes/Plotter.h"
#include "MnvPlotter.h"
#include "TStyle.h"

//  KEY: TGraph	TpiEst;1	Graph
//  KEY: TGraphErrors	TpiEst_Errs;1	Graph
//  KEY: TGraphErrors	TpiEst_Errs_Fit;1	Graph
//  KEY: TGraph	FitCurve;1	Graph

/*
int PlotTogether(TH1* h1, std::string label1, TH1* h2, std::string label2,
                 std::string tag, double ymax = -1, bool do_log_scale = false,
                 bool do_fit = false, std::string ylabel = "") {
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
    h2->GetYaxis()->SetTitle(ylabel.c_str());
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
    h1->GetYaxis()->SetTitle(ylabel.c_str());
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
*/

int plot_tpi_estimator_fit() {
  // Infiles
  TFile fin("2024-09-25_tpi_est_plot/tpi_from_range_fit.root","READ");
  std::cout << "Reading input from " << fin.GetName() << std::endl;

  TGraphErrors* data = (TGraphErrors*)fin.Get("TpiEst_Errs");
  TGraphErrors* fit = (TGraphErrors*)fin.Get("TpiEst_Errs_Fit");

  PlotUtils::MnvPlotter mnv_plotter(kCCNuPionIncStyle);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TCanvas cF("c4", "c4");
  cF.Update();

  data->Draw("AP");
  fit->Draw("PSAME");

  cF.Print("tpi_estimator_fit.png");

  //PlotUtils::MnvH1D* AaronXsecGENIE = (PlotUtils::MnvH1D*)fin1.Get("q2_xsec");
  //Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
  //                  do_cov_area_norm, include_stat,
  //                  util.m_signal_definition);
  //Plot_CrossSection(plot_info, reco_var->m_hists.m_cross_section,
  //                  m_mc_cross_section, ".", -1, false, do_bin_width_norm);

  return 0;
}

#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <tuple>
#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TList.h"
#include "TPad.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TArrayD.h"
#include "PlotUtils/MnvH1D.h"

void PlotTH1_1(TH1* h, std::string tag, double ymax = -1, bool do_log_scale = false, bool do_fit = false) {
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

void fix_q2_plot(){
  TFile fin("/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/DataXSecInputs_2023-02-21.root", "READ");

  PlotUtils::MnvH1D* old_hist = (PlotUtils::MnvH1D*)fin.Get("selection_data_q2");

  // The current binning scheme in MeV^2 and GeV^2
  //{0, 0.025e6, 0.05e6, 0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6, 0.7e6, 1.0e6, 1.3e6, 2.0e6, 3.0e6};
  //{0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.3, 2.0, 3.0};

  // GET CURRENT BIN EDGES
  const TArrayD old_bins_array = *(old_hist->GetXaxis()->GetXbins());
  const int n_old_bins = old_bins_array.GetSize();
  std::cout << "original bins\n";
  for (int i = 0; i < n_old_bins; i++) {
    std::cout << i << "  " << old_bins_array[i] << "  " << old_hist->GetBinContent(i) <<"\n";
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
  new_bins_array.Set(n_old_bins+1); // increase number of bins by 1
  int n_new_bins = new_bins_array.GetSize();
  for (int i = n_new_bins-1; i >= 2; i--) {
    new_bins_array[i] = new_bins_array[i-1]*1.e-6;
  }
  new_bins_array[1] = 6.e-3; // <-- THIS DETERMINES HOW BIG FIRST BIN APPEARS
  new_bins_array[0] = 0.;


  // MAKE NEW HIST WITH NEW BINS
  std::string new_name = Form("%s_%s",old_hist->GetName(), "1");
  PlotUtils::MnvH1D* new_hist = new PlotUtils::MnvH1D(new_name.c_str(), old_hist->GetTitle(), 
                                                      new_bins_array.GetSize()-1, new_bins_array.GetArray() );

  // COPY OLD BIN CONTENT TO NEW HIST
  // WARNING: THIS IS A PLOTTING HACK. THIS HIST'S 0TH AND 1ST BINS ARE NO
  // LONGER ACCURATE.
  for (int i = 0; i < n_old_bins; i++) {
    int new_bin_idx = i+1;
    new_hist->SetBinContent(new_bin_idx, old_hist->GetBinContent(i));
  }
  std:: cout << "final bins and bin content\n";
  for (int i = 0; i < n_new_bins; i++) {
    std::cout << i << "  " << new_bins_array[i] << "  " << new_hist->GetBinContent(i) << "\n";
  }
  std::cout << "\n";
  

  // PLOT NEW SITUATION
  tag = "Fixed";
  PlotTH1_1(new_hist, tag, ymax, do_log_scale, do_fit);

  std::cout << "DONE\n";
}

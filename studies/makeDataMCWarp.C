#ifndef makeDataMCWarp_C
#define makeDataMCWarp_C

#include "PlotUtils/MnvH1D.h"
#include "TFile.h"
#include "TH1D.h"

void makeDataMCWarp() {
  // Get In File
  TFile infile("rootfiles/MCXSecInputs_20190801.root");

  // Get POT

  // Get Data Hist

  // Get MC Hist

  // POT scale
  PlotUtils::MnvH1D* h_mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  PlotUtils::MnvH1D* h_data_pot = (PlotUtils::MnvH1D*)fin.Get("data_pot");
  float mc_pot = h_mc_pot->GetBinContent(1);
  float data_pot = h_data_pot->GetBinContent(1);

  // Make Ratio
  TH1* h_ratio = h_data->Clone("data_mc_ratio");
  h_ratio->Divide(h_mc);

  // Make Fit
  int polynomial_order = 1;
  TF1* fit_function;
  if (polynomial_order == 1) {
    fit_function = new TF1("pol1fit", "pol1", 0, 30000);

    h_ratio->Fit(fit_function, "REQ");
    double firstord = fit_function->GetParameter(1);
    double zerothord = fit_function->GetParameter(0);

    std::cout << "a1 = " << firstord << " a0 = " << zerothord << "\n"
    // double binuncert = ((secondord)*(find_val*find_val)) +
    //                   ((firstord)*(find_val)) + (zerothord);
    // if ((uncertainty == "bafflescraping") || (uncertainty == "oldhorn1") ||
    // (uncertainty == "horncurrent") )
    //    binuncert = binuncert*scale;
    // hfluxuncertainty->SetBinContent(k,binuncert);
    // hchisquare->SetBinContent(k,fit_function->GetChisquare());
  }
}

#endif  // makeDataMCWarp_C

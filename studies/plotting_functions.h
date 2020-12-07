#ifndef plotting_functions_h
#define plotting_functions_h

#include <iostream>
#include <sstream>

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TText.h"

#include "PlotUtils/MnvH1D.h"


#ifndef __CINT__ // Hide class Variable because it has c++11 stuff
#include "PlotUtils/MnvPlotter.h"
#include "../util/plot/myPlotStyle.h"
#endif

#include "SignalDefinition.h"
#include "Variable.h"

void PlotCutVar(const Variable* variable, const PlotUtils::MnvH1D* h_data, 
                const TObjArray& array_mc, 
                float data_pot, float mc_pot, 
                SignalDefinition signal_definition, 
                std::string title = "",
                std::string truth_cat = "",
                double ymax = -1, bool do_bin_width_norm = true) {

  // Never don't clone when plotting
  PlotUtils::MnvH1D* data = (PlotUtils::MnvH1D*)h_data->Clone("data");
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  if (ymax >0) mnv_plotter.axis_maximum = ymax;
  if (truth_cat == "Hadrons") mnv_plotter.legend_offset_x = .11;
  if (truth_cat == "FSPart")  mnv_plotter.legend_offset_x = .15;
  double pot_scale = data_pot/mc_pot;
  std::string label = Form("Breakdown_CutVar_%s_%s_PN_%s", variable->m_label.c_str(), 
                           GetSignalFileTag(signal_definition).c_str(), truth_cat.c_str());
  TCanvas cE ("c1","c1"); 

  std::string y_label = "Events";
  // bin width norm
  if ( do_bin_width_norm ) {
    data->Scale(1., "width");
    for (auto h : array)
      dynamic_cast<PlotUtils::MnvH1D*>(h)->Scale(1., "width");
    y_label = "Events / MeV";
  }

  mnv_plotter.DrawDataStackedMC(data, &array, pot_scale, "TR",  "Data", -1, -1, 
                               1001, variable->m_hists.m_xlabel.c_str(),
                               y_label.c_str());
  mnv_plotter.WritePreliminary("TL");
  mnv_plotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);
  mnv_plotter.AddHistoTitle(title.c_str());
  mnv_plotter.MultiPrint(&cE, label, "png");
}


void PlotWSidebandStacked(const Variable* variable, const PlotUtils::MnvH1D* h_data,
                          const TObjArray& array_mc,
                          float data_pot, float mc_pot,
                          SignalDefinition signal_definition,
                          std::string tag = "", double ymax = -1,
                          bool do_bin_width_norm = true) {

  // Never don't clone when plotting
  PlotUtils::MnvH1D* data = (PlotUtils::MnvH1D*)h_data->Clone("data");
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  std::cout << "Plotting " << variable->Name() << std::endl;
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  if (ymax >0) mnvPlotter.axis_maximum = ymax;

  double pot_scale = data_pot/mc_pot;
  std::string label = Form("Breakdown_WSideband_%s_%s_PN_%s", variable->m_label.c_str(), 
                           GetSignalFileTag(signal_definition).c_str(), tag.c_str());
  TCanvas cE ("c1","c1"); 

  std::string y_label = "Events";
  // bin width norm
  if ( do_bin_width_norm ) {
    data->Scale(1., "width");
    for (auto h : array)
      dynamic_cast<PlotUtils::MnvH1D*>(h)->Scale(1., "width");
    y_label = "Events / MeV";
  }

  mnvPlotter.DrawDataStackedMC(data, &array, pot_scale, "TR",  "Data", -1, -1, 
                               1001, variable->m_hists.m_xlabel.c_str(), y_label.c_str());

  double arrow_height = data->GetBinContent(data->GetMaximumBin()) *
                        data->GetNormBinWidth()/data->GetBinWidth(data->GetMaximumBin());
  double arrow_location = (signal_definition == kOnePi) ? 1500 : 1800;
  mnvPlotter.AddCutArrow(arrow_location, 0.0, arrow_height, 200., "R");
  mnvPlotter.WritePreliminary("TL");
  mnvPlotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);
  mnvPlotter.AddHistoTitle(tag.c_str());
  mnvPlotter.MultiPrint(&cE, label, "png");
}


void PlotFittedW(const Variable* variable, 
                 const CVUniverse& universe,
                 const PlotUtils::HistWrapper<CVUniverse>& loW_fit, 
                 const PlotUtils::HistWrapper<CVUniverse>& midW_fit, 
                 const PlotUtils::HistWrapper<CVUniverse>& hiW_fit, 
                 float data_pot, float mc_pot, 
                 SignalDefinition signal_definition, 
                 bool do_prefit = false,
                 std::string tag = "",
                 double ymax = -1,
                 bool do_bin_width_norm = true) {
  // Setup
  std::cout << "Plotting " << variable->Name() << std::endl;
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);

  double pot_scale = data_pot/mc_pot;
  std::string fit_str = do_prefit ? "Pre" : "Post";
  std::string label = Form("%sWFit_%s_%s_PN_%s", fit_str.c_str(), 
                           variable->m_label.c_str(), 
                           GetSignalFileTag(signal_definition).c_str(),
                           tag.c_str());
  TCanvas cE ("c1","c1"); 

  // Never don't clone when plotting
  PlotUtils::MnvH1D* h_data = (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_data->Clone("data");
  PlotUtils::MnvH1D* h_sig  = (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_sig.univHist(&universe)->Clone("sig");
  PlotUtils::MnvH1D* h_loW  = (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_loW.univHist(&universe)->Clone("loW");
  PlotUtils::MnvH1D* h_midW = (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_midW.univHist(&universe)->Clone("midW");
  PlotUtils::MnvH1D* h_hiW  = (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_hiW.univHist(&universe)->Clone("hiW");

  // Apply fit
  if(do_prefit) {
    ;
  }
  else {
    h_loW ->Scale(loW_fit.univHist(&universe) ->GetBinContent(1));
    h_midW->Scale(midW_fit.univHist(&universe)->GetBinContent(1));
    h_hiW ->Scale(hiW_fit.univHist(&universe) ->GetBinContent(1));
  }

  std::string y_label = "Events";

  // bin width norm
  if ( do_bin_width_norm ) {
    h_data->Scale(1., "width");
    h_sig->Scale(1.,  "width");
    h_loW->Scale(1.,  "width");
    h_midW->Scale(1., "width");
    h_hiW->Scale(1.,  "width");
    y_label = "Events / MeV";
  }

  if( ymax<0 ) ymax = h_data->GetBinContent(h_data->GetMaximumBin()) * 1.6;
  if( ymax>0 ) mnvPlotter.axis_maximum = ymax;

  // Prepare stack
    std::string legend_name = GetTruthClassification_LegendLabel(kWSideband_Signal);
    h_sig->SetTitle(legend_name.c_str());

    legend_name = GetTruthClassification_LegendLabel(kWSideband_Low);
    h_loW->SetTitle(legend_name.c_str());

    legend_name = GetTruthClassification_LegendLabel(kWSideband_Mid);
    h_midW->SetTitle(legend_name.c_str());

    legend_name = GetTruthClassification_LegendLabel(kWSideband_High);
    h_hiW->SetTitle(legend_name.c_str());

    SetHistColorScheme(h_sig,  int(kWSideband_Signal), sidebands::kWSideband_ColorScheme);
    SetHistColorScheme(h_loW,  int(kWSideband_Low),    sidebands::kWSideband_ColorScheme);
    SetHistColorScheme(h_midW, int(kWSideband_Mid),    sidebands::kWSideband_ColorScheme);
    SetHistColorScheme(h_hiW,  int(kWSideband_High),   sidebands::kWSideband_ColorScheme);

    TObjArray* array = new TObjArray();
    array->Add(h_sig);
    array->Add(h_loW);
    array->Add(h_midW);
    array->Add(h_hiW);

  // Draw
  mnvPlotter.DrawDataStackedMC(h_data, array, pot_scale, "TR",  "Data", -1, -1, 
                               1001, variable->m_hists.m_xlabel.c_str(), y_label.c_str());

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
                    const TObjArray& array_mc, 
                    float data_pot, float mc_pot, 
                    SignalDefinition signal_definition, 
                    std::string tag = "",
                    double ymax = -1,
                    bool draw_arrow = false, bool do_bin_width_norm = true) {

  // Never don't clone when plotting
  PlotUtils::MnvH1D* data = (PlotUtils::MnvH1D*)h_data->Clone("data");
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  if (ymax >0) mnvPlotter.axis_maximum = ymax;
  if(tag == "FSP") mnvPlotter.legend_offset_x = .15;
  if(tag == "Hadrons") mnvPlotter.legend_offset_x = .06;
  if(tag == "Wexp") mnvPlotter.legend_offset_x = .06;
  mnvPlotter.legend_text_size = 0.0405;

  double pot_scale = data_pot/mc_pot;
  std::string label = Form("Breakdown_BG_%s_%s_%s",
                           GetSignalFileTag(signal_definition).c_str(), 
                           variable->m_label.c_str(), tag.c_str());

  std::string y_label = "Events";

  // bin width norm
  if ( do_bin_width_norm ) {
    data->Scale(1., "width");
    for (auto h : array)
      dynamic_cast<PlotUtils::MnvH1D*>(h)->Scale(1., "width");
    y_label = "Events / MeV";
  }

  TCanvas cE ("c1","c1"); 
  mnvPlotter.DrawDataStackedMC(data, &array, pot_scale, "TR",  "Data", -1, -1, 
                               1001, variable->m_hists.m_xlabel.c_str(), y_label.c_str());

  if (draw_arrow) {
    //double arrow_height = data->GetBinContent(data->GetMaximumBin()) *
    //                      data->GetNormBinWidth()/data->GetBinWidth(data->GetMaximumBin());
    double arrow_height = 250;
    double arrow_location = (signal_definition == kOnePi) ? 1400 : 1800;
    mnvPlotter.AddCutArrow(arrow_location, 0.0, arrow_height, 200., "L");
  }

  mnvPlotter.WritePreliminary("TL");
  mnvPlotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);
  mnvPlotter.AddHistoTitle("Background");
  mnvPlotter.MultiPrint(&cE, label, "png");
  //mnvPlotter.MultiPrint(&cE, label, "eps");
}

void PlotTH1_1(TH1* h1, std::string tag, double ymax = -1, bool do_log_scale = false, bool do_fit = false) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  //TH1::SetDefaultSumw2();

  TCanvas cF ("c4","c4"); 

  h1->SetTitle(tag.c_str());
  h1->Draw("HIST");

  if(do_log_scale)
    cF.SetLogy();

  cF.Update();

  if (ymax > 0)
    h1->GetYaxis()->SetRangeUser(0,ymax);

  if(do_fit){
    // gaussian fit
    h1->Fit("gaus","0");
    h1->Draw("");
    cF.Update();
    TF1 *fit1 = (TF1*)h1->GetFunction("gaus");
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
  //cF.Print(Form("%s.eps", tag.c_str()));
}

#endif

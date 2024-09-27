
#ifndef PrintBreakDown_C
#define PrintBreakDown_C

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

#include "TLatex.h"
#include "TPaveText.h"
#include "includes/MacroUtil.h"
#include "includes/Plotter.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "plotting_functions.h"
#include "xsec/makeCrossSectionMCInputs.C"  // GetAnalysisVariables
/*
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
*/

//==============================================================================
// Main
//==============================================================================
void PrintBreakDown() {
  // Infiles

  //  std::string macro = "Sideband";
  std::string macro = "Background";
//  std::string macro = "DataSelection";
  //std::string macro = "AfterWcut";
//  std::string macro = "BeforeWcut";
  //TFile fin("DataSelection_Breakdown_20240118_ALL_mixed_p4.root", "READ");
  //TFile fin(Form("%s_Breakdown_ALL.root", macro.c_str()), "READ");
  TFile fin("BackgroundUntrackedOnePion_20240729_me1A.root", "READ");
  //  TFile fin1("DataXSecInputs_20231103_ME1A_mixed_Sys_P4.root", "READ");
  TFile fin1("DataXSecInputs_20240622_ALL_mixed_newtpibinning_noSys_p4.root", "READ"); 
  cout << "Reading input from " << fin.GetName() << endl;
  std::vector<string> variables = { "q2","ptmu", "wexp", "enu", "thetamu_deg", "mixthetapi_deg", 
                                   "pmu", "mixtpi", "pzmu", "wexp_fit"};
  std::vector<string> BdTypes = {"FSP",  "Int", "Hadrons", "Npi", "Npi0",
                                 "Npip", "WSB", "Msn",     "WBG"};
  std::string units;
  std::string label;

  //  PlotUtils::MnvH1D* h_data_pot = (PlotUtils::MnvH1D*)fin.Get("data_pot");
  //  PlotUtils::MnvH1D* h_mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  double data_pot = 8.97852e+19;
  //  double mc_pot = h_mc_pot->GetBinContent(1);
  double pot_scale = 1.;  // data_pot/mc_pot;
  std::vector<int> idxs = {6, 4, 6, 4, 3, 4, 4, 3, 4};
  if (macro == "Background"){
  for (int v = 0; v < (int)variables.size(); v++) {
    if (variables[v] == "q2") units = "MeV^{2}";
    if (variables[v] == "enu") units = "MeV^{2}";
    if (variables[v] == "ptmu") units = "MeV";
    if (variables[v] == "wexp") units = "MeV";
    if (variables[v] == "mixtpi") units = "MeV";
    if (variables[v] == "pzmu") units = "MeV";
    if (variables[v] == "wexp_fit") units = "MeV";
    if (variables[v] == "pmu") units = "MeV";
    if (variables[v] == "thetamu_deg") units = "deg";
    if (variables[v] == "mixthetapi_deg") units = "deg";

    if (variables[v] == "q2") label = "Q^{2}";
    if (variables[v] == "enu") label = "E_{#nu}";
    if (variables[v] == "ptmu") label = "p^{T}_{#mu}";
    if (variables[v] == "wexp") label = "W_{exp}";
    if (variables[v] == "mixtpi") label = "T_{#pi}";
    if (variables[v] == "pzmu") label = "p^{z}_{#mu}";
    if (variables[v] == "wexp_fit") label = "W_{exp}";
    if (variables[v] == "pmu") label = "p_{#mu}";
    if (variables[v] == "thetamu_deg") label = "#theta_{#mu}";
    if (variables[v] == "mixthetapi_deg") label = "#theta_{#pi}";

    for (int t = 0; t < (int)BdTypes.size(); t++) {
      TCanvas cE("c1", "c1");
      auto hs = new THStack("hs", Form("Breakdown %s", macro.c_str()));
      //  TObjArray* stack = new TObjArray();
      PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
      auto legend = new TLegend(0.6, 0.7, 0.84, 0.9);
      legend->SetFillStyle(0);
      legend->SetLineColor(0);
      legend->SetBorderSize(0);
      legend->SetTextFont(62);

      for (int i = 0; i < idxs[t]; i++) {
        TObject* obj = fin.Get(
            Form("%s_%s_%d", variables[v].c_str(), BdTypes[t].c_str(), i));
        TH1* h = dynamic_cast<TH1*>(obj);
        h->Scale(1., "width");
        //    stack->Add(obj);
        legend->AddEntry(h, h->GetTitle(), "f");
        hs->Add(h);
      }

      auto Title = new TPaveText(0.6, 0.7, 0.84, 0.9);
      Title->Clear();
      Title->AddText("Background");
      Title->SetShadowColor(0);
      Title->SetLineColor(0);
      Title->SetFillColor(0);
      PlotUtils::MnvH1D* data = NULL;
      //      hs->SetMinimum(0.0001);
      //    hs->SetTitle(Form("%s", macro.c_str()));
      double maxyaxis = 0.;// data->GetBinContent(data->GetMaximumBin()) * 1.6;
      if (variables[v] == "enu")
	      maxyaxis = 5.;
      if (variables[v] == "mixthetapi_deg")
	      maxyaxis = 300.;
      if (variables[v] == "mixtpi")
	      maxyaxis = 200.;
      if (variables[v] == "pmu")
	      maxyaxis = 5.;
      if (variables[v] == "ptmu")
	      maxyaxis = 30;
      if (variables[v] == "pzmu")
	      maxyaxis = 5.;
      if (variables[v] == "q2")
	      maxyaxis = 18.e-3;
      if (variables[v] == "thetamu_deg")
	      maxyaxis = 2000;
      if (variables[v] == "wexp")
	      maxyaxis = 70.;
      if (variables[v] == "wexp_fit")
	      maxyaxis = 70.;

      hs->SetMaximum(maxyaxis);
      hs->Draw("HIST");
      legend->Draw();
      //    cE.SetTitle(Form("%s", macro.c_str()));
      hs->GetXaxis()->SetTitle(Form("%s %s", label.c_str(), units.c_str()));
      hs->GetYaxis()->SetTitle(Form("Events/%s", units.c_str()));
      hs->GetXaxis()->CenterTitle();
      //      cE.SetLogy();
      Title->Draw();
      cE.Update();
      //      cE.BuildLegend(0.6,0.7,0.84,0.9);
      cE.Print(Form("BackgroundUntrackedOnePion_%s_%s_%s.png", macro.c_str(),
                    variables[v].c_str(), BdTypes[t].c_str()));
      //    mnvPlotter.DrawStackedMC(stack, 1., "TR", -1, -1,
      //                             1001, "",
      //                             variables[v].c_str());
    }
  }
  }



  if (macro == "AfterWcut" || macro == "BeforeWcut"){

  PlotUtils::MnvH1D* h_data_POT = (PlotUtils::MnvH1D*)fin1.Get("data_pot"); 
  PlotUtils::MnvH1D* h_mc_POT = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float data_pot = h_data_POT->GetBinContent(1);
  float mc_pot = h_mc_POT->GetBinContent(1);

  for (int v = 0; v < (int)variables.size(); v++) {
    if (variables[v] == "q2") units = "MeV^{2}";
    if (variables[v] == "ptmu") units = "MeV";
    if (variables[v] == "wexp") units = "MeV";
    if (variables[v] == "mixtpi") units = "MeV";
    if (variables[v] == "pzmu") units = "MeV";
    if (variables[v] == "wexp_fit") units = "MeV";
    if (variables[v] == "pmu") units = "MeV";
    if (variables[v] == "thetamu_deg") units = "deg";
    if (variables[v] == "mixthetapi_deg") units = "deg";
    if (variables[v] == "enu") units = "MeV";

    if (variables[v] == "q2") label = "Q^{2}";
    if (variables[v] == "ptmu") label = "p^{T}_{#mu}";
    if (variables[v] == "wexp") label = "W_{exp}";
    if (variables[v] == "mixtpi") label = "T_{#pi}";
    if (variables[v] == "pzmu") label = "p^{z}_{#mu}";
    if (variables[v] == "wexp_fit") label = "W_{exp}";
    if (variables[v] == "pmu") label = "p_{#mu}";
    if (variables[v] == "thetamu_deg") label = "#theta_{#mu}";
    if (variables[v] == "mixthetapi_deg") label = "#theta_{#pi}";
    if (variables[v] == "enu") label = "E_{#nu}";

    for (int t = 0; t < (int)BdTypes.size(); t++) {
      TCanvas cE("c1", "c1");
      cE.Update();
      auto hs = new THStack("hs", Form("Breakdown %s", macro.c_str()));
      //  TObjArray* stack = new TObjArray();
      PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
      auto legend = new TLegend(0.6, 0.7, 0.84, 0.9);
      legend->SetFillStyle(0);
      legend->SetLineColor(0);
      legend->SetBorderSize(0);
      legend->SetTextFont(62);
      
      PlotUtils::MnvH1D* dummy_data = nullptr;
      if (macro == "AfterWcut"){
	 dummy_data = (PlotUtils::MnvH1D*)fin1.Get(
                   Form("selection_data_%s", variables[v].c_str()));
      }
      else if (macro == "BeforeWcut"){
         dummy_data = (PlotUtils::MnvH1D*)fin1.Get(
		      Form("wsideband_data_%s", variables[v].c_str()));
      }
      else {
        std::cerr << "Error: Invalid dummy_data pointer." << std::endl;
        continue;
      }

      if (!dummy_data){
        std::cerr << "Error: Invalid dummy_data pointer." << std::endl;
        continue;
      }
      TH1D* data = (TH1D*)dummy_data->Clone(Form("data_%s", variables[v].c_str()));
      data->Scale(1.,"width");
      data->SetMarkerStyle(20);
      data->SetMarkerSize(1);
      data->SetLineWidth(1);
      data->SetLineStyle(1);
      data->SetLineColor(1);

//      legend->AddEntry(data,"Data","ple");


      for (int i = 0; i < idxs[t]; i++) {
        std::cout << variables[v] << "_" << BdTypes[t] << "_" << i << "\n";
        TObject* obj = fin.Get(
            Form("%s_%s_%d", variables[v].c_str(), BdTypes[t].c_str(), i));
        TH1* h = dynamic_cast<TH1*>(obj);
        h->Scale(data_pot/mc_pot, "width");
        //    stack->Add(obj);
        legend->AddEntry(h, h->GetTitle(), "f");
        hs->Add(h);
      }
      auto Title = new TPaveText(0.6, 0.7, 0.84, 0.9);
      Title->Clear();
      Title->AddText("Data Selection");
      Title->SetShadowColor(0);
      Title->SetLineColor(0);
      Title->SetFillColor(0);

      TLatex *latex = new TLatex( 0.5, 0.96, Form("%s", macro.c_str()));
//      AddToTmp( latex ); 

      latex->SetNDC();
      latex->SetTextSize(0.05);
      latex->SetTextColor(1);
      latex->SetTextFont(62);
      latex->SetTextAlign(22);
      latex->SetTextAngle(0);
      latex->Draw();




      //      hs->SetMinimum(0.0001);
      hs->SetTitle(Form("%s", macro.c_str()));
      double maxyaxis = 0.;// data->GetBinContent(data->GetMaximumBin()) * 1.6;
      if (variables[v] == "enu")
	      maxyaxis = 30.;
      if (variables[v] == "mixthetapi_deg")
	      maxyaxis = 2000.;
      if (variables[v] == "mixtpi")
	      maxyaxis = 1200.;
      if (variables[v] == "pmu")
	      maxyaxis = 30.;
      if (variables[v] == "ptmu")
	      maxyaxis = 200;
      if (variables[v] == "pzmu")
	      maxyaxis = 30.;
      if (variables[v] == "q2")
	      maxyaxis = 0.3;
      if (variables[v] == "thetamu_deg")
	      maxyaxis = 16000;
      if (variables[v] == "wexp")
	      maxyaxis = 350.;
      if (variables[v] == "wexp_fit")
	      maxyaxis = 350.;

      hs->SetMaximum(maxyaxis);
      hs->Draw("HIST");
//      data->Draw("SAME E1 X0");
      legend->Draw();
      //cE.SetTitle(Form("%s", macro.c_str()));
      hs->GetXaxis()->SetTitle(Form("%s %s", label.c_str(), units.c_str()));
      hs->GetYaxis()->SetTitle(Form("Events/%s", units.c_str()));
      hs->GetXaxis()->CenterTitle();
      //if (variables[v] == "q2")      cE.SetLogx();
      //Title->Draw();
      latex->Draw();
      cE.Update();
      //      cE.BuildLegend(0.6,0.7,0.84,0.9);
      cE.Print(Form("Bd_%s_%s_%s.png", macro.c_str(),
                    variables[v].c_str(), BdTypes[t].c_str()));
      //    mnvPlotter.DrawStackedMC(stack, 1., "TR", -1, -1,
      //                             1001, "",
      //                             variables[v].c_str());
    }
  }
  }


  if (macro == "DataSelection"){

  PlotUtils::MnvH1D* h_data_POT = (PlotUtils::MnvH1D*)fin1.Get("data_pot"); 
  PlotUtils::MnvH1D* h_mc_POT = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float data_pot = h_data_POT->GetBinContent(1);
  float mc_pot = h_mc_POT->GetBinContent(1);

  for (int v = 0; v < (int)variables.size(); v++) {
    if (variables[v] == "q2") units = "MeV^{2}";
    if (variables[v] == "ptmu") units = "MeV";
    if (variables[v] == "wexp") units = "MeV";
    if (variables[v] == "mixtpi") units = "MeV";
    if (variables[v] == "pzmu") units = "MeV";
    if (variables[v] == "wexp_fit") units = "MeV";
    if (variables[v] == "pmu") units = "MeV";

    if (variables[v] == "q2") label = "Q^{2}";
    if (variables[v] == "ptmu") label = "p^{T}_{#mu}";
    if (variables[v] == "wexp") label = "W_{exp}";
    if (variables[v] == "mixtpi") label = "T_{#pi}";
    if (variables[v] == "pzmu") label = "p^{z}_{#mu}";
    if (variables[v] == "wexp_fit") label = "W_{exp}";
    if (variables[v] == "pmu") label = "p_{#mu}";

    for (int t = 0; t < (int)BdTypes.size(); t++) {
      TCanvas cE("c1", "c1");
      auto hs = new THStack("hs", Form("Breakdown %s", macro.c_str()));
      //  TObjArray* stack = new TObjArray();
      PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
      auto legend = new TLegend(0.6, 0.7, 0.84, 0.9);
      legend->SetFillStyle(0);
      legend->SetLineColor(0);
      legend->SetBorderSize(0);
      legend->SetTextFont(62);

      PlotUtils::MnvH1D* dummy_data = (PlotUtils::MnvH1D*)fin1.Get(
                   Form("selection_data_%s", variables[v].c_str()));
      if (!dummy_data){
        std::cerr << "Error: Invalid dummy_data pointer." << std::endl;
        continue;
      }
      TH1D* data = (TH1D*)dummy_data->Clone(Form("data_%s", variables[v].c_str()));
      data->Scale(1.,"width");
      data->SetMarkerStyle(20);
      data->SetMarkerSize(1);
      data->SetLineWidth(1);
      data->SetLineStyle(1);
      data->SetLineColor(1);

      legend->AddEntry(data,"Data","ple");


      for (int i = 0; i < idxs[t]; i++) {
        std::cout << variables[v] << "_" << BdTypes[t] << "_" << i << "\n";
        TObject* obj = fin.Get(
            Form("%s_%s_%d", variables[v].c_str(), BdTypes[t].c_str(), i));
        TH1* h = dynamic_cast<TH1*>(obj);
        h->Scale(data_pot/mc_pot, "width");
        //    stack->Add(obj);
        legend->AddEntry(h, h->GetTitle(), "f");
        hs->Add(h);
      }
      auto Title = new TPaveText(0.6, 0.7, 0.84, 0.9);
      Title->Clear();
      Title->AddText("Data Selection");
      Title->SetShadowColor(0);
      Title->SetLineColor(0);
      Title->SetFillColor(0);
      //      hs->SetMinimum(0.0001);
      //    hs->SetTitle(Form("%s", macro.c_str()));
      double maxyaxis = data->GetBinContent(data->GetMaximumBin()) * 1.6;

      hs->SetMaximum(maxyaxis);
      hs->Draw("HIST");
      data->Draw("SAME E1 X0");
      legend->Draw();
      //    cE.SetTitle(Form("%s", macro.c_str()));
      hs->GetXaxis()->SetTitle(Form("%s %s", label.c_str(), units.c_str()));
      hs->GetYaxis()->SetTitle(Form("Events/%s", units.c_str()));
      hs->GetXaxis()->CenterTitle();
      if (variables[v] == "q2")      cE.SetLogx();
      Title->Draw();
      cE.Update();
      //      cE.BuildLegend(0.6,0.7,0.84,0.9);
      cE.Print(Form("Bd_%s_%s_%s.png", macro.c_str(),
                    variables[v].c_str(), BdTypes[t].c_str()));
      //    mnvPlotter.DrawStackedMC(stack, 1., "TR", -1, -1,
      //                             1001, "",
      //                             variables[v].c_str());
    }
  }
  }
  // PLOT W Sideband Fit
  if (macro == "Sideband") {
    PlotUtils::MnvH1D* loW_fit_wgt = (PlotUtils::MnvH1D*)fin.Get("loW_fit_wgt");
    PlotUtils::MnvH1D* midW_fit_wgt =
        (PlotUtils::MnvH1D*)fin.Get("midW_fit_wgt");
    PlotUtils::MnvH1D* hiW_fit_wgt = (PlotUtils::MnvH1D*)fin.Get("hiW_fit_wgt");

    std::vector<string> preposfit = {"prefit", "posfit"};
    double ymax = -1;
    for (int p = 0; p < (int)preposfit.size(); p++) {
      std::string fit = preposfit[p];
      for (int v = 0; v < (int)variables.size(); v++) {
        std::string name = variables[v];
        if (name == "q2") units = "MeV^{2}";
        if (name == "ptmu") units = "MeV";
        if (name == "wexp") units = "MeV";
        if (name == "mixtpi") units = "MeV";
        if (name == "pzmu") units = "MeV";
        if (name == "wexp_fit") units = "MeV";
        if (name == "pmu") units = "MeV";

        if (name == "q2") label = "Q^{2}";
        if (name == "ptmu") label = "p^{T}_{#mu}";
        if (name == "wexp") label = "W_{exp}";
        if (name == "mixtpi") label = "T_{#pi}";
        if (name == "pzmu") label = "p^{z}_{#mu}";
        if (name == "wexp_fit") label = "W_{exp}";
        if (name == "pmu") label = "p_{#mu}";
        PlotUtils::MnvH1D* h_sig =
            (PlotUtils::MnvH1D*)fin.Get(Form("%s_sig_WSideband", name.c_str()));
        PlotUtils::MnvH1D* h_loW =
            (PlotUtils::MnvH1D*)fin.Get(Form("%s_loW_WSideband", name.c_str()));
        PlotUtils::MnvH1D* h_midW = (PlotUtils::MnvH1D*)fin.Get(
            Form("%s_midW_WSideband", name.c_str()));
        PlotUtils::MnvH1D* h_hiW =
            (PlotUtils::MnvH1D*)fin.Get(Form("%s_hiW_WSideband", name.c_str()));
        PlotUtils::MnvH1D* h_data = (PlotUtils::MnvH1D*)fin.Get(
            Form("%s_data_WSideband", name.c_str()));
        PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
        TCanvas cE("c1", "c1");
        TObjArray* stack = new TObjArray();
        h_sig->SetTitle("Signal");
        h_loW->SetTitle("Low W");
        h_midW->SetTitle("Mid W");
        h_hiW->SetTitle("High W");
        h_sig->GetYaxis()->SetTitle(Form("Events/%s", units.c_str()));

        h_sig->Scale(1., "width");
        h_loW->Scale(1., "width");
        h_midW->Scale(1., "width");
        h_hiW->Scale(1., "width");
        h_data->Scale(1., "width");

        if (fit == "posfit") {
          h_loW->Scale(loW_fit_wgt->GetBinContent(1));
          h_midW->Scale(midW_fit_wgt->GetBinContent(1));
          h_hiW->Scale(hiW_fit_wgt->GetBinContent(1));
        }
        std::string xlabel = label + " " + units;
        stack->Add(h_sig);
        stack->Add(h_loW);
        stack->Add(h_midW);
        stack->Add(h_hiW);
        mnvPlotter.DrawDataStackedMC(h_data, stack, pot_scale, "TR", "Data", 2,
                                     1, 3001, xlabel.c_str(),
                                     Form("Events/%s", units.c_str()));
        //    mnvPlotter.DrawStackedMC(stack, 1.0, "TR", 2, 1, 3001,
        //              Form("%s %s", name.c_str(), units.c_str()));
        mnvPlotter.AddHistoTitle(Form("SidebandRegion %s", fit.c_str()), 0.05);

        std::string plotname = "Sideband_" + fit + "_" + name;
        mnvPlotter.MultiPrint(&cE, plotname, "png");
      }
    }
    // TODO plot pre/postfit
    /*    for (auto var : variables) {
          tag = "SidebandRegion";
          bool do_prefit = true;
          bool do_bin_width_norm = true;
          CVUniverse* universe = util.m_error_bands.at("cv").at(0);
          PlotFittedW(var, *universe, var->m_hists.m_bg_loW,
       var->m_hists.m_bg_midW, var->m_hists.m_bg_hiW, util.m_data_pot,
       util.m_mc_pot, util.m_signal_definition, do_prefit, tag, ymax,
                      do_bin_width_norm);
          do_prefit = false;
          PlotFittedW(var, *universe, var->m_hists.m_bg_loW,
       var->m_hists.m_bg_midW, var->m_hists.m_bg_hiW, util.m_data_pot,
       util.m_mc_pot, util.m_signal_definition, do_prefit, tag, ymax,
                      do_bin_width_norm);
        }*/
  }
}

#endif  // plotCrossSectionFromFile_C

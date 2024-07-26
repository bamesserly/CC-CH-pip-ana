
#ifndef RatioPlots_C
#define RatioPlots_C

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
#include "includes/Plotter.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "plotting_functions.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/FluxReweighter.h"
//==============================================================================
//Main
//==============================================================================
void Validations() {
  // Infiles

  TFile fp4("DataXSecInputs_20240617_me1A_mixed_breakdown_noSys_p4.root", "READ");
  TFile fp6("DataXSecInputs_20240617_me1A_mixed_breakdown_noSys_p6.root", "READ");
//  TFile fp6("MCXSecInputs_20240417_me1P_AaronSigDef_PEMA_nosys_p4.root", "READ");
//TFile fp6("../../aux/DataXSecInputs_20231127_ME1A_mixed_noSys_p4.root", "READ");

  std::vector<std::string> variables = {"q2", "mixtpi", "thetamu_deg", "enu", "wexp", "ptmu", "pmu", "pzmu"};	 

  PlotUtils::MnvH1D* mc_pot_p4 = (PlotUtils::MnvH1D*)fp4.Get("mc_pot");
  PlotUtils::MnvH1D* data_pot_p4 = (PlotUtils::MnvH1D*)fp4.Get("data_pot");
  PlotUtils::MnvH1D* mc_pot_p6 = (PlotUtils::MnvH1D*)fp6.Get("mc_pot");
  PlotUtils::MnvH1D* data_pot_p6 = (PlotUtils::MnvH1D*)fp6.Get("data_pot");
  double mc_scale = mc_pot_p6->GetBinContent(1)/mc_pot_p4->GetBinContent(1);
  double data_scale = data_pot_p4->GetBinContent(1)/data_pot_p6->GetBinContent(1);

  for (int i = 0; i < (int)variables.size(); i++){
    std::string var = variables[i];
    std::string Aaronvar;
    if (var == "q2")
      Aaronvar = "q2";
    else if (var == "thetamu_deg")
      Aaronvar = "muon_theta";
    else if (var == "thetapi_deg")
      Aaronvar = "pion_theta";
    else if (var == "mixtpi")
      Aaronvar = "pion_ekin";
    else if (var == "wexp"){
      Aaronvar = "W";
    }
    else if (var == "pmu"){
      Aaronvar = "muon_p";
    }
    else if (var == "pzmu"){
      Aaronvar = "muon_pz";
    }
    else if (var == "ptmu")
      Aaronvar = "muon_pt";

    std::string xlabel = "";
    if (var == "q2"){
      xlabel = "Q^{2} (MeV^{2})";
    }
    if (var == "mixtpi"){
      xlabel = "T_{#pi} (MeV)";
    }
    if (var == "thetamu_deg"){
      xlabel = "#theta_{#mu} (deg)";
    }
    if (var == "enu"){
      xlabel = "E_{#nu} (MeV)";
    }
    if (var == "wexp"){
      xlabel = "W_{exp} (MeV)";
    }
    if (var == "ptmu"){
      xlabel = "p^{T}_{#mu} (MeV)";
    }
    if (var == "pzmu"){
      xlabel = "p^{z}_{#mu} (MeV)";
    }
    if (var == "pmu"){
      xlabel = "p_{#mu} (MeV)";
    }

    PlotUtils::MnvH1D* SelMCP4 = (PlotUtils::MnvH1D*)fp4.Get(Form("selection_mc_%s", var.c_str()));
    PlotUtils::MnvH1D* SelMCP6 = (PlotUtils::MnvH1D*)fp6.Get(Form("selection_mc_%s", var.c_str()));
    PlotUtils::MnvH1D* XsecDataP4 = (PlotUtils::MnvH1D*)fp4.Get(Form("cross_section_%s", var.c_str()));
    PlotUtils::MnvH1D* XsecDataP6 = (PlotUtils::MnvH1D*)fp6.Get(Form("cross_section_%s", var.c_str()));
    PlotUtils::MnvH1D* effnumP4 = (PlotUtils::MnvH1D*)fp4.Get(Form("effnum_%s", var.c_str()));
    PlotUtils::MnvH1D* effnumP6 = (PlotUtils::MnvH1D*)fp6.Get(Form("effnum_%s", var.c_str()));
    PlotUtils::MnvH1D* XsecMCP4 = (PlotUtils::MnvH1D*)fp4.Get(Form("mc_cross_section_%s", var.c_str()));
    PlotUtils::MnvH1D* XsecMCP6 = (PlotUtils::MnvH1D*)fp6.Get(Form("mc_cross_section_%s", var.c_str()));
    PlotUtils::MnvH1D* SelDataP4 = (PlotUtils::MnvH1D*)fp4.Get(Form("selection_data_%s", var.c_str()));
    PlotUtils::MnvH1D* SelDataP6 = (PlotUtils::MnvH1D*)fp6.Get(Form("selection_data_%s", var.c_str()));
    PlotUtils::MnvH1D* effdenP4 = (PlotUtils::MnvH1D*)fp4.Get(Form("effden_%s_true", var.c_str()));
    PlotUtils::MnvH1D* effdenP6 = (PlotUtils::MnvH1D*)fp6.Get(Form("effden_%s_true", var.c_str()));
    PlotUtils::MnvH1D* bgsubdatap4 = (PlotUtils::MnvH1D*)fp4.Get(Form("bg_subbed_data_%s", var.c_str()));
    PlotUtils::MnvH1D* bgsubdatap6 = (PlotUtils::MnvH1D*)fp6.Get(Form("bg_subbed_data_%s", var.c_str()));
    PlotUtils::MnvH1D* bgp4 = (PlotUtils::MnvH1D*)fp4.Get(Form("bg_%s", var.c_str()));
    PlotUtils::MnvH1D* bgp6 = (PlotUtils::MnvH1D*)fp6.Get(Form("bg_%s", var.c_str()));


    PlotRatio(SelMCP6, SelMCP4, var, 1., "p4p6SelMC", false,       
             "P6/P4", xlabel); 
    PlotRatio(XsecDataP6, XsecDataP4, var, 1., "p4p6XsecData",  true, 
             "P6/P4", xlabel); 
    PlotRatio(SelDataP6, SelDataP4, var, data_scale, "p4p6SelData", true, 
             "P6/P4", xlabel); 
    PlotRatio(XsecMCP6, XsecMCP4, var, 1., "p4p6XsecMC",  true,       
             "P6/P4", xlabel); 
    PlotRatio(SelMCP6, SelMCP4, var, mc_scale, "p4p6SelMC", true,       
             "P6/P4", xlabel); 
    PlotRatio(effnumP6, effnumP4, var, mc_scale, "p4p6effnum", true,       
             "P6/P4", xlabel); 
    PlotRatio(effdenP6, effdenP4, var, mc_scale, "p4p6effden", true,       
             "P6/P4", xlabel); 
    PlotRatio(bgsubdatap6, bgsubdatap4, var, data_scale, "p4p6bgsubdata",  true,       
             "P6/P4", xlabel); 
    PlotRatio(bgp6, bgp4, var, mc_scale, "p4p6bg",  true,       
             "P6/P4", xlabel); 
  }
}
#endif  

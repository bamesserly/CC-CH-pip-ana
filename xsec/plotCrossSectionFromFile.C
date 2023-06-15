#ifndef plotCrossSectionFromFile_C
#define plotCrossSectionFromFile_C

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

//==============================================================================
// Main
//==============================================================================
void plotCrossSectionFromFile(int signal_definition_int = 0,
                              int plot_errors = 0) {
  // Infiles
  TFile fin("DataXSecInputs_20230611_AaronSigDef_AllAaronCuts.root", "READ");
  cout << "Reading input from " << fin.GetName() << endl;

  TFile finCCPi("/minerva/data/users/abercell/hists/Macro/GridOneLoop_MENU1PI_MinosMatched_plastic_Merged_NewdEdXCal_MinervaME1ABCDEFGLMNOP_Data_Merged_NewdEdXCal_Tracker_MinervaME1ABCDEFGLMNOP_MC.root", "READ");
  TFile fAaronSigBenMacro("DataXSecInputs_20230611_AaronSigDef_AllAaronCuts.root", "READ");

  //    TFile
  // Set up macro utility object...which gets the list of systematics for us...
  // which we need in order to read in HistWrappers...which we don't need at
  // this point...indeed we only need MnvH1Ds...so that's a TODO: write a
  // function that only loads in MnvH1D's from a file, not HWs.
  // Most of these options aren't used in this script. TODO make a CTOR that
  // doesn't require them.
  // INPUT TUPLES
  // Don't actually use the MC chain, only load it to indirectly access it's
  // systematics
  const std::string plist = "ME1A";
  std::string data_file_list = GetPlaylistFile(plist, false, false);
  std::string mc_file_list = GetPlaylistFile(plist, true, false);
  //std::string data_file_list = GetTestPlaylist(false);
  //std::string mc_file_list = GetTestPlaylist(true);

  // Macro Utility
  const std::string macro("PlotCrossSectionFromFile");
  bool do_truth = false, is_grid = false, do_systematics = true;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  CCPi::MacroUtil utilASig(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  // Get POT from file, not from any chain
  SetPOT(fin, util);
  SetPOT(fAaronSigBenMacro, utilASig);
  util.PrintMacroConfiguration(macro);

  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  for (auto var : variables) {
    if (var->Name() == sidebands::kFitVarString) {
      var->m_hists.m_stacked_wsideband = StackedHistogram<WSidebandType>(
          var->m_hists.m_label, var->m_hists.m_xlabel,
          var->m_hists.m_bins_array, kNWSidebandTypes,
          sidebands::kWSideband_ColorScheme);
    }
    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    var->LoadDataHistsFromFile(fin);

    // if(var->Name() == "ptmu")
    //  PrintUniverseContent(var->m_hists.m_cross_section);
  }

  {  // remove unwanted variables
    ContainerEraser::erase_if(variables, [](Variable* v) {
      return v->Name() == "tpi_mbr";  // || v->Name() == "wexp_fit";
    });

    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "tpi" || v->Name() == "enu"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "thetapi_deg" || v->Name() == "thetamu_deg"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "q2" || v->Name() == "wexp"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "ptmu" || v->Name() == "pzmu"; });
  }

  if (true){
    int q2nbins = CCPi::GetBinning("q2").GetSize()-1;
    PlotUtils::MnvH1D *MC_POT_Aaron = (PlotUtils::MnvH1D*)finCCPi.Get("h_MC_POT");
    PlotUtils::MnvH1D *data_POT_Aaron = (PlotUtils::MnvH1D*)finCCPi.Get("h_Data_POT");
    double MC_POT_A = MC_POT_Aaron->GetBinContent(1);
    double data_POT_A = data_POT_Aaron->GetBinContent(1);
    TFile fAaronxSec("/minerva/data/users/abercell/hists/xsec/xsec_new_jeffrey_flux_MENU1PI_plastic_MinervaME1ABCDEFGLMNOP.root", "READ");
    TFile fAaronBGs("/minerva/data/users/abercell/hists/xsec_inputs/Merge_BkgdSub_Unfold_MENU1PI_POTNorm_plastic_MinervaME1ABCDEFGLMNOP.root", "READ");
    TFile fAaronSideBands("/minerva/data/users/abercell/hists/Sideband/FittedDists/W_Sideband_MENU1PI_plastic_FittedDists.root", "READ");
    PlotUtils::MnvH1D *q2_xsec_Aaron_paper = new PlotUtils::MnvH1D("q2_xsec_Aaron_paper","q2_xsec_Aaron_paper", q2nbins, CCPi::GetBinning("q2").GetArray());

    // Cross Section reported on Aaron's paper
    PlotUtils::MnvH1D* ASigBenXSecMCq2 = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("mc_cross_section_q2");
    PlotUtils::MnvH1D* ASigBenXSecdataq2 = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("cross_section_q2");
    std::vector<double> xsec;
    xsec.push_back(90.26e-42); 
    xsec.push_back(104.81e-42);
    xsec.push_back(219.20e-42);
    xsec.push_back(402.76e-42);
    xsec.push_back(344.86e-42);
    xsec.push_back(272.35e-42);
    xsec.push_back(208.85e-42);
    xsec.push_back(279.50e-42);
    xsec.push_back(212.70e-42);
    xsec.push_back(98.25e-42);
    xsec.push_back(62.76e-42);
    xsec.push_back(20.81e-42);

    for (int i = 1; i <= q2nbins; ++i){
      q2_xsec_Aaron_paper->SetBinContent(i, xsec[i-1]);
    }

    PlotRatio(ASigBenXSecMCq2, q2_xsec_Aaron_paper, "q2", 1., "mc", false, true, "BenMacro/Aaron'sPaper", "Q^{2} MeV");
    PlotRatio(ASigBenXSecMCq2, q2_xsec_Aaron_paper, "q2", 1., "mc", false, true, "BenMacro/Aaron'sPaper", "Q^{2} MeV");
    PlotRatio(ASigBenXSecdataq2, q2_xsec_Aaron_paper, "q2", 1., "data", false, true, "BenMacro/Aaron'sPaper", "Q^{2} MeV");
    std::string niter = "2";

    int SBnbins = CCPi::GetBinning("wexp_fit").GetSize()-1;  
    const Double_t* SBbins = CCPi::GetBinning("wexp_fit").GetArray();

    PlotUtils::MnvH1D *Aaron_SB_data_prefit_aux = (PlotUtils::MnvH1D*)fAaronSideBands.Get("h_W_Sideband_Prefit_plastic_data");
    PlotUtils::MnvH1D *Aaron_SB_true_Low_prefit_aux = (PlotUtils::MnvH1D*)fAaronSideBands.Get("h_W_Sideband_Prefit_plastic_true_W_channel_lt_14_bkgd");
    PlotUtils::MnvH1D *Aaron_SB_true_Mid_prefit_aux = (PlotUtils::MnvH1D*)fAaronSideBands.Get("h_W_Sideband_Prefit_plastic_true_W_channel_14_18");
    PlotUtils::MnvH1D *Aaron_SB_true_Hi_prefit_aux = (PlotUtils::MnvH1D*)fAaronSideBands.Get("h_W_Sideband_Prefit_plastic_true_W_channel_gt_18");
    PlotUtils::MnvH1D *Aaron_SB_true_Low_posfit_aux = (PlotUtils::MnvH1D*)fAaronSideBands.Get("h_W_Sideband_AfterFit_plastic_true_W_channel_lt_14_bkgd");
    PlotUtils::MnvH1D *Aaron_SB_true_Mid_posfit_aux = (PlotUtils::MnvH1D*)fAaronSideBands.Get("h_W_Sideband_AfterFit_plastic_true_W_channel_14_18");
    PlotUtils::MnvH1D *Aaron_SB_true_Hi_posfit_aux = (PlotUtils::MnvH1D*)fAaronSideBands.Get("h_W_Sideband_AfterFit_plastic_true_W_channel_gt_18");
    PlotUtils::MnvH1D *Aaron_SB_Signal_aux = (PlotUtils::MnvH1D*)fAaronSideBands.Get("h_W_Sideband_Prefit_plastic_true_W_channel_signal");

    PlotUtils::MnvH1D *Aaron_SB_data_prefit = new PlotUtils::MnvH1D("W_Sideband_Aaron_prefit_data", "W_Sideband_Aaron_prefit_data", SBnbins, CCPi::GetBinning("wexp_fit").GetArray());
    PlotUtils::MnvH1D *Aaron_SB_true_Low_prefit = new PlotUtils::MnvH1D("W_Sideband_Aaron_prefit_true_Low", "W_Sideband_Aaron_prefit_true_Low", SBnbins, CCPi::GetBinning("wexp_fit").GetArray());
    PlotUtils::MnvH1D *Aaron_SB_true_Mid_prefit = new PlotUtils::MnvH1D("W_Sideband_Aaron_prefit_true_Mid", "W_Sideband_Aaron_prefit_true_Mid", SBnbins, CCPi::GetBinning("wexp_fit").GetArray());
    PlotUtils::MnvH1D *Aaron_SB_true_Hi_prefit = new PlotUtils::MnvH1D("W_Sideband_Aaron_prefit_true_Hi", "W_Sideband_Aaron_prefit_true_Hi", SBnbins, CCPi::GetBinning("wexp_fit").GetArray());
    PlotUtils::MnvH1D *Aaron_SB_true_Low_posfit = new PlotUtils::MnvH1D("W_Sideband_Aaron_posfit_true_Low", "W_Sideband_Aaron_posfit_true_Low", SBnbins, CCPi::GetBinning("wexp_fit").GetArray());
    PlotUtils::MnvH1D *Aaron_SB_true_Mid_posfit = new PlotUtils::MnvH1D("W_Sideband_Aaron_posfit_true_Mid", "W_Sideband_Aaron_posfit_true_Mid", SBnbins, CCPi::GetBinning("wexp_fit").GetArray());
    PlotUtils::MnvH1D *Aaron_SB_true_Hi_posfit = new PlotUtils::MnvH1D("W_Sideband_Aaron_posfit_true_Hi", "W_Sideband_Aaron_posfit_true_Hi", SBnbins, CCPi::GetBinning("wexp_fit").GetArray());
    PlotUtils::MnvH1D *Aaron_SB_Signal = new PlotUtils::MnvH1D("W_Sideband_Aaron_prefit_true_signal", "W_Sideband_Aaron_prefit_true_signal", SBnbins, CCPi::GetBinning("wexp_fit").GetArray());

    for (int i = 1; i <= SBnbins; ++i){
      double WSB_data, WSB_prefit_low, WSB_posfit_low, WSB_prefit_mid,
             WSB_posfit_mid, WSB_prefit_Hi, WSB_posfit_Hi, WSB_Signal;
      WSB_data       = Aaron_SB_data_prefit_aux->GetBinContent(i);
      WSB_prefit_low = Aaron_SB_true_Low_prefit_aux->GetBinContent(i);
      WSB_posfit_low = Aaron_SB_true_Low_posfit_aux->GetBinContent(i);
      WSB_prefit_mid = Aaron_SB_true_Mid_prefit_aux->GetBinContent(i);
      WSB_posfit_mid = Aaron_SB_true_Mid_posfit_aux->GetBinContent(i);
      WSB_prefit_Hi  = Aaron_SB_true_Hi_prefit_aux->GetBinContent(i);
      WSB_posfit_Hi  = Aaron_SB_true_Hi_posfit_aux->GetBinContent(i);
      WSB_Signal     = Aaron_SB_Signal_aux->GetBinContent(i);

      Aaron_SB_data_prefit->SetBinContent(i, WSB_data);
      Aaron_SB_true_Low_prefit->SetBinContent(i, WSB_prefit_low);
      Aaron_SB_true_Low_posfit->SetBinContent(i, WSB_posfit_low);
      Aaron_SB_true_Mid_prefit->SetBinContent(i, WSB_prefit_mid);
      Aaron_SB_true_Mid_posfit->SetBinContent(i, WSB_posfit_mid);
      Aaron_SB_true_Hi_prefit->SetBinContent(i, WSB_prefit_Hi);
      Aaron_SB_true_Hi_posfit->SetBinContent(i, WSB_posfit_Hi);
      Aaron_SB_Signal->SetBinContent(i, WSB_Signal);
    }

    PlotUtils::MnvH1D* ASigBenWSB_data = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("wsidebandfit_data_wexp_fit");
    PlotUtils::MnvH1D* ASigBenWSB_Low_prefit = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("wsidebandfit_loW_wexp_fit");
    PlotUtils::MnvH1D* ASigBenWSB_Mid_prefit = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("wsidebandfit_midW_wexp_fit");
    PlotUtils::MnvH1D* ASigBenWSB_Hi_prefit = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("wsidebandfit_hiW_wexp_fit");
    PlotUtils::MnvH1D* ASigBenWSB_Signal = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("wsidebandfit_sig_wexp_fit");
    PlotUtils::MnvH1D* ASigBen_Low_wgt = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("loW_fit_wgt");
    PlotUtils::MnvH1D* ASigBen_Mid_wgt = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("midW_fit_wgt");
    PlotUtils::MnvH1D* ASigBen_Hi_wgt = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get("hiW_fit_wgt");
    PlotUtils::MnvH1D *ASigBenWSB_Low_posfit = ASigBenWSB_Low_prefit->Clone(uniq());
    PlotUtils::MnvH1D *ASigBenWSB_Mid_posfit = ASigBenWSB_Mid_prefit->Clone(uniq());
    PlotUtils::MnvH1D *ASigBenWSB_Hi_posfit = ASigBenWSB_Hi_prefit->Clone(uniq());
    ASigBenWSB_Low_posfit->Scale(ASigBen_Low_wgt->GetBinContent(1)); 
    ASigBenWSB_Mid_posfit->Scale(ASigBen_Mid_wgt->GetBinContent(1)); 
    ASigBenWSB_Hi_posfit->Scale(ASigBen_Hi_wgt->GetBinContent(1)); 

    PlotUtils::MnvH1D* BSigBenWSB_data = (PlotUtils::MnvH1D*)fin.Get("wsidebandfit_data_wexp_fit");
    PlotUtils::MnvH1D* BSigBenWSB_Low_prefit = (PlotUtils::MnvH1D*)fin.Get("wsidebandfit_loW_wexp_fit");
    PlotUtils::MnvH1D* BSigBenWSB_Mid_prefit = (PlotUtils::MnvH1D*)fin.Get("wsidebandfit_midW_wexp_fit");
    PlotUtils::MnvH1D* BSigBenWSB_Hi_prefit = (PlotUtils::MnvH1D*)fin.Get("wsidebandfit_hiW_wexp_fit");
    PlotUtils::MnvH1D* BSigBenWSB_Signal = (PlotUtils::MnvH1D*)fin.Get("wsidebandfit_sig_wexp_fit");
    PlotUtils::MnvH1D* BSigBen_Low_wgt = (PlotUtils::MnvH1D*)fin.Get("loW_fit_wgt");
    PlotUtils::MnvH1D* BSigBen_Mid_wgt = (PlotUtils::MnvH1D*)fin.Get("midW_fit_wgt");
    PlotUtils::MnvH1D* BSigBen_Hi_wgt = (PlotUtils::MnvH1D*)fin.Get("hiW_fit_wgt");
    PlotUtils::MnvH1D *BSigBenWSB_Low_posfit = BSigBenWSB_Low_prefit->Clone(uniq());
    PlotUtils::MnvH1D *BSigBenWSB_Mid_posfit = BSigBenWSB_Mid_prefit->Clone(uniq());
    PlotUtils::MnvH1D *BSigBenWSB_Hi_posfit = BSigBenWSB_Hi_prefit->Clone(uniq());
    BSigBenWSB_Low_posfit->Scale(BSigBen_Low_wgt->GetBinContent(1)); 
    BSigBenWSB_Mid_posfit->Scale(BSigBen_Mid_wgt->GetBinContent(1)); 
    BSigBenWSB_Hi_posfit->Scale(BSigBen_Hi_wgt->GetBinContent(1)); 

    // Plotting Sidebands 
    // Plotting the comparison with Aaron file
    PlotRatio(ASigBenWSB_data, Aaron_SB_data_prefit,
              "W_SB", 1/(data_POT_A/utilASig.m_data_pot) ,
              "data_difMacros", false, true,"BenMacro/Aaron'sMacro", "W_{exp} MeV");
    PlotRatio(ASigBenWSB_Low_prefit, Aaron_SB_true_Low_prefit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "LowW_prefit_diffMacros", false, true,"BenMacro/Aaron'sMacro", "W_{exp} MeV");
    PlotRatio(ASigBenWSB_Low_posfit, Aaron_SB_true_Low_posfit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "LowW_posfit_diffMacros", false, true,"BenMacro/Aaron'sMacro", "W_{exp} MeV");
    PlotRatio(ASigBenWSB_Mid_prefit, Aaron_SB_true_Mid_prefit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "MidW_prefit_diffMacros", false, true,"BenMacro/Aaron'sMacro", "W_{exp} MeV");
    PlotRatio(ASigBenWSB_Mid_posfit, Aaron_SB_true_Mid_posfit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "MidW_posfit_diffMacros", false, true,"BenMacro/Aaron'sMacro", "W_{exp} MeV");
    PlotRatio(ASigBenWSB_Hi_prefit, Aaron_SB_true_Hi_prefit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "HiW_prefit_diffMacro", false, true,"BenMacro/Aaron'sMacro", "W_{exp} MeV");
    PlotRatio(ASigBenWSB_Hi_posfit, Aaron_SB_true_Hi_posfit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "HiW_posfit_diffMacro", false, true,"BenMacro/Aaron'sMacro", "W_{exp} MeV");
    PlotFittedW(ASigBenWSB_Signal, ASigBenWSB_Low_prefit, ASigBenWSB_Mid_prefit,
                ASigBenWSB_Hi_prefit, ASigBenWSB_data, utilASig.m_data_pot, utilASig.m_mc_pot,
                utilASig.m_signal_definition, true, "BenPrefit_ASig"); 
    PlotFittedW(ASigBenWSB_Signal, ASigBenWSB_Low_posfit, ASigBenWSB_Mid_posfit,
                ASigBenWSB_Hi_posfit, ASigBenWSB_data, utilASig.m_data_pot, utilASig.m_mc_pot,
                utilASig.m_signal_definition, false, "BenPosfit_ASig"); 
    PlotFittedW(Aaron_SB_Signal, Aaron_SB_true_Low_prefit, Aaron_SB_true_Mid_prefit,
                Aaron_SB_true_Hi_prefit, Aaron_SB_data_prefit, 1., 1.,
                utilASig.m_signal_definition, true, "AaronPrefit_ASig"); 
    PlotFittedW(Aaron_SB_Signal, Aaron_SB_true_Low_posfit, Aaron_SB_true_Mid_posfit,
                Aaron_SB_true_Hi_posfit, Aaron_SB_data_prefit, 1., 1.,
                utilASig.m_signal_definition, false, "AaronPosfit_ASig"); 

    //Comparing Ben's Macro AaronSignal deffinition and Ben's Signal Deffinition
    PlotRatio(BSigBenWSB_data, ASigBenWSB_data,
              "W_SB", 1/(data_POT_A/utilASig.m_data_pot) ,
              "data_BenMacro", false, true, "BenSignal/AaronSignal", "W_{exp} MeV");
    PlotRatio(BSigBenWSB_Low_prefit, ASigBenWSB_Low_prefit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "LowW_prefit_BenMacro", false, true,"BenSignal/AaronSignal", "W_{exp} MeV");
    PlotRatio(BSigBenWSB_Low_posfit, ASigBenWSB_Low_posfit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "LowW_posfit_BenMacro", false, true,"BenSignal/AaronSignal", "W_{exp} MeV");
    PlotRatio(BSigBenWSB_Mid_prefit, ASigBenWSB_Mid_prefit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "MidW_prefit_BenMacro", false, true,"BenSignal/AaronSignal", "W_{exp} MeV");
    PlotRatio(BSigBenWSB_Mid_posfit, ASigBenWSB_Mid_posfit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "MidW_posfit_BenMacro", false, true,"BenSignal/AaronSignal", "W_{exp} MeV");
    PlotRatio(BSigBenWSB_Hi_prefit, ASigBenWSB_Hi_prefit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "HiW_prefit_BenMacro", false, true,"BenSignal/AaronSignal", "W_{exp} MeV");
    PlotRatio(BSigBenWSB_Hi_posfit, ASigBenWSB_Hi_posfit,
              "W_SB", 1/(data_POT_A/utilASig.m_mc_pot) ,
              "HiW_posfit_BenMacro", false, true,"BenSignal/AaronSignal", "W_{exp} MeV");
    PlotFittedW(BSigBenWSB_Signal, BSigBenWSB_Low_prefit, BSigBenWSB_Mid_prefit,
                BSigBenWSB_Hi_prefit, BSigBenWSB_data, util.m_data_pot, util.m_mc_pot,
                util.m_signal_definition, true, "BenPrefit_BenSig"); 
    PlotFittedW(BSigBenWSB_Signal, BSigBenWSB_Low_posfit, BSigBenWSB_Mid_posfit,
                BSigBenWSB_Hi_posfit, BSigBenWSB_data, util.m_data_pot, util.m_mc_pot,
                util.m_signal_definition, false, "BenPosfit_BenSig"); 


      for (auto v : variables) {
      std::string var = v->Name();
      std::string Aaronvar;
      int nbins = CCPi::GetBinning(var).GetSize()-1;  
      const Double_t* bins = CCPi::GetBinning(var).GetArray();
      if (var == "q2")
        Aaronvar = "q2";
      else if (var == "thetamu_deg")
        Aaronvar = "muon_theta";
      else if (var == "thetapi_deg")
        Aaronvar = "pion_theta";
      else if (var == "tpi")
        Aaronvar = "pion_ekin";
      else if (var == "wexp"){
        Aaronvar = "W";
        niter = "10";
      }
      else if (var == "pmu"){
        Aaronvar = "muon_p";
        niter = "1";
      }
      else if (var == "pzmu"){
        Aaronvar = "muon_pz";
        niter = "1";
      }
      else if (var == "ptmu")
        Aaronvar = "muon_pt";
      else continue;
      // Getting Aaron's histograms
      PlotUtils::MnvH1D *EffDen_Aaron_aux = (PlotUtils::MnvH1D*)finCCPi.Get(Form("h_%s_plastic_EFF_DEN_pi_channel_mc",Aaronvar.c_str()));
      PlotUtils::MnvH1D *EffNum_Aaron_aux = (PlotUtils::MnvH1D*)finCCPi.Get(Form("h_%s_plastic_EFF_NUM_pi_channel_mc", Aaronvar.c_str()));
      PlotUtils::MnvH1D *xsec_Aaron_Aux = (PlotUtils::MnvH1D*)fAaronxSec.Get(Form("h_%s_plastic_pi_channel_mc_xsec_nucleon", Aaronvar.c_str()));
      PlotUtils::MnvH1D *xsec_Aaron_data_Aux = (PlotUtils::MnvH1D*)fAaronxSec.Get(Form("h_%s_plastic_data_xsec_nucleon", Aaronvar.c_str()));
      PlotUtils::MnvH1D *BGs_Aaron_Data_aux = (PlotUtils::MnvH1D*)fAaronBGs.Get(Form("h_%s_plastic_BkgdSubData", Aaronvar.c_str()));
      PlotUtils::MnvH1D *Sel_Aaron_Data_aux = (PlotUtils::MnvH1D*)fAaronBGs.Get(Form("h_%s_plastic_data", Aaronvar.c_str()));
      PlotUtils::MnvH1D *Sel_Aaron_MC_aux = (PlotUtils::MnvH1D*)fAaronBGs.Get(Form("h_%s_plastic_pi_channel_mc", Aaronvar.c_str()));
      PlotUtils::MnvH1D *Unfold_Aaron_Data_aux = (PlotUtils::MnvH1D*)fAaronBGs.Get(Form("h_%s_plastic_dataUnfold%s", Aaronvar.c_str(), niter.c_str()));
      PlotUtils::MnvH1D *Eff_Aaron_aux = (PlotUtils::MnvH1D*)fAaronBGs.Get(Form("h_%s_plastic_pi_channel_efficiency", Aaronvar.c_str()));
      PlotUtils::MnvH1D *Aaron_EFFCorr_data_aux = (PlotUtils::MnvH1D*)fAaronBGs.Get(Form("h_%s_plastic_data_eff_corrected", Aaronvar.c_str()));

      PlotUtils::MnvH1D *EffDen_Aaron = new PlotUtils::MnvH1D(Form("%s_truth_sig_Aaron", var.c_str()),Form("%s_truth_sig_Aaron", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *EffNum_Aaron = new PlotUtils::MnvH1D(Form("%s_BGs_Aaron", var.c_str()),Form("%s_BGs_Aaron", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *xsec_Aaron = new PlotUtils::MnvH1D(Form("%s_xsec_Aaron_ALL", var.c_str()),Form("%s_xsec_Aaron_ALL", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *xsec_Aaron_data = new PlotUtils::MnvH1D(Form("%s_xsec_Aaron_data_ALL", var.c_str()),Form("%s_xsec_Aaron_data", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *BGs_Aaron_Data = new PlotUtils::MnvH1D(Form("%s_BGs_Aaron_Data", var.c_str()),Form("%s_BGs_Aaron_Data", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *Sel_Aaron_Data = new PlotUtils::MnvH1D(Form("%s_Sel_Aaron_Data", var.c_str()),Form("%s_Sel_Aaron_Data", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *Sel_Aaron_MC = new PlotUtils::MnvH1D(Form("%s_Sel_Aaron_MC", var.c_str()),Form("%s_Sel_Aaron_MC", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *Unfold_Aaron_Data = new PlotUtils::MnvH1D(Form("%s_Unfold_Aaron_Data", var.c_str()),Form("%s_Unfold_Aaron_Data", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *Eff_Aaron = new PlotUtils::MnvH1D(Form("%s_Eff_Aaron", var.c_str()),Form("%s_Eff_Aaron", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());
      PlotUtils::MnvH1D *Aaron_EFFCorr_data = new PlotUtils::MnvH1D(Form("%s_Eff_corrected_Aaron", var.c_str()),Form("%s_Eff_corrected_Aaron", var.c_str()),nbins, CCPi::GetBinning(var).GetArray());

      for (int i = 1; i <= nbins; ++i){
        double EffNum, BGs_data, EffDen, xsec_mc, xsec_data, 
               Sel_data, Sel_mc, Unfold_data, Eff, EFFCorr_data;
        EffNum = EffNum_Aaron_aux->GetBinContent(i);
        BGs_data = BGs_Aaron_Data_aux->GetBinContent(i);
        EffDen = EffDen_Aaron_aux->GetBinContent(i);
        xsec_mc = xsec_Aaron_Aux->GetBinContent(i);
        xsec_data = xsec_Aaron_data_Aux->GetBinContent(i);
        Sel_data = Sel_Aaron_Data_aux->GetBinContent(i);
        Sel_mc = Sel_Aaron_MC_aux->GetBinContent(i);
        Unfold_data = Unfold_Aaron_Data_aux->GetBinContent(i);
        Eff = EffNum/EffDen;
        EFFCorr_data = Aaron_EFFCorr_data_aux->GetBinContent(i);

        EffDen_Aaron->SetBinContent(i, EffDen);
        EffNum_Aaron->SetBinContent(i, EffNum);
	BGs_Aaron_Data->SetBinContent(i, BGs_data);
        xsec_Aaron->SetBinContent(i, xsec_mc);
        xsec_Aaron_data->SetBinContent(i, xsec_data);
        Sel_Aaron_Data->SetBinContent(i, Sel_data);
        Sel_Aaron_MC->SetBinContent(i, Sel_mc);
        Unfold_Aaron_Data->SetBinContent(i, Unfold_data);
        Eff_Aaron->SetBinContent(i, Eff);
        Aaron_EFFCorr_data->SetBinContent(i,EFFCorr_data);
      }

      PlotUtils::MnvH1D* ASigBenEffdenMC = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("effden_%s_true", var.c_str()));

      PlotUtils::MnvH1D* ASigBenBGsMC = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("effnum_%s_true", var.c_str()));
      PlotUtils::MnvH1D* ASigBenBGsData = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("bg_subbed_data_%s", var.c_str()));

      PlotUtils::MnvH1D* ASigBenXSecMC = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("mc_cross_section_%s", var.c_str()));
      PlotUtils::MnvH1D* ASigBenXSecdata = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("cross_section_%s", var.c_str()));
      PlotUtils::MnvH1D* ASigBenSeldata = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("selection_data_%s", var.c_str()));
      PlotUtils::MnvH1D* ASigBenSelMC = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("selection_mc_%s", var.c_str()));
      PlotUtils::MnvH1D* ASigBenUnfold = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("unfolded_%s", var.c_str()));
      PlotUtils::MnvH1D* ASigBenEff = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("efficiency_%s", var.c_str()));
      PlotUtils::MnvH1D* ASigBenEffCorr = (PlotUtils::MnvH1D*)fAaronSigBenMacro.Get(Form("efficiency_corrected_data_%s", var.c_str()));


      PlotUtils::MnvH1D* BSigBenEffdenMC = (PlotUtils::MnvH1D*)fin.Get(Form("effden_%s_true", var.c_str()));
      PlotUtils::MnvH1D* BSigBenBGsMC = (PlotUtils::MnvH1D*)fin.Get(Form("effnum_%s_true", var.c_str()));
      PlotUtils::MnvH1D* BSigBenBGsData = (PlotUtils::MnvH1D*)fin.Get(Form("bg_subbed_data_%s", var.c_str()));
      PlotUtils::MnvH1D* BSigBenXSecMC = (PlotUtils::MnvH1D*)fin.Get(Form("mc_cross_section_%s", var.c_str()));
      PlotUtils::MnvH1D* BSigBenXSecdata = (PlotUtils::MnvH1D*)fin.Get(Form("cross_section_%s", var.c_str()));
      PlotUtils::MnvH1D* BSigBenSeldata = (PlotUtils::MnvH1D*)fin.Get(Form("selection_data_%s", var.c_str()));
      PlotUtils::MnvH1D* BSigBenSelMC = (PlotUtils::MnvH1D*)fin.Get(Form("selection_mc_%s", var.c_str()));
      PlotUtils::MnvH1D* BSigBenUnfold = (PlotUtils::MnvH1D*)fin.Get(Form("unfolded_%s", var.c_str()));
      PlotUtils::MnvH1D* BSigBenEff = (PlotUtils::MnvH1D*)fin.Get(Form("efficiency_%s", var.c_str()));
      PlotUtils::MnvH1D* BSigBenEffCorr = (PlotUtils::MnvH1D*)fin.Get(Form("efficiency_corrected_data_%s", var.c_str()));

      std::cout << "Aaron MC POT = " << MC_POT_A <<"\n";
      std::cout << "Ben Sig Def MC POT = " << util.m_mc_pot <<"\n";
      std::cout << "Ratio Ben Sig Def= " << MC_POT_A/util.m_mc_pot <<"\n";
      std::cout << "Ben Macro Sig Def POT = " << utilASig.m_mc_pot <<"\n";
      std::cout << "Ratio Aaron Sig Def = " << util.m_mc_pot/utilASig.m_mc_pot <<"\n";
      PlotRatio(ASigBenXSecdata, xsec_Aaron_data, var, 1., "data_xSec_file_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(ASigBenEffdenMC, EffDen_Aaron, var, 1/(MC_POT_A/utilASig.m_mc_pot), "mc_Effden_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(ASigBenBGsMC, EffNum_Aaron, var, 1/(MC_POT_A/utilASig.m_mc_pot), "mc_EffNum_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(ASigBenXSecMC, xsec_Aaron, var, 1. , "mc_xSec_file_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");    
      PlotRatio(ASigBenBGsData, BGs_Aaron_Data, var, 1/(data_POT_A/utilASig.m_data_pot) , "BGs_data_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(ASigBenSeldata, Sel_Aaron_Data, var, 1/(data_POT_A/utilASig.m_data_pot), "Sel_data_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(ASigBenSelMC, Sel_Aaron_MC, var, 1/(data_POT_A/utilASig.m_mc_pot), "Sel_MC_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(ASigBenUnfold, Unfold_Aaron_Data, var, 1/(data_POT_A/utilASig.m_data_pot), "Unfold_data_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");
      PlotRatio(ASigBenEff, Eff_Aaron, var, 1., "Efficiency_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(ASigBenEffCorr, Aaron_EFFCorr_data, var, 1/(data_POT_A/utilASig.m_data_pot), "Efficiency_Correction_data_diffMacro", false, true,"BenMacro/Aaron'sMacro", v->m_hists.m_xlabel + " (" + v->m_units + ")");


      PlotRatio(BSigBenXSecdata, ASigBenXSecdata, var, 1., "data_xSec_file_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(BSigBenEffdenMC, ASigBenEffdenMC, var, 1/(utilASig.m_mc_pot/util.m_mc_pot), "mc_Effden_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(BSigBenBGsMC, ASigBenBGsMC, var, 1/(utilASig.m_mc_pot/util.m_mc_pot), "mc_EffNum_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(BSigBenXSecMC, ASigBenXSecMC, var, 1. , "mc_xSec_file_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");    
      PlotRatio(BSigBenBGsData, ASigBenBGsData, var, 1/(utilASig.m_data_pot/util.m_data_pot) , "BGs_data_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(BSigBenSeldata, ASigBenSeldata, var, 1/(utilASig.m_data_pot/util.m_data_pot), "Sel_data_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(BSigBenSelMC, ASigBenSelMC, var, 1/(utilASig.m_mc_pot/util.m_mc_pot), "Sel_MC_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(BSigBenUnfold, ASigBenUnfold, var, 1/(utilASig.m_data_pot/util.m_data_pot), "Unfold_data_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");
      PlotRatio(BSigBenEff, ASigBenEff, var, 1., "Efficiency_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");  
      PlotRatio(BSigBenEffCorr, ASigBenEffCorr, var, 1/(utilASig.m_data_pot/util.m_data_pot), "Efficiency_Correction_data_OnlyBenMacro", false, true,"BenSig/AaronSig", v->m_hists.m_xlabel + " (" + v->m_units + ")");

    }  
  }

  // PLOT Event Selection, BGs (error)
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    bool do_cov_area_norm = false;
    for (auto var : variables) {
      if (var->Name() ==  "wexp_fit") continue;
      std::cout << var->Name() << "\n";
      do_cov_area_norm = false;
      Plotter plot_info(var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      bool do_bin_width_norm = true, do_log_scale = false, do_bg = true;
      bool do_tuned_bg = false;

      // selection and BG error before tuning
      if (!var->m_is_true) {
        PlotVar_Selection(plot_info, -1., do_log_scale, do_bg, do_tuned_bg,
                          do_bin_width_norm);
        if (plot_errors) PlotVar_ErrorSummary(plot_info);
        if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

        // selection and BG error after tuning
        do_tuned_bg = true;
        PlotVar_Selection(plot_info, -1., do_log_scale, do_bg, do_tuned_bg,
                          do_bin_width_norm);
        if (plot_errors) PlotVar_ErrorSummary(plot_info);
        if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);
      }
    }
  }

  // PLOT Efficiency & Migration
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;

    for (auto var : variables) {
      if (var->Name() ==  "wexp_fit") continue;
      // var->LoadMCHistsFromFile(fin, util.m_error_bands);

      const Plotter plot_info(
          var, util.m_mc_pot, util.m_data_pot, do_frac_unc, do_cov_area_norm,
          include_stat, util.m_signal_definition);

      // Efficiency
      if (var->m_is_true) {
        PlotUtils::MnvH1D* eff =
            (PlotUtils::MnvH1D*)var->m_hists.m_effnum.hist->Clone(uniq());
        PlotUtils::MnvH1D* effnum =
            (PlotUtils::MnvH1D*)var->m_hists.m_effnum.hist->Clone(uniq());
        PlotUtils::MnvH1D* effden =
            (PlotUtils::MnvH1D*)var->m_hists.m_effden.hist->Clone(uniq());
        eff->Divide(effnum, effden);
        var->m_hists.m_efficiency = (PlotUtils::MnvH1D*)eff->Clone(uniq());

        double ymax = -1.;
        bool do_log_scale = false;
        bool do_bg = true;
        bool do_tuned_bg = true;
        PlotMC(eff, plot_info, Form("Efficiency_%s", var->Name().c_str()),
               0.075, "Efficiency");
        if (plot_errors) PlotEfficiency_ErrorSummary(plot_info);
      }

      // Migration
      if (!var->m_is_true) {
        PlotUtils::MnvH2D* mig =
            (PlotUtils::MnvH2D*)var->m_hists.m_migration.hist->Clone(uniq());
        PlotMigration_AbsoluteBins(mig, var->Name());
        PlotMigration_VariableBins(mig, var->Name());
      }
    }
  }

  // PLOT Background Subtraction
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    for (auto var : variables) {
      if (var->Name() ==  "wexp_fit") continue;
      if (var->m_is_true) continue;

      Plotter plot_info(var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      double ymax = -1.;
      bool do_log_scale = false;
      bool do_bg = true;
      bool do_tuned_bg = true;
      bool do_bin_width_norm = true;
      if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

      do_tuned_bg = false;
      if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

      Plot_BGSub(plot_info, ".", ymax, do_log_scale, do_bin_width_norm);

      if (plot_errors) PlotBGSub_ErrorSummary(plot_info);
    }
  }

  // PLOT W Sideband Fit
  if (true) {
    const bool do_frac_unc = true;
    const bool do_cov_area_norm = false;
    const bool include_stat = true;

    Plotter plot_info(util.m_mc_pot, util.m_data_pot,
                                     do_frac_unc, do_cov_area_norm,
                                     include_stat, util.m_signal_definition);

    PlotUtils::MnvH1D* loW_fit_wgt = (PlotUtils::MnvH1D*)fin.Get("loW_fit_wgt");
    PlotUtils::MnvH1D* midW_fit_wgt =
        (PlotUtils::MnvH1D*)fin.Get("midW_fit_wgt");
    PlotUtils::MnvH1D* hiW_fit_wgt = (PlotUtils::MnvH1D*)fin.Get("hiW_fit_wgt");

    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, loW_fit_wgt, "WSidebandFit_loW");

    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, midW_fit_wgt,
                                    "WSidebandFit_midW");
    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, hiW_fit_wgt, "WSidebandFit_hiW");

    // Wexp distribution, stacked, all cuts except for W, pre-fit
    std::string tag;
    double ymax = -1;
    Variable* var = GetVar(variables, sidebands::kFitVarString);
    PlotWSidebandStacked(var, var->m_hists.m_wsideband_data,
                         var->GetStackArray(static_cast<WSidebandType>(0)),
                         util.m_data_pot, util.m_mc_pot,
                         util.m_signal_definition, tag, ymax);

    // TODO plot pre/postfit
    // for (auto var : variables) {
    //  tag = "SidebandRegion";
    //  bool do_prefit = true;
    //  bool do_bin_width_norm = true;
    //  CVUniverse* universe = util.m_error_bands.at("cv").at(0);
    //  PlotFittedW(var, *universe, hw_loW_fit_wgt, hw_midW_fit_wgt,
    //              hw_hiW_fit_wgt, util.m_data_pot, util.m_mc_pot,
    //              util.m_signal_definition, do_prefit, tag, ymax,
    //              do_bin_width_norm);
    //  do_prefit = false;
    //  PlotFittedW(var, *universe, hw_loW_fit_wgt, hw_midW_fit_wgt,
    //              hw_hiW_fit_wgt, util.m_data_pot, util.m_mc_pot,
    //              util.m_signal_definition, do_prefit, tag, ymax,
    //              do_bin_width_norm);
    //}
  }

  // PLOT unfolded
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    const double ymax = -1.;
    const bool do_log_scale = false;
    const bool do_bin_width_norm = true;
    for (auto reco_var : variables) {
      if (reco_var->Name() ==  "wexp_fit") continue;
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      Plot_Unfolded(plot_info, reco_var->m_hists.m_unfolded,
                    true_var->m_hists.m_effnum.hist);
      if (plot_errors) PlotUnfolded_ErrorSummary(plot_info);
    }
  }

  // PLOT cross section
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    const double ymax = -1.;
    const bool do_log_scale = false;
    const bool do_bin_width_norm = true;
    for (auto reco_var : variables) {
      if (reco_var->Name() ==  "wexp_fit") continue;
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot,
                                       do_frac_unc, do_cov_area_norm,
                                       include_stat, util.m_signal_definition);

      PlotUtils::MnvH1D* m_mc_cross_section = (PlotUtils::MnvH1D*)fin.Get(
          Form("mc_cross_section_%s", reco_var->Name().c_str()));

      // std::vector<std::string> x_bands =
      // reco_var->m_hists.m_cross_section->GetVertErrorBandNames();
      // if(reco_var->Name() == "ptmu")
      //  for (auto s : x_bands) std::cout << s << "\n";

      // std::cout << reco_var->Name() << "\n";

      Plot_CrossSection(plot_info, reco_var->m_hists.m_cross_section,
                        m_mc_cross_section);
      if (plot_errors)
        PlotCrossSection_ErrorSummary(
            plot_info);  // Adds chi2 label and prints out assumed binning.

      // Various error matrices
      // TMatrixD
      // scaled_mtx(reco_var->m_hists.m_cross_section->GetTotalErrorMatrix());
      // TMatrixD scaled_mtx(m_mc_cross_section->GetStatErrorMatrix());

      // TMatrixD
      // data_stat_err_mtx(reco_var->m_hists.m_cross_section->GetStatErrorMatrix());
      // TMatrixD
      // data_tot_err_mtx(reco_var->m_hists.m_cross_section->GetTotalErrorMatrix());
      // if (plot_errors) PlotMatrix(data_stat_err_mtx, reco_var->Name(),
      // "Data_Stat_Err_Mtx"); if (plot_errors) PlotMatrix(data_tot_err_mtx,
      // reco_var->Name(), "Data_Tot_Err_Mtx");

      // TMatrixD mc_stat_err_mtx(m_mc_cross_section->GetStatErrorMatrix());
      // TMatrixD mc_tot_err_mtx(m_mc_cross_section->GetTotalErrorMatrix());
      // if (plot_errors) PlotMatrix(mc_stat_err_mtx,   reco_var->Name(),
      // "MC_Stat_Err_Mtx"); if (plot_errors) PlotMatrix(mc_tot_err_mtx,
      // reco_var->Name(), "MC_Tot_Err_Mtx");

      // TMatrixD
      // data_sys_err_mtx(reco_var->m_hists.m_cross_section->GetSysErrorMatrix());
      // TMatrixD
      // data_sys_corr_mtx(reco_var->m_hists.m_cross_section->GetSysCorrelationMatrix
      // ()); if (plot_errors) PlotMatrix(data_sys_err_mtx,  reco_var->Name(),
      // "Data_Sys_Err_Mtx"); if (plot_errors) PlotMatrix(data_sys_cor_mtx,
      // reco_var->Name(), "Data_Sys_Cor_Mtx");

      // TMatrixD mc_sys_err_mtx (m_mc_cross_section->GetSysErrorMatrix());
      // TMatrixD mc_sys_cor_mtx(m_mc_cross_section->GetSysCorrelationMatrix());
      // if (plot_errors) PlotMatrix(mc_sys_err_mtx,    reco_var->Name(),
      // "MC_Sys_Err_Mtx"); if (plot_errors) PlotMatrix(mc_sys_cor_mtx,
      // reco_var->Name(), "MC_Sys_Cor_Mtx");

      // if (plot_errors) PlotCovarianceMatrix(scaled_mtx, reco_var->Name(),
      // "Data_StatErrMtx"); for(int c = 0; c < scaled_mtx.GetNcols(); ++c ) {
      //  for(int r = 0; r < scaled_mtx.GetNrows(); ++r ) {
      //    std::cout << scaled_mtx(r, c) << "  ";
      //  }
      //  std::cout << "\n";
      //}

      // for(int c = 0; c < scaled_mtx2.GetNcols(); ++c ) {
      //  for(int r = 0; r < scaled_mtx2.GetNrows(); ++r ) {
      //    std::cout << scaled_mtx2(r, c) << "  ";
      //  }
      //  std::cout << "\n";
      //}

      // std::cout << "data (" << scaled_mtx.GetNrows() << "," <<
      // scaled_mtx.GetNcols() << ")\n"; std::cout << "mc (" <<
      // scaled_mtx2.GetNrows() << "," << scaled_mtx2.GetNcols() << ")\n";

      // PrintChi2Info(plot_info, reco_var->m_hists.m_cross_section,
      // m_mc_cross_section); // this one works
    }
  }

  //============================================================================
}

#endif  // plotCrossSectionFromFile_C

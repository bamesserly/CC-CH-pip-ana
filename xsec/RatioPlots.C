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
// Main
//==============================================================================
void RatioPlots() {
  // Infiles
  TFile fBenEve("DataXSecInputs_20240327_ALL_mixed_AaronBinning_nosys_p4.root", "READ");
  TFile fAaronGENIE("GENIEXSECEXTRACT_AarSignalDefMCME1A_q2.root", "READ");
  TFile fAaronNoLowTpiGENIE("GENIEXSECEXTRACT_AarSignalDefTpichangeMCME1A_q2.root", "READ");
  TFile fBenEveGENIE("GENIEXSECEXTRACT_BenSignalDefMCME1A_q2.root", "READ");
  TFile fAaronXsec("/minerva/data/users/abercell/hists/xsec/xsec_new_jeffrey_flux_MENU1PI_plastic_MinervaME1ABCDEFGLMNOP.root", "READ");
  TFile fAaronXsecIn("/minerva/data/users/abercell/hists/xsec_inputs/Merge_BkgdSub_Unfold_MENU1PI_POTNorm_plastic_MinervaME1ABCDEFGLMNOP.root", "READ");
//  TFile fAaronEff("/minerva/data/users/abercell/hists/Macro/GridOneLoop_MENU1PI_MinosMatched_plastic_Merged_NewdEdXCal_MinervaME1ABCDEFGLMNOP_Data_Merged_NewdEdXCal_Tracker_MinervaME1ABCDEFGLMNOP_MC.root", "READ");
  TFile fAaronEff("/minerva/data/users/abercell/hists/xsec_inputs/Merge_Eff_MENU1PI_POTNorm_plastic_MinervaME1ABCDEFGLMNOP.root", "READ");
  TFile fBenEveMAaronSD("DataXSecInputs_20240318_ALL_AaronSignalDef_nosys_p4.root", "READ");
  TFile fBenEveSDNoTpi("DataXSecInputs_20240325_ALL_mixed_NoTpiweight_nosys_p4.root", "READ");
//  TFile fBenEveSDNoTpi("DataXSecInputs_20240304_ALL_AaronSigDef_plusAaronFidVolCuts_nosys_p4.root", "READ");
  TFile fin8("/minerva/data/users/abercell/hists/Macro/GridOneLoop_MENU1PI_MinosMatched_plastic_Merged_NewdEdXCal_MinervaME1ABCDEFGLMNOP_Data_Merged_NewdEdXCal_Tracker_MinervaME1ABCDEFGLMNOP_MC.root", "READ");
  TFile fBenEveQ2("DataXSecInputs_20240325_ALL_mixed_AaronBreakdown_nosys_p4.root", "READ");
  std::vector<std::string> variables = {"q2", "mixtpi", "thetamu_deg"};	 

  PlotUtils::MnvH1D* dataPOTBenEve =
             (PlotUtils::MnvH1D*)fBenEve.Get("data_pot");
  PlotUtils::MnvH1D* MCPOTBenEve =
             (PlotUtils::MnvH1D*)fBenEve.Get("mc_pot");
  
  if (true){
    PlotUtils::MnvH1D* AaronDataSelQ2 =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get("h_q2_plastic_data");
    PlotUtils::MnvH1D* AaronBGSubQ2 =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get("h_q2_plastic_BkgdSubData");

    PlotUtils::MnvH1D* dummyBenDataSel_Q2 = (PlotUtils::MnvH1D*)fBenEveQ2.Get("selection_data_q2");
    PlotUtils::MnvH1D* dummyBenDataSel_Q2_Aaron = (PlotUtils::MnvH1D*)fBenEveQ2.Get("selection_data_q2_Aaron");
    PlotUtils::MnvH1D* dummyBenBGSub_Q2 = (PlotUtils::MnvH1D*)fBenEveQ2.Get("bg_subbed_data_q2");
    PlotUtils::MnvH1D* dummyBenBGSub_Q2_Aaron = (PlotUtils::MnvH1D*)fBenEveQ2.Get("bg_subbed_data_q2_Aaron");
    PlotUtils::MnvH1D* BenDataSel_Q2 = (PlotUtils::MnvH1D*)AaronDataSelQ2->Clone("BenDataSel");
    BenDataSel_Q2->SetTitle("Our DataSel");
    PlotUtils::MnvH1D* BenDataSel_Q2_Aaron = (PlotUtils::MnvH1D*)AaronDataSelQ2->Clone("BenDataSel");
    BenDataSel_Q2_Aaron->SetTitle("Aaron's DataSel");
    PlotUtils::MnvH1D* BenBGSub_Q2 = (PlotUtils::MnvH1D*)AaronDataSelQ2->Clone("OurBGSub");
    BenBGSub_Q2->SetTitle("Our BGSub");
    PlotUtils::MnvH1D* BenBGSub_Q2_Aaron = (PlotUtils::MnvH1D*)AaronDataSelQ2->Clone("OurBGSub");
    BenBGSub_Q2_Aaron->SetTitle("Aaron's BGSub");

    BenDataSel_Q2->Reset();
    BenDataSel_Q2_Aaron->Reset();
    BenBGSub_Q2->Reset();
    BenBGSub_Q2_Aaron->Reset();

    std::vector<PlotUtils::MnvH1D*> hvecDataSelQ2;
    std::vector<PlotUtils::MnvH1D*> hvecBgSubQ2;

    for (int i = 1; i <= BenDataSel_Q2->GetNbinsX(); ++i){
      BenDataSel_Q2->SetBinContent(i,dummyBenDataSel_Q2->GetBinContent(i));
      BenDataSel_Q2_Aaron->SetBinContent(i,dummyBenDataSel_Q2_Aaron->GetBinContent(i));
      BenBGSub_Q2->SetBinContent(i,dummyBenBGSub_Q2->GetBinContent(i));
      BenBGSub_Q2_Aaron->SetBinContent(i,dummyBenBGSub_Q2_Aaron->GetBinContent(i));
    } 

    hvecDataSelQ2.push_back(BenDataSel_Q2);
    hvecDataSelQ2.push_back(BenDataSel_Q2_Aaron);
    hvecBgSubQ2.push_back(BenBGSub_Q2);
    hvecBgSubQ2.push_back(BenBGSub_Q2_Aaron);

    PlotRatioVec(hvecDataSelQ2, AaronDataSelQ2, "q2", 1.,"DataSelection_Bkd",false, 
             "BenSig/Aaronsig", "Q^{2} (GeV^{2})"); 
    PlotRatioVec(hvecDataSelQ2, AaronDataSelQ2, "q2", 1.,"DataSelection_Bkd",false, 
             "BenSig/Aaronsig", "Q^{2} (GeV^{2})"); 
    PlotRatioVec(hvecBgSubQ2, AaronBGSubQ2, "q2", 1.,"BGSub_Bkd",false, 
             "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  }


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
    PlotUtils::MnvH1D* AaronUnfoldData =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get(Form("h_%s_plastic_dataUnfold2", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronCorrectedData =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get(Form("h_%s_plastic_data_eff_corrected", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronEffDen =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get(Form("h_%s_plastic_pi_channel_truth_sig", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronXsecMC =
        (PlotUtils::MnvH1D*)fAaronXsec.Get(Form("h_%s_plastic_pi_channel_mc_xsec_nucleon", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronXsecData =
        (PlotUtils::MnvH1D*)fAaronXsec.Get(Form("h_%s_plastic_data_xsec_nucleon", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronEffNum =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get(Form("h_%s_plastic_pi_channel_signal", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronEff =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get(Form("h_%s_plastic_pi_channel_efficiency", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronNewEff =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get(Form("h_%s_plastic_pi_channel_efficiency", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronDataSel =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get(Form("h_%s_plastic_data", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronBGSub =
        (PlotUtils::MnvH1D*)fAaronXsecIn.Get(Form("h_%s_plastic_BkgdSubData", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronEffOtherFile =
        (PlotUtils::MnvH1D*)fAaronEff.Get(Form("h_%s_plastic_pi_channel_efficiency", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronEffNumOtherFile =
        (PlotUtils::MnvH1D*)fAaronEff.Get(Form("h_%s_plastic_EFF_NUM_pi_channel", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronEffDenOtherFile =
  //               (PlotUtils::MnvH1D*)fAaronXsec.Get("h_%s_plastic_pi_channel_truth_sig_ME1A");
        (PlotUtils::MnvH1D*)fAaronEff.Get(Form("h_%s_plastic_EFF_DEN_pi_channel", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronFluxNormalizer =
        (PlotUtils::MnvH1D*)fAaronXsec.Get(Form("h_%s_plastic_data_xsec_nucleon", Aaronvar.c_str()));
    PlotUtils::MnvH1D* AaronRevUnfoldData =
        (PlotUtils::MnvH1D*)fAaronXsec.Get(Form("h_%s_plastic_data_xsec_nucleon", Aaronvar.c_str()));
    PlotUtils::MnvH1D* h_data_POT =
        (PlotUtils::MnvH1D*)fAaronXsec.Get("h_Data_POT");
    PlotUtils::MnvH1D* h_MC_POT =
        (PlotUtils::MnvH1D*)fAaronXsec.Get("h_MC_POT");
    PlotUtils::MnvH1D* AaronNoNormXsecData =
        (PlotUtils::MnvH1D*)fAaronXsec.Get(Form("h_%s_plastic_data_xsec_nucleon", Aaronvar.c_str()));
    std::cout << "PASS 1\n";

    PlotUtils::MnvH1D* dummyBenXsecData = (PlotUtils::MnvH1D*)fBenEve.Get(Form("cross_section_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenEffCorrData = (PlotUtils::MnvH1D*)fBenEve.Get(Form("efficiency_corrected_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenEffDen = (PlotUtils::MnvH1D*)fBenEve.Get(Form("effden_%s_true", var.c_str()));
    PlotUtils::MnvH1D* dummyBenUnfoldData = (PlotUtils::MnvH1D*)fBenEve.Get(Form("unfolded_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenXsecMC = (PlotUtils::MnvH1D*)fBenEve.Get(Form("mc_cross_section_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenEff = (PlotUtils::MnvH1D*)fBenEve.Get(Form("efficiency_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenEffNum = (PlotUtils::MnvH1D*)fBenEve.Get(Form("effnum_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenDataSel = (PlotUtils::MnvH1D*)fBenEve.Get(Form("selection_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenBGSub = (PlotUtils::MnvH1D*)fBenEve.Get(Form("bg_subbed_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenXsecAaronSigData = (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get(Form("cross_section_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenXsecNoTpiWegihtData = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get(Form("cross_section_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigCutsEffDen = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get(Form("effden_%s_true", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigCutsEffNum = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get(Form("effnum_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigCutsEff = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get(Form("efficiency_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigCutsEffCorr = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get(Form("efficiency_corrected_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigCutsUnfold = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get(Form("unfolded_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigCutsDataSel = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get(Form("selection_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigCutsBGSub = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get(Form("bg_subbed_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigEffDen = (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get(Form("effden_%s_true", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigEffNum = (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get(Form("effnum_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigEff = (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get(Form("efficiency_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigEffCorr = (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get(Form("efficiency_corrected_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigUnfold = (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get(Form("unfolded_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigDataSel = (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get(Form("selection_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dummyBenAaronSigBGSub = (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get(Form("bg_subbed_data_%s", var.c_str()));
    PlotUtils::MnvH1D* dataPOTBenAaronSig =
               (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get("data_pot");
    PlotUtils::MnvH1D* MCPOTBenAaronSig =
               (PlotUtils::MnvH1D*)fBenEveMAaronSD.Get("mc_pot");
    std::cout << "PASS 2\n";
  
    PlotUtils::MnvH1D* dataPOTBenAaronSigCuts = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get("data_pot");
    std::cout << "PASS 3\n";
    PlotUtils::MnvH1D* MCPOTBenAaronSigCuts = (PlotUtils::MnvH1D*)fBenEveSDNoTpi.Get("mc_pot");
    std::cout << "PASS 4\n";
  
    PlotUtils::MnvH1D* BenXsecData = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenXsecData_%s", var.c_str()));
    BenXsecData->SetTitle("BenXsecData");
    std::cout << "PASS 5\n";
    PlotUtils::MnvH1D* BenXsecDataCorr = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenXsecDataCorr_%s", var.c_str()));
    PlotUtils::MnvH1D* BenEffCorrData = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenXsecDataCorr_%s", var.c_str()));
    BenEffCorrData->SetTitle("BenEffCorr");
    PlotUtils::MnvH1D* BenEffDen = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenXsecData_%s", var.c_str()));
    BenEffDen->SetTitle("BenEffDen");
    PlotUtils::MnvH1D* BenUnfoldData = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenUnfold_%s", var.c_str()));
    BenUnfoldData->SetTitle("BenUnfold");
    PlotUtils::MnvH1D* BenDataSel = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenDataSel_%s", var.c_str()));
    BenDataSel->SetTitle("BenDataSel");
    PlotUtils::MnvH1D* BenBGSub = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenBGSub_%s", var.c_str()));
    BenBGSub->SetTitle("BenBGSub");
    PlotUtils::MnvH1D* BenXsecMC = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenXsecDataCorr_%s", var.c_str()));
    PlotUtils::MnvH1D* BenEff = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenEff_%s", var.c_str()));
    BenEff->SetTitle("BenEff");
    PlotUtils::MnvH1D* BenEffNum = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenEffNum_%s", var.c_str()));
    BenEffNum->SetTitle("BenEffNum");
    PlotUtils::MnvH1D* BenXsecAaronSigData = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenXsecAaronSigData_%s", var.c_str()));
    BenXsecAaronSigData->SetTitle("BenXsecAaronSigData");
    PlotUtils::MnvH1D* BenAaronSigEffDen = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigEffDen_%s", var.c_str()));
    BenAaronSigEffDen->SetTitle("BenAaronSigEffDen");
    PlotUtils::MnvH1D* BenAaronSigEffNum = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigEffNum_%s", var.c_str()));
    BenAaronSigEffNum->SetTitle("Ben w/ Aaron Sig Def");
    PlotUtils::MnvH1D* BenAaronSigEff = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigEff_%s", var.c_str()));
    BenAaronSigEff->SetTitle("BenAaronSigEff");
    PlotUtils::MnvH1D* BenAaronSigEffCorr = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigEffCorr_%s", var.c_str()));
    BenAaronSigEffCorr->SetTitle("BenAaronSigEffCorr");
    PlotUtils::MnvH1D* BenAaronSigUnfold = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigUnfold_%s", var.c_str()));
    BenAaronSigUnfold->SetTitle("Ben w/ Aaron Sig Def");
  
    PlotUtils::MnvH1D* BenAaronSigDataSel = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigDataSel_%s", var.c_str()));
    BenAaronSigDataSel->SetTitle("BenAaronSigDataSel");
    PlotUtils::MnvH1D* BenAaronSigBGSub = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigBGSub_%s", var.c_str()));
    BenAaronSigBGSub->SetTitle("BenAaronSigBGSub");
    PlotUtils::MnvH1D* BenXsecNoTpiWegihtData = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenXsecAaronSigCuts_%s", var.c_str()));
    BenXsecNoTpiWegihtData->SetTitle("Ben No T_{#pi} weight");
    PlotUtils::MnvH1D* BenAaronSigCutsEffDen = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigCutsEffDen_%s", var.c_str()));
    BenAaronSigCutsEffDen->SetTitle("BenAaronSigCutsEffDen");
    PlotUtils::MnvH1D* BenAaronSigCutsEffNum = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigCutsEffNum_%s", var.c_str()));
    BenAaronSigCutsEffNum->SetTitle("BenAaronSigCutsEffNum");
    PlotUtils::MnvH1D* BenAaronSigCutsEff = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigCutsEff_%s", var.c_str()));
    BenAaronSigCutsEff->SetTitle("B/E SDef NotpiW Eff");
    PlotUtils::MnvH1D* BenAaronSigCutsEffCorr = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigCutsEffCorr_%s", var.c_str()));
    BenAaronSigCutsEffCorr->SetTitle("BenAaronSigCutsEffCorr");
    PlotUtils::MnvH1D* BenAaronSigCutsUnfold = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigCutsUnfold_%s", var.c_str()));
    BenAaronSigCutsUnfold->SetTitle("BenAaronSigCutsUnfold");
  
    PlotUtils::MnvH1D* BenAaronSigCutsDataSel = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigCutsDataSel_%s", var.c_str()));
    BenAaronSigCutsDataSel->SetTitle("BenAaronSigCutsDataSel");
    PlotUtils::MnvH1D* BenAaronSigCutsBGSub = (PlotUtils::MnvH1D*)AaronXsecData->Clone(Form("BenAaronSigCutsBGSub_%s", var.c_str()));
    BenAaronSigCutsBGSub->SetTitle("BenAaronSigCutsBGSub");
  
    BenXsecData->Reset();
    BenXsecDataCorr->Reset();
    BenEffDen->Reset();
    BenUnfoldData->Reset();
    BenEffCorrData->Reset();
    BenXsecMC->Reset();
    BenEff->Reset();
    BenEffNum->Reset();
    BenDataSel->Reset();
    BenBGSub->Reset();
    AaronNewEff->Reset();
    BenXsecAaronSigData->Reset();
    BenXsecNoTpiWegihtData->Reset();
    BenAaronSigEffDen->Reset();
    BenAaronSigEffNum->Reset();
    BenAaronSigEff->Reset();
    BenAaronSigEffCorr->Reset();
    BenAaronSigUnfold->Reset();
    BenAaronSigDataSel->Reset();
    BenAaronSigBGSub->Reset();
    BenAaronSigCutsEffDen->Reset();
    BenAaronSigCutsEffNum->Reset();
    BenAaronSigCutsEff->Reset();
    BenAaronSigCutsEffCorr->Reset();
    BenAaronSigCutsUnfold->Reset();
    BenAaronSigCutsDataSel->Reset();
    BenAaronSigCutsBGSub->Reset();

    std::vector<PlotUtils::MnvH1D*> hvec;
    std::vector<PlotUtils::MnvH1D*> hvecEffNum;
    std::vector<PlotUtils::MnvH1D*> hvecEffDen;
    std::vector<PlotUtils::MnvH1D*> hvecEff;
    std::vector<PlotUtils::MnvH1D*> hvecEffCorr;
    std::vector<PlotUtils::MnvH1D*> hvecUnfold;
    std::vector<PlotUtils::MnvH1D*> hvecDataSel;
    std::vector<PlotUtils::MnvH1D*> hvecBGSub;
    
  
    for (int i = 1; i <= BenXsecData->GetNbinsX(); ++i){
      BenXsecData->SetBinContent(i,dummyBenXsecData->GetBinContent(i));
      BenEffCorrData->SetBinContent(i,dummyBenEffCorrData->GetBinContent(i));
      BenEffDen->SetBinContent(i,dummyBenEffDen->GetBinContent(i));
      BenDataSel->SetBinContent(i,dummyBenDataSel->GetBinContent(i));
      BenBGSub->SetBinContent(i,dummyBenBGSub->GetBinContent(i));
      BenUnfoldData->SetBinContent(i,dummyBenUnfoldData->GetBinContent(i));
      BenXsecMC->SetBinContent(i,dummyBenXsecMC->GetBinContent(i));
      BenEff->SetBinContent(i,dummyBenEff->GetBinContent(i));
      BenEffNum->SetBinContent(i,dummyBenEffNum->GetBinContent(i));
      AaronNewEff->SetBinContent(i,AaronEffNum->GetBinContent(i)/AaronEffDen->GetBinContent(i));
      BenXsecAaronSigData->SetBinContent(i, dummyBenXsecAaronSigData->GetBinContent(i));
      BenXsecNoTpiWegihtData->SetBinContent(i, dummyBenXsecNoTpiWegihtData->GetBinContent(i));
      BenAaronSigEffDen->SetBinContent(i, dummyBenAaronSigEffDen->GetBinContent(i));
      BenAaronSigEffNum->SetBinContent(i, dummyBenAaronSigEffNum->GetBinContent(i));
      BenAaronSigEff->SetBinContent(i, dummyBenAaronSigEff->GetBinContent(i));
      BenAaronSigEffCorr->SetBinContent(i, dummyBenAaronSigEffCorr->GetBinContent(i));
      BenAaronSigUnfold->SetBinContent(i, dummyBenAaronSigUnfold->GetBinContent(i));
      BenAaronSigDataSel->SetBinContent(i, dummyBenAaronSigDataSel->GetBinContent(i));
      BenAaronSigBGSub->SetBinContent(i, dummyBenAaronSigBGSub->GetBinContent(i));
  
      BenAaronSigCutsEffDen->SetBinContent(i, dummyBenAaronSigCutsEffDen->GetBinContent(i));
      BenAaronSigCutsEffNum->SetBinContent(i, dummyBenAaronSigCutsEffNum->GetBinContent(i));
      BenAaronSigCutsEff->SetBinContent(i, dummyBenAaronSigCutsEff->GetBinContent(i));
      BenAaronSigCutsEffCorr->SetBinContent(i, dummyBenAaronSigCutsEffCorr->GetBinContent(i));
      BenAaronSigCutsUnfold->SetBinContent(i, dummyBenAaronSigCutsUnfold->GetBinContent(i));
      BenAaronSigCutsDataSel->SetBinContent(i, dummyBenAaronSigCutsDataSel->GetBinContent(i));
      BenAaronSigCutsBGSub->SetBinContent(i, dummyBenAaronSigCutsBGSub->GetBinContent(i));
    } 
  
    double AaronDataPOT = h_data_POT->GetBinContent(1);
    double AaronMCPOT = h_MC_POT->GetBinContent(1);
    double DataPOTScale = dataPOTBenEve->GetBinContent(1)/AaronDataPOT; 
    double MCPOTScale = MCPOTBenEve->GetBinContent(1)/AaronMCPOT;
    double POTDataScaleASig = 
           dataPOTBenEve->GetBinContent(1)/dataPOTBenAaronSig->GetBinContent(1);
    double POTMCScaleASig =
           MCPOTBenEve->GetBinContent(1)/MCPOTBenAaronSig->GetBinContent(1);
    double POTDataScaleASigCuts =
           dataPOTBenEve->GetBinContent(1)/dataPOTBenAaronSigCuts->GetBinContent(1);
    double POTMCScaleASigCuts =
           MCPOTBenEve->GetBinContent(1)/MCPOTBenAaronSigCuts->GetBinContent(1);
    std::cout << "Data scale = " << DataPOTScale << " MC Scale = " << DataPOTScale << " AaronSigDatascale = " << POTDataScaleASig << " AaronSigMCScale = " << POTMCScaleASig << "\n";
  
    AaronEffDenOtherFile->Scale(AaronMCPOT/AaronDataPOT);
    AaronEffNumOtherFile->Scale(AaronMCPOT/AaronDataPOT);
    AaronEffDenOtherFile->Scale(MCPOTScale);
    AaronEffNumOtherFile->Scale(MCPOTScale);
    BenAaronSigEffDen->Scale(POTMCScaleASig);
    BenAaronSigEffNum->Scale(POTMCScaleASig);
    BenAaronSigCutsEffDen->Scale(POTMCScaleASigCuts);
    BenAaronSigCutsEffNum->Scale(POTMCScaleASigCuts);
    
  //  AaronEffOtherFile->Scale(AaronMCPOT/AaronDataPOT);
    AaronUnfoldData->Scale(DataPOTScale);
    AaronCorrectedData->Scale(DataPOTScale);
    AaronDataSel->Scale(DataPOTScale);
    AaronBGSub->Scale(DataPOTScale);
  
    BenAaronSigUnfold->Scale(POTDataScaleASig);
    BenAaronSigDataSel->Scale(POTDataScaleASig);
    BenAaronSigBGSub->Scale(POTDataScaleASig);
    BenAaronSigEffCorr->Scale(POTDataScaleASig);
  
    BenAaronSigCutsUnfold->Scale(POTDataScaleASigCuts);
    BenAaronSigCutsDataSel->Scale(POTDataScaleASigCuts);
    BenAaronSigCutsBGSub->Scale(POTDataScaleASigCuts);
    BenAaronSigCutsEffCorr->Scale(POTDataScaleASigCuts);
  
  
    BenEffDen->Scale(MCPOTScale);
    hvec.push_back(BenXsecData);
  //  hvec.push_back(BenXsecDataCorr);
//    hvec.push_back(BenXsecAaronSigData);
//    hvec.push_back(BenXsecNoTpiWegihtData);
    hvecEffNum.push_back(BenEffNum); 
//    hvecEffNum.push_back(BenAaronSigEffNum); 
//    hvecEffNum.push_back(BenAaronSigCutsEffNum); 
    hvecEffDen.push_back(BenEffDen); 
//    hvecEffDen.push_back(BenAaronSigEffDen); 
//    hvecEffDen.push_back(BenAaronSigCutsEffDen); 
   
  
    hvecEff.push_back(BenEff); 
//    hvecEff.push_back(BenAaronSigEff); 
//    hvecEff.push_back(BenAaronSigCutsEff); 
  
    hvecEffCorr.push_back(BenEffCorrData); 
//    hvecEffCorr.push_back(BenAaronSigEffCorr); 
  //  hvecEffCorr.push_back(BenAaronSigCutsEffCorr); 
   
    hvecUnfold.push_back(BenUnfoldData); 
//    hvecUnfold.push_back(BenAaronSigUnfold); 
  //  hvecUnfold.push_back(BenAaronSigCutsUnfold); 
  
    hvecDataSel.push_back(BenDataSel); 
//    hvecDataSel.push_back(BenAaronSigDataSel); 
  //  hvecDataSel.push_back(BenAaronSigCutsDataSel); 
  
    hvecBGSub.push_back(BenBGSub); 
//    hvecBGSub.push_back(BenAaronSigBGSub); 
  //  hvecBGSub.push_back(BenAaronSigCutsBGSub); 
  
    AaronRevUnfoldData->Multiply(AaronRevUnfoldData,AaronFluxNormalizer);
    AaronRevUnfoldData->Multiply(AaronRevUnfoldData,AaronNewEff);
    AaronRevUnfoldData->Scale(h_data_POT->GetBinContent(1));
  
    AaronNoNormXsecData->Multiply(AaronNoNormXsecData, AaronFluxNormalizer);
    AaronNoNormXsecData->Scale(h_data_POT->GetBinContent(1));

    std::string xlabel = "";
    if (var == "q2"){
      xlabel = "Q^{2} (GeV^{2})";
    }
    if (var == "mixtpi"){
      xlabel = "T_{#pi} (GeV)";
    }
    if (var == "thetamu_deg"){
      xlabel = "#theta_{#mu} (deg)";
    }
    PlotRatio(BenXsecMC, AaronXsecMC, var, 1.,"BenSigvsAaronXSecMC", false, 
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenXsecData, BenXsecAaronSigData, var, 1.,"BenMacroAaronSigXSecData", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenEffDen, AaronEffDen, var, MCPOTScale,"BenSigvsAaronEffDen", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenEffCorrData, AaronCorrectedData, var, DataPOTScale,"BenSigvsAaronEffCorrection", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenEffCorrData, AaronNoNormXsecData, var, DataPOTScale,"BenSigvsAaronNoNormEffCorr", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenEff, AaronNewEff, var, 1.,"BenSigvsAaronNewEff", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(AaronNewEff, AaronEff, var, 1.,"AaronEffCheck", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenUnfoldData, AaronUnfoldData, var, DataPOTScale,"BenSigvsAaronUnfoldData", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenUnfoldData, AaronRevUnfoldData, var, DataPOTScale,"BenSigvsAaronRevUnfoldData", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(AaronRevUnfoldData, AaronUnfoldData, var, 1.,"AaronUnfoldDataCheck", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenEff, AaronEff, var, 1.,"BenSigvsAaronEff", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(AaronEffOtherFile, AaronEff, var, 1.,"AaronEffDiffFiles", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(AaronNewEff, AaronEffOtherFile, var, 1.,"AaronNewEffvsOtherFile", false,       
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(BenEffNum, AaronEffNum, var, MCPOTScale,"BenSigvsAaronEffNum", false,
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(AaronEffNumOtherFile, AaronEffNum, var, 1.,"AaronEffNumDiffFiles", false,  
             "BenSig/Aaronsig", xlabel); 
    PlotRatio(AaronEffDenOtherFile, AaronEffDen, var, 1.,"AaronEffDenDiffFiles", false,  
             "BenSig/Aaronsig", xlabel); 
   // PlotRatio(AaronEffNum, AaronEffDen, var, 1.,"AaronEfficiency", true); 
    PlotRatioVec(hvec, AaronXsecData, var, 1.,"CrossSectionsData", true,  
             "BenSig/Aaronsig", xlabel); 
    PlotRatioVec(hvecEffNum, AaronEffNumOtherFile, var, 1.,"EffNum",false, 
             "BenSig/Aaronsig", xlabel); 
    PlotRatioVec(hvecEffDen, AaronEffDenOtherFile, var, 1.,"EffDen", true, 
             "BenSig/Aaronsig", xlabel); 
    PlotRatioVec(hvecEff, AaronEffOtherFile, var, 1.,"Eff",false, 
             "BenSig/Aaronsig", xlabel); 
    PlotRatioVec(hvecEffCorr, AaronCorrectedData, var, 1.,"EffCorr",false, 
             "BenSig/Aaronsig", xlabel); 
    PlotRatioVec(hvecUnfold, AaronUnfoldData, var, 1.,"Unfold",false, 
             "BenSig/Aaronsig", xlabel); 
    PlotRatioVec(hvecDataSel, AaronDataSel, var, 1.,"DataSelection",false, 
             "BenSig/Aaronsig", xlabel); 
    PlotRatioVec(hvecBGSub, AaronBGSub, var, 1.,"BGSub",false, 
             "BenSig/Aaronsig", xlabel);  
    std::cout << "PASS 10\n";
  }
    fBenEve.Close();
    fAaronGENIE.Close();
    fAaronNoLowTpiGENIE.Close();
    fBenEveGENIE.Close();

    fAaronXsec.Close();
    fAaronXsecIn.Close();
    fAaronEff.Close();
    fBenEveMAaronSD.Close();
    fBenEveSDNoTpi.Close();
    fin8.Close();

}
#endif  

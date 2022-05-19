#ifndef crossSectionDataFromFile_C
#define crossSectionDataFromFile_C

#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvVertErrorBand.h"
#include "PlotUtils/TargetUtils.h"
#include "TDirectory.h"
#include "TFile.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/Cuts.h"
#include "includes/MacroUtil.h"
#include "includes/SignalDefinition.h"
#include "includes/Systematics.h"                // GetSystematicUniversesMap
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/Variable.h"
#include "includes/WSidebandFitter.h"
#include "includes/common_functions.h"  // GetVar, CopyHists, WritePOT, erase_if, uniq
#include "makeCrossSectionMCInputs.C"   // GetAnalysisVariables
#include "plotting_functions.h"

void LoopAndFillData(const CCPi::MacroUtil& util,
                     std::vector<Variable*> variables) {
  // Fill data distributions.
  const bool is_mc = false;
  const bool is_truth = false;
  std::cout << "*** Starting Data Loop ***" << std::endl;
  for (Long64_t i_event = 0; i_event < util.GetDataEntries(); ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    util.m_data_universe->SetEntry(i_event);

    CCPiEvent event(is_mc, is_truth, util.m_signal_definition,
                    util.m_data_universe);

    event.m_passes_cuts =
        PassesCuts(*util.m_data_universe, event.m_reco_pion_candidate_idxs,
                   is_mc, util.m_signal_definition, event.m_is_w_sideband);

    event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
    // event.m_is_w_sideband = IsWSideband(event);

    ccpi_event::FillRecoEvent(event, variables);
  }
  std::cout << "*** Done Data ***\n\n";
}

// Do W Sideband fit in every universe, return weight hist wrappers by
// reference.
void DoWSidebandTune(CCPi::MacroUtil& util, Variable* fit_var, CVHW& loW_wgt,
                     CVHW& midW_wgt, CVHW& hiW_wgt) {
  // DO FIT
  // Fit mc to data in every universe
  std::cout << "Fitting in the variable " << sidebands::kFitVarString << "\n";
  for (auto error_band : util.m_error_bands) {
    std::vector<CVUniverse*> universes = error_band.second;
    for (auto universe : universes) {
      WSidebandFitter wsb_fitter =
          WSidebandFitter(*universe, fit_var->m_hists, util.m_pot_scale);
      wsb_fitter.Fit();
      std::cout << "low weight:" << (wsb_fitter.m_fit_scale)[kLoWParamId] << "\n";
      std::cout << "mid weight:" << (wsb_fitter.m_fit_scale)[kMidWParamId] << "\n";
      std::cout << "hi weight:" << (wsb_fitter.m_fit_scale)[kHiWParamId] << "\n";
      // Store the outputs of the fits in HistWrappers
      loW_wgt.univHist(universe)->SetBinContent(
          1, (wsb_fitter.m_fit_scale)[kLoWParamId]);
      midW_wgt.univHist(universe)->SetBinContent(
          1, (wsb_fitter.m_fit_scale)[kMidWParamId]);
      hiW_wgt.univHist(universe)->SetBinContent(
          1, (wsb_fitter.m_fit_scale)[kHiWParamId]);
    }
  }

  // SYNCH AND SAVE FIT RESULTS
  std::cout << "Synching and writing bg fit stuff\n";
  loW_wgt.SyncCVHistos();
  midW_wgt.SyncCVHistos();
  hiW_wgt.SyncCVHistos();

  loW_wgt.hist->Write("fit_param_loW");
  midW_wgt.hist->Write("fit_param_midW");
  hiW_wgt.hist->Write("fit_param_hiW");
}

// The sideband fit parameters that come out of the fit are in the form of a
// histwrapper with one bin.
// In order to perform operations (multiplication) with the weights we have to
// bin the HW's using the variable's binning.
void RebinFitParamHists(UniverseMap error_bands, const double nbins,
                        const CVHW& loW_wgt, const CVHW& midW_wgt,
                        const CVHW& hiW_wgt, CVHW& loW_wgt_rebin,
                        CVHW& midW_wgt_rebin, CVHW& hiW_wgt_rebin) {
  for (auto error_band : error_bands) {
    std::vector<CVUniverse*> universes = error_band.second;
    for (auto universe : universes) {
      // Get scales in this universe
      const double loW_scale = loW_wgt.univHist(universe)->GetBinContent(1);
      const double midW_scale = midW_wgt.univHist(universe)->GetBinContent(1);
      const double hiW_scale = hiW_wgt.univHist(universe)->GetBinContent(1);

      // loop bins
      for (int i = 1; i <= nbins; ++i) {
        // std::cout << loW_wgt_rebin.hist->GetXaxis()->GetBinLowEdge(i) <<
        // std::endl;
        loW_wgt_rebin.univHist(universe)->SetBinContent(i, loW_scale);
        midW_wgt_rebin.univHist(universe)->SetBinContent(i, midW_scale);
        hiW_wgt_rebin.univHist(universe)->SetBinContent(i, hiW_scale);
      }
    }  // universes
  }    // error bands

  // Be sure to drink your ovaltine!!
  loW_wgt_rebin.SyncCVHistos();
  midW_wgt_rebin.SyncCVHistos();
  hiW_wgt_rebin.SyncCVHistos();
}

// Apply the sideband tunes to the untuned BG
// return tuned BG
void ScaleBG(Variable* var, CCPi::MacroUtil& util, const CVHW& loW_wgt,
             const CVHW& midW_wgt, const CVHW& hiW_wgt) {
  // REBIN FIT PARAMS FOR THIS VARIABLE
  // temp HW's with same binning and error bands as variable
  PlotUtils::HistWrapper<CVUniverse> loW_wgt_rebin, midW_wgt_rebin,
      hiW_wgt_rebin;

  InitializeHW(
      var, Form("loW_fit_wgt_%s", var->Name().c_str()),
      Form("W Sideband Fit Weight -- low W -- %s", var->Name().c_str()),
      util.m_error_bands, loW_wgt_rebin);
  InitializeHW(
      var, Form("midW_fit_wgt_%s", var->Name().c_str()),
      Form("W Sideband Fit Weight -- mid W -- %s", var->Name().c_str()),
      util.m_error_bands, midW_wgt_rebin);
  InitializeHW(var, Form("hiW_fit_wgt_%s", var->Name().c_str()),
               Form("W Sideband Fit Weight -- hi W -- %s", var->Name().c_str()),
               util.m_error_bands, hiW_wgt_rebin);

  // For these MnvH1Ds ^ set each bin of each universe with the universe's fit
  // weights
  RebinFitParamHists(util.m_error_bands, var->NBins(), loW_wgt, midW_wgt,
                     hiW_wgt, loW_wgt_rebin, midW_wgt_rebin, hiW_wgt_rebin);

  // APPLY TUNE TO EACH BG COMPONENT
  // tuned bg component = clone (untuned component)
  PlotUtils::MnvH1D* tuned_bg_loW =
      (PlotUtils::MnvH1D*)var->m_hists.m_bg_loW.hist->Clone(uniq());
  PlotUtils::MnvH1D* tuned_bg_midW =
      (PlotUtils::MnvH1D*)var->m_hists.m_bg_midW.hist->Clone(uniq());
  PlotUtils::MnvH1D* tuned_bg_hiW =
      (PlotUtils::MnvH1D*)var->m_hists.m_bg_hiW.hist->Clone(uniq());

  // tuned bg component = untuned component * component wgt
  tuned_bg_loW->Multiply(tuned_bg_loW, loW_wgt_rebin.hist);
  tuned_bg_midW->Multiply(tuned_bg_midW, midW_wgt_rebin.hist);
  tuned_bg_hiW->Multiply(tuned_bg_hiW, hiW_wgt_rebin.hist);

  // SUM TUNED COMPONENTS
  // total tuned bg = sum of tuned components
  PlotUtils::MnvH1D* tuned_bg = (PlotUtils::MnvH1D*)tuned_bg_loW->Clone(uniq());
  tuned_bg->Add(tuned_bg_midW);
  tuned_bg->Add(tuned_bg_hiW);

  // WRITE TUNED BG
  tuned_bg_loW->Write(Form("tuned_bg_loW_%s", var->Name().c_str()));
  tuned_bg_midW->Write(Form("tuned_bg_midW_%s", var->Name().c_str()));
  tuned_bg_hiW->Write(Form("tuned_bg_hiW_%s", var->Name().c_str()));
  tuned_bg->Write(Form("tuned_bg_%s", var->Name().c_str()));

  //// SCALE TUNED BG TO DATA
  //  tuned_bg_loW ->Scale(util.m_pot_scale);
  //  tuned_bg_midW->Scale(util.m_pot_scale);
  //  tuned_bg_hiW ->Scale(util.m_pot_scale);
  //  tuned_bg     ->Scale(util.m_pot_scale);

  //// WRITE TUNED & SCALED BG
  //  tuned_bg_loW ->Write(Form("tuned_bg_loW_POTscaled_%s",
  //  var->Name().c_str()));
  //  tuned_bg_midW->Write(Form("tuned_bg_midW_POTscaled_%s",
  //  var->Name().c_str())); tuned_bg_hiW
  //  ->Write(Form("tuned_bg_hiW_POTscaled_%s",  var->Name().c_str())); tuned_bg
  //  ->Write(Form("tuned_bg_POTscaled_%s", var->Name().c_str()));

  var->m_hists.m_tuned_bg = tuned_bg;
}

//==============================================================================
// Main
//==============================================================================
void crossSectionDataFromFile(int signal_definition_int = 0,
                              const char* plist = "ME1A") {
  //============================================================================
  // Setup
  //============================================================================

  // I/O
  TFile fin("MCXSecInputs_20220519_NoSys.root", "READ");
  std::cout << "Reading input from " << fin.GetName() << endl;

  TFile fout("DataXSecInputs_20220519_NoSys.root", "RECREATE");
  std::cout << "Output file is " << fout.GetName() << "\n";

  std::cout << "Copying all hists from fin to fout\n";
  CopyHists(fin, fout);

  // INPUT TUPLES
  // Don't actually use the MC chain, only load it to indirectly access its
  // systematics
  std::string data_file_list = GetPlaylistFile(plist, false);
  std::string mc_file_list = GetPlaylistFile("ME1A", true);

  // Macro Utility
  const std::string macro("CrossSectionDataFromFile");
  bool do_truth = false, is_grid = false, do_systematics = true;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // POT
  SetPOT(fin, fout, util);

  // Variables and histograms -- load in MC hists from fin
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  {  // remove unwanted variables
    ContainerEraser::erase_if(
        variables, [](Variable* v) { return v->Name() == "tpi_mbr"; });
    /*
      ContainerEraser::erase_if(variables, [](Variable* v) {
          return v->Name() == "tpi" || v->Name() == "enu"; });
      ContainerEraser::erase_if(variables, [](Variable* v) {
          return v->Name() == "thetapi_deg" || v->Name() == "thetamu_deg"; });
      ContainerEraser::erase_if(variables, [](Variable* v) {
          return v->Name() == "q2" || v->Name() == "wexp"; });
      ContainerEraser::erase_if(variables, [](Variable* v) {
          return v->Name() == "ptmu" || v->Name() == "pzmu"; });
    */
  }

  for (auto v : variables) {
    std::cout << "Loading hists for variable " << v->Name() << "\n";
    v->LoadMCHistsFromFile(fin, util.m_error_bands);
    v->InitializeDataHists();
  }

  //============================================================================
  // Loop Data and Make Event Selection
  //============================================================================

  LoopAndFillData(util, variables);

  // Add empty error bands to data hists and fill their CVs
  for (auto v : variables) {
    v->m_hists.m_selection_data->ClearAllErrorBands();
    v->m_hists.m_selection_data->AddMissingErrorBandsAndFillWithCV(
        *v->m_hists.m_selection_mc.hist);
  }

  SaveDataHistsToFile(fout, variables);

  //============================================================================
  // Tune Sideband
  //============================================================================
  // Hists to store the sideband tune/fit parameters/weights
  PlotUtils::HistWrapper<CVUniverse> hw_loW_fit_wgt =
      PlotUtils::HistWrapper<CVUniverse>("h_loW_fit_wgt",
                                         "W Sideband Fit Weight -- low W", 1,
                                         0., 15., util.m_error_bands);

  PlotUtils::HistWrapper<CVUniverse> hw_midW_fit_wgt =
      PlotUtils::HistWrapper<CVUniverse>("h_midW_fit_wgt",
                                         "W Sideband Fit Weight -- mid W", 1,
                                         0., 15., util.m_error_bands);

  PlotUtils::HistWrapper<CVUniverse> hw_hiW_fit_wgt =
      PlotUtils::HistWrapper<CVUniverse>("h_hiW_fit_wgt",
                                         "W Sideband Fit Weight -- high W", 1,
                                         0., 15., util.m_error_bands);

  // Sideband tune
  // Fill the fit parameter hists (by reference)
  DoWSidebandTune(util, GetVar(variables, sidebands::kFitVarString),
                  hw_loW_fit_wgt, hw_midW_fit_wgt, hw_hiW_fit_wgt);

  // Write fit weights
  hw_loW_fit_wgt.hist->Write("loW_fit_wgt");
  hw_midW_fit_wgt.hist->Write("midW_fit_wgt");
  hw_hiW_fit_wgt.hist->Write("hiW_fit_wgt");

  //============================================================================
  // In a loop over variables...
  // 1. Scale Background
  // 2. Subtract Background
  // 3. Unfold
  // 4. Efficiency Correct
  //============================================================================
  for (auto var : variables) {
    // skip non-analysis variables
    if (var->m_is_true) continue;
    if (var->Name() == std::string("tpi_mbr")) continue;
    if (var->Name() == sidebands::kFitVarString) continue;

    const char* name = var->Name().c_str();
    std::cout << "Calculating Cross Section for " << name << "\n";

    // We'll need the true version of this variable later on. Get it now.
    Variable* true_var = GetVar(variables, var->Name() + std::string("_true"));

    //============================================================================
    // Scale BG
    // i.e. apply W sideband fit to the BG in the signal region.
    //============================================================================
    ScaleBG(var, util, hw_loW_fit_wgt, hw_midW_fit_wgt, hw_hiW_fit_wgt);

    // POT scale the tuned BG
    PlotUtils::MnvH1D* tuned_POTscaled_bg =
        (PlotUtils::MnvH1D*)var->m_hists.m_tuned_bg->Clone(uniq());
    tuned_POTscaled_bg->Scale(util.m_pot_scale);

    std::cout << "  Done BG Tune\n";

    //============================================================================
    // Subtract BG
    //============================================================================
    // (Make sure empty error bands have been added to data hists)
    var->m_hists.m_bg_subbed_data =
        (PlotUtils::MnvH1D*)var->m_hists.m_selection_data->Clone(uniq());
    var->m_hists.m_bg_subbed_data->Add(tuned_POTscaled_bg, -1);

    // Write BG Sub Data
    fout.cd();
    var->m_hists.m_bg_subbed_data->Write(Form("bg_subbed_data_%s", name));

    std::cout << "  Done BG Sub\n";

    //============================================================================
    // Unfold
    //============================================================================
    MinervaUnfold::MnvUnfold mnv_unfold;
    mnv_unfold.setUseBetterStatErrorCalc(true);
    PlotUtils::MnvH2D* migration =
        (PlotUtils::MnvH2D*)var->m_hists.m_migration.hist->Clone(uniq());
    PlotUtils::MnvH1D* bg_sub_data =
        (PlotUtils::MnvH1D*)var->m_hists.m_bg_subbed_data->Clone(uniq());
    int n_iterations = 4;
    if (var->Name() == "tpi" || var->Name() == "wexp" ||
        var->Name() == "thetapi")
      n_iterations = 10;

    mnv_unfold.UnfoldHisto(var->m_hists.m_unfolded, migration, bg_sub_data,
                           RooUnfold::kBayes, n_iterations);

    // copypasta
    // Blurgh. We want the covariance matrix induced by the unfolding,
    // but MnvUnfold will only give that back to us with a call to a
    // different version of UnfoldHisto that only takes a TH1D, and
    // not a MnvH1D (so we can't just combine it with the previous call)
    TMatrixD unfolding_cov_matrix_orig;
    TH1D* unfolded_dummy =
        new TH1D(var->m_hists.m_unfolded->GetCVHistoWithStatError());
    TH2D* migration_dummy = new TH2D(migration->GetCVHistoWithStatError());
    TH1D* reco_dummy =
        new TH1D(migration->ProjectionX()->GetCVHistoWithStatError());
    TH1D* truth_dummy =
        new TH1D(migration->ProjectionY()->GetCVHistoWithStatError());
    TH1D* bg_sub_data_dummy = new TH1D(bg_sub_data->GetCVHistoWithStatError());
    mnv_unfold.UnfoldHisto(unfolded_dummy, unfolding_cov_matrix_orig,
                           migration_dummy, reco_dummy, truth_dummy,
                           bg_sub_data_dummy, RooUnfold::kBayes, 4);

    // Add cov matrix to unfolded hist
    var->m_hists.m_unfolded->PushCovMatrix(
        Form("unfolding_cov_matrix_%s", name), unfolding_cov_matrix_orig);

    // Write unfolded
    fout.cd();
    var->m_hists.m_unfolded->Write(Form("unfolded_%s", name));

    std::cout << "  Done Unfolding\n";

    //============================================================================
    // Efficiency Correct
    //============================================================================
    // Calculate efficiency

    /*
      Delete me
      { // Somehow effnum and effden have 200 flux universes
        MnvVertErrorBand *poppedFluxErrorBand =
      true_var->m_hists.m_effnum.hist->PopVertErrorBand("Flux");
        std::vector<TH1D*> fluxUniverses = poppedFluxErrorBand->GetHists();
        fluxUniverses.resize(100);
        true_var->m_hists.m_effnum.hist->AddVertErrorBand("Flux",fluxUniverses);
      }
      { // Somehow effnum and effden have 200 flux universes
        MnvVertErrorBand *poppedFluxErrorBand =
      true_var->m_hists.m_effden.hist->PopVertErrorBand("Flux");
        std::vector<TH1D*> fluxUniverses = poppedFluxErrorBand->GetHists();
        fluxUniverses.resize(100);
        true_var->m_hists.m_effden.hist->AddVertErrorBand("Flux",fluxUniverses);
      }
    */

    var->m_hists.m_efficiency =
        (PlotUtils::MnvH1D*)true_var->m_hists.m_effnum.hist->Clone(uniq());
    var->m_hists.m_efficiency->Divide(true_var->m_hists.m_effnum.hist,
                                      true_var->m_hists.m_effden.hist);

    // if(var->Name() == "ptmu")
    //  PrintUniverseContent(true_var->m_hists.m_effnum.hist);
    // if(var->Name() == "ptmu")
    //  PrintUniverseContent(true_var->m_hists.m_effden.hist);

    PlotUtils::MnvH1D* h_efficiency_corrected_data =
        (PlotUtils::MnvH1D*)var->m_hists.m_unfolded->Clone(uniq());
    // h_efficiency_corrected_data->ClearSysErrorMatrices(); // maybe we'll
    // write a new matrix when we divide? NOPE doesn't work.

    // Efficiency correct
    h_efficiency_corrected_data->Divide(var->m_hists.m_unfolded,
                                        var->m_hists.m_efficiency);

    TMatrixD unfolding_cov_matrix_effcor =
        h_efficiency_corrected_data->GetSysErrorMatrix(
            Form("unfolding_cov_matrix_%s", name));

    {  // Check to make sure the covariance matrix got divided correctly
       // for (int i = 0; i < unfolding_cov_matrix_effcor.GetNcols(); ++i) {
       //  for (int j = 0; j < unfolding_cov_matrix_effcor.GetNrows(); ++j) {
       //    std::cout << unfolding_cov_matrix_effcor[j][i] -
       //    unfolding_cov_matrix_orig[j][i];
       //  }
       //  std::cout  << "\n";
       //}
    }

    // Write efficiency and efficiency-corrected data
    fout.cd();
    var->m_hists.m_efficiency->Write(Form("efficiency_%s", name));
    h_efficiency_corrected_data->Write(
        Form("efficiency_corrected_data_%s", name));

    std::cout << "  Done Efficiency Correcting\n";

    //============================================================================
    // Normalization -- integrated flux, targets, POT (and don't forget MC,
    // too!)
    //============================================================================
    // Init the normalization hist from the eff corr just to get the error bands
    PlotUtils::MnvH1D* h_flux_normalization =
        (PlotUtils::MnvH1D*)h_efficiency_corrected_data->Clone(
            "flux_normalization");
    h_flux_normalization->ClearAllErrorBands();
    h_flux_normalization->Reset();

    // Get the flux histo, to be integrated
    static PlotUtils::FluxReweighter* frw = new PlotUtils::FluxReweighter(
        14, CCNuPionIncConsts::kUseNueConstraint, "minervame1D1M1NWeightedAve",
        PlotUtils::FluxReweighter::gen2thin,
        PlotUtils::FluxReweighter::g4numiv6,
        CCNuPionIncConsts::kNFluxUniverses);

    fout.cd();  // FRW opens a new file and changes our current dir.

    h_flux_normalization = frw->GetIntegratedFluxReweighted(
        14, h_efficiency_corrected_data, 0., 100.);

    //{ // Truncate flux universes to 10!!
    //  MnvVertErrorBand *poppedFluxErrorBand =
    //  h_flux_normalization->PopVertErrorBand("Flux"); std::vector<TH1D*>
    //  fluxUniverses = poppedFluxErrorBand->GetHists();
    //  fluxUniverses.resize(10);
    //  h_flux_normalization->AddVertErrorBand("Flux",fluxUniverses);
    //}

    {  //// Truncate flux universes to 10!!
       // MnvVertErrorBand *poppedFluxErrorBand =
       // h_flux_normalization->PopVertErrorBand("Flux"); std::vector<TH1D*>
       // fluxUniverses = poppedFluxErrorBand->GetHists(); std::cout << "flux
       // universes: " << fluxUniverses.size() << "\n";
       ////fluxUniverses.resize(10);
       ////h_flux_normalization->AddVertErrorBand("Flux",fluxUniverses);
    }

    //// remove redundant error bands
    // h_flux_normalization->PopVertErrorBand("Flux_BeamFocus");
    // h_flux_normalization->PopVertErrorBand("ppfx1_Total");

    // Convert flux units from nu/m^2/POT to nu/cm^2/POT
    h_flux_normalization->Scale(1.0e-4);

    // Divide flux integral
    PlotUtils::MnvH1D* h_cross_section =
        (PlotUtils::MnvH1D*)h_efficiency_corrected_data->Clone(uniq());

    h_cross_section->AddMissingErrorBandsAndFillWithCV(*h_flux_normalization);

    h_cross_section->Divide(h_cross_section, h_flux_normalization);

    // targets and POT norm
    static const double apothem = 865.;
    static const double upstream = 5900.;    // ~module 25 plane 1
    static const double downstream = 8430.;  // ~module 81 plane 1

    double n_target_nucleons =
        PlotUtils::TargetUtils::Get().GetTrackerNNucleons(upstream, downstream,
                                                          false,  // isMC
                                                          apothem);

    std::cout << "  flux_integral cv = "
              << h_flux_normalization->GetBinContent(1) << "\n";
    std::cout << "  N target nucleons = " << n_target_nucleons << "\n";
    std::cout << "  data pot = " << util.m_data_pot << "\n";
    static const double data_scale =
        1.0 / (n_target_nucleons * util.m_data_pot);
    h_cross_section->Scale(data_scale);

    // Write data cross section
    fout.cd();
    h_cross_section->Write(Form("cross_section_%s", name));

    // Begin MC normalization
    PlotUtils::MnvH1D* h_mc_cross_section =
        (PlotUtils::MnvH1D*)true_var->m_hists.m_effden.hist->Clone(uniq());

    h_mc_cross_section->AddMissingErrorBandsAndFillWithCV(
        *h_flux_normalization);
    h_mc_cross_section->Divide(h_mc_cross_section, h_flux_normalization);

    std::cout << "  mc pot = " << util.m_mc_pot << "\n";
    static const double mc_scale = 1.0 / (n_target_nucleons * util.m_mc_pot);
    h_mc_cross_section->Scale(mc_scale);

    // Write mc cross section
    fout.cd();
    h_mc_cross_section->Write(Form("mc_cross_section_%s", name));

    // Set covariance matrix diagonal to zero
    // copypasta
    // Jeremy tells me that the covariance matrix has the diagonal
    // errors on it, which are already included elsewhere, so we have to
    // subtract them off before adding the unfolding covariance matrix
    // back on
    TMatrixD unfolding_cov_matrix = h_cross_section->GetSysErrorMatrix(
        Form("unfolding_cov_matrix_%s", name));
    for (int i = 0; i < unfolding_cov_matrix.GetNrows(); ++i)
      unfolding_cov_matrix(i, i) = 0;
    h_cross_section->PushCovMatrix(Form("unfolding_cov_matrix_%s", name),
                                   unfolding_cov_matrix);

    // Write scaled covariance matrix
    fout.cd();
    unfolding_cov_matrix.Write(Form("unfolding_cov_matrix_%s", name));

    std::cout << "  Done flux, targets, and POT normalization\n";
  }  // vars loop
}

// PlotUtils::MnvH1D* efficiency_numerator   =
// (PlotUtils::MnvH1D*)true_var->m_hists.m_effnum.hist->Clone(uniq());
// PlotUtils::MnvH1D* efficiency_denominator =
// (PlotUtils::MnvH1D*)true_var->m_hists.m_effden.hist->Clone(uniq());
// var->m_hists.m_efficiency = (PlotUtils::MnvH1D*)efficiency->Clone(uniq());

#endif  // crossSectionDataFromFile_C

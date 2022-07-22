#ifndef GXSEClosure_C
#define GXSEClosure_C

#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/TargetUtils.h"
#include "TFile.h"
#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/MacroUtil.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "plotting_functions.h"

//==============================================================================
// Main
// TODO function-ify the real xsec script and use those functions.
// Instead of copy-paste.
//==============================================================================
void GXSEClosure(int signal_definition_int = 0) {
  // In and outfiles
  // TFile fin("rootfiles/MCXSecInputs_20190903.root", "READ");
  TFile fin(
      "/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/"
      "MCXSecInputs_20220717_ME1A_NoSys.root",
      "READ");
  cout << "Reading input from " << fin.GetName() << endl;

  // Set up macro utility object -- which does the systematics for us
  const char* plist = "ME1A";
  std::string data_file_list = GetPlaylistFile(plist, false);
  std::string mc_file_list = GetPlaylistFile(plist, true);
  // std::string mc_file_list =
  //    "/minerva/app/users/granados/cmtuser/MINERvA101/"
  //    "MINERvA-101-Cross-Section/MCME1A.txt";
  bool do_truth = false;
  bool do_systematics = true, do_grid = false;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, do_grid, do_systematics);

  // Set POT
  PlotUtils::MnvH1D* h_mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float mc_pot = h_mc_pot->GetBinContent(1);
  // checking that pot of the input file and pot of the "official file list" are
  // the same.
  std::cout << util.m_mc_pot / mc_pot << "\n";
  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  // for (auto var : variables) var->LoadMCHistsFromFile(fin,
  // util.m_error_bands);

  for (auto var : variables) {
    if (var->Name() == std::string("tpi_mbr")) continue;
    if (var->Name() == sidebands::kFitVarString) continue;
    if (var->Name() != "pmu" || var->Name() != "pmu_true") continue;

    var->LoadMCHistsFromFile(fin, util.m_error_bands);

    if (var->m_is_true) continue;

    const char* name = var->Name().c_str();

    // We'll be needing the true version of this variable
    Variable* reco_var = GetVar(variables, var->Name());
    Variable* true_var = GetVar(variables, var->Name() + std::string("_true"));

    // Closure at the background subtraction. Step 1
    PlotUtils::MnvH1D* reco_sel_mc =
        (PlotUtils::MnvH1D*)reco_var->m_hists.m_selection_mc.hist->Clone(
            uniq());
    PlotUtils::MnvH1D* BG_untuned =
        (PlotUtils::MnvH1D*)reco_var->m_hists.m_bg.hist->Clone(uniq());
    PlotUtils::MnvH1D* signal_events_reco =
        (PlotUtils::MnvH1D*)reco_var->m_hists.m_effnum.hist->Clone(uniq());
    reco_sel_mc->Add(BG_untuned, -1);  // subtraction of the background

    PlotTogether(reco_sel_mc, "sel_mc", signal_events_reco, "signal",
                 "BG_closure");
    PlotRatio(reco_sel_mc, signal_events_reco, Form("BGClosure_%s", name), 1.,
              "", true);

    // Closure at the unfolding. Step 2
    MinervaUnfold::MnvUnfold mnv_unfold;
    PlotUtils::MnvH2D* migration =
        (PlotUtils::MnvH2D*)reco_var->m_hists.m_migration.hist->Clone(uniq());

    int n_iterations = 4;
    // if (var->Name() == "tpi" || var->Name() == "wexp" ||
    //  var->Name() == "thetapi")
    //  n_iterations = 10;

    mnv_unfold.UnfoldHisto(reco_var->m_hists.m_unfolded, migration, reco_sel_mc,
                           RooUnfold::kBayes, n_iterations);

    PlotUtils::MnvH1D* unfold =
        (PlotUtils::MnvH1D*)reco_var->m_hists.m_unfolded->Clone(uniq());
    PlotUtils::MnvH1D* true_effnum =
        (PlotUtils::MnvH1D*)true_var->m_hists.m_effnum.hist->Clone(uniq());
    PlotTogether(unfold, "unfold", true_effnum, "true_effnum",
                 "Unfolding_closure");
    PlotRatio(unfold, true_effnum, Form("UnfoldClosure_%s", name), 1., "",
              false);

    // Closure at the efficiency correction. Step 3
    reco_var->m_hists.m_efficiency =
        (PlotUtils::MnvH1D*)true_var->m_hists.m_effnum.hist->Clone(uniq());
    reco_var->m_hists.m_efficiency->Divide(true_var->m_hists.m_effnum.hist,
                                           true_var->m_hists.m_effden.hist);
    PlotUtils::MnvH1D* eff_corr = unfold;
    eff_corr->Divide(reco_var->m_hists.m_unfolded,
                     reco_var->m_hists.m_efficiency);
    PlotUtils::MnvH1D* true_effden =
        (PlotUtils::MnvH1D*)true_var->m_hists.m_effden.hist->Clone(uniq());
    PlotTogether(eff_corr, "eff_corr", true_effden, "true_effden",
                 "EffCorr_closure");
    PlotRatio(eff_corr, true_effden, Form("EffCorrClosure_%s", name), 1., "",
              false);

    // Start with "fake efficiency corrected data"
    PlotUtils::MnvH1D* h_mc_cross_section =
        (PlotUtils::MnvH1D*)true_var->m_hists.m_effden.hist->Clone(uniq());

    //========================================================================
    // Prepare normalizations -- flux, targets
    //========================================================================
    // Flux normalization
    PlotUtils::MnvH1D* h_flux_normalization =
        (PlotUtils::MnvH1D*)h_mc_cross_section->Clone("flux_normalization");
    h_flux_normalization->ClearAllErrorBands();
    h_flux_normalization->Reset();

    // Get the flux histo, to be integrated
    const bool use_hundred_universes = true;
    static PlotUtils::FluxReweighter* frw = new PlotUtils::FluxReweighter(
        14, true, PlotUtils::FluxReweighter::minervame1A,
        PlotUtils::FluxReweighter::gen2thin,
        PlotUtils::FluxReweighter::g4numiv6, 0);

    h_flux_normalization =
        frw->GetIntegratedFluxReweighted(14, h_mc_cross_section, 0., 100.);

    // remove redundant error bands from flux integral
    h_flux_normalization->PopVertErrorBand("Flux_BeamFocus");
    h_flux_normalization->PopVertErrorBand("ppfx1_Total");

    // Convert flux units from nu/m^2/POT to nu/cm^2/POT
    // Or could we already by in cm^-2?
    h_flux_normalization->Scale(1.0e-4);

    PlotUtils::MnvH1D* mc_integrate_flux =
        (PlotUtils::MnvH1D*)h_flux_normalization->Clone(uniq());

    // targets and POT norm
    static const double apothem = 850.;
    static const double upstream = 5990.;    // ~module 25 plane 1
    static const double downstream = 8340.;  // ~module 81 plane 1

    const bool is_mc = false;
    double n_target_nucleons =
        PlotUtils::TargetUtils::Get().GetTrackerNNucleons(upstream, downstream,
                                                          is_mc, apothem);

    //  double n_target_nucleons =
    //      PlotUtils::TargetUtils::Get().GetTrackerNNucleons(upstream,
    //      downstream,
    //                                                        /*isMC =*/false,
    //                                                        apothem);

    // Summarize scales
    std::cout << "  flux_integral cv = "
              << h_flux_normalization->GetBinContent(1) << "\n";
    std::cout << "  N target nucleons = " << n_target_nucleons << "\n";
    std::cout << "  mc pot = " << util.m_mc_pot << "\n";
    std::cout << "  flux_integral * N target nucleons = "
              << n_target_nucleons * h_flux_normalization->GetBinContent(1)
              << "\n";
    std::cout << "  flux_integral * N target nucleons * mc pot= "
              << n_target_nucleons * h_flux_normalization->GetBinContent(1) *
                     util.m_mc_pot
              << "\n";
    // static const double data_scale = 1.0 / 3.529606e+42;

    //========================================================================
    // do the normalizations
    //========================================================================
    // Flux
    h_mc_cross_section->AddMissingErrorBandsAndFillWithCV(
        *h_flux_normalization);
    h_mc_cross_section->Divide(h_mc_cross_section, h_flux_normalization);

    // Targets & Bin width normalization
    double normFactor = 13.569879;
    static const double mc_scale =
        1.0 / (n_target_nucleons * util.m_mc_pot * normFactor);
    // static const double mc_scale = 1.0 / 3.529606e+42;
    h_mc_cross_section->Scale(mc_scale, "width");
    true_effden->Scale(1. / normFactor, "width");
    mc_integrate_flux->Scale(n_target_nucleons * util.m_mc_pot);

    //========================================================================
    // Compare with GenieXSecExtractor
    //========================================================================
    double mc_integral = h_mc_cross_section->Integral();

    if (var->Name() == "pmu") {
      TFile fin_gxse(
          "/minerva/app/users/granados/cmtuser/MATAna/"
          "cc-ch-pip-ana/GENIEXSECEXTRACT_MCME1A.root");

      // Need to convert GXSE result from GeV --> MeV and rescale the binwidth

      PlotUtils::MnvH1D* pmu_xsec_dummy =
          (PlotUtils::MnvH1D*)fin_gxse.Get("pmu_xsec");
      assert(pmu_xsec_dummy);
      PlotUtils::MnvH1D* pmu_xsec =
          (PlotUtils::MnvH1D*)h_mc_cross_section->Clone(uniq());
      pmu_xsec->Reset();
      PlotUtils::MnvH1D* integrate_flux_dummy =
          (PlotUtils::MnvH1D*)fin_gxse.Get("reweightedflux_integrated");
      assert(integrate_flux_dummy);
      PlotUtils::MnvH1D* integrate_flux =
          (PlotUtils::MnvH1D*)h_flux_normalization->Clone(uniq());
      integrate_flux->Reset();
      PlotUtils::MnvH1D* unfolded_dummy =
          (PlotUtils::MnvH1D*)fin_gxse.Get("unfolded");
      assert(unfolded_dummy);
      unfolded_dummy->Scale(1. / 1000.);
      PlotUtils::MnvH1D* unfolded =
          (PlotUtils::MnvH1D*)true_effden->Clone(uniq());
      unfolded->Reset();
      std::cout << "xSec table\n"
                << "Bins "
                << "GXSecEx "
                << " MC \n";
      for (int i = 0; i < pmu_xsec->GetNbinsX() + 1; ++i) {
        std::cout << i << "  " << pmu_xsec_dummy->GetBinLowEdge(i) << "  "
                  << pmu_xsec_dummy->GetBinContent(i) << "  "
                  << h_mc_cross_section->GetBinContent(i) << "\n";
        pmu_xsec->SetBinContent(i, pmu_xsec_dummy->GetBinContent(i));
      }
      std::cout << "Integrated Flux table\n"
                << "Bins "
                << "GXSecEx "
                << " MC \n";
      for (int i = 0; i < integrate_flux->GetNbinsX() + 1; ++i) {
        std::cout << i << "  " << integrate_flux_dummy->GetBinLowEdge(i) << "  "
                  << integrate_flux_dummy->GetBinContent(i) << "  "
                  << mc_integrate_flux->GetBinContent(i) << "\n";
        integrate_flux->SetBinContent(i,
                                      integrate_flux_dummy->GetBinContent(i));
      }
      std::cout << "Evets table\n"
                << "Bins "
                << "GXSecEx "
                << " MC \n";
      for (int i = 0; i < unfolded->GetNbinsX() + 1; ++i) {
        std::cout << i << "  " << unfolded_dummy->GetBinLowEdge(i) << "  "
                  << unfolded_dummy->GetBinContent(i) << "  "
                  << true_effden->GetBinContent(i) << "  "
                  << true_effden->GetBinContent(i) /
                         unfolded_dummy->GetBinContent(i)
                  << "\n";
        unfolded->SetBinContent(i, unfolded_dummy->GetBinContent(i));
      }
      // Fixing effect of binwidth normalization with diferent units.
      pmu_xsec->Scale(1. / 1000.);
      // What units is the flux in?
      // Maybe we need to convert flux units from nu/cm^2/POT to nu/m^2/POT?
      // pmu_xsec->Scale(1./10000. );
      // Compare integrals
      double gxse_integral = pmu_xsec->Integral();
      std::cout << "  mc xsec integral = " << mc_integral << "\n";
      std::cout << "  gxse xsec integral = " << gxse_integral << "\n";
      std::cout << "  gxse / mc = " << gxse_integral / mc_integral << "\n";
      std::cout << "  Entries mc xsec = " << h_mc_cross_section->GetEntries()
                << "\n";
      std::cout << "  Entries gxse xsec = " << pmu_xsec_dummy->GetEntries()
                << "\n";

      // Area normalize
      //     h_mc_cross_section->Scale(1. / mc_integral);
      //     pmu_xsec->Scale(1. / gxse_integral);

      // plot on top of each other
      PlotTogether(mc_integrate_flux, "mc_flux", integrate_flux, "gxse_flux",
                   "Flux_Compare", 1.e+43, false, false, "flux*POT*#C12");
      PlotTogether(h_mc_cross_section, "mc", pmu_xsec, "gxse",
                   "gxse_compare_pmu");
      PlotTogether(true_effden, "mc_unfolded", unfolded, "gxse_Unfolded",
                   "Unfolded_compare");

      // Plot ratio
      PlotRatio(true_effden, unfolded, Form("%s", var->Name().c_str()), 1.,
                "Unfolded_compare", false);
      PlotRatio(h_mc_cross_section, pmu_xsec, Form("%s", var->Name().c_str()),
                1., "GXSEClosure", false);
    }

    //========================================================================
    // Compare EVENT RATE with GenieXSecExtractor
    //========================================================================
    // Start with "fake efficiency corrected data"
    // I've made sure that this matches what we get when we unfold and effcor
    // the bgsub.
    if (var->Name() == "pmu") {
      PlotUtils::MnvH1D* h_all_signal_true =
          (PlotUtils::MnvH1D*)true_var->m_hists.m_effden.hist->Clone(uniq());

      TFile fin_gxse(
          "/minerva/app/users/granados/cmtuser/MATAna/"
          "cc-ch-pip-ana/GENIEXSECEXTRACT_MCME1A.root");

      // Need to convert GXSE result from GeV --> MeV
      PlotUtils::MnvH1D* pmu_rate_dummy =
          (PlotUtils::MnvH1D*)fin_gxse.Get("pmu_xsec_evRate");
      PlotUtils::MnvH1D* pmu_rate =
          (PlotUtils::MnvH1D*)h_all_signal_true->Clone(uniq());
      pmu_rate->Reset();
      for (int i = 0; i < pmu_rate->GetNbinsX() + 1; ++i) {
        std::cout << i << "  " << pmu_rate_dummy->GetBinLowEdge(i) << "  "
                  << pmu_rate_dummy->GetBinContent(i) << "\n";
        pmu_rate->SetBinContent(i, pmu_rate_dummy->GetBinContent(i));
      }

      // BWN
      pmu_rate->Scale(1., "width");
      h_all_signal_true->Scale(1., "width");

      // Compare integrals
      double mc_rate_integral = h_all_signal_true->Integral();
      double gxse_rate_integral = pmu_rate->Integral();
      std::cout << "  mc rate integral = " << mc_rate_integral << "\n";
      std::cout << "  gxse rate integral = " << gxse_rate_integral << "\n";
      std::cout << "  gxse / mc = " << gxse_rate_integral / mc_rate_integral
                << "\n";
      std::cout << "  Entries mc xsec = " << h_all_signal_true->GetEntries()
                << "\n";
      std::cout << "  Entries gxse rate xsec = " << pmu_rate_dummy->GetEntries()
                << "\n";

      // Area normalize
      h_all_signal_true->Scale(1. / mc_rate_integral);
      pmu_rate->Scale(1. / gxse_rate_integral);

      // plot on top of each other
      PlotTogether(h_all_signal_true, "mc", pmu_rate, "gxse",
                   "gxse_ratecompare_pmu");

      // Plot ratio
      //        PlotRatio1(h_all_signal_true, pmu_rate,
      //        Form("GXSERateClosure_%s", name), true);

      PlotRatio(h_all_signal_true, pmu_rate, Form("GXSERateClosure_%s", name),
                1., "", true);
    }
  }

  //============================================================================
}

#endif  // GXSEClosure_C

#ifndef crossSectionClosure_C
#define crossSectionClosure_C

#include "PlotUtils/MnvH1D.h"
#include "TFile.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "makeCrossSectionMCInputs.C" // GetAnalysisVariables
#include "includes/CCPiMacroUtil.h"
#include "plotting_functions.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/TargetUtils.h"


//==============================================================================
// Main
// TODO function-ify the real xsec script and use those functions.
// Instead of copy-paste.
//==============================================================================
void crossSectionClosure(int signal_definition_int = 0) {
  // In and outfiles
    //TFile fin("rootfiles/MCXSecInputs_20190903.root", "READ");
    TFile fin("rootfiles/MCXSecInputs_20190904_ME1A.root", "READ");
    cout << "Reading input from " << fin.GetName() << endl;

  // Set up macro utility object -- which does the systematics for us
    const char* plist = "ALL";
    bool do_data = false, do_mc = false, do_truth = false;
    bool do_systematics = true, do_grid = false;
    CCPiMacroUtil util(signal_definition_int, plist, do_data, do_mc, do_truth,
                       do_systematics, do_grid);

  // Set POT
  PlotUtils::MnvH1D* h_mc_pot=(PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float mc_pot = h_mc_pot->GetBinContent(1);
  util.m_mc_pot = mc_pot;

  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables = GetAnalysisVariables(util.m_signal_definition, 
                                                          do_truth_vars);

  for (auto var : variables)
    var->LoadMCHistsFromFile(fin, util.m_error_bands);

  for (auto var : variables) {
    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    if (var->m_is_true) continue;
    if (var->Name() == std::string("tpi_mbr"))  continue;
    if (var->Name() == sidebands::kFitVarString) continue; 

    if(var->Name() != "pmu") continue;

    const char* name = var->Name().c_str();

    // We'll be needing the true version of this variable
    Variable* true_var = GetVar(variables, var->Name() + std::string("_true"));

    // BG-level closure
    // A (mc selection - untuned bg)
    // =?
    // B (signal only (reco))
      PlotUtils::MnvH1D* untuned_bg  = 
          (PlotUtils::MnvH1D*)var->m_hists.m_bg.hist->Clone(uniq());
      PlotUtils::MnvH1D* bgsub_mc    = 
          (PlotUtils::MnvH1D*)var->m_hists.m_selection_mc.hist->Clone(uniq());

      bgsub_mc->Add(untuned_bg, -1);

      PlotUtils::MnvH1D* signal_only = 
          (PlotUtils::MnvH1D*)var->m_hists.m_effnum.hist->Clone(uniq());

      PlotRatio(bgsub_mc, signal_only, Form("BGClosure_%s",name));

      // do the same thing in true vars
      PlotUtils::MnvH1D* untuned_bg_true = 
          (PlotUtils::MnvH1D*)true_var->m_hists.m_bg.hist->Clone(uniq());
      PlotUtils::MnvH1D* bgsub_true = 
          (PlotUtils::MnvH1D*)true_var->m_hists.m_selection_mc.hist->Clone(uniq());

      bgsub_true->Add(untuned_bg_true, -1);

      PlotUtils::MnvH1D* signal_only_true = 
          (PlotUtils::MnvH1D*)true_var->m_hists.m_effnum.hist->Clone(uniq());

      PlotRatio(bgsub_true, signal_only_true, Form("TrueBGClosure_%s",name));

      


    // Unfolding closure
    // C (Unfold(A) == Unfold(mc selection - untuned bg))
    // =?
    // D (signal only (truth))
      MinervaUnfold::MnvUnfold mnv_unfold;
      mnv_unfold.setUseBetterStatErrorCalc(true);

      PlotUtils::MnvH2D* migration = 
          (PlotUtils::MnvH2D*)var->m_hists.m_migration.hist->Clone(uniq());

      // unfold
      mnv_unfold.UnfoldHisto(var->m_hists.m_unfolded, migration, bgsub_mc,
                             RooUnfold::kBayes, 4);

      //compare
      PlotRatio(var->m_hists.m_unfolded, signal_only_true, Form("UnfoldingClosure_%s", name));


    // Efficiency/Norm closure
    // E (Unfold(A)/efficiency)
    // =?
    // D (efficiency denominator)
      // Calculate efficiency
      PlotUtils::MnvH1D* efficiency = 
          (PlotUtils::MnvH1D*)true_var->m_hists.m_effnum.hist->Clone(uniq());

      efficiency->Divide(true_var->m_hists.m_effnum.hist,
                         true_var->m_hists.m_effden.hist);

      PlotUtils::MnvH1D* h_efficiency_corrected_data = (PlotUtils::MnvH1D*)var->m_hists.m_unfolded->Clone(uniq());

      // Efficiency correct
      h_efficiency_corrected_data->Divide(var->m_hists.m_unfolded, efficiency);

      // Flux normalization
      PlotUtils::MnvH1D* h_flux_normalization = 
          (PlotUtils::MnvH1D*)h_efficiency_corrected_data->Clone("flux_normalization");
      h_flux_normalization->ClearAllErrorBands();
      h_flux_normalization->Reset();

      // Get the flux histo, to be integrated
      const bool use_hundred_universes = true;
      static PlotUtils::FluxReweighter* frw = 
          new PlotUtils::FluxReweighter( 14, false, "minervame1a",
                                         PlotUtils::FluxReweighter::gen2thin, 
                                         PlotUtils::FluxReweighter::g4numiv6,
                                         use_hundred_universes );

      h_flux_normalization = 
          frw->GetIntegratedFluxReweighted(14, h_efficiency_corrected_data, 0., 100.);

      // remove redundant error bands
      h_flux_normalization->PopVertErrorBand("Flux_BeamFocus");
      h_flux_normalization->PopVertErrorBand("ppfx1_Total");   

      // Convert flux units from nu/m^2/POT to nu/cm^2/POT
      h_flux_normalization->Scale( 1.0e-4 );

      // Divide flux integral
      PlotUtils::MnvH1D* h_cross_section = 
          (PlotUtils::MnvH1D*)h_efficiency_corrected_data->Clone(uniq());
      h_cross_section->AddMissingErrorBandsAndFillWithCV(*h_flux_normalization);
      h_cross_section->Divide( h_cross_section, h_flux_normalization );

      // targets and POT norm
      static const double apothem    = 865.;
      static const double upstream   = 5900.; // ~module 25 plane 1
      static const double downstream = 8430.; // ~module 81 plane 1

      double n_target_nucleons = 
          PlotUtils::TargetUtils::Get().GetTrackerNNucleons(upstream, downstream,
                                                            /*isMC =*/ false,
                                                            apothem);

      std::cout << "  flux_integral cv = "  << h_flux_normalization->GetBinContent(1) << "\n";
      std::cout << "  N target nucleons = " << n_target_nucleons                      << "\n";
      std::cout << "  mc pot = " << util.m_mc_pot << "\n";
      std::cout << "  flux_integral * N target nucleons = " << n_target_nucleons * h_flux_normalization->GetBinContent(1) << "\n";
      std::cout << "  flux_integral * N target nucleons * mc pot= " << n_target_nucleons * h_flux_normalization->GetBinContent(1) * util.m_mc_pot << "\n";
      static const double data_scale = 1.0 / ( n_target_nucleons /* * util.m_data_pot*/);
      //static const double data_scale = 1.0 / 3.529606e+42;
      h_cross_section->Scale( data_scale, "width" );

      // Scale MC
      PlotUtils::MnvH1D* h_mc_cross_section = 
          (PlotUtils::MnvH1D*)true_var->m_hists.m_effden.hist->Clone(uniq());

      h_mc_cross_section->AddMissingErrorBandsAndFillWithCV(*h_flux_normalization);
      h_mc_cross_section->Divide( h_mc_cross_section, h_flux_normalization );

      static const double mc_scale = 1.0 / ( n_target_nucleons /* * util.m_mc_pot*/ );
      h_mc_cross_section->Scale(mc_scale, "width");

      double fake_data_integral = h_cross_section->Integral();
      double mc_integral = h_mc_cross_section->Integral();
      std::cout << "  fake data xsec integral   = " << fake_data_integral << "\n";
      std::cout << "  mc xsec integral = " << mc_integral << "\n";

      // compare
      PlotRatio(h_cross_section, h_mc_cross_section, Form("CrossSectionClosure_%s", name));

      /*
      // compare with GenieXSecExtractor results
      if(var->Name() == "pmu") {
        TFile fin_gxse("/minerva/app/users/bmesserl/cmtuser/Minerva_v21r1p1/GENIEXSecExtract/GenieXSecs_20190904_ME1A_pmu.root");

        // GeV --> MeV
        PlotUtils::MnvH1D* pmu_xsec_dummy = (PlotUtils::MnvH1D*)fin_gxse.Get("ds_dpmu_xsec");
        PlotUtils::MnvH1D* pmu_xsec = (PlotUtils::MnvH1D*)h_mc_cross_section->Clone(uniq());
        pmu_xsec->Reset();
        for(int i = 0; i < pmu_xsec->GetNbinsX()+1; ++i) {
          std::cout << i << "  " << pmu_xsec_dummy->GetBinLowEdge(i) << "  " << pmu_xsec_dummy->GetBinContent(i) << "\n";
          pmu_xsec->SetBinContent(i, pmu_xsec_dummy->GetBinContent(i));
        }

        double gxse_integral = pmu_xsec->Integral();
        double xsec_integral = h_cross_section->Integral();
        std::cout << "  mc xsec integral   = " << xsec_integral << "\n";
        std::cout << "  gxse xsec integral = " << gxse_integral << "\n";

        // Area norm comparison
        h_cross_section->Scale(1./xsec_integral);
        pmu_xsec->Scale(1./gxse_integral);
        
        PlotTogether(h_cross_section, "mc", pmu_xsec, "gxse", "gxse_compare_pmu");
        PlotRatio(h_cross_section, pmu_xsec, Form("GXSEClosure_%s", name));
      }
      */
  }

  //============================================================================
}

#endif // crossSectionClosure_C

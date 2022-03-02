// This script reads in root files containing migration matrices and it plots
// the migration matrices. The intent is to tweak binning in order to get
// approximately diagonal migration matrices. Approximately diagonal migration
// matrix is a metric for evaluating the reliability of your binning choice.
#ifndef binningStudy_C
#define binningStudy_C

#include "PlotUtils/MnvH1D.h"
#include "TFile.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "makeCrossSectionMCInputs.C" // GetAnalysisVariables
#include "includes/MacroUtil.h"
#include "plotting_functions.h"
//#include "includes/common_stuff.h" // SetBinVec

// TH1D::Rebin(int ngroup, const char* newname, double* xbins)
// A new histogram is created (you should specify newname). The parameter ngroup
// is the number of variable size bins in the created histogram. The array xbins
// must contain ngroup+1 elements that represent the low-edges of the bins. If the
// original histogram has errors stored (via Sumw2), the resulting histograms has
// new errors correctly calculated.



//==============================================================================
// Main
//==============================================================================
void binningStudy(int signal_definition_int = 0) {
  // In and outfiles
    //TFile fin("rootfiles/MCXSecInputs_20190616_FineBins.root", "READ");
    TFile fin("MCXSecInputs_20220225.root", "READ");
    cout << "Reading input from " << fin.GetName() << endl;

  // Set up macro utility object -- which does the systematics for us
    const char* plist = "ME1L";
    std::string data_file_list = GetPlaylistFile(plist, false);
    std::string mc_file_list = GetPlaylistFile(plist, true);
    bool do_data = false, do_mc = false, do_truth = false;
    bool do_systematics = true, do_grid = false;
    CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                         plist, do_truth, do_grid, do_systematics);
  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables = GetAnalysisVariables(util.m_signal_definition, 
                                                          do_truth_vars);

  ContainerEraser::erase_if(variables, [](Variable* v) { return v->Name() == "wexp_fit"; });

  for (auto var : variables) {
    if (var->m_is_true)
      continue;

    //if (var->Name()!="pmu")
    //  continue;

    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    PlotUtils::MnvH1D* selection = (PlotUtils::MnvH1D*)var->m_hists.m_selection_mc.hist->Clone(uniq());
    PlotUtils::MnvH2D* mig = (PlotUtils::MnvH2D*)var->m_hists.m_migration.hist->Clone(uniq());


    std::cout << "Plotting stuff\n";
    const bool do_frac_unc  = true;
    const bool include_stat = false;
    const bool do_cov_area_norm   = false;

    EventSelectionPlotInfo plot_info(var, util.m_mc_pot, util.m_data_pot,
        do_frac_unc, do_cov_area_norm, include_stat, util.m_signal_definition);

    //PlotDataMCWithError(eff, nullptr, plot_info, "EffWError");

    double ymax = -1.;
    bool do_log_scale = false;
    bool do_bg = true;
    bool do_tuned_bg = true;

    PlotMC(selection, plot_info, Form("Selection_%s",  var->Name().c_str()), -1., "N Events");
    PlotStatError(selection, plot_info, Form("StatError_%s",  var->Name().c_str()), -1., "stat error (%)");

    PlotMigration_VariableBins(mig, var->Name());
    PlotMigration_AbsoluteBins(mig, var->Name());

    //Plot_ErrorSummary(plot_info, selection, "Sel");

    // Plot Rebinned
    /*
    {
      std::vector<double> rebins_vec;
      SetBinVec(var->Name(), rebins_vec);
      double rebins_array[rebins_vec.size()];
      std::copy(rebins_vec.begin(), rebins_vec.end(), rebins_array);
      const int nbins = sizeof(rebins_array)/sizeof(*rebins_array) - 1 ;
      std::sort(rebins_array, rebins_array+nbins);
      PlotUtils::MnvH1D* hist = (MnvH1D*)selection->Rebin(nbins, "", rebins_array);

      // after rebinning
      //PlotMC(hist, plot_info, Form("Selection_%s_rebin",  var->Name().c_str()), -1., "N Events");
      PlotStatError(hist, plot_info, Form("StatError_%s_rebin",  var->Name().c_str()), -1., "stat error (%)");
      Plot_ErrorSummary(plot_info, hist, "SelRebin");

      // after rebinning + bin width norm
      hist->Scale(1., "width");
      //PlotMC(hist, plot_info, Form("Selection_%s_rebin_bwn",  var->Name().c_str()), -1., "N Events / MeV");
      //PlotStatError(hist, plot_info, Form("StatError_%s_rebin_bwn",  var->Name().c_str()), -1., "stat error (%)");


      //// Migration
      //{
      //  if(var->m_is_true)
      //    continue;

      //  PlotUtils::MnvH2D* mig = (PlotUtils::MnvH2D*)var->m_hists.m_migration.hist->Clone(uniq());
      //  PlotMigration_AbsoluteBins(mig, var->Name());
      //  PlotMigration_VariableBins(mig, var->Name());
      //}
    }
    */

    //PlotMC(eff,              plot_info, Form("Efficiency_%s", var->Name().c_str()), 0.2, "Efficiency");
    //PlotEfficiency_ErrorSummary(plot_info);
  }

  //============================================================================
}

#endif // binningStudy_C

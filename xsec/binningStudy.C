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
#include "plotting_functions2D.h"
#include "includes/Binning.h"
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
    TFile fin("MCXSecInputs_20231001_ME1A_untrackedpion_15thetabin_noSys.root", "READ");
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
  std::vector<Variable2D*> variables2D =
      GetAnalysisVariables2D(util.m_signal_definition, do_truth_vars);

  ContainerEraser::erase_if(variables, [](Variable* v) { return v->Name() == "wexp_fit"; });

  const bool do_frac_unc  = true;
  const bool include_stat = false;
  const bool do_cov_area_norm   = false;

  for (auto var : variables) {
    if (var->m_is_true)
      continue;

    //if (var->Name()!="pmu")
    //  continue;

    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    PlotUtils::MnvH1D* selection = (PlotUtils::MnvH1D*)var->m_hists.m_selection_mc.hist->Clone(uniq());
    PlotUtils::MnvH2D* mig = (PlotUtils::MnvH2D*)var->m_hists.m_migration.hist->Clone(uniq());


    std::cout << "Plotting stuff\n";

    const Plotter plot_info(var, util.m_mc_pot, util.m_data_pot,
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

  for (auto var2D : variables2D) {
    const EventSelectionPlotInfo2D plot_info2D(
          var2D, util.m_mc_pot, util.m_data_pot, do_frac_unc, do_cov_area_norm,
          include_stat, util.m_signal_definition);   
    var2D->LoadMCHistsFromFile(fin, util.m_error_bands);
    if (!var2D->m_is_true){
        PlotUtils::MnvH2D* mig2D =
            (PlotUtils::MnvH2D*)var2D->m_hists2D.m_migration.hist->Clone(uniq());               
      int nbinsx = var2D->NBinsX();
      int xmtxbins = nbinsx + 2;
      int nbinsy = var2D->NBinsY();
      int ymtxbinx = nbinsy + 2;
      double sumreco, sumtrue;
      int x, y;
      for (int bin = 1; bin <= nbinsy; ++bin){
        PlotUtils::MnvH2D* migration_reco = 
                          new PlotUtils::MnvH2D(Form("MigMtx_reco_%s_vs_%s_Bin%i", var2D->NameX().c_str(),
                          var2D->NameY().c_str(), bin), Form("%s_%s", var2D->NameX().c_str(),
                          var2D->NameY().c_str()), nbinsx, 0.0, (double)nbinsx,
                          nbinsx, 0.0, (double)nbinsx); 
        PlotUtils::MnvH2D* migration_true = 
                          new PlotUtils::MnvH2D(Form("MigMtx_true_%s_vs_%s_Bin%i", var2D->NameX().c_str(),
                          var2D->NameY().c_str(), bin), Form("%s_%s", var2D->NameX().c_str(),
                          var2D->NameY().c_str()), nbinsx, 0.0, (double)nbinsx,
                          nbinsx, 0.0, (double)nbinsx); 
        for (int j = 1; j <= nbinsx; ++j) { //it is going row by row in the
 					    //internal matrix
 	  y = bin*xmtxbins + j;
          for (int i = 1; i <= nbinsx; ++i){ //It is going columb by columb
 						 //in the internal matrix 
            x = bin*xmtxbins + i;
            sumreco = 0;
            sumtrue = 0; 
	    for (int c = 0; c < nbinsy; ++c){ //It is used to sum the 
						   //matrices vertically or horizontally
              sumreco += mig2D->GetBinContent(x, y + c*xmtxbins);
              sumtrue += mig2D->GetBinContent(x + c*xmtxbins, y);
            }
  	    migration_reco->SetBinContent(i, j, sumreco);
	    migration_true->SetBinContent(i, j, sumtrue);
          }        
        }
      // Here we print the resultant projected migration matrix
      PlotMigration_AbsoluteBins(migration_reco, Form("reco_%s_vs_%s_Bin%i", var2D->NameX().c_str(),var2D->NameY().c_str(), bin));
      PlotMigration_AbsoluteBins(migration_true, Form("true_%s_vs_%s_Bin%i", var2D->NameX().c_str(),var2D->NameY().c_str(), bin));
      }      
      PlotMigration2D(plot_info2D, mig2D, var2D->NameX(), var2D->NameY());
    }     
  }
  //============================================================================
}

#endif // binningStudy_C

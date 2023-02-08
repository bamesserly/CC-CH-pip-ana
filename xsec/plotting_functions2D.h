//#ifndef plotting_functions2D_h
//#define plotting_functions2D_h

#include <iostream>
#include <sstream>
#include <cmath>

#include "../includes/myPlotStyle.h"
#include "PlotUtils/MnvColors.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/HistogramUtils.h"
#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif
#include "PlotUtils/MnvVertErrorBand.h"
#include "SignalDefinition.h"
#include "Systematics.h"  // namespace systematics
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TList.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TObject.h"
//#include "TStyle.h"
#include "TText.h"
#include "Variable.h"
#include "Variable2D.h"
#include "PlotUtils/GridCanvas.h"

class Variable;
class Variable2D;
std::map< std::string, std::vector<string>> error_names;
//==============================================================================
// Container class
//==============================================================================
class EventSelectionPlotInfo2D {
 public:

  // Constructor with Var2D
  EventSelectionPlotInfo2D(Variable2D* variable2D, float mc_pot, float data_pot,
                         bool do_frac_unc, bool do_cov_area_norm,
                         bool include_stat, SignalDefinition signal_definition)
      : m_mnv_plotter(kCCNuPionIncStyle),
        m_variable2D(variable2D),
        m_mc_pot(mc_pot),
        m_data_pot(data_pot),
        m_do_frac_unc(do_frac_unc),
        m_do_cov_area_norm(do_cov_area_norm),
        m_include_stat(include_stat),
//  	m_error_names(),
        m_signal_definition(signal_definition) {
    m_do_frac_unc_str = m_do_frac_unc ? "Frac" : "Abs";
    m_do_cov_area_norm_str = m_do_cov_area_norm ? "CovAreaNorm" : "";
  }

  // Constructor without var
  EventSelectionPlotInfo2D(float mc_pot, float data_pot, bool do_frac_unc,
                         bool do_cov_area_norm, bool include_stat,
                         SignalDefinition signal_definition)
      : m_mnv_plotter(kCCNuPionIncStyle),
        m_variable2D(nullptr),
        m_mc_pot(mc_pot),
        m_data_pot(data_pot),
        m_do_frac_unc(do_frac_unc),
        m_do_cov_area_norm(do_cov_area_norm),
        m_include_stat(include_stat),
//      m_error_names(),
        m_signal_definition(signal_definition) {
    m_do_frac_unc_str = m_do_frac_unc ? "Frac" : "Abs";
    m_do_cov_area_norm_str = m_do_cov_area_norm ? "CovAreaNorm" : "";
  }

  // Members
  MnvPlotter m_mnv_plotter;
  Variable2D* m_variable2D;
  float m_mc_pot;
  float m_data_pot;
  bool m_do_frac_unc;
  bool m_do_cov_area_norm;
  bool m_include_stat;
  SignalDefinition m_signal_definition;
  std::string m_do_frac_unc_str;
  std::string m_do_cov_area_norm_str;
  // Add X label
  void SetXLabel(PlotUtils::MnvH2D* hist) {
    std::string label =
        m_variable2D->m_hists2D.m_labelX + " (" + m_variable2D->m_unitsX + ")";
    if (hist) hist->GetXaxis()->SetTitle(label.c_str());
  }

  // Add X label
  void SetXLabel(TH2* hist) {
    std::string label =
        m_variable2D->m_hists2D.m_labelX + " (" + m_variable2D->m_unitsX + ")";
    if (hist) hist->GetXaxis()->SetTitle(label.c_str());
  }

  // Add title to the MnvPlotter
  void SetTitle() {
    if (!gPad)
      throw std::runtime_error("Need a TCanvas. Please make one first.");
    std::string title = GetSignalName(m_signal_definition);
    m_mnv_plotter.AddHistoTitle(title.c_str());
  }

  void SetTitle(std::string title) {
    if (!gPad)
      throw std::runtime_error("Need a TCanvas. Please make one first.");
    m_mnv_plotter.AddHistoTitle(title.c_str());
  }
};
//-------------------------------------
//--------------Multipliers------------
//-------------------------------------
std::vector<double> GetMultipliers(TH2D* data_hist){
  std::vector<double> multipliers;
  double Max = data_hist->GetMaximum();
  double currentvalue = 0;
  for (int j = 1; j <= data_hist->GetNbinsY(); ++j){
    currentvalue = 0;
    for(int i = 1; i <= data_hist->GetNbinsX(); ++i ){
      if (data_hist->GetBinContent(i,j) > currentvalue)
          currentvalue = data_hist->GetBinContent(i,j);  
    }
    if (currentvalue == 0) multipliers.push_back(1.);
    else multipliers.push_back((double)nearbyint(Max/currentvalue));       
  }
  return multipliers;
}
/*
void AddToTmp( TObject* obj )
{
    //note: build this everytime to avoid obj unused compiler warning
    if ( 0 == obj )
    {
        Warning( "MnvPlotter::AddToTmp", "Attempting to add NULL object to garbage.");
        return;
    }

//#if DO_GARBAGE_COLLECTION

    //Error("MnvPlotter::AddToTmp", Form( "Adding object '%s' at index %d.", obj->GetName(), fTmpObjects.GetEntries() ) );

    //! Add to array if the object isn't in the array already
    if ( ! fTmpObjects.FindObject( obj ) )
        fTmpObjects.AddLast( obj );
//#endif
}
*/
std::vector<TH2D*> DrawErrorSummary2D(
        EventSelectionPlotInfo2D p,
        PlotUtils::MnvH2D* h,
   //     const std::string& legPos   = "TR",
        const bool   includeStat      = true,
        const bool   solidLinesOnly   = true,
        const double ignoreThreshold  = 0.00001,
        const bool covAreaNormalize = false,
        const std::string& errorGroupName  = "",
        const bool  asfrac   = false,
        const std::string &Ytitle = "",
        bool ignoreUngrouped = false,
	const double axis_minimum = MnvHist::AutoAxisLimit,
        const double axis_maximum = MnvHist::AutoAxisLimit
        )
{
    if ( ! h )
    {
	TH2D* noHist = NULL;
        std::vector<TH2D*> noHists;
	noHists.push_back(noHist);
        Error("DrawErrorSummary2D", "You passed me a NULL MnvH2D.  Nothing to do.");
        return noHists;
    }

    //set max digits to the default, since we almost never want scientific notation for errors.
    //restore the setting before returning.
    const int oldMaxDigits = TGaxis::GetMaxDigits();
    int axis_max_digits = 3;
    TGaxis::SetMaxDigits( axis_max_digits );

    // Store the pieces for a legend
    vector<TH2D*>   hists;
    vector<string> names;
    vector<string> opts;
    std::string stat_error_name = "Statistical";
    bool useDifferentLineStyles = !solidLinesOnly;

    //! Get the total error and apply styles
    TH2D *hTotalErr = (TH2D*)h->GetTotalError( includeStat, asfrac, covAreaNormalize ).Clone( Form("h_total_err_errSum_%d", __LINE__) );
    //AddToTmp( hTotalErr );

    //Error("DrawDataMCErrorSummary", Form("Total Err Max: %.2f",hTotalErr->GetMaximum() ) );

    p.m_mnv_plotter.ApplyNextLineStyle( hTotalErr, true, useDifferentLineStyles );

    //respect max/min setting the user may have used
    if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
        hTotalErr->SetMinimum( 0. );
    else
        hTotalErr->SetMinimum( axis_minimum );

    if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
        hTotalErr->SetMaximum( p.m_mnv_plotter.headroom * hTotalErr->GetMaximum() );
    else
        hTotalErr->SetMaximum( axis_maximum );

    hTotalErr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
    if (!asfrac )hTotalErr->GetYaxis()->SetTitle( Ytitle.c_str() );
    if (errorGroupName == "") {
        if (!asfrac) hTotalErr->Scale( hTotalErr->GetBinWidth(1), "width" );
//        hTotalErr->Draw( "HIST" );
        const string totalName = ( includeStat ? "Total Uncertainty" : "Total Sys. Uncertainty" );
        hists.push_back( hTotalErr );
        names.push_back( totalName );
        opts.push_back( "l" );
    }

    if ( includeStat && errorGroupName == "")
    {
        TH2D *statErr = (TH2D*)h->GetStatError(asfrac).Clone( Form("this_stat_err_%d", __LINE__) );
        //AddToTmp( statErr );

        statErr->SetLineColor( 12 );//dark gray
        statErr->SetLineStyle( 2 ); //dashed
        statErr->SetLineWidth( 3 );
        statErr->Draw("HIST SAME");
        hists.push_back( statErr );
        names.push_back( stat_error_name );
        opts.push_back( "l" );
    }

    TH2D *hTmpErr = (TH2D*)hTotalErr->Clone( Form("h_tmp_err_errSum_%d", __LINE__) );
    hTmpErr->Reset();
    map<string,TH2D*> errGroupHists;

    // plot each of the fractional contributions from all the errors.
    // first, we make a list of error bands to plot...
    bool drawn_already = (errorGroupName == "") ? true : false;
    vector<string> errNames = h->GetVertErrorBandNames();
    vector<string> otherNames = h->GetLatErrorBandNames();
    errNames.insert( errNames.end(), otherNames.begin(), otherNames.end() );
    otherNames = h->GetUncorrErrorNames();
    errNames.insert( errNames.end(), otherNames.begin(), otherNames.end() );
//    PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
    for ( vector<string>::const_iterator it_name = errNames.begin();
            it_name != errNames.end();
            ++it_name)
    {
      TH2D * hErr = NULL;

      if (h->HasVertErrorBand(*it_name))
        hErr = dynamic_cast<TH2D*>(h->GetVertErrorBand(*it_name)->GetErrorBand( asfrac, covAreaNormalize ).Clone( Form("tmp_vertError_%s", (*it_name).c_str()) ));
      else if (h->HasLatErrorBand(*it_name))
        hErr = dynamic_cast<TH2D*>(h->GetLatErrorBand(*it_name)->GetErrorBand( asfrac, covAreaNormalize ).Clone( Form("tmp_latError_%s", (*it_name).c_str()) ));
      else
        throw std::runtime_error( Form("MnvPlotter::DrawErrorSummary(): Couldn't determine error band type for error name '%s'", (*it_name).c_str()) );

      //is this histogram part of a group?
      bool inGroup = false;
      for (MnvPlotter::ErrorSummaryGroupMap::const_iterator itGroup = p.m_mnv_plotter.error_summary_group_map.begin();
              itGroup != p.m_mnv_plotter.error_summary_group_map.end(); ++itGroup ) {
        const string& errName           = itGroup->first;
        const vector<string>& histNames = itGroup->second;

        //if this histogram is not in the group we're considering, skip to the next one
        if ( find( histNames.begin(), histNames.end(), *it_name) == histNames.end() )
          continue;
        //std::cout << " MnvPlotter found " << errName << " " << *it_name << std::endl;
        // if plotting the errors from only one group,
        // we don't want to "sub-group" them any further.
        // therefore we don't do any adding of histograms.
        if (errorGroupName==errName) {
          inGroup=true;
          break;
        }
        // otherwise, if no group was specifically chosen,
        // then (since this error band is already known to be in this group)
        // we need to add it into the histogram for the group
        // (or create that histogram if it doesn't exist yet).
        else if (errorGroupName == "") {
          map<string,TH2D*>::iterator itGroupHist = errGroupHists.find(errName);

          if ( errGroupHists.end() == itGroupHist ) {
            errGroupHists[ errName ] = hErr;
          }
          else {
            MnvHist::AddInQuadrature( itGroupHist->second, hErr );
            delete hErr;
          }

          inGroup = true;
          break;
        }
      }

      // if we haven't selected a group whose constituents we want to see...
      if ( errorGroupName=="" ) {
        // we never want to show the individual plots
        // when a histogram was included in a group:
        // the sums will be drawn later.
        if (inGroup)
          continue;

        // ... when we are ignoring ungrouped errors,
        // then NOTHING gets drawn here:
        // the grouped errors were added to the group histogram
        // above (and so we don't want them here),
        // and the ungrouped ones we're ignoring altogether.
        if (ignoreUngrouped)
          continue;
      }
      // if we DID select a group to draw the constituents of,
      // and this histogram isn't part of it,
      // then we also don't want to see it.
      else if ( errorGroupName != "" && !inGroup)
        continue;

      //AddToTmp( hErr );

      // don't draw any errors that have no bins above threshold
      if ( 0 < ignoreThreshold && hErr->GetBinContent( hErr->GetMaximumBin() ) < ignoreThreshold )
        continue;

      p.m_mnv_plotter.ApplyNextLineStyle( hErr, false, useDifferentLineStyles);

      hErr->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());

      map<string,int>::const_iterator itErrCol = p.m_mnv_plotter.error_color_map.find( *it_name );
      if ( p.m_mnv_plotter.error_color_map.end() != itErrCol )
        hErr->SetLineColor( itErrCol->second );

      hErr->SetLineWidth( p.m_mnv_plotter.mc_line_width );
      if (drawn_already) {
        hErr->Draw( "HIST SAME" );
      }
      else {
        drawn_already = true;
        hTmpErr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
        if (!asfrac) hTmpErr->GetYaxis()->SetTitle( Ytitle.c_str() );
        if (!asfrac) hTmpErr->Scale(hTmpErr->GetBinWidth(1),"width");
        hTmpErr->Draw("HIST");
	
        ////respect max/min setting the user may have used
        if ( MnvHist::IsAutoAxisLimit( axis_minimum ) )
	  hErr->SetMinimum( 0. );
        else
	  hErr->SetMinimum( axis_minimum );
	
	//if (inGroup) {
	//hErr->SetMaximum( headroom * axis_maximum_group );
	//}
        if ( MnvHist::IsAutoAxisLimit( axis_maximum ) )
	  hErr->SetMaximum( p.m_mnv_plotter.headroom * hTotalErr->GetMaximum() );
        else
	  hErr->SetMaximum( axis_maximum);

        hErr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
        if (!asfrac) hErr->GetYaxis()->SetTitle( Ytitle.c_str() );
        //if (!asfrac) hErr->Scale(hErr->GetBinWidth(1),"width");
        hErr->Draw("HIST SAME");
      }
      hists.push_back( hErr );
      names.push_back( *it_name );
      opts.push_back( "l" );
    }

    //add error groups
    for ( map<string,TH2D*>::iterator itGroup = errGroupHists.begin(); itGroup != errGroupHists.end(); ++itGroup )
    {
        //   std::cout << "  (plot " << h->GetName() << ") drawing error group '" << itGroup->first << "'" << std::endl;
        const string& name = itGroup->first;
        TH2D* hist = itGroup->second;

        if ( 0 < ignoreThreshold && hist->GetBinContent( hist->GetMaximumBin() ) < ignoreThreshold )
        {
            //     std::cout << "    (... ignored because its maximum (" << hist->GetBinContent( hist->GetMaximumBin() ) << ") was below threshold)" << std::endl;
            continue;
        }

        p.m_mnv_plotter.ApplyNextLineStyle( hist, false, useDifferentLineStyles);

        map<string,int>::const_iterator itErrCol = p.m_mnv_plotter.error_color_map.find( name );
        if ( p.m_mnv_plotter.error_color_map.end() != itErrCol )
            hist->SetLineColor( itErrCol->second );

        hist->SetLineWidth( p.m_mnv_plotter.mc_line_width );
        if (!asfrac)hist->Scale(hist->GetBinWidth(1),"width");
        hist->Draw( "HIST SAME" );
        hists.push_back( hist );
        names.push_back( name );
        opts.push_back( "l" );
    }
    std::vector<std::string> vecNames;
    if (errorGroupName == ""){
      error_names.insert(pair<std::string, std::vector<std::string>>("Totals", names));
      vecNames = error_names["Totals"]; 
    }
    else{
      error_names.insert(pair<std::string, std::vector<std::string>>(errorGroupName, names));
      vecNames = error_names[errorGroupName];
    }
/*  std::string legPos = "BR";
    if ( legPos != "N" )
    {
        size_t legendWidth = p.m_mnv_plotter.GetLegendWidthInLetters( names );
        double x1,y1,x2,y2;
        p.m_mnv_plotter.DecodeLegendPosition( x1, y1, x2, y2, legPos, hists.size(), legendWidth );
        p.m_mnv_plotter.AddPlotLegend( hists, names, opts, x1, y1, x2-x1, y2-y1, p.m_mnv_plotter.legend_text_size );
    }*/

    gPad->RedrawAxis();
    gPad->Update();

    TGaxis::SetMaxDigits( oldMaxDigits );

    return hists;
}



//==============================================================================
// Some Systematics General Functions
//==============================================================================
void SetErrorGroups2D(MnvPlotter& mnv_plotter) {
  mnv_plotter.error_summary_group_map.clear();
  mnv_plotter.error_summary_group_map["Flux"].push_back("Flux");
  mnv_plotter.error_summary_group_map["NonResPi"].push_back("GENIE_Rvn1pi");
  mnv_plotter.error_summary_group_map["NonResPi"].push_back("GENIE_Rvp1pi");
  mnv_plotter.error_summary_group_map["NonResPi"].push_back("GENIE_Rvn2pi");
  mnv_plotter.error_summary_group_map["NonResPi"].push_back("GENIE_Rvp2pi");
  mnv_plotter.error_summary_group_map["2p2h"].push_back("Low_Recoil_2p2h_Tune");
  mnv_plotter.error_summary_group_map["LowQ2Pi"].push_back("LowQ2Pi");
  mnv_plotter.error_summary_group_map["Muon"].push_back("Muon_Energy_MINOS");
  mnv_plotter.error_summary_group_map["Muon"].push_back("Muon_Energy_MINERvA");
  mnv_plotter.error_summary_group_map["Muon"].push_back(
      "Muon_Energy_Resolution");
  mnv_plotter.error_summary_group_map["Muon"].push_back(
      "MINOS_Reconstruction_Efficiency");
  mnv_plotter.error_summary_group_map["Muon"].push_back("MuonAngleXResolution");
  mnv_plotter.error_summary_group_map["Muon"].push_back("MuonAngleYResolution");
  mnv_plotter.error_summary_group_map["Muon"].push_back("MuonResolution");
  mnv_plotter.error_summary_group_map["PhysicsModel"].push_back(
      "MichelEfficiency");
  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_D2_MaRES");
  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_EP_MvRES");
  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_D2_NormCCRES");
  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_MaCCQE");
  mnv_plotter.error_summary_group_map["PhysicsModel"].push_back(
      "Target_Mass_CH");
  mnv_plotter.error_summary_group_map["PhysicsModel"].push_back(
      "Target_Mass_C");
  mnv_plotter.error_summary_group_map["PhysicsModel"].push_back(
      "Target_Mass_Fe");
  mnv_plotter.error_summary_group_map["PhysicsModel"].push_back(
      "Target_Mass_H2O");
  mnv_plotter.error_summary_group_map["PhysicsModel"].push_back(
      "Target_Mass_Pb");
  mnv_plotter.error_summary_group_map["Response"].push_back("response_em");
  mnv_plotter.error_summary_group_map["Response"].push_back("response_meson");
  mnv_plotter.error_summary_group_map["Response"].push_back("response_other");
  mnv_plotter.error_summary_group_map["Response"].push_back("response_proton");
  mnv_plotter.error_summary_group_map["Diffractive"].push_back(
      "DiffractiveModelUnc");
  mnv_plotter.error_summary_group_map["Diffractive"].push_back(
      "CoherentPiUnc_C");
  mnv_plotter.error_summary_group_map["Diffractive"].push_back(
      "CoherentPiUnc_CH");
  mnv_plotter.error_summary_group_map["Diffractive"].push_back(
      "CoherentPiUnc_Fe");
  mnv_plotter.error_summary_group_map["Diffractive"].push_back(
      "CoherentPiUnc_H2O");
  mnv_plotter.error_summary_group_map["Diffractive"].push_back(
      "CoherentPiUnc_Pb");
  // for(auto g : systematics::kGenieSystematics_FSI)
  //  mnv_plotter.error_summary_group_map["Genie_FSI"].push_back(g);

  for (auto g : systematics::kGenieSystematics_FSI_nucleons)
    mnv_plotter.error_summary_group_map["Genie_FSI_nucleons"].push_back(g);

  for (auto g : systematics::kGenieSystematics_FSI_pions)
    mnv_plotter.error_summary_group_map["Genie_FSI_pions"].push_back(g);

  for (auto g : systematics::kGenieSystematics_InteractionModel)
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(g);

  mnv_plotter.error_summary_group_map["Detector"].push_back("EmuRangeCurve");
  mnv_plotter.error_summary_group_map["Detector"].push_back("Birks");
  mnv_plotter.error_summary_group_map["Detector"].push_back("BetheBloch");
  mnv_plotter.error_summary_group_map["Detector"].push_back("Mass");
  mnv_plotter.error_summary_group_map["Detector"].push_back("PartResp");
  mnv_plotter.error_summary_group_map["Detector"].push_back("TrackAngle");
  mnv_plotter.error_summary_group_map["Detector"].push_back("BeamAngle");
  mnv_plotter.error_summary_group_map["Detector"].push_back("NodeCutEff");
  mnv_plotter.error_summary_group_map["Detector"].push_back("BeamAngleX");
  mnv_plotter.error_summary_group_map["Detector"].push_back("BeamAngleY");

  mnv_plotter.error_summary_group_map["RPA"].push_back("RPA_LowQ2");
  mnv_plotter.error_summary_group_map["RPA"].push_back("RPA_HighQ2");

  //-- define colors of the standard errors
  mnv_plotter.error_color_map.clear();

  /*
    //Systematic color scheme
    mnv_plotter.error_color_map["Flux"]                   = kYellow-3;
    mnv_plotter.error_color_map["Genie_FSI"]              = kGreen+2;
    mnv_plotter.error_color_map["Genie_InteractionModel"] = kPink+2;
    mnv_plotter.error_color_map["Detector"]               = kCyan+2;
    mnv_plotter.error_color_map["RPA"]                    = kOrange+2;
    mnv_plotter.error_color_map["NonResPi"]               = kRed+2;
    mnv_plotter.error_color_map["Low_Recoil_2p2h_Tune"]   = kViolet+2;
  */
}

void Plot_ErrorGroup2D(EventSelectionPlotInfo2D p, PlotUtils::MnvH2D* h,
                     std::string error_group_name, std::string tag,
                     double ignore_threshold = 0., double ymax = -1.) {
  PlotUtils::GridCanvas* gridCanvas =new PlotUtils::GridCanvas(Form("2D_ErrorSummary_%s_%s_vs_%s_%s", tag.c_str(), p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str(), error_group_name.c_str()), 4, p.m_variable2D->NBinsY()/4 + 1, 1400, 600);

  p.m_mnv_plotter.good_colors = MnvColors::GetColors(MnvColors::k36Palette);

  // Clone hist
  PlotUtils::MnvH2D* hist = (PlotUtils::MnvH2D*)h->Clone("hist");

  // X label
  p.SetXLabel(hist);

  // For groups MnvPlotter sets the ymax to headroom*axis_maximum_group
  // Or for summary, ymax = axis_maximum
  if (ymax != -1.) {
    p.m_mnv_plotter.axis_maximum = ymax;
    p.m_mnv_plotter.axis_maximum_group = ymax;
    p.m_mnv_plotter.headroom = 1.;
  }
 
  const char* legend_position = error_group_name == "" ? "N" : "TR";
  std::cout << " antes de entrar a DrawErrorSummary2D \n";
  std::vector<TH2D*> ErrorHists = DrawErrorSummary2D(p,
        hist, p.m_include_stat, true, ignore_threshold,
        p.m_do_cov_area_norm, error_group_name, p.m_do_frac_unc);
  std::cout << " despues de entrar a DrawErrorSummary2D \n";
/*  double mx[11];
  std::cout << "Antes de GetMultipliers\n";
  std::vector<double> mult = GetMultipliers(ErrorHists[0]);
  std::cout << "Despues de GetMultipliers\n";  
  std::cout << "Antes del ciclo de multipliers\n";
  for (int i = 0; i < 11; ++i){
    if (mult[i] > 1000) mx[i] = 1.;
    else mx[i] = mult[i];
  }*/  
 // std::cout << "Despues del ciclo de multipliers\n";
  std::vector<double> mult = GetMultipliers(ErrorHists[1]);
  double mx[p.m_variable2D->NBinsY()];
  std::vector<int> bins;
  for (int i = 0; i < p.m_variable2D->NBinsY(); ++i){
    bins.push_back(i+1);
    if (mult[i] > 1000.) mx[i] = 1.;
    else mx[i] = mult[i];
  }
  bins.push_back(p.m_variable2D->NBinsY()+1);

  auto legend = new TLegend(0.80,0.20,0.98,0.33);
  gridCanvas->SetRightMargin(0.1);
  gridCanvas->SetLeftMargin(0.1);
  gridCanvas->ResetPads();
  std::vector<std::string> names;
  names = error_names[error_group_name];
  for (int i = 0; i < (int)names.size(); ++i)
    std::cout << "Group name " << names[i] << "\n";
  gridCanvas->SetRightMargin(0.01);
  gridCanvas->SetLeftMargin(0.1);
  gridCanvas->ResetPads();
  std::cout << "Despues de modificar grdCanvas\n";
//  std::vector<std::string> names;
  std::cout << "Antes del primer DrawOneHist\n";
  gridCanvas->DrawOneHist( ErrorHists[0] , "HIST", false, bins, false, mx);
  std::cout << "Despues del primer DrawOneHist";
  legend->SetNColumns(2);
  if (error_group_name == ""){
    names = error_names["Totals"];
    legend->AddEntry(ErrorHists[0], names[0].c_str(), "l");
  }
  else {
    names = error_names[error_group_name];
    legend->AddEntry(ErrorHists[0], names[0].c_str(), "l");
  }
  // Draw
  for (int i = 1; i < (int)ErrorHists.size(); ++i) {
    if (error_group_name == ""){
      legend->AddEntry(ErrorHists[i], names[i].c_str(), "l");
    }
    else 
      legend->AddEntry(ErrorHists[0], names[i].c_str(), "l");     
    gridCanvas->DrawOneHist( ErrorHists[i], "SAME HIST", false, bins, false, mx);
  }
  gridCanvas->DrawBinRanges(ErrorHists[0], 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  gridCanvas->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  gridCanvas->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("%s %s vs %s %s", tag.c_str(), p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str(), 
             GetSignalName(p.m_signal_definition).c_str()));
  std::string yaxis = "Fractional Uncertainty ";
  gridCanvas->SetYTitle(yaxis.c_str());
  legend->Draw();
  gridCanvas->ResetPads();
  gridCanvas->Draw();


  gridCanvas->Print(Form("2D_ErrorSummary_%s_%s_vs_%s_%s.png", tag.c_str(),
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   error_group_name.c_str()));
}

//-------------------------------------
//------Event Selection Plots----------
//-------------------------------------

void PlotVar_Selection2D(EventSelectionPlotInfo2D p, bool do_bg = true,
                         bool do_log_scale = false, bool do_tuned_bg = false,
			 bool do_bin_width_norm = true){

  std::cout << "Plotting Selection 2D" << p.m_variable2D->NameX() << "_vs_" << p.m_variable2D->NameY() << std::endl;
    PlotUtils::GridCanvas* StackedPlot =new PlotUtils::GridCanvas(Form("2D_Sel_%s_vs_%s",p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), 4, p.m_variable2D->NBinsY()/4 + 1, 1400, 600);
  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable2D->m_hists2D.m_selection_data);
  assert(p.m_variable2D->m_hists2D.m_selection_mc.hist);

  // Get Hists
  TH2D* data = nullptr;
  if (!p.m_variable2D->m_is_true)
    data = (TH2D*)p.m_variable2D->m_hists2D.m_selection_data->Clone(
        "data");
  TH2D* mc =
      (TH2D*)p.m_variable2D->m_hists2D.m_selection_mc.hist->Clone(
          "mc");
  TH2D* mc_error =
      new TH2D(p.m_variable2D->m_hists2D.m_selection_mc.hist->GetCVHistoWithStatError());

  // Background
  TH2D* tmp_bg = nullptr;
  if (do_bg) {
    if (do_tuned_bg)
      tmp_bg = (TH2D*)(p.m_variable2D->m_hists2D.m_tuned_bg)
                   ->Clone("bg_tmp");
    else
      tmp_bg = (TH2D*)(p.m_variable2D->m_hists2D.m_bg.hist)
                   ->Clone("bg_tmp");
  } else
    tmp_bg = NULL;

  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.8);
  data->SetMarkerColor(1);
  data->SetLineWidth(1);
  data->SetLineStyle(1);
  data->SetLineColor(1);
  tmp_bg->SetFillColor(46);
  tmp_bg->SetFillStyle(3005);
  tmp_bg->SetLineColor(46);
  tmp_bg->SetLineWidth(1);

  mc->SetFillColor(0);
  mc->SetLineColor(2);
  mc->SetLineStyle(1);
  mc->SetLineWidth(3);

  mc_error->SetFillColor(kRed-10);
  mc_error->SetFillStyle(1001);
  mc_error->SetMarkerStyle(0);

  // Log Scale
/*  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }*/

  // Y-axis limit
//  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  // p.SetXLabel(p.m_variable->m_hists2D.m_selection_mc.hist);
  p.SetXLabel(mc_error);
  p.SetXLabel(mc);
  // p.SetXLabel(data);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm)
    pot_scale = 1.;
  else
    pot_scale = p.m_data_pot / p.m_mc_pot;

  // Bin Width Normalization
  if (do_bin_width_norm) {
    if (tmp_bg) tmp_bg->Scale(pot_scale, "width");
    if (data) data->Scale(1., "width");
    mc->Scale(pot_scale, "width");
    mc_error->Scale(pot_scale, "width");
    // Y label
//  std::string yaxis = "N Events / " + p.m_variable2D->m_unitsY;
//  mc->GetYaxis()->SetTitle(yaxis.c_str());
  }
  double mx[p.m_variable2D->NBinsY()];
  std::vector<double> mult = GetMultipliers(data);
  std::vector<int> bins;
  for (int i = 0; i < p.m_variable2D->NBinsY(); ++i){
    bins.push_back(i+1);
    if (mult[i] > 1000) mx[i] = 1.;
    else mx[i] = mult[i];
  }  
  bins.push_back(p.m_variable2D->NBinsY()+1);

  auto legend = new TLegend(0.80,0.20,0.98,0.33);
  legend->AddEntry(mc, "MC", "l");
  legend->AddEntry(tmp_bg, "BG", "f");
  legend->AddEntry(data, "Data", "p");

  StackedPlot->SetRightMargin(0.01);
  StackedPlot->SetLeftMargin(0.1);
  StackedPlot->ResetPads();

  // Draw
  StackedPlot->DrawOneHist( mc_error, "E2", false, bins, false, mx);
  StackedPlot->DrawOneHist( mc, "SAME HIST", false, bins, false, mx);
  StackedPlot->DrawOneHist( tmp_bg, "SAME HIST", false, bins, false, mx);
  StackedPlot->DrawOneHist(data, "SAME E1", false, bins, false, mx);
  StackedPlot->DrawBinRanges(data, 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  StackedPlot->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  StackedPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("Selection %s vs %s %s", p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str(), 
             GetSignalName(p.m_signal_definition).c_str()));
  std::string yaxis = "N Events / " + p.m_variable2D->m_unitsY;
  StackedPlot->SetYTitle(yaxis.c_str());
  StackedPlot->ResetPads();
  legend->Draw();
  StackedPlot->Draw();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bg_str = do_tuned_bg ? "_tunedBG" : "_untunedBG";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  StackedPlot->Print(Form("2D_Sel_%s_vs_%s_%s_%s.png",
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(), bg_str.c_str(),
 			   bwn_str.c_str()));

}

void PlotVar_ErrorSummary2D(EventSelectionPlotInfo2D p) {
  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable2D->m_hists2D.m_selection_mc.hist);

  SetErrorGroups2D(p.m_mnv_plotter);

  PlotUtils::MnvH2D* sel =
      (PlotUtils::MnvH2D*)p.m_variable2D->m_hists2D.m_selection_mc.hist->Clone(
          uniq());
  Plot_ErrorGroup2D(p, sel, "", "Sel", 0.0, 0.35);
//  Plot_ErrorGroup(p, sel, "LEGENDONLY", "Sel", 0.0, 0.2);
//  Plot_ErrorGroup2D(p, sel, "2p2h", "Sel", 0.0, 0.01);
//  Plot_ErrorGroup2D(p, sel, "Detector", "Sel", 0.0, 0.15);
//  Plot_ErrorGroup2D(p, sel, "Flux", "Sel", 0.0, 0.15);
//  Plot_ErrorGroup2D(p, sel, "Genie_FSI_nucleons", "Sel", 0.004, 0.06);
//  Plot_ErrorGroup2D(p, sel, "Genie_FSI_pions", "Sel", 0.01, 0.1);
//  Plot_ErrorGroup2D(p, sel, "Genie_InteractionModel", "Sel", 0.01, 0.2);
//  Plot_ErrorGroup2D(p, sel, "Muon", "Sel", 0.0, 0.14);
//  Plot_ErrorGroup2D(p, sel, "NonResPi", "Sel", 0.0, 0.06);
//  Plot_ErrorGroup2D(p, sel, "RPA", "Sel", 0.0, 0.012);
  //  Plot_ErrorGroup(p, sel, "Michel", "Sel", 0.0, 0.15);
//  Plot_ErrorGroup2D(p, sel, "GENIE", "Sel", 0.0, 0.30);
//  Plot_ErrorGroup2D(p, sel, "Target", "Sel", 0.0, 0.15);
//  Plot_ErrorGroup2D(p, sel, "Response", "Sel", 0.0, 0.05);
//  Plot_ErrorGroup2D(p, sel, "Diffractive", "Sel", 0.0, 0.15);
//  Plot_ErrorGroup2D(p, sel, "PhysicsModel", "Sel", 0.0, 0.15);
}




//==============================================================================
// Background-Subtracted
//==============================================================================
void Plot_BGSub2D(EventSelectionPlotInfo2D p, std::string outdir = ".",
                double ymax = -1, bool do_log_scale = false,
                bool do_bin_width_norm = true) {
  std::cout << "Plotting BG-subtracted Data " << p.m_variable2D->NameX()
	    << " vs " << p.m_variable2D->NameY()
            << std::endl;

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable2D->m_hists2D.m_bg_subbed_data);
  assert(p.m_variable2D->m_hists2D.m_effnum.hist);
  
  PlotUtils::GridCanvas* StackedPlot =new PlotUtils::GridCanvas(Form("2D_BG_%s_vs_%s",p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), 4, p.m_variable2D->NBinsY()/4 + 1, 1400, 600);
  StackedPlot->SetTitle(Form("Background Subtracted %s", GetSignalName(p.m_signal_definition).c_str()));


  // Get Hists
  TH2D* bg_sub_data_w_tot_error =
      new TH2D(p.m_variable2D->m_hists2D.m_bg_subbed_data->GetCVHistoWithError());
  TH2D* bg_sub_data_w_stat_error = new TH2D(
      p.m_variable2D->m_hists2D.m_bg_subbed_data->GetCVHistoWithStatError());
  TH2D* effnum_cv =
      new TH2D(p.m_variable2D->m_hists2D.m_effnum.hist->GetCVHistoWithStatError());
  TH2D* effnum_w_stat_error =
      new TH2D(p.m_variable2D->m_hists2D.m_effnum.hist->GetCVHistoWithStatError());

  bg_sub_data_w_tot_error->SetMarkerStyle(20);
  bg_sub_data_w_tot_error->SetMarkerSize(0.8);
  bg_sub_data_w_tot_error->SetMarkerColor(1);
  bg_sub_data_w_tot_error->SetLineWidth(1);
  bg_sub_data_w_tot_error->SetLineStyle(1);
  bg_sub_data_w_tot_error->SetLineColor(1);
  bg_sub_data_w_stat_error->SetMarkerStyle(20);
  bg_sub_data_w_stat_error->SetMarkerSize(0.8);
  bg_sub_data_w_stat_error->SetMarkerColor(1);
  bg_sub_data_w_stat_error->SetLineWidth(1);
  bg_sub_data_w_stat_error->SetLineStyle(1);
  bg_sub_data_w_stat_error->SetLineColor(1);

  effnum_cv->SetFillColor(0);
  effnum_cv->SetLineColor(2);
  effnum_cv->SetLineStyle(1);
  effnum_cv->SetLineWidth(3);

  effnum_w_stat_error->SetFillColor(kRed-10);
  effnum_w_stat_error->SetFillStyle(1001);
  effnum_w_stat_error->SetMarkerStyle(0);

  // Log Scale
/*  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
*/
  // Y-axis range
//if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  p.SetXLabel(bg_sub_data_w_tot_error);
  p.SetXLabel(bg_sub_data_w_stat_error);
  p.SetXLabel(effnum_w_stat_error);

  // Overall Normalization
  double pot_scale = -99.;
  // PlotUtils::MnvH1D* tmp_bg = nullptr;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  std::string yaxis = "N Events";  
  if (do_bin_width_norm) {
    bg_sub_data_w_tot_error->Scale(1., "width");
    bg_sub_data_w_stat_error->Scale(1., "width");
    effnum_w_stat_error->Scale(pot_scale, "width");
    effnum_cv->Scale(pot_scale, "width");
    // Y label
    yaxis = yaxis + " / " + p.m_variable2D->m_unitsY;
//  effnum_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }
  // Draw
//  const bool use_hist_titles = false;
//  p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_tot_error,
//                                          effnum_w_stat_error, pot_scale, "TR",
//                                          use_hist_titles);

  // p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_stat_error,
  // effnum_w_stat_error, pot_scale, "TR");
  // p.m_mnv_plotter.DrawDataMC(bg_sub_data_w_tot_error, effnum_w_stat_error,
  // pot_scale, "TR");
  // p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_error, effnum,
  // pot_scale, "TR",
  //                                        use_hist_titles, NULL, NULL,
  //                                        p.m_do_cov_area_norm,
  //                                        p.m_include_stat);

  // POT info
//if (p.m_do_cov_area_norm)
//  p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);
//else
//  p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);

  // Label re: error bars
//p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.7, 0.75, 0.03);

  // Plot Title
//  p.m_mnv_plotter.title_size = 0.04;
//  StackedPlot->SetTitleSize(10.);
  double mx[p.m_variable2D->NBinsY()];
  std::vector<double> mult = GetMultipliers(bg_sub_data_w_stat_error);
  std::vector<int> bins;
  for (int i = 0; i < p.m_variable2D->NBinsY(); ++i){
    bins.push_back(i+1);
    if (mult[i] > 1000) mx[i] = 1.;
    else mx[i] = mult[i];
  }
  bins.push_back(p.m_variable2D->NBinsY()+1);

//  std::vector<int> bins = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

  auto legend = new TLegend(0.80,0.20,0.98,0.33);
  legend->AddEntry(effnum_cv, "MC", "l");
//legend->AddEntry(tmp_bg, "BG", "f");
  legend->AddEntry(bg_sub_data_w_stat_error, "Data", "p");

  StackedPlot->SetRightMargin(0.01);
  StackedPlot->SetLeftMargin(0.1);
  StackedPlot->ResetPads();

  // Draw
  StackedPlot->DrawOneHist( effnum_w_stat_error, "E2", false, bins, false, mx);
  StackedPlot->DrawOneHist( effnum_cv, "SAME HIST", false, bins, false, mx);
  StackedPlot->DrawOneHist(bg_sub_data_w_stat_error, "SAME E1", false, bins, false, mx);
  StackedPlot->DrawBinRanges(bg_sub_data_w_stat_error, 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  StackedPlot->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  StackedPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("Background Subtracted %s vs %s %s", p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str(),
             GetSignalName(p.m_signal_definition).c_str()));
  StackedPlot->SetYTitle(yaxis.c_str());
  StackedPlot->SetTitle(Form("Background Subtracted %s", GetSignalName(p.m_signal_definition).c_str()));
  StackedPlot->ResetPads();
  legend->Draw();
  StackedPlot->SetHistTexts();
  StackedPlot->Draw();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  StackedPlot->Print(Form("2D_BG_%s_vs_%s_%s.png",
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));

}

void PlotMigration2D(EventSelectionPlotInfo2D p, PlotUtils::MnvH2D* hist, std::string nameX, std::string nameY) {
  TGaxis::SetExponentOffset(-0.035, -0.048, "x");
  TH2D* htmp = GetHistWithUnderOverFlow(hist);
  TH2D* htmp2 = RowNormalize(htmp);
  htmp2->GetXaxis()->SetTitle(Form("Reconstructed %s per %s bin", p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_hists2D.m_xlabelY.c_str()));
  htmp2->GetYaxis()->SetTitle(Form("Truth %s per %s bin", p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_hists2D.m_xlabelY.c_str()));
  TCanvas c("c1", "c1");
  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  mnv_plotter.SetWhiteRainbowPalette();
//  mnv_plotter.SetRedHeatPalette();
  bool draw_as_matrix = false;
  bool coarseContours = false;
  bool includeFlows = true;
  bool no_text = true;
  gStyle->SetHistMinimumZero(kFALSE);
  mnv_plotter.DrawNormalizedMigrationHistogram(htmp2, draw_as_matrix, coarseContours, 
                                               includeFlows, no_text);
  // gStyle->SetPaintTextFormat("2.0f");
  // htmp2->SetMarkerSize(2);
  // htmp2->Draw("colz text");
  p.SetTitle("Migration matrix");
  c.Update();
  c.Print(Form("2D_Migration_VarBins_%s_vs_%s.png", nameX.c_str(), nameY.c_str()));
  // c.SetLogz();
  // c.Update();
  // c.Print("WMigrationMatrix_Wbins_logz.png");
  TGaxis::SetExponentOffset(0, 0, "x");
}

void Plot_Unfolded2D(EventSelectionPlotInfo2D p, MnvH2D* data, MnvH2D* mc,
                   std::string outdir = ".", double ymax = -1,
                   bool do_log_scale = false, bool do_bin_width_norm = true) {
  std::cout << "Plotting Unfolded " << p.m_variable2D->NameX() << " vs " << p.m_variable2D->NameY()<< std::endl;
  PlotUtils::GridCanvas* GridPlot =new PlotUtils::GridCanvas(Form("2D_Unfold_%s_vs_%s",p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), 4, p.m_variable2D->NBinsY()/4 + 1, 1400, 600);
  GridPlot->SetTitle(Form("Unfolding %s", GetSignalName(p.m_signal_definition).c_str()));

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

//  TCanvas canvas("c1", "c1");

  // Get Hists
  PlotUtils::MnvH2D* unfolded = (PlotUtils::MnvH2D*)data->Clone("unfolded");
  PlotUtils::MnvH2D* effnum_true = (PlotUtils::MnvH2D*)mc->Clone("effnum_true");
  TH2D* unfolded_w_tot_error = new TH2D(unfolded->GetCVHistoWithError());
  TH2D* unfolded_w_stat_error = new TH2D(unfolded->GetCVHistoWithStatError());
  TH2D* effnum_cv =
      new TH2D(effnum_true->GetCVHistoWithStatError());
  TH2D* effnum_w_stat_error =
      new TH2D(effnum_true->GetCVHistoWithStatError());

  unfolded_w_tot_error->SetMarkerStyle(20);
  unfolded_w_tot_error->SetMarkerSize(0.8);
  unfolded_w_tot_error->SetMarkerColor(1);
  unfolded_w_tot_error->SetLineWidth(1);
  unfolded_w_tot_error->SetLineStyle(1);
  unfolded_w_tot_error->SetLineColor(1);
  unfolded_w_stat_error->SetMarkerStyle(20);
  unfolded_w_stat_error->SetMarkerSize(0.8);
  unfolded_w_stat_error->SetMarkerColor(1);
  unfolded_w_stat_error->SetLineWidth(1);
  unfolded_w_stat_error->SetLineStyle(1);
  unfolded_w_stat_error->SetLineColor(1);
  effnum_cv->SetFillColor(0);
  effnum_cv->SetLineColor(2);
  effnum_cv->SetLineStyle(1);
  effnum_cv->SetLineWidth(3);

  effnum_w_stat_error->SetFillColor(kRed-10);
  effnum_w_stat_error->SetFillStyle(1001);
  effnum_w_stat_error->SetMarkerStyle(0);


  // Log Scale
/*  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
*/
  // Y-axis range
//  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  p.SetXLabel(effnum_w_stat_error);
  p.SetXLabel(unfolded_w_tot_error);
  p.SetXLabel(unfolded_w_stat_error);
  p.SetXLabel(effnum_cv);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  std::string yaxis = "N Events";
  if (do_bin_width_norm) {
    unfolded_w_tot_error->Scale(1., "width");
    unfolded_w_stat_error->Scale(1., "width");
    effnum_w_stat_error->Scale(pot_scale, "width");
    effnum_cv->Scale(pot_scale, "width");

    // Y label
    yaxis = yaxis + " / " + p.m_variable2D->m_unitsY;
  }

  // Draw
  double mx[p.m_variable2D->NBinsY()];
  std::vector<double> mult = GetMultipliers(unfolded_w_stat_error);
  std::vector<int> bins;
  for (int i = 0; i < p.m_variable2D->NBinsY(); ++i){
    bins.push_back(i+1);
    if (mult[i] > 1000) mx[i] = 1.;
    else mx[i] = mult[i];
  }
  bins.push_back(p.m_variable2D->NBinsY()+1);

//  std::vector<int> bins = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

  auto legend = new TLegend(0.80,0.20,0.98,0.33);
  legend->AddEntry(effnum_cv, "MC", "l");
//legend->AddEntry(tmp_bg, "BG", "f");
  legend->AddEntry(unfolded_w_stat_error, "Data", "p");

  GridPlot->SetRightMargin(0.01);
  GridPlot->SetLeftMargin(0.1);
  GridPlot->ResetPads();

  // Draw
  GridPlot->DrawOneHist( effnum_w_stat_error, "E2", false, bins, false, mx);
  GridPlot->DrawOneHist( effnum_cv, "SAME HIST", false, bins, false, mx);
  GridPlot->DrawOneHist(unfolded_w_stat_error, "SAME E1", false, bins, false, mx);
  GridPlot->DrawBinRanges(unfolded_w_stat_error, 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  GridPlot->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  GridPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("Unfolded %s vs %s %s", p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str(),
             GetSignalName(p.m_signal_definition).c_str()));
  GridPlot->SetYTitle(yaxis.c_str());
//GridPlot->SetTitle(Form("Background Subtracted %s", GetSignalName(p.m_signal_definition).c_str()));
  GridPlot->ResetPads();
  legend->Draw();
  GridPlot->Draw();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  GridPlot->Print(Form("2D_Unfolded_%s_vs_%s_%s.png",
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));


}
void PlotMC2D(EventSelectionPlotInfo2D p, MnvH2D* mc, std::string ylabel = "",
                   std::string outdir = ".", double ymax = -1,
                   bool do_log_scale = false, bool do_bin_width_norm = true) {
  std::cout << ylabel + " 2D " << p.m_variable2D->NameX() << " vs " << p.m_variable2D->NameY()<< std::endl;
  PlotUtils::GridCanvas* GridPlot =new PlotUtils::GridCanvas(Form("2D_%s_%s_vs_%s", ylabel.c_str(), p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), 4, p.m_variable2D->NBinsY()/4 + 1, 1400, 600);
  GridPlot->SetTitle(Form("%s %s", ylabel.c_str(), GetSignalName(p.m_signal_definition).c_str()));
  // Make sure we remembered to load the source histos from the input file.
  assert(mc);
//  TCanvas canvas("c1", "c1");
  // Get Hists
  PlotUtils::MnvH2D* hist = (PlotUtils::MnvH2D*)mc->Clone("hist");
  TH2D* hist_w_cv =
      new TH2D(hist->GetCVHistoWithError());
  TH2D* hist_w_tot_error =
      new TH2D(hist->GetCVHistoWithError());

  hist_w_tot_error->SetFillColor(kRed-10);
  hist_w_tot_error->SetFillStyle(1001);
  hist_w_tot_error->SetMarkerStyle(0);
 
  hist_w_cv->SetFillColor(0);
  hist_w_cv->SetLineColor(2);
  hist_w_cv->SetLineStyle(1);
  hist_w_cv->SetLineWidth(3);


  // Log Scale
/*  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
*/
  // Y-axis range
//  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  p.SetXLabel(hist_w_tot_error);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  std::string yaxis = ylabel;
  if (do_bin_width_norm) {
    hist_w_tot_error->Scale(pot_scale, "width");
    hist_w_cv->Scale(pot_scale, "width");
    // Y label
    yaxis = yaxis + " / " + p.m_variable2D->m_unitsY;
  }

  // Draw
//  std::cout << "Number of bins " << p.m_variable2D->NBinsY() <<"\n";
  double mx[p.m_variable2D->NBinsY()];
  std::vector<double> mult = GetMultipliers(hist_w_cv);
  std::vector<int> bins;
  for (int i = 0; i < p.m_variable2D->NBinsY(); ++i){
    bins.push_back(i+1);
//    std::cout << "Pasa 1\n";
    if (mult[i] > 1000.) mx[i] = 1.;
    else mx[i] = mult[i];
  }
  bins.push_back(p.m_variable2D->NBinsY()+1);
//  std::cout << "Pasa 2 " << sizeof(mx)/sizeof(double) <<"\n";
  auto legend = new TLegend(0.80,0.20,0.98,0.33);
  legend->AddEntry(hist_w_cv, "MC", "l");
  
  GridPlot->SetRightMargin(0.01);
  GridPlot->SetLeftMargin(0.1);
  GridPlot->ResetPads();

  // Draw
  GridPlot->DrawOneHist( hist_w_tot_error, "E2", false, bins, false, mx);
  GridPlot->DrawOneHist( hist_w_cv, "SAME HIST", false, bins, false, mx);
  GridPlot->DrawBinRanges(hist_w_cv, 2, bins, Form("%s (%s)",
                          p.m_variable2D->m_hists2D.m_xlabelY.c_str(),
                          p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  GridPlot->DrawMultipliers( sizeof(mx)/sizeof(double) , mx);
  GridPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
                      p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(ylabel + " " + p.m_variable2D->m_hists2D.m_xlabelX.c_str() 
             + " vs " + p.m_variable2D->m_hists2D.m_xlabelY.c_str() + " " 
             + GetSignalName(p.m_signal_definition));
  GridPlot->SetYTitle(yaxis.c_str());
  GridPlot->SetTitle(Form("%s %s", ylabel.c_str(), GetSignalName(p.m_signal_definition).c_str()));
//  GridPlot->ResetPads();
  legend->Draw();
  GridPlot->Draw();
  GridPlot->Update();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  GridPlot->Print(Form("2D_%s_%s_vs_%s_%s.png", ylabel.c_str(),
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));
}

void Plot_CrossSection2D(EventSelectionPlotInfo2D p, MnvH2D* data, MnvH2D* mc,
                       std::string outdir = ".", double ymax = -1,
                       bool do_log_scale = false,
                       bool do_bin_width_norm = true) {

  std::cout << "Plotting Cross Section " << p.m_variable2D->NameX() << " vs " << p.m_variable2D->NameY()<< std::endl;
  PlotUtils::GridCanvas* GridPlot =new PlotUtils::GridCanvas(Form("2D_CrossSection_%s_vs_%s",p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), 4, p.m_variable2D->NBinsY()/4 + 1, 1400, 600);
  GridPlot->SetTitle(Form("Cross Section %s", GetSignalName(p.m_signal_definition).c_str()));

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

//  TCanvas canvas("c1", "c1");

  // Get Hists
  PlotUtils::MnvH2D* xsec_data = (PlotUtils::MnvH2D*)data->Clone("xsec_data");
  PlotUtils::MnvH2D* xsec_mc = (PlotUtils::MnvH2D*)mc->Clone("xsec_mc");
  TH2D* xsec_data_w_tot_error = new TH2D(xsec_data->GetCVHistoWithError());
  TH2D* xsec_data_w_stat_error = new TH2D(xsec_data->GetCVHistoWithStatError());
  TH2D* xsec_mc_cv =
      new TH2D(xsec_mc->GetCVHistoWithStatError());
  TH2D* xsec_mc_w_stat_error =
      new TH2D(xsec_mc->GetCVHistoWithStatError());

  xsec_data_w_tot_error->SetMarkerStyle(20);
  xsec_data_w_tot_error->SetMarkerSize(0.8);
  xsec_data_w_tot_error->SetMarkerColor(1);
  xsec_data_w_tot_error->SetLineWidth(1);
  xsec_data_w_tot_error->SetLineStyle(1);
  xsec_data_w_tot_error->SetLineColor(1);
  xsec_data_w_stat_error->SetMarkerStyle(20);
  xsec_data_w_stat_error->SetMarkerSize(0.8);
  xsec_data_w_stat_error->SetMarkerColor(1);
  xsec_data_w_stat_error->SetLineWidth(1);
  xsec_data_w_stat_error->SetLineStyle(1);
  xsec_data_w_stat_error->SetLineColor(1);
  xsec_mc_cv->SetFillColor(0);
  xsec_mc_cv->SetLineColor(2);
  xsec_mc_cv->SetLineStyle(1);
  xsec_mc_cv->SetLineWidth(3);

  xsec_mc_w_stat_error->SetFillColor(kRed-10);
  xsec_mc_w_stat_error->SetFillStyle(1001);
  xsec_mc_w_stat_error->SetMarkerStyle(0);


  // Log Scale
/*  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
*/
  // Y-axis range
//  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  p.SetXLabel(xsec_mc_w_stat_error);
  p.SetXLabel(xsec_data_w_stat_error);
  p.SetXLabel(xsec_mc_cv);
  p.SetXLabel(xsec_mc_w_stat_error);
  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  std::string yaxis = "N Events";
  if (do_bin_width_norm) {
    xsec_data_w_tot_error->Scale(1., "width");
    xsec_data_w_stat_error->Scale(1., "width");
    xsec_mc_w_stat_error->Scale(1., "width");
    xsec_mc_cv->Scale(1., "width");

    // Y label
    yaxis = yaxis + " / " + p.m_variable2D->m_unitsY;
  }

  // Draw
  double mx[p.m_variable2D->NBinsY()];
  std::vector<double> mult = GetMultipliers(xsec_data_w_stat_error);
  std::vector<int> bins;
  for (int i = 0; i < p.m_variable2D->NBinsY(); ++i){
    bins.push_back(i+1);
    if (mult[i] > 1000) mx[i] = 1.;
    else mx[i] = mult[i];
  }
  bins.push_back(p.m_variable2D->NBinsY()+1);

 // std::vector<int> bins = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

  auto legend = new TLegend(0.80,0.20,0.98,0.33);
  legend->AddEntry(xsec_mc_cv, "MC", "l");
//legend->AddEntry(tmp_bg, "BG", "f");
  legend->AddEntry(xsec_data_w_stat_error, "Data", "p");

  GridPlot->SetRightMargin(0.01);
  GridPlot->SetLeftMargin(0.1);
  GridPlot->ResetPads();

  // Draw
  GridPlot->DrawOneHist( xsec_mc_w_stat_error, "SAME E2", false, bins, false, mx);
  GridPlot->DrawOneHist( xsec_mc_cv, "SAME HIST", false, bins, false, mx);
  GridPlot->DrawOneHist(xsec_data_w_tot_error, "SAME E1", false, bins, false, mx);
  GridPlot->DrawBinRanges(xsec_data_w_tot_error, 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  GridPlot->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  GridPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("Cross Section %s vs %s %s", p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str(),
             GetSignalName(p.m_signal_definition).c_str()));
  GridPlot->SetYTitle(yaxis.c_str());
//GridPlot->SetTitle(Form("Background Subtracted %s", GetSignalName(p.m_signal_definition).c_str()));
  GridPlot->ResetPads();
  legend->Draw();
  GridPlot->Draw();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  GridPlot->Print(Form("2D_CrossSection_%s_vs_%s_%s.png",
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));

}


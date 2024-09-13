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
#include "plotting_functions.h"

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
void SetErrorGroups2D(MnvPlotter& mnv_plotter);
std::vector<double> GetMultipliers(TH2D* data_hist){
//  std::cout << "Inside GetMultipliers\n";
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
std::vector<double> GetMultipliers(std::vector <TH2D*> hvec){
//  std::cout << "Inside GetMultipliers(vec)\n";
  std::vector<double> multipliers;
  double Max = 0;
  for (int k = 0; k < (int)hvec.size(); ++k){
    for (int j = 1; j <= (int)hvec[0]->GetNbinsY(); ++j){ 
      for (int i = 1; i <= (int)hvec[0]->GetNbinsX(); ++i){
        if (Max < hvec[k]->GetBinContent(i, j)) Max = hvec[k]->GetBinContent(i, j);
      }
    }
  }
//  std::cout << "Max " << Max << "\n";

  double currentvalue = 0;
//  std::cout << "hvec size = " << hvec.size() << "\n";
  for (int j = 1; j <= hvec[0]->GetNbinsY(); ++j){
    currentvalue = 0;
    for (int k = 0; k < (int)hvec.size(); ++k){
      for(int i = 1; i <= hvec[k]->GetNbinsX(); ++i ){
//        std::cout << "Instance i = " << i << " j = "  << j << 
//		" k = " << k << "\n" ;
        if (hvec[k]->GetBinContent(i,j) > currentvalue)
          currentvalue = hvec[k]->GetBinContent(i,j);  
//        std::cout << "Current value after = " << currentvalue << "\n";
      }
    }
    if (currentvalue == 0) multipliers.push_back(1.);
    else multipliers.push_back((double)nearbyint(Max/currentvalue));       
//    std::cout << "Multiplier " << j-1 << " = " << multipliers[j-1] << "\n";
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
void Plot2D(EventSelectionPlotInfo2D p, std::vector<TH2D*> hvect,
                   std::vector<std::string> names, std::string label = "", 
                   std::string ylabel = "", std::string outdir = ".", double ymax = -1,
                   bool do_log_scale = false, bool do_bin_width_norm = true) {
  std::cout << label + " 2D " << p.m_variable2D->NameX() << " vs " << p.m_variable2D->NameY()<< std::endl;
  int nxpads = 4;
  int nypads = p.m_variable2D->NBinsY()/3;
  if (p.m_variable2D->NBinsY() < 4){
      nxpads = p.m_variable2D->NBinsY();
      nypads = 1;
  }
  if (p.m_variable2D->NBinsY() > 4 && p.m_variable2D->NBinsY() < 8 )
      nypads = 2;
  PlotUtils::GridCanvas* GridPlot =new PlotUtils::GridCanvas(Form("2D_%s_%s_vs_%s", ylabel.c_str(), p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), nxpads, nypads, 1400, 600);
  GridPlot->SetTitle(Form("%s", ylabel.c_str()));
  // Make sure we remembered to load the source histos from the input file.
//  TCanvas canvas("c1", "c1");
  // Get Hists
  TH2D* hist = (TH2D*)hvect[0]->Clone("hist");
 

  // Y-axis range
//  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  p.SetXLabel(hist);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  std::string yaxis = ylabel;
/*  if (do_bin_width_norm) {
    hist->Scale(1., "width");
    // Y labe
    yaxis = yaxis + " / " + p.m_variable2D->m_unitsY;
  }*/

  // Draw
  //std::cout << "Number of bins " << p.m_variable2D->NBinsY() <<"\n";
  double mx[p.m_variable2D->NBinsY()];
  std::vector<double> mult = GetMultipliers(hvect);
  std::vector<int> bins;
  for (int i = 0; i < p.m_variable2D->NBinsY(); ++i){
    bins.push_back(i+1);
//    std::cout << "Pasa 1\n";
    if (mult[i] > 1000.) mx[i] = 1.;
    else mx[i] = mult[i];
  }
  bins.push_back(p.m_variable2D->NBinsY()+1);
//  std::cout << "Pasa 2 " << sizeof(mx)/sizeof(double) <<"\n";
  auto legend = new TLegend(0.88,0.9,1.,0.5);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->AddEntry(hist, Form("%s", names[0].c_str()), "l");
  
  GridPlot->SetRightMargin(0.12);
  GridPlot->SetLeftMargin(0.07);
  hist->GetXaxis()->SetLabelSize(0.03);
  GridPlot->SetInterpadSpace(0.008);
  //  GridPlot->ResetPads();
  // Draw
  GridPlot->DrawOneHist( hist, "E2", false, bins, false, mx);
  for (int i = 1; i < (int)hvect.size(); ++i){
    legend->AddEntry(hvect[i], Form("%s", names[i].c_str()), "l");
    GridPlot->DrawOneHist( hvect[i], "HIST SAME", false, bins, false, mx);    
  }
  GridPlot->DrawBinRanges(hist, 2, bins, Form("%s (%s)",
                          p.m_variable2D->m_hists2D.m_xlabelY.c_str(),
                          p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  GridPlot->DrawMultipliers( sizeof(mx)/sizeof(double) , mx);
  GridPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
                      p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(label + " " + p.m_variable2D->m_hists2D.m_xlabelX.c_str() 
             + " vs " + p.m_variable2D->m_hists2D.m_xlabelY.c_str());
  GridPlot->SetYTitle(yaxis.c_str());
  GridPlot->SetTitle(Form("%ss", label.c_str()));
//  GridPlot->ResetPads();
  GridPlot->Draw();
  legend->Draw();
  GridPlot->Update();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  GridPlot->Print(Form("2D_ErrorSummary_%s_%s_vs_%s_%s.png", label.c_str(),
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));
}

int AddInQuadrature2D(TH2D* a, const TH2D * b) {

  if( a->GetNbinsX() != b->GetNbinsX() || a->GetNbinsY() != b->GetNbinsY()) {
    Error("AddInQuadrature2D", "Histogram axes not compatible.  Doing nothing." );
    return 1;
  }
    
  int firstBinx = 0;
  int firstBiny = 0;
  int lastBinx  = a->GetNbinsX()+1;
  int lastBiny  = a->GetNbinsY()+1;
  for(int iBiny = firstBiny; iBiny <= lastBiny; ++iBiny){
    for( int iBinx = firstBinx; iBinx <= lastBinx; ++iBinx )
    {
      const double aVal = a->GetBinContent(iBinx,iBiny);
      const double bVal = b->GetBinContent(iBinx,iBiny);
      const double aErr = a->GetBinError(iBinx,iBiny);
      const double bErr = b->GetBinError(iBinx,iBiny);
  
      const double val = sqrt( aVal*aVal + bVal*bVal );
      const double err = ( 1E-8 < fabs(val) ) ? sqrt( pow(aVal*aErr,2) + pow(bVal*bErr,2) ) / val : 0.;
  
      a->SetBinContent(iBinx, iBiny, val);
      a->SetBinError(iBinx, iBiny, err);
    }
  }
  return 1;
}

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
    PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
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
            AddInQuadrature2D( itGroupHist->second, hErr );
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
void SetErrorGroups2D(MnvPlotter& mnv_plotter, bool is_subgroups = false) {
  mnv_plotter.error_summary_group_map.clear();
  if (!is_subgroups){
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
    mnv_plotter.error_summary_group_map["Pion_Reconstruction"].push_back("MichelEfficiency");
    mnv_plotter.error_summary_group_map["Flux"].push_back("Flux");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_CH");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_C");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_Fe");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_H2O");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_Pb");
    mnv_plotter.error_summary_group_map["Others"].push_back("response_em");
    mnv_plotter.error_summary_group_map["Others"].push_back("response_meson");
    mnv_plotter.error_summary_group_map["Others"].push_back("response_other");
    mnv_plotter.error_summary_group_map["Others"].push_back("response_proton");
    mnv_plotter.error_summary_group_map["Others"].push_back("GEANT_Proton");
    mnv_plotter.error_summary_group_map["Others"].push_back("GEANT_Pion");
    mnv_plotter.error_summary_group_map["Others"].push_back("GEANT_Neutron");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("GENIE_D2_NormCCRES");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "DiffractiveModelUnc");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "CoherentPiUnc_C");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "CoherentPiUnc_CH");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "CoherentPiUnc_Fe");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "CoherentPiUnc_H2O");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "CoherentPiUnc_Pb");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("GENIE_Rvn1pi");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("GENIE_Rvp1pi");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("GENIE_Rvn2pi");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("GENIE_Rvp2pi");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "RPA_LowQ2");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "RPA_HighQ2");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "Low_Recoil_2p2h_Tune");
    for (auto g : systematics::kGenieSystematics_InteractionModel)
      mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(g);
    for (auto g : systematics::kGenieSystematics_FSI_nucleons)
      mnv_plotter.error_summary_group_map["GENIE_FSI"].push_back(g);

    for (auto g : systematics::kGenieSystematics_FSI_pions)
      mnv_plotter.error_summary_group_map["GENIE_FSI"].push_back(g);
    mnv_plotter.error_summary_group_map["Pion_Reconstruction"].push_back("TrackAngle");
    mnv_plotter.error_summary_group_map["Muon"].push_back("BeamAngle");
    mnv_plotter.error_summary_group_map["Muon"].push_back("BeamAngleX");
    mnv_plotter.error_summary_group_map["Muon"].push_back("BeamAngleY");
    mnv_plotter.error_summary_group_map["Pion_Reconstruction"].push_back("Birks");
    mnv_plotter.error_summary_group_map["Pion_Reconstruction"].push_back("BetheBloch");
    mnv_plotter.error_summary_group_map["Pion_Reconstruction"].push_back("Mass");
    mnv_plotter.error_summary_group_map["Pion_Reconstruction"].push_back("NodeCutEff");  
  }
  if (is_subgroups){
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("NonRESPi");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("MnvTunes");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("CCQE");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("RESPi");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "Coherent-Diffractive");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back(
        "DIS-Hadronization");
    mnv_plotter.error_summary_group_map["Cross_Section_Models"].push_back("Elastic");

    mnv_plotter.error_summary_group_map["NonRESPi"].push_back("GENIE_Rvn1pi");
    mnv_plotter.error_summary_group_map["NonRESPi"].push_back("GENIE_Rvp1pi");
    mnv_plotter.error_summary_group_map["NonRESPi"].push_back("GENIE_Rvn2pi");
    mnv_plotter.error_summary_group_map["NonRESPi"].push_back("GENIE_Rvp2pi");
    mnv_plotter.error_summary_group_map["MnvTunes"].push_back(
        "Low_Recoil_2p2h_Tune");
    mnv_plotter.error_summary_group_map["MnvTunes"].push_back(
        "RPA_LowQ2");
    mnv_plotter.error_summary_group_map["MnvTunes"].push_back(
        "RPA_HighQ2");
    mnv_plotter.error_summary_group_map["CCQE"].push_back("GENIE_MaCCQE");
    mnv_plotter.error_summary_group_map["CCQE"].push_back("GENIE_VecFFCCQEshape");
    mnv_plotter.error_summary_group_map["CCQE"].push_back("GENIE_CCQEPauliSupViaKF");
    mnv_plotter.error_summary_group_map["RESPi"].push_back("GENIE_D2_MaRES");
    mnv_plotter.error_summary_group_map["RESPi"].push_back("GENIE_EP_MvRES");
    mnv_plotter.error_summary_group_map["RESPi"].push_back("GENIE_NormNCRES");
    mnv_plotter.error_summary_group_map["RESPi"].push_back("GENIE_Theta_Delta2Npi");
    mnv_plotter.error_summary_group_map["RESPi"].push_back("GENIE_D2_NormCCRES");
    mnv_plotter.error_summary_group_map["Coherent-Diffractive"].push_back("DiffractiveModelUnc");
    mnv_plotter.error_summary_group_map["Coherent-Diffractive"].push_back(
        "CoherentPiUnc_CH");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back("GENIE_AhtBY");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back("GENIE_BhtBY");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back("GENIE_CV1uBY");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back("GENIE_CV2uBY");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back("GENIE_NormDISCC");
    mnv_plotter.error_summary_group_map["Elastic"].push_back("GENIE_MaNCEL");
    mnv_plotter.error_summary_group_map["Elastic"].push_back("GENIE_EtaNCEL");



  }
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

TH2D* AddingGroups(MnvPlotter& mnv_plotter, PlotUtils::MnvH2D& h,
                  std::string subgroup) {
  SetErrorGroups2D(mnv_plotter, true);
  std::vector<std::string> group =
      mnv_plotter.error_summary_group_map[subgroup];
  TH2D* ErrGroup = (TH2D*)h.Clone(Form("%s", subgroup.c_str()));
  ErrGroup->Reset();
  for (int i = 0; i < (int)group.size(); ++i) {
    TH2D* hErr = dynamic_cast<TH2D*>(
        h.GetVertErrorBand(group[i])->GetErrorBand(true, false).Clone(uniq()));
    AddInQuadrature2D(ErrGroup, hErr);
//    std::cout << "Error " << group[i] << " " << hErr->GetBinContent(2,2)  << " Errorgroup " << ErrGroup->GetBinContent(2,2) << "\n"; 
  }
//  std::vector <TH2D*> ErrGroupVec;
//std::cout << "ErrorGroup " << subgroup << " " << ErrGroup->GetBinContent(2,2) << "\n";
  return ErrGroup;
//  ErrGroupVec.push_back(ErrGroup);
//  h.AddVertErrorBand( subgroup, ErrGroupVec);
/*  for (int j = 1; j < (int)ErrGroup->GetNbinsY(); ++j){ 
    for(int i = 1; i < (int)ErrGroup->GetNbinsY(); ++i)
      h.GetVertErrorBand(subgroup)->GetErrorBand(false, false).SetBinContent(i, j, ErrGroup->GetBinContent(i, j));
  }*/
}

void NiceColors(EventSelectionPlotInfo2D p, std::vector<TH2D*> hvec, 
		std::vector<std::string> names){

//  p.m_mnv_plotter.good_colors = MnvColors::GetColors(MnvColors::k36Palette);
  if(hvec.size() != names.size()){
    std::cout << "The size of the vectors are not compatible\n";
    return;
  }
  std::string alias = "Cross Section Models"; 

  for (int i = 0; i < (int)hvec.size(); ++i){
//    if (names[i] == "NonRESPi") alias = "Sidebands";
//    else alias = names[i];
//    map<string,int>::const_iterator itErrCol = p.m_mnv_plotter.error_color_map.find(alias);
//    if ( p.m_mnv_plotter.error_color_map.end() != itErrCol )
//      std::cout << "Colors = " << itErrCol->second << "\n";
        hvec[i]->SetLineColor(206+(i*3));
        hvec[i]->SetLineWidth(p.m_mnv_plotter.mc_line_width);
    
  }
}

void Plot_ErrorGroup2D(EventSelectionPlotInfo2D p, PlotUtils::MnvH2D* h,
                     std::string error_group_name, std::string tag,
                     double ignore_threshold = 0., double ymax = -1.) {
  //std::cout << "Error Group name = " << error_group_name << "\n";
  int nxpads = 4;
  int nypads = p.m_variable2D->NBinsY()/3;
  if (p.m_variable2D->NBinsY() < 4){
      nxpads = p.m_variable2D->NBinsY();
      nypads = 1;
  }
  if (p.m_variable2D->NBinsY() > 4 && p.m_variable2D->NBinsY() < 8 )
      nypads = 2;
  PlotUtils::GridCanvas* gridCanvas = new PlotUtils::GridCanvas(Form("2D_ErrorSummary_%s_%s_vs_%s_%s", tag.c_str(), p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str(), error_group_name.c_str()), nxpads, nypads, 1400, 600);
  std::cout << "N of pads = " << gridCanvas->GetPads().size() << "\n";
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
  std::vector<TH2D*> ErrorHists = DrawErrorSummary2D(p,
        hist, p.m_include_stat, true, ignore_threshold,
        p.m_do_cov_area_norm, error_group_name, p.m_do_frac_unc);
/*  double mx[11];
  std::cout << "Antes de GetMultipliers\n";
  std::vector<double> mult = GetMultipliers(ErrorHists[0]);
  std::cout << "Despues de GetMultipliers\n";  
  std::cout << "Antes del ciclo de multipliers\n";
  for (int i = 0; i < 11; ++i){
    if (mult[i] > 1000) mx[i] = 1.;
    else mx[i] = mult[i];
  }*/  
  std::vector<double> mult = GetMultipliers(ErrorHists);
  double mx[p.m_variable2D->NBinsY()];
  std::vector<int> bins;
  for (int i = 0; i < p.m_variable2D->NBinsY(); ++i){
    bins.push_back(i+1);
    if (mult[i] > 1000.) mx[i] = 1.;
    else mx[i] = mult[i];
  }
  bins.push_back(p.m_variable2D->NBinsY()+1);

  auto legend = new TLegend(0.88,0.9,1.,0.5);
  std::vector<std::string> names;
  names = error_names[error_group_name];
  if (error_group_name != ""){
    Plot2D(p, ErrorHists, names, Form("%s_%s", tag.c_str(), 
           error_group_name.c_str()),  "Fractional Uncertainty");
  }
  else {
    gridCanvas->SetRightMargin(0.12);
    gridCanvas->SetLeftMargin(0.07);
//    gridCanvas->ResetPads();
    ErrorHists[0]->GetXaxis()->SetLabelSize(0.03);
    gridCanvas->SetInterpadSpace(0.008);
    std::vector<std::string> names;

  gridCanvas->DrawOneHist( ErrorHists[0] , "HIST", false, bins, false, mx);
  legend->SetNColumns(1);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
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
  //gridCanvas->SetYLimits(ignore_threshold, 1.);
  gridCanvas->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  gridCanvas->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("%s %s vs %s", tag.c_str(), p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str()));
  std::string yaxis = "Fractional Uncertainty ";
  gridCanvas->SetYTitle(yaxis.c_str());
//  gridCanvas->ResetPads();
  gridCanvas->Draw();
  legend->Draw();

  gridCanvas->Print(Form("2D_ErrorSummary_%s_%s_vs_%s_%s.png", tag.c_str(),
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   error_group_name.c_str()));
  }
}

//-------------------------------------
//------Event Selection Plots----------
//-------------------------------------

void PlotVar_Selection2D(EventSelectionPlotInfo2D p, bool do_bg = true,
                         bool do_log_scale = false, bool do_tuned_bg = false,
			 bool do_bin_width_norm = true){

  std::cout << "Plotting Selection 2D" << p.m_variable2D->NameX() << "_vs_" << p.m_variable2D->NameY() << std::endl;
  int nxpads = 4;
  int nypads = p.m_variable2D->NBinsY()/3;
  if (p.m_variable2D->NBinsY() < 4){
      nxpads = p.m_variable2D->NBinsY();
      nypads = 1;
  }
  if (p.m_variable2D->NBinsY() > 4 && p.m_variable2D->NBinsY() < 8 )
      nypads = 2;
    PlotUtils::GridCanvas* StackedPlot =new PlotUtils::GridCanvas(Form("2D_Sel_%s_vs_%s",p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), nxpads, nypads, 1400, 600);
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
      new TH2D(p.m_variable2D->m_hists2D.m_selection_mc.hist->GetCVHistoWithError());

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
  tmp_bg->SetFillStyle(3000);
  tmp_bg->SetFillColor(46);
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

  auto legend = new TLegend(0.88,0.9,1.,0.7);
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->AddEntry(mc, "MC", "l");
  legend->AddEntry(tmp_bg, "BG", "f");
  legend->AddEntry(data, "Data", "p");

  mc_error->GetXaxis()->SetLabelSize(0.03);
  StackedPlot->SetTopMargin(0.1);
  StackedPlot->SetRightMargin(0.12);
  StackedPlot->SetLeftMargin(0.05);
  StackedPlot->SetInterpadSpace(0.008);
  StackedPlot->ResetPads();

  // Draw
  StackedPlot->DrawOneHist( mc_error, "E2", false, bins, false, mx);
  StackedPlot->DrawOneHist( mc, "SAME HIST", false, bins, false, mx);
  StackedPlot->DrawOneHist( tmp_bg, "SAME HIST", false, bins, false, mx);
  StackedPlot->DrawOneHist(data, "SAME E1", false, bins, false, mx);
  StackedPlot->DrawBinRanges(data, 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  StackedPlot->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  StackedPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("Selection %s vs %s", p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str()));
  std::string yaxis = "N Events / " + p.m_variable2D->m_unitsY;
  StackedPlot->SetYTitle(yaxis.c_str());
  StackedPlot->SetTitleSize(20);
  StackedPlot->ResetPads();
  StackedPlot->Draw();
  legend->Draw();

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

  SetErrorGroups2D(p.m_mnv_plotter, false);

  PlotUtils::MnvH2D* sel =
      (PlotUtils::MnvH2D*)p.m_variable2D->m_hists2D.m_selection_mc.hist->Clone(
          uniq());
  Plot_ErrorGroup2D(p, sel, "", "Sel", 0.0, 0.35);
//  Plot_ErrorGroup(p, sel, "LEGENDONLY", "Sel", 0.0, 0.2);
  Plot_ErrorGroup2D(p, sel, "2p2h", "Sel", 0.0, 0.01);
  Plot_ErrorGroup2D(p, sel, "Detector", "Sel", 0.0, 0.15);
  Plot_ErrorGroup2D(p, sel, "Flux", "Sel", 0.0, 0.15);
  Plot_ErrorGroup2D(p, sel, "Genie_FSI_nucleons", "Sel", 0.004, 0.06);
  Plot_ErrorGroup2D(p, sel, "Genie_FSI_pions", "Sel", 0.01, 0.1);
  Plot_ErrorGroup2D(p, sel, "Genie_InteractionModel", "Sel", 0.01, 0.2);
  Plot_ErrorGroup2D(p, sel, "Muon", "Sel", 0.0, 0.14);
  Plot_ErrorGroup2D(p, sel, "NonResPi", "Sel", 0.0, 0.06);
  Plot_ErrorGroup2D(p, sel, "RPA", "Sel", 0.0, 0.012);
  //  Plot_ErrorGroup(p, sel, "Michel", "Sel", 0.0, 0.15);
  Plot_ErrorGroup2D(p, sel, "GENIE", "Sel", 0.0, 0.30);
//  Plot_ErrorGroup2D(p, sel, "Target", "Sel", 0.0, 0.15);
  Plot_ErrorGroup2D(p, sel, "Response", "Sel", 0.0, 0.05);
  Plot_ErrorGroup2D(p, sel, "Diffractive", "Sel", 0.0, 0.15);
  Plot_ErrorGroup2D(p, sel, "PhysicsModel", "Sel", 0.0, 0.15);
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
  int nxpads = 4;
  int nypads = p.m_variable2D->NBinsY()/3;
  if (p.m_variable2D->NBinsY() < 4){
      nxpads = p.m_variable2D->NBinsY();
      nypads = 1;
  }
  if (p.m_variable2D->NBinsY() > 4 && p.m_variable2D->NBinsY() < 8 )
      nypads = 2;
  
  PlotUtils::GridCanvas* StackedPlot =new PlotUtils::GridCanvas(Form("2D_BG_%s_vs_%s",p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), nxpads, nypads, 1400, 600);
  StackedPlot->SetTitle("Background Subtracted");


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

  auto legend = new TLegend(0.92,0.65,0.98,0.78);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->AddEntry(effnum_cv, "MC", "l");
//legend->AddEntry(tmp_bg, "BG", "f");
  legend->AddEntry(bg_sub_data_w_stat_error, "Data", "p");

  StackedPlot->SetRightMargin(0.12);
  StackedPlot->SetLeftMargin(0.07);
  StackedPlot->SetInterpadSpace(0.008);
  StackedPlot->ResetPads();
  effnum_w_stat_error->GetXaxis()->SetLabelSize(0.03);
  // Draw
  StackedPlot->DrawOneHist( effnum_w_stat_error, "E2", false, bins, false, mx);
  StackedPlot->DrawOneHist( effnum_cv, "SAME HIST", false, bins, false, mx);
  StackedPlot->DrawOneHist(bg_sub_data_w_stat_error, "SAME E1", false, bins, false, mx);
  StackedPlot->DrawBinRanges(bg_sub_data_w_stat_error, 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  StackedPlot->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  StackedPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("Background Subtracted %s vs %s", p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str()));
  StackedPlot->SetYTitle(yaxis.c_str());
  StackedPlot->SetTitle("Background Subtracted");
  StackedPlot->ResetPads();
  StackedPlot->SetHistTexts();
  StackedPlot->Draw();
  legend->Draw();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  StackedPlot->Print(Form("2D_BG_%s_vs_%s_%s.png",
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));

}

// Error Summary BG
void PlotBG_ErrorSummary2D(EventSelectionPlotInfo2D p, bool do_tuned = false) {
  SetErrorGroups2D(p.m_mnv_plotter, false);
  PlotUtils::MnvH2D* bg;

  if (do_tuned)
    bg =
        (PlotUtils::MnvH2D*)p.m_variable2D->m_hists2D.m_tuned_bg->Clone("tuned_bg");
  else
    bg = (PlotUtils::MnvH2D*)p.m_variable2D->m_hists2D.m_bg.hist->Clone("bg");

  std::string name = p.m_variable2D->NameX() + "_vs_" + p.m_variable2D->NameY();
  std::string tuned_str = do_tuned ? "BGTuned" : "BGUntuned";

  // name, ignore threshold, ymax
//  Plot_ErrorGroup(p, bg, "LEGENDONLY", tuned_str, 0.0);
  Plot_ErrorGroup2D(p, bg, "", tuned_str, 0.0);
//  Plot_ErrorGroup(p, bg, "2p2h", tuned_str, 0.0, 0.05);
//  Plot_ErrorGroup(p, bg, "Detector", tuned_str, detector_threshold,
//                  detector_ymax);
//  Plot_ErrorGroup(p, bg, "Flux", tuned_str, 0.0, 0.15);
//  Plot_ErrorGroup(p, bg, "Genie_FSI_pions", tuned_str, FSI_threshold, FSI_ymax);
//  Plot_ErrorGroup(p, bg, "Genie_FSI_nucleons", tuned_str, FSI_threshold,
//                  FSI_ymax);
//  Plot_ErrorGroup(p, bg, "Genie_InteractionModel", tuned_str, Int_threshold,
//                  Int_ymax);
//  Plot_ErrorGroup(p, bg, "NonResPi", tuned_str, 0.0, 0.1);
//  Plot_ErrorGroup(p, bg, "RPA", tuned_str, 0.0, 0.1);
//  Plot_ErrorGroup(p, bg, "Target", tuned_str, 0.0, 0.15);
//  Plot_ErrorGroup(p, bg, "Response", tuned_str, 0.0, 0.30);
//  Plot_ErrorGroup(p, bg, "Diffractive", tuned_str, 0.0, 0.05);
//  Plot_ErrorGroup(p, bg, "PhysicsModel", tuned_str, 0.0, 0.05);
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
  std::cout << "Plotting Unfoldeding " << p.m_variable2D->NameX() << " vs " << p.m_variable2D->NameY()<< std::endl;
  int nxpads = 4;
  int nypads = p.m_variable2D->NBinsY()/3;
  if (p.m_variable2D->NBinsY() < 4){
      nxpads = p.m_variable2D->NBinsY();
      nypads = 1;
  }
  if (p.m_variable2D->NBinsY() > 4 && p.m_variable2D->NBinsY() < 8 )
      nypads = 2;
  PlotUtils::GridCanvas* GridPlot =new PlotUtils::GridCanvas(Form("2D_Unfold_%s_vs_%s",p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), nxpads, nypads, 1400, 600);
  GridPlot->SetTitle("Unfolding");

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

  auto legend = new TLegend(0.92,0.65,0.98,0.78);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->AddEntry(effnum_cv, "MC", "l");
//legend->AddEntry(tmp_bg, "BG", "f");
  legend->AddEntry(unfolded_w_stat_error, "Data", "p");

  GridPlot->SetRightMargin(0.12);
  GridPlot->SetLeftMargin(0.07);
  GridPlot->SetInterpadSpace(0.008);
  GridPlot->ResetPads();
  effnum_w_stat_error->GetXaxis()->SetLabelSize(0.03);
  // Draw
  GridPlot->DrawOneHist( effnum_w_stat_error, "E2", false, bins, false, mx);
  GridPlot->DrawOneHist( effnum_cv, "SAME HIST", false, bins, false, mx);
  GridPlot->DrawOneHist(unfolded_w_stat_error, "SAME E1", false, bins, false, mx);
  GridPlot->DrawBinRanges(unfolded_w_stat_error, 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  GridPlot->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  GridPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("Unfolded %s vs %s", p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str()));
  GridPlot->SetYTitle(yaxis.c_str());
//GridPlot->SetTitle(Form("Background Subtracted %s", GetSignalName(p.m_signal_definition).c_str()));
  GridPlot->ResetPads();
  GridPlot->Draw();
  legend->Draw();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  GridPlot->Print(Form("2D_Unfolded_%s_vs_%s_%s.png",
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));


}
void PlotUnfolded_ErrorSummary2D(EventSelectionPlotInfo2D p) {
  SetErrorGroups2D(p.m_mnv_plotter, false);
  std::cout << "In errorSummary2D Unfolded \n";
  PlotUtils::MnvH2D* unf =
      (PlotUtils::MnvH2D*)p.m_variable2D->m_hists2D.m_unfolded->Clone(uniq());

//  Plot_ErrorGroup2D(p, unf, "LEGENDONLY", "Unfolded", 0.0);
  Plot_ErrorGroup2D(p, unf, "", "Unfolded", 0.0);
//  Plot_ErrorGroup2D(p, unf, "2p2h", "Unfolded", 0.0, 0.02);
//  Plot_ErrorGroup2D(p, unf, "Detector", "Unfolded", 0.0, 0.08);
//  Plot_ErrorGroup2D(p, unf, "Flux", "Unfolded", 0.0, 0.06);
//  Plot_ErrorGroup2D(p, unf, "Genie_FSI_nucleons", "Unfolded", 0.01, 0.1);
//  Plot_ErrorGroup2D(p, unf, "Genie_FSI_pions", "Unfolded", 0.01, 0.1);
//  Plot_ErrorGroup2D(p, unf, "Genie_InteractionModel", "Unfolded", 0.01, 0.1);
//  Plot_ErrorGroup2D(p, unf, "NonResPi", "Unfolded", 0.0, 0.1);
//  Plot_ErrorGroup2D(p, unf, "RPA", "Unfolded", 0.0, 0.02);
//  Plot_ErrorGroup2D(p, unf, "Target", "Unfolded", 0.0, 0.1);
//  Plot_ErrorGroup2D(p, unf, "Response", "Unfolded", 0.0, 0.24);
//  Plot_ErrorGroup2D(p, unf, "Diffractive", "Unfolded", 0.0, 0.02);
//  Plot_ErrorGroup2D(p, unf, "PhysicsModel", "Unfolded", 0.0, 0.02);
//  Plot_ErrorGroup2D(p, unf, "Muon", "Unfolded", 0.0, 0.14);
}
void PlotMC2D(EventSelectionPlotInfo2D p, MnvH2D* mc, std::string ylabel = "",
                   std::string outdir = ".", double ymax = -1,
                   bool do_log_scale = false, bool do_bin_width_norm = true) {
  std::cout << ylabel + " 2D " << p.m_variable2D->NameX() << " vs " << p.m_variable2D->NameY()<< std::endl;
  int nxpads = 4;
  int nypads = p.m_variable2D->NBinsY()/4;
  if (p.m_variable2D->NBinsY() < 4){
      nxpads = p.m_variable2D->NBinsY();
      nypads = 1;
  }
  if (p.m_variable2D->NBinsY() > 4 && p.m_variable2D->NBinsY() < 8 )
      nypads = 2;
  PlotUtils::GridCanvas* GridPlot =new PlotUtils::GridCanvas(Form("2D_%s_%s_vs_%s", ylabel.c_str(), p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), nxpads, nypads, 1400, 600);
  GridPlot->SetTitle(Form("%s", ylabel.c_str()));
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
  auto legend = new TLegend(0.92,0.65,0.98,0.78);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->AddEntry(hist_w_cv, "MC", "l");
  
  GridPlot->SetRightMargin(0.12);
  GridPlot->SetLeftMargin(0.07);
  GridPlot->SetInterpadSpace(0.008);
  GridPlot->ResetPads();
  hist_w_tot_error->GetXaxis()->SetLabelSize(0.03);
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
             + " vs " + p.m_variable2D->m_hists2D.m_xlabelY.c_str());
  GridPlot->SetYTitle(yaxis.c_str());
  GridPlot->SetTitle(Form("%s", ylabel.c_str()));
//  GridPlot->ResetPads();
  GridPlot->Draw();
  legend->Draw();
  GridPlot->Update();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  GridPlot->Print(Form("2D_%s_%s_vs_%s_%s.png", ylabel.c_str(),
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));
}

void PlotEfficiency_ErrorSummary2D(EventSelectionPlotInfo2D p) {
  SetErrorGroups2D(p.m_mnv_plotter, false);
  std::cout << "In errorSummary2D Efficiency \n";
  PlotUtils::MnvH2D* eff =
      (PlotUtils::MnvH2D*)p.m_variable2D->m_hists2D.m_efficiency->Clone(uniq());

//  Plot_ErrorGroup2D(p, eff, "LEGENDONLY", "Unfolded", 0.0);
  Plot_ErrorGroup2D(p, eff, "", "Efficiency", 0.0);
//  Plot_ErrorGroup2D(p, eff, "2p2h", "Unfolded", 0.0, 0.02);
//  Plot_ErrorGroup2D(p, eff, "Detector", "Unfolded", 0.0, 0.08);
//  Plot_ErrorGroup2D(p, eff, "Flux", "Unfolded", 0.0, 0.06);
//  Plot_ErrorGroup2D(p, eff, "Genie_FSI_nucleons", "Unfolded", 0.01, 0.1);
//  Plot_ErrorGroup2D(p, eff, "Genie_FSI_pions", "Unfolded", 0.01, 0.1);
//  Plot_ErrorGroup2D(p, eff, "Genie_InteractionModel", "Unfolded", 0.01, 0.1);
//  Plot_ErrorGroup2D(p, eff, "NonResPi", "Unfolded", 0.0, 0.1);
//  Plot_ErrorGroup2D(p, eff, "RPA", "Unfolded", 0.0, 0.02);
//  Plot_ErrorGroup2D(p, eff, "Target", "Unfolded", 0.0, 0.1);
//  Plot_ErrorGroup2D(p, eff, "Response", "Unfolded", 0.0, 0.24);
//  Plot_ErrorGroup2D(p, eff, "Diffractive", "Unfolded", 0.0, 0.02);
//  Plot_ErrorGroup2D(p, eff, "PhysicsModel", "Unfolded", 0.0, 0.02);
//  Plot_ErrorGroup2D(p, eff, "Muon", "Unfolded", 0.0, 0.14);
}

void Plot_CrossSection2D(EventSelectionPlotInfo2D p, MnvH2D* data, MnvH2D* mc,
                       std::string outdir = ".", double ymax = -1,
                       bool do_log_scale = false,
                       bool do_bin_width_norm = true) {

  int nxpads = 4;
  int nypads = p.m_variable2D->NBinsY()/4;
  if (p.m_variable2D->NBinsY() < 4){
      nxpads = p.m_variable2D->NBinsY();
      nypads = 1;
  }
  if (p.m_variable2D->NBinsY() > 4 && p.m_variable2D->NBinsY() < 8 )
      nypads = 2;
  std::cout << "Plotting Cross Section " << p.m_variable2D->NameX() << " vs " << p.m_variable2D->NameY()<< std::endl;
  PlotUtils::GridCanvas* GridPlot =new PlotUtils::GridCanvas(Form("2D_CrossSection_%s_vs_%s",p.m_variable2D->NameX().c_str(), p.m_variable2D->NameY().c_str()), nxpads, nypads, 1400, 600);
  GridPlot->SetTitle("Cross Section");
  //GridPlot->SetInterpadSpace(0.0005);
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
  std::string yaxis = "";
  if (do_bin_width_norm) {
    xsec_data_w_tot_error->Scale(1.e42, "width");
    xsec_data_w_stat_error->Scale(1.e42, "width");
    xsec_mc_w_stat_error->Scale(1.e42, "width");
    xsec_mc_cv->Scale(1.e42, "width");
    std::string div = Form("#frac{d^{2}#sigma}{d%sd%s}",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_hists2D.m_xlabelY.c_str());
    // Y label
    yaxis = div + " (10^{-42} cm^{2}/" + p.m_variable2D->m_unitsY + "*"
	    + p.m_variable2D->m_unitsY + "/nucleon)";
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

  auto legend = new TLegend(0.88,0.65,0.98,0.78);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->AddEntry(xsec_mc_cv, "MC", "l");
//legend->AddEntry(tmp_bg, "BG", "f");
  legend->AddEntry(xsec_data_w_stat_error, "Data", "p");

  GridPlot->SetRightMargin(0.12);
  GridPlot->SetLeftMargin(0.07);
  GridPlot->SetTopMargin(0.1);
  GridPlot->SetInterpadSpace(0.008);
  GridPlot->ResetPads();
  xsec_mc_w_stat_error->GetXaxis()->SetLabelSize(0.03);
  // Draw
  GridPlot->DrawOneHist( xsec_mc_w_stat_error, "SAME E2", false, bins, false, mx);
  GridPlot->DrawOneHist( xsec_mc_cv, "SAME HIST", false, bins, false, mx);
  GridPlot->DrawOneHist(xsec_data_w_tot_error, "SAME E1", false, bins, false, mx);
  GridPlot->DrawBinRanges(xsec_data_w_tot_error, 2, bins, Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelY.c_str(), p.m_variable2D->m_unitsY.c_str()), 0.03, ".2f", 2);
  GridPlot->DrawMultipliers(sizeof(mx)/sizeof(double), mx);
  GridPlot->SetXTitle(Form("%s (%s)",p.m_variable2D->m_hists2D.m_xlabelX.c_str(), p.m_variable2D->m_unitsX.c_str()));
  p.SetTitle(Form("Cross Section %s vs %s", p.m_variable2D->m_hists2D.m_xlabelX.c_str(),
             p.m_variable2D->m_hists2D.m_xlabelY.c_str()));
  GridPlot->SetYTitle(yaxis.c_str());
  GridPlot->SetTitleSize(20);
//  GridPlot->SetTitle(Form("Background Subtracted %s", GetSignalName(p.m_signal_definition).c_str()));
//  GridPlot->ResetPads();
  GridPlot->Draw();
  legend->Draw();
  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  GridPlot->Print(Form("2D_CrossSection_%s_vs_%s_%s.png",
			   p.m_variable2D->NameX().c_str(),
                           p.m_variable2D->NameY().c_str(),
 			   bwn_str.c_str()));

}

void PlotCrossSection_ErrorSummary2D(EventSelectionPlotInfo2D p) {
  SetErrorGroups2D(p.m_mnv_plotter, false);
  std::cout << "In errorSummary2D Cross Section \n";
  PlotUtils::MnvH2D* xsec =
      (PlotUtils::MnvH2D*)p.m_variable2D->m_hists2D.m_cross_section->Clone(uniq());

//  Plot_ErrorGroup2D(p, xsec, "LEGENDONLY", "Unfolded", 0.0);
  Plot_ErrorGroup2D(p, xsec, "", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "Flux", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "Pion_Reconstruction", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "GENIE_FSI", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "Muon", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "Others", "CrossSection", 0.0);

  std::vector<TH2D*> XsecModelsvec;
  XsecModelsvec.push_back(AddingGroups(p.m_mnv_plotter, *xsec, "NonRESPi"));
  XsecModelsvec.push_back(AddingGroups(p.m_mnv_plotter, *xsec, "RESPi"));
  XsecModelsvec.push_back(AddingGroups(p.m_mnv_plotter, *xsec, "CCQE"));
  XsecModelsvec.push_back(AddingGroups(p.m_mnv_plotter, *xsec, "Coherent-Diffractive"));
  XsecModelsvec.push_back(AddingGroups(p.m_mnv_plotter, *xsec, "DIS-Hadronization"));
  XsecModelsvec.push_back(AddingGroups(p.m_mnv_plotter, *xsec, "MnvTunes"));
  XsecModelsvec.push_back(AddingGroups(p.m_mnv_plotter, *xsec, "Elastic"));

//  AddingGroups(p.m_mnv_plotter, *xsec, "NonRESPi");
//  AddingGroups(p.m_mnv_plotter, *xsec, "RESPi");
//  AddingGroups(p.m_mnv_plotter, *xsec, "CCQE");
//  AddingGroups(p.m_mnv_plotter, *xsec, "Coherent-Diffractive");
//  AddingGroups(p.m_mnv_plotter, *xsec, "DIS-Hadronization");
//  AddingGroups(p.m_mnv_plotter, *xsec, "MnvTunes");
//  AddingGroups(p.m_mnv_plotter, *xsec, "Elastic");

  SetErrorGroups2D(p.m_mnv_plotter, true);
  std::vector<std::string> names{"NonRESPi", "RESPi", "CCQE",
                                 "Coherent-Diffractive", "DIS-Hadronization",
                                 "MnvTunes", "Elastic"};
//  Plot_ErrorGroup2D(p, xsec, "Cross_Section_Models", "CrossSection", 0.0);
  NiceColors(p, XsecModelsvec, names); 

  Plot2D(p, XsecModelsvec, names, "CrossSection_Cross_Section_Models", 
         "Fractional Uncertainty");
  Plot_ErrorGroup2D(p, xsec, "RESPi", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "NonRESPi", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "CCQE", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "Coherent-Diffractive", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "DIS-Hadronization", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "MnvTunes", "CrossSection", 0.0);
  Plot_ErrorGroup2D(p, xsec, "Elastic", "CrossSection", 0.0);

}


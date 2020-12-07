#ifndef WSidebandFitter_cxx
#define WSidebandFitter_cxx

#include "WSidebandFitter.h"
#include <iostream>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TF1.h"
#include <iomanip> // setprecision


//==============================================================================
// Helper Functions
//==============================================================================
//you take ownership of this hist and must delete it
TH1D* SumHists(TObjArray* hists) {
  TH1D* total = nullptr;
  for (int itH = 0; itH != hists->GetEntries(); ++itH) {
    if (itH == 0) {
      total = new TH1D( *(TH1D*)(*hists)[itH]->Clone("my_tmp_name") );
    }
    else {
      total->Add( (TH1D*)(*hists)[itH] );
    }
  }
  return total;
}


// TFractionFitter. This is an alternative fit to Minimizer.
// Currently this method is not used.
// The result of this fit are fractions that each MC source should contribute
// to the total. 
// e.g. before fractions (0.33/0.33/0.33) --> sums to 1
//      after fractions  (0.35/0.28/0.36) --> sums to 1
// To convert to a scale:
// scale = new fraction/old fraction
void FractionFit(TH1D* data_hist, TObjArray* mc_templates,
                 std::vector<double>& weights, std::vector<double>& errs,
                 bool be_verbose,
                 double& midscale, double& hiscale) {

  if(be_verbose) std::cout<<"   Entering CalculateFitWeight"<<std::endl;

  
  //PlotTH1_1(data_hist, "data");

  for (int w = 0; w != mc_templates->GetEntries(); ++w) {
    std::string tag = std::to_string(w);
    TH1D* h_temp = new TH1D();
    h_temp = (TH1D*)mc_templates->At(w)->Clone("h_temp");;
    //PlotTH1_1(h_temp,tag);
    delete h_temp;
  }

  //cleanup
  weights.clear();
  errs.clear();
  
  //is this a onepi analysis?
  bool do1PiAna = true;

  //fit limits 
  int lower_limit = 10;
  int upper_limit = 23;
  
  //Create a TFractionFitter and fit the templates
  if(be_verbose) std::cout<<"    Calling TFractionFitter"<<std::endl;
  TFractionFitter* fit = new TFractionFitter( data_hist, mc_templates, "V");
  if(be_verbose) std::cout<<"    Leaving TFractionFitter"<<std::endl;
  fit->SetRangeX(lower_limit, upper_limit);

  //add the total MC hist
  TH1D* hist_mc = SumHists(mc_templates);
  double total_integral = hist_mc->Integral(lower_limit,upper_limit);

  // Constrain components
  // Fix signal and Wgenie < 1.4
  // Float 1.4 < Wgenie < 1.8 and Wgenie > 1.8
  for (int w = 0; w != mc_templates->GetEntries(); ++w) {
    std::cout<<setprecision(3);
    TH1D* tmp_hist = (TH1D*)(*mc_templates)[w];
    double initial_frac = tmp_hist->Integral(lower_limit,upper_limit)/total_integral;
    if(be_verbose)
      std::cout<<"initial " << w << " fraction: "<<initial_frac  << std::endl;
    switch (w) {
      case kSigParamId:
        fit->Constrain(w+1,initial_frac*0.99999, initial_frac*1.00001);
        break;
      case kLoWParamId:
        fit->Constrain(w+1,initial_frac*0.99999, initial_frac*1.00001);
        break;
      case kMidWParamId:
      case kHiWParamId:
        fit->Constrain(w+1, 0., 1.);
        break;
    }
  }

  // Do indices start at 1? is data reserved for 0?

  //fit->Constrain(kLoWParamId+1,  0.,1.);
  //fit->Constrain(kMidWParamId+1, 0.,1.);
  //fit->Constrain(kHiWParamId+1,  0.,1.);

  //fit->Constrain(kLoWParamId+1,  0.05,0.2);
  fit->Constrain(kMidWParamId+1, 0.,1.);
  fit->Constrain(kHiWParamId+1,  0.,1.);
  //fit->Constrain(kSigParamId+1,  0.0,1.0);
  int fit_status = fit->Fit();

  if (be_verbose) std::cout<<"    CV fit status "<< fit_status <<" and chi2/ndf "<<fit->GetChisquare()<<"/"<<fit->GetNDF()<< std::endl;

  //fill the weights
  double fit_frac, fit_err, N_raw, error;
  for (int w = 0; w != mc_templates->GetEntries(); ++w) {
    fit->GetResult(w, fit_frac, fit_err);
    N_raw = ((TH1D*)(*mc_templates)[w])->Integral(lower_limit, upper_limit);
    
    if (N_raw == 0.0) {
      std::cout<<"WARNING: N is 0 for template "<<w<<std::endl;
      weights.push_back(0.0);
      errs.push_back(0.0);
    }

    TH1D* tmp_hist = (TH1D*)(*mc_templates)[w];
    double initial_frac =  tmp_hist->Integral(lower_limit,upper_limit)/total_integral;
    
    //Helpful output
    if(be_verbose) {
      switch (w) {
        case kLoWParamId:
          std::cout<<"fit low fraction: "<<fit_frac<<" +/- "<<fit_err << std::endl;
          std::cout<<"fit low scale: "<<fit_frac/initial_frac<< std::endl;
          break;
        case kMidWParamId:
          std::cout<<"fit midw fraction: "<<fit_frac<<" +/- "<<fit_err << std::endl;
          std::cout<<"fit midw scale: "<<fit_frac/initial_frac<< std::endl;
          midscale = fit_frac/initial_frac;
          break;
        case kHiWParamId:
          std::cout<<"fit hiw fraction: "<<fit_frac<<" +/- "<<fit_err << std::endl;
          std::cout<<"fit hiw scale: "<<fit_frac/initial_frac<< std::endl;
          hiscale = fit_frac/initial_frac;
          break;
        case kSigParamId:
          std::cout<<"fit sig fraction: "<<fit_frac<<" +/- "<<fit_err << std::endl;
          std::cout<<"fit sig scale: "<<fit_frac/initial_frac<< std::endl;
          break;
      }
    }
    
    error = sqrt( pow(fit_err/N_raw,2.0) + pow(fit_frac,2.0)/pow(N_raw,3.0) );
    
    weights.push_back(fit_frac/N_raw);
    errs.push_back(error);
  }

  TH1F* result = (TH1F*) fit->GetPlot();
  data_hist->Draw("Ep");
  result->Draw("same");

  if(be_verbose) std::cout<<"   Leaving CalculateFitWeight"<<std::endl;
    
  delete hist_mc;
  delete fit;  
    
  return;
}    
//==============================================================================

TH1D* WSidebandFitter::m_data = NULL;  
TH1D* WSidebandFitter::m_sig  = NULL;
TH1D* WSidebandFitter::m_hiW  = NULL;
TH1D* WSidebandFitter::m_midW  = NULL;
TH1D* WSidebandFitter::m_loW  = NULL;
double WSidebandFitter::m_pot_scale = 1.0;
double WSidebandFitter::m_chi2 = 1.0;
int WSidebandFitter::m_ndf = 1;


WSidebandFitter::WSidebandFitter() {} //default constructor


WSidebandFitter::WSidebandFitter(const CVUniverse& universe, Histograms hists,
                                 const double pot_scale) {
  m_data = new TH1D(*hists.m_wsidebandfit_data);
  m_sig  = new TH1D(*hists.m_wsidebandfit_sig.univHist(&universe)); 
  m_hiW  = new TH1D(*hists.m_wsidebandfit_hiW.univHist(&universe));
  m_midW = new TH1D(*hists.m_wsidebandfit_midW.univHist(&universe));
  m_loW  = new TH1D(*hists.m_wsidebandfit_loW.univHist(&universe));
  m_pot_scale = pot_scale;
}


WSidebandFitter::~WSidebandFitter() {} //destructor


double WSidebandFitter::MinimizerFunc(const double* par) {
  double loW_scale  = par[kLoWParamId];
  double midW_scale = par[kMidWParamId];
  double hiW_scale  = par[kHiWParamId];

  TH1D* h_temp_data = (TH1D*)m_data->Clone("data_temp");
  TH1D* h_temp_sig  = (TH1D*)m_sig->Clone("sig_temp");
  TH1D* h_temp_hiW  = (TH1D*)m_hiW->Clone("hiW_temp");
  TH1D* h_temp_midW = (TH1D*)m_midW->Clone("midW_temp");
  TH1D* h_temp_loW  = (TH1D*)m_loW->Clone("loW_temp");

  h_temp_sig->Scale(m_pot_scale);
  h_temp_hiW->Scale(m_pot_scale*hiW_scale);
  h_temp_midW->Scale(m_pot_scale*midW_scale);
  h_temp_loW->Scale(m_pot_scale*loW_scale);

  double tot_data = 0.0;
  double tot_sig  = 0.0;
  double tot_hiW  = 0.0;
  double tot_midW  = 0.0;
  double tot_loW  = 0.0;

  double chi2 = 0.0;
  int ndf = 0;
  int nbins = h_temp_data->GetNbinsX();// + 1;
  
  for (int i = 0; i < nbins; ++i) {
    double data_entries = h_temp_data->GetBinContent(i);
    double error        = h_temp_data->GetBinError(i);
    double sig_entries  = h_temp_sig ->GetBinContent(i);
    double hiW_entries  = h_temp_hiW ->GetBinContent(i);
    double midW_entries = h_temp_midW->GetBinContent(i);
    double loW_entries  = h_temp_loW ->GetBinContent(i);
    double tot_mc_entries = sig_entries + hiW_entries + midW_entries + loW_entries;
    double mc_error = sqrt(fabs(tot_mc_entries));
    

    //if (data_entries == 0) continue;
    if (data_entries < 5) continue;

    tot_data += data_entries;
    tot_sig  += sig_entries;
    tot_hiW  += hiW_entries;
    tot_midW += midW_entries;
    tot_loW  += loW_entries;

    chi2 += pow((sig_entries + hiW_entries + midW_entries + loW_entries - data_entries),2) / (pow(error, 2) + pow(mc_error,2));
    ndf++;

  }
  double tot_chi2 = pow((tot_sig + tot_hiW + tot_midW + tot_loW - tot_data) / sqrt(tot_data), 2);
  m_chi2 = chi2;
  m_ndf  = ndf;
  //std::cout << m_chi2 << " " << m_ndf << " " << tot_chi2 << " " << hiW_scale << " " << loW_scale << "\n";
  return chi2;
}


void WSidebandFitter::Fit() {
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  ROOT::Math::Functor func(&MinimizerFunc, 3);
  min->SetFunction(func);

  min->SetMaxFunctionCalls(10000000);// Alex 10M
  min->SetMaxIterations(100000);// Alex 0.1M
  min->SetTolerance(0.01);// Alex 0.1micro
  min->SetPrintLevel(2);
  //min->SetErrorDef(1);
  //min->SetStrategy(2);
  //min->SetLimitedVariable(0, "usplastic", 1.0, 0.01, 0.1, 3.0);
  //min->SetLimitedVariable(1, "dsplastic", 1.0, 0.01, 0.1, 3.0);
  min->SetVariable(kLoWParamId,  "loW", 1.0, 0.001/*, 0.1, 3.0*/);
  min->SetVariable(kMidWParamId, "midW", 1.0, 0.001/*, 0.1, 3.0*/);
  min->SetVariable(kHiWParamId,  "hiW", 1.0, 0.001/*, 0.1, 3.0*/);
  min->SetVariableLimits(kLoWParamId, 1., 10);
  min->SetVariableLimits(kMidWParamId, 0., 10.);
  min->SetVariableLimits(kHiWParamId, 0., 10.);

  min->Minimize();

  cout << "==============================" << endl;
  cout << "chi2 = " << m_chi2 << " ndf = " << m_ndf << endl;
  cout << "chi2/ndf = " << m_chi2/m_ndf    << endl;
  cout << "==============================" << endl;

  m_fit_min = min->MinValue();
  m_fit_scale[kLoWParamId]      = min->X()[kLoWParamId];
  m_fit_scale[kMidWParamId]     = min->X()[kMidWParamId];
  m_fit_scale[kHiWParamId]      = min->X()[kHiWParamId];
  m_fit_scale_err[kLoWParamId]  = min->Errors()[kLoWParamId];
  m_fit_scale_err[kMidWParamId] = min->Errors()[kMidWParamId];
  m_fit_scale_err[kHiWParamId]  = min->Errors()[kHiWParamId];
}

#endif // WSidebandFitter_cxx


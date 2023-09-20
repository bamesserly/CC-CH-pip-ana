#ifndef Variable_cxx
#define Variable_cxx

#include "Variable.h"
#include "TruthMatching.h"
#include <stdlib.h> // exit


// CTOR -- default
template <class T>
Variable<T>::Variable()
  : m_pointer_to_GetValue(&T::GetDummyVar),
    m_label(),
    m_units(),
    m_hists(),
    m_is_true(false)
{}


// CTOR -- uniform binning
template <class T>
Variable<T>::Variable(const std::string label, const std::string xlabel,
                   const std::string units,
                   const int nbins, const double xmin, const double xmax, 
                   PointerToTFunction p,
                   const bool is_true)
  : m_pointer_to_GetValue(p),
    m_label(label),
    m_units(units),
    m_hists(m_label, xlabel, nbins, xmin, xmax),
    m_is_true(is_true)
{}


// CTOR -- variable binning
template <class T>
Variable<T>::Variable(const std::string label, const std::string xlabel, 
                   const std::string units,
                   const TArrayD& bins_array,
                   PointerToTFunction p, 
                   const bool is_true)
  : m_label(label),
    m_units(units),
    m_pointer_to_GetValue(p),
    m_hists(m_label, xlabel, bins_array),
    m_is_true(is_true)
{}


// GetValue defines this variable
template <class T>
double Variable<T>::GetValue (const T& t, const std::string& s) const { 
  return m_pointer_to_GetValue(t, s); 
}


// Histogram Initialization
template<class T>
template<class U>
void Variable<T>::InitializeAllHists(U systematic_univs, U systematic_univs_truth) {
  m_hists.InitializeAllHists(systematic_univs, systematic_univs_truth);
}


template<class T>
template<class U>
void Variable<T>::InitializeSidebandHists(U systematic_univs) {
  m_hists.InitializeSidebandHists(systematic_univs);
}


template<class T>
void Variable<T>::InitializeStackedHists(){
  m_hists.InitializeStackedHists();
}

template<class T>
void Variable<T>::InitializeDataHists(){
  m_hists.InitializeDataHists();
}


// Histogram Access
template <class T>
template <class U>
PlotUtils::MnvH1D* Variable<T>::GetStackComponentHist(U type) const {
  std::map<T, PlotUtils::MnvH1D*> stack_map = m_hists.GetStackMap(type);
  return stack_map[type];
}


// Save with the object names that hists were initialized with
template<class T>
void Variable<T>::WriteMCHists(TFile& fout) const {
  fout.cd();
  m_hists.m_selection_mc.hist -> Write();
  m_hists.m_bg.hist           -> Write();
  m_hists.m_bg_loW.hist       -> Write();
  m_hists.m_bg_midW.hist      -> Write();
  m_hists.m_bg_hiW.hist       -> Write();
  m_hists.m_effnum.hist       -> Write();
  m_hists.m_effden.hist       -> Write();
  if (Name() == sidebands::kFitVarString) {
    m_hists.m_wsidebandfit_sig.hist ->Write();
    m_hists.m_wsidebandfit_loW.hist ->Write();
    m_hists.m_wsidebandfit_midW.hist->Write();
    m_hists.m_wsidebandfit_hiW.hist ->Write();
    for(auto i : m_hists.m_stacked_wsideband.m_hist_map) {
      std::string legend_name = GetTruthClassification_LegendLabel(i.first);
      std::string short_name  = GetTruthClassification_Name(i.first);
      std::cout << Form("%s_%s", m_label.c_str(), short_name.c_str()) << "\n";
      i.second->Write(Form("%s_%s", m_label.c_str(), short_name.c_str()));
    }
  }
  if (!m_is_true && Name() != sidebands::kFitVarString)
    m_hists.m_migration.hist->Write();
}


template<class T>
void Variable<T>::LoadDataHistsFromFile(TFile& fin) {
  m_hists.LoadDataHistsFromFile(fin);
}


template<class T>
void Variable<T>::LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands) {
  m_hists.LoadMCHistsFromFile(fin, error_bands, m_is_true);
}


template <class T>
template <class U>
TObjArray Variable<T>::GetStackArray(U type) const {
  if(std::is_same<U, WType>::value)                     return m_hists.m_stacked_w.m_hist_array;
  else if(std::is_same<U, SignalBackgroundType>::value) return m_hists.m_stacked_sigbg.m_hist_array;
  else if(std::is_same<U, WBackgroundType>::value)      return m_hists.m_stacked_wbg.m_hist_array;
  else if(std::is_same<U, MesonBackgroundType>::value)  return m_hists.m_stacked_mesonbg.m_hist_array;
  else if(std::is_same<U, HadronType>::value)           return m_hists.m_stacked_hadron.m_hist_array;
  else if(std::is_same<U, FSParticleType>::value)       return m_hists.m_stacked_fspart.m_hist_array;
  else if(std::is_same<U, ChannelType>::value)          return m_hists.m_stacked_channel.m_hist_array;
  else if(std::is_same<U, NPionsType>::value)           return m_hists.m_stacked_npi.m_hist_array;
  else if(std::is_same<U, NPi0Type>::value)             return m_hists.m_stacked_npi0.m_hist_array;
  else if(std::is_same<U, NPipType>::value)             return m_hists.m_stacked_npip.m_hist_array;
  else if(std::is_same<U, WSidebandType>::value)        return m_hists.m_stacked_wsideband.m_hist_array;
  else if(std::is_same<U, CoherentType>::value)         return m_hists.m_stacked_coherent.m_hist_array;
  else{
    std::cerr << "GetStackArray: Unknown truth type.\n";
    std::exit(1);
  }
}

#endif // Variable_cxx

#ifndef StackedHistogram_cxx
#define StackedHistogram_cxx

#include "StackedHistogram.h"

// CTOR -- default
template <typename T>
StackedHistogram<T>::StackedHistogram()
  : m_label(),
    m_xlabel(),
    m_bins_array(0),
    m_nhists(0),
    m_color_scheme(0),
    m_hist_map(),
    m_hist_array()
{}


// CTOR -- uniform binning
template <typename T>
StackedHistogram<T>::StackedHistogram(std::string label, std::string xlabel, 
                                      int nbins, double xmin, double xmax,
                                      int nhists, int color_scheme)
  : m_label(label),
    m_xlabel(xlabel),
    m_bins_array(MakeUniformBinArray(nbins, xmin, xmax)),
    m_nhists(nhists),
    m_color_scheme(color_scheme),
    m_hist_map(),
    m_hist_array()
{
  Initialize();
}


// CTOR -- variable binning
template <typename T>
StackedHistogram<T>::StackedHistogram(std::string label, std::string xlabel, 
                                      const TArrayD& bins_array,
                                      int nhists, int color_scheme)
  : m_label(label),
    m_xlabel(xlabel),
    m_bins_array(GetSortedArray(bins_array)),
    m_nhists(nhists),
    m_color_scheme(color_scheme),
    m_hist_map(),
    m_hist_array()
{
  Initialize();
}


// Add component hists to stack
template <typename T>
void StackedHistogram<T>::Initialize() {
  for (int i=0; i != m_nhists; ++i){
    T type = static_cast<T>(i);
    PlotUtils::MnvH1D* component_hist = MakeStackComponentHist(type);
    m_hist_map[type] = component_hist;
    m_hist_array.Add(component_hist);
  }
}


// Make component hist
template <typename T>
PlotUtils::MnvH1D* StackedHistogram<T>::MakeStackComponentHist(const T type) const  {
  std::string legend_name = GetTruthClassification_LegendLabel(type);
  std::string short_name  = GetTruthClassification_Name(type);

  PlotUtils::MnvH1D* hist = nullptr;

  hist = new PlotUtils::MnvH1D( 
      uniq(), Form("%s_%s", m_label.c_str(), short_name.c_str()), 
      NBins(), m_bins_array.GetArray());
  
  SetHistColorScheme(hist, int(type), m_color_scheme);

  hist->SetTitle(legend_name.c_str());
  return hist;
}


#endif // StackedHistogram_cxx

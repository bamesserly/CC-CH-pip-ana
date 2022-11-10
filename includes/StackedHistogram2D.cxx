#ifndef StackedHistogram2D_cxx
#define StackedHistogram2D_cxx

#include "StackedHistogram.h"
#include "StackedHistogram2D.h"
#include "utilities.h"    // uniq

// CTOR -- default
template <typename T>
StackedHistogram2D<T>::StackedHistogram2D()
  : m_labelX(),
    m_labelY(),
    m_xlabelX(),
    m_xlabelY(),
    m_bins_arrayX(0),
    m_bins_arrayY(0),
    m_nhists(0),
    m_color_scheme(0),
    m_hist_map(),
    m_hist_array()
{}


// CTOR -- uniform binning
template <typename T>
StackedHistogram2D<T>::StackedHistogram2D(std::string labelX, std::string labelY,
		                          std::string xlabelX, std::string xlabelY,
                     			  int nbinsX, double xminX, double xmaxX,
                     			  int nbinsY, double xminY, double xmaxY,
                                          int nhists, int color_scheme)
  : m_labelX(labelX),
    m_labelY(labelY),
    m_xlabelX(xlabelX),
    m_xlabelY(xlabelY),
    m_bins_arrayX(MakeUniformBinArray(nbinsX, xminX, xmaxX)),
    m_bins_arrayY(MakeUniformBinArray(nbinsY, xminY, xmaxY)),
    m_nhists(nhists),
    m_color_scheme(color_scheme),
    m_hist_map(),
    m_hist_array()
{
  Initialize();
}


// CTOR -- variable binning
template <typename T>
StackedHistogram2D<T>::StackedHistogram2D(std::string labelX, std::string labelY,
                     			  std::string xlabelX, std::string xlabelY,
                 			  const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
                                          int nhists, int color_scheme)
  : m_labelX(labelX),
    m_labelY(labelY),
    m_xlabelX(xlabelX),
    m_xlabelY(xlabelY),
    m_bins_arrayX(GetSortedArray(bins_arrayX)),
    m_bins_arrayY(GetSortedArray(bins_arrayY)),
    m_nhists(nhists),
    m_color_scheme(color_scheme),
    m_hist_map(),
    m_hist_array()
{
  Initialize();
}


// Add component hists to stack
template <typename T>
void StackedHistogram2D<T>::Initialize() {
  for (int i=0; i != m_nhists; ++i){
    T type = static_cast<T>(i);
    PlotUtils::MnvH2D* component_hist = MakeStackComponentHist(type);
    m_hist_map[type] = component_hist;
    m_hist_array.Add(component_hist);
  }
}


// Make component hist
template <typename T>
PlotUtils::MnvH2D* StackedHistogram2D<T>::MakeStackComponentHist(const T type) const  {
  std::string legend_name = GetTruthClassification_LegendLabel(type);
  std::string short_name  = GetTruthClassification_Name(type);

  PlotUtils::MnvH2D* hist = nullptr;

  hist = new PlotUtils::MnvH2D( 
      uniq(), Form("%s_vs_%s_%s", m_labelX.c_str(), m_labelY.c_str(), short_name.c_str()), 
      NBinsX(), m_bins_arrayX.GetArray(), NBinsY(), m_bins_arrayY.GetArray());
  
//  SetHistColorScheme(hist, int(type), m_color_scheme);

  hist->SetTitle(legend_name.c_str());
  return hist;
}


#endif // StackedHistogram_cxx

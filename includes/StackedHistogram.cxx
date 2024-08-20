#ifndef StackedHistogram_cxx
#define StackedHistogram_cxx

#include "StackedHistogram.h"

#include "utilities.h"  // uniq

// CTOR -- default
template <typename T>
StackedHistogram<T>::StackedHistogram()
    : m_label(),
      m_xlabel(),
      m_bins_array(0),
      m_nhists(0),
      m_color_scheme(0),
      m_hist_map(),
      m_hist_array() {}

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
      m_hist_array() {
  Initialize();
}

// CTOR -- variable binning
template <typename T>
StackedHistogram<T>::StackedHistogram(std::string label, std::string xlabel,
                                      const TArrayD& bins_array, int nhists,
                                      int color_scheme)
    : m_label(label),
      m_xlabel(xlabel),
      m_bins_array(GetSortedArray(bins_array)),
      m_nhists(nhists),
      m_color_scheme(color_scheme),
      m_hist_map(),
      m_hist_array() {
  Initialize();
}

// COPY
template <typename T>
StackedHistogram<T>::StackedHistogram(const StackedHistogram& h)
    : m_label(h.m_label),
      m_xlabel(h.m_xlabel),
      m_bins_array(h.m_bins_array),
      m_nhists(h.m_nhists),
      m_color_scheme(h.m_color_scheme),
      m_hist_map(h.m_hist_map),
      m_hist_array(h.m_hist_array) {}

// COPY ASSIGNMENT OPERATOR
template <typename T>
StackedHistogram<T>& StackedHistogram<T>::operator=(const StackedHistogram& h) {
  if (this != &h) {
    m_label = h.m_label;
    m_xlabel = h.m_xlabel;
    m_bins_array = h.m_bins_array;
    m_nhists = h.m_nhists;
    m_color_scheme = h.m_color_scheme;
    m_hist_map = h.m_hist_map;
    m_hist_array = h.m_hist_array;
  }
  return *this;
}

// Add component hists to stack
template <typename T>
void StackedHistogram<T>::Initialize() {
  for (int i = 0; i != m_nhists; ++i) {
    T type = static_cast<T>(i);
    PlotUtils::MnvH1D* component_hist = MakeStackComponentHist(type);
    m_hist_map[type] = component_hist;
    m_hist_array.Add(component_hist);
  }
}

// Make component hist
template <typename T>
PlotUtils::MnvH1D* StackedHistogram<T>::MakeStackComponentHist(
    const T type) const {
  std::string legend_name = GetTruthClassification_LegendLabel(type);
  std::string short_name = GetTruthClassification_Name(type);

  PlotUtils::MnvH1D* hist = nullptr;

  hist = new PlotUtils::MnvH1D(
      uniq(), Form("%s_%s", m_label.c_str(), short_name.c_str()), NBins(),
      m_bins_array.GetArray());

  SetHistColorScheme(hist, int(type), m_color_scheme);

  hist->SetTitle(legend_name.c_str());
  return hist;
}

template <typename T>
void StackedHistogram<T>::LoadStackedFromFile(TFile& fin,
                                              UniverseMap& error_bands) {
  for (int i = 0; i != m_nhists; ++i) {
    T type = static_cast<T>(i);
    std::string legend_name = GetTruthClassification_LegendLabel(type);
    std::string short_name = GetTruthClassification_Name(type);
    std::string name = Form("%s_%s", m_label.c_str(), short_name.c_str());
    PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)fin.Get(name.c_str());

    // Hist not found, quit
    if (hist == 0) {
      std::cout << name << " hist doesn't exist. skipping...\n";
      return;
    }

    // Binning Checks
    TArrayD bins_array = *(hist->GetXaxis()->GetXbins());
    // Source histo has uniform binning
    if (bins_array.GetSize() == 0) {
      bins_array.Reset();
      bins_array = MakeUniformBinArray(hist->GetXaxis()->GetNbins(),
                                       hist->GetXaxis()->GetXmin(),
                                       hist->GetXaxis()->GetXmax());
      hist = dynamic_cast<PlotUtils::MnvH1D*>(hist->Rebin(
          bins_array.GetSize() - 1, hist->GetName(), bins_array.GetArray()));
    }

    // Compare source and destination binnings
    for (int i = 0; i < NBins(); ++i) {
      if (m_bins_array[i] != bins_array[i]) {
        std::cout << "WARNING! Binning mismatch for " << m_label << "\n";
        std::cout << "Aligning output binning to match source binning\n";
        m_bins_array = bins_array;
        break;
      }
    }

    SetHistColorScheme(hist, int(type), m_color_scheme);
    m_hist_map[type] = hist;
    m_hist_array.Add(hist);
  }  // loop over components
}

#endif  // StackedHistogram_cxx

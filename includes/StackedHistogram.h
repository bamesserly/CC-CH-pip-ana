#ifndef StackedHistogram_h
#define StackedHistogram_h

#include "TArrayD.h"
#include "Constants.h" // CCNuPionIncPlotting, SetHistColorScheme
#include "PlotUtils/MnvH1D.h"
#include "Binning.h" // MakeUniformBinArray

template<typename T>
class StackedHistogram {
  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    StackedHistogram();

    StackedHistogram(std::string label, std::string xlabel, 
                     int nbins, double xmin, double xmax,
                     int nhists, int color_scheme = 0);

    StackedHistogram(std::string label, std::string xlabel, 
                     const TArrayD& bins_array,
                     int nhists, int color_scheme = 0);


    //==========================================================================
    // Data Members
    //==========================================================================
    std::string m_label;
    std::string m_xlabel;
    TArrayD m_bins_array;
    int m_color_scheme;
    int m_nhists;

    std::map<T, PlotUtils::MnvH1D*> m_hist_map;
    TObjArray m_hist_array;


    //==========================================================================
    // Functions
    //==========================================================================
    int NBins()   const { return m_bins_array.GetSize()-1; }
    double XMin() const { return m_bins_array[0]; }
    double XMax() const { return m_bins_array[NBins()]; }
    void Initialize();
    PlotUtils::MnvH1D* MakeStackComponentHist(const T type) const;
};


// Template member functions need to be available in the header.
#include "StackedHistogram.cxx"


#endif

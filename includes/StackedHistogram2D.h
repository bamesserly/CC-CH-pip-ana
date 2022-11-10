#ifndef StackedHistogram2D_h
#define StackedHistogram2D_h

#include "TArrayD.h"
#include "Constants.h" // CCNuPionIncPlotting, SetHistColorScheme
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "Binning.h" // MakeUniformBinArray

template<typename T>
class StackedHistogram2D {
  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    StackedHistogram2D();

    StackedHistogram2D(std::string labelX, std::string labelY, 
    		     std::string xlabelX, std::string xlabelY,
                     int nbinsX, double xminX, double xmaxX,
                     int nbinsY, double xminY, double xmaxY,
                     int nhists, int color_scheme = 0);

    StackedHistogram2D(std::string labelX, std::string labelY,
                     std::string xlabelX, std::string xlabelY,
                     const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
                     int nhists, int color_scheme = 0);


    //==========================================================================
    // Data Members
    //==========================================================================
    std::string m_labelX;
    std::string m_labelY;
    std::string m_xlabelX;
    std::string m_xlabelY;
    TArrayD m_bins_arrayX;
    TArrayD m_bins_arrayY;
    int m_color_scheme;
    int m_nhists;

    std::map<T, PlotUtils::MnvH2D*> m_hist_map;
    TObjArray m_hist_array;


    //==========================================================================
    // Functions
    //==========================================================================
    int NBinsX()   const { return m_bins_arrayX.GetSize()-1; }
    int NBinsY()   const { return m_bins_arrayY.GetSize()-1; }
    double XMinX() const { return m_bins_arrayX[0]; }
    double XMinY() const { return m_bins_arrayY[0]; }
    double XMaxX() const { return m_bins_arrayX[NBinsX()]; }
    double XMaxY() const { return m_bins_arrayY[NBinsY()]; }
    void Initialize();
    PlotUtils::MnvH2D* MakeStackComponentHist(const T type) const;
};


// Template member functions need to be available in the header.
#include "StackedHistogram.cxx"
#include "StackedHistogram2D.cxx"


#endif

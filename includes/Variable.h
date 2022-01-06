#ifndef Variable_h
#define Variable_h

#include <functional>
#include "TArrayD.h"
#include "CVUniverse.h"
#include "Histograms.h"


class Variable {
  protected:
    typedef std::function<double(const CVUniverse&)> PointerToCVUniverseFunction;
    PointerToCVUniverseFunction m_pointer_to_GetValue;

  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    Variable();

    Variable(const std::string label, const std::string xaxis, 
             const std::string units,
             const int nbins, const double xmin, const double xmax,
             PointerToCVUniverseFunction p = &CVUniverse::GetDummyVar,
             const bool is_true = false);

    Variable(const std::string label, const std::string xaxis, 
             const std::string units,
             const TArrayD& bins_array,
             PointerToCVUniverseFunction p = &CVUniverse::GetDummyVar,
             const bool is_true = false);

    //==========================================================================
    // Data members
    //==========================================================================
    std::string m_label;
    std::string m_units;
    Histograms  m_hists;
    bool m_is_true;
    // also, a pointer to CV universe function for Getting value (private)

    //==========================================================================
    // Functions
    //==========================================================================
    // Access
    std::string Name() const { return m_label; }
    int NBins() const { return m_hists.NBins(); }
    int XMin() const  { return m_hists.XMin(); }
    int XMax() const  { return m_hists.XMax(); }

    // Get the variable's value
    virtual double GetValue (const CVUniverse& universe, const int hadron_ID = -1) const;

    // Histogram Initialization
    template<typename T>
    void InitializeAllHists(T systematic_univs, T systematic_univs_truth);
    template<typename T>
    void InitializeSidebandHists(T systematic_univs);
    void InitializeStackedHists();
    void InitializeDataHists();

    // Write and Load MC Hists to/from a file
    void WriteMCHists(TFile& fout) const;
    void LoadDataHistsFromFile(TFile& fin);
    void LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands);


    // Histogram Access
    template <typename T>
    PlotUtils::MnvH1D* GetStackComponentHist(T type) const;

    template <typename T>
    TObjArray GetStackArray(T type) const ;

    // Get Histograms from File
    //void GetMCHists(TFile& fin);
};


#include "Variable.cxx"


#endif // Variable_h

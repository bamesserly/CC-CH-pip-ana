#ifndef Variable_h
#define Variable_h

#include <functional>
#include "TArrayD.h"
#include "Histograms.h"

template<class T>
class Variable {
  protected:
    typedef std::function<double(const T&, const std::string&)> PointerToTFunction;
    PointerToTFunction m_pointer_to_GetValue;

  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    Variable();

    Variable(const std::string label, const std::string xaxis, 
             const std::string units,
             const int nbins, const double xmin, const double xmax,
             PointerToTFunction p = &T::GetDummyVar,
             const bool is_true = false);

    Variable(const std::string label, const std::string xaxis, 
             const std::string units,
             const TArrayD& bins_array,
             PointerToTFunction p = &T::GetDummyVar,
             const bool is_true = false);

    //==========================================================================
    // Data members
    //==========================================================================
    // [a pointer to CV universe function for getting value (private)]
    std::string m_label;
    std::string m_units;
    Histograms  m_hists;
    bool m_is_true;

    //==========================================================================
    // Functions
    //==========================================================================
    // Access
    std::string Name() const { return m_label; }
    int NBins() const { return m_hists.NBins(); }
    int XMin() const  { return m_hists.XMin(); }
    int XMax() const  { return m_hists.XMax(); }

    // Get the variable's value
    virtual double GetValue (const T& t, const std::string& s) const;

    // Histogram Initialization
    template<typename U>
    void InitializeAllHists(U systematic_univs, U systematic_univs_truth);
    template<typename U>
    void InitializeSidebandHists(U systematic_univs);
    void InitializeStackedHists();
    void InitializeDataHists();

    // Write and Load MC Hists to/from a file
    void WriteMCHists(TFile& fout) const;
    void LoadDataHistsFromFile(TFile& fin);
    void LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands);


    // Histogram Access
    template <typename U>
    PlotUtils::MnvH1D* GetStackComponentHist(U type) const;

    template <typename U>
    TObjArray GetStackArray(U type) const ;

    // Get Histograms from File
    //void GetMCHists(TFile& fin);
};

#include "Variable.cxx"

#endif // Variable_h

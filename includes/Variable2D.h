#ifndef Variable2D_h
#define Variable2D_h

#include <functional>
#include "TArrayD.h"
#include "CVUniverse.h"
#include "Histograms.h"
#include "Histograms2D.h"
#include "Variable.h"
#include "Constants.h"
#include "Binning.h"
class Variable;

class Variable2D {
  protected:
    typedef std::function<double(const CVUniverse&)> PointerToCVUniverseFunction;
    PointerToCVUniverseFunction m_pointer_to_GetValueX;
    PointerToCVUniverseFunction m_pointer_to_GetValueY;

  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    Variable2D();

    Variable2D(const std::string labelX, const std::string labelY, 
	       const std::string xaxisX, const std::string xaxisY,
	       const std::string unitsX, const std::string unitsY,  
               const int nbinsX, const double xminX, const double xmaxX,
               const int nbinsY, const double xminY, const double xmaxY,
 	       PointerToCVUniverseFunction px = &CVUniverse::GetDummyVar,
	       PointerToCVUniverseFunction py = &CVUniverse::GetDummyVar, 
	       const bool is_true = false);

    Variable2D(const std::string labelX, const std::string labelY,
               const std::string xaxisX, const std::string xaxisY,
	       const std::string unitsX, const std::string unitsY,
	       const TArrayD& bins_arrayX, const TArrayD& bins_arrayY,
   	       PointerToCVUniverseFunction px = &CVUniverse::GetDummyVar,
	       PointerToCVUniverseFunction py = &CVUniverse::GetDummyVar,
	       const bool is_true = false);

    Variable2D(const Variable*,
	       const Variable*);

    Variable2D(const std::string name,
	       const Variable*,
               const Variable*);

    //==========================================================================
    // Data members
    //==========================================================================
    std::string m_labelX;
    std::string m_unitsX;
    Histograms2D  m_hists2D;
    std::string m_labelY;
    std::string m_unitsY;
    MResponse m_response;
    bool m_is_true;
    // also, a pointer to CV universe function for Getting value (private)

    //==========================================================================
    // Functions
    //==========================================================================
    // Access
    std::string NameX() const { return m_labelX; }
    std::string NameY() const { return m_labelY; }
    int NBinsX() const { return m_hists2D.NBinsX(); }
    int NBinsY() const { return m_hists2D.NBinsY(); }
    int XMinX() const  { return m_hists2D.XMinX(); }
    int XMinY() const  { return m_hists2D.XMinY(); }
    int XMaxX() const  { return m_hists2D.XMaxX(); }
    int XMaxY() const  { return m_hists2D.XMaxY(); }

    // Get the variable's value
    virtual double GetValueX (const CVUniverse& universe, const int hadron_ID = -1) const;
    virtual double GetValueY (const CVUniverse& universe, const int hadron_ID = -1) const;

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
    PlotUtils::MnvH2D* GetStackComponentHist(T type) const;

    template <typename T>
    TObjArray GetStackArray(T type) const ;

    // Get Histograms from File
    //void GetMCHists(TFile& fin);
};


#include "Variable2D.cxx"


#endif // Variable_h

#ifndef Variable_h
#define Variable_h

#include <functional>

#include "CVUniverse.h"
#include "Histograms.h"
#include "TArrayD.h"
#include "CCPiEvent.h"

class VariableBase {
 public:
  //==========================================================================
  // Constructors
  //==========================================================================
  VariableBase();

  VariableBase(const std::string label, const std::string xaxis,
               const std::string units, const int nbins, const double xmin,
               const double xmax, const bool is_true = false);

  VariableBase(const std::string label, const std::string xaxis,
               const std::string units, const TArrayD& bins_array,
               const bool is_true = false);

  //==========================================================================
  // Data members
  //==========================================================================
  // [a pointer to CV universe function for getting value (private)]
  std::string m_label;
  std::string m_units;
  Histograms m_hists;
  bool m_is_true;

  //==========================================================================
  // Functions
  //==========================================================================
  // Access
  std::string Name() const { return m_label; }
  int NBins() const { return m_hists.NBins(); }
  int XMin() const { return m_hists.XMin(); }
  int XMax() const { return m_hists.XMax(); }

  // Get the variable's value
  virtual double GetValue(const CCPiEvent&) const = 0;

  // Histogram Initialization
  template <class U>
  void InitializeAllHists(U systematic_univs, U systematic_univs_truth);
  template <class U>
  void InitializeSidebandHists(U systematic_univs);
  void InitializeStackedHists();
  void InitializeDataHists();

  // Write and Load MC Hists to/from a file
  void WriteMCHists(TFile& fout) const;
  void LoadDataHistsFromFile(TFile& fin);
  void LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands);

  // Histogram Access
  template <class U>
  PlotUtils::MnvH1D* GetStackComponentHist(U type) const;

  template <class U>
  TObjArray GetStackArray(U type) const;

  // Get Histograms from File
  // void GetMCHists(TFile& fin);
};

class CVUVariable : public VariableBase {
 protected:
  typedef std::function<double(const CVUniverse&)> PointerToCVUFunction;
  PointerToCVUFunction m_pointer_to_GetValue;

 public:
  CVUVariable();

  CVUVariable(const std::string label, const std::string xaxis,
              const std::string units, const int nbins, const double xmin,
              const double xmax,
              PointerToCVUFunction p = &CVUniverse::GetDummyVar,
              const bool is_true = false);

  CVUVariable(const std::string label, const std::string xaxis,
              const std::string units, const TArrayD& bins_array,
              PointerToCVUFunction p = &CVUniverse::GetDummyVar,
              const bool is_true = false);

  double GetValue(const CCPiEvent& e) const;
};

class EventVariable : public VariableBase {
 protected:
  typedef std::function<double(const CCPiEvent&)> PointerToEventFunction;
  PointerToEventFunction m_pointer_to_GetValue;

 public:
  EventVariable();

  EventVariable(const std::string label, const std::string xaxis,
                const std::string units, const int nbins, const double xmin,
                const double xmax,
                PointerToEventFunction p = &CCPiEvent::GetDummyVar,
                const bool is_true = false);

  EventVariable(const std::string label, const std::string xaxis,
                const std::string units, const TArrayD& bins_array,
                PointerToEventFunction p = &CCPiEvent::GetDummyVar,
                const bool is_true = false);

  double GetValue(const CCPiEvent& e) const;
};

#include "Variable.cxx"

#endif  // Variable_h

//==============================================================================
// Container class that holds a lot of histograms for a specific variable
// Also knows how to initialize those histograms
//==============================================================================
#ifndef Histograms_h
#define Histograms_h

#include "Binning.h"  // MakeUniformBinArray
#include "CVUniverse.h"
#include "Constants.h"  // typedefs MH1D, CVHW, CVH2DW
#include "StackedHistogram.h"
#include "TArrayD.h"
#include "TFile.h"
#include "TruthMatching.h"
#include "utilities.h"  // uniq

class Histograms {
 public:
  //==========================================================================
  // Constructors
  //==========================================================================
  Histograms();

  Histograms(const std::string label, const std::string xlabel, const int nbins,
             const double xmin, const double xmax);

  Histograms(const std::string label, const std::string xlabel,
             const TArrayD& bins_array);

  // Histograms(const Histograms&);

  //==========================================================================
  // Data Members
  //==========================================================================
  // Basic Data Members
  std::string m_label;
  std::string m_xlabel;
  TArrayD m_bins_array;

  // Cross Section Pipeline
  CVHW m_bg;
  CVHW m_bg_hiW;
  CVHW m_bg_loW;
  CVHW m_bg_midW;
  MH1D* m_bg_subbed_data;
  MH1D* m_cross_section;
  CVHW m_effden;
  MH1D* m_efficiency;
  CVHW m_effnum;
  CVH2DW m_migration;
  MH1D* m_selection_data;
  MH1D* m_selection_data_tracked;
  MH1D* m_selection_data_untracked;
  MH1D* m_selection_data_mixed;
  CVHW m_selection_mc;
  CVHW m_selection_mc_tracked;
  CVHW m_selection_mc_untracked;
  CVHW m_selection_mc_mixed;
  CVHW m_selection_mc_no_tpi_weight;
  CVHW m_selection_mc_tracked_no_tpi_weight;
  CVHW m_selection_mc_untracked_no_tpi_weight;
  CVHW m_selection_mc_mixed_no_tpi_weight;
  MH1D* m_tuned_bg;
  MH1D* m_unfolded;
  MH1D* m_wsideband_data;
  MH1D* m_noWcut_data;
  MH1D* m_wsidebandfit_data;
  CVHW m_wsidebandfit_hiW;
  CVHW m_wsidebandfit_loW;
  CVHW m_wsidebandfit_midW;
  CVHW m_wsidebandfit_sig;
  MH1D* m_noWcut;

  // Stacked Histograms
  StackedHistogram<ChannelType> m_stacked_channel;
  StackedHistogram<CoherentType> m_stacked_coherent;
  StackedHistogram<FSParticleType> m_stacked_fspart;
  StackedHistogram<HadronType> m_stacked_hadron;
  StackedHistogram<MesonBackgroundType> m_stacked_mesonbg;
  StackedHistogram<NPi0Type> m_stacked_npi0;
  StackedHistogram<NPionsType> m_stacked_npi;
  StackedHistogram<NPipType> m_stacked_npip;
  StackedHistogram<SignalBackgroundType> m_stacked_sigbg;
  StackedHistogram<WType> m_stacked_w;
  StackedHistogram<WBackgroundType> m_stacked_wbg;
  StackedHistogram<WSidebandType> m_stacked_wsideband;
  StackedHistogram<PionRecoType> m_stacked_pionreco;

  //==========================================================================
  // Functions
  //==========================================================================
  int NBins() const { return m_bins_array.GetSize() - 1; }
  double XMin() const { return m_bins_array[0]; }
  double XMax() const { return m_bins_array[NBins()]; }
  void PrintBinning() const {
    for (int i = 0; i <= NBins(); ++i) std::cout << m_bins_array[i] << " ";
    std::cout << "\n";
  }

  // Histogram Initialization
  template <typename T>
  void InitializeAllHists(T systematic_univs, T systematic_univs_truth);
  template <typename T>
  void InitializeSelectionHists(T systematic_univs, T systematic_univs_truth);
  template <typename T>
  void InitializeSidebandHists(T systematic_univs);
  template <typename T>
  void InitializeMigrationHist(T systematic_univs);
  void InitializeDataHists();
  void InitializeStackedHists();

  // Stack Map Access
  std::map<WType, MH1D*> GetStackMap(WType type) const;
  std::map<SignalBackgroundType, MH1D*> GetStackMap(
      SignalBackgroundType type) const;
  std::map<WBackgroundType, MH1D*> GetStackMap(WBackgroundType type) const;
  std::map<MesonBackgroundType, MH1D*> GetStackMap(
      MesonBackgroundType type) const;
  std::map<HadronType, MH1D*> GetStackMap(HadronType type) const;
  std::map<FSParticleType, MH1D*> GetStackMap(FSParticleType type) const;
  std::map<ChannelType, MH1D*> GetStackMap(ChannelType type) const;
  std::map<NPionsType, MH1D*> GetStackMap(NPionsType type) const;
  std::map<NPi0Type, MH1D*> GetStackMap(NPi0Type type) const;
  std::map<NPipType, MH1D*> GetStackMap(NPipType type) const;
  std::map<WSidebandType, MH1D*> GetStackMap(WSidebandType type) const;
  std::map<CoherentType, MH1D*> GetStackMap(CoherentType type) const;
  std::map<PionRecoType, MH1D*> GetStackMap(PionRecoType type) const;

  // Load MC hists from file
  void LoadDataHistsFromFile(TFile& fin);
  void LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands, bool is_true);
  CVHW LoadHWFromFile(TFile& fin, UniverseMap& error_bands, std::string name);
  CVH2DW LoadH2DWFromFile(TFile& fin, UniverseMap& error_bands,
                          std::string name);
};

// Template member functions need to be available in the header.
#include "Histograms.cxx"

#endif  // Histograms_h

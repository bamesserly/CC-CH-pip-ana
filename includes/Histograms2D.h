//==============================================================================
// Container class that holds a lot of histograms for a specific variable
// Also knows how to initialize those histograms
//==============================================================================
#ifndef Histograms2D_h
#define Histograms2D_h

#include "TFile.h"
#include "TArrayD.h"
#include "utilities.h" // uniq
#include "StackedHistogram.h"
#include "StackedHistogram2D.h"
#include "TruthMatching.h"
#include "CVUniverse.h"
#include "Constants.h" // typedefs MH1D, MH2D, CVHW, CVH2DW
#include "Binning.h" // MakeUniformBinArray

class Histograms2D {
  public:
    //==========================================================================
    // Constructors
    //==========================================================================
    Histograms2D();

    Histograms2D(const std::string labelX, const std::string labelY, 
	         const std::string xlabelX, const std::string xlabelY,
                 const int nbinsX, const int nbinsY, const double xminX,
	         const double xminY, const double xmaxX, const double xmaxY);

    Histograms2D(const std::string labelX, const std::string labelY,
                 const std::string xlabelX, const std::string xlabelY,
                 const TArrayD& bins_arrayX, const TArrayD& bins_arrayY);


    //==========================================================================
    // Data Members
    //==========================================================================
    // Basic Data Members
    std::string m_labelX;
    std::string m_labelY;
    std::string m_xlabelX;
    std::string m_xlabelY;
    TArrayD m_bins_arrayX;
    TArrayD m_bins_arrayY;

    // Histograms2D -- Event Selection
    MH2D* m_selection_data;  // DATA after cuts
    CVH2DW m_selection_mc;     // MC   after cuts, with systematics
    CVH2DW m_bg;               // BACKGROUND
    CVH2DW m_bg_loW;           // BACKGROUND
    CVH2DW m_bg_midW;          // BACKGROUND
    CVH2DW m_bg_hiW;           // BACKGROUND
    CVH2DW m_effnum;           // EFF num
    CVH2DW m_effden;           // EFF den

    // Histograms2D -- Later-stage Cross Section Calculation
    MH2D* m_tuned_bg;
    MH2D* m_bg_subbed_data;
    MH2D* m_efficiency;
    MH2D* m_unfolded;
    MH2D* m_cross_section;

    // Migration
    CVH2DW m_migration;
    CVH2DW m_migration_reco;
    CVH2DW m_migration_true;
    MResponse m_response;
    MH2D* m_responseMigration;
    MH2D* m_responseTrue;
    MH2D* m_responseReco;
 
    // Histograms2D -- Sidebands
    MH2D* m_wsidebandfit_data;
    CVH2DW  m_wsidebandfit_sig;
    CVH2DW  m_wsidebandfit_loW;
    CVH2DW  m_wsidebandfit_midW;
    CVH2DW  m_wsidebandfit_hiW;

    MH2D* m_wsideband_data;

    // Stacked Histograms2D -- Cut Studies
    StackedHistogram2D <WType>                m_stacked_w;
    StackedHistogram2D <SignalBackgroundType> m_stacked_sigbg;
    StackedHistogram2D <WBackgroundType>      m_stacked_wbg;
    StackedHistogram2D <MesonBackgroundType>  m_stacked_mesonbg;
    StackedHistogram2D <HadronType>           m_stacked_hadron;
    StackedHistogram2D <FSParticleType>       m_stacked_fspart;
    StackedHistogram2D <ChannelType>          m_stacked_channel;
    StackedHistogram2D <NPionsType>           m_stacked_npi;
    StackedHistogram2D <NPi0Type>             m_stacked_npi0;
    StackedHistogram2D <NPipType>             m_stacked_npip;
    StackedHistogram2D <CoherentType>         m_stacked_coherent;

    StackedHistogram2D <WSidebandType>        m_stacked_wsideband;


    //==========================================================================
    // Functions
    //==========================================================================
    int NBinsX()   const { return m_bins_arrayX.GetSize()-1; }
    int NBinsY()   const { return m_bins_arrayY.GetSize()-1; }
    double XMinX() const { return m_bins_arrayX[0]; }
    double XMinY() const { return m_bins_arrayY[0]; }
    double XMaxX() const { return m_bins_arrayX[NBinsX()]; }
    double XMaxY() const { return m_bins_arrayY[NBinsY()]; }
    void PrintBinningX() const {
      for(int i = 0; i <= NBinsX(); ++i) std::cout << m_bins_arrayX[i] << " ";
      std::cout << "\n";
    }
    void PrintBinningY() const {
      for(int i = 0; i <= NBinsY(); ++i) std::cout << m_bins_arrayY[i] << " ";
      std::cout << "\n";
    }


    // Histogram Initialization
    template<typename T> void InitializeAllHists(T systematic_univs,       T systematic_univs_truth);
    template<typename T> void InitializeSelectionHists(T systematic_univs, T systematic_univs_truth);
    template<typename T> void InitializeSidebandHists(T systematic_univs);
    template<typename T> void InitializeMigrationHist(T systematic_univs);
    void InitializeDataHists();
    void InitializeStackedHists();

    // Stack Map Access
    std::map<WType,                MH2D*> GetStackMap (WType type)                const;
    std::map<SignalBackgroundType, MH2D*> GetStackMap (SignalBackgroundType type) const;
    std::map<WBackgroundType,      MH2D*> GetStackMap (WBackgroundType type)      const;
    std::map<MesonBackgroundType,  MH2D*> GetStackMap (MesonBackgroundType type)  const;
    std::map<HadronType,           MH2D*> GetStackMap (HadronType type)           const;
    std::map<FSParticleType,       MH2D*> GetStackMap (FSParticleType type)       const;
    std::map<ChannelType,          MH2D*> GetStackMap (ChannelType type)          const;
    std::map<NPionsType,           MH2D*> GetStackMap (NPionsType type)           const;
    std::map<NPi0Type,             MH2D*> GetStackMap (NPi0Type type)             const;
    std::map<NPipType,             MH2D*> GetStackMap (NPipType type)             const;
    std::map<WSidebandType,        MH2D*> GetStackMap (WSidebandType type)        const;
    std::map<CoherentType,         MH2D*> GetStackMap (CoherentType type)         const;

    // Load MC hists from file
    void LoadDataHistsFromFile(TFile& fin);
    void LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands, bool is_true);
//  CVHW LoadHWFromFile(TFile& fin, UniverseMap& error_bands, std::string name);
    CVH2DW LoadH2DWFromFile(TFile& fin, UniverseMap& error_bands, std::string name);
};


// Template member functions need to be available in the header.
#include "Histograms2D.cxx"


#endif // Histograms_h

#ifndef Histograms2D_cxx
#define Histograms2D_cxx

#include <algorithm>
#include "Histograms.h"
#include "Histograms2D.h"
// CTOR -- default
Histograms2D::Histograms2D() 
  : m_labelX(),          m_xlabelX(),
    m_labelY(),          m_xlabelY(),
    m_bins_arrayX(0),
    m_bins_arrayY(0),
    m_selection_data(),      m_selection_mc(),        m_bg(),
    m_bg_loW(),              m_bg_midW(),             m_bg_hiW(),
    m_effnum(),              m_effden(),
    m_stacked_w(),           m_stacked_hadron(),
    m_stacked_fspart(),      m_stacked_channel(),
    m_stacked_npi(),         m_stacked_npi0(),        m_stacked_npip(),
    m_stacked_sigbg(),       m_stacked_wbg(),         m_stacked_mesonbg(),
    m_stacked_coherent(),
    m_wsidebandfit_data(),   m_wsidebandfit_sig(),
    m_wsidebandfit_loW(),    m_wsidebandfit_midW(),   m_wsidebandfit_hiW(),
    m_stacked_wsideband(),   m_wsideband_data(),
    m_tuned_bg(),            m_bg_subbed_data(),
    m_migration(),           m_unfolded(),            m_cross_section(),
    m_migration_reco(),      m_migration_true(),      m_response(), 
    m_responseTrue(),           m_responseReco(),     m_responseMigration()
{}


// CTOR -- uniform binning
Histograms2D::Histograms2D(const std::string labelX, const std::string labelY,
                 	   const std::string xlabelX, const std::string xlabelY,
                 	   const int nbinsX, const int nbinsY, const double xminX,
                 	   const double xminY, const double xmaxX, const double xmaxY)
  : m_labelX(labelX),          m_xlabelX(xlabelX), 
    m_labelY(labelY),          m_xlabelY(xlabelY),
    m_bins_arrayX(MakeUniformBinArray(nbinsX, xminX, xmaxX)),
    m_bins_arrayY(MakeUniformBinArray(nbinsY, xminY, xmaxY)),
    m_selection_data(),      m_selection_mc(),        m_bg(),
    m_bg_loW(),              m_bg_midW(),             m_bg_hiW(),
    m_effnum(),              m_effden(),
    m_stacked_w(),           m_stacked_hadron(),
    m_stacked_fspart(),      m_stacked_channel(),
    m_stacked_npi(),         m_stacked_npi0(),        m_stacked_npip(),
    m_stacked_sigbg(),       m_stacked_wbg(),         m_stacked_mesonbg(),
    m_stacked_coherent(),
    m_wsidebandfit_data(),   m_wsidebandfit_sig(),
    m_wsidebandfit_loW(),    m_wsidebandfit_midW(),   m_wsidebandfit_hiW(),
    m_stacked_wsideband(),   m_wsideband_data(),
    m_tuned_bg(),            m_bg_subbed_data(),
    m_migration(),           m_unfolded(),            m_cross_section(),
    m_migration_reco(),      m_migration_true(),      m_response(), 
    m_responseTrue(),        m_responseReco(),        m_responseMigration()
{}


// CTOR -- variable binning
Histograms2D::Histograms2D(const std::string labelX, const std::string labelY,
                           const std::string xlabelX, const std::string xlabelY,
                           const TArrayD& bins_arrayX, const TArrayD& bins_arrayY)
  : m_labelX(labelX),          m_xlabelX(xlabelX),
    m_labelY(labelY),          m_xlabelY(xlabelY),
    m_bins_arrayX(GetSortedArray(bins_arrayX)),
    m_bins_arrayY(GetSortedArray(bins_arrayY)),
    m_selection_data(),         m_selection_mc(),        m_bg(),
    m_bg_loW(),                 m_bg_midW(),             m_bg_hiW(),
    m_effnum(),                 m_effden(),
    m_stacked_w(),              m_stacked_hadron(),
    m_stacked_fspart(),         m_stacked_channel(),
    m_stacked_npi(),            m_stacked_npi0(),        m_stacked_npip(),
    m_stacked_sigbg(),          m_stacked_wbg(),         m_stacked_mesonbg(),
    m_stacked_coherent(),
    m_wsidebandfit_data(),      m_wsidebandfit_sig(),
    m_wsidebandfit_loW(),       m_wsidebandfit_midW(),   m_wsidebandfit_hiW(),
    m_stacked_wsideband(),      m_wsideband_data(),
    m_tuned_bg(),               m_bg_subbed_data(),
    m_migration(),              m_unfolded(),            m_cross_section(),
    m_migration_reco(),         m_migration_true(),      m_response(),
    m_responseTrue(),           m_responseReco(),        m_responseMigration()
{}


// Load hists from file
  void Histograms2D::LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands,
                                       bool is_true) {
    m_selection_mc = LoadH2DWFromFile(fin, error_bands, "selection_mc");

    m_bg           = LoadH2DWFromFile(fin, error_bands, "bg");
    m_bg_loW       = LoadH2DWFromFile(fin, error_bands, "bg_loW");
    m_bg_midW      = LoadH2DWFromFile(fin, error_bands, "bg_midW");
    m_bg_hiW       = LoadH2DWFromFile(fin, error_bands, "bg_hiW");
    m_tuned_bg     = (PlotUtils::MnvH2D*)fin.Get(Form("tuned_bg_%s_vs_%s",
                                                      m_labelX.c_str(), m_labelY.c_str()));

    m_effnum       = LoadH2DWFromFile(fin, error_bands, "effnum");
    m_effden       = LoadH2DWFromFile(fin, error_bands, "effden");
   
    if (m_labelX == sidebands::kFitVarString) {
      m_wsidebandfit_sig   = LoadH2DWFromFile(fin, error_bands, "wsidebandfit_sig");
      m_wsidebandfit_loW   = LoadH2DWFromFile(fin, error_bands, "wsidebandfit_loW");
      m_wsidebandfit_midW  = LoadH2DWFromFile(fin, error_bands, "wsidebandfit_midW");
      m_wsidebandfit_hiW   = LoadH2DWFromFile(fin, error_bands, "wsidebandfit_hiW");
    }
   
    if (!is_true && m_labelX != sidebands::kFitVarString) {
      m_migration = LoadH2DWFromFile(fin, error_bands, "migration");
      m_migration_reco = LoadH2DWFromFile(fin, error_bands, "migration_reco");
      m_migration_true = LoadH2DWFromFile(fin, error_bands, "migration_truth");
//      m_responseMigration = (PlotUtils::MnvH2D*)fin.Get(Form("Migration2d_%s_vs_%s_migration",m_labelX.c_str(), m_labelY.c_str()));
//      m_responseMigration = LoadH2DWFromFile(fin, error_bands, "Migration2d")
      assert(m_migration.hist);
      assert(m_migration_reco.hist);
      assert(m_migration_true.hist);
    }
  }


/*  CVHW Histograms2D::LoadHWFromFile(TFile& fin, UniverseMap& error_bands,
                                  std::string name) {
    const bool do_erase_bands = false;
    PlotUtils::MnvH2D* hist = 
        (PlotUtils::MnvH2D*)fin.Get(Form("%s_%s_vs_%s", name.c_str(), m_labelX.c_str(), m_labelY.c_str()));

    if (hist == 0) {
      std::cout << name << " hist doesn't exist. skipping...\n";
      return CVH2DW();
    }
    else {
      TArrayD bins_arrayX = *(hist->GetXaxis()->GetXbins());
      TArrayD bins_arrayY = *(hist->GetYaxis()->GetXbins());
      // Source histo has uniform binning
      if (bins_arrayX.GetSize() == 0 && bins_arrayY.GetSize() == 0) {
        bins_arrayX.Reset();
        bins_arrayY.Reset();
        bins_arrayX = MakeUniformBinArray(hist->GetXaxis()->GetNbins(),
                                         hist->GetXaxis()->GetXmin(),
                                         hist->GetXaxis()->GetXmax());
        bins_arrayY = MakeUniformBinArray(hist->GetYaxis()->GetNbins(),
                                         hist->GetYaxis()->GetXmin(),
                                         hist->GetYaxis()->GetXmax());
        hist = dynamic_cast<PlotUtils::MnvH2D*>( hist->Rebin(bins_array.GetSize()-1,
                                                             hist->GetName(),
                                                             bins_array.GetArray()) );
      }

      // Compare source and destination binnings
      for(int i = 0; i < NBins(); ++i) {
        if (m_bins_array[i] != bins_array[i]) {
          std::cout << "WARNING! Binning mismatch for " << m_label << "\n";
          std::cout << "Aligning output binning to match source binning\n";
          m_bins_array = bins_array;
          break;
        }
      }

      return CVHW(hist, error_bands, do_erase_bands);
    }
  }*/


  CVH2DW Histograms2D::LoadH2DWFromFile(TFile& fin, UniverseMap& error_bands,
                                      std::string name) {
    const bool do_erase_bands = false;
    PlotUtils::MnvH2D* hist = 
        (PlotUtils::MnvH2D*)fin.Get(Form("%s_%s_vs_%s", name.c_str(), m_labelX.c_str(), m_labelY.c_str()));
    assert(hist);
    return CVH2DW(hist, error_bands, do_erase_bands);
  }

  void Histograms2D::LoadDataHistsFromFile(TFile& fin) {
    m_selection_data = (PlotUtils::MnvH2D*)fin.Get(Form("selection_data_%s_vs_%s", m_labelX.c_str(), m_labelY.c_str()));
    m_bg_subbed_data = (PlotUtils::MnvH2D*)fin.Get(Form("bg_subbed_data_%s_vs_%s", m_labelX.c_str(), m_labelY.c_str()));
    m_unfolded       = (PlotUtils::MnvH2D*)fin.Get(Form("unfolded_%s_vs_%s", m_labelX.c_str(), m_labelY.c_str()));
    m_cross_section  = (PlotUtils::MnvH2D*)fin.Get(Form("D_cross_section_%s_vs_%s", m_labelX.c_str(), m_labelY.c_str()));
  }

 
// Initialize Hists
  template<typename T>
  void Histograms2D::InitializeAllHists(T systematic_univs,
                                      T systematic_univs_truth) {
    // Event Section Analysis
    InitializeSelectionHists(systematic_univs, systematic_univs_truth);

    // Migration Matrix
    InitializeMigrationHist(systematic_univs);

    // Sidebands
//    InitializeSidebandHists(systematic_univs);

    // Data
    InitializeDataHists();

    // Event Selection Stacked
//    InitializeStackedHists();
  }


  template<typename T>
  void Histograms2D::InitializeSelectionHists(T systematic_univs,
                                            T systematic_univs_truth) {
    const Double_t* binsX = m_bins_arrayX.GetArray();
    const Double_t* binsY = m_bins_arrayY.GetArray();
    const char* labelX = m_labelX.c_str();
    const char* labelY = m_labelY.c_str();

    MH2D* selection_mc = new MH2D(Form("selection_mc_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY); 

    MH2D* bg      = new MH2D(Form("bg_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    MH2D* bg_loW  = new MH2D(Form("bg_loW_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    MH2D* bg_midW = new MH2D(Form("bg_midW_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    MH2D* bg_hiW  = new MH2D(Form("bg_hiW_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    MH2D* effnum  = new MH2D(Form("effnum_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    MH2D* effden  = new MH2D(Form("effden_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    const bool clear_bands = true;
    m_selection_mc = CVH2DW(selection_mc, systematic_univs,       clear_bands);
    m_bg           = CVH2DW(bg,           systematic_univs,       clear_bands);
    m_bg_loW       = CVH2DW(bg_loW,       systematic_univs,       clear_bands);
    m_bg_midW      = CVH2DW(bg_midW,      systematic_univs,       clear_bands);
    m_bg_hiW       = CVH2DW(bg_hiW,       systematic_univs,       clear_bands);
    m_effnum       = CVH2DW(effnum,       systematic_univs,       clear_bands);
    m_effden       = CVH2DW(effden,       systematic_univs_truth, clear_bands);

    delete selection_mc;
    delete bg;
    delete bg_loW;
    delete bg_midW;
    delete bg_hiW;
    delete effnum;
    delete effden;
  }


  template<typename T>
  void Histograms2D::InitializeSidebandHists(T systematic_univs) {
    const Double_t* binsX = m_bins_arrayX.GetArray();
    const Double_t* binsY = m_bins_arrayY.GetArray();
    const char* labelX = m_labelX.c_str();
    const char* labelY = m_labelY.c_str();
    MH2D* wsidebandfit_sig  = new MH2D(Form("wsidebandfit_sig_%s_vs_%s",  labelX,   labelY),
                                       Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                       NBinsY(),     binsY);

    MH2D* wsidebandfit_loW  = new MH2D(Form("wsidebandfit_loW_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    MH2D* wsidebandfit_midW = new MH2D(Form("wsidebandfit_midW_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    MH2D* wsidebandfit_hiW  = new MH2D(Form("wsidebandfit_hiW_%s_vs_%s",  labelX,   labelY),
                                  Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                  NBinsY(),     binsY);

    const bool clear_bands = true;
    m_wsidebandfit_sig  = CVH2DW(wsidebandfit_sig,  systematic_univs, clear_bands);
    m_wsidebandfit_loW  = CVH2DW(wsidebandfit_loW,  systematic_univs, clear_bands);
    m_wsidebandfit_midW = CVH2DW(wsidebandfit_midW, systematic_univs, clear_bands);
    m_wsidebandfit_hiW  = CVH2DW(wsidebandfit_hiW,  systematic_univs, clear_bands);

    delete wsidebandfit_sig;
    delete wsidebandfit_loW;
    delete wsidebandfit_midW;
    delete wsidebandfit_hiW;
  }


  void Histograms2D::InitializeStackedHists() {
    // Event Selection Stacked
    m_stacked_w       = StackedHistogram2D<WType>(
                	        m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX,
				m_bins_arrayY, int(kNWTypes));

    m_stacked_sigbg   = StackedHistogram2D<SignalBackgroundType>(
        	                m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNSignalBackgroundTypes));

    m_stacked_wbg     = StackedHistogram2D<WBackgroundType>(
                        	m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNWBackgroundTypes), 4);

    m_stacked_mesonbg = StackedHistogram2D<MesonBackgroundType>(
        	                m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNMesonBackgroundTypes));

    m_stacked_hadron  = StackedHistogram2D<HadronType>(
                        	m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNHadronTypes), 3);

    m_stacked_fspart  = StackedHistogram2D<FSParticleType>(
        	                m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNFSParticleTypes));

    m_stacked_channel = StackedHistogram2D<ChannelType>(
                        	m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNChannelTypes));

    m_stacked_npi     = StackedHistogram2D<NPionsType>(
                        	m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNNPionsTypes), 2);

    m_stacked_npi0    = StackedHistogram2D<NPi0Type>(
                        	m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNNPi0Types), 2);

    m_stacked_npip    = StackedHistogram2D<NPipType>(
                        	m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNNPipTypes), 2);

    m_stacked_coherent = StackedHistogram2D<CoherentType>(
                         	m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, int(kNCoherentTypes), 5);

    // Sideband Stacked
    m_stacked_wsideband = StackedHistogram2D<WSidebandType>(
                          	m_labelX, m_labelY, m_xlabelX, m_xlabelY, m_bins_arrayX, 
				m_bins_arrayY, kNWSidebandTypes,
                              	sidebands::kWSideband_ColorScheme );
  }

  void Histograms2D::InitializeDataHists() {
    const Double_t* binsX = m_bins_arrayX.GetArray();
    const Double_t* binsY = m_bins_arrayY.GetArray();
    const char* labelX = m_labelX.c_str();
    const char* labelY = m_labelY.c_str();
    m_selection_data    = new MH2D(Form("selection_mc_%s_vs_%s", labelX, labelY),
                                   Form("%s_%s", labelX, labelY), NBinsX(), binsX,
				   NBinsY(), binsY);

    m_wsidebandfit_data = new MH2D(Form("wsidebandfit_data_%s_vs_%s", labelX, labelY),
                                   Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                   NBinsY(), binsY);

    m_wsideband_data    = new MH2D(Form("wsideband_data_%s_vs_%s", labelX, labelY),
                                   Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                   NBinsY(), binsY);

    m_bg_subbed_data    = new MH2D(Form("bg_subbed_data_%s_vs_%s", labelX, labelY),
                                   Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                   NBinsY(), binsY);

    m_unfolded          = new MH2D(Form("unfolded_%s_vs_%s", labelX, labelY),
                                   Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                   NBinsY(), binsY);

    m_cross_section     = new MH2D(Form("cross_section_%s_vs_%s", labelX, labelY),
                                   Form("%s_%s", labelX, labelY), NBinsX(), binsX,
                                   NBinsY(), binsY);
  }


  template<typename T>
  void Histograms2D::InitializeMigrationHist(T systematic_univs) {
    const Double_t* binsX = m_bins_arrayX.GetArray();
    const Double_t* binsY = m_bins_arrayY.GetArray();
    const char* labelX = m_labelX.c_str();
    const char* labelY = m_labelY.c_str();
    const Int_t migrationBinsX = (NBinsX()+2) * (NBinsY()+2); 
    const Int_t migrationBinsY = (NBinsX()+2) * (NBinsY()+2);

    PlotUtils::MnvH2D* migration = new PlotUtils::MnvH2D(Form("migration_%s_vs_%s", 
                                         labelX, labelY), Form("%s_%s", labelX, labelY),
                                         migrationBinsX, 0.0, Double_t(migrationBinsX),
 					 migrationBinsY, 0.0, Double_t(migrationBinsY));
    PlotUtils::MnvH2D* migration_reco = new PlotUtils::MnvH2D(Form("migration_reco_%s_vs_%s", labelX, 							      labelY), Form("%s_%s", labelX, labelY),
							 NBinsX(), binsX, NBinsY(), binsY);
    PlotUtils::MnvH2D* migration_true = new PlotUtils::MnvH2D(Form("migration_truth_%s_vs_%s", labelX,
							 labelY), Form("%s_%s", labelX, labelY),
							 NBinsX(), binsX, NBinsY(), binsY);

  /*  m_response = 0; new PlotUtils::MnvH2D(Form("migration2D_%s_vs_%s", labelX,
                                       labelY), Form("%s_%s", labelX, labelY),
                                       NBinsX(), binsX, NBinsY(), binsY);
    m_responseReco = 0; new PlotUtils::MnvH2D(Form("migration2D_%s_vs_%s_reco", labelX,
                                       labelY), Form("%s_%s", labelX, labelY),
                                       NBinsX(), binsX, NBinsY(), binsY);
    m_responseTrue = 0; new PlotUtils::MnvH2D(Form("migration2D_%s_vs_%s_true", labelX,
                                       labelY), Form("%s_%s", labelX, labelY),
                                       NBinsX(), binsX, NBinsY(), binsY);*/
//    m_response.Setup(Form("Migration2d_%s_vs_%s", labelX, labelY),
//               Form("Migration2d %s vs %s", labelX, labelY),
//               NBinsX(), binsX, NBinsY(), binsY, systematic_univs);
    const bool clear_bands = true;
    m_migration  = CVH2DW(migration, systematic_univs, clear_bands);
    m_migration_reco = CVH2DW(migration_reco, systematic_univs, clear_bands);
    m_migration_true = CVH2DW(migration_true, systematic_univs, clear_bands);
    delete migration;
    delete migration_reco;
    delete migration_true;
  }


//Accessor functions
  std::map<WType, MH2D*> Histograms2D::GetStackMap (WType type ) const {
    return m_stacked_w.m_hist_map;
  }


  std::map<SignalBackgroundType, MH2D*> Histograms2D::GetStackMap (SignalBackgroundType type ) const {
    return m_stacked_sigbg.m_hist_map;
  }


  std::map<WBackgroundType, MH2D*> Histograms2D::GetStackMap (WBackgroundType type ) const {
    return m_stacked_wbg.m_hist_map;
  }


  std::map<MesonBackgroundType, MH2D*> Histograms2D::GetStackMap (MesonBackgroundType type ) const {
    return m_stacked_mesonbg.m_hist_map;
  }


  std::map<HadronType, MH2D*> Histograms2D::GetStackMap (HadronType type ) const {
    return m_stacked_hadron.m_hist_map;
  }


  std::map<FSParticleType, MH2D*> Histograms2D::GetStackMap (FSParticleType type ) const {
    return m_stacked_fspart.m_hist_map;
  }


  std::map<ChannelType, MH2D*> Histograms2D::GetStackMap (ChannelType type ) const {
    return m_stacked_channel.m_hist_map;
  }


  std::map<NPionsType, MH2D*> Histograms2D::GetStackMap (NPionsType type ) const {
    return m_stacked_npi.m_hist_map;
  }


  std::map<NPi0Type, MH2D*> Histograms2D::GetStackMap (NPi0Type type ) const {
    return m_stacked_npi0.m_hist_map;
  }


  std::map<NPipType, MH2D*> Histograms2D::GetStackMap (NPipType type ) const {
    return m_stacked_npip.m_hist_map;
  }


  std::map<WSidebandType, MH2D*> Histograms2D::GetStackMap (WSidebandType type ) const {
    return m_stacked_wsideband.m_hist_map;
  }


  std::map<CoherentType, MH2D*> Histograms2D::GetStackMap (CoherentType type) const {
    return m_stacked_coherent.m_hist_map;
  }


#endif // Histograms_cxx

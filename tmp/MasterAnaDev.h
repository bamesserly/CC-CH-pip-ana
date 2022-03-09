#ifndef MASTERANADEV_H
#define MASTERANADEV_H 1

//Inheritance
#include "AnaUtils/MinervaAnalysisTool.h"
#include "AnaUtils/MCTrack.h"

//***Class Tools***//
class TRandom3;
class INuclearTargetTool;
class IMinervaObjectAssociator;
class IIDAnchoredBlobCreator;
class IAnalysisToolUtils;
class ICalorimetryUtils;
class IMichelTool;
class IImprovedMichelTool;
class IParticleMakerTool;
class IRecoObjectTimeTool;
class IHitTaggerTool;
class IPrimaryBlobProngTool;
class IProtonUtils;
class IMinervaCoordSysTool;
class IVertexEnergyStudyTool;
class IConeUtilsTool;
class Cone;
class IMinervaMathTool;
class IProngClassificationTool;
class IAnchoredTrackFormation;
class IMasterAnaDevRecoUtils;
class INeutronBlobRecoTool;
class INeutronBlobUtils;
class IAnaBlobUtils;
class IRecoilReconstructor;
class IPathLengthTool;
class ILikelihoodParticleTool;
class IMCTrackTool;
//***veto wall**//
class IGetVetoEfficiency;
class ITrackLinearPropagator;


//for Machine Learning
class IEnergyCorrectionTool;
class IGeomUtilSvc;
class IMLVFTool;
//-- from ROOT
class TFile;
class TTree;
//class DeDetector;
//Forward declarations
#include "Event/MinervaEventFwd.h"
namespace Minerva {
class NuclearTarget;
class DeDetector;
//add for vetowall
class DeCryoTarget;
class DeVetoDetector;

}


namespace MCPartStatus
{
  enum t_type
  {
    kPrimary,
    kFinal
  };
}

namespace MCPartType
{
  enum t_type
  {
    kAny,
    kn,
    klown,
    kmidn,
    khighn,
    kp,
    kpi,
    kmeson,
    kmu,
    kem,
    kXtalk,
    kOther
  };
}


class MasterAnaDev : public MinervaAnalysisTool {

 public:

  // Standard constructor
  MasterAnaDev( const std::string& type, const std::string& name, const IInterface* parent );

  // Destructor
  ~MasterAnaDev(){};

  StatusCode initialize();
  StatusCode finalize();

  //mandatory method
  StatusCode reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth = NULL ) const;

  //mandatory
  StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, std::vector<Minerva::NeutrinoInt*>& nuInts ) const;

  //mandatory
  StatusCode tagTruth( Minerva::GenMinInteraction* truth) const;

protected:

    //! Initialize the map from z position to plane
    StatusCode fillPlaneLowZMap();

    //! Get the module/plane at the z position if there is one at this z position
    bool getVtxPlane( double z, int & mod, int & plane ) const;


    //! Given just a module number, return the ModuleID
    Minerva::ModuleID getModuleID( int mod ) const;

    /*!
      @brief Get the target of this event after applying a Z cut that the point is within N(M) planes US(DS) of the target.
      @return A pointer to the NuclearTarget for which point event passes the z cuts
      @param[in] p The Point you are checking against target z positions
      @param[out] targetCode targetID*1000 + targetZ.  (targetZ of 0 means unkown material)
     */
    const Minerva::NuclearTarget* getTarget( const Gaudi::XYZPoint& p, int& targetCode ) const;

    //! Overload of getTarget for GenMinInteraction
    const Minerva::NuclearTarget* getTarget( const Minerva::GenMinInteraction* truth, int& targetCode ) const;

    /*!
      @brief Get the target in which this XYZ point sits.
      @return A pointer to the NuclearTarget for which point event passes the z cuts
      @param[in] p The Point you are checking against target z positions
      @param[out] targetCode targetID*1000 + targetZ.  (targetZ of 0 means unkown material)
     */
    const Minerva::NuclearTarget* getExactTarget( const Gaudi::XYZPoint& p, int& targetCode ) const;

    //! Overload of getExactTarget for GenMinInteraction
    const Minerva::NuclearTarget* getExactTarget( const Minerva::GenMinInteraction* truth, int& targetCode ) const;

    /*!
      @brief Adjust the vtx to be the point where the muon intersects the zCenter of a valid nuclear target
      @return true if vtx has been associated with a nuclear target
      @param[in,out] vtx The vertex you are adjusting
      @param[in] primMuon The muon you are using to project
     */
    bool adjustVertexIntoPassiveTarget( SmartRef<Minerva::Vertex>& vtx, const SmartRef<Minerva::Particle>& muonPart ) const;

      /*!
      @brief If this xy point were in passive targets 1,2..5 what would the Z of the material be?
      @return true if everything worked
      @param[out] intVec A vector of what Z would be in the passive target
      @param[in] x X position of the point
      @param[in] y Y position of the point
     */
    bool getReferenceZ( std::vector<int>& intVec, double x, double y ) const;

    /*!
      @brief If this xy point were in passive targets 1,2..5 what would be the smallest distance to the division of materials?
      @return true if everything worked
      @param[out] dVec A vector of what the distance to the closet division between materials would be in the passive target
      @param[in] x X position of the point
      @param[in] y Y position of the point
     */
    bool getReferenceDistToDivision( std::vector<double>& dVec, double x, double y ) const;

    //! What is the distance from this z point to each of the passive targets
    bool getDistToTarget( std::vector<double>& dVec, double z ) const;
     /*! @brief Vertex must be closer than this many planes DS of the NuclearTarget
        @property{DSVertexPlanes,3}*/
    int    m_nplanes_DSVertexLimit;

    TRandom3*                 m_randomGen;
    unsigned long int         m_randomSeed;
    /*! @brief Vertex must be closer than this many planes US of the NuclearTarget
        @property{USVertexPlanes,1} */
    int    m_nplanes_USVertexLimit;

    /*! @brief Vertex must be closer than this to the upstream edge of the NuclearTarget
        @property{USVertexCut, 75*mm}*/
    double m_USVertexCut;

    /*! @brief Vertex must be closer than this to the downstream edge of the NuclearTarget
        @property{DSVertexCut, 75*mm}*/
    double m_DSVertexCut;

    /*! @brief Use the nplanes away from target criteria instead of z distance from target.
        @property{UsePlaneZCut, false}*/
    bool m_usePlanesZCut;

   /*! Lowest module number for US boundary of fiducial volume (inclusive)
        @property{ USFiducialMod, -2 }*/
    int m_usFiducialMod;

    /*! Highest module number for DS boundary of fiducial volume (inclusive)
        @property{ DSFiducialMod, 80 }*/
    int  m_dsFiducialMod;

    /*! Lowest module number for US boundary of analyzable volume (inclusive)
        @property{ USAnalyzableMod, -3 }
             */
    int m_usAnalyzableMod;

    /*! Highest module number for DS boundary of analyzable volume (inclusive)
        @property{ DSAnalyzableMod, 114 }
             */
    int m_dsAnalyzableMod;

 private:

      // Fiducial volume for pions studies
      double m_fidHexApothem;
      double m_fidUpStreamZ;
      double m_fidDownStreamZ;

      IGeomUtilSvc*       m_geomUtilSvc;

      //Muon q/p cut
      double              m_qOverpChargeCut;

      double m_beamAngleBias;

      //Binding Energy
      double              m_nuCCQEBindingEnergyMeV;

      //Muon time window parameters
      double              m_lowerTimeWindow;
      double              m_upperTimeWindow;

    std::vector<double> m_nuVtxBlobRadii;

    /*! @brief Vertex blob radii for antineutrino CCQE
    @property{VertexBlobRadius_Antinu, 100. * CLHEP::mm...}
    */
    std::vector<double> m_antinuVtxBlobRadii;

    bool m_doVtxBlobSystematic;

    /*! @brief Vertex blob radii for neutrino CCQE
    @property{VertexBlobRadii_Nu, 300. * CLHEP::mm...}
    */
      //faiza
    /*! @brief Binding energy for neutrino CCQE
      @property{CCQEBindingE_Nu, 34 * CLHEP::MeV}
     */
    double m_nuCCQEBindingE;

    /*! @brief Binding energy for antineutrino CCQE
      @property{CCQEBindingE_Antinu, 30 * CLHEP::MeV}
     */
    double m_antinuCCQEBindingE;
    /*! @brief Do you want to use the CaloTuning from the CC inclusives in tracker analysis for the tracker events?
      If false, NukeCC's tuning for the tracker will be used.
      @property{ UseCCTrackerTuning, false }
     */
    bool m_useCCTrackerTuning;


      //Primary tag names
      std::string         m_primaryHadron;
      std::string         m_primaryProton;
      std::string         m_secondaryProtons;

      //Particle scores
      double              m_minMuonScore;
      double              m_minProtonScore;

      //Isolated EM Blobs parameters

      //Maximum number of isolated blobs
      double              m_maxIsoBlobs;

      //Prong, Blob and Cluster colors
      int                 m_hadronProngsColor;
      int                 m_vertexBlobProngColor;
      int                 m_isolatedBlobProngColor;
      int                 m_dispersedEnergyColor;
      int                 m_protonProngColor;
      int                 m_coneEnergyColor;
      int                 m_neutronBlobColor;

      //Make Short Tracks
      bool                m_makeShortTracks;
      bool                m_runVertexAnchoredShortTracker;
      bool                m_runVertexEnergyStudyTool;

      bool                m_doRecursiveVertexAnchoredShortTracking;
      bool                m_doRecursiveVESTool;

      //Make Iso/Fuzz EM Blobs
      bool		  m_makeIsoMuonEMBlobs;
      bool                m_makeFuzzMuonEMBlobs;

      //Make Vertex Blobs
      unsigned int        m_numSearchRadii;
      double              m_searchStepSize;

      double              m_maxSearchDistance;
      double              m_maxStartingDistance;
      double              m_maxAllowedSearchGap;

      bool                m_makeFilamentStyleVtxBlobProngs;
      bool                m_makeMultipleRadiiStyleVtxBlobProngs;
      double              m_maxSeparationBlobVertex;

      //Truth Michel Parameters
      double              m_michel_upstreamZ;
      double              m_michel_downstreamZ;

      bool                m_fillTruthTG4ProtonCone;

      //Proton Conditions
      bool                m_useOdMatchProtons;
      double              m_maxProtonChi2;
      double              m_ProtonZLow;
      double              m_ProtonZHigh;
      double              m_ProtonApothem;

      //muonIsPlausible
      bool                m_useMuonIsPlausible;
      bool                m_doNeutron;
      bool                m_doNeutronInTrackerOnly;
      double              m_muonAxisCylinderRadius;

      //Neutron Variables
      size_t              m_maxRecoCluster;

      //for Machine Learning
      bool m_getLatticeEnergies;
////////////////////////////////////////////////////////RECOIL////////////////////
      IIDAnchoredBlobCreator*    m_vtxBlobCreator;   ///< Tool used to make vertex blob to blind the vertex activuty


     //Nuclear targets
     bool m_addMyNuclearTargets;


     //MLVFTool for Machine Learning
      IMLVFTool*                 m_MLVFTool;
      std::string                m_MLVFToolAlias;

      //ParticleMaker Tool variables
      IParticleMakerTool*       m_particleMaker;
      std::string               m_particleMakerAlias;

      // Tool to perform functions common to CC analyses
      IAnalysisToolUtils*        m_ccUtils;

      //Hit Tagger Tool
      IHitTaggerTool*           m_hitTaggerTool;
      std::string               m_hitTaggerToolAlias;

      //Calorimetry Utils Tool
      ICalorimetryUtils*        m_caloUtils;
      std::string               m_caloUtilsAlias;

      IRecoilReconstructor*      m_recoilReconstructor;      ///< Reconstructs the recoil system energy

      //PrimaryBlobProng Tool
      IPrimaryBlobProngTool*    m_primaryBlobProngTool;
      std::string               m_primaryBlobProngToolAlias;

      //ProtonUtils Tool
      IProtonUtils*             m_protonUtils;
      std::string               m_protonUtilsAlias;

      //Michel Tool
      IMichelTool*              m_michelTool;
      std::string               m_michelToolAlias;

      Minerva::DeDetector *     m_InnerDetector;	

      IMCTrackTool*             m_MCTrackTool;
 
     //iImprovedMichelTool
      IImprovedMichelTool*      m_improvedmichelTool;
      std::string               m_improvedmichelToolAlias;


      //< Tool to organize and use NuclearTargets
      INuclearTargetTool       *m_nukeTool;
      std::string         m_nukeToolAlias;


      //MinervaObjectAssociator
      IMinervaObjectAssociator* m_objectAssociator;
      std::string               m_minObjAssocAlias;

      IIDAnchoredBlobCreator*   m_vtxBlobCreatorTool;
      std::string               m_vtxBlobCreatorAlias;

      IRecoObjectTimeTool*      m_recoObjTimeTool;
      std::string               m_recoObjTimeToolAlias;

      IMinervaCoordSysTool*     m_minCoordSysTool;
      std::string               m_minCoordSysToolAlias;

      //VertexEnergyStudyTool
      IVertexEnergyStudyTool*   m_vertexEnergyStudyTool;
      std::string               m_vtxEngStudyToolAlias;

      IBlobCreatorUtils*        m_blobCreatorUtils;
      std::string               m_blobCreatorUtilsAlias;

      IConeUtilsTool*           m_coneUtilsTool;
      std::string               m_coneUtilsToolAlias;

      static IMinervaMathTool*  m_mathTool;
      std::string               m_mathToolAlias;

      IPathLengthTool*          m_pathLengthTool;
      std::string               m_pathLengthToolAlias;

      IProngClassificationTool* m_prongIntersectionTool;
      std::string               m_prongIntersectionToolAlias;

      IAnchoredTrackFormation*  m_anchoredShortTracker;
      std::string               m_shortTrackToolAlias;

      ILikelihoodParticleTool*  m_LikelihoodPIDTool;
      std::string               m_LPIDToolAlias;

      IMasterAnaDevRecoUtils*   m_masterAnaDevRecoUtils;
      std::string               m_masterAnaDevRecoUtilsAlias;

      INeutronBlobRecoTool*     m_neutronBlobRecoTool;
      std::string               m_neutronBlobRecoToolAlias;
      INeutronBlobUtils*        m_neutronBlobUtils;
      std::string               m_neutronBlobUtilsAlias;
      IAnaBlobUtils*            m_anaBlobUtils;
      std::string               m_anaBlobUtilsAlias;
     //Tool used to correct the muon momentum for the passive target material
      IEnergyCorrectionTool*     m_energyCorrectionTool;
      //---------------------
      //Stores ID detector
      const Minerva::DeDetector* m_idDet;
      //DeDetector* m_idDet;
  //---------------------
  // Vetowall
  //---------------------
      bool m_VetoWallIsImportant;

      std::string m_propagateToVetoName;

      std::string m_getVetoEffName;
      Minerva::DeVetoDetector*        m_VetoDet;
      std::string m_timeToolName;
      IRecoObjectTimeTool*       m_timeTool;
      std::vector<std::pair <int,int> >discrPairToPaddle(unsigned int minervaID) const;

      IPlexModel*   m_plexModel;
      ITrackLinearPropagator*    m_propagateToVeto;
      IGetVetoEfficiency*     m_getVetoEff;
      bool m_hasVeto;
       int m_NumVetoPaddles;
       int m_VetoPaddlePMT2;
       int m_VetoPaddlePMT1;
       int m_totalVETOPMTS;
       double m_PanelOffset;
      //---------------------
       // CryoTarget
       //---------------------


      bool m_hasCryo;
      Minerva::DeCryoTarget*           m_CryoDet;



      //---------------------
      //Helper functions
      //---------------------

    
//--------------------------------------------------------------------faiza
    IAnalysisToolUtils*        m_ccUtilsWideTimeWindow;    ///< Tool to perform functions common to CC analyses with a wide time window for systematics

    //! Get a pass control tag <pass>_<ana_signature>_<tgt>
    inline std::string getPassTag( const std::string& tgt ) const { return std::string( "pass_" + m_anaSignature + "_" + tgt ); };
    //! Get a pass control tag <pass>_<ana_signature>
    inline std::string getPassTag( )                        const { return std::string( "pass_" + m_anaSignature ); };
    //! Select the appropriate channel and use RecoilUtils to get the recoil energy
    std::string getRecoilChannel(const int muonCharge, const int targetID, const int targetZ ) const;
    double getRecoilEnergy( const Minerva::PhysicsEvent* event, const SmartRef<Minerva::Prong>& muonProng, const std::string& channel, Minerva::NeutrinoInt *eventHyp = 0 ) const;
     double getvisibleenergy(const Minerva::PhysicsEvent* event, const SmartRef<Minerva::Prong>& muon) const;
  StatusCode getBestParticle( SmartRef<Minerva::Prong> prong, SmartRef<Minerva::Particle>& particle, Minerva::Particle::ID partType) const;

      void calcMomentumByRangeVars(Minerva::PhysicsEvent* event) const;
      void tagMichelElectrons(Minerva::PhysicsEvent* event) const;
      void fillIsoProngInfo(Minerva::PhysicsEvent* event, double muon_time) const;
      Minerva::IDClusterVect getHadronIDClusters(const Minerva::PhysicsEvent* event, SmartRef<Minerva::Prong> muonProng, SmartRef<Minerva::Prong> excludeProng = (Minerva::Prong*)NULL ) const;
      double calcOpeningAngle(Gaudi::LorentzVector& vec1, Gaudi::LorentzVector& vec2) const;
  //! Declare the branches used for CCQE blobs
    void declareBlobBranches( const std::string& prefix ); 

      // truthIsPlausible - a pure virtual method
      bool truthIsPlausible( const Minerva::PhysicsEvent *event ) const;

      bool getProtonProngs( Minerva::ProngVect prongs, Minerva::ProngVect& secondaryProtonProngs, SmartRef<Minerva::Prong> &protonProng, SmartRef<Minerva::Particle> &protonPart ) const;

      void TruthFittedMichel(Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth) const;

      bool tagMichels( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth ) const;
      bool ImprovedtagMichels( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth ) const;

      bool hasTruthMichelElectron( Minerva::PhysicsEvent* event = NULL ) const;

      void tagProtonProngTruth( Minerva::PhysicsEvent *event, SmartRef<Minerva::Prong> protonProng ) const;
      void tagSecondaryProtonProngTruth( Minerva::PhysicsEvent *event, Minerva::ProngVect protonProng ) const;

      std::vector<double> tagMichelProngTruth( Minerva::Prong michelProng ) const;


 
//! Get blob-analyzable clusters and mark them as unused
    SmartRefVector<Minerva::IDCluster> getBlobbingClusters( Minerva::PhysicsEvent* event, const SmartRef<Minerva::Prong>& primMuonProng ) const;

    //! Remove Used clusters from the input param, return the unused clusters
    SmartRefVector<Minerva::IDCluster> removeUsedClusters( SmartRefVector<Minerva::IDCluster>& clusters ) const;

    //! Fill the event branches for a blob using its clusters
    //! Returns the energy of this blob's nucl, tracker, ecal clusters (ccqe recoil regions)
    double fillBlobBranches( Minerva::PhysicsEvent* event, SmartRefVector<Minerva::IDCluster>& clusters, const std::string& prefix, int color = 0 ) const;
    double fillBlobBranches( Minerva::PhysicsEvent* event, SmartRefVector<Minerva::IDCluster>& idclusters, SmartRefVector<Minerva::ODCluster>& odclusters, const std::string& prefix, int color = 0 ) const;

     //! Split clusters by ID subdet
        void separateClustersBySubdet(
        const SmartRefVector<Minerva::IDCluster>& clusters,
        SmartRefVector<Minerva::IDCluster>& nucl,
        SmartRefVector<Minerva::IDCluster>& tracker,
        SmartRefVector<Minerva::IDCluster>& ecal,
        SmartRefVector<Minerva::IDCluster>& hcal
        ) const; 

//! Get clusters created by this particle type - Stolen from CCQEAntiNuTool
         template<typename ClusterType, typename DigitType>
         SmartRefVector<ClusterType> getClustersCreatedBy( SmartRefVector<ClusterType> clusters, MCPartType::t_type primType, MCPartType::t_type fsType = MCPartType::kAny ) const;
//
              //! Get particle type that created this cluster - Stolen from CCQEAntiNuTool
         template<typename ClusterType, typename DigitType>
         MCPartType::t_type getParticleCreatedBy( SmartRef<ClusterType> cluster, MCPartStatus::t_type status ) const;
//
    
     bool createParticles( Minerva::PhysicsEvent* event, Minerva::ProngVect& hadronProngs ) const;
      //Nuclear targets: List if Nuclear targets
          std::vector<Minerva::NuclearTarget*> m_targets;

      std::vector<double> tagBlobTruth( Minerva::IDBlob blob ) const;

      std::vector<double> tagDigitsTruth( SmartRefVector<Minerva::IDDigit> digits ) const;
      bool m_correctMuonEnergyToTarget;

    /*! @brief Do you want to use the CaloTuning from the CC inclusives in tracker analysis for the tracker events?
 *       If false, NukeCC's tuning for the tracker will be used.
 *             @property{ UseCCTrackerTuning, false }
 *                  */
      bool findLongestMinervaTrack( Minerva::Prong* deProng, SmartRef<Minerva::Track> &longestMinervaTrack ) const;
//by Zubair
//!< Map from z-position of front edge of a plane to <mod,plane>
    std::map< double,std::pair<int,int> >       m_planeLowZMap;


     Minerva::ModuleID m_usFiducialModID, m_dsFiducialModID;
     Minerva::ModuleID m_usAnalyzableModID, m_dsAnalyzableModID;
//by Oscar
   bool m_useDSCal;
 //Declare DB File related variables
   bool m_useDNN;
   std::string m_getMLPredFilename;
   std::string filename;
   TFile *DBPredFile;
   TTree *dbPred;

      bool findLongestMinervaTrack( Minerva::ProngVect deProngVec, SmartRef<Minerva::Track> &longestMinervaTrack ) const;

      bool findLongestMinervaTrack( const Minerva::PhysicsEvent *event, SmartRef<Minerva::Track> &longestMinervaTrack ) const;

  protected:

    SmartRefVector<Minerva::IDCluster> getAnalyzableNonMuIDClusters(
        const Minerva::PhysicsEvent* event,
        const SmartRef<Minerva::Prong>& muon,
        bool useWideWindow = false
        ) const;
    SmartRefVector<Minerva::ODCluster> getAnalyzableNonMuODClusters(
        const Minerva::PhysicsEvent* event,
        const SmartRef<Minerva::Prong>& muon,
        bool useWideWindow = false
        ) const;



      bool adjustDNNVertexIntoCenterOfSegment( SmartRef<Minerva::Vertex>& vtx, double zCenter, const SmartRef<Minerva::Particle> & muonPart ) const;
//      void calcMomentumByRangeVars(Minerva::PhysicsEvent* event) const;
//      void tagMichelElectrons(Minerva::PhysicsEvent* event) const;
//      void getParticleRecoilE(const Minerva::PhysicsEvent* event,SmartRef<Prong> muonProng, double& proton, double& neutron, double& pion,   double& kaon, double& em, double& other,  double& xtalk) const;
      int GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane ) const;
void lonerVertex( Minerva::PhysicsEvent* event, Minerva::Track *muTrack, SmartRefVector<Minerva::IDCluster> allClusters, double targetZEnd ) const;
};

#endif // MASTERANADEV_H

#include "MasterAnaDev.h"
#include "GiGaCnv/IGiGaGeomCnvSvc.h"

//Tools - General and Analysis Specific
#include "MinervaUtils/IMinervaObjectAssociator.h"
#include "RecInterfaces/IRecoObjectTimeTool.h"
#include "RecInterfaces/IFiducialPointTool.h"
#include "GeoUtils/IMinervaCoordSysTool.h"
#include "MinervaUtils/IMinervaMathTool.h"
#include "GeoUtils/NuclearTarget.h"
#include "GeoUtils/INuclearTargetTool.h"
#include "MasterAnaDev/IMasterAnaDevRecoUtils.h"
#include "BlobFormation/IBlobCreatorUtils.h"
#include "BlobFormation/IBlobFitTool.h"
#include "EnergyRecTools/IEnergyCorrectionTool.h"
#include "BlobFormation/INeutronBlobRecoTool.h"
#include "BlobFormation/INeutronBlobUtils.h"
#include "MasterAnaDev/IAnaBlobUtils.h"
#include "MasterAnaDev/IElectronUtils.h"
#include "MasterAnaDev/IRegisterBranches.h"
#include "BlobFormation/NeutronDataStructDef.h"
#include "MasterAnaDev/INeutronInelasticReweightTool.h"

#include <GeoUtils/IPathLengthTool.h>
#include "ParticleMaker/ILikelihoodParticleTool.h"

#include "EventRecInterfaces/IPrimaryBlobProngTool.h"
#include "ProngMaker/IMichelTool.h"
#include "EnergyRecTools/ICalorimetryUtils.h"
#include "AnaUtils/IProtonUtils.h"
#include "AnaUtils/AnaFilterTags.h"
#include "AnaUtils/IMuonUtils.h"
#include "AnaUtils/IAnalysisToolUtils.h"
#include "AnaUtils/NuMIVars.h"
#include "MinervaDet/IGeomUtilSvc.h"
#include "MinervaDet/DeDetector.h"
#include "MinervaDet/DePlane.h"
#include "DetDesc/Material.h"
#include "AnaUtils/IMCTrackTool.h"
#include "AnaUtils/MCTrack.h"

#include "BlobFormation/IIDAnchoredBlobCreator.h"
#include "AnaUtils/IPhysicsCalculator.h"
#include "TruthMatcher/ITruthMatcher.h"
#include "ParticleMaker/IParticleMakerTool.h"
#include "RecUtils/ParticleExtraDataDefs.h"
#include "RecUtils/VertexExtraDataDefs.h"

#include "MinervaUtils/IHitTaggerTool.h"

#include "Event/DAQHeader.h"
#include "Event/PhysicsEvent.h"
#include "Event/TimeSlice.h"
#include "Event/MCIDDigit.h"
#include "Event/MCODDigit.h"

#include "Event/GenMinHeader.h"

#include "ProngMaker/IImprovedMichelTool.h" //ImprovedMichelTool
#include "MinervaUtils/MinervaObjectSort.h"

#include <TFile.h>
#include <TTree.h>

#include <TString.h>
#include <TRandom3.h>
#include <TMath.h>
#include <set>
#include <numeric>
//From GaudiException
#include "GaudiKernel/PhysicalConstants.h"
#include <AnaUtils/IMLVFTool.h>
#include "EnergyRecTools/IRecoilReconstructor.h"
#include "EnergyRecTools/IEnergyCorrectionTool.h"
#include "RecUtils/ITrackAddClusters.h"
#include "RecInterfaces/ITrackPropagator.h"
#include "AnaUtils/IRecoilUtils.h"
//--------------------------------------------------------------
//Vetowall tools
#include "VetoDet/DeVetoWall.h"
#include "VetoDet/DeVetoDetector.h"

#include "CalTools/IGetVetoEfficiency.h"
#include "RecInterfaces/ITrackLinearPropagator.h"
#include "Event/VetoDigit.h"
#include "Kernel/ChannelID.h"
#include <Kernel/ChannelID.h>
#include "Kernel/DiscrPairID.h"
#include <Plex/IPlexModel.h>
#include "Event/MCVetoDigit.h"
#include "Plex/IPlexModel.h"
//Cryotarget
#include "MinervaDet/DeCryoTarget.h"
// ElasticFilterTool
#include "TrackShortPatRec/IElasticFilterTool.h"


DECLARE_TOOL_FACTORY( MasterAnaDev );
// m_mathTool is STATIC so it has to be initialized
// somewhere outside of class scope.
IMinervaMathTool* MasterAnaDev::m_mathTool = NULL;

using namespace std;
namespace
{
  //cc inclusive blob
  const string RECOILBLOB_PREFIX = "recoil";
  int MAX_TRACKER_MODULE; //!< No events downstream of this module will be recorded
  const int FIRST_TRACKER_MODULE = 23; //!< Lowest module number in the tracker

  //these numbers were calculated as half the distance between adjacent plane's zCenter
  const double PLANE_WIDTH    = 18.5  * CLHEP::mm;    //!< longitudinal distance from start to end of a DePlane
  const double MOD_GAP        = 6.235 * CLHEP::mm;    //!< Air gap between plane2-1 of adjacent modules
  const double MOD_GAP_NUKE   = 5.235 * CLHEP::mm;    //!< Air gap between plane2-1 of adjacent modules in nuclear region
  const double PLANE_GAP      = 1.98  * CLHEP::mm;    //!< Air gap between plane1-2 in the same module
  //Oficial -- DSCAL
  const double GAP_PLANE1_ECAL  = 3.84  * CLHEP::mm;
  const double GAP_MOD85        = 0.43  * CLHEP::mm;
  const double ECAL_HALF_WIDTH  = 11.21 * CLHEP::mm;
  const double HCAL_HALF_WIDTH  = 34.63 * CLHEP::mm;
  const double GAP_ECAL_HIGH    = 20.46 * CLHEP::mm;
  const double HCAL_HIGH        = 47.31 * CLHEP::mm;
  const double GAP_LEAD_94      = 5.82  * CLHEP::mm;

  //see http://www.rapidtables.com/web/color/RGB_Color.htm
  const int m_muonColor          = 0x32CD32;  //lime green
  const int m_vtxBlobColor       = 0x9932CC;  //dark orchid
  const int m_muonFuzzColor      = 0xB5FF20;  //medium spring green
  const int m_isoBlobColor       = 0x1E90FF;  //dodger blue
  const int m_dispColor          = 0xB22222;  //fire brick
  const int m_emShowerColor      = 0xFF6347;  //tomato
  const int m_generalShowerColor = 0xFF4500;  //orange red




  //ccqe blobs
  const string VTXBLOB_PREFIX = "vtx";
  const string MUONFUZZBLOB_PREFIX = "mufuzz";
  const string ISOBLOB_PREFIX = "iso";
  const string DISPBLOB_PREFIX = "disp";
  const string NUEFUZZBLOB_PREFIX = "nuefuzz";

  //if MC, add e by particle type
  std::map<MCPartType::t_type, string> PART_TYPE_MAP;
}
/////////////////////////////////////////////////////////////////////////
//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MasterAnaDev::MasterAnaDev(const std::string& type, const std::string& name, const IInterface* parent )
  : MinervaAnalysisTool( type, name, parent ) {

  declareInterface<IInteractionHypothesis>(this);

  m_anaSignature = "MasterAnaDev";
  m_hypMeths.push_back( "MasterAnaDev" );

  declareProperty( "RandomSeed", m_randomSeed = 0 );
  declareProperty("HypothesisMethods", m_hypMeths  );

  //-- get the DSCal fiducial and analizable modules
  declareProperty("UseDSCal", m_useDSCal = false);

  declareProperty("NuclearTargetToolAlias", m_nukeToolAlias  = "NukeCCInclusiveNukeTool" );
  declareProperty("AddMyNuclearTargets", m_addMyNuclearTargets = true );
  declareProperty("UsePlaneZCut",   m_usePlanesZCut = true );
  declareProperty("DSVertexPlanes", m_nplanes_DSVertexLimit = 2 ); //allow vertices for 2 planes DS
  declareProperty("USVertexPlanes", m_nplanes_USVertexLimit = 1 ); //allow vertices in the plane just US
  declareProperty("USVertexZCut", m_USVertexCut = 75.0*CLHEP::mm );
  declareProperty("DSVertexZCut", m_DSVertexCut = 75.0*CLHEP::mm );

  // by default target1 - module 80 are fiducial
  declareProperty( "USFiducialModule",     m_usFiducialMod = -2 );
  declareProperty( "DSFiducialModule",     m_dsFiducialMod );//DSCAL
  //declareProperty( "DSFiducialModule",     m_dsFiducialMod = 80 );

  // ... and everything is analyzable except first 2 modules and HCAL
  declareProperty( "USAnalyzableModule",   m_usAnalyzableMod = -3 );
  declareProperty( "DSAnalyzableModule",   m_dsAnalyzableMod );//DSCAL
  //declareProperty( "DSAnalyzableModule",   m_dsAnalyzableMod = 94 );

  //! Time Window parameters for cluster selections
  declareProperty( "LowerTimeWindow",m_lowerTimeWindow = -20 * CLHEP::ns );
  declareProperty( "UpperTimeWindow",m_upperTimeWindow =  35 * CLHEP::ns );

  declareProperty("CCQEBindingE_Nu",     m_nuCCQEBindingE       = 34. * CLHEP::MeV);
  declareProperty("CCQEBindingE_Antinu", m_antinuCCQEBindingE   = 30. * CLHEP::MeV);
  declareProperty( "UseCCTrackerTuning", m_useCCTrackerTuning   = false );

  //! Particles Scores
  declareProperty( "MinMuonScore", m_minMuonScore  = 0.9 );
  declareProperty( "MinProtonScore", m_minProtonScore  = 0.05 );
  declareProperty( "BeamAngleBias", m_beamAngleBias = 0.006*CLHEP::radian );
  //! Muon q/p cut
  declareProperty( "qOverpChargeCut", m_qOverpChargeCut = 0. );
  //! Binding Energy
  declareProperty( "NuCCQEBindingEnergyMeV", m_nuCCQEBindingEnergyMeV = 34.0 * CLHEP::MeV );

  //! Making Short Tracks
  declareProperty( "MakeShortTracks",                           m_makeShortTracks                         = true );
  declareProperty( "MakeShortTracksvtxAnchoredShortTracker",   m_runVertexAnchoredShortTracker           = true );
  declareProperty( "DoRecursiveVertexAnchoredShortTracking",    m_doRecursiveVertexAnchoredShortTracking  = true );
  declareProperty( "MakeShortTracksvtxEnergyStudyTool",        m_runVertexEnergyStudyTool                = true );
  declareProperty( "DoRecursiveVESTool",                        m_doRecursiveVESTool                      = true );

  //! Make EM and Fuzz Blobs for the Muon
  declareProperty( "MakeIsoMuonEMBlobs",                        m_makeIsoMuonEMBlobs              = true );
  declareProperty( "MakeFuzzMuonEMBlobs",                       m_makeFuzzMuonEMBlobs             = true );

  //systematics for vertex blob radius
  //  //index zero is the CV
  vector<double> default_nuVtxBlobRadii;
  default_nuVtxBlobRadii.push_back( 80. * CLHEP::mm );
  vector<double> default_antinuVtxBlobRadii;
  default_antinuVtxBlobRadii.push_back( 80. * CLHEP::mm );

  declareProperty("VertexBlobRadii_Nu", m_nuVtxBlobRadii         = default_nuVtxBlobRadii );
  declareProperty("VertexBlobRadii_Antinu", m_antinuVtxBlobRadii = default_antinuVtxBlobRadii );
  declareProperty("DoVertexBlobSystematic", m_doVtxBlobSystematic = true );


  //! Primary Tag names
  declareProperty( "PrimaryHadronTagName",   m_primaryHadron     = "PrimaryHadron" );
  declareProperty( "PrimaryProtonTagName",   m_primaryProton     = "PrimaryProton" );
  declareProperty( "SecondaryProtonTagName", m_secondaryProtons  = "SecondaryProtons" );

  //! Maximum number of isolated blobs
  declareProperty( "MaxIsoBlobs", m_maxIsoBlobs = 2);

  //! For making Vertex Blobs
  declareProperty( "SearchStepSize",                    m_searchStepSize                       = 300.0 * CLHEP::mm );
  declareProperty( "NumSearchRadii",                    m_numSearchRadii                       = 1 );
  declareProperty( "FStyleMaxSearchDistance",          m_maxSearchDistance                    = 300.0 * CLHEP::mm ); //mm
  declareProperty( "FStyleMaxStartingDistance",        m_maxStartingDistance                  = 300.0 * CLHEP::mm ); //mm
  declareProperty( "FStyleMaxAllowedSearchGap",        m_maxAllowedSearchGap                  = 1e6 * CLHEP::mm );   //mm

  //! For making Vertex Blob Prongs
  declareProperty( "MakeFilamentStyleVtxBlobProngs",       m_makeFilamentStyleVtxBlobProngs        = false );
  declareProperty( "MakeMultipleRadiiStyleVtxBlobProngs",  m_makeMultipleRadiiStyleVtxBlobProngs   = true );
  declareProperty( "MaxAllowedSeparationBlobVertex",       m_maxSeparationBlobVertex               = 300.0 * CLHEP::mm );

  //! Loose Fiducial truth Michel electron cut
  declareProperty( "MichelDownstreamZ",          m_michel_downstreamZ          = 8524.19 * CLHEP::mm ); //mm // Module 83?
  declareProperty( "MichelUpstreamZ",            m_michel_upstreamZ            = 5900.91 * CLHEP::mm ); //mm //Module 25?

  //! Cone Parameters
  declareProperty("FillTruthTG4ProtonCone",     m_fillTruthTG4ProtonCone       = false); //mm

  //! Use OD Match prongs for finding protons
  declareProperty( "UseODMatchProtons",       m_useOdMatchProtons       = false);
  declareProperty( "MaxProtonChi2",           m_maxProtonChi2           = 50.  );
  declareProperty( "ProtonEndPointZLow",      m_ProtonZLow              = 4000.0  );  //-- mm
  declareProperty( "ProtonEndPointZHigh",     m_ProtonZHigh             = 10000.0 );  //-- mm
  declareProperty( "ProtonEndPointApothem",   m_ProtonApothem           = 1200  );  //-- mm


  //! Prong, Blob and Cluster Colors
  declareProperty( "IsolatedBlobProngColor",      m_isolatedBlobProngColor       = 0xFF0000 ); // pure red
  declareProperty( "DispersedEnergyColor",        m_dispersedEnergyColor         = 0x00CCFF ); // cyan
  declareProperty( "VertexBlobProngColor",        m_vertexBlobProngColor         = 0x996600 ); // brown
  declareProperty( "ProtonProngColor",            m_protonProngColor             = 0x990000 ); // maroon
  declareProperty( "HadronProngsColor",           m_hadronProngsColor            = 0xFFFF00);  // pure yellow
  declareProperty( "ConeEnergyColor",             m_coneEnergyColor              = 0xFF00FF ); // magenta
  declareProperty( "NeutronBlobColor",            m_neutronBlobColor             = 0x00FF00 ); //lime

  //! Tools and their Aliases
  declareProperty( "MinervaObjAssocAlias",        m_minObjAssocAlias           = "MASTERANADEVMinObjAssociator" );
  declareProperty( "RecoObjTimeToolAlias",        m_recoObjTimeToolAlias       = "MASTERANADEVRecoObjTimeTool" );
  declareProperty( "MinervaCoordSysToolAlias",    m_minCoordSysToolAlias       = "MASTERANADEVMinCoordSysTool" );
  declareProperty( "MinervaMathToolAlias",        m_mathToolAlias              = "MASTERANADEVMinMathTool" );

  declareProperty( "MASTERANADEVRecoUtilsAlias",  m_masterAnaDevRecoUtilsAlias = "MASTERANADEVRecoUtilsTool" );
  declareProperty( "BlobCreatorUtilsAlias",       m_blobCreatorUtilsAlias      = "MASTERANADEVBlobCreatorUtils" );
  declareProperty( "PrimaryBlobProngToolAlias",   m_primaryBlobProngToolAlias  = "MASTERANADEVPrimaryBlobProngTool" );
  declareProperty( "MichelToolAlias",             m_michelToolAlias            = "MASTERANADEVMichelFinderTool" );
  declareProperty( "ImprovedMichelToolAlias",     m_improvedmichelToolAlias    = "MASTERANADEVImprovedMichelFinderTool" );
  declareProperty( "CalorimetryUtilsAlias",       m_caloUtilsAlias             = "MASTERANADEVCaloUtils" );
  declareProperty( "ProtonUtilsAlias",            m_protonUtilsAlias           = "MASTERANADEVProtonUtils" );
  declareProperty( "ParticleToolAlias",           m_particleMakerAlias         = "MASTERANADEVParticleMaker" );

  declareProperty( "HitTaggerToolAlias",          m_hitTaggerToolAlias         = "MASTERANADEVHitTaggerTool" );
  declareProperty( "ElasticFilterToolAlias",      m_elasticFilterToolAlias     = "ElasticFilterTool" );

  declareProperty( "NeutronBlobRecoToolAlias",    m_neutronBlobRecoToolAlias   = "NeutronBlobRecoTool" );
  declareProperty( "NeutronBlobUtilsAlias",       m_neutronBlobUtilsAlias      = "NeutronBlobUtils"    );
  declareProperty( "AnaBlobUtilsAlias",           m_anaBlobUtilsAlias          = "AnaBlobUtils"        );
  declareProperty( "BlobFitToolAlias",            m_blobFitToolAlias           = "BlobFitTool"        );


  // analysis-specific options
  declareProperty( "UseMuonIsPlausible",          m_useMuonIsPlausible         = true );
  declareProperty( "DoNeutronBranches",           m_doNeutron                  = true );
  declareProperty( "DoNeutronInTrackerOnly",      m_doNeutronInTrackerOnly     = false);
  declareProperty( "MuonCylindarCut",             m_muonAxisCylinderRadius     = 100.0 * CLHEP::mm );

  declareProperty( "MaxNeutronClusters",          m_maxRecoCluster             = 1000 );
  declareProperty( "ConsiderNCprotons",           m_considerNCprotons          = false );


  //-- get the lattice energies vector? - default to no
  declareProperty( "UseLatticeEnergies", m_getLatticeEnergies = false );
  //-- get the string for ML prediction file
  declareProperty("MLPredFile", m_getMLPredFilename = "MLPredFile.root" );
  declareProperty("UseDNNVtx", m_useDNN = false );


  PART_TYPE_MAP.clear();
  PART_TYPE_MAP[MCPartType::klown] = "lown";
  PART_TYPE_MAP[MCPartType::kmidn] = "midn";
  PART_TYPE_MAP[MCPartType::khighn] = "highn";
  PART_TYPE_MAP[MCPartType::kp] = "p";
  PART_TYPE_MAP[MCPartType::kmeson] = "meson";
  PART_TYPE_MAP[MCPartType::kmu] = "mu";
  PART_TYPE_MAP[MCPartType::kem] = "em";
  PART_TYPE_MAP[MCPartType::kXtalk] = "xtalk";
  PART_TYPE_MAP[MCPartType::kOther] = "other";
  //Vetowall properties
  declareProperty( "VetoEfficiencyName"       ,  m_getVetoEffName         = "GetVetoEfficiency" );
  declareProperty("TrackPropagatorToVeto",       m_propagateToVetoName     = "TrackLinearPropagator");
  declareProperty( "RecoObjectTimeToolName"   , m_timeToolName            = "RecoObjectTimeTool" );
  declareProperty( "PanelOffset",           m_PanelOffset         = 55.0 * CLHEP::mm   );
  declareProperty( "VetoWallIsImportant",          m_VetoWallIsImportant             = true );
  declareProperty("NumVetoPaddles"                  , m_NumVetoPaddles = 6);
  declareProperty("VetoPaddlePMT1"                  , m_VetoPaddlePMT1 = 1);
  declareProperty("VetoPaddlePMT2"                  , m_VetoPaddlePMT2 = 2);
  declareProperty("TotalNumberofPMTsforVETOWALL"    , m_totalVETOPMTS = 24);

}

//=============================================================================
// Initialize
//=============================================================================

StatusCode MasterAnaDev::initialize() {

  debug()<<"Entering MasterAnaDev::initialize() ....."<<endmsg;


  //Selecting fiducial and analizable modules by Oscar
  debug()<<"DSCal region enabled? "<<m_useDSCal<<endmsg;
  m_useDSCal == true ? m_dsFiducialMod = 114 : m_dsFiducialMod = 80;
  m_useDSCal == true ? m_dsAnalyzableMod = 114 : m_dsAnalyzableMod = 94;
  m_useDSCal == true ? MAX_TRACKER_MODULE = 114 : MAX_TRACKER_MODULE = 84;


  StatusCode sc = this->MinervaAnalysisTool::initialize();
  if( sc.isFailure() ) { return Error( "Failed to initialize!", sc ); }

  m_idDet = svc<IGeomUtilSvc>("GeomUtilSvc", true)->getIDDet();
  if( 0 == m_idDet )
    {
      fatal() << "Failed to retrieve InnerDetector" << endmsg;
      return StatusCode::FAILURE;
    }


  m_geomUtilSvc = svc<IGeomUtilSvc>("GeomUtilSvc");

  m_randomGen = new TRandom3( m_randomSeed );
  if ( m_randomGen == NULL ) return StatusCode::FAILURE;

  m_fidHexApothem  = 850.0*CLHEP::mm;
  m_fidUpStreamZ   = 5990.0*CLHEP::mm;   // ~middle of module 27, plane 1
  m_fidDownStreamZ = 8340.0*CLHEP::mm;   // ~middle of module 79, plane 1

  //for Machine Learning Zubair
  //Get DB File prediction
  if( m_useDNN ){
    filename   = m_getMLPredFilename;
    debug() << "ML Prediction file: " << filename << endmsg;
    debug() << "m_useDNN = " << m_useDNN << endmsg;
    DBPredFile = new TFile(filename.c_str());
    if(DBPredFile == NULL){fatal()<<"Failed to fail ML file ="<<  filename  <<"CHECK OPTION FILE!!" <<endmsg; }
    dbPred     = (TTree*)DBPredFile->Get("MasterAnaDev");
    if(dbPred == NULL){fatal()<<"Failed to fail Get correct Tree Check option for NukeCC-> CHECK OPTION FILE!! " <<endmsg; }
  }



  //convert module numbers to moduleID's
  m_usFiducialModID = getModuleID( m_usFiducialMod );
  m_dsFiducialModID = getModuleID( m_dsFiducialMod );
  m_usAnalyzableModID = getModuleID( m_usAnalyzableMod );
  m_dsAnalyzableModID = getModuleID( m_dsAnalyzableMod );
  //NuclearTargetTool to create the NuclearTargets.
  // This way, I add the NuclearTargets I want.
  try
    {
      m_nukeTool = tool<INuclearTargetTool>("NuclearTargetTool", m_nukeToolAlias);//will be shared
      if( m_addMyNuclearTargets )
        {
          debug() << "telling " << m_nukeToolAlias << " to add NuclearTargets" << endmsg;

          // Add passive targets
          m_nukeTool->addAllPassiveNuclearTargets();
        }
      else
        debug() << "Trusting NuclearTargetTool to add all the targets I need." << endmsg;
      debug() << "Locking NuclearTargetTool " << m_nukeToolAlias << endmsg;
      m_nukeTool->lock();

      m_targets = m_nukeTool->getNuclearTargets();
    }
  catch( GaudiException& e )
    {
      error() << "Could not obtain NuclearTargetTool" << endmsg;
      return StatusCode::FAILURE;
    }

  if( fillPlaneLowZMap().isFailure() ) return StatusCode::FAILURE;

  //Get MinervaObjectAssociator
  try{
    m_objectAssociator = tool<IMinervaObjectAssociator>("MinervaObjectAssociator", m_minObjAssocAlias);
  } catch( GaudiException& e ) {
    error() << "Could not obtain tool: " << m_minObjAssocAlias << endmsg;
    return StatusCode::FAILURE;
  }

  //Get RecoObjectTime Tool
  try{
    m_recoObjTimeTool = tool<IRecoObjectTimeTool>("RecoObjectTimeTool", m_recoObjTimeToolAlias);
  } catch(GaudiException& e){
    error()<<"Could not obtain tool: " << m_recoObjTimeToolAlias << endmsg;
    return StatusCode::FAILURE;
  }

  //Get MinervaCoordSys Tool
  try{
    m_minCoordSysTool = tool<IMinervaCoordSysTool>("MinervaCoordSysTool", m_minCoordSysToolAlias);
  } catch( GaudiException& e ) {
    error() << "Could not obtain tool: " << m_minCoordSysToolAlias << endmsg;
    return StatusCode::FAILURE;
  }
  try{
    m_vtxBlobCreator = tool<IIDAnchoredBlobCreator>("VertexBlobCreator", "NukeCCVertexBlobCreator");
    debug() << "Got VertexBlobCreator named: " << m_vtxBlobCreator->name() << endmsg;
  }
  catch( GaudiException& e )
    {
      error() << "Could not obtain tool: VertexBlobCreator" << endmsg;
      return StatusCode::FAILURE;
    }
  //////////////////////////////////////////////////////////////////////////


  //Get MinervaMath Tool
  try{
    m_mathTool = tool<IMinervaMathTool>("MinervaMathTool", m_mathToolAlias);
  } catch( GaudiException& e ) {
    error() << "Could not obtain tool: " << m_mathToolAlias << endmsg;
    return StatusCode::FAILURE;
  }

  // Path length tool for momentum by range
  try {
        m_pathLengthTool = tool<IPathLengthTool>("PathLengthTool","PathLengthTool" );
                 } catch( GaudiException& e ){
                       error() << "Could not obtain tool: " << "PathLengthTool" << endmsg;
                             return StatusCode::FAILURE;
                                 }
  try {
      m_LikelihoodPIDTool = tool<ILikelihoodParticleTool>("LikelihoodPIDTool", m_LPIDToolAlias);
    } catch( GaudiException& e ){
      error() << "Could not obtain LikelihoodPIDTool: " << m_LPIDToolAlias << endmsg;
      return StatusCode::FAILURE;
    }

  //Get MASTERANADEVRecoUtils Tool
  try{
    m_masterAnaDevRecoUtils = tool<IMasterAnaDevRecoUtils>("MasterAnaDevRecoUtils", m_masterAnaDevRecoUtilsAlias);
  } catch(GaudiException& e) {
    error() << "Could not obtain tool: " << m_masterAnaDevRecoUtilsAlias <<endmsg;
    return StatusCode::FAILURE;
  }

  //Get BlobCreatorUtm_masterAnaDevRecoUtilsls Tool
  try{
    m_blobCreatorUtils = tool<IBlobCreatorUtils>("BlobCreatorUtils", m_blobCreatorUtilsAlias);
  } catch( GaudiException& e ){
    error() << "Could not obtain BlobCreatorUtils: " << m_blobCreatorUtilsAlias << endmsg;
    return StatusCode::FAILURE;
  }

  //Get PrimaryBlobProng Tool
  try{
    m_primaryBlobProngTool = tool<IPrimaryBlobProngTool>("PrimaryBlobProngTool", m_primaryBlobProngToolAlias);
  } catch( GaudiException& e ) {
    error() << "Could not obtain tool: " << m_primaryBlobProngToolAlias << endmsg;
    return StatusCode::FAILURE;
  }

  //for Muon Energy Correction by Huma
  try
    {
      m_energyCorrectionTool = tool<IEnergyCorrectionTool>("EnergyCorrectionTool");
      debug() << " Got EnergyCorrectionTool named: " << m_energyCorrectionTool->name() << endmsg;
    }
  catch( GaudiException& e )
    {
      error() << "Could not obtain EnergyCorrectionTool!" << endmsg;
      return StatusCode::FAILURE;
    }

  //Get Michel Tool
  try{
    m_michelTool = tool<IMichelTool>("MichelTool", m_michelToolAlias);
  } catch(GaudiException& e){
    error()<<"Could not obtain tool: " << m_michelToolAlias << endmsg;
    return StatusCode::FAILURE;
  }

  try {
    m_MCTrackTool = tool<IMCTrackTool>("MCTrackTool");
  } catch (GaudiException& e) {
    error() << "Could not obtain tool: MCTrackTool" << endmsg;
    return StatusCode::FAILURE;
  }

  try {
    m_InnerDetector = getDet<Minerva::DeDetector>("/dd/Structure/Minerva/Detector/InnerDetector");
  } catch(GaudiException& e) {
    error() << "Could not obtain Inner Detector" << endmsg;
    return StatusCode::FAILURE;
  }

  //Get Improved Michel Tool //ImprovedMichelTool
  try{
    m_improvedmichelTool = tool<IImprovedMichelTool>("ImprovedMichelTool");
  } catch (GaudiException& e) {
    error() << "Could not obtain tool: ImprovedMichelTool" << endmsg;
    return StatusCode::FAILURE;
  }

  //Get Calorimetry Tool
  try{
    m_caloUtils = tool<ICalorimetryUtils>("CalorimetryUtils", m_caloUtilsAlias);
  } catch( GaudiException& e ) {
    error() << "Could not obtain tool: " << m_caloUtilsAlias << endmsg;
    return StatusCode::FAILURE;
  }

  //Get ProtonUtils Tool
  try{
    m_protonUtils = tool<IProtonUtils>("ProtonUtils", m_protonUtilsAlias);
  } catch( GaudiException& e ){
    error() << "Could not obtain tool: " << m_protonUtilsAlias << endmsg;
    return StatusCode::FAILURE;
  }

  //Get ParticleMaker Tool
  try{
    m_particleMaker = tool<IParticleMakerTool>("ParticleMakerTool", m_particleMakerAlias);
  } catch( GaudiException& e ){
    error() << "Could not obtain tool: " << m_particleMakerAlias << endmsg;
    return StatusCode::FAILURE;
  }

  //Get HitTagger Tool
  try{
    m_hitTaggerTool = tool<IHitTaggerTool>("HitTaggerTool", m_hitTaggerToolAlias);
  } catch( GaudiException& e ) {
    error() << "Could not obtain tool: " << m_hitTaggerToolAlias << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_doNeutron )
    {
      //Get Neutron Blob tool
      try{
        m_neutronBlobRecoTool = tool<INeutronBlobRecoTool>(m_neutronBlobRecoToolAlias);
      } catch (GaudiException& e ) {
        error() << "Could not obtain tool: "<< m_neutronBlobRecoToolAlias<< endmsg;
        return StatusCode::FAILURE;
      }
      try{
        m_neutronBlobUtils = tool<INeutronBlobUtils>(m_neutronBlobUtilsAlias);
      } catch (GaudiException& e ) {
        error() << "Could not obtain tool: "<< m_neutronBlobUtilsAlias<< endmsg;
        return StatusCode::FAILURE;
      }


      try {
        m_anaBlobUtils = tool<IAnaBlobUtils>(m_anaBlobUtilsAlias);
      } catch( GaudiException& e ) {
        error() << "Could not get tool: "<< m_anaBlobUtilsAlias<< endmsg;
        return StatusCode::FAILURE;
      }
      try {
        m_blobFitTool = tool<IBlobFitTool>("BlobFitTool", m_blobFitToolAlias);
      } catch( GaudiException& e ) {
        error() << "Could not get tool: "<< m_blobFitToolAlias<< endmsg;
        return StatusCode::FAILURE;
      }
    }

  try
  {
    m_neutronInelasticReweight = tool<INeutronInelasticReweightTool>("NeutronInelasticReweightTool");
  }
  catch(const GaudiException& e)
  {
    error() << "Failed to get a handle to a NeutronInelasticReweightTool tool with "
            << "exception " << e.what() << endmsg;
    return e.code();
  }

  //added tool to get PMT status from DB//christian

  //  m_VetoDet = NULL;
  try{
    m_VetoDet = getDet<Minerva::DeVetoDetector>("/dd/Structure/Minerva/VetoDetector");
  }catch(GaudiException& e){
    error() << "Failed to retrieve Veto Wall Geometry. Perhaps this is the Frozen Detector?" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hasVeto = m_geomUtilSvc->hasVeto();

  if( !m_hasVeto ) {return StatusCode::SUCCESS;}
  else{warning() << "Failed to get Veto geometry" << endmsg;}

  try {
    m_getVetoEff = tool<IGetVetoEfficiency>( m_getVetoEffName );
  } catch(GaudiException& e){
    error() << "Could not obtain tool: " << m_getVetoEffName << endmsg;
    return StatusCode::FAILURE;
  }

  try{
    m_timeTool = tool<IRecoObjectTimeTool>( m_timeToolName );
  }catch(GaudiException& e){
    error()<<"Could not obtain tool: " << m_timeToolName << endmsg;
    return StatusCode::FAILURE;
  }

  //added to extrapolate track back to veto wall//christian
  try {
    m_propagateToVeto = tool<ITrackLinearPropagator>( m_propagateToVetoName );
  } catch(GaudiException& e){
    error() << "Could not obtain tool: " << m_propagateToVetoName << endmsg;
    return StatusCode::FAILURE;
  }

  try{ m_plexModel = tool<IPlexModel>("MinervaPlexModel");
  } catch(GaudiException& e){
    error() <<"Could not obtain a Plex tool" << endmsg;
    return StatusCode::FAILURE;
  }

  //added tool CryoTarget info//christian

  try{
    m_CryoDet = getDet<Minerva::DeCryoTarget>("/dd/Structure/Minerva/CryoTarget");
  } catch(...){
    warning() << "Failed to retrieve Cryo Target Geometry. Perhaps this is Frozen Detector data?" << endmsg;
  }
  m_hasCryo = m_geomUtilSvc->hasCryo();
  if( !m_hasCryo ) {warning() << "Failed to get Cryo-Target geometry" << endmsg;}
  else{}

  try{ m_electronUtils = tool<IElectronUtils>("ElectronUtils");
  } catch(GaudiException& e){
    error() <<"Could not obtain ElectronUtils tool" << endmsg;
    return StatusCode::FAILURE;
  }

  m_ccUtils              = tool<IAnalysisToolUtils>( "CCInclusiveUtils" );
  m_ccUtilsWideTimeWindow= tool<IAnalysisToolUtils>( "CCInclusiveUtils", "CCUtilsWideWindow" );
  m_recoilReconstructor  = tool<IRecoilReconstructor>( "RecoilReconstructor" );
  m_caloUtils            = tool<ICalorimetryUtils>( "CalorimetryUtils", "RecoilCalorimetryUtils" );  //same one as RecoilCalorimetryUtils
  //=========================================


  //for Machine Learning
  if( m_getLatticeEnergies ){
    try {
      m_MLVFTool = tool<IMLVFTool>("MLVFTool", m_MLVFToolAlias);
      debug() << " Got MLVF Sample Prep Tool named: " << m_MLVFTool->name() << endmsg;
    }
    catch( GaudiException& e ) {
      error() << "Could not obtain MLVF Sample Prep tool: " << m_MLVFToolAlias << endmsg;
      return StatusCode::FAILURE;
    }
  } // get lattice energies

  //if not using many vtx blob radii, keep only the first one (CV)
  if( !m_doVtxBlobSystematic && m_antinuVtxBlobRadii.size() > 1 )
    {
      m_nuVtxBlobRadii.erase(     m_nuVtxBlobRadii.begin()+1,     m_nuVtxBlobRadii.end() );
      m_antinuVtxBlobRadii.erase( m_antinuVtxBlobRadii.begin()+1, m_antinuVtxBlobRadii.end() );
    }

  ///// Get ElasticFilterTool
  try{
    m_elasticFilterTool = tool<IElasticFilterTool>("ElasticFilterTool", m_elasticFilterToolAlias);
  } catch( GaudiException& e ) {
    error() << "Could not obtain tool: " << m_elasticFilterToolAlias << endmsg;
    return StatusCode::FAILURE;
  }


  // grouped with branch fill calls
  //---------------------------------------------------------------------------
  //! Declare common ntuple branches that will be present in your AnaTuple
  //---------------------------------------------------------------------------
  declareCommonPhysicsAnaBranches();
  declareNuMIBranches();
  declareMinosMuonBranches();
  declareGenieWeightBranches();
  //declareVertexActivityStudyBranches();
  declareMatchedMichelBranches();
  declareSystematicShiftsBranches();
  declareParticleResponseBranches();
  declareModifiedParticleResponseBranches();
  // For HadronReweightTool
  declareHadronReweightBranches();
  declareInelasticReweightBranches();
  //declareBoolEventBranch( getPassTag() );                     //< isCC && is_fiducial && targetCode>0
  //---------------------------------------
  //! Declare branches for own ntuples
  //---------------------------------------
  declareBoolTruthBranch( getPassTag() );                   //< isCC && is_fiducial && targetCode>0

  //! Event - Truth
  declareBoolTruthBranch( "is_fiducial" );
  declareIntTruthBranch( "reco_is_fiducial", -1 );
  declareIntTruthBranch( "reco_is_minos_match", -1 );
  declareIntTruthBranch( "reco_has_muon_charge", -1 );
  declareIntTruthBranch( "reco_muon_is_minos_match_track", -1 );
  declareIntTruthBranch( "reco_muon_is_minos_match_stub", -1 );
  declareIntTruthBranch( "reco_pass_MasterAnaDev_precuts", 0 );
  declareIntTruthBranch( "reco_has_single_proton", 0 );

  declareBoolTruthBranch( "reco_hasGoodObjects");
  declareBoolTruthBranch( "reco_isGoodVertex");
  declareBoolTruthBranch( "reco_isFidVol_smeared");

  declareIntTruthBranch("N_pip",        -1 );
  declareIntTruthBranch("N_chargedK",   -1 );
  declareIntTruthBranch("N_lambda",     -1 );
  declareIntTruthBranch("N_K0",         -1 );
  declareIntTruthBranch("N_sigma",      -1 );
  declareIntTruthBranch("N_pim",        -1 );
  declareIntTruthBranch("N_pi0",        -1 );

  declareDoubleTruthBranch("muon_px", 0 );
  declareDoubleTruthBranch("muon_py", 0 );
  declareDoubleTruthBranch("muon_pz", 0 );
  declareDoubleTruthBranch("muon_E", -9.0 );

  declareContainerDoubleTruthBranch("pi_theta_wrtbeam",  20, -9.0 );
  declareContainerDoubleTruthBranch("pi_E",  20, -9.0 );
  declareContainerDoubleTruthBranch("pi_px",  20, -9.0 );
  declareContainerDoubleTruthBranch("pi_py",  20, -9.0 );
  declareContainerDoubleTruthBranch("pi_pz",  20, -9.0 );

  declareContainerIntTruthBranch("pi_charge", 20, 0 );
  // declareIntTruthBranch( "reco_has_michel_electron", 0 );
  declareBoolTruthBranch( "reco_has_muon" );                  //< Is there a muon in this PhysicsEvent?
  //! Event - general reco
  declareBoolEventBranch( getPassTag() );                     //< isCC && is_fiducial && targetCode>0
  declareIntEventBranch( "multiplicity", -1 );

  //for active targets, record Z of the material at that vtx_XY in each passive target
  declareContainerIntBranch( m_hypMeths, "ref_targZ", 5, -1 );
  declareContainerIntBranch( m_hypMeths, "ANN_ref_targZ", 5, -1 );
  //for active targets,record distance in XY to the nearest edge of a material if this XY were in each passive target
  declareContainerDoubleBranch( m_hypMeths, "ref_dist_to_division", 5, -1 );
  declareContainerDoubleBranch( m_hypMeths, "ANN_ref_dist_to_division", 5, -1 );
  //< distance in Z from vertex to center of the passive targets
  declareContainerDoubleBranch( m_hypMeths, "ref_dist_to_target", 5, -1 );
  declareContainerDoubleBranch( m_hypMeths, "ANN_ref_dist_to_target", 5, -1 );
  //< How far is the vertex from the nearest edge of the NuclearTarget in Z?
  declareDoubleBranch( m_hypMeths, "target_zDist", -1000. );
  declareDoubleBranch( m_hypMeths, "ANN_target_zDist", -1000. );
  //< How far is the vertex from the nearest edge of a DeTargetSection in the NuclearTarget in XY?
  declareDoubleBranch( m_hypMeths, "target_dist_to_division", -1000.);
  declareDoubleBranch( m_hypMeths, "ANN_target_dist_to_division", -1000.);
  declareIntBranch( m_hypMeths, "in_fiducial_area", 0 );//< Is the vertex in the region of interest in XY
  declareIntBranch( m_hypMeths, "ANN_in_fiducial_area", 0 );//< Is the vertex in the region of interest in XY
  declareIntBranch( m_hypMeths, "ANN_in_analyzable_area", 0 );//< Is the vertex in the region of interest in XY
  declareIntBranch( m_hypMeths, "vtx_module", -999 );//< What module is this vertex in (if appropriate)
  declareIntBranch( m_hypMeths, "vtx_plane",  -999 );//< What plane is this vertex in (if appropriate)
  declareIntTruthBranch( "vtx_module", -999 );     //< What module is the true vertex in
  declareIntTruthBranch( "vtx_plane",  -999 );     //< What plane is the true vertex in
  //< How far is the vertex from the nearest edge of the NuclearTarget in Z?
  declareDoubleTruthBranch( "target_zDist", -1000. );
  //< How far is the vertex from the nearest edge of a DeTargetSection in the NuclearTarget in XY?
  declareDoubleTruthBranch( "target_dist_to_division", -1000.);

  declareDoubleTruthBranch( "muon_theta",  -999.);

  declareIntEventBranch( "recoil_EInc", -1 );
  declareIntTruthBranch( "target_code", 0); //< This is targetID*1000 + targetZ.  Occasionally useful
  declareIntTruthBranch( "targetID", 0 );   //< True event is in nuclear target with this ID
  declareIntTruthBranch( "targetZ", 0 );    //< True event is in a nuclear target material with this Z
  declareIntTruthBranch( "in_fiducial_area", 0); //< Is the true vertex in the region of interest in XY
  declareBoolTruthBranch( "reco_has_int_vtx" );
  declareIntEventBranch( "has_single_proton", -1 );
  declareIntEventBranch( "pass_MasterAnaDev_precuts", -1 );
  declareIntEventBranch( "n_iso_trk_prongs", -1 );
  declareIntEventBranch( "n_anchored_short_trk_prongs", -1 );
  declareIntEventBranch( "n_anchored_long_trk_prongs", -1 );
  declareContainerIntEventBranch( "event_vertex_types" );
  declareContainerIntEventBranch( "event_in_time_vertex_types" );
  declareContainerDoubleEventBranch( "event_vertex_time_diff" );
  declareContainerDoubleEventBranch( "event_track_time_diff" );
  declareContainerDoubleEventBranch( "event_tracks_energy" );
  declareContainerDoubleEventBranch( "event_in_time_tracks_energy" );
  declareContainerDoubleEventBranch( "event_track_hit_energy" );
  declareContainerDoubleEventBranch( "event_track_hit_time" );
  declareContainerIntEventBranch( "event_track_nhits" );
  declareIntEventBranch( "ID_Hits_Post_Long", -1);
  declareIntEventBranch( "slice_n_hits", -1);
  declareContainerDoubleEventBranch( "slice_hit_time");
  declareContainerDoubleEventBranch( "slice_hit_energy");
  declareContainerDoubleEventBranch("all_event_start_vertex_time");
  declareContainerDoubleEventBranch("all_event_start_vertex_time_minos_match");
  declareContainerIntEventBranch("all_event_start_vertex_fv_minos_match");
  declareIntEventBranch( "has_interaction_vertex", -1);
  declareIntEventBranch( "n_minos_matches", -1);
  declareContainerIntEventBranch("recoil_lower_time_limit");
  declareContainerIntEventBranch("recoil_upper_time_limit");
  declareContainerDoubleEventBranch("recoil_summed_energy");
  declareContainerDoubleEventBranch("recoil_summed_energy_edge");
  declareContainerDoubleEventBranch("recoil_data_fraction");
  declareIntEventBranch( "event_used_start_vertices",-1);
  declareIntEventBranch( "event_unused_start_vertices",-1);
  declareContainerDoubleEventBranch( "event_extra_track_PID");//will store pid values from tracks not in the interaction vertex.
  declareBoolEventBranch( "isMinosMatchTrack");
  //! Event - reco Muon
  declareIntEventBranch( "muon_is_minos_match_track", -1 );
  declareIntEventBranch( "muon_is_minos_match_stub", -1 );
  declareIntEventBranch( "isMinosMatchTrack", -1 );
  declareIntEventBranch( "isMinosMatchStub", -1 );
  declareDoubleEventBranch( "muon_minerva_trk_chi2PerDoF", -1 );
  declareContainerDoubleEventBranch( "muon_theta_allNodes" );
  declareContainerDoubleEventBranch( "muon_thetaX_allNodes" );
  declareContainerDoubleEventBranch( "muon_thetaY_allNodes" );

  //! Event - Muon Iso Blobs/Fuzz
  declareIntEventBranch("n_muon_iso_blobs", -1 );
  declareDoubleEventBranch("muon_iso_blobs_energy", -9999.0);
  declareIntEventBranch("muon_fuzz_energy", -9999.0);

  //! Event - Blobs & energies
  declareIntEventBranch( "n_nonvtx_iso_blobs", -1 );
  declareIntEventBranch( "n_nonvtx_iso_blobs_all", -1 );
  declareDoubleEventBranch( "nonvtx_iso_blobs_energy" , -9999.0 );
  declareDoubleEventBranch( "nonvtx_iso_blobs_energy_all" , -9999.0 );
  declareDoubleEventBranch( "vtx_blobs_energy" , -9999.0 );
  declareDoubleEventBranch( "dis_id_energy", -9999.0 );
  declareDoubleEventBranch( "vtx_iso_blobs_energy_outside_radius", -9999.0 );
  declareDoubleEventBranch( "recoil_energy_nonmuon_nonvtx0mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_nonvtx50mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_nonvtx100mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_nonvtx150mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_nonvtx200mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_nonvtx250mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_nonvtx300mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_vtx0mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_vtx50mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_vtx100mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_vtx150mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_vtx200mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_vtx250mm", -9999.0);
  declareDoubleEventBranch( "recoil_energy_nonmuon_vtx300mm", -9999.0);
  declareContainerDoubleEventBranch( "vtx_blobs_vtx_energy_in_prong" );
  declareContainerDoubleEventBranch( "vtx_blobs_iso_energy_in_prong" );
  declareContainerDoubleEventBranch( "vtx_blobs_iso_energy_clusters_outside_radius_in_prong" );
  declareContainerDoubleEventBranch( "vtx_blobs_vtx_distance_in_prong" );
  declareContainerDoubleEventBranch( "vtx_blobs_iso_distance_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_energy_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_distance_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_start_position_x_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_start_position_y_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_start_position_z_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_time_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_time_difference_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_lowest_module_x_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_lowest_module_u_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_lowest_module_v_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_highest_module_x_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_highest_module_u_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_highest_module_v_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_earliest_hit_time_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_latest_hit_time_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_highest_hit_energy_in_prong" );
  declareContainerIntEventBranch( "nonvtx_iso_blobs_n_hits_in_prong" );
  declareContainerIntEventBranch( "nonvtx_iso_blobs_particle_pdg_in_prong" );
  declareContainerIntEventBranch( "nonvtx_iso_blobs_primary_particle_pdg_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_matched_energy_fraction_in_prong" );
  declareContainerDoubleEventBranch( "nonvtx_iso_blobs_data_energy_fraction_in_prong" );


  //! Event - Off-proton track cluster info
  declareContainerIntEventBranch( "clusters_found_at_end_proton_prong" );
  declareContainerDoubleEventBranch( "clusters_found_at_end_proton_prong_max_distance" );
  declareContainerIntEventBranch( "number_clusters_at_end_proton_prong" );
  declareContainerDoubleEventBranch( "visE_clusters_at_end_proton_prong" );
  declareContainerDoubleEventBranch( "calibE_clusters_at_end_proton_prong" );


  //declareIntEventBranch( "truth_has_michel_electron", 0 );
  declareIntEventBranch( "truth_improved_michel_electron", 0 );
  //declareContainerDoubleEventBranch( "truth_has_michel_from_pion_plus_momentum" );
  //declareContainerDoubleEventBranch( "truth_has_michel_from_pion_minus_momentum" );

  //! Event - Improved Michel Tool variables //ImprovedMichelTool
  declareIntEventBranch("improved_nmichel");
  declareContainerIntEventBranch("improved_michel_vertex_type");
  declareContainerIntEventBranch("improved_michel_in_vertex_point");

  declareContainerIntEventBranch("improved_michel_match_vec");
  declareContainerDoubleEventBranch("improved_michel_tvec");
  declareContainerDoubleEventBranch("improved_michel_tdiff_vec");
  declareContainerDoubleEventBranch("improved_michel_xvec");
  declareContainerDoubleEventBranch("improved_michel_yvec");
  declareContainerDoubleEventBranch("improved_michel_uvec");
  declareContainerDoubleEventBranch("improved_michel_vvec");
  declareContainerDoubleEventBranch("improved_michel_zvec");
  declareContainerDoubleEventBranch("improved_michel_dist_vec");
  declareContainerDoubleEventBranch("improved_michel_evis_vec");
  declareContainerDoubleEventBranch("improved_michel_ecalo_vec");
  declareContainerIntEventBranch("improved_michel_view_vec");
  declareContainerDoubleEventBranch("improved_michel_hit_charges");
  declareContainerDoubleEventBranch("improved_michel_hit_times");
  declareContainerDoubleEventBranch("improved_michel_hit_time_diff_vtx");
  declareContainerDoubleEventBranch("improved_michel_hit_time_diff_cluster");
  declareContainerIntEventBranch("improved_michel_ndigits");

  declareContainerIntEventBranch("improved_michel_matched_pdg");
  declareContainerIntEventBranch("improved_michel_matched_primary_pdg");
  declareContainerDoubleEventBranch("improved_michel_matched_energy_fraction");
  declareContainerDoubleEventBranch("improved_michel_data_energy_fraction");

  //! Event - Proton truth matcher
  declareIntEventBranch("proton_prong_PDG", -1);
  declareIntEventBranch("proton_prong_traj_ID", -1);
  declareContainerDoubleEventBranch("proton_prong_4p", 4, -1.0);
  declareContainerDoubleEventBranch("proton_prong_tpos", 4, -1.0);

  //for Machine Learning
  declareDoubleEventBranch( "muon_trackVertexTime", -9999. );  //needed by ML m_MLVFTool

  declareContainerIntEventBranch("sec_protons_prong_PDG");
  declareContainerIntEventBranch("sec_protons_prong_traj_ID");
  declareContainerDoubleEventBranch("seco_protons_prong_4p_E");
  declareContainerDoubleEventBranch("seco_protons_prong_4p_px");
  declareContainerDoubleEventBranch("seco_protons_prong_4p_py");
  declareContainerDoubleEventBranch("seco_protons_prong_4p_pz");
  declareContainerDoubleEventBranch("proton_prong_tpos_x");
  declareContainerDoubleEventBranch("proton_prong_tpos_y");
  declareContainerDoubleEventBranch("proton_prong_tpos_z");
  declareContainerDoubleEventBranch("proton_prong_tpos_t");

  // Momentum-by-range hadron variables
  declareContainerIntEventBranch("hadron_ntrack",                    10, -1 );
  declareContainerDoubleEventBranch("hadron_track_length",           10, -1.);
  declareContainerDoubleEventBranch("hadron_track_length_area",      10, -1.);
  declareContainerDoubleEventBranch("hadron_track_length_area_trkr", 10, -1.);
  declareContainerDoubleEventBranch("hadron_track_length_area_ecal", 10, -1.);


  //==========================================
  //add blob branches
  //==========================================

  // Removing the following branch declarations as there is no information saved to them,   
  // and the code related is modifying cluster histories and affecting other reconstruction. 
  // Can be added back in by uncommenting if they are used more carefully than at present.
  // -David Last 01/20/2022

  /*
  declareBlobBranches( VTXBLOB_PREFIX );
  declareBlobBranches( MUONFUZZBLOB_PREFIX );
  declareBlobBranches( ISOBLOB_PREFIX );
  declareBlobBranches( DISPBLOB_PREFIX );
  */

  declareBlobBranches( RECOILBLOB_PREFIX );
  declareBlobBranches( NUEFUZZBLOB_PREFIX );

  declareDoubleTruthBranch("visibleE", 0 );


  //for each vtx radius, we record ccqe-recoil
  for( unsigned int i = 0; i != m_nuVtxBlobRadii.size(); ++i )
    {
      string suffix = (0==i) ? "" : Form("_alt%02d",i);
      declareDoubleEventBranch( "blob_ccqe_recoil_E" + suffix, 0. );
    }

  //and the cc inclusive blob braches
  declareContainerDoubleEventBranch( "vtx_blob_radius", m_nuVtxBlobRadii.size(), 0. );

  /////////////////////////////////////////////////////////////////////////////////////////
  //-----------------------------------------------------
  //! Veto Wall Branches
  //-----------------------------------------------------
  declareIntEventBranch("VetoWall_NumberMatchToVeto", 0);//number of tracks in vtx options
  declareContainerDoubleEventBranch("VetoWall_muon_extrapVetoXY", 4,     -999);
  declareContainerIntEventBranch("VetoWall_PMTStatusMAP", 24,     0);//adding pandel configuration
  declareIntEventBranch("VetoWall_event_IsVeto", -9999);
  declareIntEventBranch("VetoWall_sixPush", -1);
  declareContainerIntEventBranch("VetoWall_sixPush_paddle_wall1");
  declareContainerIntEventBranch("VetoWall_sixPush_paddle_wall2");
  declareBoolEventBranch("VetoWall_extrapnowalls");
  declareBoolEventBranch("VetoWall_extrapwall1");
  declareBoolEventBranch("VetoWall_extrapwall2");
  declareBoolEventBranch("VetoWall_extrapbothwalls");
  declareBoolEventBranch("VetoWall_MuonTrkMatchToVETOwalloff");


  //-----------------------------------------------------
  //! Neutrino CCQE Interaction Hypothesis branches
  //-----------------------------------------------------

  //! Neutrino Interaction - Muon
  declareDoubleBranch( m_hypMeths, "muon_theta",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "muon_E",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "muon_P",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "muon_Px",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "muon_Py",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "muon_Pz",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "muon_T",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "muon_score",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "Q2_CCQE",  -9999.0 ); //Q2 for CCQE analysis
  declareDoubleBranch( m_hypMeths, "Q2_Inclusive",  -9999.0 ); //Q2 for Inclusive analysis
  declareDoubleBranch( m_hypMeths, "Q2_wide_window", -9999.0);
  declareDoubleBranch( m_hypMeths, "muon_theta_biasUp", 0.0);
  declareDoubleBranch( m_hypMeths, "muon_theta_biasDown", 0.0);
  declareDoubleBranch( m_hypMeths, "muon_qp", 99.0);
  declareDoubleBranch( m_hypMeths, "muon_qpqpe", 99.0);


  declareIntBranch( m_hypMeths, "target_code", 0);   //< This is targetID*1000 + targetZ.  Occasionally useful
  declareIntBranch( m_hypMeths, "ANN_target_code", 0);   //< This is targetID*1000 + targetZ.  Occasionally useful
  declareIntBranch( m_hypMeths, "targetZ", 0 ); //< Event was reconstructed into a nuclear target material with this Z
  declareIntBranch( m_hypMeths, "ANN_targetZ", 0 ); //< Event was reconstructed into a nuclear target material with this Z
  declareIntBranch( m_hypMeths, "targetID", 0 );   //< Event was reconstructed in nuclear target with this ID
  declareIntBranch( m_hypMeths, "ANN_targetID", 0 );   //< Event was reconstructed in nuclear target with this ID
  //! Neutrino Interaction - Primary Proton
  declareDoubleBranch( m_hypMeths, "proton_E_fromdEdx",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_P_fromdEdx",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_T_fromdEdx",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_Px_fromdEdx",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_Py_fromdEdx",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_Pz_fromdEdx",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_theta_fromdEdx",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_calib_energy",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_score", -9999.0);
  declareDoubleBranch( m_hypMeths, "proton_score1", -9999.0);
  declareDoubleBranch( m_hypMeths, "proton_score2", -9999.0);

  declareDoubleEventBranch( "proton_track_length", -9999. );
  declareDoubleEventBranch( "proton_track_endx", -9999. );
  declareDoubleEventBranch( "proton_track_endy", -9999. );
  declareDoubleEventBranch( "proton_track_endz", -9999. );


  declareDoubleBranch( m_hypMeths, "proton_startPointX", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_startPointY", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_startPointZ", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_endPointX", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_endPointY", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_endPointZ", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_theta", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_thetaX", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_thetaY", -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_phi", -9999.0 );
  declareIntBranch( m_hypMeths, "proton_patternRec", -1);


  //! Neutrino Interaction - Secondary Protons
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_E_fromdEdx" );
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_P_fromdEdx" );
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_T_fromdEdx" );
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_Px_fromdEdx" );
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_Py_fromdEdx" );
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_Pz_fromdEdx" );
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_theta_fromdEdx" );
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_T_fromCalo" );
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_proton_scores");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_proton_scores1");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_proton_scores2");
  declareContainerIntBranch(m_hypMeths, "sec_protons_patternRec");
  //! ESC proton selection
  declareIntBranch( m_hypMeths, "proton_prongType");
  declareContainerDoubleBranch( m_hypMeths, "proton_nodes_nodesNormE");
  declareContainerDoubleBranch( m_hypMeths, "proton_nodes_E");

  declareContainerIntBranch( m_hypMeths, "sec_protons_prongType");
  declareContainerIntBranch( m_hypMeths, "sec_protons_nodes_index");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_nodes_nodesNormE");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_nodes_E");

  //! Node Cuts
  const int MAX_TRACKS = 10;
  declareContainerIntBranch(m_hypMeths,"pion_nNodes",MAX_TRACKS, -1);
  declareContainerDoubleBranch(m_hypMeths, "pion_lastnode_Q0", MAX_TRACKS, -1.);
  declareContainerDoubleBranch(m_hypMeths, "pion_lastnode_Q1", MAX_TRACKS, -1.);
  declareContainerDoubleBranch(m_hypMeths, "pion_lastnode_Q2", MAX_TRACKS, -1.);
  declareContainerDoubleBranch(m_hypMeths, "pion_lastnode_Q3", MAX_TRACKS, -1.);
  declareContainerDoubleBranch(m_hypMeths, "pion_lastnode_Q4", MAX_TRACKS, -1.);
  declareContainerDoubleBranch(m_hypMeths, "pion_lastnode_Q5", MAX_TRACKS, -1.);

  //! Neutrino Interaction - Primary Pion
  declareDoubleBranch( m_hypMeths, "pion_score", -9999.0);
  declareDoubleBranch( m_hypMeths, "pion_score1", -9999.0);
  declareDoubleBranch( m_hypMeths, "pion_score2", -9999.0);
  declareContainerDoubleBranch( m_hypMeths, "pion_E",                       10, 0);
  declareContainerDoubleBranch( m_hypMeths, "pion_T",                       10, 0);
  declareContainerDoubleBranch( m_hypMeths, "pion_P",                       10, 0);
  declareContainerDoubleBranch( m_hypMeths, "pion_Px",                      10, 0);
  declareContainerDoubleBranch( m_hypMeths, "pion_Py",                      10, 0);
  declareContainerDoubleBranch( m_hypMeths, "pion_Pz",                      10, 0);
  declareContainerDoubleBranch( m_hypMeths, "pion_theta",                   10, 0);
  declareContainerDoubleBranch( m_hypMeths, "pion_phi",                     10, 0);
  declareContainerDoubleBranch( m_hypMeths, "hadron_piFit_score1",         10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "hadron_piFit_scoreLLR", 10, -1.0);
  declareContainerIntBranch( m_hypMeths, "hadron_isTracker",   10, -1);
  declareContainerIntBranch( m_hypMeths, "hadron_isSideECAL",  10, -1);
  declareContainerIntBranch( m_hypMeths, "hadron_isExiting",   10, -1);
  declareContainerIntBranch( m_hypMeths, "hadron_isODMatch",   10, -1);
  declareContainerIntBranch( m_hypMeths, "hadron_isForked",    10, -1);
  declareIntBranch( m_hypMeths, "hadron_number", 0);
  declareContainerIntBranch( m_hypMeths, "hadron_tm_destructCode",       10, -3);
  declareContainerIntBranch( m_hypMeths, "hadron_tm_PDGCode",            10, 0);
  declareContainerIntBranch( m_hypMeths, "hadron_tm_trackID",            10, -1);
  declareContainerDoubleBranch( m_hypMeths, "hadron_tm_beginKE",                10, -1.0 );
  declareContainerDoubleBranch( m_hypMeths, "hadron_tm_fraction",               10, -1.0 );
  declareContainerDoubleBranch( m_hypMeths, "hadron_tm_otherE",                 10, -1.0 );
  declareContainerDoubleBranch( m_hypMeths, "hadron_tm_endpt_fraction",         10, -1.0 );
  declareContainerIntBranch( m_hypMeths, "hadron_tm_isTruthMatched",     10, 1);

  // Aaron's passive-corrected variables
  declareContainerDoubleBranch( m_hypMeths, "hadron_pion_E_corr",        10, -99.9);
  declareContainerDoubleBranch( m_hypMeths, "hadron_pion_p_corr",        10, -99.9);
  declareContainerDoubleBranch( m_hypMeths, "hadron_pion_px_corr",       10, -99.9);
  declareContainerDoubleBranch( m_hypMeths, "hadron_pion_py_corr",       10, -99.9);
  declareContainerDoubleBranch( m_hypMeths, "hadron_pion_pz_corr",       10, -99.9);
  declareContainerDoubleBranch( m_hypMeths, "hadron_pion_E_recoil_corr", 10, -99.9);

  //! Brandon Michel Properties
  declareContainerIntBranch( m_hypMeths,    "hadron_endMichel_category",           10, -1);
  declareContainerDoubleBranch( m_hypMeths, "hadron_endMichel_slice_energy",       10, -1.0);
  declareContainerIntBranch( m_hypMeths,    "hadron_endMichel_ndigits",            10, -1);
  declareContainerDoubleBranch( m_hypMeths, "hadron_endMichel_energy",             10, -1.0);

 //! Ben Analysis Michel Electron Branches
 // "misc"
 declareContainerIntEventBranch("matched_michel_code"   ,MAX_TRACKS, -1);

   // "matched_michel"
   declareContainerIntEventBranch("matched_michel_avg_idx",MAX_TRACKS, -1);
   declareContainerIntEventBranch("matched_michel_end_idx",MAX_TRACKS, -1);
   declareContainerIntEventBranch("matched_michel_ov_idx" ,MAX_TRACKS, -1);

   declareContainerDoubleEventBranch("matched_michel_fit_residual", MAX_TRACKS, -9999);

   declareContainerIntEventBranch("matched_michel_avg_module"      , MAX_TRACKS, -9999);
   declareContainerIntEventBranch("matched_michel_end_module1"      , MAX_TRACKS, -9999);
   declareContainerIntEventBranch("matched_michel_end_module2"      , MAX_TRACKS, -9999);

   // Michel truth matching
   declareContainerIntEventBranch    ( "has_michel_primary_particle_trackID");
   declareContainerIntEventBranch    ( "true_michel_pdg");
   declareContainerDoubleEventBranch ( "true_michel_x1");
   declareContainerDoubleEventBranch ( "true_michel_y1");
   declareContainerDoubleEventBranch ( "true_michel_z1");
   declareContainerDoubleEventBranch ( "true_michel_x2");
   declareContainerDoubleEventBranch ( "true_michel_y2");
   declareContainerDoubleEventBranch ( "true_michel_z2");
   declareContainerIntEventBranch    ( "michel_parent_PDG");
   declareContainerIntEventBranch    ( "michel_from_decay");
   declareContainerIntEventBranch    ( "michel_parent_trackID");
   declareContainerIntEventBranch    ( "has_michel_primary_particle_PDG");
   declareContainerDoubleEventBranch ( "has_michel_primary_particle_p");
   declareContainerDoubleEventBranch ( "has_michel_primary_particle_px");
   declareContainerDoubleEventBranch ( "has_michel_primary_particle_py");
   declareContainerDoubleEventBranch ( "has_michel_primary_particle_pz");
   declareContainerDoubleEventBranch ( "has_michel_primary_particle_E");

  //! Neutrino Interaction - Secondary Pions
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_pion_scores");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_pion_scores1");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_pion_scores2");

  //! dEdXTool Variations to Primary Pion
  declareContainerDoubleBranch( m_hypMeths, "piFit_score1_Mass_biasUp",         10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "piFit_score1_Mass_biasDown",       10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "piFit_score1_BetheBloch_biasUp",   10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "piFit_score1_BetheBloch_biasDown", 10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "piFit_score1_Birks",               10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "pion_E_Mass_biasUp", 10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "pion_E_Mass_biasDown", 10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "pion_E_BetheBloch_biasUp", 10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "pion_E_BetheBloch_biasDown", 10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "pion_E_Birks", 10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "pion_theta_biasUp", 10, -1.0);
  declareContainerDoubleBranch( m_hypMeths, "pion_theta_biasDown", 10, -1.0);

  // Prong cuts
  declareContainerDoubleEventBranch("primary_prong_separation", 10, -1.);

  declareContainerDoubleEventBranch("iso_prong_energy",     30, -1.);
  declareContainerDoubleEventBranch("iso_prong_separation", 30, -1.);
  declareContainerIntEventBranch("iso_prong_nhits",         30, -1);
  declareContainerDoubleEventBranch("iso_prong_x",          30, -1);
  declareContainerDoubleEventBranch("iso_prong_y",          30, -1);
  declareContainerDoubleEventBranch("iso_prong_z",          30, -1);
  declareIntEventBranch( "iso_prongs_count",                -1);
  declareIntEventBranch( "iso_prongs_blob_count",           -1);
  declareDoubleEventBranch( "iso_prongs_total_energy",      -1);

  //! Neutrino Interaction - Neutrino Energy
  declareDoubleBranch( m_hypMeths, "enu_muon",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "enu_proton",  -9999.0 );

  //! dEdXTool Variations to Primary Proton
  declareDoubleBranch( m_hypMeths, "proton_score1_Mass_biasUp",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_score1_Mass_biasDown",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_score1_BetheBloch_biasUp",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_score1_BetheBloch_biasDown",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_score1_MEU_biasUp",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_score1_MEU_biasDown",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_score1_Birks_bias",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_E_Mass_biasUp",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_E_Mass_biasDown",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_E_BetheBloch_biasUp",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_E_BetheBloch_biasDown",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_E_MEU_biasUp",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_E_MEU_biasDown",  -9999.0 );
  declareDoubleBranch( m_hypMeths, "proton_E_Birks_bias",  -9999.0 );

  //! dEdXTool Variations to Secondary Protons
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_score1_Mass_biasUp");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_score1_Mass_biasDown");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_score1_BetheBloch_biasUp");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_score1_BetheBloch_biasDown");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_score1_MEU_biasUp");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_score1_MEU_biasDown");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_score1_Birks_bias");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_E_Mass_biasUp");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_E_Mass_biasDown");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_E_BetheBloch_biasUp");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_E_BetheBloch_biasDown");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_E_MEU_biasUp");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_E_MEU_biasDown");
  declareContainerDoubleBranch( m_hypMeths, "sec_protons_E_Birks_bias");

  //! Adding recoil energy branches
  declareDoubleBranch( m_hypMeths, "recoil_E", -1.);     //< Calorimetric recoil energy
  declareDoubleBranch( m_hypMeths, "ANN_recoil_E", -1.);     //< Calorimetric recoil energy
  declareDoubleBranch( m_hypMeths, "visible_E", -1.);    //< Recoil Energy with no calorimetric correction at all
  declareDoubleBranch( m_hypMeths, "ANN_visible_E", -1.);    //< Recoil Energy with no calorimetric correction at all
  declareDoubleBranch( m_hypMeths, "recoil_passivecorrected", -1.); //< Recoil Energy with a passive matrl correction
  declareDoubleBranch( m_hypMeths, "recoil_E_wide_window", -1.); //< Calorimetric recoil energy using an extra wide time window
  declareDoubleBranch( m_hypMeths, "hadron_recoil_CCInc",         0.0);// Recoil energy CCInclusive Spline correction 
  declareDoubleBranch( m_hypMeths, "hadron_recoil",               0.0);// Recoil Energy CCPionInc Spline correction 
  declareDoubleBranch( m_hypMeths, "hadron_recoil_default",       0.0);// Recoil Energy Default Spline correction
  declareDoubleBranch( m_hypMeths, "hadron_recoil_two_track",     0.0);// Recoil Energy with CCPi Spline correction
  //-------------------------------------------------------------------------------------

  //!neutron blob branches
  if (m_doNeutron )
    {
      std::vector<std::string>::iterator it;
      debug()<<"Declaring Blob Branches"<<endmsg;
      declareContainerDoubleBranch( m_hypMeths, "HadronE" );

      debug()<<"BlobIntBranches Size = "<<m_neutronBlobUtils->GetIntDataMemberNames().size()<<endmsg;
      for( uint i = 0; i< m_neutronBlobUtils->GetIntDataMemberNames().size(); i++ ) declareContainerIntBranch( m_hypMeths, m_neutronBlobUtils->GetIntDataMemberNames()[i] );
      for( uint i = 0; i< m_neutronBlobUtils->GetDoubleDataMemberNames().size(); i++ ) declareContainerDoubleBranch( m_hypMeths, m_neutronBlobUtils->GetDoubleDataMemberNames()[i] );
    }
  declareIntBranch( m_hypMeths, "EvtHasNBlobTracks", 0 );
  declareIntBranch( m_hypMeths, "EvtHasNBlobIncTracks", 0 );


  declareIntEventBranch(  "NneutronClusters0", 0);
  declareIntEventBranch(  "NneutronClusters1", 0);
  declareIntEventBranch(  "NneutronClusters2", 0);
  declareIntEventBranch(  "NneutronClusters3", 0);
  declareIntEventBranch(  "PassNeutronMaxClusterCut", 0);
  declareIntEventBranch(  "NeutronMaxClusters", m_maxRecoCluster);
  //for Machine Learning
  declareIntEventBranch( "n_indices" );
  declareContainerIntEventBranch("latticeEnergyIndices", "n_indices");
  declareContainerDoubleEventBranch("latticeNormEnergySums", "n_indices");
  declareContainerDoubleEventBranch("latticeRelativeTimes", "n_indices");

  
  declareContainerIntEventBranch("ANN_segments", 2, -999);
  declareContainerDoubleEventBranch("ANN_plane_probs", 2, -1.0);
  declareContainerIntEventBranch("ANN_vtx_modules");
  declareContainerIntEventBranch("ANN_vtx_planes");

  declareContainerDoubleEventBranch("ANN_vtx");
  declareContainerDoubleEventBranch("muon_corrected_p");

  //organize branches don't care default value.
  std::vector<const IRegisterBranches*> tools_registering_branches;
  tools_registering_branches.push_back(dynamic_cast<IRegisterBranches*>(m_electronUtils));

  std::vector<std::string> branchnames;
  for (auto it_tool = tools_registering_branches.cbegin();it_tool!=tools_registering_branches.cend();++it_tool){
    if (!(*it_tool)) continue;
    branchnames = (*it_tool)->IntEventBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareIntEventBranch(*it_branches,-1);
    }
    branchnames = (*it_tool)->DoubleEventBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareDoubleEventBranch(*it_branches,-999);
    }
    branchnames = (*it_tool)->IntProngBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareIntProngBranch(*it_branches,-1);
    }
    branchnames = (*it_tool)->DoubleProngBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareDoubleProngBranch(*it_branches,-999);
    }
    branchnames = (*it_tool)->IntTruthBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareIntTruthBranch(*it_branches,-1);
    }
    branchnames = (*it_tool)->DoubleTruthBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareDoubleTruthBranch(*it_branches,-999);
    }
    branchnames = (*it_tool)->ContainerIntEventBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareContainerIntEventBranch(*it_branches);
    }
    branchnames = (*it_tool)->ContainerDoubleEventBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareContainerDoubleEventBranch(*it_branches);
    }
    branchnames = (*it_tool)->ContainerIntProngBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareContainerIntProngBranch(*it_branches,4,-1);
    }
    branchnames = (*it_tool)->ContainerDoubleProngBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareContainerDoubleProngBranch(*it_branches,4,-999);
    }
    branchnames = (*it_tool)->ContainerIntTruthBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareContainerIntTruthBranch(*it_branches);
    }
    branchnames = (*it_tool)->DoubleTruthBranches();
    for (auto it_branches = branchnames.cbegin();it_branches!=branchnames.cend();++it_branches ){
      declareContainerDoubleTruthBranch(*it_branches);
    }
  }

  debug()<<"Exiting MasterAnaDev::initialize() ....."<<endmsg;

  return sc;
}


//=============================================================================
// reconstructEvent()
//=============================================================================
StatusCode MasterAnaDev::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth ) const {

  info() << "Entering MasterAnaDev::reconstructEvent() now ....." << endmsg;
  StatusCode sc;

  if (truth) {
    truth->filtertaglist()->setOrAddFilterTag( "reco_hasGoodObjects", false );
    truth->filtertaglist()->setOrAddFilterTag( "reco_isGoodVertex", false );
    truth->filtertaglist()->setOrAddFilterTag( "reco_isFidVol_smeared", false );
  }

  event->filtertaglist()->setOrAddFilterTag( "isMinosMatchTrack", false );

  if( getDAQHeader()->isMCTrigger() ) {
    debug() << truthData() << " Truth Process Type " << truth->processType() << endmsg;
  }
  //========================================================================//
  //for Machine Learning
  Minerva::DAQHeader* header;
  if( exist<Minerva::DAQHeader>( Minerva::DAQHeaderLocation::Default ))
    header = get<Minerva::DAQHeader>( Minerva::DAQHeaderLocation::Default );
  else{
    error() << "Could not locate DAQHeader" << endmsg;
    return StatusCode::FAILURE;
  }
  //
  ulonglong gpstime = header->gpsTime();
  //========================================================================//

  //! Skip events with a truth pointer whose processType is insufficiently defined
  //! The Truth pointer will be NULL while processing data
  //-----------------------------------------------------------------------------------
  if( truth && ( truth->processType()==Minerva::GenMinInteraction::kNoProcess || truth->processType()==Minerva::GenMinInteraction::kUnknownProcess ) ) {
    debug() << "Truth pointer points to an insufficiently defined process - no sense in proceeding !" << endmsg;
    counter("refused unknown truth interaction")++;
    return StatusCode::SUCCESS;
  }
  //=========================================================//

  //for Machine Learning
  Minerva::GenMinHeader *genMinHeader = 0;
  if( truth ){
    if( exist<Minerva::GenMinHeader>( Minerva::GenMinHeaderLocation::Default ) ){
      genMinHeader = get<Minerva::GenMinHeader>( Minerva::GenMinHeaderLocation::Default );
    }
    else {
      error() << "No GenMinHeader in this event!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  //----------------------------------------------------------------------
  //! CUT : Does this event carry a reco object with a BadObject flag?
  //----------------------------------------------------------------------
  if( event->filtertaglist()->isFilterTagTrue( AnaFilterTags::BadObject() ) ) {
    error() << "Found an event flagged with a BadObject! Refusing to analyze ..." << endmsg;
    counter("REFUSED_A_BADOBJECT") += 1;
    return StatusCode::SUCCESS; // Things are bad, but we didn't crash.
  }
  if (truth) truth->filtertaglist()->setOrAddFilterTag( "reco_hasGoodObjects", true);
  //--------------------------------------------------------------------
  //! NOTE: WE DO NOT CUT ON INTERACTION VERTICES ANYMORE. IF YOU WANT
  //! THE ANALYSIS TO BEHAVE AS BEFORE YOU NEED TO CUT ON HAS_INTERACTION_VERTEX
  //--------------------------------------------------------------------


  //count vertices of the start variety used and unused
  VertexVect unused_vtx = event->select<Minerva::Vertex>( "Unused", "StartPoint" );
  VertexVect used_vtx = event->select<Minerva::Vertex>( "Used", "StartPoint" );

  event->setIntData("event_unused_start_vertices",unused_vtx.size());
  event->setIntData("event_used_start_vertices",used_vtx.size());
  //get the number of minos matches in the event
  //get primary prongs. Get output of stubs+matches
  int count_matches = 0;
  std::vector<double> allvtxtimes;
  std::vector<double> allvtxtimes_minos_match;
  std::vector<int> allvtxfv_minos_match;
  for( Minerva::VertexVect::iterator itvtx=used_vtx.begin(); itvtx!=used_vtx.end(); ++itvtx ) {
    debug() << "  Vertex Position          = " << (*itvtx)->position() << endmsg;
    allvtxtimes.push_back((*itvtx)->getTime());
    ProngVect primProngs; ProngVect unattachedProngs;
    m_objectAssociator->getProngs_fromSourceVertex( primProngs, unattachedProngs, *itvtx );

    ProngVect minosTrackProngs, minosStubProngs;
    minosTrackProngs.clear();
    minosStubProngs.clear();

    int isFV = m_minCoordSysTool->inFiducial((*itvtx)->position().X(),(*itvtx)->position().Y(),(*itvtx)->position().Z(),850.0,6127,8183)?1:0;


    for( ProngVect::const_iterator itProng = primProngs.begin(); itProng != primProngs.end(); ++itProng ) {
      if( (*itProng)->MinosTrack() )     minosTrackProngs.push_back( *itProng );
      else if( (*itProng)->MinosStub() ) minosStubProngs.push_back(  *itProng );
    }


    if( !minosTrackProngs.empty() ) {
      allvtxtimes_minos_match.push_back((*itvtx)->getTime());
      allvtxfv_minos_match.push_back(isFV);
      debug() << " Found a vertex with " << minosTrackProngs.size() << " MinosTrack Prong(s)..." << endmsg;
      ++count_matches;
    } else if( !minosStubProngs.empty() ) {
      allvtxtimes_minos_match.push_back((*itvtx)->getTime());
      allvtxfv_minos_match.push_back(isFV);
      debug() << " Found a MinosStub Prong..." << endmsg;
      ++count_matches;

    } // end check for MinosTrack/Stub prongs
  } //end loop over "used" vertices. There never seem to be unused vertices

  event->setContainerDoubleData("all_event_start_vertex_time", allvtxtimes);
  event->setContainerDoubleData("all_event_start_vertex_time_minos_match",allvtxtimes_minos_match);
  event->setContainerIntData("all_event_start_vertex_fv_minos_match",allvtxfv_minos_match);
  event->setIntData("n_minos_matches",count_matches);


  //---------------------------------------------------
  //! Let's Get all slice hits and energies
  //---------------------------------------------------
  std::cout << "Doing slice variables" << std::endl;
  Minerva::IDClusterVect sliceclusters = event->select<Minerva::IDCluster>("All","!LowActivity&!XTalkCandidate");
  std::vector<double> slice_hit_time;
  std::vector<double> slice_hit_energy;
  int slice_n_hits = 0;
  for( Minerva::IDClusterVect::iterator itsc=sliceclusters.begin(); itsc!=sliceclusters.end(); itsc++ ) {
    SmartRefVector<Minerva::IDDigit> digits = (*itsc)->digits();
    for ( SmartRefVector<Minerva::IDDigit>::const_iterator itDig = digits.begin(); itDig != digits.end(); ++itDig ) {
      slice_hit_time.push_back((*itDig)->time());
      slice_hit_energy.push_back((*itDig)->normEnergy());
      slice_n_hits++;
    }
  }

  // If a slice has no hits, skip it!
  if ( slice_n_hits == 0 ) {
    info() << " * Skipping slice with no hits! " << endmsg;
    return StatusCode::SUCCESS;
  }

  // sort out clusters in this slice for use in ElasticFilterTool
  sc = StatusCode::SUCCESS;
  if ( sliceclusters.size() != 0 ) {
    std::vector<double> zRange;
    std::vector<double> xyRange;
    std::vector<double> xySeed;
    // need overloaded version of sortClusters to handle the OD
    sc = m_elasticFilterTool->sortClusters ( sliceclusters, zRange, xyRange, xySeed );
    if ( sc == StatusCode::FAILURE ) {
      warning() << " reconstructEvent() did not provide clusters to sort!" << endmsg;
    }
    else {
      debug() << " slice has " << sliceclusters.size() << " clusters." << endmsg;
      debug() << "  zRange for slice: " << zRange << endmsg;
    }
  }

  event->setContainerDoubleData("slice_hit_time",slice_hit_time);
  event->setContainerDoubleData("slice_hit_energy",slice_hit_energy);
  event->setIntData("slice_n_hits",slice_n_hits);

  std::cout << "This slice has " << slice_n_hits << " hits" << std::endl;
  std::cout << "This slice has " << std::accumulate(slice_hit_energy.begin(),slice_hit_energy.end(),0) << " total hit energy" << std::endl;

  //-----------------------------------------------------------------
  //! Look at muon properties (charge, score, MINOS-match etc.)
  //! CUT: Skip if no muons present in event
  //-----------------------------------------------------------------
  //  SmartRef<Minerva::Particle> primMuon;
  SmartRef<Minerva::Prong> muonProng;
  SmartRef<Minerva::Particle> muonPart;
  bool is_minos_track = false, is_minos_stub = false;
  counter("before back exiting tracks")++;
  bool has_muon = MuonUtils->findMuonProng( event, muonProng, muonPart);
  bool has_electron = false, has_lepton = has_muon, has_NCproton = false;
  event->filtertaglist()->setOrAddFilterTag( "has_muon", has_muon );
  truth && truth->filtertaglist()->setOrAddFilterTag( "reco_has_muon", has_muon );

  if( has_muon ){
    event->setIntData("HasNoBackExitingTracks",0); //Not Electron event
    //-----------------------------------------------------
    //! Use virtual bool function truthIsPlausible( ) to
    //! look for plausible muon
    //-----------------------------------------------------
    muonProng->setIntData( "look_for_plausibility", 1 );

    //--------------------------------------
    //! CUT: Skip events with no muon
    //--------------------------------------
    if( truth ) {
      truth->filtertaglist()->setOrAddFilterTag( "reco_minos_match", has_muon );
      truth->setIntData( "reco_is_minos_match", (int)has_muon );
    }
    counter("has_muon")+=has_muon;
    //--------------------------------------
    //! CUT: Skip muons with no charge
    //--------------------------------------
    int muon_charge = 0;
    MuonUtils->muonCharge( muonProng, muon_charge, m_qOverpChargeCut );
    if( truth ) {
      truth->filtertaglist()->setOrAddFilterTag( "reco_muon_charge", ( muon_charge!=0 ) );
      truth->setIntData( "reco_has_muon_charge", (int)( muon_charge!=0 ) );
    }
    counter("muon_charge") += ( muon_charge!=0 );
    debug() << "Muon Particle Score: " << muonPart->score() << endmsg;
    //--------------------------------------------------
    //! CUT: Only keep muons with a certain score
    //! Find out if muons are MINOS-matched or not
    //--------------------------------------------------
    if (muonPart->score() >= m_minMuonScore) {
      muonProng->filtertaglist()->setOrAddFilterTag( AnaFilterTags::PrimaryMuon(), true );

      if (muonProng->MinosTrack()) is_minos_track = true;
      if (muonProng->MinosStub())  is_minos_stub  = true;
      if (is_minos_stub && is_minos_track) counter("MuonHasMinosStubAndTrack")++;
      if (!is_minos_stub && !is_minos_track) counter("MuonIsNotMinosMatched")++;

      event->filtertaglist()->setOrAddFilterTag("isMinosMatchTrack", is_minos_track );
      event->filtertaglist()->setOrAddFilterTag("isMinosMatchStub", is_minos_stub );
      event->setIntData("isMinosMatchTrack", is_minos_track );
      event->setIntData("isMinosMatchStub", is_minos_stub );
      if (truth) {
        truth->filtertaglist()->setOrAddFilterTag( "reco_isMinosMatchTrack", is_minos_track );
        truth->filtertaglist()->setOrAddFilterTag( "reco_isMinosMatchStub", is_minos_stub );
        truth->setIntData("reco_muon_is_minos_match_track", (int)is_minos_track);
        truth->setIntData("reco_muon_is_minos_match_stub", (int)is_minos_stub);
      }
    } else {
      debug()<<"Muon prong does not pass score cut"<<endmsg;
      counter("muon not passing cut")++;
      return StatusCode::SUCCESS;
    }


  } //  if ( has_muon )

  if ( !has_muon ) {
    counter("consider nue")++;
    has_electron = m_electronUtils->findElectronProng( event, muonProng, muonPart );
    has_lepton = has_muon || has_electron;
  
    if (has_electron) {
      // this next line has a dubious literal translation! assume the goal is to try the steps below for events that lack a muon but contain an electron
      muonProng->filtertaglist()->setOrAddFilterTag( AnaFilterTags::PrimaryMuon(), true );
      counter("nue candidate")++;
      m_electronUtils->PassNuERecoCuts(event);
    }
  }

  // bail here and return success if we are not looking for NC protons
  if (!has_lepton && !m_considerNCprotons) return StatusCode::SUCCESS;

  if ( !has_lepton && m_considerNCprotons ) {
    // THIS IS WHERE WE MUST LOOK FOR PROTON PRONG AS A LONE TRACK
    //has_NCproton = APPROPRIATE FUNCTION CALL FROM ELASTICFILTERTOOL
    debug() << "Event did not have candidate lepton.  Considering Neutral Current proton activity..." << endmsg;
    debug() << " NOT ACTUALLY LOOKING FOR PROTONS YET! " << endmsg;
    // has_NCproton = true; // if we find it
    counter("NCproton_candidate")++;

    // consider neutral current
    bool has_vtx = ( event->hasInteractionVertex() );
    debug() << " No likely muon or electron.  Boolean interaction vertex status: " << has_vtx << endmsg;

    double bareSignalEnergy=-999.9, bareSignalShape=-999.9, averageNormEnergyPerHit=-999.9;
    // AMM: mystery segfault here is fixed by skipping time slices with no hits 
    sc = m_elasticFilterTool->getBareSignalMetrics( bareSignalEnergy, bareSignalShape, averageNormEnergyPerHit );
    debug() << " ElasticFilterTool: Bare signal Energy " << bareSignalEnergy << endmsg;
    debug() << "                  : Bare signal Shape " << bareSignalShape << endmsg;
    debug() << "                  : Average Norm Energy per Hit " << averageNormEnergyPerHit << endmsg;

    ////////return StatusCode::SUCCESS;  // we need to bail on this event until further structure is developed

    // MasterAnaDevRecoUtils has function runDetachedShortTracker( untracked clusters, XYZPoint, short track vector ) to try,
    //  after checking for existence of any existing proton tracks.

    markEvent( event );
    //-------------------------------
    //! Interpret event
    //-------------------------------
    std::vector<Minerva::NeutrinoInt*> nuInt;
    interpretEvent( event, truth, nuInt );

    //---------------------------------------------------
    //! Add interaction hypothesis to physics event
    //---------------------------------------------------
    sc = addInteractionHyp( event, nuInt );
    return sc;
  }


  if( muonProng ){
    Minerva::Track *muTrack = muonProng->minervaTracks()[0];
    double muTrackTime =  m_recoObjTimeTool->trackVertexTime(muTrack);
    event->setDoubleData("muon_trackVertexTime", muTrackTime );
  }

  //m_lowRecoilUtils->SaveRecoilBranches(event,"CCNuE_Nu_Tracker",muonProng);

  //--------------------------------------------------
  //! CUT: Check if there is an interaction vertex
  //--------------------------------------------------
  bool has_vtx = ( event->hasInteractionVertex() );
  truth && truth->filtertaglist()->setOrAddFilterTag( "reco_has_int_vtx", has_vtx );
  counter("has_vertex") += has_vtx;
  //if( !event->hasInteractionVertex() ){
  if( !has_vtx){
    debug() << "Event doesn't have an interaction vertex" << endmsg;
    std::cout << "NO INTERACTION VERTEX" << std::endl;
    event->setIntData("has_interaction_vertex",0);

    //Check if at least one start vertex is in FV. Don't save rock muons.
    int fv_count = 0;
    for(size_t fv_e=0;fv_e<allvtxfv_minos_match.size();fv_e++){
      fv_count+=allvtxfv_minos_match[fv_e];
    }

    //if(fv_count==0)
    counter("refused bad vertex")++;
    return StatusCode::SUCCESS; //none of the minos matches orginate from vertices in the FV or there were no minos matches. Both cases we don't care about.

    markEvent( event );
    //-------------------------------
    //! Interpret event
    //-------------------------------
    std::vector<Minerva::NeutrinoInt*> nuInt;
    interpretEvent( event, truth, nuInt );

    //---------------------------------------------------
    //! Add interaction hypothesis to physics event
    //---------------------------------------------------
    sc = addInteractionHyp( event, nuInt );
    return sc;

  }
  else{
    event->setIntData("has_interaction_vertex",1);
  }
  if (truth) truth->filtertaglist()->setOrAddFilterTag( "reco_isGoodVertex", true);
  //---------------------------------------------------------------
  //! Get the interaction vertex
  //! If you can't get the vertex, tag the event as a BadObject
  //---------------------------------------------------------------
  SmartRef< Minerva::Vertex > vertexOfPhysEvent = (event->interactionVertex()).target();
  if( !vertexOfPhysEvent ) {
    bool pass = true; std::string tag = "BadObject";
    event->filtertaglist()->addFilterTag(tag,pass);
    event->setIntData("NullVertex",1);
    error() << "This vertex is NULL! Flag this event as bad!" << endmsg;
    return StatusCode::SUCCESS;
  }

  //-----------------------------------------------------------------------------------

  //for Machine Learnig
  if (m_getLatticeEnergies) {
  std::vector<int> latticeEnergyIndices;
  std::vector<double> latticeNormEnergySums;
  std::vector<double> latticeRelativeTimes;

    double signature_time = event->getDoubleData("muon_trackVertexTime");
    debug() << "filling the branches for the lattice energies" << endmsg;
    m_MLVFTool->getLatticeValues(event,
                                 signature_time,
                                 latticeEnergyIndices,
                                 latticeRelativeTimes,
                                 latticeNormEnergySums);

    
  event->setIntData("n_indices", latticeEnergyIndices.size() );
  event->setContainerIntData("latticeEnergyIndices", latticeEnergyIndices);
  event->setContainerDoubleData("latticeNormEnergySums", latticeNormEnergySums);
  event->setContainerDoubleData("latticeRelativeTimes", latticeRelativeTimes);
  
   }

///////////////////////Veto code shifted after image file code

  Minerva::VetoDigits* vdigits = 0;
  try {
    vdigits = get<Minerva::VetoDigits>(Minerva::VetoDigitLocation::Default);
  } catch(...) {
    error() << "Invalid VetoDigitLocation " << endmsg;
    return StatusCode::SUCCESS;
  }

  //This fills a map which associates discriminatorPairs to the time when they had their 6th push
  std::map<Minerva::DiscrPairID,double> timeOfSixthPush;
  std::vector<std::pair <int,int> > WallPad;
  std::vector<int> wall1Paddles, wall2Paddles;
  double earliestSixPushTime = 9999999.9;
  bool sixPushes = false;
  int sixPush_matchtrack = 0;
  for ( Minerva::VetoDigits::iterator v = vdigits->begin(); v != vdigits->end(); ++v) {

    Minerva::ChannelID c( (*v)->key() );
    Minerva::DiscrPairID discrID = m_plexModel->getDiscrPair( c );
    if ( (*v)->discrFired() == 1 && c.hit() == 5 ) {
      sixPushes = true;
      if( ((*v)->rawTime()) < earliestSixPushTime) {
        earliestSixPushTime = (*v)->rawTime();
        WallPad = discrPairToPaddle(discrID.minervaID());
      }
      if ( timeOfSixthPush.find(discrID) == timeOfSixthPush.end() ) {
        timeOfSixthPush[discrID] = (*v)->rawTime();
      }
      else if ( (*v)->rawTime() < timeOfSixthPush[discrID] ) {
        timeOfSixthPush[discrID] = (*v)->rawTime();
      }
    }//if discrFired
  }//loop over vetodigits


  //Will tell us what the earliest time is that 6 pushes from the Veto Wall occured and on what paddle
  if (sixPushes) {
    debug() << "Found a sixth push at time: " << earliestSixPushTime << endmsg;

    //Fill Paddle vectors with the dead paddles from both walls
    for( std::vector< std::pair<int,int> >::iterator it = WallPad.begin(); it != WallPad.end(); it++ ){
      info() << "veto (wall,paddle) in saturated discriminator = (" << it->first << "," << it->second << ")" << endmsg;
      if(it->first==1){wall1Paddles.push_back(it->second) ;}
      if(it->first==2){wall2Paddles.push_back(it->second) ;}
    }

    event->setContainerIntData("VetoWall_sixPush_paddle_wall1", wall1Paddles);
    event->setContainerIntData("VetoWall_sixPush_paddle_wall2", wall2Paddles);
  }


  //////////////////////////////////////////////////////////
  //End sixpush
  //////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////
  // make map if pmt are on or off
  ////////////////////////////////////////////////////////////////////
  typedef std::map<int, bool>  Paddle_StatusOFFMap;
  typedef std::map< std::pair<int,int>, bool>  Paddle_StatusOFFMap_1;


  Paddle_StatusOFFMap paddle_statusMap;
  Paddle_StatusOFFMap_1 paddle_statusMap_1;
  std::vector<const Minerva::DeVetoWall*>::const_iterator itWall;
  const std::vector<const Minerva::DeVetoWall*>  walls = m_VetoDet->getDeVetoWalls();
  unsigned int wall_number;
  debug()<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endmsg;
  debug()<<"Filling Veto status Map = " << endmsg;
  debug()<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endmsg;
  std::vector<int> MapVector;
  int counterVETOMAP[m_totalVETOPMTS];
  int input1=999;
  int input2=999;
  int count=0;
  for(itWall = walls.begin(); itWall !=walls.end(); ++itWall){
    wall_number = ((*itWall)->getWallID()).wall();
    for (int j = 1; j<=m_NumVetoPaddles;j++){
      count++;
      //Gets the ID of each pmt on the VETOwall
      Minerva::VetoPMTID pmtid_1( Minerva::DetectorID::Veto, wall_number, j, m_VetoPaddlePMT1 );
      Minerva::VetoPMTID pmtid_2( Minerva::DetectorID::Veto, wall_number, j, m_VetoPaddlePMT2 );
      debug()<<" Geting PMT ID -  For Wall_number: "<< wall_number << "  Paddle:" << j <<endmsg;
      int veto_1=-999;
      int veto_2=-999;
      //for a given PMT ID - check the pmt status returns true if status is good

      //veto_1 = m_getVetoEff->getPMTStatus_hw( pmtid_1,  gpstime       );
      //veto_2 = m_getVetoEff->getPMTStatus_hw( pmtid_2,  gpstime       );

      if(veto_1 == -999 ){warning()<<"Failure:: getPMTStatus_hw() = " << veto_1 << " for PMTID1 = "<<pmtid_1<<endmsg; }
      if(veto_2 == -999 ){warning()<<"Failure:: getPMTStatus_hw() = " << veto_2 << " for PMTID2 = "<<pmtid_2<<endmsg; }

      debug()<< "Geting PMT Status for pmtid_1( "<< pmtid_1<<" ) = " << veto_1 << endmsg;
      debug()<< "Geting PMT Status for pmtid_2( "<< pmtid_2<<" ) = " << veto_2 << endmsg;

      if(veto_1==1){input1=1;}
      else{input1= -1;}
      if(veto_2==1){input2=1;}
      else{input2= -1;}
      //Veto_status_Map[pmtid_1]=std::make_pair (gpstime,veto_1);
      //Veto_status_Map[pmtid_2]=std::make_pair (gpstime,veto_2);
      debug()<< "counter wall_number*j + 1 = "<< count - 1<< endmsg;
      counterVETOMAP[count - 1]= input1;
      debug()<< "counter wall_number*j + 11 = "<< count + 11<< endmsg;
      counterVETOMAP[count + 11 ]= input2;
      if(veto_1==false || veto_2==false){
        paddle_statusMap[j]=false;
        paddle_statusMap_1[std::make_pair(wall_number,j)]=false;

      }
      else{paddle_statusMap[j]=true;
        paddle_statusMap_1[std::make_pair(wall_number,j)]=true;
      }

      debug()<<"wall_number = "<<wall_number << "j ="<< j <<"  paddle_statusMap_1[make_pair(wall_number,j)] = "<<  paddle_statusMap_1[std::make_pair(wall_number,j)]<<endmsg;

    }



  }

  debug()<<"filling vector MAP"<< endmsg;
  for(int i =0; i<m_totalVETOPMTS;i++ ){
    MapVector.push_back(counterVETOMAP[i]);
    debug()<<"intputs of Map index= "<< i << " Input =  " <<counterVETOMAP[i]<<endmsg;
  }


  event->setContainerIntData("VetoWall_PMTStatusMAP", MapVector);

  //---------------------------------------------------------------------------
  //! CUT: Check if vertex is inside fiducial volume
  //! according to SignalApothem, SignalUpstreamZ, and SignalDownstreamZ.
  //---------------------------------------------------------------------------

  SmartRef< Minerva::Vertex > interactionVertex = event->interactionVertex();
  bool is_fiducial = eventVertexIsFiducial( event );
  Gaudi::XYZPoint vtx_position = interactionVertex->position();
    if (vtx_position.z() != vtx_position.z()) {
      warning()<<"NaN vertex!"<<endmsg;
      return StatusCode::SUCCESS;
    }


  if( truth ) {
    truth->filtertaglist()->setOrAddFilterTag( "reco_fiducial", is_fiducial );
    truth->setIntData( "reco_is_fiducial", (int)is_fiducial );

    Gaudi::XYZPoint vtx_smear(vtx_position.x() + m_randomGen->Gaus(0.0,0.91*CLHEP::mm), vtx_position.y() + m_randomGen->Gaus(0.0,1.25*CLHEP::mm), vtx_position.z() + m_randomGen->Gaus(0.0, 10.0*CLHEP::mm));
    if ( FiducialPointTool->isFiducial(vtx_smear, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ ) ) {
      truth->filtertaglist()->setOrAddFilterTag( "reco_isFidVol_smeared", true );
    }
    else truth->filtertaglist()->setOrAddFilterTag( "reco_isFidVol_smeared", false );
  }
  counter("is_fiducial") += is_fiducial;
  if( !FiducialPointTool->isFiducial(vtx_position, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ ) ) {
    debug() << "Event vertex is outside fiducial volume" << endmsg;
 //   return StatusCode::SUCCESS;
  }



  //==================================================================//
  //-------------------------------------------------------------------------
  // for Machine Learning:  RUN MACHINE LEARNING PREDICTION
  //-------------------------------------------------------------------------
  vector<int> ANN_segments;
  vector<double> ANN_plane_probs;
  int DNN_segment;
  double DNN_plane_probs;

  if ( has_muon && muonProng && m_useDNN ){
    std::string Cuts;
    if( truth ) Cuts = Form("run==%d && subrun==%d && gate==%d && slice_number==%d", genMinHeader->RunNumber(), genMinHeader->SubRunNumber(),truth->NumEventInFile(), event->sliceNumbers()[0]);
    else Cuts = Form("run==%d && subrun==%d && gate==%d && slice_number==%d", header->runNumber(), header->subRunNumber(),header->gateNumber(), event->physicsEventNumber());

    // Draw with option goff and generate three variables
    debug() << "DNN Cuts = " << Cuts << endmsg;
    //first segment predictions
    Int_t n = dbPred->Draw("ANN_best_segment:ANN_best_prediction",Cuts.c_str(),"goff");
    debug() << "DNN The arrays' dimension is " << n << endmsg;
    if( n == 1 ){
      // Retrieve variables 1 and 2
      DNN_segment  = (int)*dbPred->GetVal(0);
      DNN_plane_probs = *dbPred->GetVal(1);
      info() << "Found prediction: DNN segment, DNN plane probs = " << DNN_segment << " " << DNN_plane_probs << endmsg;
      ANN_segments.push_back(DNN_segment);
      ANN_plane_probs.push_back(DNN_plane_probs);
    }
    else{
      if( truth ) error() << "No ML Prediction for this event: " << genMinHeader->RunNumber() << ", " << genMinHeader->SubRunNumber() << ", " << truth->NumEventInFile() << ", " << event->sliceNumbers()[0] <<  endmsg;
      else error() << "No ML Prediction for this event: " <<  header->runNumber() << ", " << header->subRunNumber() << ", " << header->gateNumber() << ", " << event->physicsEventNumber() <<  endmsg;
    //trun this off   
    return StatusCode::SUCCESS;
    }
    //second segment predictions
    Int_t m = dbPred->Draw("ANN_segment_2:ANN_prob_2",Cuts.c_str(),"goff");
    if( m == 1 ){
      // Retrieve variables 1 and 2
      DNN_segment  = (int)*dbPred->GetVal(0);
      DNN_plane_probs = *dbPred->GetVal(1);
      info() << "Found prediction: DNN segment, DNN plane probs = " << DNN_segment << " " << DNN_plane_probs << endmsg;
      ANN_segments.push_back(DNN_segment);
      ANN_plane_probs.push_back(DNN_plane_probs);
    }
    debug() << "DNN The arrays' dimension is " << n << endmsg;
    event->setContainerIntData("ANN_segments", ANN_segments );
    event->setContainerDoubleData("ANN_plane_probs", ANN_plane_probs );
  }

  //--------------------------------------------------------------------------------------------------------------------
  // Find chi2 of MINERvA part (not MINOS) of muon track (for comparisons with J. Wolcott's chiPerDoF for electrons)
  //--------------------------------------------------------------------------------------------------------------------
  if(muonProng) {
    SmartRef<Minerva::Track> muonTrack;
    bool foundMinervaPart = findLongestMinervaTrack( muonProng, muonTrack );
    if( foundMinervaPart && muonTrack->nNodes() ) {
      verbose() <<"Chi2 of the MINERvA part of the muon track is " << muonTrack->chi2PerDoF() <<" Nodes: "<< muonTrack->nNodes() << endmsg;
      event->setDoubleData( "muon_minerva_trk_chi2PerDoF", muonTrack->chi2PerDoF() );
    }

    //Get muontime
    double muonTime = m_recoObjTimeTool->trackVertexTime(muonProng->minervaTracks()[0]);
    //---------------------------------------------------
    //! Store all node angles to fix angular biases due to hadronic energy near the vertex (docdb 12047)
    //---------------------------------------------------
    std::vector<double> thetaNodes, thetaXNodes, thetaYNodes;
    for(uint i=0; i<muonTrack->nNodes(); i++){
      const Minerva::Node* node=muonTrack->nodes()[i];
      double ax=node->state().ax();
      double ay=node->state().ay();
      double theta=m_minCoordSysTool->thetaWRTBeam( ax, ay, 1 );
      double thetaX=m_minCoordSysTool->thetaXWRTBeam( ax, ay, 1 );
      double thetaY=m_minCoordSysTool->thetaYWRTBeam( ay, 1 );
      thetaNodes.push_back(theta);
      thetaXNodes.push_back(thetaX);
      thetaYNodes.push_back(thetaY);
    }

    event->setContainerDoubleData("muon_theta_allNodes", thetaNodes);
    event->setContainerDoubleData("muon_thetaX_allNodes", thetaXNodes);
    event->setContainerDoubleData("muon_thetaY_allNodes", thetaYNodes);

    //----------------------------------------------------------------------------------------
    //! Add PID to prongs not associated with the interaction vertex, but are in time
    //----------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------
    //! Veto wall: muon track  extrapolation to Wall checks
    //----------------------------------------------------------------------------------------
    // Veto wall Varibles
    //----------------------------------------------------------------------------------------

    debug()<<"Starting Veto wall Checks" << endmsg;
    bool found_VetoMatchedTrk = false;

    // Variables for extracting information from the prong/event
    int     wallHit=0;
    double  xPos[2];         //Xpos of track on veto wall
    double  yPos[2];         //YPos of track on veto wall

    // -------------------------------------------------------------------
    // Now we're going to check if the muon Prong was matched in time to the veto wall
    if((muonProng)->Veto()) {

      found_VetoMatchedTrk = true;
      debug() << "FOUND A MUON PRONG THAT WAS MATCHED IN TIME TO THE VETO WALL" << endmsg;

    }//if muonProng is veto matched in time

    // -------------------------------------------------------------------
    // Lets get the veto wall extrapolated location from the prong, regardless of a match in time
    // Fill arrays with information stored in the prong, about the veto track extrapolation
    // WallHit: 0 to 3.  0 = no walls, 1 = only Wall 1, 2 = only Wall 2, 3 = both Walls

    if ( (muonProng) -> hasIntData("WallHit"))    wallHit    = (muonProng) -> getIntData("WallHit");
    if ( (muonProng) -> hasDoubleData("WallOneXPosition")) xPos[0] = (muonProng) -> getDoubleData("WallOneXPosition");
    if ( (muonProng) -> hasDoubleData("WallOneYPosition")) yPos[0] = (muonProng) -> getDoubleData("WallOneYPosition");
    if ( (muonProng) -> hasDoubleData("WallTwoXPosition")) xPos[1] = (muonProng) -> getDoubleData("WallTwoXPosition");
    if ( (muonProng) -> hasDoubleData("WallTwoYPosition")) yPos[1] = (muonProng) -> getDoubleData("WallTwoYPosition");

    // if this track has wallHit = 0 it didn't land on any walls, so skip (in which case, how do we get the filter set?)
    if (wallHit == 0) event->filtertaglist()->setOrAddFilterTag("VetoWall_extrapnowalls", true);

    if (wallHit == 1) {
      event->filtertaglist()->setOrAddFilterTag("VetoWall_extrapwall1", true);
      event->setDoubleData("trk_XPosW1",xPos[0]);
      event->setDoubleData("trk_YPosW1",yPos[0]);
    } // if only wall 1 is hit

    else if(wallHit == 2){
      event->filtertaglist()->setOrAddFilterTag("VetoWall_extrapwall2", true);
      event->setDoubleData("trk_XPosW2",xPos[1]);
      event->setDoubleData("trk_YPosW2",yPos[1]);
    } // if only wall 2 is hit

    else if (wallHit == 3) {
      event->filtertaglist()->setOrAddFilterTag("VetoWall_extrapbothwalls", true);
      event->setDoubleData("trk_XPosW1",xPos[0]);
      event->setDoubleData("trk_YPosW1",yPos[0]);
      event->setDoubleData("trk_XPosW2",xPos[1]);
      event->setDoubleData("trk_YPosW2",yPos[1]);
    } // if both walls are hit


    // -------------------------------------------------------------------
    // Now we're going to check if the muon track time was after a sixth Push from a discriminator on the veto wall
    if (sixPushes) {

      double earlytracktime = 99999999.9;

      SmartRefVector<Minerva::Track> minervaTracks = (muonProng)->minervaTracks();
      for(SmartRefVector<Minerva::Track>::iterator itTrack=minervaTracks.begin(); itTrack!=minervaTracks.end(); itTrack++) {
        double time = m_timeTool->trackVertexTime(*itTrack);
        debug() << "Track time is " << time << endmsg;
        if (time < earlytracktime) { earlytracktime = time;
          debug() << "Updating earliest track time to " << earlytracktime << endmsg; }
      }//end loop over tracks in muonProng

      debug() << "Correcting earlytracktime for time of flight" << endmsg;

      //first veto wall is ~3m away from front face of MINERvA, gives a 1ns time difference, any more detail is smaller than our timing resolution
      earlytracktime -= (3/CLHEP::c_light);

      debug() << "The earliest track time (TOF corrected) = " << earlytracktime << " and the earliest 6 push time = " << earliestSixPushTime << endmsg;

      if (earlytracktime > earliestSixPushTime) {
        debug() << "This event likely happened in the veto wall dead time, let's check the paddle extrap" << endmsg;
        bool sixPushes_1 = false;
        bool sixPushes_2 = false;

        if(muonProng->hasIntData("WallOnePaddleExtrap")){

          if ((wall1Paddles[0] == (muonProng->getIntData("WallOnePaddleExtrap"))) || (wall1Paddles[1] == (muonProng->getIntData("WallOnePaddleExtrap"))) ){
            sixPushes_1 = true;
          }//muon paddle wall 1 extrap = dead paddle on wall 1
        }//muon extrapd to a paddle on wall 1

        if(muonProng->hasIntData("WallTwoPaddleExtrap")){
          if ((wall2Paddles[0] == (muonProng->getIntData("WallTwoPaddleExtrap"))) || (wall2Paddles[1] == (muonProng->getIntData("WallTwoPaddleExtrap"))) ){
            sixPushes_2 = true;
          }//muon paddle wall 2 extrap = dead paddle on wall 2
        }//muon extrapd to a paddle on wall 2

        if(sixPushes_1 == true && sixPushes_2 == true){sixPush_matchtrack=1;}
        else{ sixPush_matchtrack=0;}

        event->setIntData("VetoWall_sixPush", sixPush_matchtrack);

      }//if a sixth push occured before the muon tracktime

    }//if we have sixPushes from the veto wall this gate
    ////////////////////////////////////////////////////////////////
    debug()<<"Enter Vetowall Check if powered off"<< endmsg;
    ////////////////////////////////////////////////////////////////
    const int NWalls = 2;
    double wall_z = 0.;
    bool   retry = false;
    bool   overlap = false;
    int paddleDiff = 0;
    bool extrapVeto = false;
    double YPosVeto[NWalls + 1];
    double XPosVeto[NWalls + 1];
    std::vector<int> overlap_paddlevector;
    std::vector<double> XandYextrapolate_VETO;
    Paddle_StatusOFFMap_1 paddle_statusMap_extra;
    std::vector < std::pair<int,int> > extrapolate_Veto;

    // ******************************************
    // Zero out quantities here. This way something sensible gets written to the TES even if we don't have a candidate track.

    for(int i = 0; i< NWalls + 1; ++i) {
      YPosVeto[i]          = -9999.9;
      XPosVeto[i]          = -9999.9;
    }


    int countwallveto = 0;
    for(itWall = walls.begin(); itWall !=walls.end(); ++itWall){

      Gaudi::XYZPoint positive_point_in_wall;
      Gaudi::XYZPoint negative_point_in_wall;

      Minerva::PaddleID PaddleEx;
      Minerva::PaddleID SecondPaddleEx;

      wall_number = ((*itWall)->getWallID()).wall();
      wall_z = (*itWall)->getZCenter();

      debug() << endmsg;
      debug() << "Veto Wall " << wall_number << " Z: " << wall_z << endmsg;

      m_propagateToVeto->propagate(muonTrack, wall_z + m_PanelOffset, positive_point_in_wall);
      m_propagateToVeto->propagate(muonTrack, wall_z - m_PanelOffset, negative_point_in_wall);

      // ******************************************
      // First query the geometry to see if the point is any paddle, try positive offset paddle first.
      try{
        debug() <<"Trying Paddle Ex" <<endmsg;
        PaddleEx = (*itWall)->getPaddleID(positive_point_in_wall);
        overlap = true;
      }
      catch(MinervaException e) {
        debug() << "Point outside Veto Wall, Trying again" << endmsg;
        retry = true;
      } // end of positive offset test.


      // ******************************************
      // Try to find the paddle on the negative offset paddle.
      if( retry ){
        try{
          PaddleEx = (*itWall)->getPaddleID(negative_point_in_wall);
        }
        catch(MinervaException e) {


          debug() << "Point outside Veto Wall, Moving on to the next wall" << endmsg;
          debug() << " .. pos point x, y = " <<  positive_point_in_wall.x() << "," << positive_point_in_wall.y()<< endmsg;
          debug() << " .. neg point x, y = " <<  negative_point_in_wall.x() << "," << negative_point_in_wall.y()<< endmsg;
          continue;

        }

      }//END OF RETRY

      if(overlap){
        try{
          SecondPaddleEx = (*itWall)->getPaddleID(negative_point_in_wall);
          debug() << "Overlap found!" << endmsg;
        }
        catch(MinervaException e) {
          overlap = false;
        }
      } // end of if overlap

      ////////////////////////////////////////
      // check  overlap
      ////////////////////////////////////////
      if( overlap ) {
        paddleDiff = PaddleEx.paddle() - SecondPaddleEx.paddle();

        // if in overlap and have wrong difference between paddle numbers, exit this prong!
        if( -1 != paddleDiff && 1 != paddleDiff ) {
          warning() << "Difference between overlaping paddles is not 1!!!" << endmsg;
          warning() << "Diference between overlaping paddles is " <<  paddleDiff << endmsg;
          continue;
        }

        // after we've exited the prong if bad delta paddle - so now we have an overlap situation, where delta paddle is valid

      }

      overlap_paddlevector.push_back (PaddleEx.paddle());
      extrapolate_Veto.push_back(std::make_pair(wall_number,PaddleEx.paddle()));

      debug() << "Track extrapolate to Paddle: " << PaddleEx.paddle() << endmsg;

      if(overlap){ debug() << "Track also extrapolate to Second Paddle: " << SecondPaddleEx.paddle() << endmsg;
        overlap_paddlevector.push_back (SecondPaddleEx.paddle());
        extrapolate_Veto.push_back(std::make_pair(wall_number,SecondPaddleEx.paddle()));

      }
      // ******************************************
      // now have either a single extrap or overlap extrap track
      // if overlap is true, found a point in pos and neg z of the same wall
      // if overlap is false and retry is false, found a point only in pos z
      // if overlap is false and retry is true, found a point only in neg z

      // get the x, y position of the extrapolated track
      if ( (overlap) || (!overlap && !retry) ) {// if overlap is false and retry is false, found a point only in pos z
        XPosVeto[wall_number] = positive_point_in_wall.x();
        YPosVeto[wall_number] = positive_point_in_wall.y();

        XandYextrapolate_VETO.push_back(positive_point_in_wall.x());
        XandYextrapolate_VETO.push_back(positive_point_in_wall.y());

      }
      else if (!overlap && retry) {    // if overlap is false and retry is true, found a point only in neg z
        XPosVeto[wall_number] = negative_point_in_wall.x();
        YPosVeto[wall_number] = negative_point_in_wall.y();

        XandYextrapolate_VETO.push_back(negative_point_in_wall.x());
        XandYextrapolate_VETO.push_back(negative_point_in_wall.y());

      }


      debug() << "Track at Veto Wall X:" << XPosVeto[wall_number] << endmsg;
      debug() << "Track at Veto Wall Y:" << YPosVeto[wall_number] << endmsg;

      countwallveto++;
      extrapVeto = true;

    }//end of itwall loop





    event->setIntData("VetoWall_NumberMatchToVeto", countwallveto);

    debug()<<"Passed extrapolating to VetoWall"<<endmsg;



    if(extrapVeto==true){
      if (countwallveto==2){
        event->setContainerDoubleData("VetoWall_muon_extrapVetoXY", XandYextrapolate_VETO);
        debug()<<"Found 2  extrapolated!!"<<endmsg;
      }

      if (countwallveto==1){
        XandYextrapolate_VETO.push_back(-999);
        XandYextrapolate_VETO.push_back(-999);
        event->setContainerDoubleData("VetoWall_muon_extrapVetoXY", XandYextrapolate_VETO);
        debug()<<"Found 1  extrapolated!!"<<endmsg;
      }



      bool pmtOfforOn=true;

      for(std::vector< std::pair<int,int> >::const_iterator it = extrapolate_Veto.begin(); it != extrapolate_Veto.end();++it) {

        //Using the list of Paddle ©that were extrapolated too and see what corrpening Paddle_map of veto status.
        pmtOfforOn  =  paddle_statusMap_1[std::make_pair(it->first,it->second)];
        debug()<<"it->first " <<it->first <<"it->second = " << it->first<<"   paddle_statusMap_1[std::make_pair(it->first,it->first)]= "<<paddle_statusMap_1[std::make_pair(it->first,it->second)] <<endmsg;

        if(pmtOfforOn==false){continue;}

      }

      if(pmtOfforOn==false){
        debug()<<"Matched to paddle that has a PMT marked OFF"<<endmsg;
        event->filtertaglist()->setOrAddFilterTag("VetoWall_MuonTrkMatchToVETOwalloff", true);
      }
      else{
        // didn't match to map a paddle that was off
        event->filtertaglist()->setOrAddFilterTag("VetoWall_MuonTrkMatchToVETOwalloff", false);
      }
    }//end of extrapVeto==true

    else{
      debug()<<"Found 0  extrapolated!!"<<endmsg;
      XandYextrapolate_VETO.push_back(-999);
      XandYextrapolate_VETO.push_back(-999);
      XandYextrapolate_VETO.push_back(-999);
      XandYextrapolate_VETO.push_back(-999);
      event->setContainerDoubleData("VetoWall_muon_extrapVetoXY", XandYextrapolate_VETO);
      event->filtertaglist()->setOrAddFilterTag("VetoWall_MuonTrkMatchToVETOwalloff", false);
    }
    event->setIntData("VetoWall_event_IsVeto", (int)found_VetoMatchedTrk);

    //----------------------------------------------------------------------------------------
    //! END of Vetowall Checks
    //----------------------------------------------------------------------------------------

    Minerva::ProngVect extraProngs = m_masterAnaDevRecoUtils->getInTimeTrackPID(event,muonTime);
    std::vector<double> extraPID;
    std::vector<Minerva::Particle::ID> particleHypotheses;
    particleHypotheses.push_back(Minerva::Particle::Pion);
    particleHypotheses.push_back(Minerva::Particle::Proton);

    for( Minerva::ProngVect::iterator itProng = extraProngs.begin(); itProng!= extraProngs.end();++itProng){
      debug()<<"Make particles from extra prongs"<<endmsg;
      IParticleMakerTool::NameAliasListType toolsToUse;
      toolsToUse.push_back( std::make_pair("dEdXTool", "dEdXTool") );  // name, alias pair

      StatusCode sc_dEdX = m_particleMaker->makeParticles( *itProng, particleHypotheses, toolsToUse );
      if (!sc_dEdX) {
        debug()<<"Could not make particles in extra prong!"<<endmsg;
      }
    }//End loop over prongs to create particles
    Minerva::ProngVect extrasecondaryProtonProngs;
    SmartRef<Minerva::Prong> extraprotonProng;
    SmartRef<Minerva::Particle> extraprotonPart;
    bool hasextraProton = getProtonProngs( extraProngs, extrasecondaryProtonProngs, extraprotonProng, extraprotonPart );
    //now get PID scores for each prong (only save the equivalent of score1 which is used for protons associated with the primary vtx
    //TODO: add systematic variations into the flow....
    if(hasextraProton){
      extraPID.push_back(extraprotonPart->getDoubleData(ParticleExtraDataDefs::dEdXScore1()));
      for(Minerva::ProngVect::iterator itProng = extrasecondaryProtonProngs.begin(); itProng!= extrasecondaryProtonProngs.end();++itProng){
        SmartRef<Minerva::Particle> extraprotonpart = (*itProng)->bestParticle();
        extraPID.push_back(extraprotonpart->getDoubleData(ParticleExtraDataDefs::dEdXScore1()));
      }
    }//end if has protons
    event->setContainerDoubleData("event_extra_track_PID",extraPID);

    //-------------------------------------------------------------------------------------
    //! Create vertex-anchored short track Prongs, refit vertex
    //! Does it make sense to refit the Vtx reliably with such short tracks ? Not sure
    //-------------------------------------------------------------------------------------
    int n_anchored_long_trk_prongs = event->primaryProngs().size() - 1;
    event->setIntData("n_anchored_long_trk_prongs", n_anchored_long_trk_prongs );

    //-------------------------------------------------------------------------------
    //! Do you want to obtain the vector of ID Clusters from the PhysEvent first ?
    //-------------------------------------------------------------------------------
    Minerva::IDClusterVect idClusterVector = m_masterAnaDevRecoUtils->getInTimeAnalyzableIDClusters( event, muonProng );
    if( !idClusterVector.size() ) {
      debug() << "Size of vector containing ID Clusters is " << idClusterVector.size() << " - nothing to make short Tracks from !" << endmsg;
    }
    event->setIntData("ID_Hits_Post_Long", idClusterVector.size());

    //-------------------------------------------------------------------------------------------
    //! A vector for holding all the newly made short tracks (from the clusters near the Vtx)
    //-------------------------------------------------------------------------------------------
    std::vector<Minerva::Track*> *deShortTrackVector = new std::vector<Minerva::Track*>;

    //----------------------------------------------------------------------------------------------------
    //! Make short tracks either with the VertexAnchoredShortTracker or VertexEnergyStudyTool or both
    //----------------------------------------------------------------------------------------------------
    if( m_makeShortTracks && idClusterVector.size() != 0) {
      sc = m_masterAnaDevRecoUtils->runShortTrackers( event, deShortTrackVector, m_runVertexAnchoredShortTracker, m_runVertexEnergyStudyTool );
      if( sc.isFailure() ) return sc;
    }

    //---------------------------------------------------------------------
    //! Delete the vector created, don't need it anymore
    //! Do you need to delete the track pointers (commented out line) ?
    //---------------------------------------------------------------------
    // for (std::vector<Minerva::Track*>::iterator itTrk =deShortTrackVector->begin(); itTrk != deShortTrackVector->end(); ++itTrk) { delete *itTrk; }
    deShortTrackVector->clear();
    delete deShortTrackVector;

    int n_anchored_short_trk_prongs = event->primaryProngs().size() - n_anchored_long_trk_prongs - 1;
    event->setIntData("n_anchored_short_trk_prongs", n_anchored_short_trk_prongs);

    //---------------------------------------------------------------------
    //! CUT: check if the number of outgoing tracks == 2 (multiplicity)
    //---------------------------------------------------------------------
    int multiplicity = vertexOfPhysEvent->getNOutgoingTracks();
    event->setIntData( "multiplicity",  multiplicity );

    debug()<<"Multiplicity = "<<multiplicity<<endmsg;
    SmartRefVector<Minerva::Track> Tracks = vertexOfPhysEvent->getOutgoingTracks();
    debug()<<"Outgoing tracks size =" << Tracks.size() << endmsg;
    for(SmartRefVector<Minerva::Track>::iterator itTrack = Tracks.begin(); itTrack != Tracks.end() ; ++itTrack) {
      verbose() << "Track address = "<< *itTrack << " Track Pat Rec History = " << (*itTrack)->patRecHistory()
                << " Track type = " << (*itTrack)->type() << " Track nodes = " << (*itTrack)->nNodes()
                << " Track VisE = " << (*itTrack)->visibleEnergy() << " Track theta = " << (*itTrack)->theta()
                << " Track phi = " << (*itTrack)->phi() << " Track first cluster time = " << (*itTrack)->firstNode()->idcluster()->time()
                << " Track time = " << (*itTrack)->time() << endmsg;
    }

    int fitConverged = 0;
    event->interactionVertex()->getIntData( VertexExtraDataDefs::FitConverged(), fitConverged );


    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Adjust the vertex position and muon momentum to the center of a passive target where applicable
    //if( has_vtx && has_lepton && muonProng && muonPart && !m_useDNN )
    if( has_vtx && has_lepton && muonProng && muonPart )
      {
        int targetCode = 0;
        getTarget( event->interactionVertex()->position(), targetCode );
        int targetID = (targetCode - targetCode%1000) / 1000;
        if( 1 <= targetCode && targetID <= 5 )
          {
            // If the vertex fit did not converge (or vertex has only 1 track),
            // then change the vertex position at the position where the muon intersects the zCenter of the relevant passive target.
            if( !fitConverged )
              {
                debug() << " This is a passive target event without a good vertex fit.  I will adjust the vertex position to the center of the target." << endmsg;
                SmartRef<Minerva::Vertex> primVtx = event->interactionVertex();
                adjustVertexIntoPassiveTarget( primVtx, muonPart );
              }


            // Correct muon momentum for muon energy loss in passive targets if event is in passive target sample.
            debug() << " This is a passive target event, I will adjust muon energy" << endmsg;
            double vertexZ = event->interactionVertex()->position().z();
            Gaudi::LorentzVector muonfourVec;
            bool corrOK = m_energyCorrectionTool->getCorrectedEnergy( muonProng, muonPart, vertexZ, muonfourVec).isSuccess();
            counter( Form("EMuCorrOK_Tgt%d", targetID ) ) += corrOK;
            if( corrOK  )
              {
                debug() << "  Before Momentum = " << muonPart->momentumVec() << endmsg;
                debug() << "  After  Momentum = " << muonfourVec << endmsg;
                debug() << "  Delta           = " << ( muonfourVec - muonPart->momentumVec() ) << endmsg;
                vector<double> corrP;
                corrP.push_back( muonfourVec.Px() );
                corrP.push_back( muonfourVec.Py() );
                corrP.push_back( muonfourVec.Pz() );
                corrP.push_back( muonfourVec.E() );
                //event->setContainerDoubleData( "muon_corrected_p", corrP );
                //if( m_correctMuonEnergyToTarget )
                //  muonPart->setMomentumVec(muonfourVec);
              }
          }
      }

    //==============================for Machine Learning====================================//
    // Adjust the vertex position and muon momentum to the center of a segment where applicable
    if( has_vtx && has_muon && muonProng && muonPart && m_useDNN )
      {
        for (int targetid=1; targetid < 6; targetid++ ){
          Minerva::NuclearTarget* tgt = m_nukeTool->getNuclearTarget( targetid );
          verbose() << "***** targetID, z start, z center, z end = " << tgt->getTargetID() << " " << tgt->getZStart() << " " << tgt->getZCenter() << " " << tgt->getZEnd() << endmsg;
        }
        debug() << "Adjusting vertex position based on DNN prediction. " << endmsg;
        debug() << "Make sure that m_useDNN is true? " << m_useDNN << endmsg;
        int targetID = -999;
        int target_id = -999;

        vector<int> DNN_vtx_modules, DNN_vtx_planes;
        if( DNN_vtx_modules.size() > 0 ) DNN_vtx_modules.clear();
        if( DNN_vtx_planes.size() > 0 ) DNN_vtx_planes.clear();

        if( 1 <= ANN_segments[0] && ANN_segments[0] <= 174 ) //modified from 214 to 174 takes tracker region into consideration
          {
            error() << "ANN_segments = " << ANN_segments << endmsg;
            int DNN_vtx_plane;
            int DNN_vtx_module;
            double DNN_z_center;
            targetID = GetTargetFromSegment( ANN_segments[0], DNN_vtx_module, DNN_vtx_plane );
            error() << "first error targetID="<< targetID<< "ANN_segments=" << ANN_segments << "DNN_vtx_module = " << DNN_vtx_module << "DNN_vtx_plane=" << DNN_vtx_plane << endmsg;
            DNN_vtx_modules.push_back(DNN_vtx_module);
            DNN_vtx_planes.push_back(DNN_vtx_plane);
            // If the vertex fit did not converge (or vertex has only 1 track),
            // then change the vertex position at the position where the muon intersects the zCenter of the relevant passive target.
            if( targetID == -1 ){
              debug() << "NukeCC getPlaneID and get the z center" << endmsg;
              vector<const Minerva::DePlane*> planes = m_idDet->getDePlanes();
              for( vector<const Minerva::DePlane*>::const_iterator itPlane = planes.begin(); itPlane != planes.end(); ++itPlane ){

                if( DNN_vtx_module == (*itPlane)->getPlaneID().module() && DNN_vtx_plane == (int)( (*itPlane)->getPlaneID().plane()) ){
                  DNN_z_center = (*itPlane)->getZCenter();
                  debug() << "plane, module, z center = " << DNN_vtx_module << ", " << DNN_vtx_plane << ", " << DNN_z_center << endmsg;
                }
              }
            }
            else {
              debug() << "NukeCC getNuclearTarget and get the z center of target " << targetID << endmsg;
              Minerva::NuclearTarget* target = m_nukeTool->getNuclearTarget( targetID );
              if( 0 == target ){
                warning() << "Cannot project to target with targetID = " << targetID << " because the targetID does not exist" << endmsg;
                return 0;
              }
              //set the offset for the target positions here:
              if( targetID == 1 ) DNN_z_center = target->getZCenter() - 2.3;
              if( targetID == 2 ) DNN_z_center = target->getZCenter() - 3.5;
              if( targetID == 3 ) DNN_z_center = target->getZCenter();
              if( targetID == 4 ) DNN_z_center = target->getZCenter() - 3.4;
              if( targetID == 5 ) DNN_z_center = target->getZCenter() - 3.05;
              if( targetID == 6 ) DNN_z_center = 5305.3;
              debug() << "DNN z center = " << DNN_z_center << endmsg;
            }

            error() << "second error targetID="<< targetID<< "ANN_segments=" << ANN_segments << "DNN_vtx_module = " << DNN_vtx_module << "DNN_vtx_plane=" << DNN_vtx_plane << endmsg;
            //SmartRef<Minerva::Vertex> primVtx = event->interactionVertex();
            vector<double> MLVtx;
	    SmartRef<Minerva::Vertex> primVtx = new Minerva::Vertex(event->interactionVertex()->position());
            adjustDNNVertexIntoCenterOfSegment( primVtx, DNN_z_center, muonPart );
            //save the the second highest probability segment
            target_id = GetTargetFromSegment( DNN_segment, DNN_vtx_module, DNN_vtx_plane );
            DNN_vtx_modules.push_back(DNN_vtx_module);
            DNN_vtx_planes.push_back(DNN_vtx_plane);

            MLVtx.push_back(primVtx->position().x());
            MLVtx.push_back(primVtx->position().y());
            MLVtx.push_back(primVtx->position().z());

            // error() << "DNN_segment=" << DNN_segment << "DNN_vtx_module = " << DNN_vtx_module << "DNN_vtx_plane=" << DNN_vtx_plane << endmsg;
            event->setContainerIntData("ANN_vtx_modules", DNN_vtx_modules);
            event->setContainerIntData("ANN_vtx_planes", DNN_vtx_planes);
            
            event->setContainerDoubleData("ANN_vtx", MLVtx);
            warning() << " ANN_vtx filled " << endmsg;
          //}

        //Correct muon momentum for muon energy loss in passive targets if event is in passive target sample.
        //for Muon Energy correction by Huma
        debug() << " This is a machine learning predicted event, I will adjust muon energy" << endmsg;
        double vertexZ = primVtx->position().z();
        //double vertexZ = event->interactionVertex()->position().z();
        Gaudi::LorentzVector muonfourVec;
        bool corrOK = m_energyCorrectionTool->getCorrectedEnergy( muonProng, muonPart, vertexZ, muonfourVec).isSuccess();
        counter( Form("EMuCorrOK_Tgt%d", targetID ) ) += corrOK;
        if( corrOK  )
          {
            debug() << "  Before Momentum = " << muonPart->momentumVec() << endmsg;
            debug() << "  After  Momentum = " << muonfourVec << endmsg;
            debug() << "  Delta           = " << ( muonfourVec - muonPart->momentumVec() ) << endmsg;
            vector<double> corrP;
            corrP.push_back( muonfourVec.Px() );
            corrP.push_back( muonfourVec.Py() );
            corrP.push_back( muonfourVec.Pz() );
            corrP.push_back( muonfourVec.E() );
            event->setContainerDoubleData( "muon_corrected_p", corrP );
            //if( m_correctMuonEnergyToTarget )
            //  muonPart->setMomentumVec(muonfourVec);
         }
      }
    }


  //==========================================================================//
  //---------------------------------------------------
  //! Start using dEdX Tool
  //---------------------------------------------------
  debug() <<"Making particle hypotheses" << endmsg;
  ProngVect hadronProngs;

  //--------------------------------------------
  //! Create Particles and get hadronProngs
  //--------------------------------------------
  if( !createParticles( event, hadronProngs ) )
    info()<<"Could not create Particles"<<endmsg;

  //! Look and Tag with the ImprovedMichelTool
  //---------------------------------------------------
  bool hasImprovedMichels = has_lepton && ImprovedtagMichels( event, truth );
  event->filtertaglist()->setOrAddFilterTag("hasNoImprovedMichelElectrons", !hasImprovedMichels );
  verbose() << "Finished tagging Improved Michels " << endmsg;
  tagMichelElectrons( event );

  Minerva::ProngVect primary_prongs = event->primaryProngs();
  Minerva::ProngVect pionProngs;
  pionProngs.clear();
  for(unsigned int prong = 0; prong <primary_prongs.size(); ++prong) {
    if( primary_prongs[prong]->filtertaglist()->filterTagExists("PrimaryHadron")
        &&
        !primary_prongs[prong]->filtertaglist()->filterTagExists("PrimaryMuon") ) {
      pionProngs.push_back( primary_prongs[prong] );
    }
  }
  fillMatchedMichelBranches( event, hadronProngs );
  //---------------------------------------------------
  //! Let's count vertex types
  //---------------------------------------------------
  std::vector<int> vtx_types_all;
  std::vector<int> vtx_types_intime;
  std::vector<double> vtx_time_diff;

  const SmartRefVector<Minerva::Vertex> tmp_vertices = event->select<Minerva::Vertex>();
  for ( SmartRefVector<Minerva::Vertex>::const_iterator itVtx = tmp_vertices.begin(); itVtx != tmp_vertices.end(); ++itVtx ) {
    double vertex_time = (*itVtx)->getTime();//m_recoObjTimeTool->trackVertexTime((*itVtx));
    double time_diff = vertex_time - muonTime;
    vtx_time_diff.push_back(time_diff);
    if( time_diff > m_lowerTimeWindow && time_diff < m_upperTimeWindow )vtx_types_intime.push_back((*itVtx)->type());
    vtx_types_all.push_back((*itVtx)->type());
  }
  event->setContainerIntData("event_vertex_types", vtx_types_all);
  event->setContainerIntData("event_in_time_vertex_types", vtx_types_intime);
  event->setContainerDoubleData("event_vertex_time_diff", vtx_time_diff);
  //-----------------------------------------------------
  //! Let's count tracks
  //-----------------------------------------------------
  std::vector<double> track_energy_all;
  std::vector<double> track_energy_intime;
  std::vector<double> track_time_diff;
  std::vector<double >track_hit_time;
  std::vector<double >track_hit_energy;
  std::vector<int > track_nhits;
  const SmartRefVector<Minerva::Track> tmp_tracks = event->select<Minerva::Track>();
  for ( SmartRefVector<Minerva::Track>::const_iterator itTrack = tmp_tracks.begin(); itTrack != tmp_tracks.end(); ++itTrack ) {
    double track_time = m_recoObjTimeTool->trackVertexTime((*itTrack));
    double time_diff =  track_time - muonTime;
    track_time_diff.push_back(time_diff);
    if( time_diff > m_lowerTimeWindow && time_diff < m_upperTimeWindow )track_energy_intime.push_back((*itTrack)->type());
    track_energy_all.push_back((*itTrack)->type());

    //now lets had track by track hit time and charge information
    int n_hits = 0;
    //get clusters to get digits
    std::vector<Minerva::IDCluster*> clusters = (*itTrack)->idclusters();
    for(uint i=0;i<clusters.size();i++){
      SmartRefVector<Minerva::IDDigit> digits = clusters[i]->digits();
      for ( SmartRefVector<Minerva::IDDigit>::const_iterator itDig = digits.begin(); itDig != digits.end(); ++itDig ) {
        track_hit_time.push_back((*itDig)->time());
        track_hit_energy.push_back((*itDig)->normEnergy());
        n_hits++;
      }
    }//end cluster
    track_nhits.push_back(n_hits);
  }//end tracks
  event->setContainerDoubleData("event_tracks_energy",track_energy_all);
  event->setContainerDoubleData("event_in_time_tracks_energy",track_energy_intime);
  event->setContainerDoubleData("event_track_time_diff", track_time_diff);
  event->setContainerDoubleData("event_track_hit_time", track_hit_time);
  event->setContainerDoubleData("event_track_hit_energy", track_hit_energy);
  event->setContainerIntData("event_track_nhits",track_nhits);

  //--------------------------------------------------------------
  //! Momentum-by-range
  //--------------------------------------------------------------

  calcMomentumByRangeVars(event);

  const double muon_time = m_recoObjTimeTool->trackVertexTime(muonProng->minervaTracks()[0]);
  fillIsoProngInfo(event, muon_time);

  //-------------------------------------------------------
  //! Let's Get recoil values for various windows
  //-------------------------------------------------------
  std::cout << "Doing recoil variables" << std::endl;
  std::vector<double> recoil_summed_energy;
  std::vector<double> recoil_summed_energy_edge;
  std::vector<double> recoil_data_fraction;
  std::vector<int> recoil_lower_time_limit;
  std::vector<int> recoil_upper_time_limit;
  Minerva::IDClusterVect  muonClusters = muonProng->getAllIDClusters();
  Minerva::IDClusterVect mod_slice_clusters;
  for( Minerva::IDClusterVect::iterator itsc=sliceclusters.begin(); itsc!=sliceclusters.end(); itsc++ ) {
    bool nonMuonCluster = true;
    for( Minerva::IDClusterVect::iterator muit = muonClusters.begin(); muit!=muonClusters.end();muit++){
      if((*muit)==(*itsc)){
        nonMuonCluster = false;
        break;
      }
    }
    if(nonMuonCluster) mod_slice_clusters.push_back(*itsc);
  }

  for(int low=-20;low!=0;low++){
    for(int high=35;high!=0;high--){
      Minerva::IDClusterVect selectedclusters;
      Minerva::IDClusterVect selectedclustersedge;
      for( Minerva::IDClusterVect::iterator itsc=mod_slice_clusters.begin(); itsc!=mod_slice_clusters.end(); itsc++ ) {
        double time_diff = (*itsc)->time()-muonTime;
        if(time_diff>low && time_diff <high){
          selectedclusters.push_back(*itsc);
          //Special check for containment. 3 strip (<=3 or >=124)
          SmartRefVector<Minerva::IDDigit> clusDigits = (*itsc)->digits();
          SmartRefVector<Minerva::IDDigit>::iterator itDig;
          for(itDig = clusDigits.begin();itDig!=clusDigits.end();itDig++){
            if((*itDig)->strip()<=3 || (*itDig)->strip()>=124) selectedclustersedge.push_back(*itsc);
          }
        }
      }
      std::vector<const Minerva::IDCluster*> idclusters;
      idclusters.insert(idclusters.begin(), selectedclusters.begin(), selectedclusters.end());
      if(truth)recoil_data_fraction.push_back(TruthMatcher->getDataFraction(idclusters));
      else recoil_data_fraction.push_back(-1);
      //now apply constants
      double rec_energy = m_caloUtils->applyCalConsts( selectedclusters, "Default", false );
      double rec_energy_edge = m_caloUtils->applyCalConsts( selectedclustersedge, "Default", false );
      //now store
      recoil_summed_energy.push_back(rec_energy);
      recoil_lower_time_limit.push_back(low);
      recoil_upper_time_limit.push_back(high);
      recoil_summed_energy_edge.push_back(rec_energy_edge);


    }
  }
  //now store
  event->setContainerDoubleData("recoil_summed_energy",recoil_summed_energy);
  event->setContainerIntData("recoil_lower_time_limit",recoil_lower_time_limit);
  event->setContainerIntData("recoil_upper_time_limit",recoil_upper_time_limit);
  event->setContainerDoubleData("recoil_data_fraction",recoil_data_fraction);
  event->setContainerDoubleData("recoil_summed_energy_edge",recoil_summed_energy_edge);


  //=======================================
  // Blobbing Section
  //------------------
  // Create temporary blobs and prongs.
  //--
  // Get quanties from them but do not
  // save them in event store.
  //--
  // Used, Unused have a funny meaning.
  // All blob-analyzable clusters
  // are marked as Unused although some
  // of them are on tracks.
  //=======================================
  if( muonProng )
    {
      //filling nue style fuzz clusters.
      SmartRefVector<Minerva::IDCluster> NuEidClusters = m_electronUtils->NuEStyleLeptonFuzzClusters(event,muonProng);
      fillBlobBranches(event,NuEidClusters,NUEFUZZBLOB_PREFIX);

      // Removing the following section as there is no information saved from the block, but the `makeIsolatedIDBlobs` line affects the
      // histories of some clusters to Used affecting later blobbing. If it becomes desirable to add this information back to the output,
      // one must be careful to understand how changes to the objects in TES affect later parts of the code. There are efforts underway to
      // map this reconstruction path for reference.
      // -David Last 01/20/2022

      /*
      // get clusters for blobbing, which are all marked Unused
      SmartRefVector<Minerva::IDCluster> blobbingClusters = getBlobbingClusters( event, muonProng );

      // Blobber 3: Isolated Blobs
      vector<Minerva::IDBlob*> isoBlobs; //will have to delete the prongs/blobs
      m_primaryBlobProngTool->makeIsolatedIDBlobs( blobbingClusters, &isoBlobs );

      SmartRefVector<Minerva::IDCluster> isoBlobbedClusters = removeUsedClusters( blobbingClusters );
      // "Blobber" 4: remaining dispersed clusters in the QE recoil zone
      SmartRefVector<Minerva::IDCluster> dispBlobbedClusters = blobbingClusters;

      ##### The lines below were previously after line 2720 -David Last 01/20/2022 #####

      //color the clusters not used for specific blobs as general recoil
      SmartRefVector<Minerva::IDCluster> allBlobbedClusters = dispBlobbedClusters;
      allBlobbedClusters.insert( allBlobbedClusters.end(), isoBlobbedClusters.begin(), isoBlobbedClusters.end() );
      */

      // "Blobber" 5: the clusters used for the cc inclusive recoil energy
      //
      //get the clusters used for recoil
      SmartRefVector<Minerva::IDCluster> idClusters = getAnalyzableNonMuIDClusters( event, muonProng );
      SmartRefVector<Minerva::ODCluster> odClusters = getAnalyzableNonMuODClusters( event, muonProng );
      fillBlobBranches( event, idClusters, odClusters, RECOILBLOB_PREFIX );
      event->setIntData( "blob_" + RECOILBLOB_PREFIX + "_nBlobs", 1 );

    }//done blobbing



  //===============================================================================================


  //----------------------------------------------------------------
  //! Get a single proton coming out from the interaction vertex
  //! Now, let's look for the most energetic proton...
  //----------------------------------------------------------------
  std::cout << "On to protons " << std::endl;
  info()<<"Looking for a proton"<<endmsg;
  Minerva::ProngVect secondaryProtonProngs;
  SmartRef<Minerva::Prong> protonProng;
  SmartRef<Minerva::Particle> protonPart;

  bool hasProton = getProtonProngs( hadronProngs, secondaryProtonProngs, protonProng, protonPart );
  event->filtertaglist()->setOrAddFilterTag("hasProton", hasProton );
  event->setIntData("has_proton",hasProton);
  if( hasProton ) {
    info()<<"A proton was found"<<endmsg;
    counter("has_Proton")+=1;
    if( truth ) tagProtonProngTruth( event, protonProng );
  }
  debug() <<"Number of secondary proton prongs = "<< secondaryProtonProngs.size() << endmsg;

  //----------------------------------------------------
  // Assigning vectors for storing info from clusters
  // found at the end of the proton track
  //----------------------------------------------------
  std::vector<int> clustersFoundAtEndOfProtonProng;
  std::vector<double> maxDistanceClusterFromEndOfProtonProng;
  std::vector<int> numOfClustersAtEndOfProtonProng;
  std::vector<double> visEnergyOfClustersAtEndOfProtonProng;
  std::vector<double> calibEnergyOfClustersAtEndOfProtonProng;

  if( hasProton ) {
    //Make a Cone at the end of the proton track
    bool foundClusters = false;
    double max_cluster_distance = -9999.;
    Minerva::IDClusterVect idClustersInsideCone;
    double totalVisEAtProtonTrackEnd = 0;
    sc = m_masterAnaDevRecoUtils->makeConeAtEndProtonTrack( event, protonProng, muonProng, foundClusters, max_cluster_distance, idClustersInsideCone, m_fillTruthTG4ProtonCone );

    verbose() << "Number of clusters inside cone in MasterAnaDev.cpp (hasProton=true) " << idClustersInsideCone.size() << endmsg;
    verbose() << "Value of foundClusters (hasProton=true) " << foundClusters << endmsg;

    //----------------------------------------------------------------------------------
    // Find total visible energy of clusters present at end of Primary Proton track
    //----------------------------------------------------------------------------------
    for( Minerva::IDClusterVect::iterator itClust = idClustersInsideCone.begin(); itClust != idClustersInsideCone.end(); ++itClust ) {
      totalVisEAtProtonTrackEnd += (*itClust)->energy();
    }

    //-----------------------------------------------------------------------------------
    // Find total calibrated energy of clusters present at end of Primary Proton track
    //-----------------------------------------------------------------------------------
    double calibEAtProtonTrackEnd = m_caloUtils->applyCalConsts( idClustersInsideCone, "Default", false );

    //-------------------------------------------
    //Storing variables for the Primary Proton
    //-------------------------------------------
    clustersFoundAtEndOfProtonProng.push_back( (int)foundClusters );
    maxDistanceClusterFromEndOfProtonProng.push_back( max_cluster_distance );
    numOfClustersAtEndOfProtonProng.push_back( idClustersInsideCone.size() );
    visEnergyOfClustersAtEndOfProtonProng.push_back( totalVisEAtProtonTrackEnd );
    calibEnergyOfClustersAtEndOfProtonProng.push_back( calibEAtProtonTrackEnd );

    //-------------------------------------
    // Fill Primary Proton track length
    //-------------------------------------
    SmartRef<Minerva::Vertex> BackVertex;
    m_objectAssociator->getVertex_fromTrackBack( BackVertex, protonProng->minervaTracks()[0] );
    double protonTrackLength = -9999.;
    if( BackVertex ) {
      const Gaudi::XYZPoint protonTrackEndPos = BackVertex->position();
      protonTrackLength = m_mathTool->distanceBetween( vertexOfPhysEvent->position(), protonTrackEndPos );
      event->setDoubleData("proton_track_length", protonTrackLength);
      event->setDoubleData("proton_track_endx", protonTrackEndPos.x());
      event->setDoubleData("proton_track_endy", protonTrackEndPos.y());
      event->setDoubleData("proton_track_endz", protonTrackEndPos.z());
    }

    //----------------------------------------------------
    // Color the clusters belonging to the cone
    //----------------------------------------------------
    for( Minerva::IDClusterVect::iterator itClust = idClustersInsideCone.begin(); itClust != idClustersInsideCone.end(); ++itClust ) {
      if( 0 <= m_coneEnergyColor ) m_hitTaggerTool->applyColorTag( *itClust, m_coneEnergyColor );
    }

  } //End of if condition for hasProton

  if( secondaryProtonProngs.size()>0 ) {
    info() << "Secondary proton(s) found" << endmsg;
    for( ProngVect::iterator secProtonProng = secondaryProtonProngs.begin(); secProtonProng != secondaryProtonProngs.end(); ++secProtonProng ) {
      bool foundClusters = false;
      double max_cluster_distance = -9999.;
      Minerva::IDClusterVect idClustersInsideCone;
      double totalVisEAtProtonTrackEnd = 0;
      sc = m_masterAnaDevRecoUtils->makeConeAtEndProtonTrack( event, *secProtonProng, muonProng, foundClusters, max_cluster_distance, idClustersInsideCone, false );

      verbose() << "Number of clusters inside cone in MasterAnaDev.cpp (secondaryProtonProngs.size()>0) " << idClustersInsideCone.size() << endmsg;
      verbose() << "Value of foundClusters (secondaryProtonProngs.size()>0) " << foundClusters << endmsg;

      //-----------------------------------------------------------------------------------
      // Find total visible energy of clusters present at end of Secondary Proton track
      //-----------------------------------------------------------------------------------
      for( Minerva::IDClusterVect::iterator itClust = idClustersInsideCone.begin(); itClust != idClustersInsideCone.end(); ++itClust ) {
        totalVisEAtProtonTrackEnd += (*itClust)->energy();
      }

      //-------------------------------------------------------------------------------------
      // Find total calibrated energy of clusters present at end of Secondary Proton track
      //-------------------------------------------------------------------------------------
      double calibEAtProtonTrackEnd = m_caloUtils->applyCalConsts( idClustersInsideCone, "Default", false );

      //---------------------------------------------
      //Storing variables for each Secondary Proton
      //---------------------------------------------
      clustersFoundAtEndOfProtonProng.push_back( (int)foundClusters );
      maxDistanceClusterFromEndOfProtonProng.push_back( max_cluster_distance );
      numOfClustersAtEndOfProtonProng.push_back( idClustersInsideCone.size() );
      visEnergyOfClustersAtEndOfProtonProng.push_back( totalVisEAtProtonTrackEnd );
      calibEnergyOfClustersAtEndOfProtonProng.push_back( calibEAtProtonTrackEnd );

      //----------------------------------------------------
      // Color the clusters belonging to the cone
      //----------------------------------------------------
      for( Minerva::IDClusterVect::iterator itClust = idClustersInsideCone.begin(); itClust != idClustersInsideCone.end(); ++itClust ) {
        if( 0 <= m_coneEnergyColor ) m_hitTaggerTool->applyColorTag( *itClust, m_coneEnergyColor );
      }


    } //End of for loop over secondaryProtonProngs

    if(truth) tagSecondaryProtonProngTruth( event, secondaryProtonProngs);//tag secondaries!!

  } //End of if condition for secondaryProtonProngs

  //---------------------------------------------------------------
  //Storing variables for Primary+Secondary Protons in the event
  //---------------------------------------------------------------
  event->setContainerIntData("clusters_found_at_end_proton_prong", clustersFoundAtEndOfProtonProng);
  event->setContainerDoubleData("clusters_found_at_end_proton_prong_max_distance", maxDistanceClusterFromEndOfProtonProng);
  event->setContainerIntData("number_clusters_at_end_proton_prong", numOfClustersAtEndOfProtonProng);
  event->setContainerDoubleData("visE_clusters_at_end_proton_prong", visEnergyOfClustersAtEndOfProtonProng);
  event->setContainerDoubleData("calibE_clusters_at_end_proton_prong", calibEnergyOfClustersAtEndOfProtonProng);


  //----------------------------------------------------------------------------------------------------------------
  //! Blobbing Stuff
  //! JO has chosen to make VertexBlobProngs first, then add EM Blobs to MuonProng, then make IsoBlobProngs
  //! Reason for this is to pick out the recoil "near the Vertex" first, attach "muon-related brems/deltas"
  //! to the MuonProng second, then determine the "non-Vertex" recoil.
  //! Kenyi - please let me know if you disagree with this !
  //----------------------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------
  //! Create Vertex Blob Prongs
  //! User can choose to make VtxBlobProngs either "Filament Style" or "Multiple Radii Style"
  // V.V.I. - The vector of clusters 'idClusterVector' contains clusters "in time" with the muon
  // V.V.I. - Is this okay ? What if there is an out-of-time overlay activity right near the vertex
  // V.V.I. - INVESTIGATE !
  //----------------------------------------------------------------------------------------------------

  // Minerva::ProngVect newVtxBlobProngs;
  // if( (RecoilUtils->createVtxBlobProngs( event, newVtxBlobProngs )).isFailure() )
  // return StatusCode::FAILURE;

  debug() << "Creating vertex blob prongs" << endmsg;
  std::vector<Minerva::Prong*> *newVtxBlobProngs = new std::vector<Minerva::Prong*>;
  if( m_makeFilamentStyleVtxBlobProngs ) {
    sc = m_primaryBlobProngTool->makeFilamentStyleVtxBlobProngs( event, idClusterVector, vertexOfPhysEvent,
                                                                 m_maxSearchDistance, m_maxStartingDistance, m_maxAllowedSearchGap,
                                                                 m_maxSeparationBlobVertex, newVtxBlobProngs );

  } else if( m_makeMultipleRadiiStyleVtxBlobProngs ) {
    sc = m_primaryBlobProngTool->makeMultipleRadiiStyleVtxBlobProngs( event, idClusterVector, vertexOfPhysEvent,
                                                                      m_searchStepSize, m_numSearchRadii,
                                                                      m_maxSeparationBlobVertex, newVtxBlobProngs );

  } else {
    info() << "You have not specified a preference for making VtxBlobProngs, please specify one if you would like Vertex Blob Prongs !" << endmsg;
  }

  if( sc.isFailure() ) return sc;

  debug() << "Found " << newVtxBlobProngs->size() << " new Vertex Blob Prongs ....." << endmsg;

  //-----------------------------------------------------------------------
  //! You can add other conditions here, for examining the VtxBlobProngs
  //! before you start adding them to the evtMgr, user has control !
  //-----------------------------------------------------------------------
  if( newVtxBlobProngs->size() ) addObject( event, *newVtxBlobProngs );

  //---------------------------------------------------------------
  //! If there are Vertex Prongs, always promote them to primary
  //! V.V.I.- INVESTIGATE WHY THIS IS CAUSING CONFUSION
  //---------------------------------------------------------------
  /*Comment this part temporarily
    for( std::vector<Minerva::Prong*>::iterator p=newVtxBlobProngs->begin(); p!=newVtxBlobProngs->end(); ++p ) {
    event->promoteProngToPrimary( *p );
    }
  */

  //-----------------------------------------------------------------
  //! Color the clusters belonging to the Vertex Blob Prong
  //-----------------------------------------------------------------
  for( std::vector<Minerva::Prong*>::iterator itProng =newVtxBlobProngs->begin(); itProng != newVtxBlobProngs->end(); ++itProng ) {
    if( 0 <= m_vertexBlobProngColor ) m_hitTaggerTool->applyColorTag( *itProng, m_vertexBlobProngColor );
  }

  //-----------------------------------------------------------
  //! Delete vector newVtxBlobProngs but not their elements
  //-----------------------------------------------------------
  delete newVtxBlobProngs;

  //------------------------------------------------------------------------------
  //! Create Isolated 3D EM Blobs and add it to muonProng
  //! These attach the muon-related activity (brems/deltas) to the muon track
  //------------------------------------------------------------------------------
  if ( m_makeIsoMuonEMBlobs ) {
    debug()<<"Create Isolated 3D EM Blobs and add it to muon Prong"<<endmsg;
    m_masterAnaDevRecoUtils->makeIsoEMMuonBlobs( event, muonProng );
  }
  if ( m_makeFuzzMuonEMBlobs ) {
    debug()<<"Create Fuzz EM Blobs and add it to muon Prong"<<endmsg;
    m_masterAnaDevRecoUtils->makeFuzzEMMuonBlobs( event, muonProng );
  }

  //Temporary check of number of tracks in muon prong
  Minerva::TrackVect muontracks = muonProng->minervaTracks();
  debug()<<" tracks in muon prong = " << muontracks.size() << endmsg;

  //!Create primary Blob Prongs - NOT SURE IF WE SHOULD RUN THIS GENERAL METHOD NOW !
  // m_primaryBlobProngTool->makePrimaryBlobProngs( event );

  //------------------------------------------------------------------------------------------------------------------------------
  //! Create Shower Blob Prongs from the remaining "Unused" clusters
  //! They become the non-Vertex recoil Prong
  //! V.V.I. - We have noticed events (10200/69/521, 10203/54/601) where there is activity from overlay in the time slice.
  //! However this is out-of-time (i.e. -20, 35 ns) with the muon time.
  //! This activity will not be blobbed if we stick to "in time" clusters only.
  //! Hence we are resorting to collecting "ALL UNUSED" clusters (except XTalk & LowAct) in the event, for blobbing purposes.
  //! For event selection, we cut only on the blobs that are in time with the neutrino interaction we are interested in.
  //------------------------------------------------------------------------------------------------------------------------------
  debug() << "Create Shower Blob Prongs from remaining \"unused\" clusters" << endmsg;

  std::vector<Minerva::Prong*> *newIsoIdProngs = new std::vector<Minerva::Prong*>;

  Minerva::IDClusterVect  unusedClustsForBlobs = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

  sc = m_primaryBlobProngTool->makeShowerBlobProngs( event, unusedClustsForBlobs, newIsoIdProngs );
  if( sc.isFailure() ) return sc;

  //-----------------------------------------------------------------
  //! Color the clusters belonging to the non-Vertex recoil Prong
  //-----------------------------------------------------------------
  for( std::vector<Minerva::Prong*>::iterator itProng = newIsoIdProngs->begin(); itProng != newIsoIdProngs->end(); ++itProng ) {
    if( 0 <= m_isolatedBlobProngColor ) m_hitTaggerTool->applyColorTag( *itProng, m_isolatedBlobProngColor );

    debug() << "Prong Type " << (*itProng)->type() << endmsg;

    SmartRefVector<Minerva::IDCluster> prongClusts = (*itProng)->getAllIDClusters();
    for( SmartRefVector<Minerva::IDCluster>::iterator itClust = prongClusts.begin(); itClust != prongClusts.end(); ++itClust ) {
      debug() << "Prong Cluster Time " << (*itClust)->time() << endmsg;
    }
  }

  // ProngVect newIsoIdProngs;
  // if( (RecoilUtils->createIsoBlobProngs( event, newIsoIdProngs )).isFailure() )
  // return StatusCode::FAILURE;

  debug() << "Found " << newIsoIdProngs->size() << " new Isolated ID Blob Prongs ....." << endmsg;

  //--------------------------------------------------------------------------
  //! You can add other conditions here, for examining the ShowerBlobProngs
  //! before you start adding them to the evtMgr, user has control !
  //--------------------------------------------------------------------------
  if( newIsoIdProngs->size() ) addObject( event, *newIsoIdProngs );

  //Delete vector newIsoIdProngs but not their elements
  delete newIsoIdProngs;

  //!@todo cut on number of isolated blobs (<= 2)?
  //! I think the actual cut should be performed on the isolated blobs outside
  //! a certain (spheric?) region from the vertex to still be defined (30cm radius?)
  int n_nonvtx_iso_blobs_all = 0, n_nonvtx_iso_blobs = 0;
  double nonvtx_iso_blobs_energy_all = 0., nonvtx_iso_blobs_energy = 0.;
  double vtx_blobs_energy = 0.;
  double vtxisoblob_clusE_outsideradius = 0.;
  //  double muonTime = m_recoObjTimeTool->trackVertexTime(muonProng->minervaTracks()[0]);

  std::vector<double> vtxBlobEnergyInVtxBlobProng;
  std::vector<double> isoBlobEnergyInVtxBlobProng;
  std::vector<double> isoBlobDistanceInVtxBlobProng;
  std::vector<double> nonVtxIsoBlobEnergyInProng;
  std::vector<double> nonVtxIsoBlobDistanceInProng;
  std::vector<double> isoBlobEnergyOutsideRadiusInVtxBlobProng;
  std::vector<double> nonVtxIsoBlobStartXInProng;
  std::vector<double> nonVtxIsoBlobStartYInProng;
  std::vector<double> nonVtxIsoBlobStartZInProng;
  std::vector<double> nonVtxIsoBlobTimeInProng;
  std::vector<double> nonVtxIsoBlobTimeDiffInProng;
  std::vector<double> nonVtxIsoBlobLowestModuleXInProng;
  std::vector<double> nonVtxIsoBlobLowestModuleUInProng;
  std::vector<double> nonVtxIsoBlobLowestModuleVInProng;
  std::vector<double> nonVtxIsoBlobHighestModuleXInProng;
  std::vector<double> nonVtxIsoBlobHighestModuleUInProng;
  std::vector<double> nonVtxIsoBlobHighestModuleVInProng;
  std::vector<double> nonVtxIsoBlobEarliestHitTimeInProng;
  std::vector<double> nonVtxIsoBlobLatestHitTimeInProng;
  std::vector<double> nonVtxIsoBlobHighestHitEnergyInProng;
  std::vector<int> nonVtxIsoBlobNHitsInProng;
  std::vector<int> nonVtxIsoBlobParticlePDGInProng;
  std::vector<int> nonVtxIsoBlobPrimaryParticlePDGInProng;
  std::vector<double> nonVtxIsoBlobMatchedEnergyFractionInProng;
  std::vector<double> nonVtxIsoBlobDataEnergyFractionInProng;

  Minerva::ProngVect prongs = event->select<Minerva::Prong>( "Used:Unused", "All" );
  debug() << "Event has now " << prongs.size() << " prongs" << endmsg;

  for( Minerva::ProngVect::iterator p=prongs.begin(); p!=prongs.end(); p++ ) {
    debug() << "Prong Creation Signature = "<< (*p)->creationSignature() << endmsg;
    //Get calibrated vertex energy
    if( ((*p)->creationSignature()== PrimaryBlobProngExtraDataDefs::InteractionVtxBlobProngName()) ||
        ((*p)->creationSignature()== PrimaryBlobProngExtraDataDefs::StartPointVtxBlobProngName()) ){

      double calE = m_caloUtils->applyCalConsts( *p, "Default", false );
      vtx_blobs_energy += calE;
      debug() << "Found a " << (*p)->creationSignature() << ", nblobs: " << (*p)->idblobs().size() << ", energy: " << calE << " MeV" << endmsg;

      // loop over blobs in prong
      Minerva::IDBlobVect idBlobs = (*p)->idblobs();
      for( Minerva::IDBlobVect::iterator itBlob=idBlobs.begin(); itBlob!=idBlobs.end(); itBlob++ ){
        if( (*itBlob)->patRecHistory() == Minerva::IDBlob::VertexBlobPatRec ) {
          double vtxE = m_caloUtils->applyCalConsts( *itBlob, "Default", false );
          vtxBlobEnergyInVtxBlobProng.push_back(vtxE);
        }
        else if( (*itBlob)->patRecHistory() == Minerva::IDBlob::IsolatedIDBlobPatRec ) {
          double isoE = m_caloUtils->applyCalConsts( *itBlob, "Default", false );
          isoBlobEnergyInVtxBlobProng.push_back(isoE);

          const Gaudi::XYZPoint blob_pos = (*itBlob)->startPoint();
          const Gaudi::XYZPoint vtx_pos = vertexOfPhysEvent->position();

          double distance = m_mathTool->distanceBetween(blob_pos, vtx_pos);
          isoBlobDistanceInVtxBlobProng.push_back(distance);
          /*
            const Gaudi::XYZPoint blob_pos_us = m_blobCreatorUtils->calcUpstreamBlobStartPoint(*itBlob);
            const Gaudi::XYZPoint blob_pos_ds = m_blobCreatorUtils->calcDownstreamBlobStartPoint(*itBlob);
            const Gaudi::XYZPoint vtx_pos = vertexOfPhysEvent->position();
            double distance1 = m_mathTool->distanceBetween(blob_pos_us, vtx_pos);
            double distance2 = m_mathTool->distanceBetween(blob_pos_ds, vtx_pos);
            double distance = (distance1<distance2)? distance1 : distance2;
            isoBlobDistanceInVtxBlobProng.push_back(distance);
            debug()<<"Iso vertex blob US (x,y,z) = "<< blob_pos_us.x() << ", "<< blob_pos_us.y()<<", "<< blob_pos_us.z()<< endmsg;
            debug()<<"Iso vertex blob DS (x,y,z) = "<< blob_pos_ds.x() << ", "<< blob_pos_ds.y()<<", "<< blob_pos_ds.z()<< endmsg;
            debug()<<"vertex position (x,y,z) = "<< vtx_pos.x() << ", "<< vtx_pos.y()<<", "<< vtx_pos.z()<< endmsg;
            debug()<<"distance iso vertex blob US from vertex="<< distance1<<endmsg;
            debug()<<"distance iso vertex blob DS from vertex ="<< distance2<<endmsg;
            debug()<<"distance iso vertex blob minimum from vertex="<< distance<<endmsg;
            debug()<<"distance iso vertex blob default start point ="<< (*itBlob)->startPoint()<<endmsg;
          */

          debug()<< "Vertex Iso Blob centroid position X,Z = " << ((*itBlob)->energyCentroidXZ()).first<< ", " << ((*itBlob)->energyCentroidXZ()).second <<endmsg;
          debug()<< "Vertex Iso Blob centroid position U,Z = " << (*itBlob)->energyCentroidUZ().first<< ", " << (*itBlob)->energyCentroidXZ().second <<endmsg;
          debug()<< "Vertex Iso Blob centroid position V,Z = " << (*itBlob)->energyCentroidVZ().first<< ", " << (*itBlob)->energyCentroidXZ().second <<endmsg;
          //! Also, calculate the energy of the clusters from these blobs that are outside
          //! the sphere of radius: m_searchStepSize
          SmartRefVector<Minerva::IDCluster> vtxisoblobClusters = (*itBlob)->clusters();
          SmartRefVector<Minerva::IDCluster> vtxisoblobClusters_outside;
          SmartRefVector<Minerva::IDCluster>::iterator itCluster = vtxisoblobClusters.begin();
          double minclus_distance = 1e6;
          for (; itCluster!=vtxisoblobClusters.end(); itCluster++ ) {
            double dist = m_mathTool->twoDDistance(*itCluster,vertexOfPhysEvent->position());
            if (minclus_distance>dist) minclus_distance = dist;
            if (dist> m_searchStepSize){
              vtxisoblobClusters_outside.push_back(*itCluster);
            }
          }
          debug()<<"Distance vtx iso blob to vertex = "<< distance<< endmsg;
          debug()<<"Distance iso vertex blob min cluster = " << minclus_distance << endmsg;
          double clusE_outside = m_caloUtils->applyCalConsts( vtxisoblobClusters_outside, "Default", false );
          isoBlobEnergyOutsideRadiusInVtxBlobProng.push_back(clusE_outside);
          vtxisoblob_clusE_outsideradius += clusE_outside;
        }
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------
    //! Get calibrated non-vertex energy for ALL the isolated blobs
    //! Here we have blobbed up "ALL UNUSED" clusters (except XTalk & LowAct) in the event, in time or out of time w.r.t. muon
    //! So activity from overlay will also appear as blobs, if they are in the same time slice.
    //------------------------------------------------------------------------------------------------------------------------------
    if( (*p)->creationSignature() == PrimaryBlobProngExtraDataDefs::IsolatedIDBlobProngName() ) {
      n_nonvtx_iso_blobs_all += (*p)->idblobs().size();
      double calE = m_caloUtils->applyCalConsts( *p, "Default", false );
      nonvtx_iso_blobs_energy_all += calE;

      //------------------------------------------------------------------------------------------------------------------
      //! Loop over Isolated Blobs in prong that are in time with the neutrino interaction that we are interested in.
      //! Get their calibrated energy.
      //! For event selection purposes, we cut ONLY on these blobs
      //------------------------------------------------------------------------------------------------------------------
      Minerva::IDBlobVect isoBlobs = (*p)->idblobs();
      for( Minerva::IDBlobVect::iterator itBlob=isoBlobs.begin(); itBlob!=isoBlobs.end(); itBlob++ ) {
        double time_diff = (*itBlob)->time() - muonTime;
        // get calibrated energy of isolated blob prong
        //@channel = "Defaul"
        //@applyPolyLine = false
        if( time_diff > m_lowerTimeWindow && time_diff < m_upperTimeWindow ) {
          n_nonvtx_iso_blobs++;
          double calE_inTime = m_caloUtils->applyCalConsts( *itBlob, "Default", false );
          nonvtx_iso_blobs_energy += calE_inTime;
          nonVtxIsoBlobEnergyInProng.push_back(calE_inTime);

          const Gaudi::XYZPoint blob_pos = (*itBlob)->startPoint();
          const Gaudi::XYZPoint vtx_pos = vertexOfPhysEvent->position();
          double distance = m_mathTool->distanceBetween(blob_pos, vtx_pos);
          /*
            const Gaudi::XYZPoint blob_pos_us = m_blobCreatorUtils->calcUpstreamBlobStartPoint(*itBlob);
            const Gaudi::XYZPoint blob_pos_ds = m_blobCreatorUtils->calcDownstreamBlobStartPoint(*itBlob);

            double distance1 = m_mathTool->distanceBetween(blob_pos_us, vtx_pos);
            double distance2 = m_mathTool->distanceBetween(blob_pos_ds, vtx_pos);
            double distance = (distance1<distance2)? distance1 : distance2;
          */
          //double distance = sqrt( pow(blob_pos.x()-vtx_pos.x(),2) + pow(blob_pos.y()-vtx_pos.y(),2) + pow(blob_pos.z()-vtx_pos.z(),2) );
          nonVtxIsoBlobDistanceInProng.push_back(distance);
          nonVtxIsoBlobStartXInProng.push_back(blob_pos.X());
          nonVtxIsoBlobStartYInProng.push_back(blob_pos.Y());
          nonVtxIsoBlobStartZInProng.push_back(blob_pos.Z());
          nonVtxIsoBlobTimeInProng.push_back((*itBlob)->time());
          nonVtxIsoBlobTimeDiffInProng.push_back(time_diff);
          nonVtxIsoBlobLowestModuleXInProng.push_back((*itBlob)->moduleLowX());
          nonVtxIsoBlobLowestModuleUInProng.push_back((*itBlob)->moduleLowU());
          nonVtxIsoBlobLowestModuleVInProng.push_back((*itBlob)->moduleLowV());
          nonVtxIsoBlobHighestModuleXInProng.push_back((*itBlob)->moduleHighX());
          nonVtxIsoBlobHighestModuleUInProng.push_back((*itBlob)->moduleHighU());
          nonVtxIsoBlobHighestModuleVInProng.push_back((*itBlob)->moduleHighV());
          verbose()<<"Distance non-vtx iso blob-vertex = "<< distance<<endmsg;

          Minerva::IDDigitVect blobhits = (*itBlob)->getAllDigits();
          std::vector<const Minerva::IDCluster*> blobidclusters;
          SmartRefVector<Minerva::IDCluster> clustersblob = (*itBlob)->clusters();
          blobidclusters.insert(blobidclusters.begin(), clustersblob.begin(), clustersblob.end());

          double earlyhit = 1e20;
          double latehit = -1e20;
          double highesthitenergy = 0;
          int nblobhits = 0;
          for(Minerva::IDDigitVect::iterator itBlobHits=blobhits.begin(); itBlobHits!=blobhits.end();itBlobHits++){
            if( (*itBlobHits)->calTime() < earlyhit) earlyhit =(*itBlobHits)->time();
            if( (*itBlobHits)->calTime() > latehit)  latehit =(*itBlobHits)->time();
            if( (*itBlobHits)->normEnergy() > highesthitenergy) highesthitenergy=(*itBlobHits)->normEnergy();
            nblobhits+=1;
          }
          nonVtxIsoBlobEarliestHitTimeInProng.push_back(earlyhit);
          nonVtxIsoBlobLatestHitTimeInProng.push_back(latehit);
          nonVtxIsoBlobHighestHitEnergyInProng.push_back(highesthitenergy);
          nonVtxIsoBlobNHitsInProng.push_back(nblobhits);

          //Match truth
          Minerva::IDBlob myblob = (*(*itBlob).data());
          std::vector<double> bestPart = tagBlobTruth(myblob);
          double data_fraction = TruthMatcher->getDataFraction(blobidclusters);
          nonVtxIsoBlobParticlePDGInProng.push_back(bestPart[0]);
          nonVtxIsoBlobPrimaryParticlePDGInProng.push_back(bestPart[1]);
          nonVtxIsoBlobMatchedEnergyFractionInProng.push_back(bestPart[2]);
          nonVtxIsoBlobDataEnergyFractionInProng.push_back(data_fraction);

        }
      }
    }
  } //End of for loop over Prongs

  verbose()<<"nonvtx_iso_blobs_energy = " << nonvtx_iso_blobs_energy << endmsg;
  verbose()<<"nonvtx_iso_blobs_energy_all = " << nonvtx_iso_blobs_energy_all << endmsg;

  event->setDoubleData("vtx_blobs_energy", vtx_blobs_energy);
  event->setDoubleData("vtx_iso_blobs_energy_outside_radius", vtxisoblob_clusE_outsideradius);
  event->setContainerDoubleData("vtx_blobs_vtx_energy_in_prong", vtxBlobEnergyInVtxBlobProng);
  event->setContainerDoubleData("vtx_blobs_iso_energy_in_prong", isoBlobEnergyInVtxBlobProng);
  event->setContainerDoubleData("vtx_blobs_iso_energy_clusters_outside_radius_in_prong", isoBlobEnergyOutsideRadiusInVtxBlobProng);
  event->setContainerDoubleData("vtx_blobs_iso_distance_in_prong", isoBlobDistanceInVtxBlobProng);
  event->setIntData("n_nonvtx_iso_blobs", n_nonvtx_iso_blobs);
  event->setIntData("n_nonvtx_iso_blobs_all", n_nonvtx_iso_blobs_all);
  event->setDoubleData("nonvtx_iso_blobs_energy", nonvtx_iso_blobs_energy);
  event->setDoubleData("nonvtx_iso_blobs_energy_all", nonvtx_iso_blobs_energy_all);
  event->setContainerDoubleData("nonvtx_iso_blobs_energy_in_prong", nonVtxIsoBlobEnergyInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_distance_in_prong", nonVtxIsoBlobDistanceInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_start_position_x_in_prong",nonVtxIsoBlobStartXInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_start_position_y_in_prong",nonVtxIsoBlobStartYInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_start_position_z_in_prong",nonVtxIsoBlobStartZInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_time_in_prong",nonVtxIsoBlobTimeInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_time_difference_in_prong",nonVtxIsoBlobTimeDiffInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_lowest_module_x_in_prong",nonVtxIsoBlobLowestModuleXInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_lowest_module_u_in_prong",nonVtxIsoBlobLowestModuleUInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_lowest_module_v_in_prong",nonVtxIsoBlobLowestModuleVInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_highest_module_x_in_prong",nonVtxIsoBlobHighestModuleXInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_highest_module_u_in_prong",nonVtxIsoBlobHighestModuleUInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_highest_module_v_in_prong",nonVtxIsoBlobHighestModuleVInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_earliest_hit_time_in_prong",nonVtxIsoBlobEarliestHitTimeInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_latest_hit_time_in_prong",nonVtxIsoBlobLatestHitTimeInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_highest_hit_energy_in_prong",nonVtxIsoBlobHighestHitEnergyInProng);
  event->setContainerIntData("nonvtx_iso_blobs_n_hits_in_prong",nonVtxIsoBlobNHitsInProng);
  event->setContainerIntData("nonvtx_iso_blobs_particle_pdg_in_prong",nonVtxIsoBlobParticlePDGInProng);
  event->setContainerIntData("nonvtx_iso_blobs_primary_particle_pdg_in_prong",nonVtxIsoBlobPrimaryParticlePDGInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_matched_energy_fraction_in_prong",nonVtxIsoBlobMatchedEnergyFractionInProng);
  event->setContainerDoubleData("nonvtx_iso_blobs_data_energy_fraction_in_prong",nonVtxIsoBlobDataEnergyFractionInProng);

  //-------------------------------------------------------------------------------------
  //! Look at the remaining ID clusters in the Tracker/ECAL regions (dispersed energy)
  //! in the [-20,35] ns time window
  //-------------------------------------------------------------------------------------
  Minerva::IDClusterVect unusedIDClusters = m_masterAnaDevRecoUtils->getInTimeAnalyzableIDClusters( event, muonProng );
  double iddisEnergy = m_caloUtils->applyCalConsts( unusedIDClusters, "Default", false );

  event->setDoubleData("dis_id_energy", iddisEnergy);

  //-------------------------------------------------------------------------------------
  //! Look at all ID clusters in the Tracker/ECAL regions
  //! in the [-20,35] ns time window which are not associated with the muon
  //! and do not occur inside various vertex regions
  //-------------------------------------------------------------------------------------

  Minerva::IDClusterVect nonMuonNonVtxClusters0mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonNonVtxIDClusters( event, muonProng, 0);
  Minerva::IDClusterVect nonMuonNonVtxClusters50mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonNonVtxIDClusters( event, muonProng, 50);
  Minerva::IDClusterVect nonMuonNonVtxClusters100mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonNonVtxIDClusters( event, muonProng, 100);
  Minerva::IDClusterVect nonMuonNonVtxClusters150mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonNonVtxIDClusters( event, muonProng, 150);
  Minerva::IDClusterVect nonMuonNonVtxClusters200mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonNonVtxIDClusters( event, muonProng, 200);
  Minerva::IDClusterVect nonMuonNonVtxClusters250mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonNonVtxIDClusters( event, muonProng, 250);
  Minerva::IDClusterVect nonMuonNonVtxClusters300mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonNonVtxIDClusters( event, muonProng, 300);

  double nonMuonNonVtxEnergy0mm = m_caloUtils->applyCalConsts( nonMuonNonVtxClusters0mm, "Default", false );
  double nonMuonNonVtxEnergy50mm = m_caloUtils->applyCalConsts( nonMuonNonVtxClusters50mm, "Default", false );
  double nonMuonNonVtxEnergy100mm = m_caloUtils->applyCalConsts( nonMuonNonVtxClusters100mm, "Default", false );
  double nonMuonNonVtxEnergy150mm = m_caloUtils->applyCalConsts( nonMuonNonVtxClusters150mm, "Default", false );
  double nonMuonNonVtxEnergy200mm = m_caloUtils->applyCalConsts( nonMuonNonVtxClusters200mm, "Default", false );
  double nonMuonNonVtxEnergy250mm = m_caloUtils->applyCalConsts( nonMuonNonVtxClusters250mm, "Default", false );
  double nonMuonNonVtxEnergy300mm = m_caloUtils->applyCalConsts( nonMuonNonVtxClusters300mm, "Default", false );

  event->setDoubleData("recoil_energy_nonmuon_nonvtx0mm",nonMuonNonVtxEnergy0mm);
  event->setDoubleData("recoil_energy_nonmuon_nonvtx50mm",nonMuonNonVtxEnergy50mm);
  event->setDoubleData("recoil_energy_nonmuon_nonvtx100mm",nonMuonNonVtxEnergy100mm);
  event->setDoubleData("recoil_energy_nonmuon_nonvtx150mm",nonMuonNonVtxEnergy150mm);
  event->setDoubleData("recoil_energy_nonmuon_nonvtx200mm",nonMuonNonVtxEnergy200mm);
  event->setDoubleData("recoil_energy_nonmuon_nonvtx250mm",nonMuonNonVtxEnergy250mm);
  event->setDoubleData("recoil_energy_nonmuon_nonvtx300mm",nonMuonNonVtxEnergy300mm);


  Minerva::IDClusterVect nonMuonVtxClusters0mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonVtxIDClusters( event, muonProng, 0);
  Minerva::IDClusterVect nonMuonVtxClusters50mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonVtxIDClusters( event, muonProng, 50);
  Minerva::IDClusterVect nonMuonVtxClusters100mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonVtxIDClusters( event, muonProng, 100);
  Minerva::IDClusterVect nonMuonVtxClusters150mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonVtxIDClusters( event, muonProng, 150);
  Minerva::IDClusterVect nonMuonVtxClusters200mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonVtxIDClusters( event, muonProng, 200);
  Minerva::IDClusterVect nonMuonVtxClusters250mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonVtxIDClusters( event, muonProng, 250);
  Minerva::IDClusterVect nonMuonVtxClusters300mm = m_masterAnaDevRecoUtils->getInTimeAnalyzableNonMuonVtxIDClusters( event, muonProng, 300);
  debug()<<"first try: nonMuonVtxClusters0mm size: "<<nonMuonVtxClusters0mm.size()<<endmsg;

  double nonMuonVtxEnergy0mm = m_caloUtils->applyCalConsts( nonMuonVtxClusters0mm, "Default", false );
  double nonMuonVtxEnergy50mm = m_caloUtils->applyCalConsts( nonMuonVtxClusters50mm, "Default", false );
  double nonMuonVtxEnergy100mm = m_caloUtils->applyCalConsts( nonMuonVtxClusters100mm, "Default", false );
  double nonMuonVtxEnergy150mm = m_caloUtils->applyCalConsts( nonMuonVtxClusters150mm, "Default", false );
  double nonMuonVtxEnergy200mm = m_caloUtils->applyCalConsts( nonMuonVtxClusters200mm, "Default", false );
  double nonMuonVtxEnergy250mm = m_caloUtils->applyCalConsts( nonMuonVtxClusters250mm, "Default", false );
  double nonMuonVtxEnergy300mm = m_caloUtils->applyCalConsts( nonMuonVtxClusters300mm, "Default", false );

  event->setDoubleData("recoil_energy_nonmuon_vtx0mm",nonMuonVtxEnergy0mm);
  event->setDoubleData("recoil_energy_nonmuon_vtx50mm",nonMuonVtxEnergy50mm);
  event->setDoubleData("recoil_energy_nonmuon_vtx100mm",nonMuonVtxEnergy100mm);
  event->setDoubleData("recoil_energy_nonmuon_vtx150mm",nonMuonVtxEnergy150mm);
  event->setDoubleData("recoil_energy_nonmuon_vtx200mm",nonMuonVtxEnergy200mm);
  event->setDoubleData("recoil_energy_nonmuon_vtx250mm",nonMuonVtxEnergy250mm);
  event->setDoubleData("recoil_energy_nonmuon_vtx300mm",nonMuonVtxEnergy300mm);

  //---------------------------------------------------
  // Color the clusters in the dispersed energy
  //---------------------------------------------------
  for( Minerva::IDClusterVect::iterator itClust = unusedIDClusters.begin(); itClust != unusedIDClusters.end(); ++itClust ) {
    if( 0 <= m_dispersedEnergyColor ) m_hitTaggerTool->applyColorTag( *itClust, m_dispersedEnergyColor );
  }

  //---------------------------------------------------
  //! use the neutron blobs for computation
  //---------------------------------------------------
  SmartRef<Minerva::Particle> neutronPart = m_neutronBlobUtils->computeNeutronFromHydrogen( muonPart );

  if( m_doNeutron )
    {
      debug()<<"Do Neutron and clear neutron data map"<<endmsg;
      m_neutronBlobUtils->ClearMap();
      if(truth)
        {
          m_anaBlobUtils->ResetTraj();
          m_anaBlobUtils->SetTraj();
        }


      //Get all clusters not on a track
      debug()<<"Create untrackedClusters Vector"<<endmsg;
      std::set< SmartRef<Minerva::IDCluster> > trackedClusterSets;
      //Add  clusters associated with the tracks on vertex and hadronProngs
      SmartRefVector<Minerva::Track> vertexTracks = event->interactionVertex()->getTracks();
      debug()<<"Start Looping"<<endmsg;
      for( SmartRefVector<Minerva::Track>::iterator itTrack = vertexTracks.begin(); itTrack != vertexTracks.end();++itTrack )
        {
          std::vector< Minerva::IDCluster* > clustersOnTrack = (*itTrack)->idclusters();
          for (std::vector< Minerva::IDCluster* >::iterator itC = clustersOnTrack.begin(); itC != clustersOnTrack.end(); ++itC )  trackedClusterSets.insert( *itC );
        }
      //Add clusters on the hadronProngs
      for (SmartRefVector<Minerva::Prong>::iterator itProng = hadronProngs.begin(); itProng != hadronProngs.end(); ++itProng)
        {
          SmartRefVector<Minerva::IDCluster> allIDClusters = (*itProng)->getAllIDClusters();
          for ( SmartRefVector<Minerva::IDCluster>::iterator itC = allIDClusters.begin(); itC != allIDClusters.end(); ++itC ) trackedClusterSets.insert( *itC );
        }
      Minerva::IDClusterVect trackedClusters( trackedClusterSets.begin(), trackedClusterSets.end() );
      Minerva::IDClusterVect untrackedClusters;

      for ( Minerva::IDClusterVect::iterator itC = nonMuonNonVtxClusters0mm.begin(); itC != nonMuonNonVtxClusters0mm.end(); ++itC )
        {
          bool duplicate = std::find( trackedClusters.begin(), trackedClusters.end() , *itC ) != trackedClusters.end();
          if (!duplicate) untrackedClusters.push_back( *itC );
        }



      //for( SmartRefVector<Minerva::Track>::iterator itTrack = vertexTracks.begin(); itTrack != vertexTracks.end(); ++itTrack )
      //{
      //  //loop through each track attached to the vertex and remove duplicate from untrackedClusters
      //  Minerva::IDClusterVect::iterator itCluster;
      //  std::vector< Minerva::IDCluster* > clustersOnTrack = (*itTrack)->idclusters();
      //  debug()<<"Cluster and iterator created"<<endmsg;
      //  for( itCluster = untrackedClusters.begin(); itCluster != untrackedClusters.end(); )
      //  {
      //    bool duplicate = std::find( clustersOnTrack.begin(), clustersOnTrack.end(), *itCluster) != clustersOnTrack.end() ;
      //    (duplicate)? debug()<<"Duplicate"<<endmsg : debug()<<"Not duplicate"<<endmsg;
      //    if (duplicate) untrackedClusters.erase( itCluster );
      //    else ++itCluster;
      //  }
      //}
      debug()<<"clusters are untracked"<<endmsg;

      //Get all untracked clusters that satisfy some criteria
      //far away from muon
      //must be in tracker
      //trackable or heavyionizing
      SmartRefVector<Minerva::IDCluster> nonMuonClusters;
      bool cleanTrack = true;
      if(cleanTrack)
        {
          for( SmartRefVector<Minerva::IDCluster>::iterator itC = untrackedClusters.begin(); itC != untrackedClusters.end(); ++itC )
            {
              bool inTracker =  ( (*itC)->subdet() == Minerva::IDCluster::Tracker);
              if ( !inTracker &&  m_doNeutronInTrackerOnly) continue;

              bool isTrackable = ( (*itC)->type() == Minerva::IDCluster::Trackable );
              bool isHeavyIonizing = ( (*itC)->type() == Minerva::IDCluster::HeavyIonizing );
              bool closeToMuon = m_neutronBlobUtils->isClusterCloseToMuon(muonPart, (*itC), m_muonAxisCylinderRadius);

              bool saveCluster = (!closeToMuon && isTrackable ) || isHeavyIonizing;
              if (saveCluster) nonMuonClusters.push_back( *itC );
            }
        }
      else nonMuonClusters = untrackedClusters;

      //Clear vector and pass on the event if there are too many clusters
      if( nonMuonClusters.size() < m_maxRecoCluster ) event->setIntData("PassNeutronMaxClusterCut", 1 );
      else nonMuonClusters.clear();


      std::vector<Minerva::MinervaRecoBase::History> nonMuonClusterHistories;
      for( unsigned int i = 0 ; i < nonMuonClusters.size(); i++)
        {
          nonMuonClusterHistories.push_back( nonMuonClusters[i] -> history() );
          nonMuonClusters[i]->setHistory(  Minerva::MinervaRecoBase::Unused );
        }


      SmartRefVector<Minerva::IDBlob> AllIDBlobs = m_neutronBlobRecoTool->findTrackableCandidates( nonMuonClusters, 100, 0.0 ); //form all blobs
      SmartRefVector<Minerva::IDBlob> AllIDContiguousBlobs = m_neutronBlobRecoTool->combineContiguous2DBlobs( AllIDBlobs ); //form all blobs

      int AllIDBlobsNClusters = 0, AllIDContiguousBlobsNClusters = 0;
      for( SmartRefVector<Minerva::IDBlob>::iterator itall = AllIDBlobs.begin(); itall != AllIDBlobs.end(); ++itall ) AllIDBlobsNClusters+= (*itall)->nclusters();
      for( SmartRefVector<Minerva::IDBlob>::iterator itall = AllIDContiguousBlobs.begin(); itall != AllIDContiguousBlobs.end(); ++itall ) AllIDContiguousBlobsNClusters+= (*itall)->nclusters();


      event->setIntData("NneutronClusters0",untrackedClusters.size());
      event->setIntData("NneutronClusters1",nonMuonClusters.size() );
      event->setIntData("NneutronClusters2",AllIDBlobsNClusters  );
      event->setIntData("NneutronClusters3",AllIDContiguousBlobsNClusters );


      debug()<<"nonMuonNonVtxClusters0mm size: "<<nonMuonNonVtxClusters0mm.size()<<endmsg;
      debug()<<"nonMuonClusters size: "<<nonMuonClusters.size()<<endmsg;
      debug()<<"AllIDBlobs size: "<<AllIDBlobs.size()<<endmsg;
      debug()<<"AllIDBlobs clusters size: "<<AllIDBlobsNClusters<<endmsg;
      debug()<<"AllIDContinugousBlobs size: "<<AllIDContiguousBlobs.size()<<endmsg;
      debug()<<"AllIDContiguousBlobs clusters size: "<<AllIDContiguousBlobsNClusters<<endmsg;

      debug()<<"Resetting History"<<endmsg;
      for( unsigned int i = 0 ; i < nonMuonClusters.size(); i++)
        {
          nonMuonClusters[i]->setHistory(  nonMuonClusterHistories[i] );
        }



      std::vector<NeutronDataStruct::BlobInfo> BlobsInfo; //create blobinfo data container
      debug()<<"Created BlobInfo"<<endmsg;
      Gaudi::XYZTVector neutron_vtx( vertexOfPhysEvent->position().X(), vertexOfPhysEvent->position().Y(), vertexOfPhysEvent->position().Z(), vertexOfPhysEvent->getTime());
      debug()<<"neutron_vtx set at ("<<neutron_vtx.X() <<", "<< neutron_vtx.Y() <<", " << neutron_vtx.Z()<<")"<<endmsg;
      for( SmartRefVector<Minerva::IDBlob>::iterator it = AllIDContiguousBlobs.begin(); it!=AllIDContiguousBlobs.end(); it++ )
        {
          debug()<<"Fill data for blob "<<it-AllIDContiguousBlobs.begin()<<endmsg;
          m_hitTaggerTool->applyColorTag( *it, m_neutronBlobColor );

          bool blobIsUpstream = (*it)->position().Z() < vertexOfPhysEvent->position().Z();

          NeutronDataStruct::BlobInfo data;
          data.id = it -  AllIDContiguousBlobs.begin();
          m_neutronBlobUtils->FillBlobDataIndependent( data, *it );
          m_neutronBlobUtils->FillBlobDataVtx( data, *it, neutron_vtx );
          m_neutronBlobUtils->FillBlobDataPart(data, *it, neutronPart );
          m_neutronBlobUtils->FillBlobDataMuon(data, *it, muonPart );

          double ax,ay, residual;
          Gaudi::XYZPoint upstream, downstream;
          bool fitok = m_blobFitTool->Fit( (*it)->clusters(), ax, ay, upstream, downstream, residual );
          if( fitok )
            {
              blobIsUpstream = (downstream- vertexOfPhysEvent->position() ).R() < (upstream - vertexOfPhysEvent->position() ).R();
              double time_begin =std::min(data.blobBegPos.T(), data.blobEndPos.T() );
              double time_end = std::max(data.blobBegPos.T(), data.blobEndPos.T() );

              if( blobIsUpstream )
                {
                  (*it)->setPosition( downstream );
                  data.blobPos.SetXYZT( downstream.X(), downstream.Y(), downstream.Z(), time_begin );
                  data.blobBegPos.SetXYZT( downstream.X(), downstream.Y(), downstream.Z(), time_begin );
                  data.blobEndPos.SetXYZT( upstream.X(), upstream.Y(), upstream.Z(), time_begin );
                }
              else
                {
                  (*it)->setPosition( upstream );
                  data.blobPos.SetXYZT( upstream.X(), upstream.Y(), upstream.Z(), time_begin );
                  data.blobBegPos.SetXYZT( upstream.X(), upstream.Y(), upstream.Z(), time_begin );
                  data.blobEndPos.SetXYZT( downstream.X(), downstream.Y(), downstream.Z(), time_begin );
                }
              data.blobDirection = data.blobEndPos - data.blobBegPos;
              data.flightPath = data.blobBegPos - data.vtx;
            }



          debug()<<"GetBlobTrack"<<endmsg;
          m_anaBlobUtils->GetBlobTrack( data );
          info()<<"GetBlobTrackScore"<<endmsg;
          m_anaBlobUtils->GetBlobTrackScore( data );

          //Fill MC truth Info
          if( truth )
            {
              debug()<<"Start filling blob truth"<<endmsg;
              m_anaBlobUtils->FillTruth(data, *it);
              debug()<<"Finished with truth"<<endmsg;
            }
          BlobsInfo.push_back( data );

          debug()<<"Fill BlobsInfo to Data members"<<endmsg;
          //if(truth) m_anaBlobUtils->ResetTraj();
        }
      m_neutronBlobUtils->FillData( BlobsInfo, haveNeutrinoMC() ); //Fill reco related blobInfo branches and truth if haveNeutrinoMC
      debug()<<"NBlobs Before Filling: "<<AllIDContiguousBlobs.size()<<endmsg;
      debug()<<"NBlobs In Container: "<<BlobsInfo.size()<<endmsg;
      debug()<<"NBlobs After Filling: "<<endmsg;

      //m_anaBlobUtils->GetInclusiveBlobTracks( BlobsInfo );
      //---------------------------------------------------
      //! end computing neutron blobs
      //---------------------------------------------------
    }



  //--------------------------------------------------
  //! If we got to this line, set pass precuts to true
  //--------------------------------------------------
  bool pass_MasterAnaDev_precuts = true;
  event->setIntData("pass_MasterAnaDev_precuts",(int)pass_MasterAnaDev_precuts);
  event->filtertaglist()->setOrAddFilterTag("pass_MasterAnaDev_precuts", pass_MasterAnaDev_precuts );

  //---------------------------------------------------
  //! Fill common physics ana branches
  //---------------------------------------------------
  fillCommonPhysicsAnaBranches( event );
  //---------------------------------------------------
  //! Fill NuMI branches
  //---------------------------------------------------
  fillNuMIBranches(event);

  //---------------------------------------------------
  //! Fill Particle Reponse branches
  // Create a vector of all clusters associated
  // to the Non-vtx Isolated Blobs and the Non-vtx
  // dispersed blobs and use it as the input
  // for the fill function
  //---------------------------------------------------
  Minerva::IDClusterVect nonVtxIDClusters;

  //First, attach all clusters associated with the non-vtx isolated blobs
  Minerva::ProngVect isoProngs = event->select<Minerva::Prong>( "Used:Unused", "All" );
  debug() << "Event has now " << prongs.size() << " prongs" << endmsg;

  for( Minerva::ProngVect::iterator p=prongs.begin(); p!=prongs.end(); p++ ) {
    if( (*p)->creationSignature() == PrimaryBlobProngExtraDataDefs::IsolatedIDBlobProngName() ){
      Minerva::IDClusterVect clusters = (*p)->getAllIDClusters();
      for(Minerva::IDClusterVect::iterator itCluster=clusters.begin(); itCluster!=clusters.end(); itCluster++)
        {
          nonVtxIDClusters.push_back( *itCluster );
        }
    }
  }

  //Now, attach dispersed clusters
  for(Minerva::IDClusterVect::iterator itCluster=unusedIDClusters.begin(); itCluster!=unusedIDClusters.end(); itCluster++)
    nonVtxIDClusters.push_back( *itCluster );

  // Finally, Fill Particle Responses
  fillParticleResponseBranches(event,nonVtxIDClusters,"cal");

  //---------------------------------------------------
  //! New Particle response systematics
  //---------------------------------------------------
  if( haveNeutrinoMC() ) {
     info() << "Filling ParticleResponse branches" << endmsg;
     // Removing re-declaration of the muonProng object. -David Last 01/21/2022
     //SmartRef<Minerva::Prong> muonProng = (Minerva::Prong*)NULL;
     Minerva::ProngVect primaryProngs =  event->primaryProngs();
     const double mintimeCut = event->time() - 20.0*CLHEP::ns;
     const double maxtimeCut = event->time() + 35.0*CLHEP::ns;
     Minerva::IDClusterVect allProngIDClusters;
     std::vector<Minerva::IDClusterVect> container_idclusters;
     std::vector<Minerva::ODClusterVect> container_odclusters;
     int iClusters = 0;
     int idClusters = 0;
     // Removing re-declaration of the muonProng object. 
     // I have left the commented out block below until we know whether we wish to have these branches for both muon and electron events or just muon events.
     // However, the code below doesn't select muons only since in the electron case `muonProng` is set to have the `PrimaryMuon` tag set as true.
     // -David Last 01/21/2022
     /*
     for(unsigned int prong = 0; prong < primaryProngs.size(); prong++) {
       if( primaryProngs[prong]->filtertaglist()->filterTagExists("PrimaryMuon") ) {
         muonProng    = primaryProngs[prong];
         info()<< "Muon Prong Size " << muonProng->getAllIDClusters().size()<<endmsg;
         break;
       }
     }
     */
     const Minerva::IDClusterVect& muonClustersID    = muonProng->getAllIDClusters();
     for(unsigned int prong = 0; prong < primaryProngs.size(); prong++) {
       if( primaryProngs[prong]->filtertaglist()->filterTagExists("PrimaryPion") ) {
         Minerva::IDClusterVect tmp_idclusters = primaryProngs[prong]->getAllIDClusters();
         Minerva::IDClusterVect idclusters;
         Minerva::ODClusterVect odclusters;
         idClusters += tmp_idclusters.size();
         Minerva::IDClusterVect::iterator i_tmp_idclusters = tmp_idclusters.begin();
         for( ; i_tmp_idclusters != tmp_idclusters.end(); ++i_tmp_idclusters ){
           //skip if not in time window
           if( (*i_tmp_idclusters)->time() < mintimeCut ||  maxtimeCut < (*i_tmp_idclusters)->time() ) continue;
           //skip if it's a muon cluster
           if( std::find( muonClustersID.begin(), muonClustersID.end(), (*i_tmp_idclusters) ) != muonClustersID.end() ) continue;
           //blegh.  This is because the short tracker will use a cluster twice.
           // This should remove duplicate clusters in all clusters and prong clusters.
           //This is very hacky
           if( std::find( allProngIDClusters.begin(), allProngIDClusters.end(), (*i_tmp_idclusters) ) != allProngIDClusters.end() )  continue;
           //skip xtalk
           if( (*i_tmp_idclusters)->type() == 5 ) continue;
           allProngIDClusters.push_back( *i_tmp_idclusters );
           idclusters.push_back( *i_tmp_idclusters );
           ++iClusters;
         }
         container_idclusters.push_back(idclusters);
         container_odclusters.push_back(odclusters);
         tmp_idclusters.clear();
         idclusters.clear();
         odclusters.clear();
       }
     }
     info() << "Container size " << container_idclusters.size() << " " << container_odclusters.size() << endmsg;
     std::sort( allProngIDClusters.begin(), allProngIDClusters.end() );
     Minerva::IDClusterVect::iterator duplicate = std::adjacent_find( allProngIDClusters.begin(), allProngIDClusters.end() );
     if( duplicate != allProngIDClusters.end() ) info()<<"Found a cluster duplicate in part response"<<endmsg;
     Minerva::IDClusterVect hadClustersID = getHadronIDClusters( event, muonProng );
     Minerva::ODClusterVect hadClustersOD;
     if(iClusters>hadClustersID.size()) info()<<"Number of pion prong clusters > hadClusters"<<endmsg;
     info()<<"Total ID Prong clusters " <<idClusters<<endmsg;
     info()<<"Total Container clusters "<<iClusters<<endmsg;
     info()<<"Total HadronID clusters " <<hadClustersID.size()<<endmsg;
     fillModifiedParticleResponseBranches(event,container_idclusters,container_odclusters,hadClustersID,hadClustersOD);
     info() << "End filling ParticleResponse branches" << endmsg;
   }
  } // if(muonProng)

  markEvent( event );
  //-------------------------------
  //! Interpret event
  //-------------------------------
  std::vector<Minerva::NeutrinoInt*> nuInt;
  interpretEvent( event, truth, nuInt );


  //-------------------------------
  //! Mark event
  //-------------------------------
  //markEvent( event );

  debug() << "   Found " << nuInt.size() << " interpretations." << endmsg;
  if( !nuInt.empty() )
    {
      // Add control filtertags to the event
      for( vector<Minerva::NeutrinoInt*>::const_iterator it = nuInt.begin(); it != nuInt.end(); ++it )
        {
          //there should only be one NeutrinoInt.  choose the one with the best score to set the event hadronic energy.
          double hypRecoilE = 0.;
          if( (*it)->getDoubleData( AnaFilterTags::RecoilEnergy(), hypRecoilE ) )
            event->setDoubleData( AnaFilterTags::RecoilEnergy(), hypRecoilE );

          if( (*it)->filtertaglist()->isFilterTagTrue( "is_neutrino" ) )
            event->filtertaglist()->setOrAddFilterTag( "is_neutrino", true );
          else if( (*it)->filtertaglist()->isFilterTagTrue( "is_antineutrino" ) )
            event->filtertaglist()->setOrAddFilterTag( "is_antineutrino", true );

          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag() ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag(), true );
          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("Tgt1") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("Tgt1"), true );
          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("Tgt2") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("Tgt2"), true );
          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("Tgt3") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("Tgt3"), true );
          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("Tgt4") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("Tgt4"), true );
          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("Tgt5") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("Tgt5"), true );
          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("Tracker") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("Tracker"), true );

          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("C") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("C"), true );
          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("Fe") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("Fe"), true );
          if( (*it)->filtertaglist()->isFilterTagTrue( getPassTag("Pb") ) )
            event->filtertaglist()->setOrAddFilterTag( getPassTag("Pb"), true );
        }

      fillCCInclusiveBranches( event );
      fillNuMIBranches( event );
    }
  //---------------------------------------------------
  //! Add interaction hypothesis to physics event
  //---------------------------------------------------
  sc = addInteractionHyp( event, nuInt );

  info() << "Exiting MasterAnaDev::reconstructEvent() now ....." << endmsg;

  return sc;
} //End of reconstructEvent( )


//=============================================================================
// interpretEvent() --
//=============================================================================
StatusCode MasterAnaDev::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction* truth, std::vector<Minerva::NeutrinoInt*>& nuInt ) const {

  info() << "Entering MasterAnaDev::interpretEvent() now ....." << endmsg;
  if( truth ){
    debug() << "This event has a matched MC interaction" << endmsg;
  }

  //! Construct neutrino hypotheses
  Minerva::NeutrinoInt* ccqeHyp = new Minerva::NeutrinoInt( "MasterAnaDev" );
  nuInt.push_back( ccqeHyp );
  ccqeHyp->setNeutrinoFlavor( Minerva::NeutrinoInt::MuonFlavor );
  ccqeHyp->setInteractionCurrent( Minerva::NeutrinoInt::ChargedCurrent );
  ccqeHyp->setInteractionType( Minerva::NeutrinoInt::QuasiElastic );

  if(!event->hasInteractionVertex()){
    //if(!has_vtx){  //! Construct neutrino hypothesis
    //! Construct neutrino hypothesis
    return StatusCode::SUCCESS; // have to do this for events without interaction vertices because basically all the objects below are not valid. The code will crash and burn.
  }

  SmartRef<Minerva::Prong> muonProng;
  SmartRef<Minerva::Particle> muonPart;
  ProngVect hadronProngs;
  bool has_muon = MuonUtils->findMuonProng( event, muonProng, muonPart);
  bool has_electron = !has_muon && m_electronUtils->findElectronProng(event,muonProng,muonPart);
  bool has_lepton = true;
  if (!has_muon && !has_electron) has_lepton = false; // CODE WAS return StatusCode::SUCCESS;
  // MUST PROTECT BELOW FOR !has_lepton
  SmartRef<Minerva::Particle> neutronPart;
  if (has_lepton) neutronPart = m_neutronBlobUtils->computeNeutronFromHydrogen( muonPart );
  if (!has_lepton) {
    ccqeHyp->setInteractionCurrent( Minerva::NeutrinoInt::NeutralCurrent );
    ccqeHyp->setInteractionType( Minerva::NeutrinoInt::Elastic );
  }

  //check for missing in_fiducial_area branch
  SmartRef<Minerva::Vertex> vtx = event->interactionVertex();
  //Gaudi::XYZPoint vtxPos = vtx->position();

  //! Do the fiducial cut on the projected point (XY only this time)
  bool isFiducialXY   = FiducialPointTool->isFiducial(vtx, m_signalApothem, m_usFiducialModID, m_dsFiducialModID);
  bool isAnalyzableXY = FiducialPointTool->isFiducial(vtx, m_analyzeApothem,m_usAnalyzableModID,m_dsAnalyzableModID);


  if( ! isAnalyzableXY )
    {
      debug() << "    Projected vertex is not in the analyzable volume." << endmsg;
      // if( ! force )
      return StatusCode::SUCCESS;
    }


  int muon_charge = 0;
  if ( has_muon && MuonUtils->muonCharge( muonProng, muon_charge, m_qOverpChargeCut ) == StatusCode::FAILURE )
    {
      debug() << "Could not get a muon charge.  Mark the event with a low score and do your best." << endmsg;
    } else if (has_electron && event->hasDoubleData(NuMIVar::horn_curr))
    {
      muon_charge = event->getDoubleData(NuMIVar::horn_curr) >0 ? 1:-1;
    }


  double muon_P,muon_Px, muon_Py, muon_Pz, muon_E, muon_T, muon_theta, muon_score = -1.0, muon_theta_biasUp, muon_theta_biasDown, muon_qp = -9999.9, muon_qpqpe = -9999.9;

  if ( has_lepton ) {
    int muon_minosTrackQuality = 0;
    //!Get muon momentum, energy and angle
    Gaudi::LorentzVector muon_4p = muonPart->momentumVec();
    muon_P = muonPart->momentumVec().P();
    muon_Px = muonPart->momentumVec().px();
    muon_Py = muonPart->momentumVec().py();
    muon_Pz = muonPart->momentumVec().pz();
    muon_E = muonPart->momentumVec().E();
    muon_T = muon_E - MinervaUnits::M_mu;
    muon_theta = MuonUtils->muTheta( muonPart );
    muon_score = muonPart->score();
    muon_theta_biasUp   = m_minCoordSysTool->thetaWRTBeam(muon_4p,m_beamAngleBias) - muon_theta;
    muon_theta_biasDown = m_minCoordSysTool->thetaWRTBeam(muon_4p, -1.0*m_beamAngleBias) - muon_theta;
    muon_qpqpe = MuonUtils->minosQPQPE(muonProng);
    if (muonProng->MinosTrack() && !muonProng->minosTracks().empty()) {
      SmartRef<Minerva::MinosRecoTrack> minosTrack = muonProng->minosTracks()[0];
      muon_qp = minosTrack->qp();
      muon_minosTrackQuality = minosTrack->quality();
    }
    //--------------------------------------------------------------
    //! Hadronic Recoil
    //--------------------------------------------------------------
    double hadron_recoilE_CCInc = RecoilUtils->calcRecoilEFromClusters( event, muonProng, "CCInclusive", -1);
    double hadron_recoilE = RecoilUtils->calcRecoilEFromClusters( event, muonProng, "CCNuPionInc");
    double hadron_recoilE_default = RecoilUtils->calcRecoilEFromClusters( event, muonProng, "Default" );
    double hadron_recoilE_two_track = RecoilUtils->calcRecoilEFromClusters( event, muonProng, "NukeCCPion_TwoTrack_Nu_Tracker");

    ccqeHyp->setDoubleData("hadron_recoil",               hadron_recoilE);
    ccqeHyp->setDoubleData("hadron_recoil_default",       hadron_recoilE_default);
    ccqeHyp->setDoubleData("hadron_recoil_two_track",     hadron_recoilE_two_track);
    ccqeHyp->setDoubleData("hadron_recoil_CCInc",         hadron_recoilE_CCInc);
  }

  //--------------------------------------------------------------
  //! Get MCTracks vector
  //--------------------------------------------------------------
  std::vector<MCTrack> MCTracksVector;
    if (truth)
      MCTracksVector = m_MCTrackTool->getMCTracks(truth);

  //-----------------------------------------------
  //! Get proton kinematics stuff from dEdx info
  //-----------------------------------------------
  double proton_P = -9999., proton_E = -9999., proton_T = -9999., proton_theta = -9999., proton_Px = -9999.,  proton_Py = -9999.,  proton_Pz = -9999.;
  ProngVect primaryProngs = event->primaryProngs();
  SmartRef<Minerva::Prong> hadronProng;
  SmartRef<Minerva::Prong> pionProng;
  SmartRef<Minerva::Prong> protonProng;
  SmartRef<Minerva::Particle> protonPart;
  SmartRef<Minerva::Particle> pionPart;
  double totalCalibEPrimaryProton   = -9999.;

  ProngVect secondaryProtonProngs;
  std::vector<int> sec_protons_prongType; // 0: ok, 1: kinked, 2: forked
  std::vector<double> sec_protons_E;
  std::vector<double> sec_protons_P;
  std::vector<double> sec_protons_T;
  std::vector<double> sec_protons_Px;
  std::vector<double> sec_protons_Py;
  std::vector<double> sec_protons_Pz;
  std::vector<double> sec_protons_T_Calo;
  std::vector<double> sec_protons_theta;
  std::vector<double> sec_protons_proton_scores;
  std::vector<double> sec_protons_proton_scores1;
  std::vector<double> sec_protons_proton_scores2;
  std::vector<double> sec_protons_pion_scores;
  std::vector<double> sec_protons_pion_scores1;
  std::vector<double> sec_protons_pion_scores2;
  std::vector<int> sec_protons_nodes_index;
  std::vector<double> sec_protons_nodes_nodesNormE;
  std::vector<double> sec_protons_nodes_E;

  //------------------------------------------------------------------------
  //! dEdx Tool variations parameters for Primary and Secondary Protons
  //------------------------------------------------------------------------
  double proton_score1_Mass_Up=-9999.9, proton_score1_Mass_Down=-9999.9, proton_score1_BetheBloch_Up=-9999.9, proton_score1_BetheBloch_Down=-9999.9;
  double proton_score1_MEU_Up=-9999.9, proton_score1_MEU_Down=-9999.9, proton_score1_Birks=-9999.9;
  double proton_E_Mass_Up=-9999.9, proton_E_Mass_Down=-9999.9, proton_E_BetheBloch_Up=-9999.9, proton_E_BetheBloch_Down=-9999.9;
  double proton_E_MEU_Up=-9999.9, proton_E_MEU_Down=-9999.9, proton_E_Birks=-9999.9;
  int proton_patternRec = -1;

  std::vector<double> sec_protons_score1_Mass_Up, sec_protons_score1_Mass_Down, sec_protons_score1_BetheBloch_Up, sec_protons_score1_BetheBloch_Down;
  std::vector<double> sec_protons_score1_MEU_Up, sec_protons_score1_MEU_Down, sec_protons_score1_Birks;
  std::vector<double> sec_protons_E_Mass_Up, sec_protons_E_Mass_Down, sec_protons_E_BetheBloch_Up, sec_protons_E_BetheBloch_Down;
  std::vector<double> sec_protons_E_MEU_Up, sec_protons_E_MEU_Down, sec_protons_E_Birks;
  std::vector<int> sec_protons_patternRec;

  debug()<< "Getting most energetic (Primary) proton and Secondary proton prongs" << endmsg;
  bool protonFound = false;
  bool secondaryProtonsFound = false;

  for (ProngVect::iterator itProng = primaryProngs.begin(); itProng != primaryProngs.end(); ++itProng) {

      bool isPrimaryMuon = false, isPrimaryHadron = false;

      if ( (*itProng)->filtertaglist()->filterTagExists("PrimaryMuon") ) {
        (*itProng)->filtertaglist()->checkFilterTag( "PrimaryMuon", isPrimaryMuon );
      }

      if ( (*itProng)->filtertaglist()->filterTagExists("PrimaryHadron") ) {
        (*itProng)->filtertaglist()->checkFilterTag( "PrimaryHadron", isPrimaryHadron );
      }

     // if (isPrimaryMuon && !isPrimaryHadron) muonProng = *itProng;
      if (isPrimaryHadron && !isPrimaryMuon) hadronProngs.push_back(*itProng);
      if (isPrimaryMuon && isPrimaryHadron ) {
        warning()<<"Prong is two primary particles!"<<endmsg;
     }
  }
     
  for (ProngVect::iterator itProngs = primaryProngs.begin(); itProngs != primaryProngs.end(); ++itProngs) {
    if( (*itProngs)->filtertaglist()->filterTagExists( m_primaryProton ) ) {
      protonProng = (*itProngs);


      totalCalibEPrimaryProton = m_caloUtils->applyCalConsts( protonProng, "Default", false, true );
      debug() << "Total calib energy of primary proton " << totalCalibEPrimaryProton << endmsg;

      SmartRefVector<Minerva::Particle> pionparticles = (*itProngs)->particlesWithIDCode(Minerva::Particle::Pion, "dEdX");
      if (pionparticles.size()>0) pionPart = pionparticles[0];
      protonPart  = (*itProngs)->bestParticle();
      if (protonPart) protonFound = true;
    }
    else if( (*itProngs)->filtertaglist()->filterTagExists( m_secondaryProtons ) ) {
      secondaryProtonProngs.push_back( *itProngs );
      secondaryProtonsFound = true;
    }
  }

  double proton_score = -9999., proton_score1 = -9999., proton_score2 = -9999.;
  double pion_score = -9999., pion_score1 = -9999., pion_score2 = -9999.;

  std::vector<double> protonE;
  double proton_startPointX = -9999., proton_startPointY = -9999., proton_startPointZ = -9999.;
  double proton_endPointX = -9999., proton_endPointY = -9999., proton_endPointZ = -9999.;
  double proton_thetaX = -9999., proton_thetaY = -9999., proton_phi = -9999.;
  int proton_prongType = -9999;

  std::vector<double> proton_nodes_nodesNormE;
  std::vector<double> proton_nodes_E;
  if (protonFound) {
    debug() << "Proton prongs gotten"<<endmsg;

    proton_score  = protonPart->score();
    proton_score1 = protonPart->getDoubleData(ParticleExtraDataDefs::dEdXScore1());
    proton_score2 = protonPart->getDoubleData(ParticleExtraDataDefs::dEdXScore2());

    // Set proton prong type, kinked=1, forked=2, ok = 0
    if ( protonProng->Kinked() ) proton_prongType = 1;
    else if ( protonProng->Forked() ) proton_prongType = 2;
    else proton_prongType = 0;

    debug() << "Proton range score1: "<< proton_score1 << endmsg;

    proton_P = protonPart->momentumVec().P();
    proton_E = protonPart->momentumVec().E();
    proton_T = proton_E - MinervaUnits::M_p;
    proton_Px= protonPart->momentumVec().px();
    proton_Py= protonPart->momentumVec().py();
    proton_Pz= protonPart->momentumVec().pz();
    proton_patternRec = protonProng->minervaTracks()[0]->patRecHistory();


    //-----------------------------------------------------------
    //-----------------------------------------------------------
    //double theta = protonProng->minervaTracks().front()->theta();

    proton_thetaX  = m_minCoordSysTool->thetaXWRTBeam(protonPart->momentumVec());
    proton_thetaY  = m_minCoordSysTool->thetaYWRTBeam(protonPart->momentumVec());
    proton_phi     = m_minCoordSysTool->phiWRTBeam(protonPart->momentumVec());

    proton_endPointX   = (protonProng->minervaTracks().back())->lastState().x();
    proton_endPointY   = (protonProng->minervaTracks().back())->lastState().y();
    proton_endPointZ   = (protonProng->minervaTracks().back())->lastState().z();

    proton_startPointX = (protonProng->minervaTracks().back())->firstState().x();
    proton_startPointY = (protonProng->minervaTracks().back())->firstState().y();
    proton_startPointZ = (protonProng->minervaTracks().back())->firstState().z();

    proton_nodes_nodesNormE = m_masterAnaDevRecoUtils->getNormalizedLastNodesEnergy( protonProng->minervaTracks().back(), true );
    proton_nodes_E = m_masterAnaDevRecoUtils->getNormalizedLastNodesEnergy( protonProng->minervaTracks().back(), false );

    //------------------------------------------------------
    // dEdx Tool variations for Primary Proton
    //------------------------------------------------------
    if( protonPart->hasDoubleData("score1_Mass_Up") ) {
      proton_score1_Mass_Up = protonPart->getDoubleData("score1_Mass_Up") - proton_score1;
      proton_E_Mass_Up      = protonPart->getDoubleData("FitE_Mass_Up") - proton_E;
    }
    if( protonPart->hasDoubleData("score1_Mass_Down") ) {
      proton_score1_Mass_Down = protonPart->getDoubleData("score1_Mass_Down") - proton_score1;
      proton_E_Mass_Down      = protonPart->getDoubleData("FitE_Mass_Down") - proton_E;
    }
    if( protonPart->hasDoubleData("score1_BetheBloch_Up") ) {
      proton_score1_BetheBloch_Up = protonPart->getDoubleData("score1_BetheBloch_Up") - proton_score1;
      proton_E_BetheBloch_Up      = protonPart->getDoubleData("FitE_BetheBloch_Up") - proton_E;
    }
    if( protonPart->hasDoubleData("score1_BetheBloch_Down") ) {
      proton_score1_BetheBloch_Down = protonPart->getDoubleData("score1_BetheBloch_Down") - proton_score1;
      proton_E_BetheBloch_Down      = protonPart->getDoubleData("FitE_BetheBloch_Down") - proton_E;
    }
    if( protonPart->hasDoubleData("score1_MEU_Up") ) {
      proton_score1_MEU_Up = protonPart->getDoubleData("score1_MEU_Up") - proton_score1;
      proton_E_MEU_Up      = protonPart->getDoubleData("FitE_MEU_Up") - proton_E;
    }
    if( protonPart->hasDoubleData("score1_MEU_Down") ) {
      proton_score1_MEU_Down = protonPart->getDoubleData("score1_MEU_Down") - proton_score1;
      proton_E_MEU_Down      = protonPart->getDoubleData("FitE_MEU_Down") - proton_E;
    }
    if( protonPart->hasDoubleData("score1_Birks") ) {
      proton_score1_Birks  = protonPart->getDoubleData("score1_Birks") - proton_score1;
      proton_E_Birks       = protonPart->getDoubleData("FitE_Birks") - proton_E;
    }

    //----------------------------------------
    // Rotated proton's angle w.r.t beam
    //----------------------------------------
    proton_theta  = m_minCoordSysTool->thetaWRTBeam(protonPart->momentumVec());

    //double proton_thetaX = m_minCoordSysTool->thetaXWRTBeam(protonPart->momentumVec());
    //double proton_thetaY = m_minCoordSysTool->thetaYWRTBeam(protonPart->momentumVec());

    debug()<< "Calculated from dEdx: proton_E, proton_T, proton_P, proton_theta"<<endmsg;

  } else debug() << " Could not get Primary proton prong" << endmsg;

  if (secondaryProtonsFound) {
    int index = 0;
    for(ProngVect::iterator secprongs = secondaryProtonProngs.begin(); secprongs != secondaryProtonProngs.end(); ++secprongs) {
      double secproton_proton_score = -9999., secproton_proton_score1 = -9999., secproton_proton_score2 = -999.;
      double totalCalibESecondaryProton = -9999.;

      totalCalibESecondaryProton = m_caloUtils->applyCalConsts( (*secprongs), "Default", false, true );
      debug() << "Total calib energy of secondary proton " << totalCalibESecondaryProton << endmsg;

      SmartRef<Minerva::Particle> secprongPart = (*secprongs)->bestParticle();

      secproton_proton_score  = secprongPart->score();
      secproton_proton_score1 = secprongPart->getDoubleData(ParticleExtraDataDefs::dEdXScore1());
      secproton_proton_score2 = secprongPart->getDoubleData(ParticleExtraDataDefs::dEdXScore2());

      sec_protons_proton_scores.push_back( secproton_proton_score );
      sec_protons_proton_scores1.push_back( secproton_proton_score1 );
      sec_protons_proton_scores2.push_back( secproton_proton_score2 );

      sec_protons_E.push_back( secprongPart->momentumVec().E() );
      sec_protons_P.push_back( secprongPart->momentumVec().P() );
      sec_protons_T.push_back( secprongPart->momentumVec().E() - MinervaUnits::M_p );
      sec_protons_Px.push_back( secprongPart->momentumVec().px() );
      sec_protons_Py.push_back( secprongPart->momentumVec().py() );
      sec_protons_Pz.push_back( secprongPart->momentumVec().pz() );
      sec_protons_T_Calo.push_back( totalCalibESecondaryProton );

      sec_protons_patternRec.push_back(protonProng->minervaTracks()[0]->patRecHistory());
      //Add px,py,pz here
      int prongtype = -9999;
      if ((*secprongs)->Kinked() ) prongtype = 1;
      else if ((*secprongs)->Forked() ) prongtype = 2;
      else prongtype = 0;

      sec_protons_prongType.push_back( prongtype );


      //------------------------------------------------
      // dEdx Tool variations for Secondary Protons
      //------------------------------------------------
      if( secprongPart->hasDoubleData("score1_Mass_Up") ) {
        sec_protons_score1_Mass_Up.push_back( secprongPart->getDoubleData("score1_Mass_Up") - proton_score1 );
        sec_protons_E_Mass_Up.push_back( secprongPart->getDoubleData("FitE_Mass_Up") - proton_E );
      }
      if( secprongPart->hasDoubleData("score1_Mass_Down") ) {
        sec_protons_score1_Mass_Down.push_back( secprongPart->getDoubleData("score1_Mass_Down") - proton_score1 );
        sec_protons_E_Mass_Down.push_back( secprongPart->getDoubleData("FitE_Mass_Down") - proton_E );
      }
      if( secprongPart->hasDoubleData("score1_BetheBloch_Up") ) {
        sec_protons_score1_BetheBloch_Up.push_back( secprongPart->getDoubleData("score1_BetheBloch_Up") - proton_score1 );
        sec_protons_E_BetheBloch_Up.push_back( secprongPart->getDoubleData("FitE_BetheBloch_Up") - proton_E );

        debug() << "Proton Score1 BB Up Variation " << secprongPart->getDoubleData("score1_BetheBloch_Up") << endmsg;
      }
      if( secprongPart->hasDoubleData("score1_BetheBloch_Down") ) {
        sec_protons_score1_BetheBloch_Down.push_back( secprongPart->getDoubleData("score1_BetheBloch_Down") - proton_score1 );
        sec_protons_E_BetheBloch_Down.push_back( secprongPart->getDoubleData("FitE_BetheBloch_Down") - proton_E );

        debug() << "Proton Score1 BB Down Variation " << secprongPart->getDoubleData("score1_BetheBloch_Down") << endmsg;
      }
      if( secprongPart->hasDoubleData("score1_MEU_Up") ) {
        sec_protons_score1_MEU_Up.push_back( secprongPart->getDoubleData("score1_MEU_Up") - proton_score1 );
        sec_protons_E_MEU_Up.push_back( secprongPart->getDoubleData("FitE_MEU_Up") - proton_E );
      }
      if( secprongPart->hasDoubleData("score1_MEU_Down") ) {
        sec_protons_score1_MEU_Down.push_back( secprongPart->getDoubleData("score1_MEU_Down") - proton_score1 );
        sec_protons_E_MEU_Down.push_back( secprongPart->getDoubleData("FitE_MEU_Down") - proton_E );
      }
      if( secprongPart->hasDoubleData("score1_Birks") ) {
        sec_protons_score1_Birks.push_back( secprongPart->getDoubleData("score1_Birks") - proton_score1 );
        sec_protons_E_Birks.push_back( secprongPart->getDoubleData("FitE_Birks") - proton_E );
      }

      //----------------------------------------
      // ESC proton variables
      //----------------------------------------
      std::vector<double> thisNodesNormE = m_masterAnaDevRecoUtils->getNormalizedLastNodesEnergy( (*secprongs)->minervaTracks().back(), true );
      std::vector<double>  thisNodesE=m_masterAnaDevRecoUtils->getNormalizedLastNodesEnergy( (*secprongs)->minervaTracks().back(), false );

      sec_protons_nodes_index.insert( sec_protons_nodes_index.end(), thisNodesE.size(), index );
      sec_protons_nodes_nodesNormE.insert(sec_protons_nodes_nodesNormE.end(), thisNodesNormE.begin(), thisNodesNormE.end() );
      sec_protons_nodes_E.insert( sec_protons_nodes_E.end(), thisNodesE.begin(), thisNodesE.end() );
      //----------------------------------------
      // Rotated proton's angle w.r.t beam
      //----------------------------------------
      sec_protons_theta.push_back( m_minCoordSysTool->thetaWRTBeam(secprongPart->momentumVec()) );
      index++;

    } //End of for loop over Secondary proton prongs

  } else verbose() << "No secondary protons were found" << endmsg;

  //--------------------------------------------------------
  //! Set pion scores
  //--------------------------------------------------------
  //Primary prong pion particle

  std::vector<double> hadron_piFit_scoreLLR;
  std::vector<int> hadron_isTracker, hadron_isSideECAL;
  std::vector<int> hadron_isExiting, hadron_isODMatch, hadron_isForked;
  std::vector<int> hadron_endMichel_category, hadron_endMichel_ndigits;
  std::vector<double> hadron_endMichel_slice_energy, hadron_endMichel_energy;
  std::vector<int> hadron_tm_PDGCode, hadron_tm_trackID, hadron_tm_destructCode, hadron_tm_isTruthMatched;
  std::vector<double> hadron_tm_beginKE, hadron_tm_fraction, hadron_tm_otherE;

  std::vector<double> pion_E;
  std::vector<double> pion_T;
  std::vector<double> pion_P;
  std::vector<double> pion_Px;
  std::vector<double> pion_Py;
  std::vector<double> pion_Pz;
  std::vector<double> pion_theta;
  std::vector<double> pion_phi;

  // Passive-corrected pion kinematics
  std::vector<double> hadron_pion_E_corr,  hadron_pion_p_corr;
  std::vector<double> hadron_pion_px_corr, hadron_pion_py_corr;
  std::vector<double> hadron_pion_pz_corr, hadron_pion_E_recoil_corr;

  // Xiango's Node cut
  std::vector<int>    pion_nNodes;
  std::vector<double> pion_lastnode_Q0;
  std::vector<double> pion_lastnode_Q1;
  std::vector<double> pion_lastnode_Q2;
  std::vector<double> pion_lastnode_Q3;
  std::vector<double> pion_lastnode_Q4;
  std::vector<double> pion_lastnode_Q5;

  //dEdXTool biasing
  std::vector<double> hadron_piFit_score1_Mass_Up, hadron_piFit_score1_Mass_Down, hadron_piFit_score1_BetheBloch_Up, hadron_piFit_score1_BetheBloch_Down, hadron_piFit_score1;
  std::vector<double> hadron_piFit_score1_MEU_Up, hadron_piFit_score1_MEU_Down, hadron_piFit_score1_Birks;
  std::vector<double> hadron_pion_E_Mass_Up, hadron_pion_E_Mass_Down, hadron_pion_E_BetheBloch_Up, hadron_pion_E_BetheBloch_Down;
  std::vector<double> hadron_pion_E_Birks;
  std::vector<double> hadron_pion_theta_biasUp, hadron_pion_theta_biasDown;

  int iter;
  for (iter = 0; iter < 10; iter++){
    pion_E.push_back(-1.);pion_T.push_back(-9999.);pion_P.push_back(-1.);pion_Px.push_back(-1.);pion_Py.push_back(-1.);pion_Pz.push_back(-1.);pion_theta.push_back(-9.);pion_phi.push_back(-9.);
    hadron_piFit_scoreLLR.push_back(-101.0);
    hadron_isTracker.push_back(-1); hadron_isSideECAL.push_back(-1);
    hadron_isExiting.push_back(-1); hadron_isODMatch.push_back(-1); hadron_isForked.push_back(-1);
    hadron_endMichel_category.push_back(0); hadron_endMichel_slice_energy.push_back(-1.0);hadron_endMichel_ndigits.push_back(-1);hadron_endMichel_energy.push_back(-1.0);
    hadron_pion_theta_biasUp.push_back(-9.0);hadron_pion_theta_biasDown.push_back(-9.0);
    hadron_piFit_score1_Mass_Up.push_back(-999.9); hadron_piFit_score1_Mass_Down.push_back(-999.9); hadron_piFit_score1_BetheBloch_Up.push_back(-999.9);
    hadron_piFit_score1_BetheBloch_Down.push_back(-999.9); hadron_piFit_score1_MEU_Up.push_back(-999.9); hadron_piFit_score1_MEU_Down.push_back(-999.9);
    hadron_piFit_score1_Birks.push_back(-999.9); hadron_piFit_score1.push_back(-1.0);
    hadron_pion_E_Mass_Up.push_back(-999.9); hadron_pion_E_Mass_Down.push_back(-999.9); hadron_pion_E_BetheBloch_Up.push_back(-999.9); hadron_pion_E_BetheBloch_Down.push_back(-999.9);
    hadron_pion_E_Birks.push_back(-999.9);
    hadron_tm_PDGCode.push_back(0);hadron_tm_destructCode.push_back(-3);hadron_tm_beginKE.push_back(-1.0);hadron_tm_trackID.push_back(-1);hadron_tm_isTruthMatched.push_back(-1); hadron_tm_fraction.push_back(-1.0); hadron_tm_otherE.push_back(-1.0);

    //Xianguo node cut.
    pion_nNodes.push_back(-1);
    pion_lastnode_Q0.push_back(-1.);
    pion_lastnode_Q1.push_back(-1.);
    pion_lastnode_Q2.push_back(-1.);
    pion_lastnode_Q3.push_back(-1.);
    pion_lastnode_Q4.push_back(-1.);
    pion_lastnode_Q5.push_back(-1.);

    // Passive-corrected pion kinematics
    hadron_pion_E_corr.push_back(-1.0); hadron_pion_p_corr.push_back(-1.0);
    hadron_pion_px_corr.push_back(-1.0); hadron_pion_py_corr.push_back(-1.0);
    hadron_pion_pz_corr.push_back(-1.0); hadron_pion_E_recoil_corr.push_back(-1.0);
  }

  int hadron_number = hadronProngs.size();
  if (hadron_number > 10) hadron_number = 10;
  iter=0;
  for (ProngVect::iterator itProngs = hadronProngs.begin(); itProngs != hadronProngs.end() && iter <10; ++itProngs, ++iter) {
    debug()<<"iter: "<<iter<<endmsg;
    std::cout << "iter = " << iter << "\n";
    SmartRef<Minerva::Particle> pionPart1 = (Minerva::Particle*)NULL;
    SmartRef<Minerva::Prong> hadronProng = *itProngs;
    SmartRefVector<Minerva::Track> hadronTracks = hadronProng->minervaTracks();
    if (!hadronProng) {
        error()<<"Should not be missing this prong!"<<endmsg;
        continue;
    }
    getBestParticle(hadronProng, pionPart1, Minerva::Particle::Pion);
    getBestParticle(hadronProng, protonPart, Minerva::Particle::Proton);
    hadron_isTracker[iter]   = hadronProng->Tracker();
    hadron_isSideECAL[iter]  = hadronProng->SideEcal();
    hadron_isExiting[iter]   = hadronProng->ExitID();
    hadron_isODMatch[iter]   = hadronProng->OdMatch();
    hadron_isForked[iter]    = hadronProng->Forked();

    if (hadronProng->hasIntData("hasEndpointMichel")){
      hadron_endMichel_category[iter] = hadronProng->getIntData("endMichel_category");
      hadron_endMichel_slice_energy[iter] = hadronProng->getDoubleData("endMichel_slice_energy");
      hadron_endMichel_ndigits[iter] = hadronProng->getIntData("endMichel_ndigits");
      hadron_endMichel_energy[iter] = hadronProng->getDoubleData("endMichel_energy");
    }

    if ( pionPart1 ) {

      Gaudi::LorentzVector pion_4p = pionPart1->momentumVec();
      double p_mag = sqrt( pow(pion_4p.px(),2.0) + pow(pion_4p.py(),2.0) + pow(pion_4p.pz(),2.0) );
      pion_E[iter]  = pion_4p.E();
      pion_T[iter]  = pion_E[iter]-MinervaUnits::M_pion;
      pion_P[iter]  = p_mag;
      pion_Px[iter] = pion_4p.px();
      pion_Py[iter] = pion_4p.py();
      pion_Pz[iter] = pion_4p.pz();
      Gaudi::LorentzVector tmp_pion_4p(pion_Px[iter], pion_Py[iter], pion_Pz[iter], pion_E[iter]);
      pion_theta[iter] = m_minCoordSysTool->thetaWRTBeam(tmp_pion_4p);
      pion_phi[iter] = m_minCoordSysTool->phiWRTBeam(tmp_pion_4p);
      hadron_pion_theta_biasUp[iter] = m_minCoordSysTool->thetaWRTBeam(pion_4p, m_beamAngleBias) - pion_theta[iter];
      hadron_pion_theta_biasDown[iter] = m_minCoordSysTool->thetaWRTBeam(pion_4p, -1.0*m_beamAngleBias) - pion_theta[iter];

       // Kinematics: passive material-corrected
       Gaudi::LorentzVector pion_4p_corr;

       double vtx_z = event->interactionVertex()->position().z();
       StatusCode sc = m_energyCorrectionTool->getCorrectedEnergy(hadronProng, pionPart1, vtx_z, pion_4p_corr);
       if(!sc) pion_4p_corr = pion_4p;

            hadron_pion_E_recoil_corr[iter] =
                m_caloUtils->applyCalConsts(hadronProng,"Default",true,true);

            double e_cal_correction = pion_4p_corr.E();
            double p_cal_correction = sqrt(pow(e_cal_correction,2) - pow(pionPart1->mass(),2));

            hadron_pion_E_corr[iter]  = sqrt(p_cal_correction*p_cal_correction + MinervaUnits::M_pion*MinervaUnits::M_pion);
            hadron_pion_px_corr[iter] = p_cal_correction*sin(pion_theta[iter])*cos(pion_phi[iter]);
            hadron_pion_py_corr[iter] = p_cal_correction*sin(pion_theta[iter])*sin(pion_phi[iter]);
            hadron_pion_pz_corr[iter] = p_cal_correction*cos(pion_theta[iter]);
            hadron_pion_p_corr[iter]  = sqrt(hadron_pion_px_corr[iter]*hadron_pion_px_corr[iter] + hadron_pion_py_corr[iter]*hadron_pion_py_corr[iter] + hadron_pion_pz_corr[iter]*hadron_pion_pz_corr[iter]);

      //score1
      if ( pionPart1->hasDoubleData("score1") ) {
        double sc1 = pionPart1->getDoubleData(ParticleExtraDataDefs::dEdXScore1());
        double Epi = pion_E[iter];

        hadron_piFit_score1[iter]   = sc1;

        //dEdXTool biasing
        if (pionPart1->hasDoubleData("score1_Mass_Up")) {
          hadron_piFit_score1_Mass_Up[iter] = pionPart1->getDoubleData("score1_Mass_Up") - sc1;
          hadron_pion_E_Mass_Up[iter]       = pionPart1->getDoubleData("FitE_Mass_Up") - Epi;
            std::cout <<"Pass hadron pion Mass Up " << hadron_pion_E_Mass_Up[iter]  << "\n" ;
        }

        if (pionPart1->hasDoubleData("score1_Mass_Down")) {
          hadron_piFit_score1_Mass_Down[iter] = pionPart1->getDoubleData("score1_Mass_Down") - sc1;
          hadron_pion_E_Mass_Down[iter]       = pionPart1->getDoubleData("FitE_Mass_Down") - Epi;
          std::cout <<"Pass hadron pion Mass Down " << hadron_pion_E_Mass_Down[iter]  << "\n" ;
        }

        if (pionPart1->hasDoubleData("score1_Birks")) {
          hadron_piFit_score1_Birks[iter] = pionPart1->getDoubleData("score1_Birks") - sc1;
          hadron_pion_E_Birks[iter]       = pionPart1->getDoubleData("FitE_Birks") - Epi;
          std::cout <<"Pass hadron pion E Birks " << hadron_pion_E_Birks[iter]  << "\n" ;
        }
        if (pionPart1->hasDoubleData("score1_BetheBloch_Up")) {
          hadron_piFit_score1_BetheBloch_Up[iter] = pionPart1->getDoubleData("score1_BetheBloch_Up") - sc1;
          hadron_pion_E_BetheBloch_Up[iter]       = pionPart1->getDoubleData("FitE_BetheBloch_Up") - Epi;
          std::cout <<"Pass hadron pion Mass Up " << hadron_pion_E_BetheBloch_Up[iter]  << "\n" ;
        }
        if (pionPart1->hasDoubleData("score1_BetheBloch_Down")) {
          hadron_piFit_score1_BetheBloch_Down[iter] = pionPart1->getDoubleData("score1_BetheBloch_Down") - sc1;
          hadron_pion_E_BetheBloch_Down[iter]       = pionPart1->getDoubleData("FitE_BetheBloch_Down") - Epi;
        std::cout <<"Pass hadron pion Mass Up " << hadron_pion_E_BetheBloch_Down[iter]  << "\n" ;      
        }
      }
      else {
        hadron_piFit_score1[iter]      = 0.0;
        warning()<<"No score1 data available!"<<endmsg;
      }

      //LLR Tool
      std::vector<Minerva::Particle::ID> particleHypotheses;
      particleHypotheses.push_back( Minerva::Particle::Pion );
      Minerva::ParticleVect particlesLikelihood;
      m_LikelihoodPIDTool->makeParticles( hadronProng , particlesLikelihood, particleHypotheses );
      if(!particlesLikelihood.empty()){
        hadron_piFit_scoreLLR[iter] = particlesLikelihood[0]->score();
      }
      else {
        warning()<<"No LLR score available"<<endmsg;
        std::cout << "No LLR = " << hadron_piFit_scoreLLR[iter] << "\t" << iter << "\n";
      }
      if (pionPart1->hasDoubleData("score1_Mass_Up")) hadron_piFit_score1_Mass_Up[iter] = pionPart1->getDoubleData("score1_Mass_Up") - hadron_piFit_score1[iter];
      }
      // Xianguo's node cut variables
      Minerva::Track *tmppiontrack = hadronTracks[0];
      const int nNodes = tmppiontrack->nNodes();
      const int nMax = 6;
      double pEnergies[nMax];
      for(int jj=0; jj<nMax; jj++){
        pEnergies[jj] = -999;

        const int iread = nNodes - 1 - jj;
        if(iread>=0){
          Minerva::Node * node = tmppiontrack->nodes()[iread];
          //http://nusoft.fnal.gov/minerva/minervadat/software_doxygen/HEAD/MINERVA/LikelihoodPIDTool_8cpp-source.html
          pEnergies[jj]=  node->idcluster()->energy() / sqrt( 1 + pow( node->state().ax(),2 ) + pow( node->state().ay(),2 ) ) ;
        }
      }
      pion_nNodes[iter]= nNodes;
      pion_lastnode_Q0[iter] = pEnergies[0];
      pion_lastnode_Q1[iter] = pEnergies[1];
      pion_lastnode_Q2[iter] = pEnergies[2];
      pion_lastnode_Q3[iter] = pEnergies[3];
      pion_lastnode_Q4[iter] = pEnergies[4];
      pion_lastnode_Q5[iter] = pEnergies[5];
      // End xiango's stuff

      //=============================================================
      //! Truth matching (MC ONLY beyond this point)
      //=============================================================
      //const Minerva::TG4Trajectory* endpt_trajectory = NULL;
      const Minerva::TG4Trajectory* hadron_trajectory = NULL;
      double fraction = -1.0;
      double other_energy = -1.0;
      if (!haveNeutrinoMC()) continue;

      //===========================================================
      // Truth match the whole prong
      //===========================================================
      StatusCode sc = TruthMatcher->getPrimaryTG4Trajectory( hadronProng, hadron_trajectory, fraction, other_energy, false);

      if (sc == StatusCode::FAILURE) {
        warning()<<"Truth Matching failed"<<endmsg;
        hadron_tm_isTruthMatched[iter] = 0;
        hadron_tm_fraction[iter] = -1.0;
        hadron_tm_otherE[iter] = other_energy;
        continue;
      }

      hadron_tm_isTruthMatched[iter] = 1;
      hadron_tm_fraction[iter]       = fraction;
      hadron_tm_otherE[iter]         = other_energy;

      hadron_tm_PDGCode[iter]        = hadron_trajectory->GetPDGCode();
      hadron_tm_trackID[iter]        = hadron_trajectory->GetTrackId();
      hadron_tm_beginKE[iter]        = hadron_trajectory->GetInitialMomentum().E() - hadron_trajectory->GetInitialMomentum().M();

      //=============================================================
      //! Find matching MCTrack and get Cliff's topology
      //=============================================================
        int topology = -2;
        int first_inelastic = 0;
        double first_inelastic_endp = -9.9;
        double first_inelastic_distance = 0.0;
        if ( MCTracksVector.size() > 0 ) {
          debug()<<"Have MCTracks"<<endmsg;
          topology = -1;
          bool foundMCT = false;
          for (unsigned int i_mct = 0; i_mct < MCTracksVector.size(); ++i_mct) {
            debug()<<"MCTrack "<<i_mct<<endmsg;

            // If this track is a child or not primary track skip it
            if (MCTracksVector[i_mct].GetParentTrackId() != -1) continue;

            // MC track's PDG and energy
            int mct_pdg = MCTracksVector[i_mct].GetPDGCode();
            double mct_E = MCTracksVector[i_mct].GetInitialMomentum().E();

            // if found parent mctrack
            // It looks to me like we're verifying that the info we found so
            // far agrees with the stuff we found in the previous truthmatching
            // block.
            if (mct_pdg == hadron_tm_PDGCode[iter] && mct_E == hadron_trajectory->GetInitialMomentum().E() ) {
              foundMCT = true;

              // Find the first daughters that are not a small angle elastic scatter
              bool isElastic = true;
              int parent_idx = i_mct;
              double mct_p;
              std::vector<int> DaughtersVector;
              debug()<<"  Matched to TM track. Entering isElastic loop"<<endmsg;
              int counter_isElastic = 0;
              while (isElastic && counter_isElastic < 50) {
                counter_isElastic++;
                mct_p = MCTracksVector[parent_idx].GetFinalMomentum().P();
                first_inelastic_distance += sqrt(pow(MCTracksVector[parent_idx].GetInitialPosition().x()-MCTracksVector[parent_idx].GetFinalPosition().x(),2.0)
                    +pow(MCTracksVector[parent_idx].GetInitialPosition().y()-MCTracksVector[parent_idx].GetFinalPosition().y(),2.0)
                    +pow(MCTracksVector[parent_idx].GetInitialPosition().z()-MCTracksVector[parent_idx].GetFinalPosition().z(),2.0));
                DaughtersVector = MCTracksVector[parent_idx].GetDaughterIndices();
                debug()<<"    Daughter size is "<<DaughtersVector.size()<<endmsg;
                if (DaughtersVector.empty()) isElastic = false;
                else {
                  for (std::vector<int>::iterator itD = DaughtersVector.begin(); itD != DaughtersVector.end(); ++itD) {
                    int d_idx = *itD;
                    debug()<<"      Daughter PDG: "<<MCTracksVector[d_idx].GetPDGCode()<<" Process: "<<MCTracksVector[d_idx].GetProcessName()<<endmsg;

                    //Elastic scatter
                    if ( MCTracksVector[d_idx].GetPDGCode() == mct_pdg && MCTracksVector[d_idx].GetProcessName() == "IsTheMotherParticle" ) {
                      parent_idx = d_idx; //update parent_idx to get new daughters

                      //it's elastic, but is it a small scattering angle (<15 deg)?
                      Gaudi::LorentzVector tmp_daught_4p(MCTracksVector[d_idx].GetInitialMomentum());
                      Gaudi::LorentzVector tmp_parent_4p(MCTracksVector[parent_idx].GetFinalMomentum());
                      double scat_angle_rad = calcOpeningAngle( tmp_daught_4p, tmp_parent_4p );
                      if (scat_angle_rad < 0.2618) isElastic = true;
                      else isElastic = false;

                      debug()<<"      Daughter angle is "<<scat_angle_rad<<endmsg;
                      break;
                    }

                    //Inelastic scatter
                    else if (MCTracksVector[d_idx].GetProcessName() != "IsTheMotherParticle" && MCTracksVector[d_idx].GetProcessName() != "hadElastic") {
                      isElastic = false;
                      debug()<<"      Not an elastic scatter"<<endmsg;
                      break;
                    }

                    //Keep looking
                    else {
                      debug()<<"      Other"<<endmsg;
                    }
                  }
                }
              }//while loop search for inelastic or high-angle elastic scatter

              // Probably one of those damn elastic interactions that changes
              // charge.
              // Default to using the original primary particle
              if (counter_isElastic>=50) {
                parent_idx = i_mct;
                DaughtersVector = MCTracksVector[i_mct].GetDaughterIndices();
                mct_p = MCTracksVector[i_mct].GetFinalMomentum().P();
                first_inelastic_distance = sqrt(pow(MCTracksVector[i_mct].GetInitialPosition().x()-MCTracksVector[i_mct].GetFinalPosition().x(),2.0)
                    +pow(MCTracksVector[i_mct].GetInitialPosition().y()-MCTracksVector[i_mct].GetFinalPosition().y(),2.0)
                    +pow(MCTracksVector[i_mct].GetInitialPosition().z()-MCTracksVector[i_mct].GetFinalPosition().z(),2.0));
              }

              debug()<<"  Exited isElastic loop"<<endmsg;

              //final momentum - after all small angle elastic scatters
//              hadron_tm_endMomentum[iter] = mct_p;

              //now get information about the first hadronic interaction for reweighting, starting with the daughter particle the previous loops ends with
              counter_isElastic = 0;
              int tmp_part_idx = parent_idx;
              while (first_inelastic==0 && counter_isElastic < 50) {
                counter_isElastic++;
                first_inelastic_endp = MCTracksVector[tmp_part_idx].GetFinalMomentum().P();

                //already counted the distance for the initial particle in the previous sum
                if (counter_isElastic!=0) first_inelastic_distance +=
                    sqrt(
                      pow(MCTracksVector[tmp_part_idx].GetInitialPosition().x()-MCTracksVector[tmp_part_idx].GetFinalPosition().x(),2.0) +
                      pow(MCTracksVector[tmp_part_idx].GetInitialPosition().y()-MCTracksVector[tmp_part_idx].GetFinalPosition().y(),2.0) +
                      pow(MCTracksVector[tmp_part_idx].GetInitialPosition().z()-MCTracksVector[tmp_part_idx].GetFinalPosition().z(),2.0));
                std::vector<int> tmp_dvec = MCTracksVector[tmp_part_idx].GetDaughterIndices();
                if (tmp_dvec.empty()) {
                  first_inelastic = 0;
                  break;
                }
                else {
                  isElastic = false;
                  bool isDecay = false;
                  int particle_counter = 0;
                  for (std::vector<int>::iterator itD = tmp_dvec.begin(); itD != tmp_dvec.end(); ++itD) {
                    int d_idx = *itD;
                    if (MCTracksVector[d_idx].GetProcessName() == "IsTheMotherParticle") {
                      tmp_part_idx = d_idx;
                      isElastic = true;
                      break;
                    }
                    else if (MCTracksVector[d_idx].GetProcessName() == "Decay") {
                      isDecay = true;
                      break;
                    }
                    if (MCTracksVector[d_idx].GetPDGCode() == mct_pdg) particle_counter++;
                  }
                  if (isElastic) continue;
                  else if (isDecay) {
                    first_inelastic=2;
                    break;
                  }
                  else if (particle_counter > 0) first_inelastic = 4;
                  else first_inelastic = 1;
                }
              }
              if (first_inelastic!=2 && first_inelastic_endp  < 25.0*CLHEP::MeV) first_inelastic = 0;

              //event topologies
              bool isDecay = false;
              bool isHighAngleElastic = false;

              bool isChargedPion = false;
              bool isProton = false;
              bool isNeutron = false;
              bool isPi0 = false;
              for (unsigned int n = 0; n < DaughtersVector.size(); ++n) {
                int daught_index = DaughtersVector[n];

                if (MCTracksVector[daught_index].GetProcessName() == "Decay") isDecay = true;
                else if (MCTracksVector[daught_index].GetProcessName() == "IsTheMotherParticle") isHighAngleElastic = true;
                else if ( MCTracksVector[daught_index].GetPDGCode() == 111) isPi0 = true;
                else if ( abs(MCTracksVector[daught_index].GetPDGCode()) == 2212 ) isProton = true;
                else if ( abs(MCTracksVector[daught_index].GetPDGCode()) == 2112 ) isNeutron = true;
                else if ( abs(MCTracksVector[daught_index].GetPDGCode()) == 211 ) isChargedPion = true;

                if (mct_pdg == 211) debug()<<"  Mother ("<<mct_p<<") Daughter ("<<MCTracksVector[daught_index].GetInitialMomentum().P()<<","
                  <<MCTracksVector[daught_index].GetPDGCode()<<") "
                    <<n<<" created by: "<<MCTracksVector[daught_index].GetProcessName()<<endmsg;

              }

              //NI: 0   AIF: 1   DIF: 2   QEX: 3   K: 4  Elastic-K: 5
              if (abs(mct_pdg) == 211) {
                if (mct_p < 25.0*CLHEP::MeV)                  topology = 0;
                else if (isDecay)                             topology = 2;
                else if (isHighAngleElastic)                  topology = 5;
                else if (isChargedPion)                       topology = 4;
                else if (isPi0)                               topology = 3;
                else                                          topology = 1;
              }
              else if (abs(mct_pdg) == 2212) {
                if (mct_p < 100.0*CLHEP::MeV)            topology = 0;
                else if (isHighAngleElastic)             topology = 5;
                else if (isProton)                       topology = 4;
                else if (isNeutron)                      topology = 3;
                else                                     topology = 1;
              }
              debug()<<"  Topology "<<(mct_pdg)<<": "<<topology<<endmsg;
              break; //break out of loop over MCTracksVector
            }//if found parent mctrack
          }//loop over mctracks vector to find the parent mctrack

          if (!foundMCT) {
            counter("MISSED_MCT")++;
          }
          else {
            counter("MISSED_MCT")+= 0;
          }
        }//if there are mc tracks

        // How did this prong truly meet its fate?
	hadron_tm_destructCode[iter] = topology;
//        hadron_tm_firstInelasticCode[iter] = first_inelastic;
//        hadron_tm_firstInelasticEndP[iter] = first_inelastic_endp;
//        hadron_tm_firstInelasticDistance[iter] = first_inelastic_distance;

//        TG4Trajectory* hadron_daughter_trajectory = getMostEnergeticDaughter(hadron_trajectory);
//        if (hadron_daughter_trajectory) {     
//          hadron_tm_PDGCodeDaughter[iter] = hadron_daughter_trajectory->GetPDGCode();
//        } 


  }//Loop over hadron prongs

  if( pionPart) {
    if ( !pionPart->isMultiMass() ){
      pion_score  = pionPart->score();
      pion_score1 = pionPart->getDoubleData(ParticleExtraDataDefs::dEdXScore1());
      pion_score2 = pionPart->getDoubleData(ParticleExtraDataDefs::dEdXScore2());
    }
  }

  if (secondaryProtonsFound) {
    for(ProngVect::iterator secprongs = secondaryProtonProngs.begin(); secprongs != secondaryProtonProngs.end(); ++secprongs) {
      double secproton_pion_score = -9999., secproton_pion_score1 = -9999., secproton_pion_score2 = -999.;

      SmartRefVector<Minerva::Particle> secpionparticles = (*secprongs)->particlesWithIDCode(Minerva::Particle::Pion, "dEdX");
      if (secpionparticles.size()>0) {
        SmartRef<Minerva::Particle> secprongPart = secpionparticles[0];
        if( !secprongPart->isMultiMass() ) {
          secproton_pion_score  = secprongPart->score();
          secproton_pion_score1 = secprongPart->getDoubleData(ParticleExtraDataDefs::dEdXScore1());
          secproton_pion_score2 = secprongPart->getDoubleData(ParticleExtraDataDefs::dEdXScore2());

          sec_protons_pion_scores.push_back( secproton_pion_score );
          sec_protons_pion_scores1.push_back( secproton_pion_score1 );
          sec_protons_pion_scores2.push_back( secproton_pion_score2 );
        }
      }
    }
  }

  //--------------------------------------------------------------------------------
  //! Reconstruct Neutrino Energy based on CCQE hypotheses with muon information
  //--------------------------------------------------------------------------------
  double enu_muon = -9999., enu_proton = -9999., q2_muon = -9999., q2_proton = -9999.;
  if (has_muon) {
    enu_muon = PhysicsCalculator->nuEnergyCCQE( muon_E, muon_P, muon_theta, muon_charge, m_nuCCQEBindingEnergyMeV );
    q2_muon  = PhysicsCalculator->qSquaredCCQE( muon_E, muon_P, muon_theta, muon_charge, m_nuCCQEBindingEnergyMeV );
  }

  //--------------------------------------------------------
  //! Reconstruct Neutrino Energy using proton kinematics
  //--------------------------------------------------------
  if (protonFound && has_muon)
    enu_proton = m_protonUtils->nuEnergyCCQE(proton_theta,proton_T,m_nuCCQEBindingEnergyMeV);

  if (protonFound && !has_lepton) {
    enu_proton = m_protonUtils->nuEnergyNCE(proton_theta,proton_T,m_nuCCQEBindingEnergyMeV);
    q2_proton = m_protonUtils->qSquaredCCQE(proton_T,m_nuCCQEBindingEnergyMeV,"elastic");
    debug() << "Proton (NCE hypothesis) E_nu = " << enu_proton << " and Q^2 = " << q2_proton << endmsg;
  }


  //Stuff for NukeCCInclusive-style nuclear targets
  // Check to see if the vertex is close to a nuclear target
  //const auto& vtxPos = event->interactionVertex()->position();
  //SmartRef<Minerva::Vertex> vtx = event->interactionVertex();
  Gaudi::XYZPoint vtxPos = event->interactionVertex()->position();
  int targetCode = 0;
  const Minerva::NuclearTarget *targ = getExactTarget( vtxPos, targetCode );

  if( targ )
    debug() << "   This is a passive target event with TargetCode = " << targetCode << endmsg;

  //Use these to fill the containers
  vector<int> intVec;
  vector<double> doubleVec;

  //ccqeHyp->filtertaglist()->setOrAddFilterTag("is_cc", isCC);
  ccqeHyp->setIntData( "in_fiducial_area",   isFiducialXY );
  ccqeHyp->setIntData( "in_analyzable_area", isAnalyzableXY );

  //--------------------------
  //set things specific to Nuke
  int targetID = (targetCode - targetCode%1000) / 1000;
  int targetZ  = targetCode % 1000;
  ccqeHyp->setIntData("target_code",targetCode );
  ccqeHyp->setIntData("targetID",   targetID   );
  ccqeHyp->setIntData("targetZ",    targetZ    );

  //TODO: Do we actually use ref_targZ?
  // store the relation of the vertex to all passive targets
  getReferenceZ( intVec, vtxPos.x(), vtxPos.y() );
  ccqeHyp->setContainerIntData( "ref_targZ", intVec );
  getReferenceDistToDivision( doubleVec, vtxPos.x(), vtxPos.y() );
  ccqeHyp->setContainerDoubleData( "ref_dist_to_division", doubleVec );
  getDistToTarget( doubleVec, vtxPos.z() );
  ccqeHyp->setContainerDoubleData( "ref_dist_to_target", doubleVec );

  // cut on these for passive target event selection
  if( targ )
    {
      vector<double> orig_short_vtx; //what was the vertex before adustment to target?
      if( ccqeHyp->getContainerDoubleData("orig_short_vtx", orig_short_vtx ) )
        ccqeHyp->setDoubleData("target_zDist", targ->distanceToTarget( orig_short_vtx[2] ) );
      else
        ccqeHyp->setDoubleData("target_zDist", targ->distanceToTarget( vtxPos.z() ) );

      ccqeHyp->setDoubleData("target_dist_to_division", targ->distanceToDivision( vtxPos.x(), vtxPos.y() ) );
    }

  //store the module/plane of the vertex (-999 for events in passive material)
  int vtxMod = -999, vtxPlane = -999;
  bool inPlane = getVtxPlane( vtxPos.z(), vtxMod, vtxPlane );


  if ( has_lepton ) {
    ccqeHyp->filtertaglist()->setOrAddFilterTag( "in_plane", inPlane );
    ccqeHyp->setIntData( "vtx_module", vtxMod );
    ccqeHyp->setIntData( "vtx_plane", vtxPlane );
   }
////////******** Machine Learning stuff: TargetID, TargetZ filled here separately **********////////////

  if(event->hasContainerDoubleData("ANN_vtx"))
  {

  const Gaudi::XYZPoint vtxPosANN(event->getContainerDoubleData("ANN_vtx")[0], event->getContainerDoubleData("ANN_vtx")[1], event->getContainerDoubleData("ANN_vtx")[2]);

  bool ANNisFiducialXY   = FiducialPointTool->isFiducial(vtxPosANN, m_signalApothem, m_usFiducialModID, m_dsFiducialModID);
  bool ANNisAnalyzableXY = FiducialPointTool->isFiducial(vtxPosANN, m_analyzeApothem,m_usAnalyzableModID,m_dsAnalyzableModID);
  ccqeHyp->setIntData( "ANN_in_fiducial_area",   ANNisFiducialXY );
  ccqeHyp->setIntData( "ANN_in_analyzable_area",   ANNisAnalyzableXY );
  
  int ANNtargetCode = 0;
  const Minerva::NuclearTarget *targ = getExactTarget( vtxPosANN, ANNtargetCode );

  //set things specific to Nuke
  int ANNtargetID = (ANNtargetCode - ANNtargetCode%1000) / 1000;
  int ANNtargetZ  = ANNtargetCode % 1000;
  ccqeHyp->setIntData("ANN_target_code",ANNtargetCode );
  ccqeHyp->setIntData("ANN_targetID",   ANNtargetID   );
  ccqeHyp->setIntData("ANN_targetZ",    ANNtargetZ    );

  if( targ )
    debug() << "   This is a passive target event with TargetCode = " << targetCode << endmsg;
  
  getReferenceZ( intVec, vtxPosANN.x(), vtxPosANN.y() );
  ccqeHyp->setContainerIntData( "ANN_ref_targZ", intVec );
  getReferenceDistToDivision( doubleVec, vtxPosANN.x(), vtxPosANN.y() );
  ccqeHyp->setContainerDoubleData( "ANN_ref_dist_to_division", doubleVec );
  getDistToTarget( doubleVec, vtxPosANN.z() );
  ccqeHyp->setContainerDoubleData( "ANN_ref_dist_to_target", doubleVec );

  if( targ )
    {
      ccqeHyp->setDoubleData("ANN_target_zDist", targ->distanceToTarget( event->getContainerDoubleData("ANN_vtx")[2] ));
      ccqeHyp->setDoubleData("ANN_target_dist_to_division", targ->distanceToDivision(event->getContainerDoubleData("ANN_vtx")[0], event->getContainerDoubleData("ANN_vtx")[1]));
    }

////////////////////Bring recoil stuff for ML here /////////////////////

  if ( has_lepton ) {
    //this isn't the best way to get the mass of the struck nucleon, but it's fine for now
    const double nucleonMass = ( (MinervaUnits::M_p + MinervaUnits::M_n ) / 2.0);

    //! Get the recoil energy from the recoil clusters.
    string ANNrecoilChannel = getRecoilChannel( muon_charge, ANNtargetID, ANNtargetZ );
    double ANNrecoilE = getRecoilEnergy( event, muonProng, ANNrecoilChannel, ccqeHyp );
    double visibleE = getvisibleenergy(event, muonProng);

    SmartRefVector<Minerva::IDCluster> idClusters = getAnalyzableNonMuIDClusters( event, muonProng, true );
    SmartRefVector<Minerva::ODCluster> odClusters = getAnalyzableNonMuODClusters( event, muonProng, true );

    double recoilE_passive = m_caloUtils->applyCalConsts( idClusters, odClusters , "Default" , false);
    //double recoilE_passive = m_caloUtils->applyCalConsts( idClusters, odClusters , recoilChannel , true);
    double CorrectedMuonEnergy;
    if(event->hasContainerDoubleData("muon_corrected_p"))
      {
    CorrectedMuonEnergy = event->getContainerDoubleData("muon_corrected_p")[3];
    }  
    //double ANNtotalE = muon_E + ANNrecoilE;
    double ANNtotalE = CorrectedMuonEnergy + ANNrecoilE;

    //double bindingE   = ( muon_charge < 0 ) ? m_nuCCQEBindingE : m_antinuCCQEBindingE;
    //double qsquaredQE = PhysicsCalculator->qSquaredCCQE( muonPart->momentumVec().E(), muonPart->momentumVec().P(), muon_theta, muon_charge, m_nuCCQEBindingEnergyMeV );

    // Calculate W & Bjorken x
    //double qsquared   = PhysicsCalculator->qSquared( totalE, muonE, mu_theta ); original
    double ANNqsquared   = PhysicsCalculator->qSquared( ANNtotalE, muon_E, muon_theta );
    double ANNwmass = PhysicsCalculator->W( nucleonMass, ANNrecoilE, ANNqsquared );
    double xBj = PhysicsCalculator->xBjorken( ANNqsquared, nucleonMass, ANNrecoilE );

    //TODO: Does this clobber the numbers that CCQENu came up with on its own?
    //      I recommend we think about putting this in a separate CCInclusive hypothesis.
    ccqeHyp->setEnergy( ANNtotalE );
    ccqeHyp->setDoubleData( AnaFilterTags::RecoilEnergy(), ANNrecoilE );
    ccqeHyp->setDoubleData("ANN_neutrino_E", ANNtotalE );
    ccqeHyp->setDoubleData("ANN_recoil_E", ANNrecoilE );
    ccqeHyp->setDoubleData("ANN_visible_E",visibleE);
    ccqeHyp->setDoubleData("recoil_passivecorrected",recoilE_passive);
    
    } // if ( has_lepton )

  }

///////////////////////////////////////****** Machine Learning Stuff ends here *********////////////////////////////////////////////////////////////////


  if ( has_lepton ) {
    ccqeHyp->setNeutrinoHelicity( PhysicsCalculator->getHelicity( muon_charge ) );
    ccqeHyp->setLeptonEnergy( muonPart->momentumVec() );
    /////ccqeHyp->setEnergy( enu_muon ); //should I also use enu_proton somewhow?
    ccqeHyp->setVertex( event->interactionVertex() );

    ccqeHyp->setDoubleData( "muon_P", muon_P );
    ccqeHyp->setDoubleData( "muon_Px", muon_Px );
    ccqeHyp->setDoubleData( "muon_Py", muon_Py );
    ccqeHyp->setDoubleData( "muon_Pz", muon_Pz );
    ccqeHyp->setDoubleData( "muon_E", muon_E );
    ccqeHyp->setDoubleData( "muon_T", muon_T );
    ccqeHyp->setDoubleData( "muon_theta", muon_theta );
    ccqeHyp->setDoubleData( "muon_score", muon_score );
    ccqeHyp->setDoubleData("muon_theta_biasUp",muon_theta_biasUp);
    ccqeHyp->setDoubleData("muon_theta_biasDown",muon_theta_biasDown);
    ccqeHyp->setDoubleData("muon_qp",muon_qp );
    ccqeHyp->setDoubleData("muon_qpqpe",muon_qpqpe);
    ccqeHyp->setDoubleData( "enu_muon", enu_muon );
    ccqeHyp->setDoubleData( "proton_P_fromdEdx", proton_P );
    ccqeHyp->setDoubleData( "proton_E_fromdEdx", proton_E );
    ccqeHyp->setDoubleData( "proton_T_fromdEdx", proton_T );
    ccqeHyp->setDoubleData( "proton_Px_fromdEdx", proton_Px );
    ccqeHyp->setDoubleData( "proton_Py_fromdEdx", proton_Py );
    ccqeHyp->setDoubleData( "proton_Pz_fromdEdx", proton_Pz );
    ccqeHyp->setDoubleData( "proton_theta_fromdEdx", proton_theta );
    ccqeHyp->setDoubleData( "proton_calib_energy", totalCalibEPrimaryProton );
    ccqeHyp->setIntData( "proton_patternRec",proton_patternRec);
    ccqeHyp->setDoubleData( "enu_proton", enu_proton );
    ccqeHyp->setDoubleData( "proton_score", proton_score );
    ccqeHyp->setDoubleData( "proton_score1", proton_score1 );
    ccqeHyp->setDoubleData( "proton_score2", proton_score2 );
    ccqeHyp->setDoubleData( "pion_score", pion_score );
    ccqeHyp->setDoubleData( "pion_score1", pion_score1 );
    ccqeHyp->setDoubleData( "pion_score2", pion_score2 );
    ccqeHyp->setContainerDoubleData( "pion_E", pion_E );
    ccqeHyp->setContainerDoubleData( "pion_T", pion_T );
    ccqeHyp->setContainerDoubleData( "pion_P", pion_P );
    ccqeHyp->setContainerDoubleData( "pion_Px", pion_Px );
    ccqeHyp->setContainerDoubleData( "pion_Py", pion_Py );
    ccqeHyp->setContainerDoubleData( "pion_Pz", pion_Pz );
    ccqeHyp->setContainerDoubleData( "pion_theta", pion_theta );
    ccqeHyp->setContainerDoubleData( "pion_phi", pion_phi );
    ccqeHyp->setContainerDoubleData("hadron_pion_E_corr",  hadron_pion_E_corr);
    ccqeHyp->setContainerDoubleData("hadron_pion_p_corr",  hadron_pion_p_corr);
    ccqeHyp->setContainerDoubleData("hadron_pion_px_corr", hadron_pion_px_corr);
    ccqeHyp->setContainerDoubleData("hadron_pion_py_corr", hadron_pion_py_corr);
    ccqeHyp->setContainerDoubleData("hadron_pion_pz_corr", hadron_pion_pz_corr);
    ccqeHyp->setContainerDoubleData("hadron_pion_E_recoil_corr", hadron_pion_E_recoil_corr);
    ccqeHyp->setContainerDoubleData( "hadron_piFit_scoreLLR", hadron_piFit_scoreLLR);
    ccqeHyp->setContainerDoubleData("hadron_piFit_score1", hadron_piFit_score1);
    ccqeHyp->setContainerIntData( "hadron_isTracker", hadron_isTracker);
    ccqeHyp->setContainerIntData( "hadron_isSideECAL", hadron_isSideECAL);
    ccqeHyp->setContainerIntData( "hadron_isForked", hadron_isForked);
    ccqeHyp->setContainerIntData( "hadron_isODMatch", hadron_isODMatch);
    ccqeHyp->setContainerIntData( "hadron_isExiting", hadron_isExiting);
    ccqeHyp->setIntData("hadron_number", hadron_number);
    ccqeHyp->setContainerIntData("hadron_tm_PDGCode",hadron_tm_PDGCode);
    ccqeHyp->setContainerIntData("hadron_tm_destructCode", hadron_tm_destructCode);
    ccqeHyp->setContainerIntData("hadron_tm_trackID", hadron_tm_trackID);
    ccqeHyp->setContainerDoubleData("hadron_tm_beginKE", hadron_tm_beginKE );
    ccqeHyp->setContainerDoubleData("hadron_tm_otherE", hadron_tm_otherE );
    ccqeHyp->setContainerIntData("hadron_tm_isTruthMatched", hadron_tm_isTruthMatched);
    ccqeHyp->setContainerDoubleData("hadron_tm_fraction", hadron_tm_fraction );
    ccqeHyp->setContainerIntData("hadron_endMichel_category", hadron_endMichel_category);
    ccqeHyp->setContainerDoubleData("hadron_endMichel_slice_energy", hadron_endMichel_slice_energy);
    ccqeHyp->setContainerIntData("hadron_endMichel_ndigits", hadron_endMichel_ndigits);
    ccqeHyp->setContainerDoubleData("hadron_endMichel_energy", hadron_endMichel_energy);
    ccqeHyp->setContainerDoubleData("piFit_score1_Mass_biasUp", hadron_piFit_score1_Mass_Up);
    ccqeHyp->setContainerDoubleData("piFit_score1_Mass_biasDown", hadron_piFit_score1_Mass_Down);
    ccqeHyp->setContainerDoubleData("piFit_score1_BetheBloch_biasUp",   hadron_piFit_score1_BetheBloch_Up);
    ccqeHyp->setContainerDoubleData("piFit_score1_BetheBloch_biasDown", hadron_piFit_score1_BetheBloch_Down);
    ccqeHyp->setContainerDoubleData("piFit_score1_Birks", hadron_piFit_score1_Birks);
    ccqeHyp->setContainerDoubleData( "pion_E_Mass_biasUp", hadron_pion_E_Mass_Up );
    ccqeHyp->setContainerDoubleData( "pion_E_Mass_biasDown", hadron_pion_E_Mass_Down );
    ccqeHyp->setContainerDoubleData( "pion_E_BetheBloch_biasUp", hadron_pion_E_BetheBloch_Up );
    ccqeHyp->setContainerDoubleData( "pion_E_BetheBloch_biasDown", hadron_pion_E_BetheBloch_Down );
    ccqeHyp->setContainerDoubleData( "pion_theta_biasUp", hadron_pion_theta_biasUp );
    ccqeHyp->setContainerDoubleData( "pion_theta_biasDown", hadron_pion_theta_biasDown );
    ccqeHyp->setContainerDoubleData( "pion_E_Birks", hadron_pion_E_Birks );
    ccqeHyp->setContainerDoubleData( "sec_protons_E_fromdEdx", sec_protons_E );
    ccqeHyp->setContainerDoubleData( "sec_protons_P_fromdEdx", sec_protons_P );
    ccqeHyp->setContainerDoubleData( "sec_protons_T_fromdEdx", sec_protons_T );
    ccqeHyp->setContainerDoubleData( "sec_protons_Px_fromdEdx", sec_protons_Px );
    ccqeHyp->setContainerDoubleData( "sec_protons_Py_fromdEdx", sec_protons_Py );
    ccqeHyp->setContainerDoubleData( "sec_protons_Pz_fromdEdx", sec_protons_Pz );
    ccqeHyp->setContainerDoubleData( "sec_protons_theta_fromdEdx", sec_protons_theta );
    ccqeHyp->setContainerDoubleData( "sec_protons_T_fromCalo", sec_protons_T_Calo );
    ccqeHyp->setContainerIntData("sec_protons_patternRec",sec_protons_patternRec);
    ccqeHyp->setContainerDoubleData( "sec_protons_proton_scores", sec_protons_proton_scores );
    ccqeHyp->setContainerDoubleData( "sec_protons_proton_scores1", sec_protons_proton_scores1 );
    ccqeHyp->setContainerDoubleData( "sec_protons_proton_scores2", sec_protons_proton_scores2 );
    ccqeHyp->setContainerDoubleData( "sec_protons_pion_scores", sec_protons_pion_scores );
    ccqeHyp->setContainerDoubleData( "sec_protons_pion_scores1", sec_protons_pion_scores1 );
    ccqeHyp->setContainerDoubleData( "sec_protons_pion_scores2", sec_protons_pion_scores2 );
    //ccqeHyp->setQ2( q2_muon );
    ccqeHyp->setDoubleData( "Q2_CCQE", q2_muon );

    //----------------------------------------------------------------
    ccqeHyp->setDoubleData("proton_endPointX",proton_endPointX);
    ccqeHyp->setDoubleData("proton_endPointY",proton_endPointY);
    ccqeHyp->setDoubleData("proton_endPointZ",proton_endPointZ);

    ccqeHyp->setDoubleData("proton_startPointX",proton_startPointX);
    ccqeHyp->setDoubleData("proton_startPointY",proton_startPointY);
    ccqeHyp->setDoubleData("proton_startPointZ",proton_startPointZ);

    ccqeHyp->setDoubleData("proton_theta",proton_theta);
    ccqeHyp->setDoubleData("proton_thetaX",proton_thetaX);
    ccqeHyp->setDoubleData("proton_thetaY",proton_thetaY);
    ccqeHyp->setDoubleData("proton_phi",proton_phi);
    //-----------------------------------------------------------------
    //---------------------------------
    //--!dE/dx biasing information
    ////------------------------------
    ccqeHyp->setDoubleData( "proton_score1_Mass_biasUp", proton_score1_Mass_Up );
    ccqeHyp->setDoubleData( "proton_score1_Mass_biasDown", proton_score1_Mass_Down );
    ccqeHyp->setDoubleData( "proton_score1_BetheBloch_biasUp", proton_score1_BetheBloch_Up );
    ccqeHyp->setDoubleData( "proton_score1_BetheBloch_biasDown", proton_score1_BetheBloch_Down );
    ccqeHyp->setDoubleData( "proton_score1_MEU_biasUp", proton_score1_MEU_Up );
    ccqeHyp->setDoubleData( "proton_score1_MEU_biasDown", proton_score1_MEU_Down );
    ccqeHyp->setDoubleData( "proton_score1_Birks_bias", proton_score1_Birks );
    ccqeHyp->setDoubleData( "proton_E_Mass_biasUp", proton_E_Mass_Up );
    ccqeHyp->setDoubleData( "proton_E_Mass_biasDown", proton_E_Mass_Down );
    ccqeHyp->setDoubleData( "proton_E_BetheBloch_biasUp", proton_E_BetheBloch_Up );
    ccqeHyp->setDoubleData( "proton_E_BetheBloch_biasDown", proton_E_BetheBloch_Down );
    ccqeHyp->setDoubleData( "proton_E_MEU_biasUp", proton_E_MEU_Up );
    ccqeHyp->setDoubleData( "proton_E_MEU_biasDown", proton_E_MEU_Down );
    ccqeHyp->setDoubleData( "proton_E_Birks_bias", proton_E_Birks );

    ccqeHyp->setContainerDoubleData( "sec_protons_score1_Mass_biasUp", sec_protons_score1_Mass_Up );
    ccqeHyp->setContainerDoubleData( "sec_protons_score1_Mass_biasDown", sec_protons_score1_Mass_Down );
    ccqeHyp->setContainerDoubleData( "sec_protons_score1_BetheBloch_biasUp", sec_protons_score1_BetheBloch_Up );
    ccqeHyp->setContainerDoubleData( "sec_protons_score1_BetheBloch_biasDown", sec_protons_score1_BetheBloch_Down );
    ccqeHyp->setContainerDoubleData( "sec_protons_score1_MEU_biasUp", sec_protons_score1_MEU_Up );
    ccqeHyp->setContainerDoubleData( "sec_protons_score1_MEU_biasDown", sec_protons_score1_MEU_Down );
    ccqeHyp->setContainerDoubleData( "sec_protons_score1_Birks_bias", sec_protons_score1_Birks );
    ccqeHyp->setContainerDoubleData( "sec_protons_E_Mass_biasUp", sec_protons_E_Mass_Up );
    ccqeHyp->setContainerDoubleData( "sec_protons_E_Mass_biasDown", sec_protons_E_Mass_Down );
    ccqeHyp->setContainerDoubleData( "sec_protons_E_BetheBloch_biasUp", sec_protons_E_BetheBloch_Up );
    ccqeHyp->setContainerDoubleData( "sec_protons_E_BetheBloch_biasDown", sec_protons_E_BetheBloch_Down );
    ccqeHyp->setContainerDoubleData( "sec_protons_E_MEU_biasUp", sec_protons_E_MEU_Up );
    ccqeHyp->setContainerDoubleData( "sec_protons_E_MEU_biasDown", sec_protons_E_MEU_Down );
    ccqeHyp->setContainerDoubleData( "sec_protons_E_Birks_bias", sec_protons_E_Birks );

    //! Set Node Cut Branches
    ccqeHyp->setContainerIntData("pion_nNodes",pion_nNodes);
    ccqeHyp->setContainerDoubleData("pion_lastnode_Q0",pion_lastnode_Q0);
    ccqeHyp->setContainerDoubleData("pion_lastnode_Q1",pion_lastnode_Q1);
    ccqeHyp->setContainerDoubleData("pion_lastnode_Q2",pion_lastnode_Q2);
    ccqeHyp->setContainerDoubleData("pion_lastnode_Q3",pion_lastnode_Q3);
    ccqeHyp->setContainerDoubleData("pion_lastnode_Q4",pion_lastnode_Q4);
    ccqeHyp->setContainerDoubleData("pion_lastnode_Q5",pion_lastnode_Q5);

    //-------------------
    // !ESC proton variables
    //------------------
    ccqeHyp->setIntData("proton_prongType",proton_prongType);
    ccqeHyp->setContainerDoubleData("proton_nodes_E",proton_nodes_E);
    ccqeHyp->setContainerDoubleData("proton_nodes_nodesNormE",proton_nodes_nodesNormE);

    ccqeHyp->setContainerIntData( "sec_protons_prongType", sec_protons_prongType );
    ccqeHyp->setContainerIntData( "sec_protons_nodes_index", sec_protons_nodes_index );
    ccqeHyp->setContainerDoubleData( "sec_protons_nodes_nodesNormE", sec_protons_nodes_nodesNormE );
    ccqeHyp->setContainerDoubleData( "sec_protons_nodes_E", sec_protons_nodes_E );
  } // if ( has_lepton )

  if( has_lepton && m_doNeutron )
    {
      debug()<<"Filling Neutron Braches"<<endmsg;
      std::vector<double> neutron4P;
      neutron4P.push_back(neutronPart->momentumVec().X() ); neutron4P.push_back(neutronPart->momentumVec().Y() ); neutron4P.push_back(neutronPart->momentumVec().Z() );
      neutron4P.push_back(neutronPart->momentumVec().E() );
      ccqeHyp->setContainerDoubleData( "HadronE", neutron4P );
      int NTracks, NIncTracks;

      for( uint i = 0; i< m_neutronBlobUtils->GetIntDataMemberNames().size(); i++)
        {
          const std::string key = m_neutronBlobUtils->GetIntDataMemberNames()[i];
          for( uint j = 0; j<m_neutronBlobUtils->GetDataInt()[ key ].size();j++ ) {
            debug()<<"Key "<<key<<" contains "<<m_neutronBlobUtils->GetDataInt()[ key ][j]<<endmsg;

            if (key == "BlobNTracks" ) NTracks += m_neutronBlobUtils->GetDataInt()[ key ][j];
            else if (key == "BlobNIncTracks") NIncTracks += m_neutronBlobUtils->GetDataInt()[ key ][j];
          }
          ccqeHyp->setContainerIntData( key, m_neutronBlobUtils->GetDataInt()[ key ] );

        }

      for( uint i = 0; i< m_neutronBlobUtils->GetDoubleDataMemberNames().size(); i++)
        {
          const std::string key = m_neutronBlobUtils->GetDoubleDataMemberNames()[i];
          for( uint j = 0; j<m_neutronBlobUtils->GetDataDouble()[ key ].size();j++ ) debug()<<"Key "<<key<<" contains "<<m_neutronBlobUtils->GetDataDouble()[ key ][j]<<endmsg;
          ccqeHyp->setContainerDoubleData( key, m_neutronBlobUtils->GetDataDouble()[ key ] );
        }

      ccqeHyp->setIntData( "EvtHasNBlobTracks", NTracks );
      ccqeHyp->setIntData( "EvtHasNBlobIncTracks", NIncTracks );

      m_neutronBlobUtils->ClearMap();
    }

  if ( has_lepton ) {
    //this isn't the best way to get the mass of the struck nucleon, but it's fine for now
    const double nucleonMass = ( (MinervaUnits::M_p + MinervaUnits::M_n ) / 2.0);

    //! Get the recoil energy from the recoil clusters.
    string recoilChannel = getRecoilChannel( muon_charge, targetID, targetZ );
    double recoilE = getRecoilEnergy( event, muonProng, recoilChannel, ccqeHyp );
    double visibleE = getvisibleenergy(event, muonProng);

    SmartRefVector<Minerva::IDCluster> idClusters = getAnalyzableNonMuIDClusters( event, muonProng, true );
    SmartRefVector<Minerva::ODCluster> odClusters = getAnalyzableNonMuODClusters( event, muonProng, true );

    double recoilE_passive = m_caloUtils->applyCalConsts( idClusters, odClusters , "Default" , false);
    //double recoilE_passive = m_caloUtils->applyCalConsts( idClusters, odClusters , recoilChannel , true);

    double totalE = muon_E + recoilE;

    //double bindingE   = ( muon_charge < 0 ) ? m_nuCCQEBindingE : m_antinuCCQEBindingE;
    //double qsquaredQE = PhysicsCalculator->qSquaredCCQE( muonPart->momentumVec().E(), muonPart->momentumVec().P(), muon_theta, muon_charge, m_nuCCQEBindingEnergyMeV );

    // Calculate W & Bjorken x
    //double qsquared   = PhysicsCalculator->qSquared( totalE, muonE, mu_theta ); original
    double qsquared   = PhysicsCalculator->qSquared( totalE, muon_E, muon_theta );
    double wmass = PhysicsCalculator->W( nucleonMass, recoilE, qsquared );
    double xBj = PhysicsCalculator->xBjorken( qsquared, nucleonMass, recoilE );

    //TODO: Does this clobber the numbers that CCQENu came up with on its own?
    //      I recommend we think about putting this in a separate CCInclusive hypothesis.
    ccqeHyp->setEnergy( totalE );
    ccqeHyp->setW( wmass );
    ccqeHyp->setDoubleData( "Q2_Inclusive", qsquared );
    ccqeHyp->setXBjorken( xBj );
    ccqeHyp->setInelasticity( recoilE / totalE );
    ccqeHyp->setDoubleData( AnaFilterTags::RecoilEnergy(), recoilE );
    ccqeHyp->setDoubleData("recoil_E", recoilE );
    ccqeHyp->setDoubleData("visible_E",visibleE);
    ccqeHyp->setDoubleData("recoil_passivecorrected",recoilE_passive);

    //add wide time window kinematics
    {
      double recoilE_wide = ccqeHyp->getDoubleData("recoil_E_wide_window"); //we trust that it is filled in getRecoilEnergy
      double totalE_wide = muon_E + recoilE_wide;

      //Calculate E, W, x, y, q2 for the shfited recoil
      double qsquared_wide   = PhysicsCalculator->qSquared( totalE_wide, muon_E, muon_theta );
      double W_wide = PhysicsCalculator->W( nucleonMass, recoilE_wide, qsquared_wide );
      double xBj_wide = PhysicsCalculator->xBjorken( qsquared_wide, nucleonMass, recoilE_wide );

      ccqeHyp->setDoubleData( "E_wide_window", totalE_wide );
      ccqeHyp->setDoubleData( "W_wide_window", W_wide );
      ccqeHyp->setDoubleData( "Q2_wide_window", qsquared_wide );
      ccqeHyp->setDoubleData( "x_wide_window", xBj_wide );
      ccqeHyp->setDoubleData( "y_wide_window", recoilE_wide / totalE_wide );
    }
    //--------------------------------------------------------------------------------
  } // if ( has_lepton )

  if (has_muon) {
    fillMinosMuonBranches( ccqeHyp, muonProng );
    //! Fill the systematic shifts. See SystematicShiftsExtensions.cc.
    fillSystematicShiftsBranches( ccqeHyp, event );
  }

  info() << "Exiting MasterAnaDev::interpretEvent() now ....." << endmsg;

  return StatusCode::SUCCESS;

} //End of interpretEvent( )



//================================================================================
// BEGIN helper functions
//================================================================================
bool MasterAnaDev::adjustDNNVertexIntoCenterOfSegment( SmartRef<Minerva::Vertex>& vtx, double zCenter, const SmartRef<Minerva::Particle>& muonPart ) const
{
  const Minerva::Particle * muon = muonPart.target();
  Gaudi::XYZTVector muon_vtx = muon->startPos();
  Gaudi::XYZTVector muon_p   = muon->momentumVec();
  double projx = muon_vtx.x() + ( muon_p.x() / muon_p.z() ) * ( zCenter - muon_vtx.z() );
  double projy = muon_vtx.y() + ( muon_p.y() / muon_p.z() ) * ( zCenter - muon_vtx.z() );
  Gaudi::XYZTVector projVtx = Gaudi::XYZTVector( projx, projy, zCenter, muon_vtx.t() );
  vtx->setPosition( Gaudi::XYZPoint( projVtx.x(), projVtx.y(), projVtx.z() ) );
  return true;
}

bool MasterAnaDev::adjustVertexIntoPassiveTarget( SmartRef<Minerva::Vertex>& vtx, const SmartRef<Minerva::Particle>& muonPart ) const
{
  //! Do we think this event comes from a passive target?
  int targetCode = 0;
  const Minerva::NuclearTarget *targ = getTarget( vtx->position(), targetCode );
  info() << "***** targetID, z start, z center, z end, vertex pos, targetCode = " << targ->getTargetID() << " " << targ->getZStart() << " " << targ->getZCenter() << " " << targ->getZEnd() << " " << vtx->position().z() << " " << targetCode << endmsg;
  if( ! targ )
    return false;

  //! If this event comes from a target, project to the zCenter of it and update the vertex position to that point.
  Gaudi::XYZTVector projVtx(0,0,0,0);
  m_nukeTool->projectToTarget( muonPart.target(), targ, projVtx );
  vtx->setPosition( Gaudi::XYZPoint( projVtx.x(), projVtx.y(), projVtx.z() ) );
  return true;
}

//=======================================================
// Get the target this event is inside or close enough to
//=======================================================
const Minerva::NuclearTarget* MasterAnaDev::getTarget( const Gaudi::XYZPoint& p, int& targetCode ) const
{
  const Minerva::NuclearTarget *targ = 0;
  targetCode = 0;

  //! Loop over targets to find one in which the muon might have originated
  for( vector<Minerva::NuclearTarget*>::const_iterator itTarg = m_targets.begin(); itTarg != m_targets.end(); ++itTarg )
    {
      debug() << "Checking muon against target " << (*itTarg)->getName() << endmsg;

      double zLow  = m_nukeTool->getZ_nPlanesUS( (*itTarg), m_nplanes_USVertexLimit );
      double zHigh = m_nukeTool->getZ_nPlanesDS( (*itTarg), m_nplanes_DSVertexLimit );

      if( ! m_usePlanesZCut )
        {
          zLow  = (*itTarg)->getZStart() - m_USVertexCut;
          zHigh = (*itTarg)->getZEnd() + m_DSVertexCut;
        }


      //! if the muon's origin is in the vicinity of the nuclear target (Z only)
      if( FiducialPointTool->isFiducial( p, 99999999.9, zLow, zHigh ) )
        {
          debug() << "   Z position of the event is consistent with an event from the target" << endmsg;
          targ = *itTarg;

          //! See if we can get a Material at this point.
          const Material *mat = m_nukeTool->getSectionMaterial( targ, p.x(), p.y() );

          if( mat )
            debug() << " Mat Z before cast: " << mat->Z() << endl;

          //! targetCode = 1000*targetID + targetZ.  targetZ is 0 if material is unknown
          targetCode = 1000 * (*itTarg)->getTargetID();
          //When calculating target Z:
          //Use Nint to round to nearest integer Z when calculating target code
          //Mixtures of materials (such as Fe + Mn) will result in mat->Z returning
          //A non-int Z which will get truncated when C casts it to int.
          //Work around this with Nint, which rounds a double to the nearest int.

          if( mat )
            targetCode += TMath::Nint(mat->Z());

          if( mat )
            debug() << "Mat Z after cast: " << targetCode << ", " << TMath::Nint( mat->Z() ) << endl;

          break;
        }
    }//foreach(targets)

  return targ;
}

const Minerva::NuclearTarget* MasterAnaDev::getTarget( const Minerva::GenMinInteraction* truth, int& targetCode ) const
{
  return getTarget( Gaudi::XYZPoint( truth->Vtx().x(), truth->Vtx().y(), truth->Vtx().z() ), targetCode );
}

//=======================================================
// Get the exact target (point must be in material)
//=======================================================
const Minerva::NuclearTarget* MasterAnaDev::getExactTarget( const Gaudi::XYZPoint& p, int& targetCode ) const
{
  targetCode = 0;
  //When calculating target Z:
  //Use Nint to round to nearest integer Z when calculating target code
  //Mixtures of materials (such as Fe + Mn) will result in mat->Z returning
  //A non-int Z which will get truncated when C casts it to int.
  //Work around this with Nint, which rounds a double to the nearest int.

  const Minerva::NuclearTarget *targ = m_nukeTool->getNuclearTarget( p );
  if( targ )
    {
      debug() << "  Point is exactly in a nuclear target" << endmsg;

      //! targetCode = 1000*targetID + targetZ.  targetZ is 0 if material is unknown
      targetCode = 1000 * targ->getTargetID();
      if( targ->getTargetID() != 3 && 4978.95 <= p.z() && p.z() <= 4983.98 ) targetCode = 1000 * 3;
      if( 5637.53 <= p.z() && p.z() <= 5645.73 ) targetCode = 1000 * 4;

      //! See if we can get a Material at this point.
      const Material *mat = m_nukeTool->getSectionMaterial( targ, p.x(), p.y() );
      if( mat )
        targetCode += TMath::Nint( mat->Z() );
      else
        warning() << "Cannot identify a section material.  Something may be wrong." << endmsg;
    }

  return targ;
}

const Minerva::NuclearTarget* MasterAnaDev::getExactTarget( const Minerva::GenMinInteraction* truth, int& targetCode ) const
{
  return getExactTarget( Gaudi::XYZPoint( truth->Vtx().x(), truth->Vtx().y(), truth->Vtx().z() ), targetCode );
}


//=======================================================
// Get the Z of the would-be nucleus in reference targets
//=======================================================
bool MasterAnaDev::getReferenceZ( vector<int>& intVec, double x, double y ) const
{
  debug() << "  getReferenceZ..." << endmsg;
  intVec.clear();
  for( vector<Minerva::NuclearTarget*>::const_iterator itTarg = m_targets.begin(); itTarg != m_targets.end(); ++itTarg )
    {
      const Material *mat = m_nukeTool->getSectionMaterial( *itTarg, x, y );
      int refZ = 0;
      //When calculating refZ:
      //Use Nint to round to nearest integer Z when calculating target code
      //Mixtures of materials (such as Fe + Mn) will result in mat->Z returning
      //A non-int Z which will get truncated when C casts it to int.
      //Work around this with Nint, which rounds a double to the nearest int.
      if( mat )
        refZ = TMath::Nint(mat->Z());
      else
        debug() << "Could not get Z for reference target = " << (*itTarg)->getName() << " at (x,y) = ( " << x << ", " << y << " ) " << endmsg;
      intVec.push_back( refZ );
    }
  return true;
}

//================================================================
// Get the distance to the would-be division in reference targets
//================================================================
bool MasterAnaDev::getReferenceDistToDivision( vector<double>& dVec, double x, double y ) const
{
  debug() << "  getReferenceDistToDivision..." << endmsg;
  dVec.clear();

  for( vector<Minerva::NuclearTarget*>::const_iterator itTarg = m_targets.begin(); itTarg != m_targets.end(); ++itTarg )
    dVec.push_back( (*itTarg)->distanceToDivision( x, y ) );

  return true;

}

bool MasterAnaDev::getDistToTarget( vector<double>& dVec, double z ) const
{
  debug() << "  getDistToTarget..." << endmsg;
  dVec.clear();

  for( vector<Minerva::NuclearTarget*>::const_iterator itTarg = m_targets.begin(); itTarg != m_targets.end(); ++itTarg )
    dVec.push_back( z - (*itTarg)->getZCenter() );

  return true;
}


//====================================================================//

StatusCode MasterAnaDev::fillPlaneLowZMap( )
{
  m_planeLowZMap.clear();

  typedef vector<const Minerva::DePlane*> DePlanes;

  bool printPlanes = false;
  const string viewStr( "0XUV" );
  const double zStep = .00001 * CLHEP::mm;

  DePlanes planes = m_idDet->getDePlanes();
  for( DePlanes::const_iterator plane = planes.begin(); plane != planes.end(); ++plane )
    {
      const double zCenter = (*plane)->getZCenter();
      const int modno      = (*plane)->getPlaneID().module();
      const int planeno    = (*plane)->getPlaneID().plane();
      const int view       = int( (*plane)->getView() );

      // do not use
      if( modno > MAX_TRACKER_MODULE )
        break;

      debug() << "module number: " << modno <<" And plane: "<<planeno<<" with view: "<<view<<endmsg; //Oscar

      // firgure out how far this plane extends US in z before encroaching on another plane
      
/*      double deltaZ = PLANE_WIDTH / 2.;
      if( 2 == planeno )
        deltaZ += PLANE_GAP / 2.;
      else
      {
        if( modno < FIRST_TRACKER_MODULE )
          deltaZ += MOD_GAP_NUKE /2.;
        else
          deltaZ += MOD_GAP / 2.;
      }
*/
//Oscar Stuff/////

    double deltaZ;
      //Just for the Target and tracker region
      if(modno < 85)
        {
          deltaZ = PLANE_WIDTH / 2.;
          if( 2 == planeno )
            deltaZ += PLANE_GAP / 2.;
          else
            {
              if( modno < FIRST_TRACKER_MODULE )
                deltaZ += MOD_GAP_NUKE /2.;
              else
                {
                  deltaZ += MOD_GAP / 2.;

                }
            }
        }

      //Just for the ECAL By Oscar
      if(modno >84 && modno < 95)
        {
          deltaZ = ECAL_HALF_WIDTH; //11.21
          if(planeno == 1)
            deltaZ += GAP_PLANE1_ECAL; // 3.84
          if(modno == 85 && planeno == 1)
            deltaZ += GAP_MOD85; //0.43
        }
      //Just for the HCAL By Oscar
      if(modno > 94)
        {
          deltaZ = HCAL_HALF_WIDTH; //34.63
        }
      const double zLow    = zCenter - deltaZ;
      m_planeLowZMap[zLow] = std::pair<int,int>(modno, planeno);

      if(printPlanes)
        {
          Gaudi::XYZPoint testPoint(0., 0., zCenter );
          while( (*plane)->isInside(testPoint) )
            testPoint.SetCoordinates(0., 0., testPoint.z() - zStep);
          const double zStart = testPoint.z() + zStep/2.;

          testPoint.SetCoordinates(0., 0., zCenter );
          while( (*plane)->isInside(testPoint) )
            testPoint.SetCoordinates(0., 0., testPoint.z() + zStep);
          const double zStop = testPoint.z() - zStep/2.;

    //        const double zHigh   = zCenter + deltaZ; //an approximation.  there's a different gap between planes and modules

          //Just for ECAL Oscar Stuff///
          if(modno > 84 && modno < 95)
            {
              deltaZ = GAP_ECAL_HIGH; //20.46
              if(planeno == 1)
                deltaZ += GAP_PLANE1_ECAL; // 3.84
              if(modno == 85 && planeno == 1)
                deltaZ += GAP_MOD85; //0.43
              if(modno == 94 && planeno ==2)
                deltaZ += GAP_LEAD_94;//5.82
            }
          //Just for HCAL
          if(modno > 94)
            {
              deltaZ = HCAL_HIGH; //47.31
            }
        
          double zHigh   = zLow + deltaZ;
    
          std::cout
            << "Module: " << modno
            << ", Plane: " << planeno
            << ", view: " << viewStr[view]
            << Form("; zCenter: %.2fmm",(*plane)->getZCenter())
            << Form("; zRange[est]: [%.2f,%.2f]", zLow, zHigh)
            << Form("; zRange: [%.2f,%.2f]", zStart, zStop)
            << std::endl;
        }
    }

  return StatusCode::SUCCESS;
}

bool MasterAnaDev::getVtxPlane( double z, int & modno, int & planeno ) const
{
  bool foundPlane = false;

  if( z < m_planeLowZMap.begin()->first )
    {
      debug() << "MasterAnaDev::getVtxPlane: z point " << z << " is farther upstream than any plane." << endmsg;
      counter( "vtxplane_lowZ" )++;
    }
  else if( m_planeLowZMap.rbegin()->first < z )
    {
      debug() << "MasterAnaDev::getVtxPlane: z point " << z << " is farther upstream than any plane." << endmsg;
      counter( "vtxplane_highZ" )++;
    }
  else
    {
      std::map< double, std::pair<int,int> >::const_iterator itPlane = --m_planeLowZMap.lower_bound( z );
      const double zLow = itPlane->first;
      modno   = itPlane->second.first;
      planeno = itPlane->second.second;


      // firgure out how far DS this plane extends in z before encroaching on another plane, starting from the US point of the plane
 /* 
      double deltaZ = PLANE_WIDTH;
      deltaZ += PLANE_GAP / 2.;
      if( modno < FIRST_TRACKER_MODULE )
        deltaZ += MOD_GAP_NUKE /2.;
      else
        deltaZ += MOD_GAP / 2.;
*/


/////Oscar Stuff/////

    double deltaZ;
      //Just for target and tracker region (added by Oscar)
      if(modno < 85)
        {
           deltaZ = PLANE_WIDTH;
           deltaZ += PLANE_GAP / 2.;
          if( modno < FIRST_TRACKER_MODULE )
            deltaZ += MOD_GAP_NUKE /2.;
          else
            deltaZ += MOD_GAP / 2.;
        }
      //Just for ECAL
      if(modno >84 && modno < 95)
        {
          deltaZ = GAP_ECAL_HIGH; //20.46
          if(planeno == 1)
            deltaZ += GAP_PLANE1_ECAL; // 3.84
          if(modno == 85 && planeno == 1)
            deltaZ += GAP_MOD85; //0.43
          if(modno == 94 && planeno ==2)
            deltaZ += GAP_LEAD_94;//5.82
        }
      //Just for HCAL
      if(modno > 94)
        {
          deltaZ = HCAL_HIGH; //47.31
        }
      //check that this z position is inside the range of a plane
      const double zHigh = zLow + deltaZ;
      if( z < zHigh )
        foundPlane = true;
      else
        debug() << "getVtxPlane: Position within fiducial limits but not inside a plane at z = " << z << endmsg;

      debug() << Form(
                      "  getVtxPlane: z = %f - modno = %d, planeno = %d - lowz = %f, deltaz = %f, highz = %f",
                      z, modno, planeno, zLow, deltaZ, zHigh )
              << endmsg;
    }

  //set nonsense if a plane was not found
  if( ! foundPlane )
    modno = planeno = -999;

  return foundPlane;
}

// Given a module number, return the correct moduleID
Minerva::ModuleID MasterAnaDev::getModuleID( int mod ) const
{
  Minerva::SubdetID::Subdet subdet;
  if( mod < 23 )
    subdet = Minerva::SubdetID::NuclTargs;
  else if( mod < 85 )
    subdet = Minerva::SubdetID::Tracker;
  else if( mod < 95 )
    subdet = Minerva::SubdetID::ECAL;
  else
    subdet = Minerva::SubdetID::HCAL;
  return Minerva::ModuleID( Minerva::DetectorID::ID, subdet, mod );
}
///////////////////////////////////////////////////////////////////////////////////////

//  Blobbing functions
//============================================================================
//============================================================================
SmartRefVector<Minerva::IDCluster> MasterAnaDev::getBlobbingClusters( Minerva::PhysicsEvent* event, const SmartRef<Minerva::Prong>& muonProng ) const
{
  //this gets in-time clusters that are not XTalk and not on the muon prong
  SmartRefVector<Minerva::IDCluster> blobClusters = getAnalyzableNonMuIDClusters(event, muonProng);

  for( SmartRefVector<Minerva::IDCluster>::iterator i = blobClusters.begin(); i != blobClusters.end(); )
    {
      //mark all clusters as unused
      (*i)->setHistory( Minerva::MinervaRecoBase::Unused );

      //remove HCAL and LA clusters which are analyzable for recoil energy but not blobbing
      if( Minerva::IDCluster::HCAL == (*i)->subdet() || Minerva::IDCluster::LowActivity == (*i)->type() )
        i = blobClusters.erase(i);
      else
        ++i;
    }

  return blobClusters;
}

SmartRefVector<Minerva::IDCluster> MasterAnaDev::removeUsedClusters( SmartRefVector<Minerva::IDCluster>& clusters ) const
{
  SmartRefVector<Minerva::IDCluster> used;
  for( SmartRefVector<Minerva::IDCluster>::iterator i = clusters.begin(); i != clusters.end(); )
    {
      if( Minerva::MinervaRecoBase::Used == (*i)->history() )
        {
          used.push_back( *i );
          i = clusters.erase(i);
        }
      else
        ++i;
    }
  return used;
}

double MasterAnaDev::fillBlobBranches( Minerva::PhysicsEvent* event, SmartRefVector<Minerva::IDCluster>& idclusters, const string& prefix, int color ) const
{
  SmartRefVector<Minerva::ODCluster> odclusters;
  return fillBlobBranches( event, idclusters, odclusters, prefix, color );
}

double MasterAnaDev::fillBlobBranches( Minerva::PhysicsEvent* event, SmartRefVector<Minerva::IDCluster>& idclusters, SmartRefVector<Minerva::ODCluster>& odclusters, const string& prefix, int color ) const
{
  //return how much energy is in the subdets used for ccqe-recoil E
  double eInNuclTrkEcal = 0.;

  SmartRefVector<Minerva::IDCluster> nucl, tracker, ecal, hcal;
  separateClustersBySubdet( idclusters, nucl, tracker, ecal, hcal );

  event->setIntData(  string( "blob_" + prefix + "_nClus"), idclusters.size() + odclusters.size() );
  event->setIntData(  string( "blob_" + prefix + "_nClus_nucl"), nucl.size() );
  event->setIntData(  string( "blob_" + prefix + "_nClus_tracker"), tracker.size() );
  event->setIntData(  string( "blob_" + prefix + "_nClus_ecal"), ecal.size() );
  event->setIntData(  string( "blob_" + prefix + "_nClus_hcal"), hcal.size() );
  event->setIntData(  string( "blob_" + prefix + "_nClus_od"), odclusters.size() );

  {
    double nuclE = m_caloUtils->applyCalConsts( nucl );
    double trackerE = m_caloUtils->applyCalConsts( tracker );
    double ecalE = m_caloUtils->applyCalConsts( ecal );
    double hcalE = m_caloUtils->applyCalConsts( hcal );
    double odE = m_caloUtils->applyCalConsts( odclusters );
    double totE = nuclE + trackerE + ecalE + hcalE + odE;
    eInNuclTrkEcal = nuclE + trackerE + ecalE;

    debug() << "Tracker Energy =" << trackerE <<endl;
    debug() << "Hcal Energy =" << hcalE <<endl;
    debug() << "OD Energy =" << odE <<endl;
    debug() << "Ecal Energy =" << ecalE <<endl;
    debug() << "Nucl Energy =" << nuclE <<endl;
    debug() << "Total Energy =" << totE <<endl;

    event->setDoubleData( string( "blob_" + prefix + "_E"), totE );
    event->setDoubleData( string( "blob_" + prefix + "_E_nucl"),    nuclE );
    event->setDoubleData( string( "blob_" + prefix + "_E_tracker"), trackerE );
    event->setDoubleData( string( "blob_" + prefix + "_E_ecal"),    ecalE );
    event->setDoubleData( string( "blob_" + prefix + "_E_hcal"),    hcalE );
    event->setDoubleData( string( "blob_" + prefix + "_E_od"),      odE );
  }


  if( 0 != color )
    {
      for( SmartRefVector<Minerva::IDCluster>::iterator i = idclusters.begin(); i != idclusters.end(); ++i )
        m_hitTaggerTool->applyColorTag( *i, color );


      for( SmartRefVector<Minerva::ODCluster>::iterator i = odclusters.begin(); i != odclusters.end(); ++i )
        m_hitTaggerTool->applyColorTag( *i, color );
    }

  //if MC, add e by particle type
  if( haveNeutrinoMC() )
    {
      for( std::map<MCPartType::t_type, string>::iterator i = PART_TYPE_MAP.begin(); i != PART_TYPE_MAP.end(); ++i )
        {
          MCPartType::t_type type = i->first;
          string& name = i->second;

          SmartRefVector<Minerva::IDCluster> nucl_fType = getClustersCreatedBy<Minerva::IDCluster,Minerva::IDDigit>(nucl, type);
          SmartRefVector<Minerva::IDCluster> tracker_fType = getClustersCreatedBy<Minerva::IDCluster,Minerva::IDDigit>(tracker, type);
          SmartRefVector<Minerva::IDCluster> ecal_fType = getClustersCreatedBy<Minerva::IDCluster,Minerva::IDDigit>(ecal, type);
          SmartRefVector<Minerva::IDCluster> hcal_fType = getClustersCreatedBy<Minerva::IDCluster,Minerva::IDDigit>(hcal, type);
          SmartRefVector<Minerva::ODCluster> odclusters_fType = getClustersCreatedBy<Minerva::ODCluster,Minerva::ODDigit>(odclusters, type);


          double nuclE = m_caloUtils->applyCalConsts( nucl_fType );
          double trackerE = m_caloUtils->applyCalConsts( tracker_fType );
          double ecalE = m_caloUtils->applyCalConsts( ecal_fType );
          double hcalE = m_caloUtils->applyCalConsts( hcal_fType );
          double odE = m_caloUtils->applyCalConsts( odclusters_fType );
          double totE = nuclE + trackerE + ecalE + hcalE + odE;


          event->setDoubleData( string( "blob_" + prefix + "_E_" + name),    totE );
          event->setDoubleData( string( "blob_" + prefix + "_E_nucl_" + name),    nuclE );
          event->setDoubleData( string( "blob_" + prefix + "_E_tracker_" + name), trackerE );
          event->setDoubleData( string( "blob_" + prefix + "_E_ecal_" + name),    ecalE );
          event->setDoubleData( string( "blob_" + prefix + "_E_hcal_" + name),    hcalE );
          event->setDoubleData( string( "blob_" + prefix + "_E_od_" + name),      odE );
        }

    }

  return eInNuclTrkEcal;
}
/////////////////////////////////////////////////////////////////////////////////////////
void MasterAnaDev::separateClustersBySubdet(
                                            const SmartRefVector<Minerva::IDCluster>& clusters,
                                            SmartRefVector<Minerva::IDCluster>& nucl,
                                            SmartRefVector<Minerva::IDCluster>& tracker,
                                            SmartRefVector<Minerva::IDCluster>& ecal,
                                            SmartRefVector<Minerva::IDCluster>& hcal
                                            ) const
{
  for( SmartRefVector<Minerva::IDCluster>::const_iterator i = clusters.begin(); i != clusters.end(); ++i )
    {
      Minerva::IDCluster::Subdet subdet = (*i)->subdet();
      if(      Minerva::IDCluster::NuclTargs == subdet )
        nucl.push_back( *i );
      else if( Minerva::IDCluster::Tracker == subdet )
        tracker.push_back( *i );
      else if( Minerva::IDCluster::ECAL == subdet )
        ecal.push_back( *i );
      else if( Minerva::IDCluster::HCAL == subdet )
        hcal.push_back( *i );
    }
}

//Note could add radius to branch names
void MasterAnaDev::declareBlobBranches( const string& prefixBefore )
{
  //loop over vertex blob radius variations for all blobs but the recoil blob
  const unsigned int nVtxRadii = ( prefixBefore == RECOILBLOB_PREFIX ) ? 1 : m_nuVtxBlobRadii.size();
  for( unsigned int i = 0; i != nVtxRadii; ++i )
    {
      const string suffix = (0==i) ? "" : Form("_alt%02d",i);
      const string prefix = prefixBefore + suffix;

      //counts
      declareIntEventBranch(    "blob_"+prefix+"_nBlobs", 0 );
      declareIntEventBranch(    "blob_"+prefix+"_nClus", 0 );
      declareIntEventBranch(    "blob_"+prefix+"_nClus_nucl", 0 );
      declareIntEventBranch(    "blob_"+prefix+"_nClus_tracker", 0 );
      declareIntEventBranch(    "blob_"+prefix+"_nClus_ecal", 0 );
      declareIntEventBranch(    "blob_"+prefix+"_nClus_hcal", 0 );
      declareIntEventBranch(    "blob_"+prefix+"_nClus_od", 0 );

      //energy by subdet
      declareDoubleEventBranch( "blob_"+prefix+"_E", 0. );
      declareDoubleEventBranch( "blob_"+prefix+"_E_nucl", 0. );
      declareDoubleEventBranch( "blob_"+prefix+"_E_tracker", 0. );
      declareDoubleEventBranch( "blob_"+prefix+"_E_ecal", 0. );
      declareDoubleEventBranch( "blob_"+prefix+"_E_hcal", 0. );
      declareDoubleEventBranch( "blob_"+prefix+"_E_od", 0. );

      //energy by what created it
      for( std::map<MCPartType::t_type, string>::iterator i = PART_TYPE_MAP.begin(); i != PART_TYPE_MAP.end(); ++i )
        {
          declareDoubleEventBranch( "blob_"+prefix+"_E_" + i->second, 0. );
          declareDoubleEventBranch( "blob_"+prefix+"_E_nucl_" + i->second, 0. );
          declareDoubleEventBranch( "blob_"+prefix+"_E_tracker_" + i->second, 0. );
          declareDoubleEventBranch( "blob_"+prefix+"_E_ecal_" + i->second, 0. );
          declareDoubleEventBranch( "blob_"+prefix+"_E_hcal_" + i->second, 0. );
          declareDoubleEventBranch( "blob_"+prefix+"_E_od_" + i->second, 0. );
        }
    }//loop over vertex blob radii
}




//=============================================================================
// Get particle type info from truth information
//=============================================================================
template<typename ClusterType, typename DigitType>
SmartRefVector<ClusterType> MasterAnaDev::getClustersCreatedBy( SmartRefVector<ClusterType> clusters, MCPartType::t_type primType, MCPartType::t_type fsType ) const
{

  SmartRefVector<ClusterType> rval;
  for( typename SmartRefVector<ClusterType>::iterator clus = clusters.begin(); clus != clusters.end(); ++clus )
    {

      if( primType == getParticleCreatedBy<ClusterType,DigitType>( *clus, MCPartStatus::kPrimary ) )
        {
          if( MCPartType::kAny == fsType || fsType == getParticleCreatedBy<ClusterType,DigitType>( *clus, MCPartStatus::kFinal ) )
            rval.push_back(*clus);
        }

    }

  return rval;
}

template<typename ClusterType, typename DigitType>
MCPartType::t_type MasterAnaDev::getParticleCreatedBy( SmartRef<ClusterType> cluster, MCPartStatus::t_type status ) const
{
  if( status == MCPartStatus::kPrimary )
    {
      double lowneutronEnergy = 0;
      double midneutronEnergy = 0;
      double highneutronEnergy = 0;
      double protonEnergy = 0;
      double mesonEnergy = 0;
      double emEnergy = 0;
      double muEnergy = 0;
      double xtalkEnergy = 0;
      double otherEnergy = 0;

      double maxEnergy = -1.;
      MCPartType::t_type maxType = MCPartType::kOther;

      SmartRefVector<DigitType> digits = cluster->digits();
      for( typename SmartRefVector<DigitType>::iterator dig = digits.begin(); dig != digits.end(); ++dig )
        {
          const double normE = (*dig)->normEnergy();

          vector<const DigitType*> digVec;
          digVec.push_back( *dig );
          const Minerva::TG4Trajectory* constTrajectory = NULL;
          double fraction = 0., other_energy = 0.0;

          if( ! digVec.empty() )
            TruthMatcher->getPrimaryTG4Trajectory(digVec, constTrajectory, fraction, other_energy );
          if( 0 == constTrajectory )
            { //no trajectory means Xtalk
              xtalkEnergy += normE;
              if( maxEnergy < xtalkEnergy )
                {
                  maxEnergy = xtalkEnergy;
                  maxType = MCPartType::kXtalk;
                }
            }
          else
            {
              const int pdg = constTrajectory->GetPDGCode();
              switch( pdg )
                {
                case 2112: //neutron
                  {
                    double neutron_mom = constTrajectory->GetInitialMomentum().P();
                    //double neutron_ke = neutron_mom*neutron_mom/2/MinervaUnits::M_n; //<--- nonrelativistic?
                    double neutron_ke = sqrt( pow( neutron_mom,2 ) + pow( MinervaUnits::M_n, 2 ) ) - MinervaUnits::M_n;
                    if( neutron_ke < 50 * CLHEP::MeV )
                      {
                        lowneutronEnergy += normE;
                        if( maxEnergy < lowneutronEnergy )
                          {
                            maxEnergy = lowneutronEnergy;
                            maxType = MCPartType::klown;
                          }
                      }
                    else if( neutron_ke < 150 * CLHEP::MeV )
                      {
                        midneutronEnergy += normE;
                        if( maxEnergy < midneutronEnergy )
                          {
                            maxEnergy = midneutronEnergy;
                            maxType = MCPartType::kmidn;
                          }
                      }
                    else
                      {
                        highneutronEnergy += normE;
                        if( maxEnergy < highneutronEnergy )
                          {
                            maxEnergy = highneutronEnergy;
                            maxType = MCPartType::khighn;
                          }
                      }
                  }
                  break;

                case 2212: //proton
                  protonEnergy += normE;
                  if( maxEnergy < protonEnergy )
                    {
                      maxEnergy = protonEnergy;
                      maxType = MCPartType::kp;
                    }
                  break;

                case 211:  //pi+
                case -211: //pi-
                case 321:  //K+
                case -321: //K-
                case 311:  //K0
                  mesonEnergy += normE;
                  if( maxEnergy < mesonEnergy )
                    {
                      maxEnergy = mesonEnergy;
                      maxType = MCPartType::kmeson;
                    }
                  break;

                case 13:  //mu+
                case -13: //mu-
                  muEnergy += normE;
                  if( maxEnergy < muEnergy )
                    {
                      maxEnergy = muEnergy;
                      maxType = MCPartType::kmu;
                    }
                  break;

                case 11:   //e+
                case -11:  //e-
                case 111:  //pi0
                case -111: //pi0
                case 22:   //gamma
                  emEnergy += normE;
                  if( maxEnergy < emEnergy )
                    {
                      maxEnergy = emEnergy;
                      maxType = MCPartType::kem;
                    }
                  break;

                default:
                  otherEnergy += normE;
                  if( maxEnergy < otherEnergy )
                    {
                      maxEnergy = otherEnergy;
                      maxType = MCPartType::kOther;
                    }
                  break;
                } //end swtich on particle pdg
            } //end if found trajectory
        } //end loop over digits

      return maxType;

    }
  else if( status == MCPartStatus::kFinal )
    {
      vector<const ClusterType*> clusVec;
      clusVec.push_back( cluster );

      const Minerva::TG4Trajectory* constTrajectory = NULL;
      double fraction = 0., other_energy = 0.0;

      if( ! clusVec.empty() )
        TruthMatcher->getTG4Trajectory( clusVec, constTrajectory, fraction, other_energy );
      if( 0 == constTrajectory )
        return MCPartType::kXtalk;
      else
        {
          const int pdg = constTrajectory->GetPDGCode();
          switch( pdg )
            {
            case 2112: //neutron
              {
                double neutron_mom = constTrajectory->GetInitialMomentum().P();
                double neutron_ke = sqrt( pow( neutron_mom,2 ) + pow( MinervaUnits::M_n, 2 ) ) - MinervaUnits::M_n;
                if( neutron_ke < 50 * CLHEP::MeV )
                  return MCPartType::klown;
                else if( neutron_ke < 150 * CLHEP::MeV )
                  return MCPartType::kmidn;
                else
                  return MCPartType::khighn;
              }

            case 2212: //proton
              return MCPartType::kp;

            case 211:  //pi+
            case -211: //pi-
            case 321:  //K+
            case -321: //K-
            case 311:  //K0
              return MCPartType::kmeson;

            case 13:  //mu+
            case -13: //mu-
              return MCPartType::kmu;

            case 11:   //e+
            case -11:  //e-
            case 111:  //pi0
            case -111: //pi0
            case 22:   //gamma
              return MCPartType::kem;

            default:
              return MCPartType::kOther;
            } //end of pdg switch
        } //end if found trajectory
    } //end if primary/final

  warning() << "Asked for particle with unknown status.  status can be Primary or Final." << endmsg;

  return MCPartType::kOther;
}

///////////////////////////////////////////////////////////////////////////////////////////


//==============================================================================
// Momentum by range
//==============================================================================

void MasterAnaDev::calcMomentumByRangeVars(Minerva::PhysicsEvent* event) const {
  // Declare vars
  // Indices here are for each hadron
  std::vector<int>    hadron_ntrack;
  std::vector<double> hadron_track_length;
  std::vector<double> hadron_track_length_area;
  std::vector<double> hadron_track_length_area_trkr;
  std::vector<double> hadron_track_length_area_ecal;
  const int MAX_TRACKS = 10;

  for( int kl = 0; kl < MAX_TRACKS; ++kl ){
    hadron_ntrack.push_back(-1);
    hadron_track_length.push_back(-1.0);
    hadron_track_length_area.push_back(-1.0);
    hadron_track_length_area_trkr.push_back(-1.0);
    hadron_track_length_area_ecal.push_back(-1.0);
  }
  const double tracker_ecal_boundary = 8590.07; // from Ron's document!
  // Collect the hadronic prongs
  Minerva::ProngVect primary_prongs = event->primaryProngs();
  Minerva::ProngVect pion_prongs;
  pion_prongs.clear();
  for(unsigned int prong = 0; prong <primary_prongs.size(); ++prong) {
    if( primary_prongs[prong]->filtertaglist()->filterTagExists("PrimaryHadron")
        &&
        !primary_prongs[prong]->filtertaglist()->filterTagExists("PrimaryMuon") ) {
      pion_prongs.push_back( primary_prongs[prong] );
    }
  }
  info()<<"Number of Hadron Prongs: "<<pion_prongs.size()<<endmsg;

  // Loop prongs
  for(unsigned int prong = 0; prong < pion_prongs.size(); prong++) {
    double track_length = 0.0;
    double track_length_area = 0.0;
    double track_length_area_trkr = 0.0;
    double track_length_area_ecal = 0.0;
    if (prong == MAX_TRACKS) break;
    // loop prong's tracks
    const SmartRefVector<Minerva::Track>& tracks = pion_prongs[prong]->minervaTracks();
    info() << "prong ntracks: " << tracks.size() << endmsg;
    for (SmartRefVector<Minerva::Track>::const_iterator t = tracks.begin(); t != tracks.end(); ++t) {

       // loop track's nodes
       const std::vector<Minerva::Node*>& nodes = (*t)->nodes();

      for (std::vector<Minerva::Node*>::const_iterator n = nodes.begin(); n != nodes.end(); ++n) {
        std::vector<Minerva::Node*>::const_iterator current = n;
        std::vector<Minerva::Node*>::const_iterator next = n; ++next;

        if (next == nodes.end()) break;

        // start and end point of this node
        Gaudi::XYZPoint start( (*current)->state().x(),
                               (*current)->state().y(),
                               (*current)->state().z() );

        Gaudi::XYZPoint end( (*next)->state().x(),
                             (*next)->state().y(),
                             (*next)->state().z() );
        // calculate area density between two nodes
        // use path length tool to get distance
        double area_density = m_pathLengthTool->ColumnarDensity(start,end);
        // multiply by material density
        track_length      += (end - start).R();
        track_length_area += area_density;
        // ignore the case where one node is upstream and the other is downstream
        // of the boundary
        if (start.z() > tracker_ecal_boundary && end.z() > tracker_ecal_boundary) {
          track_length_area_ecal += area_density;
        } else {
          track_length_area_trkr += area_density;
        }
      } // end node loop
    } // end track loop
    info() << "end track loop, filling branches" << endmsg;

      hadron_ntrack[prong]  =  int(tracks.size()) ;

    info() << "track_length_area: " << track_length_area << endmsg;

    hadron_track_length[prong]           = track_length;
    hadron_track_length_area[prong]      = track_length_area;
    hadron_track_length_area_trkr[prong] = track_length_area_trkr;
    hadron_track_length_area_ecal[prong] = track_length_area_ecal;
  } // end prong loop

  event->setContainerIntData(   "hadron_ntrack",                 hadron_ntrack);
  event->setContainerDoubleData("hadron_track_length",           hadron_track_length);
  event->setContainerDoubleData("hadron_track_length_area",      hadron_track_length_area);
  event->setContainerDoubleData("hadron_track_length_area_trkr", hadron_track_length_area_trkr);
  event->setContainerDoubleData("hadron_track_length_area_ecal", hadron_track_length_area_ecal);
}//end of MasterAnaDev::calcMomentumByRangeVars

//==============================================================================
// Pion-Muon Opening Angle
//==============================================================================
double MasterAnaDev::calcOpeningAngle(Gaudi::LorentzVector& pion_vec,
                                     Gaudi::LorentzVector& muon_vec) const {
  double pionMuonAngle = pion_vec.px()*muon_vec.px() + pion_vec.py()*muon_vec.py() + pion_vec.pz()*muon_vec.pz();
  pionMuonAngle /= sqrt( pow(pion_vec.px(),2.0) + pow(pion_vec.py(),2.0) + pow(pion_vec.pz(),2.0) );
  pionMuonAngle /= sqrt( pow(muon_vec.px(),2.0) + pow(muon_vec.py(),2.0) + pow(muon_vec.pz(),2.0) );
  pionMuonAngle = acos(pionMuonAngle);

  return pionMuonAngle;
}

//==============================================================================
// Isolated blob info
//! Get calibrated non-vertex energy for ALL the isolated blobs.
//! Here we have blobbed up "ALL UNUSED" clusters (except XTalk & LowAct) in
//  the event, in time or out of time w.r.t. muon.
//! So activity from overlay will also appear as blobs, if they are in the same
//  time slice.
//==============================================================================
void MasterAnaDev::fillIsoProngInfo(Minerva::PhysicsEvent* event,
                                   double muon_time) const {
  const double kLowerTimeWindow = -20 * CLHEP::ns;
  const double kUpperTimeWindow =  35 * CLHEP::ns;

  warning() << "Entered fillIsoProngInfo()" << endmsg;

  // Create Shower Blob Prongs from remaining "unused" clusters
  std::vector<Minerva::Prong*> *new_iso_ID_prongs = new std::vector<Minerva::Prong*>;
  Minerva::IDClusterVect  unused_clusts_for_blobs = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
  StatusCode sc = m_primaryBlobProngTool->makeShowerBlobProngs( event, unused_clusts_for_blobs, new_iso_ID_prongs );
  if( sc.isFailure() ){
    error() << "Make Isolated Prongs FAILED" << endmsg;
    return;
  }
  //--------------------------------------------------------------------------
  //! You can add other conditions here, for examining the ShowerBlobProngs
  //! before you start adding them to the evtMgr, user has control !
  //--------------------------------------------------------------------------
  if( new_iso_ID_prongs->size() ) addObject( event, *new_iso_ID_prongs );
  //Delete vector new_iso_ID_prongs but not their elements
  delete new_iso_ID_prongs;

  // declare
  int iso_prongs_count = 0;
  int iso_prongs_blob_count = 0;
  double iso_prongs_total_energy = 0.;
  std::vector<double> iso_prong_energy;
  std::vector<double> iso_prong_separation;
  std::vector<int>    iso_prong_nhits;
  std::vector<double> iso_prong_x;
  std::vector<double> iso_prong_y;
  std::vector<double> iso_prong_z;
  for (int i = 0; i < 30; i++){
    iso_prong_energy.push_back(-1.);
    iso_prong_separation.push_back(-1.);
    iso_prong_nhits.push_back(-1);
    iso_prong_x.push_back(-1.);
    iso_prong_y.push_back(-1.);
    iso_prong_z.push_back(-1.);
  }

  // get interaction vertex
  SmartRef<Minerva::Vertex> vertex = event->interactionVertex();

  // Get ALL prongs, not just primaries.
  Minerva::ProngVect prongs = event->select<Minerva::Prong>( "Used:Unused", "All" );
  warning() << "Event has now " << prongs.size() << " prongs" << endmsg;
  int p_index = 0;
  for( Minerva::ProngVect::iterator p=prongs.begin(); p!=prongs.end() && p_index < 30; ++p, ++p_index  ) {
    warning() << "Prong Creation Signature = "<< (*p)->creationSignature() << endmsg;
    // Only look at Isolated ID prongs
    if( (*p)->creationSignature() == PrimaryBlobProngExtraDataDefs::IsolatedIDBlobProngName() ) {
      // Basic info about prong
      iso_prongs_count          += 1;
      iso_prongs_blob_count     += (*p)->idblobs().size();
      double prong_energy        = m_caloUtils->applyCalConsts( *p, "CCInclusive", false );
      iso_prongs_total_energy   += prong_energy;
      iso_prong_energy[p_index]  = prong_energy;

      warning() << "prong_energy " << prong_energy << endmsg;

      // Loop over the prong's blobs, to get position and nhits
      int prong_nhits = 0;
      double prong_separation = 9999.;
      double prong_x = 9999., prong_y = 9999., prong_z = 9999.;
      Minerva::IDBlobVect blobs = (*p)->idblobs();
      for( Minerva::IDBlobVect::iterator it_blob=blobs.begin(); it_blob!=blobs.end(); it_blob++ ) {
        // Make sure blob is in time with muon
        double time_diff = (*it_blob)->time() - muon_time;
        if( time_diff > kLowerTimeWindow && time_diff < kUpperTimeWindow ) {
          const Gaudi::XYZPoint blob_pos = (*it_blob)->startPoint();
          prong_separation = std::min((blob_pos - vertex->position()).R(),
                                      prong_separation ); // this is in mm
          prong_x = std::min(prong_x, blob_pos.X());
          prong_y = std::min(prong_y, blob_pos.Y());
          prong_z = std::min(prong_z, blob_pos.Z());
        }
        prong_nhits += (*it_blob)->getAllDigits().size();
      }
      iso_prong_x[p_index] = prong_x;
      iso_prong_y[p_index] = prong_y;
      iso_prong_z[p_index] = prong_z;
      iso_prong_nhits[p_index]      = prong_nhits;
      iso_prong_separation[p_index] = prong_separation;
    } // End must be isolated prong
  } //End of for loop over Prongs

  event->setIntData            ("iso_prongs_count",        iso_prongs_count);
  event->setIntData            ("iso_prongs_blob_count",   iso_prongs_blob_count);
  event->setDoubleData         ("iso_prongs_total_energy", iso_prongs_total_energy);
  event->setContainerDoubleData("iso_prong_energy",        iso_prong_energy);
  event->setContainerDoubleData("iso_prong_separation",    iso_prong_separation);
  event->setContainerIntData   ("iso_prong_nhits",         iso_prong_nhits);
  event->setContainerDoubleData("iso_prong_x",             iso_prong_x);
  event->setContainerDoubleData("iso_prong_y",             iso_prong_y);
  event->setContainerDoubleData("iso_prong_z",             iso_prong_z);
}

//==============================================================================
// Return the particle hypothesis with the largest score for a given particle
// type
//==============================================================================

StatusCode MasterAnaDev::getBestParticle(SmartRef<Minerva::Prong> prong,
                                        SmartRef<Minerva::Particle>& particle,
                                        Minerva::Particle::ID partType) const {
  verbose()<<"Enter MasterAnaDev::getBestParticle"<<endmsg;

  SmartRefVector<Minerva::Particle> partHypVec = prong->particles();

//get prong's best dEdX particle hypothesis
  double largest_score = -1.0;

  verbose()<<"  Prong has "<<partHypVec.size()<<" particle hypotheses"<<endmsg;
  for( SmartRefVector<Minerva::Particle>::iterator itPart = partHypVec.begin(); itPart != partHypVec.end(); ++itPart ) {
    double score1 = -9.0;
    if ( (*itPart)->isMultiMass() ) continue;
    if ( (*itPart)->hasDoubleData("score1") ) score1  = (*itPart)->getDoubleData("score1");
    verbose()<<"    score1 "<< score1 <<" part "<<(*itPart)->idcode()<<endmsg;
    if ( score1 > largest_score && (*itPart)->idcode() == partType ) {

      verbose()<<"    Updating particle"<<endmsg;
      particle = *itPart;
      largest_score = score1;

    }
  }
  verbose()<<"Exiting MasterAnaDev::getBestParticle"<<endmsg;
  std::cout << " Inside getBestParticle"<< "\n";
  return StatusCode::SUCCESS;
}

//===================================================================================================================
// Investigate if it is plausible for the true event to be a signal event
// truthIsPlausible( ) is now a pure virtual function of MinervaAnalysisTool, so you must implement it
// It is called automatically by PhysicsEventAnalysisAlg AFTER reconstructEvent() and interpretEvent() are run.
// If it returns false, the event is not included in the analysis ntuple
//===================================================================================================================
bool MasterAnaDev::truthIsPlausible( const Minerva::PhysicsEvent *physEvent ) const {

  info() << "Entering truthIsPlausible( ) now ..... " << endmsg;

  if( NULL == physEvent ) {
    debug() << "Physics Event is a NULL pointer - no sense in proceeding !" << endmsg;
    return false;
  }

  //This might need more work.
  //Look at all prongs. If ONE prong is the mc event as defined below then true else false
  if(!physEvent->hasInteractionVertex()){
    bool isPlausible = false;
    VertexVect used_vtx = physEvent->select<Minerva::Vertex>( "Used", "StartPoint" );
    for( Minerva::VertexVect::iterator itvtx=used_vtx.begin(); itvtx!=used_vtx.end(); ++itvtx ) {
      ProngVect primProngs; ProngVect unattachedProngs;
      m_objectAssociator->getProngs_fromSourceVertex( primProngs, unattachedProngs, *itvtx );
      for( ProngVect::const_iterator itProng = primProngs.begin(); itProng != primProngs.end(); ++itProng ) {
        if( (*itProng)->MinosTrack() || (*itProng)->MinosStub()) isPlausible = muonIsPlausible( (*itProng));
        if(isPlausible) return true;
      }
    }//end loop over vertices


    return false;
  }

  //------------------------------------------------------------------------------------------------------------
  //There are 2 methods available here. One uses muonIsPlausible( ), the other is a more general method.
  //The general method could be used for the muon or proton analysis, if using this tool.
  //------------------------------------------------------------------------------------------------------------
  if( m_useMuonIsPlausible ) {

    bool foundDeProng = false;
    ProngVect primaryProngs = physEvent->primaryProngs();
    if( primaryProngs.size() ) {
      for( ProngVect::iterator itProng = primaryProngs.begin(); itProng != primaryProngs.end(); ++itProng ) {
        if( (*itProng)->hasIntData("look_for_plausibility") ) {
          if( (*itProng)->getIntData("look_for_plausibility") == 1 ) {
            foundDeProng = true;
            bool is_muon_plausible = muonIsPlausible( (*itProng) );
            if( !is_muon_plausible ) {
              debug() << "This MC muon is not plausible because it has a MC energy fraction below threshold !" << endmsg;
              return false;
            } else {
              debug() << "This MC muon is plausible !" << endmsg;
              counter("muon_plausible")+=is_muon_plausible;
              return true;
            }
          }
        } else {
          debug() << "This primary Prong does not have look_for_plausibility set and is probably not a muon Prong !" << endmsg;
        }
      }
      if( !foundDeProng ) throw MinervaException( "truthIsPlausible: Did not find muon Prong amongst primary Prongs, even though you had marked it earlier !" );

    } else {
      debug() << "There are no Primary Prongs in this Physics Event !" << endmsg;
      return false;
    }

  } else {

    //------------------------------------------------------------------------------------------
    // Obtain a vector of Primary Prongs from the event
    // Should the other prongs (non-primary) be included here too ?
    // Find the longest track from all the prongs (could be muon, proton, etc.)
    // This relies on the assumption that rock muons are filtered out early on and are
    // not part of the physics event.
    //------------------------------------------------------------------------------------------
    SmartRef<Minerva::Track> deLongestTrack;
    bool foundaTrack = findLongestMinervaTrack( physEvent, deLongestTrack );
    verbose() << "Inside truthIsPlausible(): Value of foundaTrack " << foundaTrack << endmsg;
    if( foundaTrack && deLongestTrack->nNodes() ) {
      debug() <<"Inside truthIsPlausible(): Chi2 of deLongestTrack from findLongestMinervaTrack( ) is " << deLongestTrack->chi2PerDoF() <<" Nodes: "<< deLongestTrack->nNodes() << endmsg;
    } else {
      debug() <<"Did not find any tracks in this PhysicsEvent !" << endmsg;
      return false;
    }

    //----------------------------------------------------------------------------------------------------------------
    // Find all the Prongs that hold this track and verify if they all carry MC identities
    // If not, then the track probably comes from a data overlay and not MC. Hence it is not your signal event !
    //----------------------------------------------------------------------------------------------------------------
    SmartRefVector<Minerva::Prong> prongsFromLongestTrack;
    StatusCode sc = m_objectAssociator->getProngs_containingMinervaTrack( prongsFromLongestTrack, deLongestTrack );
    if( sc.isFailure() ) {
      debug() << "The method did not work correctly!" << endmsg;
      return false;
    }

    if( prongsFromLongestTrack.size() ) {
      for( SmartRefVector<Minerva::Prong>::iterator itProng = prongsFromLongestTrack.begin(); itProng != prongsFromLongestTrack.end(); ++itProng ) {
        if( !prongIsMC( *itProng ) ) return false;
      }
    } else {
      throw MinervaException( "The longest track does not belong to a Prong! Please investigate!" );
    }

  } //End of condition m_useMuonIsPlausible

  info() << "Exiting truthIsPlausible( ) now ..... " << endmsg;
  return true;

} //End of truthIsPlausible( )


//=============================================================================
// Tag truth events as signal
//=============================================================================
StatusCode MasterAnaDev::tagTruth(Minerva::GenMinInteraction* truth) const {

  debug()<<"MasterAnaDev::tagTruth"<<endmsg;
  fillInelasticReweightBranches(truth);
  fillGenieWeightBranches( truth );
  // for HadronReweightTool
  fillHadronReweightBranches( truth );
  int nChargePi=0;


  const SmartRefVector<Minerva::TG4Trajectory> pri_trajectories = truth->trajectories();
  SmartRefVector<Minerva::TG4Trajectory>::const_iterator it_mcpart;

  int nPiPlus = 0, nChargedK = 0, nLambda = 0, nK0 = 0, nSigma = 0, nPiMinus = 0, nPi0 = 0;
  std::vector<int> otherPDG, t_pi_charge;
  std::vector<double> t_pi_px, t_pi_py, t_pi_pz, t_pi_E, t_pi_theta;

  for (int iter = 0; iter < 20; iter++) {
    t_pi_px.push_back(-1.0); t_pi_py.push_back(-1.0); t_pi_pz.push_back(-1.0); t_pi_E.push_back(-1.0);
    t_pi_theta.push_back(-9.0);
    t_pi_charge.push_back(0);
  }

  for (it_mcpart = pri_trajectories.begin(); it_mcpart != pri_trajectories.end(); ++it_mcpart){
    if ( (*it_mcpart)->GetParentId() != 0 ) continue;

    Gaudi::LorentzVector tmp_4p = (*it_mcpart)->GetInitialMomentum();
    int partPDG = (*it_mcpart)->GetPDGCode();

    if ( abs(partPDG) == 211 && nChargePi < 20){
      t_pi_E[nChargePi]  = (*it_mcpart)->GetInitialMomentum().E();
      t_pi_px[nChargePi] = (*it_mcpart)->GetInitialMomentum().px();
      t_pi_py[nChargePi] = (*it_mcpart)->GetInitialMomentum().py();
      t_pi_pz[nChargePi] = (*it_mcpart)->GetInitialMomentum().pz();
      t_pi_theta[nChargePi] = m_minCoordSysTool->thetaWRTBeam(tmp_4p);
      t_pi_charge[nChargePi] = partPDG>0 ? 1 : -1;

      nChargePi++;
      if ( partPDG>0 ) {
        nPiPlus++;
      }
      else {
        nPiMinus++;
      }

    }

    else if ( partPDG == 111 && nPi0 < 20 ) {
      nPi0++;
    }

    else if ( abs(partPDG) == 321 && nChargedK < 20 ){
      nChargedK++;
    }
    else if ( (abs(partPDG) == 130 || abs(partPDG) == 310 || abs(partPDG) == 311) && nK0 < 20){
      nK0++;
    }
    else if ( abs(partPDG) == 3222 && nSigma < 20){
      nSigma++;
    }
    else if (abs(partPDG) == 3122 && nLambda < 20){
      nLambda++;
    }

  }

  Gaudi::LorentzVector t_mu4p = truth->PrimFSLepton();
  double t_mu_px = t_mu4p.px(), t_mu_py = t_mu4p.py(), t_mu_pz = t_mu4p.pz(), t_mu_E= t_mu4p.E();

  truth->setIntData("N_pip",        nPiPlus);
  truth->setIntData("N_chargedK",   nChargedK);
  truth->setIntData("N_lambda",     nLambda);
  truth->setIntData("N_K0",         nK0);
  truth->setIntData("N_sigma",      nSigma);
  truth->setIntData("N_pim",        nPiMinus);
  truth->setIntData("N_pi0",        nPi0);

  truth->setDoubleData("muon_px", t_mu_px);
  truth->setDoubleData("muon_py", t_mu_py);
  truth->setDoubleData("muon_pz", t_mu_pz);
  truth->setDoubleData("muon_E",  t_mu_E);
  truth->setDoubleData( "muon_theta",  m_minCoordSysTool->thetaWRTBeam(  truth->PrimFSLepton() ));

  truth->setContainerDoubleData("pi_px", t_pi_px);
  truth->setContainerDoubleData("pi_py", t_pi_py);
  truth->setContainerDoubleData("pi_pz", t_pi_pz);
  truth->setContainerDoubleData("pi_E",  t_pi_E);
  truth->setContainerDoubleData("pi_theta_wrtbeam",  t_pi_theta);
  truth->setContainerIntData("pi_charge", t_pi_charge);

  //True target ID based on NukeCCInclusive.  I think NukeCCInclusive macros use
  //these branches for studying signal definitions.
  //TODO: Is this what standard reconstruction needs too?
  int targetCode = 0;
  getExactTarget(truth, targetCode);

  //Unpack targetCode into a target ID (position in z) and struck nucleus (target Z).
  const int targetID = (targetCode - targetCode % 1000) / 1000;
  const int targetZ = targetCode % 1000;

  truth->setIntData("target_code", targetCode); //TODO: Do we need this?  My notes on NukeCCInclusive say no.
  truth->setIntData("targetID", targetID);
  truth->setIntData("targetZ", targetZ);

  vector<int> intVec;
  vector<double> doubleVec;
  getReferenceZ( intVec, truth->Vtx().x() , truth->Vtx().y() );
  truth->setContainerIntData( "ref_targZ", intVec );
  getReferenceDistToDivision( doubleVec, truth->Vtx().x(), truth->Vtx().y()  );
  truth->setContainerDoubleData( "ref_dist_to_division", doubleVec );
  getDistToTarget( doubleVec, truth->Vtx().z() );
  truth->setContainerDoubleData( "ref_dist_to_target", doubleVec );

  if( 0 < targetID )
    {
      truth->setDoubleData("target_zDist", m_nukeTool->getNuclearTarget(targetID)->distanceToTarget(truth->Vtx().z()));
      truth->setDoubleData("target_dist_to_division", m_nukeTool->getNuclearTarget(targetID)->distanceToDivision(truth->Vtx().x(), truth->Vtx().y()));
    }


  //! Cut 3: Is the vertex in the fiducial/analyzable region
  bool isFiducialXY   = FiducialPointTool->isFiducial( truth, m_signalApothem, m_usFiducialModID, m_dsFiducialModID );
  bool isAnalyzableXY = FiducialPointTool->isFiducial(truth,m_analyzeApothem,m_usAnalyzableModID,m_dsAnalyzableModID);
  truth->setIntData("is_fiducial_area", isFiducialXY);
  truth->setIntData("is_analyzable_area", isAnalyzableXY);

  //store the module/plane of the vertex (-999 for events in passive material)
  int vtxMod = -999, vtxPlane = -999;
  bool inPlane = getVtxPlane( truth->Vtx().z(), vtxMod, vtxPlane );
  truth->filtertaglist()->setOrAddFilterTag( "in_plane", inPlane );
  truth->setIntData( "vtx_module", vtxMod );
  truth->setIntData( "vtx_plane", vtxPlane );

  //TODO: NukeCCInclusive team, there's lots of other code about FS particles that could go here.

  //--------------------------------
  // Is it inside fiducial vol ?
  //--------------------------------
   bool is_fiducial = FiducialPointTool->isFiducial( truth, m_signalApothem, m_signalUpstreamZ, m_signalDownstreamZ );
   if( truth ) {
     truth->filtertaglist()->setOrAddFilterTag( "is_fiducial", is_fiducial );
    }

  //!@todo Tag true-CCQE events and fill variables for efficiency corrections...

   const std::vector<int> fsPDG  = truth->fSpdg();
   const std::vector<double> fsE = truth->fsParticlesE();
	
   double protonE = 0.0;
   double neutronE = 0.0;
   double pionE = 0.0;
   double pizeroE = 0.0;
   double leptonE = 0.0;
   double gammaE = 0.0;
   double otherE = 0.0;
   double missingE = 0.0;
	
   int protonN = 0;
   int neutronN = 0;
   int pionN = 0;
   int pizeroN = 0;
   int leptonN = 0;
   int gammaN = 0;
   int otherN = 0;
	
   for( unsigned int iPart = 0; iPart != fsPDG.size(); ++iPart )
     {
       const double E = fsE[iPart];
       const int pdg  = fsPDG[iPart];
	    
       if(pdg==2212){
         protonE += E - MinervaUnits::M_p;
         if(E - MinervaUnits::M_p > 5.0)protonN++;
       } else if(pdg==2112){
         neutronE += E - MinervaUnits::M_n;
         if(E - MinervaUnits::M_n > 5.0)neutronN++;
       } else if(abs(pdg)==211){
         pionE += E;
         pionN++;
       } else if(pdg==111){
         pizeroE += E;
         pizeroN++;
       } else if(abs(pdg)==13 || abs(pdg)==11){
         leptonE += E;
         leptonN++;
       } else if(fabs(pdg)==22){
         gammaE += E;
         gammaN++;
       } else if(pdg>=2000000000){
          missingE += E;
       } else if(pdg>1000000000){
         // do nothing, assume negligible energy for nucleons
          // save me the trouble of asking what that nucleon's mass was.
       } else if(pdg>=2000){
         otherE += E - MinervaUnits::M_p;
         otherN++;
         //primarily strange baryons 3122s
       } else if(pdg<=-2000){
         otherE += E + MinervaUnits::M_p;
         otherN++;
         //primarily anti-baryons -2112 and -2212s,
       } else {
         otherE += E;
         otherN++;
         //primarily kaons and eta mesons.
       }
	  } //end loop over fs particles
   truth->setDoubleData( "visibleE", protonE + pionE - pionN*MinervaUnits::M_pion + pizeroE + gammaE + otherE );

  return StatusCode::SUCCESS;
}

void MasterAnaDev::tagProtonProngTruth( Minerva::PhysicsEvent *event, SmartRef<Minerva::Prong> protonProng ) const {

  StatusCode sc;
  if ( !haveNeutrinoMC() ) return;

  const Minerva::TG4Trajectory* prong_traj = NULL;
  double prong_traj_fractionE = -1.0, prong_traj_otherE = -1.0;
  int prong_traj_PDG = -1, prong_traj_ID;
  std::vector<double> prong_traj_4p;
  std::vector<double> prong_traj_final_tpos;

  //------------------------------------------------------------------
  //! Let's tag the highest energy primary trajectory first (proton candidate)
  //------------------------------------------------------------------
  sc = TruthMatcher->getPrimaryTG4Trajectory(protonProng, prong_traj, prong_traj_fractionE, prong_traj_otherE);
  if (sc == StatusCode::FAILURE){
    debug()<<"Truth Matching Failed while getting the proton candidate trajectory"<<endmsg;
    return;
  }
  //! Set proton candidate PDG, four-momentum and position
  prong_traj_PDG = prong_traj->GetPDGCode();
  prong_traj_ID = prong_traj->GetTrackId();

  prong_traj_4p.push_back( prong_traj->GetInitialMomentum().E() );
  prong_traj_4p.push_back( prong_traj->GetInitialMomentum().px() );
  prong_traj_4p.push_back( prong_traj->GetInitialMomentum().py() );
  prong_traj_4p.push_back( prong_traj->GetInitialMomentum().pz() );

  prong_traj_final_tpos.push_back( prong_traj->GetFinalPosition().T() );
  prong_traj_final_tpos.push_back( prong_traj->GetFinalPosition().x() );
  prong_traj_final_tpos.push_back( prong_traj->GetFinalPosition().y() );
  prong_traj_final_tpos.push_back( prong_traj->GetFinalPosition().z() );

  event->setIntData("proton_prong_PDG", prong_traj_PDG);
  event->setIntData("proton_prong_traj_ID", prong_traj_ID);
  event->setContainerDoubleData("proton_prong_4p", prong_traj_4p);
  event->setContainerDoubleData("proton_prong_tpos", prong_traj_final_tpos);

}


void MasterAnaDev::tagSecondaryProtonProngTruth( Minerva::PhysicsEvent *event, Minerva::ProngVect protonProngs ) const {

  StatusCode sc;
  if ( !haveNeutrinoMC() ) return;

  std::vector<double> fractionE;
  std::vector<double> fractionOtherE;
  std::vector<int> trajpdg;
  std::vector<int> trajID;
  std::vector<double> traj_4p_E;
  std::vector<double> traj_4p_px;
  std::vector<double> traj_4p_py;
  std::vector<double> traj_4p_pz;
  std::vector<double> traj_tpos_x;
  std::vector<double> traj_tpos_y;
  std::vector<double> traj_tpos_z;
  std::vector<double> traj_tpos_t;

  for( SmartRefVector<Minerva::Prong>::iterator itProng = protonProngs.begin(); itProng != protonProngs.end(); ++itProng )  {
    const Minerva::TG4Trajectory* prong_traj = NULL;
    double prong_traj_fractionE = -1.0, prong_traj_otherE = -1.0;
    //int prong_traj_PDG = -1, prong_traj_ID;
    //------------------------------------------------------------------
    //! Let's tag the highest energy primary trajectory first (proton candidate)
    //------------------------------------------------------------------
    sc = TruthMatcher->getPrimaryTG4Trajectory(*itProng, prong_traj, prong_traj_fractionE, prong_traj_otherE);
    if (sc == StatusCode::FAILURE){
      debug()<<"Truth Matching Failed while getting the proton candidate trajectory"<<endmsg;
      return;
    }
    trajpdg.push_back(prong_traj->GetPDGCode());
    trajID.push_back(prong_traj->GetTrackId());

    traj_4p_E.push_back(prong_traj->GetInitialMomentum().E());
    traj_4p_px.push_back(prong_traj->GetInitialMomentum().px());
    traj_4p_py.push_back(prong_traj->GetInitialMomentum().py());
    traj_4p_pz.push_back(prong_traj->GetInitialMomentum().pz());

    traj_tpos_x.push_back(prong_traj->GetFinalPosition().x());
    traj_tpos_y.push_back(prong_traj->GetFinalPosition().y());
    traj_tpos_z.push_back(prong_traj->GetFinalPosition().z());
    traj_tpos_t.push_back(prong_traj->GetFinalPosition().T());

  }

  event->setContainerIntData("sec_protons_prong_PDG",trajpdg);
  event->setContainerIntData("sec_protons_prong_traj_ID",trajID);
  event->setContainerDoubleData("seco_protons_prong_4p_E",traj_4p_E);
  event->setContainerDoubleData("seco_protons_prong_4p_px",traj_4p_px);
  event->setContainerDoubleData("seco_protons_prong_4p_py",traj_4p_py);
  event->setContainerDoubleData("seco_protons_prong_4p_pz",traj_4p_pz);
  event->setContainerDoubleData("proton_prong_tpos_x",traj_tpos_x);
  event->setContainerDoubleData("proton_prong_tpos_y",traj_tpos_y);
  event->setContainerDoubleData("proton_prong_tpos_z",traj_tpos_z);
  event->setContainerDoubleData("proton_prong_tpos_t",traj_tpos_t);


}

std::vector<double> MasterAnaDev::tagMichelProngTruth( Minerva::Prong michelProng ) const {

  StatusCode sc;

  //  Minerva::Prong *tmp_prong = &michelProng;
  std::vector<const Minerva::IDCluster*> idclusters;
  SmartRefVector<Minerva::IDCluster> clustersFromProng = michelProng.getAllIDClusters();
  idclusters.insert(idclusters.begin(), clustersFromProng.begin(), clustersFromProng.end());
  const Minerva::TG4Trajectory* prong_traj = NULL;
  const Minerva::TG4Trajectory* prong_parent_traj = NULL;
  double prong_traj_fractionE = -1.0, prong_traj_otherE=-1.0;
  //  double prong_parent_traj_fractionE = -1.0, prong_parent_traj_otherE=-1.0;
  double prong_traj_PDG = -1;
  double prong_parent_traj_PDG = -1;
  std::vector<double> bestPart;
  bestPart.push_back(prong_traj_PDG);
  bestPart.push_back(prong_parent_traj_PDG);
  bestPart.push_back(prong_traj_fractionE);
  if ( !haveNeutrinoMC() ) return bestPart;
  //------------------------------------------------------------------
  //! Let's tag the highest energy primary trajectory first (michel candidate)
  //------------------------------------------------------------------
  sc = TruthMatcher->getTG4Trajectory(idclusters, prong_traj, prong_traj_fractionE, prong_traj_otherE);
  if (sc == StatusCode::FAILURE){
    debug()<<"Truth Matching Failed while getting the michel candidate trajectory"<<endmsg;
    return bestPart;
  }
  prong_parent_traj = TruthMatcher->getPrimaryTG4Trajectory(prong_traj,false);
  //! Set michel candidate PDG, four-momentum and position
  prong_traj_PDG = prong_traj->GetPDGCode();
  prong_parent_traj_PDG = prong_parent_traj->GetPDGCode();
  bestPart[0] = prong_traj_PDG;
  bestPart[1] = prong_parent_traj_PDG;
  bestPart[2] = prong_traj_fractionE;
  debug() << "This is the matched information(pdg,pdg,fraction): (" <<bestPart[0]<<","<<bestPart[1] << "," << bestPart[2] << ")" << endmsg;


  return bestPart;
}

std::vector<double> MasterAnaDev::tagBlobTruth( Minerva::IDBlob blob ) const {

  StatusCode sc;

  //  Minerva::Prong *tmp_prong = &blobProng;
  std::vector<const Minerva::IDDigit*> iddigits;
  SmartRefVector<Minerva::IDDigit> digitsFromBlob = blob.getAllDigits();
  iddigits.insert(iddigits.begin(), digitsFromBlob.begin(), digitsFromBlob.end());
  const Minerva::TG4Trajectory* blob_traj = NULL;
  const Minerva::TG4Trajectory* blob_parent_traj = NULL;
  double blob_traj_fractionE = -1.0, blob_traj_otherE=-1.0;

  double blob_traj_PDG = -1;
  double blob_parent_traj_PDG = -1;
  std::vector<double> bestPart;
  bestPart.push_back(blob_traj_PDG);
  bestPart.push_back(blob_parent_traj_PDG);
  bestPart.push_back(blob_traj_fractionE);
  if ( !haveNeutrinoMC() ) return bestPart;
  //------------------------------------------------------------------
  //! Let's tag the highest energy primary trajectory first (michel candidate)
  //------------------------------------------------------------------
  sc = TruthMatcher->getTG4Trajectory(iddigits, blob_traj, blob_traj_fractionE, blob_traj_otherE);
  if (sc == StatusCode::FAILURE){
    debug()<<"Truth Matching Failed while getting the blob candidate trajectory"<<endmsg;
    return bestPart;
  }
  blob_parent_traj = TruthMatcher->getPrimaryTG4Trajectory(blob_traj,false);
  //! Set michel candidate PDG, four-momentum and position
  blob_traj_PDG = blob_traj->GetPDGCode();
  blob_parent_traj_PDG = blob_parent_traj->GetPDGCode();
  bestPart[0] = blob_traj_PDG;
  bestPart[1] = blob_parent_traj_PDG;
  bestPart[2] = blob_traj_fractionE;
  debug() << "This is the matched information(pdg,pdg,fraction): (" <<bestPart[0]<<","<<bestPart[1] << "," << bestPart[2] << ")" << endmsg;


  return bestPart;
}

std::vector<double> MasterAnaDev::tagDigitsTruth( SmartRefVector<Minerva::IDDigit> digits ) const {

  StatusCode sc;

  //  Minerva::Prong *tmp_prong = &blobProng;
  std::vector<const Minerva::IDDigit*> iddigits;
  iddigits.insert(iddigits.begin(), digits.begin(), digits.end());
  const Minerva::TG4Trajectory* digits_traj = NULL;
  const Minerva::TG4Trajectory* digits_parent_traj = NULL;
  double digits_traj_fractionE = -1.0, digits_traj_otherE=-1.0;

  double digits_traj_PDG = -1;
  double digits_parent_traj_PDG = -1;
  std::vector<double> bestPart;
  bestPart.push_back(digits_traj_PDG);
  bestPart.push_back(digits_parent_traj_PDG);
  bestPart.push_back(digits_traj_fractionE);
  if ( !haveNeutrinoMC() ) return bestPart;
  //------------------------------------------------------------------
  //! Let's tag the highest energy primary trajectory first (michel candidate)
  //------------------------------------------------------------------
  sc = TruthMatcher->getTG4Trajectory(iddigits, digits_traj, digits_traj_fractionE, digits_traj_otherE);
  if (sc == StatusCode::FAILURE){
    debug()<<"Truth Matching Failed while getting the digits candidate trajectory"<<endmsg;
    return bestPart;
  }
  digits_parent_traj = TruthMatcher->getPrimaryTG4Trajectory(digits_traj,false);
  //! Set michel candidate PDG, four-momentum and position
  digits_traj_PDG = digits_traj->GetPDGCode();
  digits_parent_traj_PDG = digits_parent_traj->GetPDGCode();
  bestPart[0] = digits_traj_PDG;
  bestPart[1] = digits_parent_traj_PDG;
  bestPart[2] = digits_traj_fractionE;
  debug() << "This is the matched information(pdg,pdg,fraction): (" <<bestPart[0]<<","<<bestPart[1] << "," << bestPart[2] << ")" << endmsg;


  return bestPart;
}


//==================================================
// createParticles
//==================================================
bool MasterAnaDev::createParticles( Minerva::PhysicsEvent* event, Minerva::ProngVect& hadronProngs ) const
{
  bool status = true;

  ProngVect primaryProngs = event->primaryProngs();
  std::vector<Minerva::Particle::ID> particleHypotheses;
  particleHypotheses.push_back(Minerva::Particle::Pion);
  particleHypotheses.push_back(Minerva::Particle::Proton);

  for (ProngVect::iterator itProngs = primaryProngs.begin(); itProngs != primaryProngs.end(); ++itProngs) {

    // only non-muon prongs...
    if( (*itProngs)->filtertaglist()->filterTagExists( AnaFilterTags::PrimaryMuon() ) ) continue;

    hadronProngs.push_back(*itProngs);
    (*itProngs)->filtertaglist()->setOrAddFilterTag( m_primaryHadron, true );

    //Let's color hadron prongs
    m_hitTaggerTool->applyColorTag( *itProngs, m_hadronProngsColor );

    debug()<<"Make particles"<<endmsg;
    IParticleMakerTool::NameAliasListType toolsToUse;
    toolsToUse.push_back( std::make_pair("dEdXTool", "dEdXTool") );  // name, alias pair

    StatusCode sc_dEdX = m_particleMaker->makeParticles( *itProngs, particleHypotheses, toolsToUse );
    if (!sc_dEdX) {
      debug()<<"Could not make particles in prong!"<<endmsg;
    }

    //make an endpoint vertex blob for the hadron prong - currently this depends on the prong order!!!
    debug()<<"Make endpoint blob"<<endmsg;
    SmartRef<Minerva::Vertex> endpoint_vtx;
    m_objectAssociator->getVertex_fromTrackBack( endpoint_vtx, (*itProngs)->minervaTracks().back() );
    if (!endpoint_vtx) {
      warning()<<"Could not find a back vertex for this prong!"<<endmsg;
      continue;
    }

    //search for endpoint michels
    debug()<<"Search for michel"<<endmsg;
    Minerva::Prong michelProng;
    bool foundMichel = m_michelTool->findMichel( endpoint_vtx, michelProng );
    SmartRef<Minerva::Prong> smart_michelProng(&michelProng);
    double michel_mc_frac = -1.0;
    muonIsPlausible( smart_michelProng, michel_mc_frac);
    if (foundMichel) {
      (*itProngs)->setIntData("hasEndpointMichel",1);
      (*itProngs)->setIntData("endMichel_category",michelProng.getIntData("category"));
      (*itProngs)->setIntData("endMichel_ndigits",michelProng.getIntData("ndigits"));
      (*itProngs)->setDoubleData("endMichel_energy",michelProng.getDoubleData("energy"));
      (*itProngs)->setDoubleData("endMichel_slice_energy",michelProng.getDoubleData("slice_energy"));
    }
  }//end loop over primary prongs

  return status;
}

void MasterAnaDev::tagMichelElectrons(Minerva::PhysicsEvent* event) const {
  //-- Calculate and Store Michel electron data (Matching will occur in later macros)
  //-- Note, Michel index != hadron index
  verbose() << "Enter MasterAnaDev::tagMichelElectrons" << endmsg;
  // initialize variables
  int michel_count;
  const int MAX_TRACKS = 10;

    std::vector<double> michel_vtx_x(MAX_TRACKS,-9999);
    std::vector<double> michel_vtx_y(MAX_TRACKS,-9999);
    std::vector<double> michel_vtx_z(MAX_TRACKS,-9999);

    std::vector<int> matched_michel_code(MAX_TRACKS,-1);
    std::vector<int> matched_michel_idx(MAX_TRACKS,-1);
    std::vector<int> matched_michel_avg_idx(MAX_TRACKS,-1);
    std::vector<int> matched_michel_end_idx(MAX_TRACKS,-1);
    std::vector<int> matched_michel_ov_idx(MAX_TRACKS,-1);

    std::vector<double> matched_michel_cal_energy(MAX_TRACKS,-1);
    std::vector<double> matched_michel_vis_energy(MAX_TRACKS,-1);
    std::vector<double> matched_michel_time(MAX_TRACKS,-1);
    std::vector<double> matched_michel_slice(MAX_TRACKS,-1);
    std::vector<int> matched_michel_views(MAX_TRACKS,-1);

    std::vector<double> matched_michel_avg_dist(MAX_TRACKS,-1);
    std::vector<double> matched_michel_end_dist(MAX_TRACKS,-1);
    std::vector<double> matched_michel_ov_dist(MAX_TRACKS,-1);

    std::vector<double> matched_michel_position_x(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_position_y(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_position_u(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_position_v(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_position_z(MAX_TRACKS,-1);

    std::vector<int>    matched_michel_fitPass(MAX_TRACKS,-1);
    std::vector<double> matched_michel_fit_residual(MAX_TRACKS,-1);

    std::vector<double> matched_michel_fit_AX(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_fit_AY(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_fit_AU(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_fit_AV(MAX_TRACKS,-9999);

    std::vector<double> matched_michel_end_X1(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_end_Y1(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_end_U1(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_end_V1(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_end_Z1(MAX_TRACKS,-9999);

    std::vector<double> matched_michel_end_X2(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_end_Y2(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_end_U2(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_end_V2(MAX_TRACKS,-9999);
    std::vector<double> matched_michel_end_Z2(MAX_TRACKS,-9999);

    // new 2018.11.21 modules corresponding to michel position
    std::vector<int> matched_michel_avg_module( MAX_TRACKS,  -9999);
    std::vector<int> matched_michel_end_module1( MAX_TRACKS, -9999);
    std::vector<int> matched_michel_end_module2( MAX_TRACKS, -9999);

  m_improvedmichelTool->Reset();

  // 1.) Collect michel candidates -- make sure the hits look like michels
  m_improvedmichelTool->FindMichel(event->earliestSliceNumber());
  m_improvedmichelTool->ApplyQualityCuts();
  michel_count = m_improvedmichelTool->GetMichelCount(); //How many found?

  // Initialize more variables with michel_count
    std::vector<double> has_michel_cal_energy(michel_count,-1);
    std::vector<double> has_michel_vis_energy(michel_count,-1);
    std::vector<double> has_michel_time(michel_count,-1);
    std::vector<int>    has_michel_slice( michel_count, -1 );
    std::vector<double> has_michel_position_x(michel_count,-9999);
    std::vector<double> has_michel_position_y(michel_count,-9999);
    std::vector<double> has_michel_position_u(michel_count,-9999);
    std::vector<double> has_michel_position_v(michel_count,-9999);
    std::vector<double> has_michel_position_z(michel_count,-1);
    std::vector<int> has_michel_nClusters(michel_count,-1);
    std::vector<int> has_michel_views(michel_count,-1);

    std::vector<int> has_michel_primary_particle_trackID(michel_count,-9999);
    std::vector<int> has_michel_primary_particle_PDG(michel_count,0);
    std::vector<double> has_michel_primary_particle_p(michel_count,-9999);
    std::vector<double> has_michel_primary_particle_px(michel_count,-9999);
    std::vector<double> has_michel_primary_particle_py(michel_count,-9999);
    std::vector<double> has_michel_primary_particle_pz(michel_count,-9999);
    std::vector<double> has_michel_primary_particle_E(michel_count,-9999);

    std::vector<int> true_michel_pdg(michel_count, 0);
    std::vector<double> true_michel_x1(michel_count,-9999);
    std::vector<double> true_michel_y1(michel_count,-9999);
    std::vector<double> true_michel_z1(michel_count,-9999);
    std::vector<double> true_michel_x2(michel_count,-9999);
    std::vector<double> true_michel_y2(michel_count,-9999);
    std::vector<double> true_michel_z2(michel_count,-9999);
    std::vector<int> michel_from_decay(michel_count,0);
    std::vector<int> michel_parent_PDG(michel_count,0);
    std::vector<int> michel_parent_trackID(michel_count,-9999);

    std::vector<int>    has_michel_fitPass(michel_count,-1);
    std::vector<double> has_michel_fit_residual(michel_count,-1);
    std::vector<double> has_michel_fit_AX(michel_count,-9999);
    std::vector<double> has_michel_fit_AY(michel_count,-9999);
    std::vector<double> has_michel_fit_AU(michel_count,-9999);
    std::vector<double> has_michel_fit_AV(michel_count,-9999);
    std::vector<double> has_michel_end_X1(michel_count,-9999);
    std::vector<double> has_michel_end_Y1(michel_count,-9999);
    std::vector<double> has_michel_end_U1(michel_count,-9999);
    std::vector<double> has_michel_end_V1(michel_count,-9999);
    std::vector<double> has_michel_end_Z1(michel_count,-9999);
    std::vector<double> has_michel_end_X2(michel_count,-9999);
    std::vector<double> has_michel_end_Y2(michel_count,-9999);
    std::vector<double> has_michel_end_U2(michel_count,-9999);
    std::vector<double> has_michel_end_V2(michel_count,-9999);
    std::vector<double> has_michel_end_Z2(michel_count,-9999);

    // new 2018.11.21 modules corresponding to michel position
    std::vector<int> has_michel_avg_module( michel_count,  -9999);
    std::vector<int> has_michel_end_module1( michel_count, -9999);
    std::vector<int> has_michel_end_module2( michel_count, -9999);

  // 2.) Now, collect the vertices that we'll try to match to michels
  SmartRefVector<Minerva::Vertex> vertices;
  vertices.clear();
    // a.) The interaction vertex is special. It will always be the 0th index of
    // our various container branches
    SmartRef<Minerva::Vertex> int_vtx = event->interactionVertex();
    if (!int_vtx) error() << "No Interaction Vertex Found!" << endmsg;
    else{
      info()<<"Got interaction vertex for Michel matching" <<endmsg;
      vertices.push_back(int_vtx);
    }

    // b.) Next get the hadronic prongs so we can get their vertices
    Minerva::ProngVect primary_prongs = event->primaryProngs();
    Minerva::ProngVect pionProngs;
    pionProngs.clear();
    for(unsigned int prong = 0; prong <primary_prongs.size(); ++prong) {
      if( primary_prongs[prong]->filtertaglist()->filterTagExists("PrimaryHadron")
          &&
          !primary_prongs[prong]->filtertaglist()->filterTagExists("PrimaryMuon") ) {
        pionProngs.push_back( primary_prongs[prong] );
      }
    }
    info()<<"Number of Primary Prongs: "<<pionProngs.size()<<endmsg;

    //c.) ...And get the prongs' vertices
    for(unsigned int i = 0; i < pionProngs.size(); i++) {
      Gaudi::XYZPoint endpoint = (pionProngs[i]->minervaTracks().back())->lastState().position();
      SmartRef<Minerva::Vertex> end_vtx = new Minerva::Vertex(endpoint);
      vertices.push_back(end_vtx);
      if (vertices.size() == MAX_TRACKS) break;
    }

    info()<<"Vertex sizes: "<<vertices.size()<<endmsg;

    // d.) Save the vertex positions
    for (unsigned int i = 0; i < vertices.size(); i++){
      michel_vtx_x[i] = vertices[i]->position().x() ;
      michel_vtx_y[i] = vertices[i]->position().y() ;
      michel_vtx_z[i] = vertices[i]->position().z() ;
    }

  // 3.) Calculate the distances between each vertex and the michels
  m_improvedmichelTool->OptimalDistance( vertices, 500, 500 );

  // 4.) Check each vertex for a match
  for(unsigned int v = 0; v < vertices.size(); ++v) {
    bool match = m_improvedmichelTool->OptimalMatch(v, 500, 500);
    if (match){
      matched_michel_code[v]    = m_improvedmichelTool->GetMatchedMichelCode();
      matched_michel_idx[v]     = m_improvedmichelTool->GetMatchedMichelIndex();
      matched_michel_avg_idx[v] = m_improvedmichelTool->GetMatchedAvgNofMichelIndex();
      matched_michel_end_idx[v] = m_improvedmichelTool->GetMatchedEndMichelIndex();
      matched_michel_ov_idx[v]  = m_improvedmichelTool->GetMatchedOneViewMichelIndex();

      // new 2018.11.21
      // time slice of vertex v's michel.
      SmartRefVector<Minerva::IDCluster> michel_clusters = m_improvedmichelTool->GetMichelClusters(matched_michel_idx[v]);
      Minerva::TimeSlice* time_slice = getTimeSlice( michel_clusters[0] );
      int slice = time_slice ? time_slice->sliceNumber() : -1;
      matched_michel_slice[v] = slice;

      matched_michel_time[v]       = m_improvedmichelTool->GetMichelTime();
      matched_michel_cal_energy[v] = m_improvedmichelTool->GetMichelCalorimetricEnergy();
      matched_michel_vis_energy[v] = m_improvedmichelTool->GetMichelVisibleEnergy();
      matched_michel_views[v]      = m_improvedmichelTool->GetMichelViewCode();

      matched_michel_fitPass[v]      = m_improvedmichelTool->GetMichelFitPass();
      matched_michel_fit_residual[v] = m_improvedmichelTool->GetMichelFitResidual();

      matched_michel_avg_dist[v] = m_improvedmichelTool->GetMatchedAvgNofMichelDist();
      matched_michel_end_dist[v] = m_improvedmichelTool->GetMatchedEndMichelDist();
      matched_michel_ov_dist[v]  = m_improvedmichelTool->GetMatchedOneViewMichelDist();

      // new 2018.11.21
      // get the module corresponding to vertex v's michel's avg position      
      Gaudi::XYZPoint hitXYZMAvg( matched_michel_position_x[v], matched_michel_position_y[v], matched_michel_position_z[v] );
      const Minerva::DePlane* trkrHitPlaneAvg   = m_InnerDetector->getDePlane(hitXYZMAvg);
      if( !trkrHitPlaneAvg )  // Cases where we're not in a tracker module (and I'm too lazy to find the regular module
      {
        matched_michel_avg_module[v] = -6;
      }
      else
      {
        const Minerva::PlaneID& trkrHitPlaneIdAvg = trkrHitPlaneAvg->getPlaneID();
        matched_michel_avg_module[v] = trkrHitPlaneIdAvg.module();
      }

      matched_michel_position_x[v] = m_improvedmichelTool->GetMichelPositionX();
      matched_michel_position_y[v] = m_improvedmichelTool->GetMichelPositionY();
      matched_michel_position_u[v] = m_improvedmichelTool->GetMichelPositionU();
      matched_michel_position_v[v] = m_improvedmichelTool->GetMichelPositionV();
      matched_michel_position_z[v] = m_improvedmichelTool->GetMichelPositionZ();

      matched_michel_fit_AX[v] = m_improvedmichelTool->GetMichelFitAX();
      matched_michel_fit_AY[v] = m_improvedmichelTool->GetMichelFitAY();
      matched_michel_fit_AU[v] = m_improvedmichelTool->GetMichelFitAU();
      matched_michel_fit_AV[v] = m_improvedmichelTool->GetMichelFitAV();

      matched_michel_end_X1[v] = m_improvedmichelTool->GetMichelEndpointX1();
      matched_michel_end_Y1[v] = m_improvedmichelTool->GetMichelEndpointY1();
      matched_michel_end_U1[v] = m_improvedmichelTool->GetMichelEndpointU1();
      matched_michel_end_V1[v] = m_improvedmichelTool->GetMichelEndpointV1();
      matched_michel_end_Z1[v] = m_improvedmichelTool->GetMichelEndpointZ1();

      // new 2018.11.21      
      // get the first module corresponding to vertex v's michel's endpoint position     
      if( matched_michel_end_Z1[v] > 0 )
      {
        Gaudi::XYZPoint hitXYZMEnd1( matched_michel_end_X1[v], matched_michel_end_Y1[v], matched_michel_end_Z1[v] );
        const Minerva::DePlane* trkrHitPlaneEnd1   =  m_InnerDetector->getDePlane(hitXYZMEnd1);
        if( !trkrHitPlaneEnd1 )  // Cases where we're not in a tracker module (and I'm too lazy to find the regular module
        {
          matched_michel_end_module1[v] = -6;
        }
        else
        {
          const Minerva::PlaneID& trkrHitPlaneIdEnd1 = trkrHitPlaneEnd1->getPlaneID();
          matched_michel_end_module1[v] = trkrHitPlaneIdEnd1.module();
        }
      }

      matched_michel_end_X2[v] = m_improvedmichelTool->GetMichelEndpointX2();
      matched_michel_end_Y2[v] = m_improvedmichelTool->GetMichelEndpointY2();
      matched_michel_end_U2[v] = m_improvedmichelTool->GetMichelEndpointU2();
      matched_michel_end_V2[v] = m_improvedmichelTool->GetMichelEndpointV2();
      matched_michel_end_Z2[v] = m_improvedmichelTool->GetMichelEndpointZ2();

      // new 2018.11.21
      // get the second module corresponding to vertex v's michel's endpoint position

      if( matched_michel_end_Z2[v] > 0 )
      {
        Gaudi::XYZPoint hitXYZMEnd2( matched_michel_end_X2[v], matched_michel_end_Y2[v], matched_michel_end_Z2[v] );
        const Minerva::DePlane* trkrHitPlaneEnd2   = m_InnerDetector->getDePlane(hitXYZMEnd2);
        if( !trkrHitPlaneEnd2 )  // Cases where we're not in a tracker module (and I'm too lazy to find the regular module
        {
          matched_michel_end_module2[v] = -6;
        }
        else
        {
          const Minerva::PlaneID& trkrHitPlaneIdEnd2 = trkrHitPlaneEnd2->getPlaneID();
          matched_michel_end_module2[v] = trkrHitPlaneIdEnd2.module();
        }
      }
    }
  }

  // Setup for getting truth info of the michels themselves
    Minerva::TG4Trajectories* trajectories=NULL;
  if( haveNeutrinoMC() ) {
    //get list of TG4Trajectories in event in order to track lineage of Michel electron
    trajectories = get<Minerva::TG4Trajectories>( "MC/TG4Trajectories" );
  }
  //===========================================================================
  // LOOP MICHELS -- info about the michels themselves
  //===========================================================================
  for(int m = 0; m < michel_count; m++) {
    const SmartRefVector<Minerva::IDCluster>& michel_clusters =
                                    m_improvedmichelTool->GetMichelClusters(m);

    int michel_nClusters     = michel_clusters.size();
    has_michel_nClusters[m]  = michel_nClusters;

    // new 2018.11.21
    // time slice of michel cluster m
    Minerva::TimeSlice* time_slice = getTimeSlice( michel_clusters[0] );
    int slice = time_slice ? time_slice->sliceNumber() : -1;
    has_michel_slice[m] = slice;

    has_michel_cal_energy[m] = m_improvedmichelTool->GetMichelCalorimetricEnergy(m);
    has_michel_vis_energy[m] = m_improvedmichelTool->GetMichelVisibleEnergy(m);
    has_michel_time[m]       = m_improvedmichelTool->GetMichelTime(m);
    has_michel_views[m]      = m_improvedmichelTool->GetMichelViewCode(m);
    has_michel_position_x[m] = m_improvedmichelTool->GetMichelPositionX(m);
    has_michel_position_y[m] = m_improvedmichelTool->GetMichelPositionY(m);
    has_michel_position_u[m] = m_improvedmichelTool->GetMichelPositionU(m);
    has_michel_position_v[m] = m_improvedmichelTool->GetMichelPositionV(m);
    has_michel_position_z[m] = m_improvedmichelTool->GetMichelPositionZ(m);

    // new 2018.11.21. module numbers corresponding to positions
    Gaudi::XYZPoint hitXYZ( has_michel_position_x[m], has_michel_position_y[m], has_michel_position_z[m] );
    const Minerva::DePlane* hasTrkrHitPlaneAvg    =  m_InnerDetector->getDePlane(hitXYZ);
    if( !hasTrkrHitPlaneAvg )  // Cases where we're not in a tracker module
    {
      has_michel_avg_module[m] = -6;
    }
    else
    {
      const Minerva::PlaneID& hasTrkrHitPlaneIdAvg  = hasTrkrHitPlaneAvg->getPlaneID();
      has_michel_avg_module[m] = hasTrkrHitPlaneIdAvg.module();
    }

    //Endpoint Z
    has_michel_end_Z1[m]     = m_improvedmichelTool->GetMichelEndpointZ1(m);
    has_michel_end_Z2[m]     = m_improvedmichelTool->GetMichelEndpointZ2(m);
    bool fitPass = m_improvedmichelTool->GetMichelFitPass(m);
    has_michel_fitPass[m] = (int)fitPass;

    if(fitPass){
      has_michel_fit_residual[m] = m_improvedmichelTool->GetMichelFitResidual(m);

      has_michel_fit_AX[m] = m_improvedmichelTool->GetMichelFitAX(m);
      has_michel_fit_AY[m] = m_improvedmichelTool->GetMichelFitAY(m);
      has_michel_fit_AU[m] = m_improvedmichelTool->GetMichelFitAU(m);
      has_michel_fit_AV[m] = m_improvedmichelTool->GetMichelFitAV(m);

      has_michel_end_X1[m] = m_improvedmichelTool->GetMichelEndpointX1(m);
      has_michel_end_Y1[m] = m_improvedmichelTool->GetMichelEndpointY1(m);
      has_michel_end_U1[m] = m_improvedmichelTool->GetMichelEndpointU1(m);
      has_michel_end_V1[m] = m_improvedmichelTool->GetMichelEndpointV1(m);

      has_michel_end_X2[m] = m_improvedmichelTool->GetMichelEndpointX2(m);
      has_michel_end_Y2[m] = m_improvedmichelTool->GetMichelEndpointY2(m);
      has_michel_end_U2[m] = m_improvedmichelTool->GetMichelEndpointU2(m);
      has_michel_end_V2[m] = m_improvedmichelTool->GetMichelEndpointV2(m);

      // new 2018.11.21. module numbers corresponding to positions
      Gaudi::XYZPoint hitXYZEnd1( has_michel_end_X1[m], has_michel_end_Y1[m], has_michel_end_Z1[m] );
      const Minerva::DePlane* hasTrkrHitPlaneEnd1  = m_InnerDetector->getDePlane(hitXYZEnd1);
      if( !hasTrkrHitPlaneEnd1 )  // Cases where we're not in a tracker module
      {
        has_michel_end_module1[m] = -6;
      }
      else
      {
        const Minerva::PlaneID& hasTrkrHitPlaneIdEnd1= hasTrkrHitPlaneEnd1->getPlaneID();
        has_michel_end_module1[m] = hasTrkrHitPlaneIdEnd1.module();
      }

      Gaudi::XYZPoint hitXYZEnd2( has_michel_end_X2[m], has_michel_end_Y2[m], has_michel_end_Z2[m] );
      const Minerva::DePlane* hasTrkrHitPlaneEnd2   = m_InnerDetector->getDePlane(hitXYZEnd2);
      if( !hasTrkrHitPlaneEnd2 )  // Cases where we're not in a tracker module
      {
        has_michel_end_module2[m] = -6;
      }
      else
      {
        const Minerva::PlaneID& hasTrkrHitPlaneIdEnd2 = hasTrkrHitPlaneEnd2->getPlaneID();
        has_michel_end_module2[m] = hasTrkrHitPlaneIdEnd2.module();
      }
    } // end fitPass

    // Michel Truth info
    if( !trajectories ) continue;
    if( haveNeutrinoMC() ) {
      Minerva::TG4Trajectory* trueMichel = NULL;
      const Minerva::TG4Trajectory* MichelParent = NULL;
      const Minerva::TG4Trajectory* MichelPrimaryParent = NULL;
      double fraction;
      double other_energy;
      SmartRefVector<Minerva::IDCluster> michel_clusters =
                                     m_improvedmichelTool->GetMichelClusters(m);
      std::vector<const Minerva::IDCluster*> functional_michel_clusters;
      for(unsigned int c = 0; c < michel_clusters.size(); c++)
        functional_michel_clusters.push_back(michel_clusters[c]);

      // check if overlay
      if(TruthMatcher->getDataFraction(functional_michel_clusters) < 0.5) {
        //find TG4Trajectory that created the Michel leading the the Michel track
        StatusCode parent_sc =
          TruthMatcher->getTG4Trajectory(functional_michel_clusters,
                                          MichelParent,fraction,other_energy);
        if(!parent_sc) continue;

        //use final position of Michel parent to get Michel trajectory
        for(Minerva::TG4Trajectories::iterator it = trajectories->begin();
                                                it != trajectories->end(); ++it) {
          Minerva::TG4Trajectory * traject = *it;
          //if(traject->GetInitialPosition() == MichelParent->GetFinalPosition()) {
          if(traject->GetParentId() == MichelParent->GetTrackId()){
            if(abs(traject->GetPDGCode())==14) continue;
            trueMichel = traject;

            //Check if the Michel came from a decay or some kind of scattering
            if(trueMichel->GetProcessName().find("Decay") != std::string::npos)
              michel_from_decay[m] = true;

            true_michel_pdg[m] = trueMichel->GetPDGCode();
            true_michel_x1[m] = trueMichel->GetInitialPosition().X();
            true_michel_y1[m] = trueMichel->GetInitialPosition().Y();
            true_michel_z1[m] = trueMichel->GetInitialPosition().Z();
            true_michel_x2[m] = trueMichel->GetFinalPosition().X();
            true_michel_y2[m] = trueMichel->GetFinalPosition().Y();
            true_michel_z2[m] = trueMichel->GetFinalPosition().Z();
            break;
          }
        }
        michel_parent_PDG[m]                   = MichelParent->GetPDGCode();
        michel_parent_trackID[m]               = MichelParent->GetTrackId();

        MichelPrimaryParent                    = TruthMatcher->getPrimaryTG4Trajectory(MichelParent,false);
        has_michel_primary_particle_trackID[m] = MichelPrimaryParent->GetTrackId();
        has_michel_primary_particle_PDG[m]     = MichelPrimaryParent->GetPDGCode();
        has_michel_primary_particle_p[m]       = MichelPrimaryParent->GetInitialMomentum().P();
        has_michel_primary_particle_px[m]      = MichelPrimaryParent->GetInitialMomentum().px();
        has_michel_primary_particle_py[m]      = MichelPrimaryParent->GetInitialMomentum().py();
        has_michel_primary_particle_pz[m]      = MichelPrimaryParent->GetInitialMomentum().pz();
        has_michel_primary_particle_E[m]       = MichelPrimaryParent->GetInitialMomentum().E();
      }
    } // end haveNeutrinoMC
  } // end michel loop

  //============================================================================
  // Fill Branches
  //============================================================================
    event->setIntData             ( "has_michel_count",      michel_count);
    event->setContainerDoubleData ( "has_michel_cal_energy", has_michel_cal_energy);
    event->setContainerDoubleData ( "has_michel_vis_energy", has_michel_vis_energy);
    event->setContainerDoubleData ( "has_michel_time",       has_michel_time);
    event->setContainerIntData    ( "has_michel_slice",      has_michel_slice);
    event->setContainerIntData    ( "has_michel_views",      has_michel_views);
    event->setContainerIntData    ( "has_michel_nClusters",  has_michel_nClusters);

    event->setContainerDoubleData("has_michel_position_x",has_michel_position_x);
    event->setContainerDoubleData("has_michel_position_y",has_michel_position_y);
    event->setContainerDoubleData("has_michel_position_u",has_michel_position_u);
    event->setContainerDoubleData("has_michel_position_v",has_michel_position_v);
    event->setContainerDoubleData("has_michel_position_z",has_michel_position_z);

    event->setContainerIntData("has_michel_fitPass",has_michel_fitPass);
    event->setContainerDoubleData("has_michel_fit_residual",has_michel_fit_residual );
    event->setContainerDoubleData("has_michel_fit_AX",has_michel_fit_AX);
    event->setContainerDoubleData("has_michel_fit_AY",has_michel_fit_AY);
    event->setContainerDoubleData("has_michel_fit_AU",has_michel_fit_AU);
    event->setContainerDoubleData("has_michel_fit_AV",has_michel_fit_AV);
    event->setContainerDoubleData("has_michel_end_X1",has_michel_end_X1);
    event->setContainerDoubleData("has_michel_end_Y1",has_michel_end_Y1);
    event->setContainerDoubleData("has_michel_end_U1",has_michel_end_U1);
    event->setContainerDoubleData("has_michel_end_V1",has_michel_end_V1);
    event->setContainerDoubleData("has_michel_end_Z1",has_michel_end_Z1);
    event->setContainerDoubleData("has_michel_end_X2",has_michel_end_X2);
    event->setContainerDoubleData("has_michel_end_Y2",has_michel_end_Y2);
    event->setContainerDoubleData("has_michel_end_U2",has_michel_end_U2);
    event->setContainerDoubleData("has_michel_end_V2",has_michel_end_V2);
    event->setContainerDoubleData("has_michel_end_Z2",has_michel_end_Z2);

    // new 2018.11.21. module numbers corresponding to positions
    event->setContainerIntData( "has_michel_avg_module" ,has_michel_avg_module );
    event->setContainerIntData( "has_michel_end_module1",has_michel_end_module1);
    event->setContainerIntData( "has_michel_end_module2",has_michel_end_module2);

    event->setContainerIntData("matched_michel_code"   ,matched_michel_code   );
    event->setContainerIntData("matched_michel_avg_idx",matched_michel_avg_idx);
    event->setContainerIntData("matched_michel_end_idx",matched_michel_end_idx);
    event->setContainerIntData("matched_michel_ov_idx" ,matched_michel_ov_idx );

    event->setContainerDoubleData("matched_michel_fit_residual", matched_michel_fit_residual );
    event->setContainerIntData("matched_michel_end_module1",  matched_michel_end_module1);
    event->setContainerIntData("matched_michel_end_module2",  matched_michel_end_module2);

    if( haveNeutrinoMC() ){
      event->setContainerIntData("has_michel_primary_particle_trackID",has_michel_primary_particle_trackID);
      event->setContainerIntData("michel_from_decay",michel_from_decay);
      event->setContainerIntData("true_michel_pdg",true_michel_pdg);
      event->setContainerDoubleData("true_michel_x1",true_michel_x1);
      event->setContainerDoubleData("true_michel_y1",true_michel_y1);
      event->setContainerDoubleData("true_michel_z1",true_michel_z1);
      event->setContainerDoubleData("true_michel_x2",true_michel_x2);
      event->setContainerDoubleData("true_michel_y2",true_michel_y2);
      event->setContainerDoubleData("true_michel_z2",true_michel_z2);
      event->setContainerIntData("michel_parent_PDG",michel_parent_PDG);
      event->setContainerIntData("michel_parent_trackID",michel_parent_trackID);
      event->setContainerIntData("has_michel_primary_particle_PDG",has_michel_primary_particle_PDG);
      event->setContainerDoubleData("has_michel_primary_particle_p",has_michel_primary_particle_p);
      event->setContainerDoubleData("has_michel_primary_particle_px",has_michel_primary_particle_px);
      event->setContainerDoubleData("has_michel_primary_particle_py",has_michel_primary_particle_py);
      event->setContainerDoubleData("has_michel_primary_particle_pz",has_michel_primary_particle_pz);
      event->setContainerDoubleData("has_michel_primary_particle_E",has_michel_primary_particle_E);
    }
}

//=======================================================================================================
// Find proton prongs, pick the most energetic one and save the others in a vector of prongs
// -protonProng and protonPart are the Minerva::Prong and Minerva::Particle of the
//  most energetic proton found
// secondaryProtonProngs is a vector of all extra proton prongs found besides the most energetic one.
//=======================================================================================================
bool MasterAnaDev::getProtonProngs( Minerva::ProngVect prongs, Minerva::ProngVect& secondaryProtonProngs, SmartRef<Minerva::Prong> &protonProng, SmartRef<Minerva::Particle> &protonPart ) const
{
  debug() << "Enter MasterAnaDev::getProtonProngs" << endmsg;

  bool foundProton = false;
  double max_proton_E = -1.0;
  std::map< SmartRef<Minerva::Prong>, SmartRef<Minerva::Particle> > protonProngPartMap;

  //-- loop over prongs
  Minerva::ProngVect::iterator prong;
  for( prong = prongs.begin(); prong != prongs.end(); prong++ ){

    verbose() << "   prong bit-field signature = " << (*prong)->typeBitsToString() << endmsg;
    int numProtonPartsinProng  = 0;
    //bool prongHasProtonPartandPassProtonCuts = false;
    SmartRef<Minerva::Particle> prongPart;

    //-- continue only with prongs that have tracks in the inner detector
    Minerva::TrackVect tracks = (*prong)->minervaTracks();
    if( tracks.empty() ) {
      warning() << "  This prong contains an empty vector of tracks, skipping!" << endmsg;
      continue;
    } else if( (*prong)->MinosTrack() || (*prong)->MinosStub() ) {
      verbose() << "  This is a MINOS matched prong, skipping!" << endmsg;
      continue;
    }

    //-- make use odMatch prongs an option
    if( !m_useOdMatchProtons && (*prong)->OdMatch() ) {
      verbose() << "  The prong is OdMatch, skipping!" << endmsg;
      continue;
    }

    //-- get the endpoint of the prong's last track
    SmartRef<Minerva::Track> track = tracks[tracks.size() -1];
    Gaudi::XYZPoint endpoint = track->lastState().position();

    //-- check the prong's destruction vertex is within the specify volume
    if( m_minCoordSysTool->inFiducial(endpoint.x(),endpoint.y(),endpoint.z(),m_ProtonApothem,m_ProtonZLow,m_ProtonZHigh) ) {

      //-- get the prong's particles
      Minerva::ParticleVect partHypVec = (*prong)->particles();
      Minerva::ParticleVect::iterator part;

      debug() << "  Considering a prong with " << partHypVec.size() << " particle hypotheses and " << tracks.size() << " tracks." << endmsg;

      //-- loop over all particle hypotheses
      for( part = partHypVec.begin(); part != partHypVec.end(); part++ ){

        debug() << "   Found a " << (*part)->idcode() << " with signature: "
                << (*part)->methodSignature() << " and score: " << (*part)->score() << endmsg;

        //---------------------------------------------------------------------------------------
        //! for now we will only accept protons whose mass was not change during the fit.
        //! there has not been sufficient studies for us to choose the proton-pion fit over
        //! the proton-proton fit based on the score alone.
        //!--------------------------------------------------------------------------------------

        //-- proton conditional statement
        if( (*part)->idcode() == Minerva::Particle::Proton && !(*part)->isMultiMass() ) {
          bool passProtonCut = false;

          debug() << "   Found a proton particle hypothesis!" << endmsg;

          //-- proton's requirements
          if( (*part)->methodSignature().find("dEdX") != std::string::npos ) {
            if ( track->chi2PerDoF() < m_maxProtonChi2 && (*part)->score() > m_minProtonScore ) {
              //prongHasProtonPartandPassProtonCuts = true;
              passProtonCut = true;
            } else {
              debug() << "  Chi2/DoF and/or score is not consistent with a proton" << endmsg;
            }
          } else {
            debug() << "  Method signature is not consistent with a proton" << endmsg;
          }

          //-- proton's pass cut conditional
          if( passProtonCut ) {
            prongPart = (*part);
            debug() << "   Method signature and/or score is consistent with a proton!" << endmsg;
            numProtonPartsinProng++;

          } //-- if pass cut
        } //-- if proton particle
      } //-- for (parts)
    } //-- if pass z requirements
    debug() << "numProtonPartsinProng = " << numProtonPartsinProng << endmsg;

    debug() << "Adding prongs with one single proton part" <<endmsg;
    //-- Only accept prongs with one proton part
    if ( numProtonPartsinProng == 1 ){
      (*prong)->updateBestParticle( prongPart );
      protonProngPartMap[ (*prong) ] = prongPart;
      (*prong)->filtertaglist()->setOrAddFilterTag(m_secondaryProtons, true);
      if (prongPart->momentumVec().E() > max_proton_E)
        max_proton_E = prongPart->momentumVec().E();
    }
  } //-- for(prongs)

  //-----------------------------------------------------------------
  // -- Look for the most energetic proton prong
  //-----------------------------------------------------------------
  debug()<< " Look for most energetic proton prong in event " << endmsg;
  if (max_proton_E < 0 )
    {
      debug()<<" Ups, max_proton_E = " << max_proton_E << endmsg;
      return false;
    }
  for( std::map< SmartRef<Minerva::Prong>, SmartRef<Minerva::Particle> >::iterator it = protonProngPartMap.begin(); it != protonProngPartMap.end(); ++it ) {
    SmartRef<Minerva::Prong> prong = it->first;
    SmartRef<Minerva::Particle> part = it->second;
    if ( part->momentumVec().E() == max_proton_E )
      {
        protonProng = prong;
        protonPart  = part;
        foundProton = true;
        //-- color the proton prong
        m_hitTaggerTool->applyColorTag( protonProng, m_protonProngColor );
        //tag most energetic proton as the primary Proton
        protonProng->filtertaglist()->setFilterTag(m_secondaryProtons, false);
        protonProng->filtertaglist()->setOrAddFilterTag(m_primaryProton, true);
      }
    else
      secondaryProtonProngs.push_back( prong );
  }
  return foundProton;

}


//========================================
// ImprovedtagMichels
//========================================
bool MasterAnaDev::ImprovedtagMichels(Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth) const
{
  bool foundImprovedMichel = false;
  std::cout << "attempting new michel" << std::endl;
  //--tag truth michel
  if (truth){
    bool hasTruthMichelFiducial = hasTruthMichelElectron(event);
    event->setIntData("truth_improved_michel_electron",hasTruthMichelFiducial);
  }

  //--get the tagged muon prong
  Minerva::ProngVect primaryProngs = event->primaryProngs();
  SmartRef<Minerva::Prong> muonProng;

  for (ProngVect::iterator itProngs = primaryProngs.begin(); itProngs !=primaryProngs.end(); ++itProngs) {
    if( (*itProngs)->filtertaglist()->filterTagExists( AnaFilterTags::PrimaryMuon() ) ) {
      muonProng = *itProngs;
      break;
    }
  }

  //--look for michel electrons in vertices
  Minerva::Prong vtx_improvedmichelProng;
  const SmartRefVector<Minerva::Vertex> vertices = event->select<Minerva::Vertex>();
  debug()<<"Looking for Improved michels in " <<vertices.size()-1<<" vertex points" <<endmsg;
  //event->setIntData("has_n_vertex_points", vertices.size()-1);
  warning()<<"Loking for Improved michels in " <<vertices.size()-1<<" vertex points" <<endmsg;

  std::vector<int> improved_michel_vertex_type;
  std::vector<int> improved_michel_in_vertex_point;

  std::vector<double> improved_michel_tvec;
  std::vector<double> improved_michel_tdiff_vec;
  std::vector<double> improved_michel_xvec;
  std::vector<double> improved_michel_yvec;
  std::vector<double> improved_michel_uvec;
  std::vector<double> improved_michel_vvec;
  std::vector<double> improved_michel_zvec;
  std::vector<double> improved_michel_dist_vec;
  std::vector<double> improved_michel_evis_vec;
  std::vector<double> improved_michel_ecalo_vec;

  std::vector<double> improved_michel_hit_charges;
  std::vector<double> improved_michel_hit_times;
  std::vector<double> improved_michel_hit_time_diff_vtx;
  std::vector<double> improved_michel_hit_time_diff_cluster;
  std::vector<double> improved_michel_matched_energy_fraction;
  std::vector<double> improved_michel_data_energy_fraction;

  std::vector<int> improved_michel_match_vec;
  std::vector<int> improved_michel_view_vec;
  std::vector<int> improved_michel_matched_pdg;
  std::vector<int> improved_michel_matched_primary_pdg;
  std::vector<int> improved_michel_ndigits;
  double muonTime = m_recoObjTimeTool->trackVertexTime(muonProng->minervaTracks()[0]);
  const double max_distance = 300.0; //in mm

  //-- look for the Improved Michel
  //while( improved_michel ) {
  m_improvedmichelTool->Reset();
  m_improvedmichelTool->FindMichel(event->earliestSliceNumber());
  m_improvedmichelTool->ApplyQualityCuts();


  //--loop over the vertices
  for ( SmartRefVector<Minerva::Vertex>::const_iterator itVtx = vertices.begin(); itVtx != vertices.end(); ++itVtx ) {

    SmartRef<Minerva::Vertex> myvtx = (*itVtx);
    //--no need to look at the endpoint of the muon track
    if ( (*itVtx)->getIncomingTrack()==muonProng->minervaTracks()[0] ) continue;
    if ( !m_masterAnaDevRecoUtils->vertexTracksInTime(myvtx,muonTime)) continue;//if any tracks associated with the vertex are outside tracking time window don't use it

    verbose() << "This vertex is not at the end point of a Muon. Vertex type is = " << (*itVtx)->type() << endmsg;

    //bool improved_michel = true;

    Gaudi::XYZPoint pos = (*itVtx)->position();
    double x = pos.x(); double y = pos.y(); double z = pos.z();
    double t = (*itVtx)->getTime();

    //--grab vtx info
    int improved_vtx_type = -1;


    if( (*itVtx)->type() == Minerva::Vertex::StartPoint )
      improved_vtx_type = 0;
    else if ( (*itVtx)->type() == Minerva::Vertex::Kinked )
      improved_vtx_type = 1;
    else if ( (*itVtx)->type() == Minerva::Vertex::StopPoint )
      improved_vtx_type = 2;
    else
      debug()<<"Could not determine vertex point!"<<endmsg;


    bool match = m_improvedmichelTool->Match( x, y, z, max_distance, max_distance);
    if (!match){
      debug()<<"'Did NOT Found a michel electron' Said the ImprovedMichelTool" <<endmsg;
      improved_vtx_type = 99;
      continue;
    }
    debug()<<"'Found a michel electron' Said the ImprovedMichelTool" <<endmsg;

    //--reset michel variables
    double michel_t = -9999.9999;
    double michel_tdiff = -9999.9999;
    double michel_x = -9999.9999;
    double michel_y = -9999.9999;
    double michel_u = -9999.9999;
    double michel_v = -9999.9999;
    double michel_z = -9999.9999;
    double michel_dist = -9999.9999;
    double michel_ecalo = -9999.9999;
    double michel_evis  = -9999.9999;
    int michel_view = -1;
    std::cout << "FOUND ONE" << std::endl;
    std::vector<const Minerva::IDCluster*> idclusters;
    const SmartRefVector<Minerva::IDCluster> clustersfromcandidate = m_improvedmichelTool->GetMichelClusters();
    idclusters.insert(idclusters.begin(), clustersfromcandidate.begin(), clustersfromcandidate.end());
    SmartRefVector<Minerva::IDDigit> alliddigits;
    //Loop over id clusters, and extract hit information
    for( SmartRefVector<Minerva::IDCluster>::const_iterator itClus = clustersfromcandidate.begin(); itClus != clustersfromcandidate.end(); ++itClus ){
      double clus_t = (*itClus)->time();
      SmartRefVector<Minerva::IDDigit> iddigits = (*itClus)->digits();
      alliddigits.insert(alliddigits.end(),iddigits.begin(),iddigits.end());
      for( SmartRefVector<Minerva::IDDigit>::const_iterator itDigit = iddigits.begin(); itDigit != iddigits.end(); ++itDigit ){
        improved_michel_hit_charges.push_back((*itDigit)->normEnergy());
        improved_michel_hit_times.push_back((*itDigit)->calTime());
        improved_michel_hit_time_diff_vtx.push_back((*itDigit)->calTime()-t);
        improved_michel_hit_time_diff_cluster.push_back((*itDigit)->calTime()-clus_t);
      }//end loop digits
    }//end loop clusters
    //Do truth matching if this is MC
    if (truth){
      std::vector<double> bestPart = tagDigitsTruth(alliddigits);
      double data_fraction = TruthMatcher->getDataFraction(idclusters);
      improved_michel_matched_pdg.push_back(bestPart[0]);
      improved_michel_matched_primary_pdg.push_back(bestPart[1]);
      improved_michel_matched_energy_fraction.push_back(bestPart[2]);
      improved_michel_data_energy_fraction.push_back(data_fraction);
    }
    else{
      improved_michel_matched_pdg.push_back(-1);
      improved_michel_matched_primary_pdg.push_back(-1);
      improved_michel_matched_energy_fraction.push_back(-1.0);
      improved_michel_data_energy_fraction.push_back(-1.0);
    }

    michel_t = m_improvedmichelTool->GetMichelTime();
    michel_tdiff = (michel_t - t);

    michel_x = m_improvedmichelTool->GetMichelPositionX();
    michel_y = m_improvedmichelTool->GetMichelPositionY();
    michel_u = m_improvedmichelTool->GetMichelPositionU();
    michel_v = m_improvedmichelTool->GetMichelPositionV();
    michel_z = m_improvedmichelTool->GetMichelPositionZ();
    michel_dist = sqrt(pow(michel_x - x,2) + pow(michel_y - y,2) + pow(michel_z -z,2));

    michel_ecalo = m_improvedmichelTool->GetMichelCalorimetricEnergy();
    michel_evis  = m_improvedmichelTool->GetMichelVisibleEnergy();
    michel_view  = m_improvedmichelTool->GetMichelViewCode();


    //--push variable sin the vectors (ONLY if matched)
    improved_michel_vertex_type.push_back(improved_vtx_type);
    improved_michel_in_vertex_point.push_back( itVtx - vertices.begin() );

    improved_michel_tvec.push_back(michel_t);
    improved_michel_tdiff_vec.push_back(michel_tdiff);
    improved_michel_xvec.push_back(michel_x);
    improved_michel_yvec.push_back(michel_y);
    improved_michel_uvec.push_back(michel_u);
    improved_michel_vvec.push_back(michel_v);
    improved_michel_zvec.push_back(michel_z);
    improved_michel_dist_vec.push_back(michel_dist);
    improved_michel_evis_vec.push_back(michel_evis);
    improved_michel_ecalo_vec.push_back(michel_ecalo);
    improved_michel_view_vec.push_back(michel_view);
    improved_michel_match_vec.push_back(int(match));
    improved_michel_ndigits.push_back(alliddigits.size());

  } //End of loop over vertices

  event->setContainerIntData("improved_michel_vertex_type", improved_michel_vertex_type);
  event->setContainerIntData("improved_michel_in_vertex_point", improved_michel_in_vertex_point);
  event->setIntData("improved_nmichel", improved_michel_vertex_type.size());

  event->setContainerDoubleData("improved_michel_tvec", improved_michel_tvec);
  event->setContainerDoubleData("improved_michel_tdiff_vec", improved_michel_tdiff_vec);
  event->setContainerDoubleData("improved_michel_xvec", improved_michel_xvec);
  event->setContainerDoubleData("improved_michel_yvec", improved_michel_yvec);
  event->setContainerDoubleData("improved_michel_uvec", improved_michel_uvec);
  event->setContainerDoubleData("improved_michel_vvec", improved_michel_vvec);
  event->setContainerDoubleData("improved_michel_zvec", improved_michel_zvec);
  event->setContainerDoubleData("improved_michel_dist_vec", improved_michel_dist_vec);
  event->setContainerDoubleData("improved_michel_evis_vec", improved_michel_evis_vec);
  event->setContainerDoubleData("improved_michel_ecalo_vec", improved_michel_ecalo_vec);
  event->setContainerDoubleData("improved_michel_hit_charges", improved_michel_hit_charges);
  event->setContainerDoubleData("improved_michel_hit_times", improved_michel_hit_times);
  event->setContainerDoubleData("improved_michel_hit_time_diff_vtx", improved_michel_hit_time_diff_vtx);
  event->setContainerDoubleData("improved_michel_hit_time_diff_cluster", improved_michel_hit_time_diff_cluster);
  event->setContainerDoubleData("improved_michel_matched_energy_fraction",improved_michel_matched_energy_fraction);
  event->setContainerDoubleData("improved_michel_data_energy_fraction",improved_michel_data_energy_fraction);
  event->setContainerIntData("improved_michel_view_vec", improved_michel_view_vec);
  event->setContainerIntData("improved_michel_match_vec", improved_michel_match_vec); // whether this track is Michel tagged
  event->setContainerIntData("improved_michel_matched_pdg",improved_michel_matched_pdg);
  event->setContainerIntData("improved_michel_matched_primary_pdg",improved_michel_matched_primary_pdg);
  event->setContainerIntData("improved_michel_ndigits",improved_michel_ndigits);
  verbose() << "Finished setting Improved michel variables" << endmsg;

  if (improved_michel_view_vec.size()>0) foundImprovedMichel = true;

  return foundImprovedMichel;
}

//=====================================
// hasTruthMichelElectron
//=====================================
bool MasterAnaDev::hasTruthMichelElectron(Minerva::PhysicsEvent* event) const {

  debug()<<"Entering MasterAnaDev::hasTruthMichelElectron(  ) ....."<< endmsg;

  bool hasMichel = false;

  Gaudi::LorentzVector positron_pos;
  Gaudi::LorentzVector electron_pos;

  std::map<int,double> pion_plus_ids_momentum, pion_minus_ids_momentum;
  std::map<int,int> muon_plus_id_parent, muon_minus_id_parent;
  std::map<int,Gaudi::LorentzVector> electron_parent_pos, positron_parent_pos;
  std::vector<double> pion_plus_momentum, pion_minus_momentum;

  Minerva::TG4Trajectories* alltrajects = get<Minerva::TG4Trajectories>( "MC/TG4Trajectories" );
  for(Minerva::TG4Trajectories::iterator it = alltrajects->begin(); it != alltrajects->end(); ++it) {
    Minerva::TG4Trajectory* traject = *it;

    //! look for charged pions
    if( traject->GetPDGCode()==211 ){
      pion_plus_ids_momentum[traject->GetTrackId()] = (traject->GetInitialMomentum()).P();
    }
    if( traject->GetPDGCode()==-211 ){
      pion_minus_ids_momentum[traject->GetTrackId()] = (traject->GetInitialMomentum()).P();
    }

    //! look for muons from decays
    if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == -13 ) {
      muon_plus_id_parent[traject->GetTrackId()] = traject->GetParentId();
    }

    if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == 13 ) {
      muon_minus_id_parent[traject->GetTrackId()] = traject->GetParentId();
    }

    //! look for electrons from decays
    if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == -11 ) {
      positron_parent_pos[traject->GetParentId()] = traject->GetInitialPosition();
    }

    if( traject->GetProcessName() == "muMinusCaptureAtRest" && traject->GetPDGCode() == 11 ) {
      electron_parent_pos[traject->GetParentId()] = traject->GetInitialPosition();
    }

  } // End of TG4Trajectories loop

  //! loop over electrons
  for( std::map<int,Gaudi::LorentzVector>::iterator electron = electron_parent_pos.begin(); electron != electron_parent_pos.end(); ++electron ){
    // check if electron parent is a muon-
    if( muon_minus_id_parent.count( electron->first ) ){
      //check if electron is inside fiducial volume
      Gaudi::LorentzVector electron_pos = electron->second;
      Gaudi::XYZPoint p( electron_pos.x(), electron_pos.y(), electron_pos.z() );
      //Open them as they are open in NukeCCInclusive
      //if( FiducialPointTool->isFiducial( p, m_analyzeApothem, m_michel_upstreamZ, m_michel_downstreamZ ) ){
      event->setIntData("truth_has_michel_electron",1);
      hasMichel = true;
      // check if muon parent is a pion-
      if( pion_minus_ids_momentum.count( muon_minus_id_parent[electron->first] ) ){
        pion_minus_momentum.push_back( pion_minus_ids_momentum[ muon_minus_id_parent[electron->first] ]  );
        //}
      }
    }
  }

  //! loop over positrons
  for( std::map<int,Gaudi::LorentzVector>::iterator positron = positron_parent_pos.begin(); positron != positron_parent_pos.end(); ++positron ){
    // check if positron parent is a muon+
    if( muon_plus_id_parent.count( positron->first ) ){
      // check if positron is inside fiducial volume
      Gaudi::LorentzVector positron_pos = positron->second;
      Gaudi::XYZPoint p( positron_pos.x(), positron_pos.y(), positron_pos.z() );
      //open them as they are open in NukeCCInclusive
      //if( FiducialPointTool->isFiducial( p, m_analyzeApothem, m_michel_upstreamZ, m_michel_downstreamZ ) ){
      event->setIntData("truth_has_michel_electron",1);
      hasMichel = true;
      // check if muon parent is a pi+
      if( pion_plus_ids_momentum.count( muon_plus_id_parent[positron->first] ) ){
        pion_plus_momentum.push_back( pion_plus_ids_momentum[ muon_plus_id_parent[positron->first] ]  );
        //}
      }
    }
  }

  event->setContainerDoubleData("truth_has_michel_from_pion_plus_momentum", pion_plus_momentum);
  event->setContainerDoubleData("truth_has_michel_from_pion_minus_momentum", pion_minus_momentum);

  debug()<<"Leaving MasterAnaDev::hasTruthMichelElectron(  ) ....."<<endmsg;
  return hasMichel;
}

bool MasterAnaDev::findLongestMinervaTrack( Minerva::Prong* deProng, SmartRef<Minerva::Track> &longestTrackInsideProng ) const {

  verbose() << "Entering MasterAnaDev::findLongestMinervaTrack( ) ....." << endmsg;

  Minerva::TrackVect minervaTrkVec = deProng->minervaTracks();
  if( minervaTrkVec.size() ) {
    //Find the longest track inside the Prong (can be muon, proton, etc.)
    unsigned int most_nodes = 0;
    for( Minerva::TrackVect::iterator itTrk = minervaTrkVec.begin(); itTrk != minervaTrkVec.end(); ++itTrk ) {
      if( (*itTrk)->nNodes() > most_nodes ) {
        most_nodes = (*itTrk)->nNodes();
        longestTrackInsideProng = (*itTrk);
      }
    } //end of for loop over tracks

    if( longestTrackInsideProng->nNodes() )
      debug() <<"For each prong "<< longestTrackInsideProng->chi2PerDoF() <<" Nodes: "<< longestTrackInsideProng->nNodes() << endmsg;

  } else {
    debug() << "This Prong contains no MINERvA Tracks, so skipping !" << endmsg;
    return false;
  }

  verbose() << "Exiting MasterAnaDev::findLongestMinervaTrack( ) ....." << endmsg;
  return true;
}

bool MasterAnaDev::findLongestMinervaTrack( Minerva::ProngVect deProngVec, SmartRef<Minerva::Track> &longestMinervaTrack ) const {

  bool foundTrackInsideAProng = false;

  //------------------------------------------------------------------------------------------------------
  // This returns the longest track (whichever prong it belongs to) after looping over all the Prongs
  //------------------------------------------------------------------------------------------------------
  unsigned int prong_nodes = 0;
  SmartRef<Minerva::Track> insideProngTrack;

  for( ProngVect::iterator itProngs = deProngVec.begin(); itProngs != deProngVec.end(); ++itProngs ) {
    foundTrackInsideAProng = findLongestMinervaTrack( *itProngs, insideProngTrack );

    if( foundTrackInsideAProng && insideProngTrack->nNodes() ) {
      if( insideProngTrack->nNodes() > prong_nodes ) {
        prong_nodes = (insideProngTrack)->nNodes();
        longestMinervaTrack = insideProngTrack;
      }
    }

  } //End of for loop over Prongs

  if( !longestMinervaTrack->nNodes() ) {
    return false;
  } else {
    debug() << "After looping over all Prongs " << longestMinervaTrack->chi2PerDoF() <<" Nodes: "<< longestMinervaTrack->nNodes() << endmsg;
    return true;
  }
}

bool MasterAnaDev::findLongestMinervaTrack( const Minerva::PhysicsEvent *physEvent, SmartRef<Minerva::Track> &longestMinervaTrack ) const {

  //----------------------------------------------------------------------------------------
  // This returns the longest track found (whichever Primary Prong it belongs to)
  //----------------------------------------------------------------------------------------
  ProngVect primaryProngs = physEvent->primaryProngs();
  if( !primaryProngs.size() ) {
    verbose() << "There are no Primary Prongs in this physics event!" << endmsg;
    return false;
  } else {
    return findLongestMinervaTrack( primaryProngs, longestMinervaTrack );
  }

  return false;
}
//================================
// Get recoil energy from clusters
//================================
string MasterAnaDev::getRecoilChannel(
                                      const int muonCharge,
                                      const int targetID,
                                      const int targetZ
                                      ) const
{

  string channel;

  string prefix = "NukeCC_";
  if(       1 == muonCharge )
    prefix += "AntiNu_";
  else if( -1 == muonCharge )
    prefix += "Nu_";

  if(      1 == targetID ) channel = prefix + "Tgt1";
  else if( 2 == targetID ) channel = prefix + "Tgt2";
  else if( 3 == targetID ) channel = prefix + "Tgt3";
  else if( 4 == targetID ) channel = prefix + "Tgt4";
  else if( 5 == targetID ) channel = prefix + "Tgt5";
  else
    {
      if( m_useCCTrackerTuning )
        channel = ( 1 == muonCharge ) ? "AntiCCInclusive" : "CCInclusive";
      else
        channel = prefix + "Tracker";
    }


  verbose() << " I know targetZ is unused..." << targetZ << endmsg;

  //don't bother with nucleus specific for now
  if( 1 <= targetID && targetID <= 5 )
    {
      if( 6  == targetZ ) channel += "_C";
      if( 26 == targetZ ) channel += "_Fe";
      if( 82 == targetZ ) channel += "_Pb";
    }
  debug() << " channel = " << channel << endmsg;
  return channel;
}

SmartRefVector<Minerva::IDCluster> MasterAnaDev::getAnalyzableNonMuIDClusters( const Minerva::PhysicsEvent* event, const SmartRef<Minerva::Prong>& muon, bool useWideTimeWindow ) const
{
  //get in-time clusters with analyzable type
  bool includeUsed = true;
  SmartRefVector<Minerva::IDCluster> clusters =
    useWideTimeWindow ?
    m_ccUtilsWideTimeWindow->getAnalyzableIDClusters( event, includeUsed ) :
    m_ccUtils->getAnalyzableIDClusters( event, includeUsed );

  //remove muons clusters
  const SmartRefVector<Minerva::IDCluster> muClusters = muon->getAllIDClusters();
  for( SmartRefVector<Minerva::IDCluster>::iterator i = clusters.begin(); i != clusters.end(); )
    {
      if( find( muClusters.begin(), muClusters.end(), *i ) != muClusters.end() )
        i = clusters.erase(i);
      else
        ++i;
    }

  return clusters;
}

SmartRefVector<Minerva::ODCluster> MasterAnaDev::getAnalyzableNonMuODClusters( const Minerva::PhysicsEvent* event, const SmartRef<Minerva::Prong>& muon, bool useWideTimeWindow ) const
{
  //get in-time clusters with analyzable type
  bool includeUsed = true;
  SmartRefVector<Minerva::ODCluster> clusters =
    useWideTimeWindow ?
    m_ccUtilsWideTimeWindow->getAnalyzableODClusters( event, includeUsed ) :
    m_ccUtils->getAnalyzableODClusters( event, includeUsed );

  //remove muons clusters
  const SmartRefVector<Minerva::ODCluster> muClusters = muon->getAllODClusters();
  for( SmartRefVector<Minerva::ODCluster>::iterator i = clusters.begin(); i != clusters.end(); )
    {
      if( find( muClusters.begin(), muClusters.end(), *i ) != muClusters.end() )
        i = clusters.erase(i);
      else
        ++i;
    }

  return clusters;
}

double MasterAnaDev::getRecoilEnergy( const Minerva::PhysicsEvent* event, const SmartRef<Minerva::Prong>& muon, const string& channel, Minerva::NeutrinoInt *ccqeHyp ) const
{
  bool useWideWindow = true;
  SmartRefVector<Minerva::IDCluster> idClusters = getAnalyzableNonMuIDClusters( event, muon, useWideWindow );
  SmartRefVector<Minerva::ODCluster> odClusters = getAnalyzableNonMuODClusters( event, muon, useWideWindow );
  double caloEwide = m_caloUtils->applyCalConsts( idClusters, odClusters, channel );
  ccqeHyp->setDoubleData( "recoil_E_wide_window", caloEwide );

  double caloE = RecoilUtils->calcRecoilEFromClusters( event, muon, channel );
  return caloE;
}
double MasterAnaDev::getvisibleenergy(const Minerva::PhysicsEvent* event, const SmartRef<Minerva::Prong>& muon) const
{
  bool useWideWindow = true;
  SmartRefVector<Minerva::IDCluster> idClusters = getAnalyzableNonMuIDClusters( event, muon, useWideWindow );
  SmartRefVector<Minerva::ODCluster> odClusters = getAnalyzableNonMuODClusters( event, muon, useWideWindow );

  double totalE = 0.0;
  double tempE = 0.0;
  for( SmartRefVector<Minerva::IDCluster>::iterator idClus = idClusters.begin(); idClus != idClusters.end(); ++idClus ){
    SmartRefVector<Minerva::IDDigit> digits = (*idClus)->digits();

    for( SmartRefVector<Minerva::IDDigit>::iterator itDig = digits.begin(); itDig != digits.end(); ++itDig ) {
      tempE = (*itDig)->normEnergy() * (*idClus)->getDigitOwnershipFraction(*itDig);
      totalE += tempE;

    }
  }

  for( SmartRefVector<Minerva::ODCluster>::iterator odClus = odClusters.begin(); odClus != odClusters.end(); ++odClus ){
    SmartRefVector<Minerva::ODDigit> digits = (*odClus)->digits();

    for( SmartRefVector<Minerva::ODDigit>::iterator itDig = digits.begin(); itDig != digits.end(); ++itDig ) {
      tempE = (*itDig)->normEnergy() * (*odClus)->getDigitOwnershipFraction(*itDig);
      totalE += tempE;

    }

  }

  return totalE;
}

//==============================================================================
// Function used to get input for new particle response systematics
//==============================================================================
Minerva::IDClusterVect MasterAnaDev::getHadronIDClusters(
    const Minerva::PhysicsEvent* event,
    SmartRef<Minerva::Prong> muonProng,
    SmartRef<Minerva::Prong> excludeProng ) const {
  verbose()<<"Entered MasterAnaDev::getHadronIDClusters"<<endmsg;

  Minerva::IDClusterVect nonmuon_clusters_id;

  Minerva::IDClusterVect allClustersID = event->select<Minerva::IDCluster>("Used:Unused", "!XTalkCandidate");
  info()<<"Got "<<allClustersID.size()<<" Total ID Clusters"<<endmsg;

  const Minerva::IDClusterVect& muonClustersID    = muonProng->getAllIDClusters();

  const double mintimeCut = event->time() - 20.0*CLHEP::ns;
  const double maxtimeCut = event->time() + 35.0*CLHEP::ns;

  Minerva::IDClusterVect::iterator i_allClustersID = allClustersID.begin();
  for( ; i_allClustersID != allClustersID.end(); ++i_allClustersID ){
    //skip if not in time window
    if( (*i_allClustersID)->time() < mintimeCut ||  maxtimeCut < (*i_allClustersID)->time() ) continue;

    //skip if it's a muon cluster
    if( std::find( muonClustersID.begin(), muonClustersID.end(), (*i_allClustersID) ) != muonClustersID.end() ) continue;

    nonmuon_clusters_id.push_back( *i_allClustersID );
  }

  info()<<"Got "<<nonmuon_clusters_id.size()<<" ID Hadron Clusters No Muon"<<endmsg;

  //This is inefficient, but this is the only way I know how to offload this
  if(excludeProng) {
    Minerva::IDClusterVect tmp_nonmuon_clusters_id;
    const Minerva::IDClusterVect excludeClustersID = excludeProng->getAllIDClusters();
    debug()<<"Got "<<excludeClustersID.size()<<" ID Exclude Clusters"<<endmsg;

    //Minerva::IDClusterVect::iterator i_nonmuon_clusters_id = nonmuon_clusters_id.begin();
    for( uint i = 0; i < nonmuon_clusters_id.size(); ++i ){
      //skip if in hadron cluster
      if( std::find( excludeClustersID.begin(), excludeClustersID.end(), nonmuon_clusters_id[i] ) != excludeClustersID.end() )  continue;

      tmp_nonmuon_clusters_id.push_back( nonmuon_clusters_id[i] );
    }
    nonmuon_clusters_id = tmp_nonmuon_clusters_id;
  }

  debug()<<"Got "<<nonmuon_clusters_id.size()<<" ID Hadron Clusters"<<endmsg;

  verbose()<<"Exit MasterAnaDev::getHadronIDClusters"<<endmsg;
  return nonmuon_clusters_id;
}
//=================================================================================//

int MasterAnaDev::GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane ) const
{
  int targetID = -1;
  vtx_module = -999;
  vtx_plane = -999;

  if( segment == 0 ){
    vtx_module = -5;
    vtx_plane = 0;
    targetID = -1;
  }
  else if( segment == 1 ){
    vtx_module = -5;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 2 ){
    vtx_module = -5;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 3 ){
    vtx_module = -4;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 4 ){
    vtx_module = -4;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 5 ){
    vtx_module = -3;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 6 ){
    vtx_module = -3;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 7 ){
    vtx_module = -2;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 8 ){
    vtx_module = -2;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 9 ){
    vtx_module = -1;
    vtx_plane = 1;
    targetID = 1;
  }
  else if( segment == 10 ){
    vtx_module = 0;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 11 ){
    vtx_module = 0;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 12 ){
    vtx_module = 1;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 13 ){
    vtx_module = 1;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 14 ){
    vtx_module = 2;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 15 ){
    vtx_module = 2;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 16 ){
    vtx_module = 3;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 17 ){
    vtx_module = 3;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 18 ){
    vtx_module = 4;
    vtx_plane = 1;
    targetID = 2;
  }
  else if( segment == 19 ){
    vtx_module = 5;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 20 ){
    vtx_module = 5;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 21 ){
    vtx_module = 6;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 22 ){
    vtx_module = 6;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 23 ){
    vtx_module = 7;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 24 ){
    vtx_module = 7;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 25 ){
    vtx_module = 8;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 26 ){
    vtx_module = 8;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 27 ){
    vtx_module = 9;
    vtx_plane = 2;
    targetID = 3;
  }
  else if( segment == 28 ){
    vtx_module = 11;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 29 ){
    vtx_module = 11;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 30 ){
    vtx_module = 12;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 31 ){
    vtx_module = 12;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 32 ){
    vtx_module = 13;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 33 ){
    vtx_module = 13;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 34 ){
    vtx_module = 14;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 35 ){
    vtx_module = 14;
    vtx_plane = 2;
    targetID = -1;
  }
  //segment 36 is water!!
  else if( segment == 36 ){
    vtx_module = -999;
    vtx_plane = -999;
    targetID = 6;
  }
  else if( segment == 37 ){
    vtx_module = 15;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 38 ){
    vtx_module = 15;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 39 ){
    vtx_module = 16;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 40 ){
    vtx_module = 16;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 41 ){
    vtx_module = 17;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 42 ){
    vtx_module = 17;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 43 ){
    vtx_module = 18;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 44 ){
    vtx_module = 18;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 45 ){
    vtx_module = 19;
    vtx_plane = 1;
    targetID = 4;
  }
  else if( segment == 46 ){
    vtx_module = 20;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 47 ){
    vtx_module = 20;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 48 ){
    vtx_module = 21;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 49 ){
    vtx_module = 21;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 50 ){
    vtx_module = 22;
    vtx_plane = 1;
    targetID = 5;
  }
  else if( segment == 51 ){
    vtx_module = 23;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 52 ){
    vtx_module = 23;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 53 ){
    vtx_module = 24;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 54 ){
    vtx_module = 24;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 55 ){
    vtx_module = 25;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 56 ){
    vtx_module = 25;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 57 ){
    vtx_module = 26;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 58 ){
    vtx_module = 26;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 59 ){
    vtx_module = 27;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 60 ){
    vtx_module = 27;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 61 ){
    vtx_module = 28;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 62 ){
    vtx_module = 28;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 63 ){
    vtx_module = 29;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 64 ){
    vtx_module = 29;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 65 ){
    vtx_module = 30;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 66 ){
    vtx_module = 30;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 67 ){
    vtx_module = 31;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 68 ){
    vtx_module = 31;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 69 ){
    vtx_module = 32;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 70 ){
    vtx_module = 32;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 71 ){
    vtx_module = 33;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 72 ){
    vtx_module = 33;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 73 ){
    vtx_module = 34;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 74 ){
    vtx_module = 34;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 75 ){
    vtx_module = 35;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 76 ){
    vtx_module = 35;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 77 ){
    vtx_module = 36;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 78 ){
    vtx_module = 36;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 79 ){
    vtx_module = 37;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 80 ){
    vtx_module = 37;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 81 ){
    vtx_module = 38;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 82 ){
    vtx_module = 38;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 83 ){
    vtx_module = 39;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 84 ){
    vtx_module = 39;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 85 ){
    vtx_module = 40;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 86 ){
    vtx_module = 40;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 87 ){
    vtx_module = 41;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 88 ){
    vtx_module = 41;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 89 ){
    vtx_module = 42;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 90 ){
    vtx_module = 42;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 91 ){
    vtx_module = 43;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 92 ){
    vtx_module = 43;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 93 ){
    vtx_module = 44;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 94 ){
    vtx_module = 44;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 95 ){
    vtx_module = 45;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 96 ){
    vtx_module = 45;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 97 ){
    vtx_module = 46;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 98 ){
    vtx_module = 46;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 99 ){
    vtx_module = 47;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 100 ){
    vtx_module = 47;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 101 ){
    vtx_module = 48;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 102 ){
    vtx_module = 48;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 103 ){
    vtx_module = 49;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 104 ){
    vtx_module = 49;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 105 ){
    vtx_module = 50;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 106 ){
    vtx_module = 50;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 107 ){
    vtx_module = 51;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 108 ){
    vtx_module = 51;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 109 ){
    vtx_module = 52;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 110 ){
    vtx_module = 52;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 111 ){
    vtx_module = 53;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 112 ){
    vtx_module = 53;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 113 ){
    vtx_module = 54;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 114 ){
    vtx_module = 54;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 115 ){
    vtx_module = 55;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 116 ){
    vtx_module = 55;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 117 ){
    vtx_module = 56;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 118 ){
    vtx_module = 56;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 119 ){
    vtx_module = 57;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 120 ){
    vtx_module = 57;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 121 ){
    vtx_module = 58;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 122 ){
    vtx_module = 58;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 123 ){
    vtx_module = 59;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 124 ){
    vtx_module = 59;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 125 ){
    vtx_module = 60;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 126 ){
    vtx_module = 60;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 127 ){
    vtx_module = 61;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 128 ){
    vtx_module = 61;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 129 ){
    vtx_module = 62;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 130 ){
    vtx_module = 62;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 131 ){
    vtx_module = 63;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 132 ){
    vtx_module = 63;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 133 ){
    vtx_module = 64;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 134 ){
    vtx_module = 64;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 135 ){
    vtx_module = 65;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 136 ){
    vtx_module = 65;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 137 ){
    vtx_module = 66;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 138 ){
    vtx_module = 66;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 139 ){
    vtx_module = 67;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 140 ){
    vtx_module = 67;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 141 ){
    vtx_module = 68;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 142 ){
    vtx_module = 68;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 143 ){
    vtx_module = 69;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 144 ){
    vtx_module = 69;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 145 ){
    vtx_module = 70;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 146 ){
    vtx_module = 70;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 147 ){
    vtx_module = 71;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 148 ){
    vtx_module = 71;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 149 ){
    vtx_module = 72;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 150 ){
    vtx_module = 72;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 151 ){
    vtx_module = 73;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 152 ){
    vtx_module = 73;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 153 ){
    vtx_module = 74;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 154 ){
    vtx_module = 74;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 155 ){
    vtx_module = 75;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 156 ){
    vtx_module = 75;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 157 ){
    vtx_module = 76;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 158 ){
    vtx_module = 76;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 159 ){
    vtx_module = 77;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 160 ){
    vtx_module = 77;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 161 ){
    vtx_module = 78;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 162 ){
    vtx_module = 78;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 163 ){
    vtx_module = 79;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 164 ){
    vtx_module = 79;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 165 ){
    vtx_module = 80;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 166 ){
    vtx_module = 80;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 167 ){
    vtx_module = 81;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 168 ){
    vtx_module = 81;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 169 ){
    vtx_module = 82;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 170 ){
    vtx_module = 82;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 171 ){
    vtx_module = 83;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 172 ){
    vtx_module = 83;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 173 ){
    vtx_module = 84;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 174 ){
    vtx_module = 84;
    vtx_plane = 2;
    targetID = -1;
  }
  //Oscar added this part for DSCAL region
  else if( segment == 175 ){
    vtx_module = 85;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 176 ){
    vtx_module = 85;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 177 ){
    vtx_module = 86;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 178 ){
    vtx_module = 86;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 179 ){
    vtx_module = 87;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 180 ){
    vtx_module = 87;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 181 ){
    vtx_module = 88;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 182 ){
    vtx_module = 88;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 183 ){
    vtx_module = 89;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 184 ){
    vtx_module = 89;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 185 ){
    vtx_module = 90;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 186 ){
    vtx_module = 90;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 187 ){
    vtx_module = 91;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 188 ){
    vtx_module = 91;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 189 ){
    vtx_module = 92;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 190 ){
    vtx_module = 92;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 191 ){
    vtx_module = 93;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 192 ){
    vtx_module = 93;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 193 ){
    vtx_module = 94;
    vtx_plane = 1;
    targetID = -1;
  }
  else if( segment == 194){
    vtx_module = 94;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 195 ){
    vtx_module = 95;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 196 ){
    vtx_module = 96;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 197 ){
    vtx_module = 97;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 198 ){
    vtx_module = 98;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 199 ){
    vtx_module = 99;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 200 ){
    vtx_module = 100;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 201 ){
    vtx_module = 101;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 202 ){
    vtx_module = 102;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 203 ){
    vtx_module = 103;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 204 ){
    vtx_module = 104;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 205 ){
    vtx_module = 105;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 206 ){
    vtx_module = 106;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 207 ){
    vtx_module = 107;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 208 ){
    vtx_module = 108;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 209 ){
    vtx_module = 109;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 210 ){
    vtx_module = 110;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 211 ){
    vtx_module = 111;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 212 ){
    vtx_module = 112;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 213 ){
    vtx_module = 113;
    vtx_plane = 2;
    targetID = -1;
  }
  else if( segment == 214 ){
    vtx_module = 114;
    vtx_plane = 2;
    targetID = -1;
  }
  if( targetID == -999 ) warning() << "Can't find any module/plane for given segment " << segment << endmsg;
  return targetID;
}
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
//Assign the discriminator pair to the wall and paddle
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
std::vector< std::pair<int,int> > MasterAnaDev::discrPairToPaddle( unsigned int minervaID ) const{

  std::vector< std::pair<int,int> > WallPaddle;

  if( minervaID == 1224832 ){
    WallPaddle.push_back(std::make_pair(1,1));
    WallPaddle.push_back(std::make_pair(1,6));
    WallPaddle.push_back(std::make_pair(2,1));
    WallPaddle.push_back(std::make_pair(2,6));
  }
  else if( minervaID == 1224960 ){
    WallPaddle.push_back(std::make_pair(1,2));
    WallPaddle.push_back(std::make_pair(1,3));
    WallPaddle.push_back(std::make_pair(2,2));
    WallPaddle.push_back(std::make_pair(2,3));
  }
  else if( minervaID == 1224992 ){
    WallPaddle.push_back(std::make_pair(1,4));
    WallPaddle.push_back(std::make_pair(1,5));
    WallPaddle.push_back(std::make_pair(2,4));
    WallPaddle.push_back(std::make_pair(2,5));
  }

  return WallPaddle;
}
//=============================================================================
// End of function
//=============================================================================

//Serialize neutron inelastic interactions so they can be reweighted to a different neutron interaction model later.
void MasterAnaDev::declareInelasticReweightBranches()
{
  const std::string prefix = "neutronInelasticReweight";
  declareIntTruthBranch(prefix + "NPaths", 0);
  declareIntTruthBranch(prefix + "NPoints", 0);

  declareContainerIntTruthBranch(prefix + "PDG");
  declareContainerIntTruthBranch(prefix + "IntCode");
  declareContainerIntTruthBranch(prefix + "IntCodePerSegment");
  declareContainerIntTruthBranch(prefix + "TrackID");
  declareContainerIntTruthBranch(prefix + "ID");
  declareContainerIntTruthBranch(prefix + "ParentID");
  declareContainerIntTruthBranch(prefix + "NTrajPoints");
  declareContainerIntTruthBranch(prefix + "NTrajPointsSaved");
  declareContainerIntTruthBranch(prefix + "TrajPointFlag");
  declareContainerIntTruthBranch(prefix + "Nuke");

  declareContainerIntTruthBranch(prefix + "NInelasticChildren");
  declareContainerIntTruthBranch(prefix + "InelasticChildPDGs");

  declareContainerDoubleTruthBranch(prefix + "InitialE");
  declareContainerDoubleTruthBranch(prefix + "FinalE");
  declareContainerDoubleTruthBranch(prefix + "InitialSigmaE");
  declareContainerDoubleTruthBranch(prefix + "FinalSigmaE");
  declareContainerDoubleTruthBranch(prefix + "ColumnarDensity");
  declareContainerDoubleTruthBranch(prefix + "PosX");
  declareContainerDoubleTruthBranch(prefix + "PosY");
  declareContainerDoubleTruthBranch(prefix + "PosZ");
  declareContainerDoubleTruthBranch(prefix + "TimeDelta");
}

void MasterAnaDev::fillInelasticReweightBranches(Minerva::GenMinInteraction* truth) const
{
  //TODO: Remove me
  warning() << "Starting MasterAnaDev::fillInelasticReweightBranches()...\n";

  const auto rwInfo = m_neutronInelasticReweight->getReweightInfo(truth);

  const std::string prefix = "neutronInelasticReweight";
  truth->setIntData(prefix + "NPaths", rwInfo->nPaths);
  truth->setIntData(prefix + "NPoints", rwInfo->nPoints);

  truth->setContainerIntData(prefix + "PDG", rwInfo->pDG);
  truth->setContainerIntData(prefix + "IntCode", rwInfo->intCode);
  truth->setContainerIntData(prefix + "IntCodePerSegment", rwInfo->intCodePerSegment);
  truth->setContainerIntData(prefix + "TrackID", rwInfo->trackID);
  truth->setContainerIntData(prefix + "ID", rwInfo->id);
  truth->setContainerIntData(prefix + "ParentID", rwInfo->parentID);
  truth->setContainerIntData(prefix + "NTrajPoints", rwInfo->ntrajPoints);
  truth->setContainerIntData(prefix + "NTrajPointsSaved", rwInfo->ntrajPointsSaved);
  truth->setContainerIntData(prefix + "TrajPointFlag", rwInfo->trajPointFlag);
  truth->setContainerIntData(prefix + "Nuke", rwInfo->nuke);

  truth->setContainerIntData(prefix + "NInelasticChildren", rwInfo->nInelasticChildren);
  truth->setContainerIntData(prefix + "InelasticChildPDGs", rwInfo->inelasticChildPDGs);

  truth->setContainerDoubleData(prefix + "InitialE", rwInfo->initialE);
  truth->setContainerDoubleData(prefix + "FinalE", rwInfo->finalE);
  truth->setContainerDoubleData(prefix + "InitialSigmaE", rwInfo->initialSigmaE);
  truth->setContainerDoubleData(prefix + "FinalSigmaE", rwInfo->finalSigmaE);
  truth->setContainerDoubleData(prefix + "ColumnarDensity", rwInfo->columnarDensity);
  truth->setContainerDoubleData(prefix + "PosX", rwInfo->posX);
  truth->setContainerDoubleData(prefix + "PosY", rwInfo->posY);
  truth->setContainerDoubleData(prefix + "PosZ", rwInfo->posZ);
  truth->setContainerDoubleData(prefix + "TimeDelta", rwInfo->timedelta);

  //TODO: Remove me
  warning() << "Finished MasterAnaDev::fillInelasticReweightBranches()\n";
}

//=============================================================================
// Finalize --
//=============================================================================

StatusCode MasterAnaDev::finalize() {

  StatusCode sc = this->MinervaAnalysisTool::finalize();
  if( sc.isFailure() ) { return Error( "Failed to finalize!", sc ); }

  return sc;

}
//#########################################################################################
//
// End MasterAnaDev.cpp
//
//#########################################################################################

#ifndef CVUniverse_cxx
#define CVUniverse_cxx

#include <algorithm> // max_element
#include <cmath> //isfinite

#include "PlotUtils/MnvTuneSystematics.h"

#include "CVUniverse.h"
#include "common_functions.h" // FixAngle


//==============================================================================
// Constructor
//==============================================================================
CVUniverse::CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
  : PlotUtils::MinervaUniverse(chw, nsigma)
{}

//==============================================================================
// Get Branches and Calculate Quantities for the universe/event
//==============================================================================
//==============================
// Muon Variables
// A lot of these are now defined in DCVU, but I'm not yet confident that
// they're being used consistently/correctly. So I'm using mine.
//==============================
// Reco
  double CVUniverse::GetThetamuDeg() const { return ConvertRadToDeg(GetThetamu()); }
  double CVUniverse::GetPZmu() const { return GetDouble("CCNuPionInc_muon_pz"); }
  double CVUniverse::GetPTmu() const {
    return sqrt(pow(GetDouble("CCNuPionInc_muon_px"),2.0) + pow(GetDouble("CCNuPionInc_muon_py"),2.0));
  }

// True
  double CVUniverse::GetPmuTrue() const { 
    return sqrt( pow(GetDouble("truth_mu_px"),2.0) + 
        pow(GetDouble("truth_mu_py"),2.0) + 
        pow(GetDouble("truth_mu_pz"),2.0) );
  }
  double CVUniverse::GetPTmuTrue() const { 
    return sqrt( pow(GetDouble("truth_mu_px"),2.0) +
                 pow(GetDouble("truth_mu_py"),2.0) );
  }
  double CVUniverse::GetPZmuTrue() const { return GetDouble("truth_mu_pz"); }
  double CVUniverse::GetEmuTrue() const { return sqrt( pow(GetPmuTrue(),2.0) + pow(CCNuPionIncConsts::MUON_MASS,2.0) ); }
  double CVUniverse::GetThetamuTrue() const { return FixAngle(GetDouble("truth_mu_theta_wrtbeam"));}
  double CVUniverse::GetThetamuTrueDeg() const {return ConvertRadToDeg(GetThetamuTrue());}

//==============================
// Event-wide Variables
//==============================
// Reco (always MeV, radians)
  double CVUniverse::GetEnu() const { return GetEmu()+GetEhad(); }
  double CVUniverse::GetEhad() const { return GetCalRecoilEnergy(); }
  double CVUniverse::GetCalRecoilEnergy() const { return GetDouble("CCNuPionInc_hadron_recoil_CCInc"); }
  double CVUniverse::GetQ2() const { return CalcQ2(GetEnu(), GetEmu(), GetThetamu()); }
  double CVUniverse::Getq0() const { return Calcq0(GetEnu(), GetEmu()); }
  double CVUniverse::Getq3() const { return Calcq3(GetQ2(), GetEnu(), GetEmu()); }
  double CVUniverse::GetWexp() const { return CalcWexp(GetQ2(), GetEhad()); }

// True (always MeV, radians)
  double CVUniverse::GetEhadTrue() const { return GetEnuTrue() - GetEmuTrue();}
  double CVUniverse::GetWexpTrue() const { return CalcWexp( GetQ2True(), GetEhadTrue() );}
  double CVUniverse::GetWgenie() const { return GetDouble("mc_w"); }

//==============================
// Hadron (Track) Variables
//==============================
// Reco (always MeV, radians)
  double CVUniverse::GetTpi(RecoPionIdx hadron) const { 
    if (hadron == -1) {
      std::cerr << "CVU::GetTpi: pion_idx = -1.\n" 
        "In the future this will be the code for a vertex pion.\n";
      throw hadron;
    }
    //return (GetVecElem("CCNuPionInc_hadron_pion_E", hadron) 
    //          - CCNuPionIncConsts::CHARGED_PION_MASS)/0.96;
    double Epi = GetVecElem("CCNuPionInc_hadron_pion_E", hadron);
    if (Epi == 0) { 
      return GetTpiMBR(hadron); // TODO maybe do from momentum instead
                                // Not sure what this fail mode means.
      //std::cout << "CVUniverse::GetTpi: warning Epi < 0\n"; 
      //return 0;
    }
    return Epi - CCNuPionIncConsts::CHARGED_PION_MASS;
  }


  double CVUniverse::GetTpiMBR(RecoPionIdx hadron) const { 
    if (hadron == -1) {
      std::cerr << "CVU::GetTpiMBR: pion_idx = -1.\n" 
        "In the future this will be the code for a vertex pion.\n"
        "In that case, this function won't make sense.\n";
      throw hadron;
    }
    double TLA = GetVecElem("hadron_track_length_area", hadron);
    return 2.3112*TLA + 37.03;
  }


  double CVUniverse::GetThetapi(RecoPionIdx hadron) const { 
    if (hadron == -1) {
      std::cerr << "CVU::GetThetapi: pion_idx = -1.\n" 
        "In the future this will be the code for a vertex pion.\n"
        "In that case, this function won't make sense.\n";
      throw hadron;
    }
    return FixAngle(GetVecElem("CCNuPionInc_hadron_pion_theta", hadron));
  }


  double CVUniverse::GetThetapiDeg(RecoPionIdx hadron) const { 
    return ConvertRadToDeg(GetThetapi(hadron));
  }


  double CVUniverse::GetLLRScore(RecoPionIdx hadron ) const { 
    if (hadron == -1) {
      std::cout << "CVU::GetLLRScore pion_idx = -1.\n" 
        "In the future this will be the code for a vertex pion.\n"
        "In that case, this function won't make sense.\n";
      throw hadron;
    }
    else if (hadron < -1) {
      std::cerr << "CVU::GetLLRScore bogus pion_idx." << hadron << "\n"; 
      return -1;
      //throw hadron;
    }
    return GetVecElem("CCNuPionInc_hadron_piFit_scoreLLR", hadron); 
  }


  double CVUniverse::GetdEdxScore(RecoPionIdx hadron) const { 
    if (hadron == -1) {
      std::cout << "CVU::GetdEdxScore: pion_idx = -1.\n" 
        "In the future this will be the code for a vertex pion.\n"
        "In that case, this function won't make sense.\n";
      throw hadron;
    }
    return GetVecElem("CCNuPionInc_hadron_piFit_score1", hadron); 
  }

  double CVUniverse::GetEmichel(RecoPionIdx hadron) const {return GetVecElem("has_michel_cal_energy", hadron);}

  double CVUniverse::GetEnode0(RecoPionIdx hadron) const {return GetVecElem("CCNuPionInc_pion_lastnode_Q0", hadron);}
  double CVUniverse::GetEnode1(RecoPionIdx hadron) const {return GetVecElem("CCNuPionInc_pion_lastnode_Q1", hadron);}
  double CVUniverse::GetEnode01(RecoPionIdx hadron) const {return GetEnode0(hadron) + GetEnode1(hadron);}
  double CVUniverse::GetEnode2(RecoPionIdx hadron) const {return GetVecElem("CCNuPionInc_pion_lastnode_Q2", hadron);}
  double CVUniverse::GetEnode3(RecoPionIdx hadron) const {return GetVecElem("CCNuPionInc_pion_lastnode_Q3", hadron);}
  double CVUniverse::GetEnode4(RecoPionIdx hadron) const {return GetVecElem("CCNuPionInc_pion_lastnode_Q4", hadron);}
  double CVUniverse::GetEnode5(RecoPionIdx hadron) const {return GetVecElem("CCNuPionInc_pion_lastnode_Q5", hadron);}

// True (always MeV, radians)
  // Get TRUE Hadron Quantities (always MeV, radians)
  // NOTE: These 'truth_pi_*' containers all have 20 elements
  // Their non-default elements correspond to all of and only the CHARGED
  // (N.B. + & -) pions from the event's TG4Trajectory primary
  // trajectories.
  // Pion indices here do not correspond to reco hadron indices.
  // Be careful with these functions. Make sure you understand the input
  // and output.
  // input: truth index. output: may be pip or pim.
  double CVUniverse::GetTpiTrue(TruePionIdx idx) const {
    double t_pi_E = GetVecElem("truth_pi_E", idx);
    if (t_pi_E == -1.) {
      std::cerr << "CVU::GetTpiTrue: Default energy.\n"
        "Tried to access truth pion energy for a nonexistent "
        "truth pion trajectory.\n";
      throw t_pi_E;
    }
    return  t_pi_E - CCNuPionIncConsts::CHARGED_PION_MASS;
  }

  // !!!!!! Vector filled with pip and pim !!!!!!!
  std::vector<double> CVUniverse::GetTpiTrueVec() const {
    std::vector<double> ret;
    const int n_true_pions = GetNChargedPionsTrue();
    for (TruePionIdx idx = 0; idx<n_true_pions; ++idx) {
      ret.push_back( GetTpiTrue(idx) );
    }
    return ret;
  }


  // input: truth index. output: may be pip or pim.
  double CVUniverse::GetThetapiTrue(TruePionIdx idx) const { 
    double t_pi_theta = GetVecElem("truth_pi_theta_wrtbeam", idx);
    if (t_pi_theta == -9.0) {
      std::cerr << "CVU::GetThetapiTrue: Default angle.\n"
        "Tried to access truth pion angle for a nonexistent "
        "truth pion trajectory.\n";
      throw t_pi_theta;
    }
    return FixAngle(t_pi_theta);
  }

  double CVUniverse::GetThetapiTrueDeg(TruePionIdx idx) const {
    return ConvertRadToDeg(GetThetapiTrue(idx));
  }


  // input: truth index. output: may be pip or pim.
  int CVUniverse::GetPiChargeTrue(TruePionIdx idx) const {
    int t_pi_charge = GetVecElem("truth_pi_charge", idx);
    if (t_pi_charge == 0) {
      std::cerr << "CVU::GetPiChargeTrue: Default charge.\n"
        "Tried to access truth pion charge for a nonexistent "
        "truth pion trajectory.\n";
      throw t_pi_charge;
    }
    return t_pi_charge;
  }

  int CVUniverse::GetNChargedPionsTrue() const {
    return GetInt("truth_N_pip") + GetInt("truth_N_pim");
  }

//==============================
// Misc
//==============================
  double CVUniverse::GetTpiFResidual(const int hadron, const bool MBR) const {
    double T_reco  = !MBR ? GetTpi(hadron) : GetTpiMBR(hadron);
    int true_index = -1;
    true_index = GetVecElem("CCNuPionInc_hadron_tm_trackID", hadron);
    double T_true  = GetVecElem("CCNuPionInc_hadron_tm_beginKE", hadron);
    double fresid = (!std::isfinite(T_reco/T_true) || (std::abs(T_reco/T_true))>10000) ? 
      0.96 : T_reco/T_true - 1.;
    int true_pdg   = GetVecElem("CCNuPionInc_hadron_tm_PDGCode", hadron);
    return fresid;
  }

  double CVUniverse::GetLargestPrimProngSep() const {
    std::vector<double> sep_vec = GetVec<double>("prong_separation");
    return *std::max_element(sep_vec.begin(), sep_vec.end());
  }

  double CVUniverse::GetLargestIsoProngSep() const {
    std::vector<double> sep_vec = GetVec<double>("iso_prong_separation");
    return *std::max_element(sep_vec.begin(), sep_vec.end());
  }

  int CVUniverse::GetNIsoProngs() const { return GetDouble("iso_prongs_count"); }

//==============================
// New Study variables
//==============================
int CVUniverse::GetNAnchoredShortTracks() const { return GetInt("n_anchored_short_trk_prongs"); }
int CVUniverse::GetNAnchoredLongTracks() const { return GetInt("n_anchored_long_trk_prongs"); }
double CVUniverse::GetFitVtxX() const { return GetVecElem("CCNuPionInc_vtx", 0); } // cm?
double CVUniverse::GetFitVtxY() const { return GetVecElem("CCNuPionInc_vtx", 1); } // cm?
double CVUniverse::GetFitVtxZ() const { return GetVecElem("CCNuPionInc_vtx", 2); } // cm? 
int CVUniverse::GetTrackReconstructionMethod(RecoPionIdx hadron) const { return GetVecElem("CCNuPionInc_hadron_1stTrackPatRec", hadron); }
int CVUniverse::GetNNodes(RecoPionIdx hadron) const { return GetVecElem("CCNuPionInc_pion_nNodes", hadron); }

//==============================
// Dummy access for variable constructors
//==============================
  double CVUniverse::GetDummyVar()               const { return -999.; }
  double CVUniverse::GetDummyHadVar(const int x) const { return -999.; }

//==============================================================================
// Calculate Quantities(Always MeV)
//==============================================================================
  double CVUniverse::CalcQ2(const double Enu, const double Emu, const double Thetamu) const {
    double Q2 = 2.0*Enu*( Emu - sqrt( pow(Emu,2.0) - pow(CCNuPionIncConsts::MUON_MASS,2.0) )*cos(Thetamu) ) - pow(CCNuPionIncConsts::MUON_MASS,2.0);
    if ( Q2 < 0. ) Q2 = 0.0;
    return Q2;
  }


  double CVUniverse::CalcWexp(const double Q2, const double Ehad) const {
    double W = pow(CCNuPionIncConsts::PROTON_MASS,2.0) - Q2 + 
      2.0*(CCNuPionIncConsts::PROTON_MASS)*Ehad;
    W = W > 0 ? sqrt(W) : 0.0;
    return W;
  }


  double CVUniverse::Calcq0(const double Enu, const double Emu) const {return Enu-Emu;}


  double CVUniverse::Calcq3(const double Q2, const double Enu, const double Emu) const { return sqrt(Q2 + pow(Enu - Emu,2.0)); }

//==============================================================================
// Get Event Weight
//==============================================================================
  double CVUniverse::GetWeight() const {
    const bool do_warping = false;
    double wgt_flux=1.,   wgt_2p2h=1.;
    double wgt_rpa=1.,    wgt_lowq2=1.;
    double wgt_genie=1.,  wgt_mueff=1.;
    double wgt_anisodd=1.;

    // genie
    wgt_genie = GetGenieWeight();
    //if (do_warping)
    //  wgt_genie = GetGenieWarpWeight();

    // flux
    wgt_flux = GetFluxAndCVWeight();

    // rpa
    wgt_rpa = GetRPAWeight();

    // MINOS efficiency
    if (!m_is_truth && GetBool("isMinosMatchTrack"))
      wgt_mueff = GetMinosEfficiencyWeight();

    //2p2h
    wgt_2p2h = GetLowRecoil2p2hWeight();

    //// low Q2
    //if (do_warping) {
    //  double q2 = GetQ2True();
    //  q2 = q2/1000000.; // pass to function as GeV^2
    //  wgt_lowq2 = GetLowQ2PiWarpWeight(q2, CCNuPionIncShifts::kLowQ2PiChannel); 
    //}

    // aniso delta decay weight -- currently being used for warping
    if (do_warping)
      wgt_anisodd = GetVecElem("truth_genie_wgt_Theta_Delta2Npi",4);

    return wgt_genie*wgt_flux*wgt_2p2h*wgt_rpa*wgt_lowq2*wgt_mueff*wgt_anisodd;
  }

//==============================================================================
// Warping Studies Weights
//==============================================================================
  double CVUniverse::GetGenieWarpWeight() const {
    double wgt = GetVecElem("truth_genie_wgt_MaRES", 4);
    wgt = 1 + (wgt-1)*2; // double the size of the shift from 1 (e.g. 1.1 --> 1.2)
    return wgt;
  }


  double CVUniverse::GetLowQ2PiWarpWeight(double q2, std::string channel) const { 
    if(!PlotUtils::IsCCRes(*this)) 
      return 1.;
    else
      return PlotUtils::weight_lowq2pi().getWeight(q2, channel, +1);
  }


  double CVUniverse::GetAnisoDeltaDecayWarpWeight() const {
    return GetVecElem("truth_genie_wgt_Theta_Delta2Npi",4);
  }

//==============================================================================
// Get pion candidate index corresponding to highest KE pion.
// -2 --> empty input pion candidates vector
// -1 --> vertex pion
// >= 0 --> index of tracked hadron
// ** No pions with KE > 0 --> fail
// ** TODO what happens when MBR and original method disagree
//==============================================================================
int CVUniverse::GetHighestEnergyPionCandidateIndex(const std::vector<int>& pion_candidate_idxs) const {
  if (pion_candidate_idxs.empty()) {
    return CCNuPionIncConsts::kEmptyPionCandidateVector; // == -2
  }

  const int dummy_idx = -99;
  int largest_tpi_idx = dummy_idx;
  double largest_tpi  = 0.;
  for(uint iter = 0; iter < pion_candidate_idxs.size(); ++iter) {
    int current_idx = pion_candidate_idxs[iter];
    double current_tpi = -997.;
    if (current_idx == CCNuPionIncConsts::kIsVertexPion) {
      std::cerr << "GetHighestEnergyPionCandidateIndex: pion_idx = -1.\n" 
                   "In the future this will be the code for a vertex pion.\n";
      std::exit(2);
      //current_tpi = universe.GetVertexTpi(current_idx);
    }
    else{
      current_tpi = GetTpi(current_idx);
      if (current_tpi > largest_tpi) {
        largest_tpi_idx = current_idx;
        largest_tpi = current_tpi;
      }
    }
  } 
  if (largest_tpi_idx == dummy_idx) {
    //return pion_candidate_idxs[0];
    std::cerr << "GetHighestEnergyPionCandidateIndex: no pion with KE > 0!\n";
    return -3;
  }
  return largest_tpi_idx;
}

TruePionIdx CVUniverse::GetHighestEnergyTruePionIndex() const {
  std::vector<double> tpi_vec = GetTpiTrueVec(); // pip and pim
  const int n_true_pions = GetNChargedPionsTrue();
  //if (n_true_pions != tpi_vec.size()) {
  //  std::cerr << "GetHighestEnergyTruePionIndex Error: True pion size mismatch\n";
  //  std::exit(1);
  //}

  TruePionIdx reigning_idx = -1;
  double reigning_tpi = 0;
  for(TruePionIdx idx = 0; idx < n_true_pions; ++idx) {
    if (tpi_vec[idx] > reigning_tpi && GetPiChargeTrue(idx) > 0.) {
      reigning_idx = idx;
      reigning_tpi = tpi_vec[idx];
    }
  }

  // This solution below doesn't account for pion charge.
  // see cache/AlgTests for some better attempts at std::alg/lambda solutions
  // TruePionIdx reigning_idx = distance(tpi_vec, max_element(tpi_vec.begin(), tpi_vec.end()));
  
  return reigning_idx;
}

//==============================================================================
// Print arachne link
//==============================================================================
  void CVUniverse::PrintArachneLink() const {
    int link_size = 200;
    char link [link_size];
    int run    = GetInt("mc_run");
    int subrun = GetInt("mc_subrun");
    int gate   = GetInt("mc_nthEvtInFile")+1;
    int slice  = GetVecElem("slice_numbers", 0);
    sprintf(link, 
        "http://minerva05.fnal.gov/Arachne/arachne.html\?det=SIM_minerva&recoVer=v21r1p1&run=%d&subrun=%d&gate=%d&slice=%d", 
        run, subrun, gate, slice);
    //strncpy(); // FAIL
    //memcpy();  // FAIL
    std::cout << link << std::endl;
  }


#endif // CVUniverse_cxx

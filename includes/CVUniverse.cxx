#ifndef CVUniverse_cxx
#define CVUniverse_cxx

#include "CVUniverse.h"

#include <algorithm>  // max_element
#include <cmath>      //isfinite

#include "PlotUtils/MnvTuneSystematics.h"
#include "utilities.h" // FixAngle

//==============================================================================
// Constructor
//==============================================================================
CVUniverse::CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
    : PlotUtils::MinervaUniverse(chw, nsigma) {}

//==============================================================================
// Get Branches and Calculate Quantities for the universe/event
//==============================================================================
//==============================
// Muon Variables
// A lot of these are now defined in DCVU, but I'm not yet confident that
// they're being used consistently/correctly. So I'm using mine.
//==============================
// Reco
double CVUniverse::GetThetamuDeg() const {
  return ConvertRadToDeg(GetThetamu());
}
double CVUniverse::GetPZmu() const { return GetDouble("MasterAnaDev_muon_Pz"); }
double CVUniverse::GetPTmu() const {
  return sqrt(pow(GetDouble("MasterAnaDev_muon_Px"), 2.0) +
              pow(GetDouble("MasterAnaDev_muon_Py"), 2.0));
}

// True
double CVUniverse::GetPmuTrue() const {
  return sqrt(pow(GetDouble("truth_muon_px"), 2.0) +
              pow(GetDouble("truth_muon_py"), 2.0) +
              pow(GetDouble("truth_muon_pz"), 2.0));
}
double CVUniverse::GetPTmuTrue() const {
  return sqrt(pow(GetDouble("truth_muon_px"), 2.0) +
              pow(GetDouble("truth_muon_py"), 2.0));
}
double CVUniverse::GetPZmuTrue() const { return GetDouble("truth_muon_pz"); }
double CVUniverse::GetEmuTrue() const {
  return sqrt(pow(GetPmuTrue(), 2.0) + pow(CCNuPionIncConsts::MUON_MASS, 2.0));
}
double CVUniverse::GetThetamuTrue() const {
  return FixAngle(GetDouble("truth_muon_theta"));
}
double CVUniverse::GetThetamuTrueDeg() const {
  return ConvertRadToDeg(GetThetamuTrue());
}

//==============================
// Event-wide Variables
//==============================
// Reco (always MeV, radians)
double CVUniverse::GetEnu() const { return GetEmu() + GetEhad(); }
//double CVUniverse::GetEhad() const { return GetDouble("MasterAnaDev_hadron_recoil_CCInc"); }
//double CVUniverse::GetEhad() const { return GetDouble("MasterAnaDev_recoil_E"); }
// double CVUniverse::GetCalRecoilEnergy() const {
//  return GetDouble("MasterAnaDev_hadron_recoil_CCInc");
//}
double CVUniverse::GetQ2() const {
  return CalcQ2(GetEnu(), GetEmu(), GetThetamu());
}
double CVUniverse::Getq0() const { return Calcq0(GetEnu(), GetEmu()); }
double CVUniverse::Getq3() const { return Calcq3(GetQ2(), GetEnu(), GetEmu()); }
double CVUniverse::GetWexp() const { return CalcWexp(GetQ2(), GetEhad()); }

// True (always MeV, radians)
// double CVUniverse::GetEhadTrue() const { return GetEnuTrue() - GetEmuTrue();
// }
double CVUniverse::GetWexpTrue() const {
  return CalcWexp(GetQ2True(), GetEhadTrue());
}
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
  // return (GetVecElem("MasterAnaDev_hadron_pion_E", hadron)
  //          - CCNuPionIncConsts::CHARGED_PION_MASS)/0.96;
  double Epi = GetVecElem("MasterAnaDev_pion_E", hadron);
  if (Epi == 0) {
    return GetTpiMBR(hadron);  // TODO maybe do from momentum instead
                               // Not sure what this fail mode means.
    // std::cout << "CVUniverse::GetTpi: warning Epi < 0\n";
    // return 0;
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
  return 2.3112 * TLA + 37.03;
}

double CVUniverse::GetThetapi(RecoPionIdx hadron) const {
  if (hadron == -1) {
    std::cerr << "CVU::GetThetapi: pion_idx = -1.\n"
                 "In the future this will be the code for a vertex pion.\n"
                 "In that case, this function won't make sense.\n";
    throw hadron;
  }
  return FixAngle(GetVecElem("MasterAnaDev_pion_theta", hadron));
}

double CVUniverse::GetThetapiDeg(RecoPionIdx hadron) const {
  return ConvertRadToDeg(GetThetapi(hadron));
}

double CVUniverse::GetLLRScore(RecoPionIdx hadron) const {
  if (hadron == -1) {
    std::cout << "CVU::GetLLRScore pion_idx = -1.\n"
                 "In the future this will be the code for a vertex pion.\n"
                 "In that case, this function won't make sense.\n";
    throw hadron;
  } else if (hadron < -1) {
    std::cerr << "CVU::GetLLRScore bogus pion_idx." << hadron << "\n";
    return -1;
    // throw hadron;
  }
  return GetVecElem("MasterAnaDev_hadron_piFit_scoreLLR", hadron);
}

double CVUniverse::GetdEdxScore(RecoPionIdx hadron) const {
  if (hadron == -1) {
    std::cout << "CVU::GetdEdxScore: pion_idx = -1.\n"
                 "In the future this will be the code for a vertex pion.\n"
                 "In that case, this function won't make sense.\n";
    throw hadron;
  }
  return GetVecElem("MasterAnaDev_hadron_piFit_score1", hadron);
}

double CVUniverse::GetEmichel(RecoPionIdx hadron) const {
  return GetVecElem("has_michel_cal_energy", hadron);
}

double CVUniverse::GetEnode0(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q0", hadron);
}
double CVUniverse::GetEnode1(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q1", hadron);
}
double CVUniverse::GetEnode01(RecoPionIdx hadron) const {
  return GetEnode0(hadron) + GetEnode1(hadron);
}
double CVUniverse::GetEnode2(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q2", hadron);
}
double CVUniverse::GetEnode3(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q3", hadron);
}
double CVUniverse::GetEnode4(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q4", hadron);
}
double CVUniverse::GetEnode5(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q5", hadron);
}

double CVUniverse::GetPZpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_Pz", hadron);
}
double CVUniverse::GetPXpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_Px", hadron);
}
double CVUniverse::GetPYpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_Py", hadron);
}

double CVUniverse::Gett(RecoPionIdx h) const {
  ROOT::Math::PxPyPzEVector mu4v = GetMuon4V();
  return Calct(GetPXpi(h), GetPYpi(h), GetPZpi(h), GetEpi(h), mu4v.Px(),
               mu4v.Py(), mu4v.Pz(), GetEmu());
}

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
  return t_pi_E - CCNuPionIncConsts::CHARGED_PION_MASS;
}

// !!!!!! Vector filled with pip and pim !!!!!!!
std::vector<double> CVUniverse::GetTpiTrueVec() const {
  std::vector<double> ret;
  const int n_true_pions = GetNChargedPionsTrue();
  for (TruePionIdx idx = 0; idx < n_true_pions; ++idx) {
    ret.push_back(GetTpiTrue(idx));
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
// Ehad Variables
//==============================
// assumes pion
double CVUniverse::GetEpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_E", hadron);
}

// AKA GetErecoil
// If the calorimetric energy is too great, abandon trying to calculate tracked
// energy separately. This requires a coordinated effort from both functions.
double CVUniverse::GetEhad() const {
  return GetCalRecoilEnergy() + GetTrackRecoilEnergy();
}

// Untracked recoil energy
double CVUniverse::GetCalRecoilEnergy() const {
  const double ecal_nopi = GetCalRecoilEnergyNoPi_DefaultSpline();
  if (ecal_nopi > 1000)
    return GetCalRecoilEnergy_CCPiSpline();
  else
    return GetCalRecoilEnergyNoPi_Corrected(ecal_nopi);
}

// (Tracked) recoil energy, not determined from calorimetry
double CVUniverse::GetTrackRecoilEnergy() const {
  //if (GetPionCandidates().empty())
  //  std::cout << "CVU::GetETrackedRecoilEnergy WARNING: no pion candidates!\n";

  double etracks = 0.;

  const double ecal_nopi = GetCalRecoilEnergyNoPi_DefaultSpline();
  if (ecal_nopi > 1000) {
    etracks = 0;
  } else {
    for (const auto& pi_idx : GetPionCandidates()) {
      // std::cout << GetVecElem("MasterAnaDev_hadron_tm_PDGCode", pi_idx) << "
      // ";
      etracks += GetEpi(pi_idx);
    }
    // std::cout << "\n";
  }

  return etracks;
}

// Cal recoil energy minus calorimetrically-measured pion energy.
// Used to determined whether we should try to use the correction or not.
double CVUniverse::GetCalRecoilEnergyNoPi_DefaultSpline() const {
  double nopi_recoilE = GetCalRecoilEnergy_DefaultSpline();
  for (const auto& pi_idx : GetPionCandidates()) {
    nopi_recoilE -= GetCalEpi(pi_idx);
  }
  return nopi_recoilE;
}

// Apply an additive, ad hoc correction to the CalRecoilENoPi
double CVUniverse::GetCalRecoilEnergyNoPi_Corrected(
    const double ecal_nopi) const {
  double ecal_nopi_corrected = ecal_nopi;

  RecoPionIdx best_pion =
      GetHighestEnergyPionCandidateIndex(GetPionCandidates());
  if (Gett(best_pion) < 125.e3) return ecal_nopi_corrected;

  TArrayD ecal_nopi_bins = CCPi::GetBinning("ecal_nopi");

  const std::vector<double> corrections = {-0.060, -0.050, -0.210, -0.180,
                                           -0.165, -0.180, -0.180, -0.180,
                                           -0.195, -0.360, -0.400};

  // Find the correction corresponding to the nominal ehad-nopi.
  // Set the corrected value to be the (nominal - correction) (where the
  // corrections are all negative, so in effect we're shifting everything
  // UP)
  for (int i_bin = 0; i_bin < ecal_nopi_bins.GetSize() - 1; ++i_bin) {
    if (((1e3) * ecal_nopi_bins[i_bin] < ecal_nopi) &&
        (ecal_nopi < (1e3) * ecal_nopi_bins[i_bin + 1])) {
      if (ecal_nopi < 1000) {
        ecal_nopi_corrected = ecal_nopi - (1e3) * corrections[i_bin];
        break;
      }
    }
  }

  return ecal_nopi_corrected;
}

// RecoilUtils->calcRecoilEFromClusters(event, muonProng, "Default" );
double CVUniverse::GetCalRecoilEnergy_DefaultSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_default");
}  //

// Total recoil with CCPi spline correction.
// Spline measured from: CC, 1pi+, 1mu, NBaryons,
// True W exp < 1.4 GeV, thmu_true < 20 deg, 35 < tpi < 350 MeV, Minos match
// RecoilUtils->calcRecoilEFromClusters(event, muonProng,
// "NukeCCPion_TwoTrack_Nu_Tracker");
double CVUniverse::GetCalRecoilEnergy_CCPiSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_two_track");
}

// Passive-corrected pion energy
// NukeCCPion_pion_recoilE_passive =
// m_caloUtils->applyCalConsts(hadronProng,"Default",true,true);
double CVUniverse::GetCalEpi(int iProng) const {
  return GetVecElem("MasterAnaDev_hadron_pion_E_recoil_corr", iProng);
}

// Ehad CCInclusive Spline Variables
// Ehad ccinclusive splines -- doesn't account for pion
double CVUniverse::GetCalRecoilEnergy_CCIncSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_CCInc");
}
// Old-school attempt at Ehad - Epi
// Subtract dEdx-tool-pion-energy from Ehad-from-calcRecoilEFromClusters
// ("CCInclusive" splines)
double CVUniverse::GetCalRecoilEnergyNoPi_CCIncSpline() const {
  double ecal_nopi = GetCalRecoilEnergy_CCIncSpline();
  for (const auto& pi_idx : GetPionCandidates()) ecal_nopi -= GetEpi(pi_idx);
  return ecal_nopi;
}

// Ehad truth variables
double CVUniverse::GetEhadTrue() const { return GetEnuTrue() - GetElepTrue(); }

// Given a reco track, get the track's truth KE
double CVUniverse::GetTpiTrueMatched(RecoPionIdx hadron) const {
  // TruePionIdx true_index = -1;
  // TruePionIdx true_index = GetVecElem("MasterAnaDev_hadron_tm_trackID",
  // hadron); int true_pdg   = GetVecElem("MasterAnaDev_hadron_tm_PDGCode",
  // hadron); if(true_pdg != 211) std::cout << "pion mis-identified as a " <<
  // true_pdg << "!\n";
  return GetVecElem("MasterAnaDev_hadron_tm_beginKE", hadron);
}

// Given a reco track, get the track's truth E, assuming it's a pion
double CVUniverse::GetEpiTrueMatched(RecoPionIdx hadron) const {
  return GetTpiTrueMatched(hadron) + CCNuPionIncConsts::CHARGED_PION_MASS;
}

// Need a function that gets truth tracked energy, total E for pions, just KE
// for protons.
double CVUniverse::GetAllTrackEnergyTrue() const {
  double etracks = 0;
  for (const auto& pi_idx : GetPionCandidates()) {
    etracks += GetVecElem("MasterAnaDev_hadron_tm_beginKE", pi_idx);
    // std::cout << GetVecElem("MasterAnaDev_hadron_tm_PDGCode", pi_idx) << " ";
    if (abs(GetVecElem("MasterAnaDev_hadron_tm_PDGCode", pi_idx)) ==
        211)  // TODO may want to only not add pion mass when is proton or
              // neutron
      etracks += CCNuPionIncConsts::CHARGED_PION_MASS;
  }
  // std::cout << "\n";
  return etracks;
}

// Given reco pions – subtract their true matched energy from truth ehad = (enu
// - emu)
double CVUniverse::GetCalRecoilEnergyNoPiTrue() const {
  double nopi_recoilE = GetEhadTrue();
  for (const auto& pi_idx : GetPionCandidates())
    nopi_recoilE -= GetEpiTrueMatched(pi_idx);
  return nopi_recoilE;
}

//==============================
// Misc
//==============================
double CVUniverse::GetTpiFResidual(const int hadron, const bool MBR) const {
  double T_reco = !MBR ? GetTpi(hadron) : GetTpiMBR(hadron);
  int true_index = -1;
  true_index = GetVecElem("MasterAnaDev_hadron_tm_trackID", hadron);
  double T_true = GetVecElem("MasterAnaDev_hadron_tm_beginKE", hadron);
  double fresid =
      (!std::isfinite(T_reco / T_true) || (std::abs(T_reco / T_true)) > 10000)
          ? 0.96
          : T_reco / T_true - 1.;
  int true_pdg = GetVecElem("MasterAnaDev_hadron_tm_PDGCode", hadron);
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

double CVUniverse::GetWexpFResidual() const {
  return GetWexp() / GetWexpTrue() - 1.;
  // return GetWexp()/GetWgenie() - 1.;
}

//==============================
// New Study variables
//==============================
int CVUniverse::GetNAnchoredShortTracks() const {
  return GetInt("n_anchored_short_trk_prongs");
}
int CVUniverse::GetNAnchoredLongTracks() const {
  return GetInt("n_anchored_long_trk_prongs");
}
double CVUniverse::GetFitVtxX() const {
  return GetVecElem("MasterAnaDev_vtx", 0);
}  // cm?
double CVUniverse::GetFitVtxY() const {
  return GetVecElem("MasterAnaDev_vtx", 1);
}  // cm?
double CVUniverse::GetFitVtxZ() const {
  return GetVecElem("MasterAnaDev_vtx", 2);
}  // cm?
int CVUniverse::GetTrackReconstructionMethod(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_hadron_1stTrackPatRec", hadron);
}
int CVUniverse::GetNNodes(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_nNodes", hadron);
}

//==============================
// Dummy access for variable constructors
//==============================
double CVUniverse::GetDummyVar() const { return -999.; }
double CVUniverse::GetDummyHadVar(const int x) const { return -999.; }

//==============================================================================
// Calculate Quantities(Always MeV)
//==============================================================================
double CVUniverse::CalcQ2(const double Enu, const double Emu,
                          const double Thetamu) const {
  double Q2 =
      2.0 * Enu *
          (Emu - sqrt(pow(Emu, 2.0) - pow(CCNuPionIncConsts::MUON_MASS, 2.0)) *
                     cos(Thetamu)) -
      pow(CCNuPionIncConsts::MUON_MASS, 2.0);
  if (Q2 < 0.) Q2 = 0.0;
  return Q2;
}

double CVUniverse::CalcWexp(const double Q2, const double Ehad) const {
  double W = pow(CCNuPionIncConsts::PROTON_MASS, 2.0) - Q2 +
             2.0 * (CCNuPionIncConsts::PROTON_MASS)*Ehad;
  W = W > 0 ? sqrt(W) : 0.0;
  return W;
}

double CVUniverse::Calcq0(const double Enu, const double Emu) const {
  return Enu - Emu;
}

double CVUniverse::Calcq3(const double Q2, const double Enu,
                          const double Emu) const {
  return sqrt(Q2 + pow(Enu - Emu, 2.0));
}

double CVUniverse::Calct(const double pxpi, const double pypi,
                         const double pzpi, const double epi, const double pxmu,
                         const double pymu, const double pzmu,
                         const double emu) const {
  return pow((epi + emu - pzmu - pzpi), 2.0) + pow(pxpi + pxmu, 2.0) +
         pow(pypi + pymu, 2.0);
}
//==============================================================================
// Functions to make fidvol cuts
//==============================================================================
  bool CVUniverse::leftlinesCut (const double a,const double x,const double y) const {
        double b, yls, yli;
        b = a*(2*sqrt(3)/3);
        yls = (sqrt(3)/3)*x + b;
        yli = -(sqrt(3)/3)*x - b;
        if (y > yli && y < yls) return true;
        else return false;
  } 
  bool CVUniverse::rightlinesCut (const double a,const double x,const double y) const {
        double b, yls, yli;
        b = a*(2*sqrt(3)/3);
        yls = -(sqrt(3)/3)*x + b;
        yli = (sqrt(3)/3)*x - b;
        if (y > yli && y < yls) return true;
        else return false;
  }



//==============================================================================
// Get Event Weight
//==============================================================================
double CVUniverse::GetWeight() const {
  const bool do_warping = false;
  double wgt_flux = 1., wgt_2p2h = 1.;
  double wgt_rpa = 1., wgt_lowq2 = 1.;
  double wgt_genie = 1., wgt_mueff = 1.;
  double wgt_anisodd = 1.;

  // genie
  wgt_genie = GetGenieWeight();
  // if (do_warping)
  //  wgt_genie = GetGenieWarpWeight();

  // flux
  wgt_flux = GetFluxAndCVWeight();

  // rpa
  wgt_rpa = GetRPAWeight();

  // MINOS efficiency
  if (!m_is_truth && GetBool("isMinosMatchTrack"))
    wgt_mueff = GetMinosEfficiencyWeight();

  // 2p2h
  wgt_2p2h = GetLowRecoil2p2hWeight();

  //// low Q2
  // if (do_warping) {
  //  double q2 = GetQ2True();
  //  q2 = q2/1000000.; // pass to function as GeV^2
  //  wgt_lowq2 = GetLowQ2PiWarpWeight(q2, CCNuPionIncShifts::kLowQ2PiChannel);
  //}

  // aniso delta decay weight -- currently being used for warping
  if (do_warping)
    wgt_anisodd = GetVecElem("truth_genie_wgt_Theta_Delta2Npi", 4);

  return wgt_genie * wgt_flux * wgt_2p2h * wgt_rpa * wgt_lowq2 * wgt_mueff *
         wgt_anisodd;
}

//==============================================================================
// Warping Studies Weights
//==============================================================================
double CVUniverse::GetGenieWarpWeight() const {
  double wgt = GetVecElem("truth_genie_wgt_MaRES", 4);
  wgt = 1 + (wgt - 1) *
                2;  // double the size of the shift from 1 (e.g. 1.1 --> 1.2)
  return wgt;
}

double CVUniverse::GetLowQ2PiWarpWeight(double q2, std::string channel) const {
  if (!PlotUtils::IsCCRes(*this))
    return 1.;
  else
    return PlotUtils::weight_lowq2pi().getWeight(q2, channel, +1);
}

double CVUniverse::GetAnisoDeltaDecayWarpWeight() const {
  return GetVecElem("truth_genie_wgt_Theta_Delta2Npi", 4);
}

//==============================================================================
// Get pion candidate index corresponding to highest KE pion.
// -2 --> empty input pion candidates vector
// -1 --> vertex pion
// >= 0 --> index of tracked hadron
// ** No pions with KE > 0 --> fail
// ** TODO what happens when MBR and original method disagree
//==============================================================================
int CVUniverse::GetHighestEnergyPionCandidateIndex(
    const std::vector<int>& pion_candidate_idxs) const {
  if (pion_candidate_idxs.empty()) {
    return CCNuPionIncConsts::kEmptyPionCandidateVector;  // == -2
  }

  const int dummy_idx = -99;
  int largest_tpi_idx = dummy_idx;
  double largest_tpi = 0.;
  for (uint iter = 0; iter < pion_candidate_idxs.size(); ++iter) {
    int current_idx = pion_candidate_idxs[iter];
    double current_tpi = -997.;
    if (current_idx == CCNuPionIncConsts::kIsVertexPion) {
      std::cerr << "GetHighestEnergyPionCandidateIndex: pion_idx = -1.\n"
                   "In the future this will be the code for a vertex pion.\n";
      std::exit(2);
      // current_tpi = universe.GetVertexTpi(current_idx);
    } else {
      current_tpi = GetTpi(current_idx);
      if (current_tpi > largest_tpi) {
        largest_tpi_idx = current_idx;
        largest_tpi = current_tpi;
      }
    }
  }
  if (largest_tpi_idx == dummy_idx) {
    // return pion_candidate_idxs[0];
    std::cerr << "GetHighestEnergyPionCandidateIndex: no pion with KE > 0!\n";
    return -3;
  }
  return largest_tpi_idx;
}

TruePionIdx CVUniverse::GetHighestEnergyTruePionIndex() const {
  std::vector<double> tpi_vec = GetTpiTrueVec();  // pip and pim
  const int n_true_pions = GetNChargedPionsTrue();
  // if (n_true_pions != tpi_vec.size()) {
  //  std::cerr << "GetHighestEnergyTruePionIndex Error: True pion size
  //  mismatch\n"; std::exit(1);
  //}

  TruePionIdx reigning_idx = -1;
  double reigning_tpi = 0;
  for (TruePionIdx idx = 0; idx < n_true_pions; ++idx) {
    if (tpi_vec[idx] > reigning_tpi && GetPiChargeTrue(idx) > 0.) {
      reigning_idx = idx;
      reigning_tpi = tpi_vec[idx];
    }
  }

  // This solution below doesn't account for pion charge.
  // see cache/AlgTests for some better attempts at std::alg/lambda solutions
  // TruePionIdx reigning_idx = distance(tpi_vec, max_element(tpi_vec.begin(),
  // tpi_vec.end()));

  return reigning_idx;
}

//==============================================================================
// Print arachne link
//==============================================================================
void CVUniverse::PrintArachneLink() const {
  int link_size = 200;
  char link[link_size];
  int run = GetInt("mc_run");
  int subrun = GetInt("mc_subrun");
  int gate = GetInt("mc_nthEvtInFile") + 1;
  int slice = GetVecElem("slice_numbers", 0);
  sprintf(link,
          "http://minerva05.fnal.gov/Arachne/"
          "arachne.html\?det=SIM_minerva&recoVer=v21r1p1&run=%d&subrun=%d&gate="
          "%d&slice=%d",
          run, subrun, gate, slice);
  // strncpy(); // FAIL
  // memcpy();  // FAIL
  std::cout << link << std::endl;
}

//==============================================================================
// Get and set pion candidates
//==============================================================================
void CVUniverse::SetPionCandidates(std::vector<RecoPionIdx> c) {
  m_pion_candidates = c;
}

std::vector<RecoPionIdx> CVUniverse::GetPionCandidates() const {
  return m_pion_candidates;
}

#endif  // CVUniverse_cxx

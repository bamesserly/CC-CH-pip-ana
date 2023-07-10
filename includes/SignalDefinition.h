#ifndef SignalDefinition_H
#define SignalDefinition_H

#include "includes/CVUniverse.h"
//#include "includes/Constants.h"  // namespace CCNuPionIncConsts

// Count the number of signal definitions we've made so far
static int kNSigDefs = 0;

class SignalDefinition {
 public:
  // Key analysis decision as enums
  enum class PionReco {
    kTracked,
    kUntracked,
    kTrackedAndUntracked,
    kNPionRecoTypes
  };
  enum class WRegion { k1_4, k1_8, kSIS, kNoW, kNWValueTypes };
  enum class NPions { kOnePi, kNPi, kNPiTypes };
  enum class Thetamu { kTwentyDeg, kThirteenDeg, kNThetamuTypes };

  static const std::map<PionReco, double> kTpiMinValues;
  static const std::map<WRegion, double> kWMinValues;
  static const std::map<WRegion, double> kWMaxValues;
  static const std::map<NPions, int> kNPiMaxValues;
  static const std::map<Thetamu, double> kThetamuMaxValues;

  // Enum key analysis decisions are members
  const PionReco m_pion_reco;
  const WRegion m_w_region;
  const NPions m_n_pions;
  const Thetamu m_thetamu;

  // CTOR -- init key analysis decision enums and set min/max values
  SignalDefinition(const PionReco pr, const WRegion wv, const NPions np,
                   const Thetamu tm = Thetamu::kTwentyDeg)
      : m_pion_reco(pr),
        m_w_region(wv),
        m_n_pions(np),
        m_thetamu(tm),
        m_tpi_min(kTpiMinValues.at(m_pion_reco)),
        m_w_min(kWMinValues.at(m_w_region)),
        m_w_max(kWMaxValues.at(m_w_region)),
        m_n_pi_max(kNPiMaxValues.at(m_n_pions)),
        m_thetamu_max(kThetamuMaxValues.at(m_thetamu)),
        m_id(kNSigDefs){
    kNSigDefs++;
  };

  // Determined at initialization
  const double m_tpi_min;      // MeV
  const double m_w_min;        // MeV
  const double m_w_max;        // MeV
  const int m_n_pi_max;
  const double m_thetamu_max;  // rad

  // const signal definition values
  const int m_n_pi_min = 1;
  const double m_tpi_max = 350.;         // MeV
  const int m_IsoProngCutVal = 2;        // strictly fewer than
  const double m_PmuMinCutVal = 1500.;   // MeV/c
  const double m_PmuMaxCutVal = 20000.;  // MeV/c
  const double m_ZVtxMinCutVal = 5990.;  // cm
  const double m_ZVtxMaxCutVal = 8340.;  // cm
  const double m_ApothemCutVal = 850.;   // cm

  const int m_id;
};

// Signal Definition Constants
using PionReco = SignalDefinition::PionReco;
using WRegion = SignalDefinition::WRegion;
using NPions = SignalDefinition::NPions;
using Thetamu = SignalDefinition::Thetamu;

const std::map<PionReco, double> SignalDefinition::kTpiMinValues{
    {PionReco::kTracked, 35.}, {PionReco::kUntracked, 0.}, {PionReco::kTrackedAndUntracked, 0.}};

const std::map<WRegion, double> SignalDefinition::kWMinValues{
    {WRegion::k1_4, 0.}, {WRegion::k1_8, 0.}, {WRegion::kSIS, 1400.}, {WRegion::kNoW, 0.}};  // MeV

const std::map<WRegion, double> SignalDefinition::kWMaxValues{
    {WRegion::k1_4, 1400.}, {WRegion::k1_8, 1800.}, {WRegion::kSIS, 1800.}, {WRegion::kNoW, 9999.}};  // MeV

const std::map<NPions, int> SignalDefinition::kNPiMaxValues{{NPions::kOnePi, 1}, {NPions::kNPi, 99}};

const std::map<Thetamu, double> SignalDefinition::kThetamuMaxValues{
    {Thetamu::kTwentyDeg, 0.3491}, {Thetamu::kThirteenDeg, 0.226892803}};  // radians


// Make some signal definitions
static const SignalDefinition kOnePi        ( PionReco::kTrackedAndUntracked, WRegion::k1_4, NPions::kOnePi);
static const SignalDefinition kOnePiTracked ( PionReco::kTracked,             WRegion::k1_4, NPions::kOnePi);
static const SignalDefinition kOnePiNoW     ( PionReco::kTrackedAndUntracked, WRegion::kNoW, NPions::kOnePi);
static const SignalDefinition kNPi          ( PionReco::kTrackedAndUntracked, WRegion::k1_8, NPions::kNPi);
static const SignalDefinition kNPiNoW       ( PionReco::kTrackedAndUntracked, WRegion::kNoW, NPions::kNPi);
static const SignalDefinition kNuke         ( PionReco::kTracked,             WRegion::k1_4, NPions::kOnePi, Thetamu::kThirteenDeg);

// Map int to SignalDefinition to we can pass by command line
static const std::map<int,SignalDefinition> SignalDefinitionMap{{kOnePi.m_id,kOnePi},{kOnePiTracked.m_id,kOnePiTracked},{kNuke.m_id,kNuke}};

// Truth topology particle counts
// From Aaron
std::map<string, int> GetParticleTopology(
    const std::vector<int>& FS_PDG, const std::vector<double>& FS_energy, const SignalDefinition sig_def) {
  std::map<std::string, int> genie_n;

  // Overarching categories: nucleons, mesons
  genie_n["muons"] = 0;     // Muons, photons (do we want electrons...)
  genie_n["photons"] = 0;   // Photons are filled if there are electrons
  genie_n["pi_zeros"] = 0;  // Pions
  genie_n["piplus"] = 0;
  genie_n["piplus_range"] = 0;  // Pi+ that passes the kinematic cuts
  genie_n["piminus"] = 0;
  genie_n["pions"] = 0;
  genie_n["kaons"] = 0;  // Other mesons
  genie_n["charms"] = 0;
  genie_n["mesons"] = 0;
  genie_n["protons"] = 0;  // Baryons
  genie_n["neutrons"] = 0;
  genie_n["nucleons"] = 0;
  genie_n["heavy_baryons"] = 0;
  genie_n["nuclei"] = 0;
  genie_n["others"] = 0;  // others

  double tpi = 0;
  for (uint p = 0; p < FS_PDG.size(); p++) {
    // Get tpi
    tpi = (FS_energy[p] - MinervaUnits::M_pion);

    // Look at every particle's pdg
    switch (FS_PDG[p]) {
      case 13:
        genie_n["muons"]++;
        break;
      case 22:
        // Check the energy.  If below 10 MeV, it's a nuclear deexcitation,
        // which we ignore
        if (FS_energy[p] > 10) genie_n["photons"]++;
        break;
      case 11:
      case -11:
        genie_n["photons"]++;
        break;
      case 211:
        if (sig_def.m_tpi_min< tpi &&
            tpi < sig_def.m_tpi_max)
          genie_n["piplus_range"]++;
        genie_n["piplus"]++;
        genie_n["pions"]++;
        genie_n["mesons"]++;
        break;
      case -211:
        genie_n["piminus"]++;
        genie_n["pions"]++;
        genie_n["mesons"]++;
        break;
      case 111:
        genie_n["pi_zeros"]++;
        // genie_n["photons"]++;
        genie_n["pions"]++;
        genie_n["mesons"]++;
        break;
      case 130:
      case 310:
      case 311:
      case 321:
      case -130:
      case -310:
      case -311:
      case -321:
        genie_n["kaons"]++;
        genie_n["mesons"]++;
        break;
      case 411:
      case 421:
      case 431:
      case -411:
      case -421:
      case -431:
        genie_n["charms"]++;
        genie_n["mesons"]++;
        break;
      case 2212:
        genie_n["protons"]++;
        genie_n["nucleons"]++;
        break;
      case 2112:
        genie_n["neutrons"]++;
        genie_n["nucleons"]++;
        break;
      case 2000000101:  // bindino is binding energy placeholder
        break;
      default:
        if (FS_PDG[p] > 3000 && FS_PDG[p] < 5000)
          genie_n["heavy_baryons"]++;
        else if (FS_PDG[p] > 1000000000 && FS_PDG[p] < 1099999999)
          genie_n["nuclei"]++;
        else
          genie_n["others"]++;
    }
  }

  return genie_n;
}

// Given the particle topology of GetParticleTopology, is this signal?
// From Aaron
bool Is1PiPlus(const std::map<std::string, int>& particles) {
  //  if( incoming =! 14 || current =! 1 ) return false;

  // which particles are part of your signal (those you count, those you want
  // any) Don't need two variables, just here to keep bookkeeping straight
  TString counted_signal("muons piplus piplus_range pions mesons");
  TString any_baryon("protons neutrons nucleons heavy_baryons nuclei");

  // Count bkgd particles
  int n_bg = 0;
  for (auto p : particles) {
    std::string particle = p.first;
    int count = p.second;
    if (counted_signal.Contains(particle) || any_baryon.Contains(particle))
      continue;
    n_bg += count;
  }

  // SIGNAL
  // 1 mu, 1 pi+, no other mesons
  if (particles.at("muons") == 1 && particles.at("piplus") == 1 &&
      particles.at("pions") == particles.at("piplus") &&
      particles.at("mesons") == particles.at("piplus") && n_bg == 0)
    return true;  // any number of baryons, nothing else

  return false;
}

// Number of abs(pdg) == 211 true TG4Trajectories which also:
// (1) are pip, (2) satisfy a KE restriction
int NSignalPions(const CVUniverse& univ, const SignalDefinition sig_def) {
  int n_signal_pions = 0;
  int n_true_pions = univ.GetNChargedPionsTrue();
  for (TruePionIdx idx = 0; idx < n_true_pions; ++idx) {
    double t_pi = univ.GetTpiTrue(idx);
    double theta_pi = univ.GetThetapiTrue(idx);
    if (univ.GetPiChargeTrue(idx) > 0 &&
        t_pi > sig_def.m_tpi_min &&
        t_pi < sig_def.m_tpi_max
        //&& (theta_pi < 1.39626 || 1.74533 < theta_pi))
    )
      ++n_signal_pions;
  }
  return n_signal_pions;
}

int NOtherParticles(const CVUniverse& univ) {
  int n_other_particles = 0;
  n_other_particles = univ.GetInt("truth_N_chargedK") +
                      univ.GetInt("truth_N_K0") + univ.GetInt("truth_N_sigma") +
                      univ.GetInt("truth_N_lambda");
  return n_other_particles;
}

bool ZVtxIsSignal(const CVUniverse& univ, const SignalDefinition sig_def) {
  double vtx_z = univ.GetVecElem("mc_vtx", 2);
  return sig_def.m_ZVtxMinCutVal < vtx_z &&
                 vtx_z < sig_def.m_ZVtxMaxCutVal
             ? true
             : false;
}

bool XYVtxIsSignal(const CVUniverse& univ, const SignalDefinition sig_def) {
  return univ.IsInHexagon(univ.GetVecElem("mc_vtx", 0),  // x
                          univ.GetVecElem("mc_vtx", 1),  // y
                          sig_def.m_ApothemCutVal);
}

bool IsSignal(const CVUniverse& univ, SignalDefinition sig_def = kOnePi) {
  int n_signal_pions = NSignalPions(univ, sig_def);

  const std::map<std::string, int> particles = GetParticleTopology(
      univ.GetVec<int>("mc_FSPartPDG"), univ.GetVec<double>("mc_FSPartE"), sig_def);

  // TODO switch the pion multiplicity check to use the particles variable
  // TODO Is1PiPlus obviously isn't checking for any number of pions
  return univ.GetInt("mc_current") == 1              &&
      univ.GetInt("mc_incoming") == 14               &&
      univ.GetBool("truth_is_fiducial")              &&
      ZVtxIsSignal(univ, sig_def)                    &&
      XYVtxIsSignal(univ, sig_def)                   &&
      univ.GetThetalepTrue() < sig_def.m_thetamu_max &&
      sig_def.m_w_min < univ.GetWexpTrue()           &&
      univ.GetWexpTrue() < sig_def.m_w_max           &&
      particles.at("piplus_range") == 1              &&
      Is1PiPlus(particles)                           &&
      sig_def.m_PmuMinCutVal < univ.GetPmuTrue()     &&
      univ.GetPmuTrue() < sig_def.m_PmuMaxCutVal     &&
      sig_def.m_n_pi_min <= n_signal_pions           &&
      n_signal_pions <= sig_def.m_n_pi_max           &&
      univ.GetInt("truth_N_pi0") == 0                &&
      univ.GetInt("truth_N_pim") == 0;
}

std::string GetSignalName(SignalDefinition sig_def) {
  switch (sig_def.m_id) {
    case kOnePi.m_id:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X  (W < 1.4 GeV)";
    case kOnePiTracked.m_id:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X  (W < 1.4 GeV, tracked)";
    case kOnePiNoW.m_id:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X";
    case kNPi.m_id:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X  (W < 1.8 GeV)";
    case kNPiNoW.m_id:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X";
    default:
      return "UNKNOWN SIGNAL";
  }
}

std::string GetSignalFileTag(SignalDefinition sig_def) {
  switch (sig_def.m_id) {
    case kOnePi.m_id:
      return "1Pi";
    case kOnePiTracked.m_id:
      return "1PiTracked";
    case kOnePiNoW.m_id:
      return "1PiNoW";
    case kNPi.m_id:
      return "NPi";
    case kNPiNoW.m_id:
      return "NPiNoW";
    default:
      return "UNKNOWN SIGNAL";
  }
}

#endif  // Signal_H

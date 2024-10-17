#ifndef CVUniverse_H
#define CVUniverse_H

#include <TVector3.h>

#include "Binning.h"    // CCPi::GetBinning for ehad_nopi
#include "Constants.h"  // CCNuPionIncConsts, CCNuPionIncShifts, Reco/TruePionIdx
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/LowRecoilPionReco.h"
#include "PlotUtils/MinervaUniverse.h"

class CVUniverse : public PlotUtils::MinervaUniverse {
 private:
  // Pion Candidates - clear these when SetEntry is called
  std::vector<RecoPionIdx> m_pion_candidates;
  LowRecoilPion::MichelEvent<CVUniverse> m_vtx_michels;

 public:
#include "PlotUtils/LowRecoilPionFunctions.h"
#include "PlotUtils/MichelFunctions.h"
#include "PlotUtils/MuonFunctions.h"
#include "PlotUtils/RecoilEnergyFunctions.h"
#include "PlotUtils/TruthFunctions.h"
#include "PlotUtils/WeightFunctions.h"
  // CTOR
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0);

  // DTOR
  virtual ~CVUniverse(){};

  // Print arachne link
  void PrintArachneLink() const;

  // Dummy access for variable constructors
  virtual double GetDummyVar() const;
  virtual double GetDummyHadVar(const int x) const;

  bool m_passesTrackedCuts;
  bool m_passesTracklessCuts;
  bool m_passesTrackedSideband;
  bool m_passesTracklessSideband;
  bool m_passesTrackedExceptW;
  bool m_passesTracklessExceptW;

  bool m_is_signal;

  // No stale cache!
  virtual void OnNewEntry() override {
    m_pion_candidates.clear();
    m_vtx_michels = LowRecoilPion::MichelEvent<CVUniverse>();
    assert(m_vtx_michels.m_idx == -1);
    m_passesTrackedCuts = false;
    m_passesTracklessCuts = false;
    m_passesTrackedSideband = false;
    m_passesTracklessSideband = false;
    m_passesTrackedExceptW = false;
    m_passesTracklessExceptW = false;
    m_is_signal = false;
  }

  virtual bool IsVerticalOnly() const override { return true; }

  bool IsSignal() { return m_is_signal; }
  void SetIsSignal(bool is_signal) { m_is_signal = is_signal; }
  // Get and set pion candidates
  TruePionIdx GetHighestEnergyTruePionIndex() const;
  int GetHighestEnergyPionCandidateIndex(const std::vector<int>& pions) const;
  std::vector<RecoPionIdx> GetPionCandidates() const;
  void SetPionCandidates(std::vector<RecoPionIdx> c);
  void SetVtxMichels(const LowRecoilPion::MichelEvent<CVUniverse>& m) {
    m_vtx_michels = m;
  }
  LowRecoilPion::MichelEvent<CVUniverse> GetVtxMichels() const {
    return m_vtx_michels;
  }
  void SetPassesTrakedTracklessCuts(
      bool passesTrackedCuts, bool passesTracklessCuts, bool tracked_sideband,
      bool trackless_sideband, bool tracked_all_ex_w, bool trackless_all_ex_w);

  virtual double ApplyCaloTuning(const double& cal_recoil_energy) const {
    return cal_recoil_energy;
  }

  //==============================================================================
  // Analysis Variables
  //==============================================================================
  // muon
  virtual double GetPTmu() const;
  virtual double GetPXmu() const;
  virtual double GetPYmu() const;
  virtual double GetPZmu() const;
  virtual double GetThetamuDeg() const;

  // event-wide
  virtual double GetEhad() const;
  virtual double GetEnu() const;
  virtual double GetQ2() const;
  virtual double GetWexp() const;
  virtual double GetTracklessWexp() const;
  virtual double GetTrackedWexp() const;
  virtual double Getq0() const;
  virtual double GetTracklessq3() const;

  // pion
  virtual double GetALR(RecoPionIdx) const;
  virtual double GetAdlerCosTheta(RecoPionIdx) const;
  virtual double GetAdlerPhi(RecoPionIdx) const;
  virtual double GetEpi(RecoPionIdx) const;
  virtual double GetPT(RecoPionIdx) const;
  virtual double GetPXpi(RecoPionIdx) const;
  virtual double GetPYpi(RecoPionIdx) const;
  virtual double GetPZpi(RecoPionIdx) const;
  virtual double GetPpi(RecoPionIdx) const;
  virtual double GetThetapi(RecoPionIdx) const;
  virtual double GetThetapiDeg(RecoPionIdx) const;
  virtual double GetTpi(RecoPionIdx) const;
  virtual double GetTpiMBR(RecoPionIdx) const;
  virtual double GetpimuAngle(RecoPionIdx) const;
  virtual double Gett(RecoPionIdx) const;
  virtual double GetTpiTrackless() const;
  virtual double GetBestDistance() const;
  virtual double GetMixedTpi(RecoPionIdx) const;
  virtual double GetPpionCorr(RecoPionIdx hadron) const;

  //==============================================================================
  // Truth
  //==============================================================================
  virtual double GetALRTrue(TruePionIdx) const;
  virtual double GetAdlerCosThetaTrue(TruePionIdx) const;
  virtual double GetAdlerPhiTrue(TruePionIdx) const;
  virtual double GetAllTrackEnergyTrue() const;
  virtual double GetEmuTrue() const;
  virtual double GetIntVtxXTrue() const;
  virtual double GetIntVtxYTrue() const;
  virtual double GetIntVtxZTrue() const;
  virtual double GetPTTrue(TruePionIdx) const;
  virtual double GetPTmuTrue() const;
  virtual double GetPXmuTrue() const;
  virtual double GetPYmuTrue() const;
  virtual double GetPZmuTrue() const;
  virtual double GetPmuTrue() const;
  virtual double GetThetamuTrue() const;
  virtual double GetThetamuTrueDeg() const;
  virtual double GetThetapiTrue(TruePionIdx) const;
  virtual double GetThetapiTrueDeg(TruePionIdx) const;
  virtual double GetTpiTrue(TruePionIdx) const;
  virtual double GetWexpTrue() const;
  virtual double GetWgenie() const;
  virtual double GetpimuAngleTrue(TruePionIdx) const;
  virtual int GetNChargedPionsTrue() const;
  virtual int GetPiChargeTrue(TruePionIdx) const;
  virtual std::vector<double> GetTpiTrueVec() const;
  virtual double GetMixedTpiTrue(TruePionIdx) const;

  //==============================
  // Ehad (GetErecoil) Variables
  //==============================
  // ehad and related variables
  virtual double GetCalEpi(RecoPionIdx) const;
  virtual double GetCalRecoilEnergy() const;
  virtual double GetCalRecoilEnergyNoPi_Corrected(const double ecal_nopi) const;
  virtual double GetCalRecoilEnergyNoPi_DefaultSpline() const;
  virtual double GetCalRecoilEnergy_CCPiSpline() const;
  virtual double GetCalRecoilEnergy_DefaultSpline() const;
  virtual double GetNonCalRecoilEnergy() const;
  virtual double GetTrackRecoilEnergy() const;

  // ehad old variables
  virtual double GetCalRecoilEnergyNoPi_CCIncSpline() const;
  virtual double GetCalRecoilEnergy_CCIncSpline() const;

  // ehad truth variables
  virtual double GetCalRecoilEnergyNoPiTrue() const;
  virtual double GetEhadTrue() const;
  virtual double GetEpiTrueMatched(RecoPionIdx) const;
  virtual double GetTpiTrueMatched(RecoPionIdx) const;

  //==============================================================================
  // Cuts, Systematics, Studies
  //==============================================================================
  virtual bool IsInHexagon(double x, double y, double apothem) const;
  virtual bool IsInPlastic() const;
  virtual double GetEmichel(RecoPionIdx) const;
  virtual double GetEnode0(RecoPionIdx) const;
  virtual double GetEnode01(RecoPionIdx) const;
  virtual double GetEnode1(RecoPionIdx) const;
  virtual double GetEnode2(RecoPionIdx) const;
  virtual double GetEnode3(RecoPionIdx) const;
  virtual double GetEnode4(RecoPionIdx) const;
  virtual double GetEnode5(RecoPionIdx) const;
  virtual double GetFitVtxX() const;
  virtual double GetFitVtxY() const;
  virtual double GetFitVtxZ() const;
  virtual double GetLLRScore(RecoPionIdx) const;
  virtual double GetLargestIsoProngSep() const;
  virtual double GetLargestPrimProngSep() const;
  virtual double GetTpiFResidual(const int hadron) const;
  //  virtual double GetTpiFResidual(const int hadron,
  //                                 const bool MBR = false) const;
  virtual double GetWexpFResidual() const;
  virtual double GetdEdxScore(RecoPionIdx) const;
  virtual int GetNAnchoredLongTracks() const;
  virtual int GetNAnchoredShortTracks() const;
  virtual int GetNIsoProngs() const;
  virtual int GetNNodes(RecoPionIdx) const;
  virtual int GetNhadrons() const;
  virtual int GetTrackReconstructionMethod(RecoPionIdx) const;

  //==============================================================================
  // Weights
  //==============================================================================
  virtual double GetAnisoDeltaDecayWarpWeight() const;
  virtual double GetDiffractiveWeight() const;
  virtual double GetGenieWarpWeight(double p) const;
  virtual double GetLowQ2PiWeight(double q2, std::string channel) const;
  virtual double GetWeight() const;

  //==============================================================================
  // Physics Calculations
  //==============================================================================
  TVector3 AdlerAngle(int RefSystemDef, double dmumom, double dpimom,
                      TVector3 NeuDir, TVector3 MuDir, TVector3 PiDir,
                      double Enu) const;
  double CalcQ2(const double Enu, const double Emu, const double Thetamu) const;
  double CalcWexp(const double Q2, const double Ehad) const;
  double Calcq0(const double Enu, const double Emu) const;
  double Calcq3(const double Q2, const double Enu, const double Emu) const;
  double Calct(const double epi, const double emu, const double pzpi,
               const double pzmu, const double pxpi, const double pxmu,
               const double pypi, const double pymu) const;

  //==============================================================================
  // Untracked pions functions
  //==============================================================================
  // Mehreen's full fit, circa 2022-02:
  // Tpi = p + q*range + r*sqrt(range) with
  // p = -2.93015 +- 4.44962    // Yikes, BTW
  // r = 0.132851 +- 0.0199247
  // q = 3.95884  +- 0.657313
  // Mehreen's fit, Minerva week
  // KE = q*range + r*sqrt(range)
  // q = 0.210207 +- 2.38011e-3
  // r = 2.90140 +- 6.06231
  virtual double GetTpiUntracked(double michel_range) const {
    return GetTpiFromRange(michel_range);
    //    return -2.93 + 0.133 * michel_range + 3.96 * sqrt(michel_range);
    //    return 0.210207 * michel_range + 2.9014 * sqrt(michel_range);
  }

  virtual double GetEavail() const;
  virtual double GetThetapitrackless() const;
  virtual double GetThetapitracklessDeg() const;
  virtual double GetMixedThetapiDeg(RecoPionIdx) const;
  virtual double GetThetapitracklessTrue() const;
  virtual double GetThetapitracklessTrueDeg() const;
  virtual double GetMixedThetapiTrueDeg(TruePionIdx) const;
  virtual double thetaWRTBeam(double x, double y, double z) const {
    double pyp = -1.0 * sin(MinervaUnits::numi_beam_angle_rad) * z +
                 cos(MinervaUnits::numi_beam_angle_rad) * y;
    double pzp = cos(MinervaUnits::numi_beam_angle_rad) * z +
                 sin(MinervaUnits::numi_beam_angle_rad) * y;
    double denom2 = pow(x, 2) + pow(pyp, 2) + pow(pzp, 2);
    if (0. == denom2)
      return -9999.;
    else
      return acos(pzp / sqrt(denom2));
  }
  virtual int GetNTruePions() const {
    return GetInt("truth_FittedMichel_all_piontrajectory_trackID_sz");
  }
  virtual double GetTrueTpi() const {
    std::vector<double> tpi_vec = GetTpiTrueVec();  // pip and pim
    const int n_true_pions = GetNChargedPionsTrue();

    TruePionIdx reigning_idx = -1;
    double reigning_tpi = 9999;
    for (TruePionIdx idx = 0; idx < n_true_pions; ++idx) {
      if (tpi_vec[idx] < reigning_tpi && GetPiChargeTrue(idx) > 0.) {
        reigning_idx = idx;
        reigning_tpi = tpi_vec[idx];
      }
    }

    return GetTpiTrue(reigning_idx);
  }
};

#endif  // CVUniverse_H

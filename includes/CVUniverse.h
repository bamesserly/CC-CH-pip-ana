#ifndef CVUniverse_H
#define CVUniverse_H

#include "Binning.h"    // CCPi::GetBinning for ehad_nopi
#include "Constants.h"  // CCNuPionIncConsts, CCNuPionIncShifts, Reco/TruePionIdx
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/MinervaUniverse.h"

class CVUniverse : public PlotUtils::MinervaUniverse {
 private:
  // Pion Candidates - clear these when SetEntry is called
  std::vector<RecoPionIdx> m_pion_candidates;
 public:
#include "PlotUtils/MuonFunctions.h"
#include "PlotUtils/TruthFunctions.h"
#include "PlotUtils/WeightFunctions.h"
#include "PlotUtils/RecoilEnergyFunctions.h"
#include "PlotUtils/MichelFunctions.h"
  // CTOR
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0);

  // DTOR
  virtual ~CVUniverse(){};

  // No stale cache!
  virtual void OnNewEntry() override { m_pion_candidates.clear(); }

  // Get Branches and Calculate Quantities for the universe/event
  // Muon Variables
  virtual double GetPTmu() const;
  virtual double GetThetamuDeg() const;
  virtual double GetPXmu() const;
  virtual double GetPYmu() const;
  virtual double GetPZmu() const;

  virtual double GetEmuTrue() const;
  virtual double GetPmuTrue() const;
  virtual double GetPTmuTrue() const;
  virtual double GetPZmuTrue() const;
  virtual double GetThetamuTrue() const;
  virtual double GetThetamuTrueDeg() const;

  // Event-wide Variables
  virtual double GetEnu() const;
  virtual double GetQ2() const;
  virtual double Getq0() const;
  virtual double Getq3() const;
  virtual double GetWexp() const;

  virtual double GetWexpTrue() const;
  virtual double GetWgenie() const;

  // Hadron (track) Variables
  virtual double GetTpi(RecoPionIdx) const;
  virtual double GetTpiMBR(RecoPionIdx) const;
  virtual double GetThetapi(RecoPionIdx) const;
  virtual double GetThetapiDeg(RecoPionIdx) const;
  virtual double GetLLRScore(RecoPionIdx) const;
  virtual double GetdEdxScore(RecoPionIdx) const;
  virtual double GetEmichel(RecoPionIdx) const;

  virtual double GetEnode0(RecoPionIdx) const;
  virtual double GetEnode1(RecoPionIdx) const;
  virtual double GetEnode01(RecoPionIdx) const;
  virtual double GetEnode2(RecoPionIdx) const;
  virtual double GetEnode3(RecoPionIdx) const;
  virtual double GetEnode4(RecoPionIdx) const;
  virtual double GetEnode5(RecoPionIdx) const;

  virtual double GetPZpi(RecoPionIdx) const;
  virtual double GetPXpi(RecoPionIdx) const;
  virtual double GetPYpi(RecoPionIdx) const;
  virtual double GetPpi(RecoPionIdx) const;

  virtual double Gett(RecoPionIdx) const;
   
  virtual int GetNhadrons() const;

  // With these truth hadron variables, SEE the warning in the .cxx
  virtual double GetTpiTrue(TruePionIdx) const;
  virtual std::vector<double> GetTpiTrueVec() const;
  virtual double GetThetapiTrue(TruePionIdx) const;
  virtual double GetThetapiTrueDeg(TruePionIdx) const;
  virtual int GetPiChargeTrue(TruePionIdx) const;
  virtual int GetNChargedPionsTrue() const;
  virtual double GetAllTrackEnergyTrue() const;

  // Misc
  virtual double GetLargestPrimProngSep() const;
  virtual double GetLargestIsoProngSep() const;
  virtual int GetNIsoProngs() const;
  double GetTpiFResidual(const int hadron, const bool MBR = false) const;
  double GetWexpFResidual() const;

  // Some study variables
  virtual int GetNAnchoredShortTracks() const;
  virtual int GetNAnchoredLongTracks() const;
  virtual double GetFitVtxX() const;
  virtual double GetFitVtxY() const;
  virtual double GetFitVtxZ() const;
  virtual int GetTrackReconstructionMethod(RecoPionIdx) const;
  virtual int GetNNodes(RecoPionIdx) const;

  // Reco pion energy from dEdx tool.
  virtual double GetEpi(RecoPionIdx) const;

  // Ehad new variables
  virtual double GetEhad() const;
  virtual double GetCalRecoilEnergy() const;
  virtual double GetTrackRecoilEnergy() const;
  virtual double GetNonCalRecoilEnergy() const;
  virtual double GetCalRecoilEnergyNoPi_DefaultSpline() const;
  virtual double GetCalRecoilEnergyNoPi_Corrected(const double ecal_nopi) const;
  virtual double GetCalRecoilEnergy_DefaultSpline() const;
  virtual double GetCalRecoilEnergy_CCPiSpline() const;
  virtual double GetCalEpi(RecoPionIdx) const;

  // Ehad old variables
  virtual double GetCalRecoilEnergy_CCIncSpline() const;
  virtual double GetCalRecoilEnergyNoPi_CCIncSpline() const;

  // Ehad truth variables
  virtual double GetEhadTrue() const;
  virtual double GetTpiTrueMatched(RecoPionIdx) const;
  virtual double GetEpiTrueMatched(RecoPionIdx) const;
  virtual double GetCalRecoilEnergyNoPiTrue() const;

  // Dummy access for variable constructors
  virtual double GetDummyVar() const;
  virtual double GetDummyHadVar(const int x) const;

  // Calculate Quantities(Always MeV)
  double CalcQ2(const double Enu, const double Emu, const double Thetamu) const;
  double CalcWexp(const double Q2, const double Ehad) const;
  double Calcq0(const double Enu, const double Emu) const;
  double Calcq3(const double Q2, const double Enu, const double Emu) const;
  double Calct(const double epi, const double emu, const double pzpi,
               const double pzmu, const double pxpi, const double pxmu,
               const double pypi, const double pymu) const;
  // Functions to make fidvol cut
  virtual bool rightlinesCut (const double a,const double x,const double y) const;
  virtual bool leftlinesCut (const double a,const double x,const double y) const;
  // Get Weight
  virtual double GetWeight() const;

  // Warping
  virtual double GetGenieWarpWeight() const;
  virtual double GetLowQ2PiWarpWeight(double q2, std::string channel) const;
  virtual double GetAnisoDeltaDecayWarpWeight() const;

  // Get Highest Energy Pion from vector of candidates
  int GetHighestEnergyPionCandidateIndex(
      const std::vector<int>& pion_candidate_idxs) const;
  TruePionIdx GetHighestEnergyTruePionIndex() const;

  // Print arachne link
  void PrintArachneLink() const;

  // Pion Candidates
  void SetPionCandidates(std::vector<RecoPionIdx> c);
  std::vector<RecoPionIdx> GetPionCandidates() const;
};

#endif  // CVUniverse_H

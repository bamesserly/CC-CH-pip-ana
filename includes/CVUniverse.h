#ifndef CVUniverse_H
#define CVUniverse_H

#include <TVector3.h>

#include "Binning.h"    // CCPi::GetBinning for ehad_nopi
#include "Constants.h"  // CCNuPionIncConsts, CCNuPionIncShifts, Reco/TruePionIdx
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/MinervaUniverse.h"

class CVUniverse : public PlotUtils::MinervaUniverse {
 private:
  // Pion Candidates - clear these when SetEntry is called
  std::vector<RecoPionIdx> m_pion_candidates;

 public:
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

  // No stale cache!
  virtual void OnNewEntry() override { m_pion_candidates.clear(); }

  // Get and set pion candidates
  TruePionIdx GetHighestEnergyTruePionIndex() const;
  int GetHighestEnergyPionCandidateIndex(const std::vector<int>& pions) const;
  std::vector<RecoPionIdx> GetPionCandidates() const;
  void SetPionCandidates(std::vector<RecoPionIdx> c);

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
  virtual double Getq0() const;
  virtual double Getq3() const;

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

  //==============================================================================
  // Adler
  //==============================================================================
  TVector3 TejinRefSys(TVector3 Pnu, TVector3 Pmu, TVector3 var) const;
  virtual TVector3 GetPnuVecWRTB() const;
  virtual TVector3 GetPnuVecWRTBTrue() const;
  virtual TVector3 GetPpiVecWRTB(RecoPionIdx) const;
  virtual TVector3 GetPpiVecWRTBTrue(TruePionIdx) const;
  virtual double GetPXnu() const;
  virtual double GetPXnuTrue() const;
  virtual double GetPXpiTrue(TruePionIdx) const;
  virtual double GetPYnu() const;
  virtual double GetPYnuTrue() const;
  virtual double GetPYpiTrue(TruePionIdx) const;
  virtual double GetPZnu() const;
  virtual double GetPZnuTrue() const;
  virtual double GetPZpiTrue(TruePionIdx) const;
  virtual double GetPpiTrue(TruePionIdx) const;
  virtual double GetthetaZ() const;
  virtual double GetthetaZTrue() const;

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
  virtual bool leftlinesCut(const double a, const double x,
                            const double y) const;
  virtual bool rightlinesCut(const double a, const double x,
                             const double y) const;
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
  virtual double GetTpiFResidual(const int hadron,
                                 const bool MBR = false) const;
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
  virtual double GetGenieWarpWeight() const;
  virtual double GetLowQ2PiWarpWeight(double q2, std::string channel) const;
  virtual double GetWeight() const;

  //=============================================================================3
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
};

#endif  // CVUniverse_H

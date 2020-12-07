#ifndef CVUniverse_H
#define CVUniverse_H

#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/DefaultCVUniverse.h"

#include "Constants.h" // CCNuPionIncConsts, CCNuPionIncShifts, Reco/TruePionIdx

class CVUniverse : public PlotUtils::DefaultCVUniverse {
  public:
    // CTOR
    CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma=0);

    // DTOR
    virtual ~CVUniverse(){};

    // Get Branches and Calculate Quantities for the universe/event
    // Muon Variables
    virtual double GetEmu() const;
    virtual double GetPmu() const;
    virtual double GetPTmu() const;
    virtual double GetPZmu() const;
    virtual double GetThetamu() const;
    virtual double GetThetamuDeg() const;

    virtual double GetEmuTrue() const;
    virtual double GetPmuTrue() const; 
    virtual double GetPTmuTrue() const; 
    virtual double GetPZmuTrue() const; 
    virtual double GetThetamuTrue() const;
    virtual double GetThetamuTrueDeg() const;

    // Event-wide Variables
    virtual double GetEnu() const; 
    virtual double GetEhad() const; // just = GetCalRecoilEnergy
    virtual double GetCalRecoilEnergy() const;
    virtual double GetQ2() const;
    virtual double Getq0() const;
    virtual double Getq3() const;
    virtual double GetWexp() const;

    virtual double GetEnuTrue() const ;
    virtual double GetEhadTrue() const ;
    virtual double GetQ2True() const ;
    virtual double Getq0True() const ;
    virtual double Getq3True() const ;
    virtual double GetWexpTrue() const ;
    virtual double GetWgenie() const ;

    // Hadron (track) Variables
    virtual double GetTpi(RecoPionIdx hadron) const;
    virtual double GetTpiMBR(RecoPionIdx hadron) const;
    virtual double GetThetapi(RecoPionIdx hadron) const;
    virtual double GetThetapiDeg(RecoPionIdx hadron) const;
    virtual double GetLLRScore(RecoPionIdx hadron ) const;
    virtual double GetdEdxScore(RecoPionIdx hadron) const;
    virtual double GetEmichel(RecoPionIdx hadron) const;

    virtual double GetEnode0(RecoPionIdx hadron) const;
    virtual double GetEnode1(RecoPionIdx hadron) const;
    virtual double GetEnode01(RecoPionIdx hadron) const;
    virtual double GetEnode2(RecoPionIdx hadron) const;
    virtual double GetEnode3(RecoPionIdx hadron) const;
    virtual double GetEnode4(RecoPionIdx hadron) const;
    virtual double GetEnode5(RecoPionIdx hadron) const;

    // With these truth hadron variables, SEE the warning in the .cxx
    virtual double GetTpiTrue(TruePionIdx idx) const;
    virtual std::vector<double> GetTpiTrueVec() const;
    virtual double GetThetapiTrue(TruePionIdx idx) const;
    virtual double GetThetapiTrueDeg(TruePionIdx idx) const;
    virtual int GetPiChargeTrue(TruePionIdx idx) const;
    virtual int GetNChargedPionsTrue() const;

    // Misc
    virtual double GetLargestPrimProngSep() const;
    virtual double GetLargestIsoProngSep() const;
    virtual int GetNIsoProngs() const;
    double GetTpiFResidual(const int hadron, const bool MBR = false) const;

    // Some study variables
    virtual int GetNAnchoredShortTracks() const;
    virtual int GetNAnchoredLongTracks() const;
    virtual double GetFitVtxX() const;
    virtual double GetFitVtxY() const;
    virtual double GetFitVtxZ() const;
    virtual int GetTrackReconstructionMethod(RecoPionIdx hadron) const;
    virtual int GetNNodes(RecoPionIdx hadron) const;

    // Dummy access for variable constructors
    virtual double GetDummyVar() const;
    virtual double GetDummyHadVar(const int x) const;

    // Calculate Quantities(Always MeV)
    double CalcQ2(const double Enu, const double Emu, const double Thetamu) const;
    double CalcWexp(const double Q2, const double Ehad) const;
    double Calcq0(const double Enu, const double Emu) const;
    double Calcq3(const double Q2, const double Enu, const double Emu) const;

    // Get Weight
    virtual double GetWeight() const;

    // Warping
    virtual double GetGenieWarpWeight() const;
    virtual double GetLowQ2PiWarpWeight(double q2, std::string channel) const;
    virtual double GetAnisoDeltaDecayWarpWeight() const;

    // Get Highest Energy Pion from vector of candidates
    int GetHighestEnergyPionCandidateIndex(const std::vector<int>& pion_candidate_idxs) const;
    TruePionIdx GetHighestEnergyTruePionIndex() const;

    // Print arachne link
    void PrintArachneLink() const;
};

#endif // CVUniverse_H

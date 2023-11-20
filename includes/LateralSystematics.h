#ifndef LateralSystematics_H
#define LateralSystematics_H

#include <iostream>

#include "CVUniverse.h"
#include "Constants.h" // CCNuPionIncShifts 
#include "common_functions.h" // FixAngle
#include "TRandom.h" // gRandom


//TODO shift for michel electron energy
class BirksShiftUniverse : public CVUniverse {
  public:
    BirksShiftUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
      : CVUniverse(chw, nsigma)
    {}


    virtual double GetTpi(int hadron) const {
      double shift_val = GetVecElem("MasterAnaDev_pion_E_Birks", hadron);
      return shift_val+CVUniverse::GetTpi(hadron);
    }


    virtual double GetdEdxScore(int hadron) const {
      double shift_val = GetVecElem("MasterAnaDev_piFit_score1_Birks",
                                    hadron);
      return shift_val+CVUniverse::GetdEdxScore(hadron);
    }


    virtual std::string ShortName() const { return "Birks"; }
    virtual std::string LatexName() const { return "Birks"; }
};


//TODO get muon E shift from muonutils!!
class BetheBlochShiftCVUniverse : public CVUniverse {
  public:
    BetheBlochShiftCVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
      : CVUniverse(chw, nsigma)
    {}


    virtual double GetEmu() const { 
      double shift_val = m_nsigma*30.0; //30 MeV from MuonUtils
      return shift_val+CVUniverse::GetEmu();
    }


    virtual double GetTpi(int hadron) const {
      double shift_val = m_nsigma > 0 ? 
          GetVecElem("MasterAnaDev_pion_E_BetheBloch_biasUp",   hadron) :
          GetVecElem("MasterAnaDev_pion_E_BetheBloch_biasDown", hadron);
      return shift_val+CVUniverse::GetTpi(hadron);
    }


    virtual double GetdEdxScore(int hadron) const {
      double shift_val = m_nsigma > 0 ? 
          GetVecElem("MasterAnaDev_piFit_score1_BetheBloch_biasUp",   hadron) :
          GetVecElem("MasterAnaDev_piFit_score1_BetheBloch_biasDown", hadron);
      return shift_val+CVUniverse::GetdEdxScore(hadron);
    }


    virtual std::string ShortName() const { return "BetheBloch"; }
    virtual std::string LatexName() const { return "Bethe Bloch"; }
};


//TODO get muon E shift from muonutils!!
//TODO vertical mass uncertainty?
class DetectorMassShiftCVUniverse : public CVUniverse {
  public:
    DetectorMassShiftCVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
      : CVUniverse(chw, nsigma)
    {}


    virtual double GetTpi(int hadron) const {
      double shift_val = m_nsigma > 0 ?
          GetVecElem("MasterAnaDev_pion_E_Mass_biasUp",   hadron) :
          GetVecElem("MasterAnaDev_pion_E_Mass_biasDown", hadron);
      return shift_val+CVUniverse::GetTpi(hadron);
    }


    virtual double GetdEdxScore(int hadron) const {
      double shift_val = m_nsigma > 0 ?
          GetVecElem("MasterAnaDev_piFit_score1_Mass_biasUp",   hadron) :
          GetVecElem("MasterAnaDev_piFit_score1_Mass_biasDown", hadron);
      return shift_val+CVUniverse::GetdEdxScore(hadron);
    }
    virtual std::string ShortName() const { return "Mass"; }
    virtual std::string LatexName() const { return "Detector Mass"; }
};


class TrackAngleShiftCVUniverse : public CVUniverse {
  public:
    TrackAngleShiftCVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
      : CVUniverse(chw, nsigma)
    {}


    virtual double GetThetamu() const { 
      double shift_val = (CCNuPionIncShifts::muon_angle_res)*
                             (gRandom->Gaus(0.0,1.0));
      return FixAngle(shift_val + CVUniverse::GetThetamu());
    }


    virtual double GetThetapi(int hadron) const { 
      double shift_val = (CCNuPionIncShifts::pion_angle_res)*
                              (gRandom->Gaus(0.0,1.0));
      return FixAngle(shift_val + CVUniverse::GetThetapi(hadron));
    }


    virtual std::string ShortName() const { return "TrackAngle"; }
    virtual std::string LatexName() const { return "Track Angle"; }
};


class BeamAngleShiftCVUniverse : public CVUniverse {
  public:
    BeamAngleShiftCVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
      : CVUniverse(chw, nsigma)
    {}

    virtual double GetThetamu() const { 
      double shift_val = m_nsigma > 0 ? 
          GetDouble("MasterAnaDev_muon_theta_biasUp") :
          GetDouble("MasterAnaDev_muon_theta_biasDown") ;
      return FixAngle(shift_val + CVUniverse::GetThetamu());
    }
    virtual double GetThetapi(int hadron) const { 
      double shift_val = m_nsigma > 0 ?
          GetVecElem("MasterAnaDev_pion_theta_biasUp",   hadron):
          GetVecElem("MasterAnaDev_pion_theta_biasDown", hadron);
      return FixAngle(shift_val + CVUniverse::GetThetapi(hadron));
    }

    virtual std::string ShortName() const { return "BeamAngle"; }
    virtual std::string LatexName() const { return "Beam Angle"; }
};

class NodeCutEffUniverse : public CVUniverse {
 public:
  NodeCutEffUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
      : CVUniverse(chw, nsigma) {}

  virtual double GetWeight() const override {
    if (IsTruth())
      return CVUniverse::GetWeight();

    // I don't really want to put in the prong for such a small reweight
    // Leading pion
    double eff_weight = 1.;
    if (CVUniverse::GetPionCandidates().size() > 0) {
      double Tpi =
          CVUniverse::GetTpi(CVUniverse::GetHighestEnergyPionCandidateIndex(CVUniverse::GetPionCandidates())) /
          1000;

      if (Tpi < 0.15) eff_weight *= (1 + m_nsigma * 0.0015);
      if (Tpi > 0.15) eff_weight *= (1 + m_nsigma * 0.01);
    }
    return eff_weight * CVUniverse::GetWeight();
  }
  virtual std::string ShortName() const override { return "NodeCutEff"; }
  virtual std::string LatexName() const override {
    return "Node Cut Efficiency";
  }
  virtual bool IsVerticalOnly() const override { return true; }
};

#endif // LateralSystematics_H

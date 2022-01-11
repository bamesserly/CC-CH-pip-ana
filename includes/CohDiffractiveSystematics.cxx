#ifndef CohDiffractiveSystematics_CXX
#define CohDiffractiveSystematics_CXX 1

#include "CohDiffractiveSystematics.h"

#include <iostream>

#include "PlotUtils/TargetUtils.h"

std::vector<CVUniverse*> GetCohDiffractiveSystematics(
    PlotUtils::ChainWrapper* chain) {
  std::vector<CVUniverse*> ret;
  ret.push_back(new DiffractiveUniverse(chain, -1., fracDiffractiveUnc));
  ret.push_back(new DiffractiveUniverse(chain, 1., fracDiffractiveUnc));

  ret.push_back(
      new CoherentPiPlasticUniverse(chain, -1., fracCoherentPiUncTracker));
  ret.push_back(
      new CoherentPiPlasticUniverse(chain, 1., fracCoherentPiUncTracker));

  return ret;
};

std::map<std::string, std::vector<CVUniverse*> >
GetCohDiffractiveSystematicsMap(PlotUtils::ChainWrapper* chain) {
  std::map<std::string, std::vector<CVUniverse*> > ret;

  ret["DiffractiveModelUnc"].push_back(
      new DiffractiveUniverse(chain, -1., fracDiffractiveUnc));
  ret["DiffractiveModelUnc"].push_back(
      new DiffractiveUniverse(chain, 1., fracDiffractiveUnc));

  ret["CoherentPiUncTracker_CH"].push_back(
      new CoherentPiPlasticUniverse(chain, -1., fracCoherentPiUncTracker));
  ret["CoherentPiUncTracker_CH"].push_back(
      new CoherentPiPlasticUniverse(chain, 1., fracCoherentPiUncTracker));

  return ret;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Class Definitions
// Constructor
DiffractiveUniverse::DiffractiveUniverse(PlotUtils::ChainWrapper* chw,
                                         double nsigma, double fracDiffUnc)
    : CVUniverse(chw, nsigma), m_fracDiffUnc(fracDiffUnc) {}

double DiffractiveUniverse::GetDiffractiveWeight() const {
  if (GetInt("mc_intType") != 4) return 1.;
  if (!IsInPlastic()) return 1.;

  double cv_wgt = CVUniverse::GetDiffractiveWeight();
  double shift_val = m_nsigma * m_fracDiffUnc * (cv_wgt - 1);
  return cv_wgt + shift_val;
}

std::string DiffractiveUniverse::ShortName() const {
  return "DiffractiveModelUnc";
}

std::string DiffractiveUniverse::LatexName() const {
  return "Diffractive Pion Uncertainty";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class Definitions
// Constructor
CoherentPiPlasticUniverse::CoherentPiPlasticUniverse(
    PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc)
    : CVUniverse(chw, nsigma), m_fracCohPiUnc(fracCohPiUnc) {}

double CoherentPiPlasticUniverse::GetCoherentPiWeight(double thpi_true,
                                                      double tpi_true) const {
  if (GetInt("mc_intType") != 4) return 1.;

  // In Plastic?
  if (PlotUtils::TargetUtils::Get().InPassiveTargetMC(
          GetIntVtxXTrue(), GetIntVtxYTrue(), GetIntVtxZTrue(),
          GetInt("mc_targetZ")))
    return 1.;

  double shift_val = 1 + m_nsigma * m_fracCohPiUnc;
  return shift_val * CVUniverse::GetCoherentPiWeight(thpi_true, tpi_true);
}

std::string CoherentPiPlasticUniverse::ShortName() const {
  return "CoherentPiUnc_CH";
}

std::string CoherentPiPlasticUniverse::LatexName() const {
  return "Coherent Pion Scintillator Uncertainty";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif  // CohDiffractiveSystematics_CXX

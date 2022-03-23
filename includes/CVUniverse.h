#ifndef CVUniverse_H
#define CVUniverse_H

#include "Binning.h"    // CCPi::GetBinning for ehad_nopi
#include "Constants.h"  // CCNuPionIncConsts, CCNuPionIncShifts, Reco/TruePionIdx
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/MinervaUniverse.h"
#include <TVector3.h>

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
#include "UniverseFunctions/WeightFunctions.h"
#include "UniverseFunctions/TruthFunctions.h"
#include "UniverseFunctions/EhadFunctions.h"
#include "UniverseFunctions/PhysicsCalculations.h"
#include "UniverseFunctions/UtilityFunctions.h"
#include "UniverseFunctions/AnalysisFunctions.h"
  // CTOR
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0);

  // DTOR
  virtual ~CVUniverse(){};

  // No stale cache!
  virtual void OnNewEntry() override { m_pion_candidates.clear(); }

  // Dummy access for variable constructors
  virtual double GetDummyVar() const;
  virtual double GetDummyHadVar(const int x) const;

  // Print arachne link
  void PrintArachneLink() const;

  // Get/set all/best pion candidates
  void SetPionCandidates(std::vector<RecoPionIdx> c);
  std::vector<RecoPionIdx> GetPionCandidates() const;
  int GetHighestEnergyPionCandidateIndex(const std::vector<int>& pions) const;
  TruePionIdx GetHighestEnergyTruePionIndex() const;
};

#endif  // CVUniverse_H

#ifndef PionReco_H
#define PionReco_H

enum class PionRecoType
{
  kT, kUT, kNPionRecoTypes
};

//==============================================================================
// Simple Signal-Background (no BG breakdown)
//==============================================================================
//// TODO -- set this in event when we know what we're doing.
//PionRecoType GetPionRecoType(const CCPiEvent& event) {
//  return event.m_pionreco_type;
//}

std::string GetTruthClassification_LegendLabel(PionRecoType category) {
  switch (category) {
    case PionRecoType::kT:
      return "Tracked";
    case PionRecoType::kUT:
      return "Untracked";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(PionRecoType category) {
  switch (category) {
    case PionRecoType::kT:
      return "Tracked";
    case PionRecoType::kUT:
      return "Untracked";
    default:
      return "ERROR";
  }
}


#endif

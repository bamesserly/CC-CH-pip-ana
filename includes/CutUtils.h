#ifndef CutUtils_h
#define CutUtils_h

#include "Constants.h"  // enum ECuts, CCNuPionIncConsts
#include "Michel.h" // endpoint::MichelMap

// Analysis Cuts - default vector
const std::vector<ECuts> kCutsVector = {
    kNoCuts,
    kPrecuts,
    kVtx,
    kMinosMuon,
    kPmu,
    kAtLeastOneMichel,
    kAtLeastOnePionCandidate,
    kLLR,
    kNode,
    kTrackQuality,
    kWexp,
    kIsoProngs,
    kPionMult
};

// Remove W cut from cuts vector
const std::vector<ECuts> GetWSidebandCuts() {
  std::vector<ECuts> w_sideband_cuts = kCutsVector;
  w_sideband_cuts.erase(
      std::remove(w_sideband_cuts.begin(), w_sideband_cuts.end(), kWexp),
      w_sideband_cuts.end());
  return w_sideband_cuts;
}

// Gaudi tool cuts - only work when checking truth tuple
bool IsPrecut(ECuts c) {
  if (c == kNoCuts || c == kGoodObjects || c == kGoodVertex ||
      c == kFiducialVolume || c == kMinosActivity || c == kPrecuts)
    return true;
  else
    return false;
}

// Cut Names
std::string GetCutName(ECuts cut) {
  switch (cut) {
    case kNoCuts:
      return "No Cuts";

    case kGoodObjects:
      return "Good Objects";

    case kGoodVertex:
      return "Good Vertex";

    case kFiducialVolume:
      return "Fiducial Volume";

    case kMinosActivity:
      return "MINOS Activity";

    case kPrecuts:
      return "Anatool Precuts";

    case kVtx:
      return "vertex position Cut";

    case kMinosMatch:
      return "MINOS Muon";

    case kMinosCharge:
      return "MINOS Charge";

    case kMinosCoil:
      return "MINOS Coil";

    case kMinosMuon:
      return "MINOS Muon";

    case kThetaMu:
      return "Muon Angle";

    case kDeadTime:
      return "Dead Time";

    case kWexp:
      return "$W_{experimental}$";

    case kIsoBlobs:
      return "$<$1 Isolated Blobs";

    case kIsoProngs:
      return "$<$2 Isolated Prongs";

    case kIsoProngSep:
      return "Iso Prong Sep $<$ 300";

    case kNProngs:
      return "Max 2 Hadr Prongs";

    case kNPionCandidates:
      return "$\\pi$ candidate";

    case kPionCandidateQuality:
      return "Quality $\\pi$ candidate";

    case kAtLeastOneMichel:
      return "$>$= 1 Michel";

    case kAtLeastOneBrandonMichel:
      return "$>$= 1 Brandon Michel";

    case kAtLeastOneAnchoredProng:
      return "$>$= 1 Anchored Prong";

    case kAtLeastOnePionCandidateTrack:
      return "$>$= 1 Hadron Track";

    case kAtLeastOneLLRCandidate:
      return "$>$= 1 Track Pass LLR Cut";

    case kAtLeastOneNodeCandidate:
      return "$>$= 1 Track Pass Node Cut";

    case kExactlyOneEndpointMichel:
      return "== 1 Michel";

    case kNode:
      return "Node";

    case kPionMult:
      return "Pion Multiplicity";

    case kOldMichel:
      return "Old Michel Cut";

    case kdEdx:
      return "dEdx PID";

    case kLLR:
      return "LLR PID";

    case kAllCuts:
      return "Total";

    case kTrackQuality:
      return "General Track Quality";

    case kPmu:
      return "1.5 GeV $<$ Pmu $<$ 20 GeV";

    case kAtLeastOnePionCandidate:
      return "At Least One Pion";

    default:
      std::cout << "ERROR: GetCutName unknown cut!" << std::endl;
      return "";
  };
}

// Get pion candidate indexes from michel map
// (our cuts strategy enforces a 1-1 michel-pion candidate match)
std::vector<int> GetHadIdxsFromMichels(const endpoint::MichelMap endpoint_michels, const trackless::MichelEvent vtx_michels) {
  std::vector<int> ret;

  // endpoint michels
  for (auto m : endpoint_michels) ret.push_back(m.second.had_idx);

  // vertex michels
  // When m_idx is set (i.e. != -1), then we have a good vertex michel.
  // In that case, a -1 in this hadron index return vector is the code that for
  // this analysis that we have a good vertex michel.
  if (vtx_michels.m_idx != -1) ret.push_back(-1); 

  return ret;
}

#endif

#ifndef CutUtils_h
#define CutUtils_h

#include "Constants.h"    // enum ECuts, CCNuPionIncConsts

// Analysis Cuts - default vector
const std::vector<ECuts> kCutsVector = {kNoCuts,
                                        kPrecuts,
  					khasIntVtx, // AaronCuts
  					kMultiplicityCut, // AaronCuts
  					kExitingMuon, // AaronCuts
  					kBrokenRockMuonCut, // AaronCuts
  					kHadronContainment, // AaronCuts
                                        kVtx,
                                        kMinosMuon,
                                        kPmu,
                                        kThetamu,
                                        kAtLeastOneMichel,
					kGoodMomentum, // AaronCuts
                                        kAtLeastOnePionCandidateTrack,
                                        kLLR,
                                        kNode,
					kTpiCut,
                                        kWexp,
                                        kIsoProngs,
                                        kPionMult};

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

    case kMinosMuon:
      return "MINOS Muon";

    case kWexp:
      return "$W_{experimental}$";

    case kIsoProngs:
      return "$<$2 Isolated Prongs";

    case kNPionCandidates:
      return "$\\pi$ candidate";

    case kAtLeastOneMichel:
      return "$>$= 1 Michel";

    case kAtLeastOnePionCandidateTrack:
      return "$>$= 1 Hadron Track";

    case kNode:
      return "Node";

    case kPionMult:
      return "Pion Multiplicity";

    case kLLR:
      return "LLR PID";

    case kAllCuts:
      return "Total";

    case kTrackQuality:
      return "General Track Quality";

    case kPmu:
      return "1.5 GeV $<$ Pmu $<$ 20 GeV";

    case kTpiCut:
      return "35 $<$ Tpi 350 MeV";

    case kThetamu:
      return "$\\theta_\\mu$ $<$ 13";

    case kAtLeastOnePionCandidate:
      return "At Least One Pion";

    default:
      std::cout << "ERROR: GetCutName unknown cut!" << std::endl;
      return "";
  };
}

#endif

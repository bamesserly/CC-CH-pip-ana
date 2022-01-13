//==============================================================================
// I'd probably like to turn this into a class and/or break all of this up into
// smaller parts. But I'm hesitant to mess with it.
//
// How it all works:
// * PassesCuts calls PassesCut on given vector of cuts.
// * A static container of Michel objects is shlepped around, passed to each
// PassesCut.
// * The michels correspond to the pion candidates.
// * One of the cuts fills the container, and subsequent cuts remove
// michels/pions from the container.
// * A container of pion candidate indices corresponding to the michel
// container is returned by reference from PassesCuts.
//==============================================================================

#ifndef Cuts_H
#define Cuts_H

#include <tuple>
#include <vector>

#include "CVUniverse.h"
#include "Constants.h"  // enum ECuts, CCNuPionIncConsts
#include "Michel.h"
#include "SignalDefinition.h"

//==============================================================================
// USE THESE CUTS
//
// Note 0: The michel cut IS the selector of the hadron candidate. Without it,
// you cannot make further cuts on pion variables nor plot pion variables.
// Note 1: The order of precuts matters for some things:
// 1. objects 2. vertex 3. FV 4. Minos activity
//==============================================================================
const std::vector<ECuts> kCutsVector = {
    kNoCuts,
    kPrecuts,
    kVtx,
    kMinosMuon,
    kAtLeastOnePionCandidateTrack,
    kAtLeastOneMichel,
    kLLR,
    kNode,
    kWexp,  // Calling this after pion candidate determined!
    kIsoProngs,
    kPionMult,
    kPmu};

const std::vector<ECuts> GetWSidebandCuts() {
  std::vector<ECuts> w_sideband_cuts = kCutsVector;
  w_sideband_cuts.erase(
      std::remove(w_sideband_cuts.begin(), w_sideband_cuts.end(), kWexp),
      w_sideband_cuts.end());
  return w_sideband_cuts;
}

//==============================================================================
// These cuts were made in the ana tool
// So for reco data and MC these are always true.
//==============================================================================
bool IsPrecut(ECuts c) {
  if (c == kNoCuts || c == kGoodObjects || c == kGoodVertex ||
      c == kFiducialVolume || c == kMinosActivity || c == kPrecuts)
    return true;
  else
    return false;
}

//==============================================================================
// Function Declarations
//==============================================================================
// Call the cut functions -- fill-in/return the ref to the good pion candidate
bool PassesCuts(CVUniverse&, std::vector<int>& pion_candidate_idxs, bool is_mc,
                SignalDefinition, std::vector<ECuts> cuts = kCutsVector);

// also tell whether we are w sideband
bool PassesCuts(CVUniverse&, std::vector<int>& pion_candidate_idxs,
                const bool is_mc, const SignalDefinition, bool& is_w_sideband,
                std::vector<ECuts> cuts = kCutsVector);

// NEW return passes_all_cuts, is_w_sideband, and pion_candidate_indices
std::tuple<bool, bool, std::vector<int>> PassesCuts(
    CVUniverse&, const bool is_mc, const SignalDefinition,
    const std::vector<ECuts> cuts = kCutsVector);

EventCount PassedCuts(const CVUniverse&, std::vector<int>& pion_candidate_idxs,
                      bool is_mc, SignalDefinition,
                      std::vector<ECuts> cuts = kCutsVector);

bool PassesCut(const CVUniverse&, const ECuts cut, const bool is_mc,
               const SignalDefinition, endpoint::MichelMap& endpoint_michels,
               endpoint::MichelMap& vertex_michels);

// Get a vector of integers which are unique hadron prong identifiers.
// The length of this vector is the number of pion candidate tracks found.
// std::vector<int> GetPionCandidates(CVUniverse&);

// Helper
std::string GetCutName(ECuts cut);

// Cuts functions -- eventwide -- TRUTH ONLY
bool GoodObjectsCut(const CVUniverse&);
bool GoodVertexCut(const CVUniverse&);
bool FiducialVolumeCut(const CVUniverse&);
bool MinosActivityCut(const CVUniverse&);

// Cut functions -- eventwide
bool MinosMatchCut(const CVUniverse&);
bool MinosChargeCut(const CVUniverse&);
bool WexpCut(const CVUniverse&, SignalDefinition);
bool IsoProngCut(const CVUniverse&);
std::vector<int> GetQualityPionCandidateIndices(const CVUniverse&);
bool HadronQualityCuts(const CVUniverse&, const RecoPionIdx pion_candidate_idx);
bool vtxCut(const CVUniverse& univ);
bool zVertexCut(const CVUniverse& univ, const double upZ, const double downZ);
bool XYVertexCut(const CVUniverse& univ, const double a);
bool PmuCut(const CVUniverse& univ);

// Cuts functions -- on pion candidate tracks
bool LLRCut(const CVUniverse&, const RecoPionIdx pion_candidate_idx);
bool NodeCut(const CVUniverse&, const RecoPionIdx pion_candidate_idx);

// Retired
bool DeadTimeCut(const CVUniverse&);
bool MinosCoilCut(const CVUniverse&);
bool IsoBlobCut(const CVUniverse&);
bool IsoProngSepCut(const CVUniverse&);
bool ThetamuCut(const CVUniverse&);
bool BrandonMinosChargeCut(const CVUniverse&);
bool CCIncMinosChargeCut(const CVUniverse&);
bool ExactlyOneEndpointMichelCut(const CVUniverse&, SignalDefinition);
bool AtLeastOneBrandonMichelCut(const CVUniverse&);
bool AtLeastOneAnchoredProngCut(const CVUniverse&);
bool AtLeastOneNodeCandidateCut(const CVUniverse&);
bool AtLeastOneLLRCandidateCut(const CVUniverse&);

#endif

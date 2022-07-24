//==============================================================================
// This file contains all event selection cuts definitions.
//
// As well as generic PassesCut and PassesCuts functions.
//
// PassesCuts does more than just check whether events pass all cuts. Checking
// cuts is expensive. So additionally, it checks whether the event-universe is
// w sideband, and it returns the tracks deemed to be pion candidates.
//==============================================================================
#ifndef Cuts_H
#define Cuts_H

#include <tuple>
#include <vector>

#include "CVUniverse.h"
#include "Constants.h"  // enum ECuts, CCNuPionIncConsts
#include "CutUtils.h"
#include "Michel.h"
#include "SignalDefinition.h"

//==============================================================================
// Generic Pass Cut(s) Functions
//      * passes, is_sideband, pion_indices = PassesCuts()
//      * bool PassesCut(cut)
//      * PassedCuts <-- just an event counter
//==============================================================================
// PassesCuts v3 (latest and greatest)
std::tuple<bool, bool, std::vector<int>> PassesCuts(
    CVUniverse&, const bool is_mc, const SignalDefinition,
    const std::vector<ECuts> cuts = kCutsVector);

// Event Counter
EventCount PassedCuts(const CVUniverse&, std::vector<int>& pion_candidate_idxs,
                      bool is_mc, SignalDefinition,
                      std::vector<ECuts> cuts = kCutsVector);

// Passes Single, Given Cut
std::tuple<bool, endpoint::MichelMap, trackless::MichelEvent> PassesCut(
    const CVUniverse& univ, const ECuts cut, const bool is_mc,
    const SignalDefinition signal_definition);

//==============================================================================
// Cuts Definitions
//==============================================================================
// Gaudi tool cuts -- read from truth tuple.
// (won't work if we pass a reco mc universe)
bool GoodObjectsCut(const CVUniverse&);
bool GoodVertexCut(const CVUniverse&);
bool FiducialVolumeCut(const CVUniverse&);
bool MinosActivityCut(const CVUniverse&);

// Cut Definitions -- eventwide
bool MinosMatchCut(const CVUniverse&);
bool MinosChargeCut(const CVUniverse&);
bool WexpCut(const CVUniverse&, SignalDefinition);
bool IsoProngCut(const CVUniverse&);
bool vtxCut(const CVUniverse& univ);
bool zVertexCut(const CVUniverse& univ, const double upZ, const double downZ);
bool XYVertexCut(const CVUniverse& univ, const double a);
bool PmuCut(const CVUniverse& univ);

// Cuts Definitions -- exclusive, i.e. on pion candidate tracks
bool HadronQualityCuts(const CVUniverse&, const RecoPionIdx pion_candidate_idx);
bool LLRCut(const CVUniverse&, const RecoPionIdx pion_candidate_idx);
bool NodeCut(const CVUniverse&, const RecoPionIdx pion_candidate_idx);

//==============================================================================
// Helper
//==============================================================================
// Get candidate pions that pass the minimal HadronQualityCuts
std::vector<int> GetQualityPionCandidateIndices(const CVUniverse&);

// bool AtLeastOnePionCut(const CVUniverse& univ) {
//  std::tuple<> GetAllMichels();
//}

//==============================================================================
// Retiring
//==============================================================================

// PassesCuts v2 (being deprecated)
bool PassesCuts(CVUniverse&, std::vector<int>& pion_candidate_idxs, bool is_mc,
                SignalDefinition, std::vector<ECuts> cuts = kCutsVector);

// PassesCuts v1 (being deprecated)
bool PassesCuts(CVUniverse&, std::vector<int>& pion_candidate_idxs,
                const bool is_mc, const SignalDefinition, bool& is_w_sideband,
                std::vector<ECuts> cuts = kCutsVector);

// PassesCut v1 (being deprecated)
bool PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
               const SignalDefinition signal_definition,
               endpoint::MichelMap& endpoint_michels,
               endpoint::MichelMap& vtx_michels);

#endif

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
#include "Constants.h"        // enum ECuts, CCNuPionIncConsts, PassesCutsInfo
#include "CutUtils.h"         // kCutsVector
#include "Michel.h"           // endpoint::Michel, endpoint::MichelMap
#include "MichelTrackless.h"  // trackless::MichelEvent
#include "SignalDefinition.h"

//==============================================================================
// Generic Pass Cut(s) Functions
//      * PassesCutsInfo(passes, is_sideband, all_except_w, pion_idxs) =
//      PassesCuts()
//      * tuple(passes, endpoint_michels, vtx_michels) = PassesCut(cut)
//      * PassedCuts <-- just an event counter
//==============================================================================
// NEW return passes_all_cuts, is_w_sideband, and pion_candidate_idxs
// PassesCuts v3 (latest and greatest))
PassesCutsInfo PassesCuts(CVUniverse&, const bool is_mc, const SignalDefinition,
                          const std::vector<ECuts> cuts = kCutsVector);

// Event Counter
EventCount PassedCuts(const CVUniverse&, std::vector<int>& pion_candidate_idxs,
                      bool is_mc, SignalDefinition,
                      std::vector<ECuts> cuts = kCutsVector);

// Passes Single, Given Cut
// New, to be implemented.
std::tuple<bool, endpoint::MichelMap, trackless::MichelEvent<CVUniverse>>
PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
          const SignalDefinition, const endpoint::MichelMap&,
          const trackless::MichelEvent<CVUniverse>&);

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
bool ThetamuCut(const CVUniverse& univ);

// Cuts Definitions -- exclusive, i.e. on pion candidate tracks
bool HadronQualityCuts(const CVUniverse&, const RecoPionIdx pion_candidate_idx);
bool LLRCut(const CVUniverse&, const RecoPionIdx pion_candidate_idx);
bool NodeCut(const CVUniverse&, const RecoPionIdx pion_candidate_idx);
bool tpiCut(const CVUniverse&, const RecoPionIdx pion_candidate_idx);

//==============================================================================
// Helper
//==============================================================================
// Get candidate pions that pass the minimal HadronQualityCuts
std::vector<int> GetQualityPionCandidateIndices(const CVUniverse&);

std::vector<int> GetHadIdxsFromMichels(
    const endpoint::MichelMap endpoint_michels,
    const trackless::MichelEvent<CVUniverse> vtx_michels =
        trackless::MichelEvent<CVUniverse>());

// bool AtLeastOnePionCut(const CVUniverse& univ) {
//  std::tuple<> GetAllMichels();
//}

//==============================================================================
// BEING DEPRECATED
//==============================================================================

// Call the cut functions -- fill-in/return the ref to the good pion candidate
// PassesCuts v1 (being deprecated))
bool PassesCuts(CVUniverse&, std::vector<int>& pion_candidate_idxs, bool is_mc,
                SignalDefinition, std::vector<ECuts> cuts = kCutsVector);

// also tell whether we are w sideband
// PassesCuts v2 (being deprecated)
bool PassesCuts(CVUniverse&, std::vector<int>& pion_candidate_idxs,
                const bool is_mc, const SignalDefinition, bool& is_w_sideband,
                std::vector<ECuts> cuts = kCutsVector);

// PassesCut v1 (being deprecated))
bool PassesCut(const CVUniverse&, const ECuts cut, const bool is_mc,
               const SignalDefinition, endpoint::MichelMap& endpoint_michels,
               endpoint::MichelMap& vertex_michels);

#endif

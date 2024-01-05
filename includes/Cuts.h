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

#include "CCPiEvent.h"
#include "CVUniverse.h"
#include "Constants.h"  // enum ECuts, CCNuPionIncConsts, PassesCutsInfo
#include "CutUtils.h"   // kCutsVector
#include "Michel.h"     // endpoint::Michel, endpoint::MichelMap
#include "PlotUtils/LowRecoilPionCuts.h"
#include "PlotUtils/LowRecoilPionReco.h"
#include "SignalDefinition.h"

bool PassesCuts(const CCPiEvent&, const SignalDefinition&,
                const std::vector<ECuts>&);
bool PassesCut(const CCPiEvent&, const ECuts&, const SignalDefinition&);

enum class VerticalUniverseCheck { kUNCHECKED, kPASS, kFAIL };

std::tuple<bool, VerticalUniverseInfo> PassesCuts(
    const CCPiEvent& event, const SignalDefinition& signal_definition,
    const std::vector<ECuts>& cuts, const VerticalUniverseInfo& in_vert_info);

/*
// I hate these functions
tuple<bool, VerticalUniverseCheck> CheckBasicCuts(
    const CCPiEvent& event, const bool is_mc,
    const SignalDefinition& signal_definition,
    const VerticalUniverseCheck vert_checked_status) {
  if (event.m_universe->IsVerticalOnly() &&
      vert_checked_status != VerticalUniverseCheck::kUNCHECKED) {
    return {vert_checked_status == VerticalUniverseCheck::kPASS ? true : false,
            vert_checked_status};
  }

  bool passes_basic_cuts =
      PassesCuts(event, is_mc, signal_definition, GetBasicCuts());

  VerticalUniverseCheck updated_status;
  if (event.m_universe->IsVerticalOnly())
    updated_status = passes_basic_cuts ? VerticalUniverseCheck::kPASS
                                       : VerticalUniverseCheck::kFAIL;
  else
    updated_status = vert_checked_status;

  // Can't decide which of these I hate more.
  // VerticalUniverseCheck updated_status =
  //    !event.m_universe->IsVerticalOnly()
  //        ? vert_checked_status
  //        : passes_basic_cuts ? VerticalUniverseCheck::kPASS
  //                            : VerticalUniverseCheck::kFAIL;

  return {passes_basic_cuts, updated_status};
}

tuple<std::vector<PionCandidate>, std::vector<PionCandidate>, bool>
GetPionCandidates(const CCPiEvent& event, const bool is_mc,
                  const SignalDefinition& signal_definition,
                  const std::vector<PionCandidate> vert_pion_candidates,
                  const bool made_vert_pion_candidates) {
  if (event.m_universe->IsVerticalOnly() && made_vert_pion_candidates)
    return {vert_pion_candidates, vert_pion_candidates,
            made_vert_pion_candidate};

  std::vector<PionCandidate> pion_candidates =
      MakePassingPionCandidates(event, is_mc, signal_definition);

  bool updated_status;
  std::vector<PionCandidate> updated_vert_pion_candidates;
  if (event.m_universe->IsVerticalOnly()) {
    updated_status = true;
    updated_vert_pion_candidates = pion_candidates;
  } else {
    updated_status = made_vert_pion_candidate;
  }

  return {pion_candidates, updated_vert_pion_candidates, updated_status};
}
*/

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
                      bool is_mc, const SignalDefinition,
                      std::vector<ECuts> cuts = kCutsVector);

// Passes Single, Given Cut
// New, to be implemented.
std::tuple<bool, endpoint::MichelMap, LowRecoilPion::MichelEvent<CVUniverse>>
PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
          const SignalDefinition, const endpoint::MichelMap&,
          const LowRecoilPion::MichelEvent<CVUniverse>&);

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
bool WexpCut(const CVUniverse&, const SignalDefinition);
bool IsoProngCut(const CVUniverse&);
bool VtxInFiducialCut(const CVUniverse& univ);
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

std::vector<int> GetHadIdxsFromMichels(
    const endpoint::MichelMap endpoint_michels,
    const LowRecoilPion::MichelEvent<CVUniverse> vtx_michels =
        LowRecoilPion::MichelEvent<CVUniverse>());

// bool AtLeastOnePionCut(const CVUniverse& univ) {
//  std::tuple<> GetAllMichels();
//}

#endif

#ifndef CCPiEvent_cxx
#define CCPiEvent_cxx

#include "CCPiEvent.h"
#include "Cuts.h"
#include "utilities.h" // ContainerEraser
#include "PlotUtils/LowRecoilPionCuts.h"
#include "PlotUtils/LowRecoilPionReco.h"

std::tuple<UntrackedMichels, bool> GetUntrackedMichels(const CVUniverse* universe) {
  UntrackedMichels untracked_michels;
  bool pass = hasMichel::hasMichelCut(*universe, untracked_michels);
  pass = pass && BestMichelDistance2D::BestMichelDistance2DCut(*universe, untracked_michels);
  pass = pass && GetClosestMichel::GetClosestMichelCut(*universe, untracked_michels);
  return {untracked_michels, pass};
}

endpoint::MichelMap GetTrackedPionCandidates(const CVUniverse* universe) {
  const CVUniverse& univ = *universe;
  // Get michels subject to prior quality cuts
  endpoint::MichelMap michels = endpoint::GetQualityMichels(univ);

  // TODO assert that no two michels in this map have the same idx.

  // Remove michels if their associated pion track fails cuts

  // Basic hadron track quality
  ContainerEraser::erase_if(
      michels, [&univ](std::pair<int, endpoint::Michel> mm) {
        return !HadronQualityCuts(univ, mm.second.had_idx);
      });

  // LLR -- proton vs pion separation PID score
  ContainerEraser::erase_if(michels,
                            [&univ](std::pair<int, endpoint::Michel> mm) {
                              return !LLRCut(univ, mm.second.had_idx);
                            });

  // Node cut -- remove interacting pions (which have bad tpi reco)
  ContainerEraser::erase_if(michels,
                            [&univ](std::pair<int, endpoint::Michel> mm) {
                              return !NodeCut(univ, mm.second.had_idx);
                            });

  return michels;

  /*
    pass = signal_definition.m_n_pi_min <= michels.size() &&
           michels.size() <= signal_definition.m_n_pi_max;

    // Subject trackless michels to cuts
    bool pass_tracked_pion_cuts = PassesCuts(*this, m_signal_definition, GetTrackedPionCuts());
    //m_universe->m_michels = m_michels;



    // modify michels
    // the quality track, LLR, and node cuts may have removed the michels of
    // failed tracks
    case kAtLeastOnePionCandidate:
      pass = michels.size() > 0;  // || vtx_michels.m_idx != -1;
      break;

    // TODO check trackless michels and no double-counting
    case kPionMult: {
      break;
    }

    case kAtLeastOnePionCandidateTrack:
      pass = GetQualityPionCandidateIndices(univ).size() >
             0;  // ||
                 // vtx_michels.m_idx != -1;
      break;
  */
}

struct PionCandidates {
  UntrackedMichels untracked_michels;
  endpoint::MichelMap tracked_michels;
};

//==============================================================================
// CTOR
//==============================================================================
CCPiEvent::CCPiEvent(const bool is_mc, const bool is_truth,
                     const SignalDefinition signal_definition,
                     CVUniverse* universe)
    : m_is_mc(is_mc),
      m_is_truth(is_truth),
      m_signal_definition(signal_definition),
      m_universe(universe),
      m_is_signal(is_mc ? IsSignal(*universe, signal_definition) : false),
      m_w_type(is_mc ? GetWSidebandType(*universe, signal_definition,
                                        sidebands::kNWFitCategories)
                     : kNWSidebandTypes),
      m_weight(is_mc ? universe->GetWeight() : 1.) {}


//==============================================================================
// Computationally expensive reconstruction and cuts checking.
//
// This function updates the internal state of 'this' event.
//
// It is extremely not const.
//
// N.B. due to the expensive nature of these steps, quit the function as we go
// along if the event fails basic cuts (first), pion cuts (second), w cut
// (third).
//==============================================================================
std::tuple<VerticalUniverseInfo, LoopStatusCode> CCPiEvent::Process(
    const VerticalUniverseInfo& in_vert_info) {
  VerticalUniverseInfo vert_info = in_vert_info;
  bool is_vert_only = m_universe->IsVerticalOnly();

  // Basic cuts
  if (is_vert_only) {
    if (vert_info.checked_basic_cuts) {
      m_passes_cuts = vert_info.passes_basic_cuts;
    } else {
      m_passes_cuts = vert_info.passes_basic_cuts =
          PassesCuts(*this, m_signal_definition, GetBasicCuts());
      vert_info.checked_basic_cuts = true;
    }
  } else {
    m_passes_cuts = PassesCuts(*this, m_signal_definition, GetBasicCuts());
  }

  // TODO could apply this cut as a basic cut if our sig def requires it
  //kAtLeastOnePionCandidateTrack

  // Fail basic cuts? EXIT
  if (!m_passes_cuts) return {vert_info, LoopStatusCode::SKIP};

  // Build trackless michels
  //
  // In Mehreen's system, there is both a bool pass and michel objects. Maybe
  // you can't have a pass without michels, but the converse isn't necessarily
  // true.
  if (m_signal_definition.m_do_untracked_michel_reco) {
    if (is_vert_only) {
      if (vert_info.made_untracked_michels) {
        m_untracked_michels = vert_info.untracked_michels;
        m_untracked_michels_pass = vert_info.untracked_michels_pass;
      } else {
        std::tie(m_untracked_michels, m_untracked_michels_pass) = GetUntrackedMichels(m_universe);
        vert_info.untracked_michels = m_untracked_michels;
        vert_info.untracked_michels_pass = m_untracked_michels_pass;
        vert_info.made_untracked_michels = true;
      }
    } else {
      std::tie(m_untracked_michels, m_untracked_michels_pass) = GetUntrackedMichels(m_universe);
    }
  }

  // Build tracked michels
  //
  // In my system pass iff michel objects exist.
  if (m_signal_definition.m_do_tracked_michel_reco) {
    m_tracked_michels = GetTrackedPionCandidates(m_universe);
  }

  // TODO assert unique tracked pion michel indices

  // tracked only
  if(m_signal_definition.m_do_tracked_michel_reco && !m_signal_definition.m_do_untracked_michel_reco) {
    m_passes_cuts = m_tracked_michels.empty(); 
    if (!m_passes_cuts) return {vert_info, LoopStatusCode::SKIP};
  }
  // untracked only
  if(!m_signal_definition.m_do_tracked_michel_reco && m_signal_definition.m_do_untracked_michel_reco) {
    m_passes_cuts = m_untracked_michels_pass;
    if (!m_passes_cuts) return {vert_info, LoopStatusCode::SKIP};
  }

  // tracked and untracked
  assert(m_signal_definition.m_do_tracked_michel_reco && m_signal_definition.m_do_untracked_michel_reco);

  // Failed to find any pion candidates. EXIT.
  if(m_tracked_michels.empty() && !m_untracked_michels_pass) {
    m_passes_cuts = false;
    return {vert_info, LoopStatusCode::SKIP};
  }

  // We got at least one pion candidate from one or both of these methods

  // if we got pion candidates from BOTH methods, look for re-use of the same michel cluster
  if(!m_tracked_michels.empty() && m_untracked_michels_pass) {
    // Best untracked michel -- unique michel index
    // TODO check multiple passing michels, not just the best
    int unique_michel_idx_untracked = m_untracked_michels.m_nmichels[m_untracked_michels.m_idx].tuple_idx;

    // potentially multiple tracked pion candidates
    std::vector<int> unique_michel_idx_tracked;
    for (auto candidate : m_tracked_michels) {
      unique_michel_idx_tracked.push_back(candidate.first);
    }

    // is the untracked michel found among the tracked michels?
    auto it = std::find(unique_michel_idx_tracked.begin(), unique_michel_idx_tracked.end(), unique_michel_idx_untracked);

    if (it != unique_michel_idx_tracked.end()) {
      // the best untracked michel candidate was also matched to a track.
      int index = std::distance(unique_michel_idx_tracked.begin(), it);
      std::cout << "The integer " << unique_michel_idx_untracked << " was found at index: " << index << ". " << m_tracked_michels.size() << std::endl;
    } else {
      // the best untracked michel candidate was not among the michels matched to a track, this is an additional pion candidate!
      std::cout << "The integer " << unique_michel_idx_untracked<< " was not found in the vector." << ". " << m_tracked_michels.size() << std::endl;
    }
  }

  /*
  */

  /*
    // Contruct pions -- the most computationally expensive
    std::tie(m_pion_candidates, ret_vert_info) = GetPionCandidates(event,
    is_mc, signal_definition, ret_vert_info);

    // No pions -- EXIT
    if (!PassesCut(event, ECuts::kPionMult))
      return {ret_vert_info, LoopStatusCode::SKIP};

    // Check W cut and sideband
    m_passes_all_cuts_except_w = true;
    bool passes_w_cut = false;
    std::tie(passes_w_cut, m_is_w_sideband) = PassesCut(event, kWexp);
    m_passes_cuts = passes_w_cut ? true : false;

    // Fail W and not sideband -- EXIT
    if (!passes_w_cut && !m_is_w_sideband)
      return {ret_vert_info, LoopStatusCode::SKIP};
  */

  // SUCCESS
  return {vert_info, LoopStatusCode::SUCCESS};
}

double CCPiEvent::GetDummyVar() const { return -99.; }

#endif  // CCPiEvent_cxx

#ifndef CCPiEvent_cxx
#define CCPiEvent_cxx

#include "CCPiEvent.h"

#include "Cuts.h"
#include "PlotUtils/LowRecoilPionCuts.h"
#include "PlotUtils/LowRecoilPionReco.h"
#include "utilities.h"  // ContainerEraser

int n_more_untracked = 0;
int n_more_tracked = 0;
int n_same = 0;
int n_tracked_multipi = 0;
int n_untracked_multipi = 0;
int n_multipi_agree = 0;

std::tuple<UntrackedMichels, bool> GetUntrackedMichels(
    const CVUniverse* universe) {
  UntrackedMichels untracked_michels;
  bool pass = hasMichel::hasMichelCut(*universe, untracked_michels);
  pass = pass && BestMichelDistance2D::BestMichelDistance2DCut(
                     *universe, untracked_michels);
  pass = pass &&
         GetClosestMichel::GetClosestMichelCut(*universe, untracked_michels);
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
    bool pass_tracked_pion_cuts = PassesCuts(*this, m_signal_definition,
    GetTrackedPionCuts());
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

  //============================================================================
  // Basic cuts (EXIT if fail)
  //============================================================================
    bool do_check_cuts = (!vert_info.checked_basic_cuts && is_vert_only) || !is_vert_only;
    bool passes_basic_cuts = false;
    if (do_check_cuts) {
      passes_basic_cuts = PassesCuts(*this, m_signal_definition, GetBasicCuts());
      if (is_vert_only) {
        vert_info.checked_basic_cuts = true;
        vert_info.passes_basic_cuts = passes_basic_cuts;
      }
    }
    if (is_vert_only) {
      assert(vert_info.checked_basic_cuts);
      passes_basic_cuts = vert_info.passes_basic_cuts;
    }
    m_passes_cuts = passes_basic_cuts;

    // TODO could apply this cut as a basic cut if our sig def requires it
    // kAtLeastOnePionCandidateTrack

  // Fail basic cuts? EXIT
  if (!m_passes_cuts) return {vert_info, LoopStatusCode::SKIP};

  // TODO assert unique tracked pion michel indices

  //============================================================================
  // BUILD tracked michels if the SD calls for it.
  // We've already built untracked michels, if needed, outside of this function
  // -- it's the same for all universes. They're saved to the event.
  //============================================================================
    ////============================================================================
    //// Build trackless michels
    ////============================================================================
    //  // In Mehreen's system, there is both a bool pass and michel objects. We
    //  can
    //  // fail and still have michel objects.
    //  if (m_signal_definition.m_do_untracked_michel_reco) {
    //    if (is_vert_only) {
    //      if (!vert_info.made_untracked_michels) {
    //        std::tie(m_untracked_michels, m_untracked_michels_pass) =
    //        GetUntrackedMichels(m_universe); vert_info.untracked_michels =
    //        m_untracked_michels; vert_info.untracked_michels_pass =
    //        m_untracked_michels_pass; vert_info.made_untracked_michels = true;
    //      } else {
    //        m_untracked_michels = vert_info.untracked_michels;
    //        m_untracked_michels_pass = vert_info.untracked_michels_pass;
    //      }
    //    } else {
    //      std::tie(m_untracked_michels, m_untracked_michels_pass) =
    //      GetUntrackedMichels(m_universe);
    //    }
    //  }

    //==========================================================================
    // Build tracked michels
    //==========================================================================
    // In my system pass iff michel objects exist.
    if (m_signal_definition.m_do_tracked_michel_reco) {
      m_tracked_michels = GetTrackedPionCandidates(m_universe);
    }

  //==========================================================================
  // REQUIRE at least one michel? EXIT if none
  //==========================================================================
  // tracked only analysis
  if (m_signal_definition.m_do_tracked_michel_reco &&
      !m_signal_definition.m_do_untracked_michel_reco) {
    m_passes_cuts = !m_tracked_michels.empty();
  }
  // untracked only analysis
  else if (!m_signal_definition.m_do_tracked_michel_reco &&
           m_signal_definition.m_do_untracked_michel_reco) {
    assert(m_untracked_michels_pass && !m_untracked_michels.m_nmichels.empty());
    m_passes_cuts = m_untracked_michels_pass;
  }
  // tracked and untracked analysis
  else {
    assert(m_signal_definition.m_do_tracked_michel_reco &&
           m_signal_definition.m_do_untracked_michel_reco);
    m_passes_cuts = !m_tracked_michels.empty() || (m_untracked_michels_pass && !m_untracked_michels.m_nmichels.empty());
  }

  if (!m_passes_cuts) return {vert_info, LoopStatusCode::SKIP};

  // If we're here, we have a signal michel candidate

  //============================================================================
  // Begin 2023-10-03 study 
  //============================================================================
  // If we're not analyzing both SD's, then quit.
  if (!m_signal_definition.m_do_tracked_michel_reco ||
      !m_signal_definition.m_do_untracked_michel_reco) std::exit(9);

  //============================================================================
  // Let's just count who found how many michels
  // And multipi?
  //============================================================================
    int n_untracked_michels = m_untracked_michels.m_nmichels.size();
    int n_tracked_michels = m_tracked_michels.size();

    if (n_untracked_michels > n_tracked_michels) n_more_untracked ++;
    if (n_tracked_michels > n_untracked_michels) n_more_tracked ++;
    if (n_tracked_michels == n_untracked_michels) n_same++;
    if (n_untracked_michels > 2) { n_untracked_multipi++; if(m_is_mc) m_universe->PrintArachneLink(); if(m_is_mc)std::cout << m_universe->GetNChargedPionsTrue() << "\n";}
    if (n_tracked_michels > 2) n_tracked_multipi++;
    if (n_tracked_michels > 2 && n_untracked_multipi > 2) n_multipi_agree++;
  //============================================================================

  

  /*
  //============================================================================
  // If using both reco michel methods, and both found a michel:
  // check for duplicates (i.e. re-use of the same cluster)
  //============================================================================
  if (m_signal_definition.m_do_tracked_michel_reco &&
      m_signal_definition.m_do_untracked_michel_reco &&
      !m_tracked_michels.empty() && m_untracked_michels_pass &&
      !m_untracked_michels.m_nmichels.empty()) {
    //==========================================================================
    // Look for the best untracked michel among the tracked michels.
    //==========================================================================
    // Get the unique michel index of the best untracked michel
    int unique_michel_idx_untracked =
        m_untracked_michels.m_nmichels[m_untracked_michels.m_idx].tuple_idx;

    // TODO check multiple passing michels, not just the best

    // Get the unique indices for the michels matched to each of the tracked
    // pion candidates.
    std::vector<int> unique_michel_idx_tracked;
    for (auto candidate : m_tracked_michels) {
      unique_michel_idx_tracked.push_back(candidate.first);
    }

    // Search for the one untracked michel among the tracked michels
    auto it =
        std::find(unique_michel_idx_tracked.begin(),
                  unique_michel_idx_tracked.end(), unique_michel_idx_untracked);

    if (it != unique_michel_idx_tracked.end()) {
      // the best untracked michel candidate was also matched to a track.
      int index = std::distance(unique_michel_idx_tracked.begin(), it);
      std::cout << "The integer " << unique_michel_idx_untracked
                << " was found at index: " << index << ". "
                << m_tracked_michels.size() << std::endl;
    } else {
      // the best untracked michel candidate was not among the michels matched
      // to a track, this is an additional pion candidate!
      std::cout << "The integer " << unique_michel_idx_untracked
                << " was not found in the vector."
                << ". " << m_tracked_michels.size() << std::endl;
    }
  }
  */

  ////============================================================================
  //// Track Cuts
  ////============================================================================
  // if (is_vert_only) {
  //  if (vert_info.checked_track_cuts) {
  //    m_passes_cuts = vert_info.passes_track_cuts;
  //  } else {
  //    m_passes_track_cuts = vert_info.passes_track_cuts = PassesCuts(*this,
  //    m_signal_definition, GetTrackCuts()); vert_info.checked_track_cuts =
  //    true;
  //  }
  //} else {
  //  m_passes_track_cuts = PassesCuts(*this, m_signal_definition,
  //  GetTrackCuts());
  //}
  // m_passes_cuts = m_passes_cuts && m_passes_track_cuts;

  // I'll come back to this
  ////============================================================================
  //// W Cut
  ////============================================================================
  //// Note: Performing the W cut sneakily bypasses the const-ness of the input
  //// Event to this process function. Namely, while returning the pass bool, it
  //// also sets the event's `m_is_w_sideband` member.
  // if(m_signal_definition.)
  // if (is_vert_only) {
  //  if (vert_info.checked_w_cut) {
  //    m_passes_cuts = vert_info.passes_w_cut;
  //  } else {
  //    m_passes_w_cut = vert_info.passes_w_cuts = PassesCuts(*this,
  //    m_signal_definition, GetWCut()); vert_info.checked_w_cuts = true;
  //  }
  //} else {
  //  m_passes_w_cuts = PassesCuts(*this, m_signal_definition, GetWCut());
  //}
  // m_passes_cuts = m_passes_cuts && m_passes_w_cuts;
  // assert(m_wexp != -99. && "Calling the W cut failed to set
  // CCPiEvent::m_wexp");

  //// Fail W and not sideband -- EXIT
  // if (!m_passes_w_cut && !m_is_w_sideband)
  //  return {ret_vert_info, LoopStatusCode::SKIP};

  // m_is_w_sideband

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

  */

  // SUCCESS
  return {vert_info, LoopStatusCode::SUCCESS};
}

// Like FillCutVars, this function loops through cuts and calls PassesCut.
// Michel containers updated as we go, but thrown away at the end.
namespace ccpi_event{
  std::pair<EventCount, EventCount> FillCounters(const CCPiEvent& event, const EventCount& s, const EventCount& b) {
    EventCount signal = s;
    EventCount bg = b;

    bool pass = true;
    for (auto i_cut : kCutsVector) {
      if (event.m_is_truth != IsPrecut(i_cut)) continue;

      bool passes_this_cut = PassesCut(event, i_cut, event.m_signal_definition);

      pass = pass && passes_this_cut;

      if (!pass) break;
      if (!event.m_is_mc) {
        signal[i_cut] += event.m_weight;  // selected data
      } else {
        if (event.m_is_signal)
          signal[i_cut] += event.m_weight;  // selected mc signal
        else
          bg[i_cut] += event.m_weight;  // selected mc bg
      }
    }  // cuts loop
    return {signal, bg};
  }
}

double CCPiEvent::GetDummyVar() const { return -99.; }

#endif  // CCPiEvent_cxx

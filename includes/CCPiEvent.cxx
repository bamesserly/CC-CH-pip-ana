#ifndef CCPiEvent_cxx
#define CCPiEvent_cxx

#include "CCPiEvent.h"

#include "Cuts.h"
#include "PlotUtils/LowRecoilPionCuts.h"
#include "PlotUtils/LowRecoilPionReco.h"

std::tuple<TracklessMichels, bool> GetTracklessMichels(const CVUniverse* universe) {
  TracklessMichels trackless_michels;
  bool pass = hasMichel::hasMichelCut(*universe, trackless_michels);
  pass = pass && BestMichelDistance2D::BestMichelDistance2DCut(*universe, trackless_michels);
  pass = pass && GetClosestMichel::GetClosestMichelCut(*universe, trackless_michels);
  return {trackless_michels, pass};
}

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

  // Fail basic cuts -- EXIT
  if (!m_passes_cuts) return {vert_info, LoopStatusCode::SKIP};

  // Build trackless michels
  if (is_vert_only) {
    if (vert_info.made_trackless_michels) {
      m_trackless_michels = vert_info.trackless_michels;
      m_trackless_michels_pass = vert_info.trackless_michels_pass;
    } else {
      std::tie(m_trackless_michels, m_trackless_michels_pass) = GetTracklessMichels(m_universe);
      vert_info.trackless_michels = m_trackless_michels;
      vert_info.trackless_michels_pass = m_trackless_michels_pass;
      vert_info.made_trackless_michels = true;
    }
  } else {
    std::tie(m_trackless_michels, m_trackless_michels_pass) = GetTracklessMichels(m_universe);
  }

  // Build tracked michels
  m_endpoint_michels = endpoint::GetQualityMichels(m_univ);
  m_universe->m_endpoint_michels = m_endpoint_michels;

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

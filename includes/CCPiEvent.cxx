#ifndef CCPiEvent_cxx
#define CCPiEvent_cxx

#include "CCPiEvent.h"

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
// PROCESS
// Macro-level reco and cuts-checking.
// Update internal state of 'this' event. Extremely not const.
// The most computationally expensive thing we do.
//==============================================================================
std::tuple<VertUniverseInfo, LoopStatusCode> CCPiEvent::Process(const VertUniverseInfo& vert_info) {
  VertUniverseInfo ret_vert_info = vert_info;
  /*
    // Basic, event-wide, quick cuts
    std::tie(m_passes_basic_cuts, ret_vert_info) = CheckBasicCuts(event, is_mc, signal_definition, ret_vert_info);

    // Fail basic cuts -- EXIT
    if (!m_passes_basic_cuts)
      return {ret_vert_info, LoopStatusCode::SKIP};

    // Contruct pions -- the most computationally expensive
    std::tie(m_pion_candidates, ret_vert_info) = GetPionCandidates(event, is_mc, signal_definition, ret_vert_info);

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
  return {ret_vert_info, LoopStatusCode::SUCCESS};
}

double CCPiEvent::GetDummyVar() const { return -99.; }

#endif  // CCPiEvent_cxx

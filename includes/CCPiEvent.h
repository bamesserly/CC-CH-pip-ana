//==============================================================================
// CCPiEvent is a container struct holding misc info about a universe-event:
// * passes cuts
// * is signal
// * event weight
// * is mc, is truth
// * vector of candidate pion indices
// * signal definition currently being used
// * whether and what kind of w sideband it is
// And it owns a pointer to the corresponding CVUniverse.
//
// This container exists to be passed around to functions that need all this
// info: functions that fill histograms, check whether the event-universe
// passed cuts, etc.
//
// CCPiEvent has no functions of its own.
//==============================================================================
#ifndef CCPiEvent_h
#define CCPiEvent_h

#include "CVUniverse.h"
#include "SignalDefinition.h"
#include "TruthCategories/Sidebands.h" // WSidebandType

struct VertUniverseInfo {
  bool checked_basic_cuts = false;
  bool passes_basic_cuts = false;
  bool made_pion_candidates = false;
  //std::vector<PionCandidate> pion_candidates;
};

enum class LoopStatusCode {SUCCESS, SKIP};

class CCPiEvent {
 public:
  CCPiEvent(const bool is_mc, const bool is_truth,
            const SignalDefinition signal_definition, CVUniverse* universe);

  // Fixed by the constructor
  const bool m_is_mc;
  const bool m_is_truth;
  const SignalDefinition m_signal_definition;
  CVUniverse* m_universe;
  bool m_is_signal;
  WSidebandType m_w_type;
  double m_weight;

  // Process does reco and makes cuts to fill these additional variables
  std::tuple<VertUniverseInfo, LoopStatusCode> Process(const VertUniverseInfo& vert_info);
  bool m_passes_cuts;
  bool m_is_w_sideband;
  bool m_passes_all_cuts_except_w;

  virtual double GetDummyVar() const;

  // New
  //std::vector<PionCandidate> m_pion_candidates;
};


#endif  // CCPiEvent

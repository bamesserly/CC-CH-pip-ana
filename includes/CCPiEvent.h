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
#include "Michel.h"
#include "SignalDefinition.h"
#include "TruthCategories/Sidebands.h"  // WSidebandType
#include "Constants.h" // UntrackedMichels

struct VerticalUniverseInfo {
  VerticalUniverseInfo() :
    checked_basic_cuts(),
    passes_basic_cuts(),
    made_pion_candidates(),
    untracked_michels(),
    untracked_michels_pass(),
    made_untracked_michels(),
    checked_track_cuts(),
    passes_track_cuts(),
    wexp(-99.)
  {}
  bool checked_basic_cuts;
  bool passes_basic_cuts;
  bool made_pion_candidates;
  UntrackedMichels untracked_michels;
  bool untracked_michels_pass;
  bool made_untracked_michels;
  bool checked_track_cuts;
  bool passes_track_cuts;
  double wexp;
  // std::vector<PionCandidate> pion_candidates;
};

enum class LoopStatusCode { SUCCESS, SKIP };

std::tuple<UntrackedMichels, bool> GetUntrackedMichels(const CVUniverse*);
std::tuple<endpoint::MichelMap, bool> GetTrackedMichels(const CVUniverse*);

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

  // Additional, heavy-lifting function to do reco and check cuts
  std::tuple<VerticalUniverseInfo, LoopStatusCode> Process(
      const VerticalUniverseInfo& vert_info);

  // Trackless michels
  UntrackedMichels m_untracked_michels;
  bool m_untracked_michels_pass;

  // Tracked Michels
  endpoint::MichelMap m_tracked_michels;

  // cuts and sideband status
  bool m_passes_cuts;
  bool m_passes_track_cuts;
  double m_wexp; // don't call the cvu function more than needed
  mutable bool m_is_w_sideband; // cuts functions will update this as a side-effect, hence the mutable.
  bool m_passes_all_cuts_except_w;

  virtual double GetDummyVar() const;

  // New
  // std::vector<PionCandidate> m_pion_candidates;

  // GetPionVariables functions
  // double GetTpi(); ...
};

#endif  // CCPiEvent

#ifndef CCPiEvent_h
#define CCPiEvent_h

#include "CVUniverse.h"
#include "Constants.h"  // typedef RecoPionIdx, EventCount
#include "Michel.h"
#include "SignalDefinition.h"
#include "TruthCategories/Sidebands.h"         // WSidebandType
#include "TruthCategories/SignalBackground.h"  // SignalBackgroundType
#include "Variable.h"
class Variable;

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

struct CCPiEvent {
  CCPiEvent(const bool is_mc, const bool is_truth,
            const SignalDefinition signal_definition, CVUniverse* universe);

  // Fixed by the constructor
  const bool m_is_mc;
  const bool m_is_truth;
  const SignalDefinition m_signal_definition;
  CVUniverse* m_universe;
  std::vector<RecoPionIdx>
      m_reco_pion_candidate_idxs;  // initialized empty, filled by PassesCuts
  bool m_is_signal;
  double m_weight;
  WSidebandType m_w_type;

  endpoint::MichelMap m_endpoint_michels;
  vertex::MichelEvent m_vertex_michels;

  // Fixed (directly) outside of constructor -- with time-intensive functions
  bool m_passes_cuts;                     // PassesCuts
  bool m_is_w_sideband;                   // IsWSideband
  RecoPionIdx m_highest_energy_pion_idx;  // GetHighestEnergyPionCandidateIndex

};

// Helper Functions
// bool IsWSideband(CCPiEvent&);
bool PassesCuts(CCPiEvent&, bool& is_w_sideband);
bool PassesCuts(CCPiEvent&, std::vector<ECuts>);
std::tuple<bool, bool, std::vector<int>> PassesCuts(CCPiEvent&);
RecoPionIdx GetHighestEnergyPionCandidateIndex(const CCPiEvent&);
SignalBackgroundType GetSignalBackgroundType(const CCPiEvent&);

// Helper Fill Histo Functions
namespace ccpi_event {
// Xsec analysis fill functions
void FillSelected(const CCPiEvent&, const std::vector<Variable*>&);
void FillRecoEvent(const CCPiEvent&, const std::vector<Variable*>&);
void FillWSideband(const CCPiEvent&, const std::vector<Variable*>&);
void FillTruthEvent(const CCPiEvent&, const std::vector<Variable*>&);
void FillEfficiencyDenominator(const CCPiEvent&, const std::vector<Variable*>&);
void FillMigration(const CCPiEvent&, const std::vector<Variable*>&,
                   std::string name);

// Study functions
void FillWSideband_Study(CCPiEvent&, std::vector<Variable*>);
void FillCounters(const CCPiEvent&,
                  const std::pair<EventCount*, EventCount*>& counters);
void FillCutVars(CCPiEvent&, const std::vector<Variable*>&);
void FillStackedHists(const CCPiEvent&,
                      const std::vector<Variable*>&);  // all variables
void FillStackedHists(const CCPiEvent&, Variable*,
                      const double fill_value = -999.);  // Single variable
}  // namespace ccpi_event

#endif  // CCPiEvent

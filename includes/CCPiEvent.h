#ifndef CCPiEvent_h
#define CCPiEvent_h

#include "SignalDefinition.h"
#include "CVUniverse.h"
#include "TruthCategories/Sidebands.h" // WSidebandType
#include "TruthCategories/SignalBackground.h" // WSidebandType
#include "Constants.h" // typedef RecoPionIdx, EventCount

#ifndef __CINT__ // Hide class Variable because it has c++11 stuff
#include "Variable.h"
#endif
class Variable;

//==============================================================================
// Container class that knows things about an event:
// * passes cuts
// * is signal
// * event weight
// * is mc, truth
// * vector of candidate pion indices
// * signal definition currently being used
// * whether and what kind of w sideband it is
// And it owns a pointer to the corresponding CVUniverse.
// Has a bunch of helper functions that fill histograms, given a Variable.
//==============================================================================

class CCPiEvent {
  public: 
    CCPiEvent(const bool is_mc, const bool is_truth,
              const SignalDefinition signal_definition, 
              CVUniverse* universe);


    // Fixed by the constructor
    const bool m_is_mc;
    const bool m_is_truth;
    const SignalDefinition m_signal_definition;
    CVUniverse* m_universe;
    std::vector<RecoPionIdx> m_reco_pion_candidate_idxs; // initialized empty, filled by PassesCuts
    bool m_is_signal;
    double m_weight;
    WSidebandType m_w_type;


    // Fixed (directly) outside of constructor -- with time-intensive functions
    bool m_passes_cuts;   // PassesCuts
    bool m_is_w_sideband; // IsWSideband
    RecoPionIdx m_highest_energy_pion_idx; // GetHighestEnergyPionCandidateIndex
};


// Helper Functions
//bool IsWSideband(CCPiEvent&);
bool PassesCuts(CCPiEvent&, bool& is_w_sideband);
RecoPionIdx GetHighestEnergyPionCandidateIndex(const CCPiEvent&);
SignalBackgroundType GetSignalBackgroundType(const CCPiEvent&);
bool PassesCuts(CCPiEvent& e, std::vector<ECuts> cuts = kCutsVector);

// Helper Fill Histo Functions
namespace ccpi_event {
  // Xsec analysis fill functions
  void FillSelected              ( const CCPiEvent&, const std::vector<Variable*>&);
  void FillRecoEvent             ( const CCPiEvent&, const std::vector<Variable*>&);
  void FillWSideband             ( const CCPiEvent&, const std::vector<Variable*>&);
  void FillTruthEvent            ( const CCPiEvent&, const std::vector<Variable*>&);
  void FillEfficiencyDenominator ( const CCPiEvent&, const std::vector<Variable*>&);
  void FillMigration             ( const CCPiEvent&, const std::vector<Variable*>&, std::string name);

  // Study functions
  void FillWSideband_Study   (       CCPiEvent&, std::vector<Variable*>);
  void FillCounters          ( const CCPiEvent&, const std::pair<EventCount*, EventCount*>& counters);
  void FillCutVars           (       CCPiEvent&, const std::vector<Variable*>&);
  void FillStackedHists      ( const CCPiEvent&, const std::vector<Variable*>&);  // all variables
  void FillStackedHists      ( const CCPiEvent&, Variable*, const double fill_value = -999.); // Single variable
}

#endif // CCPiEvent

#ifndef Coherent_H
#define _H

#include <stdlib.h> // exit, EXIT_FAILURE

#include "../CVUniverse.h"
#include "Interaction.h" // GetChannelType
#include "SignalBackground.h"  // GetSignalBackgroundType

enum CoherentType
{ 
  kCOHERENT_S, kCOHERENT_B, kNOTCOHERENT_S, kNOTCOHERENT_B, kNCoherentTypes
};

//==============================================================================

CoherentType GetCoherentType(const CVUniverse& universe,
                             SignalDefinition signal_definition = kOnePi) {
  bool is_signal   = GetSignalBackgroundType(universe, signal_definition) == kS;
  bool is_coherent = GetChannelType(universe) == kCOHPI;
  if (is_signal)
    return is_coherent ? kCOHERENT_S : kNOTCOHERENT_S;
  else
    return is_coherent ? kCOHERENT_B : kNOTCOHERENT_B;
}

std::string GetTruthClassification_LegendLabel(CoherentType coherent_category){
  switch (coherent_category){
    case kCOHERENT_S:
      return "Signal - Coherent";
    case kNOTCOHERENT_S:
      return "Signal - Not Coherent";
    case kCOHERENT_B:
      return "BG - Coherent";
    case kNOTCOHERENT_B:
      return "BG - Not Coherent";
    default: 
      return "ERROR";
  };
}

std::string GetTruthClassification_Name(CoherentType coherent_category){
  switch (coherent_category){
    case kCOHERENT_S:
      return "Signal - Coherent";
    case kNOTCOHERENT_S:
      return "Signal - Not Coherent";
    case kCOHERENT_B:
      return "BG - Coherent";
    case kNOTCOHERENT_B:
      return "BG - Not Coherent";
    default: 
      return "ERROR";
  };
}
#endif

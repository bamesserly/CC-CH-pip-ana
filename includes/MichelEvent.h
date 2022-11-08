#ifndef MichelEvent_h
#define MichelEvent_h

//==============================================================================
// Data container class containing all the michel info needed to do a trackless
// michel analysis.
// The trackless::Michel class initializes these.
//==============================================================================

#include <vector>

namespace trackless {
class Michel;  // forward declare
struct MichelEvent {
  int m_idx = -1; // Index for Best Michel in nmichels
};
} // namespace trackless

#endif

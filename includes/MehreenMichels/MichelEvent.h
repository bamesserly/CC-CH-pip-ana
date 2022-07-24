#ifndef MichelEvent_h
#define MichelEvent_h

#include "CVUniverse.h"
#include "Cluster.h"
#include "Michel.h"

struct MichelEvent {
  int m_idx = -1;                // Index for Best Michel in nmichels
  double m_bestdist = 9999.;     // in mm
  std::vector<double> m_best2D;  // 0: XZ, 1: UZ, 2:VZ
  double m_best_XZ = 9999.;
  double m_best_UZ = 9999.;
  double m_best_VZ = 9999.;
  int m_matchtype;  // 0 NULL 1 UPVTX 2 DOWNVTX 3 UPCLUS 4 DOWNCLUS
  std::vector<Michel> m_nmichels;        // nmatched michels
  std::vector<Michel> m_ntruepiparents;  // michels with true pion parent
  std::vector<Michel>
      m_nmichelspass;  // if some distance cut is applied, we can store the
                       // michels that passed for this event in here
  double best_x = 9999.;
  double best_y = 9999.;
  double best_z = 9999.;
  double b_truex = 9999.;
  double b_truey = 9999.;
  double b_truez = 9999.;
  int bestparentpdg = -1;
  int bestparenttrackid = -1;
  int eventtype =
      0;  // 0 = null, 1 = only 1 pi+ and no other pion, 2= npi+ and other pion,
          // 3 = npi0 and no other pion, 4 = kaons in event, 5 = other
  double lowTpi = 9999.;
};
#endif

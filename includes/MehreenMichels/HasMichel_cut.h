//#include "CVUniverse.h"
//#include "Michel.h"
// Cutter written to select Michels
// PlotUtils includes
#include "PlotUtils/Cut.h"

template <class UNIVERSE, class EVENT>
class hasMichel : public PlotUtils::Cut<UNIVERSE, EVENT> {
 public:
  hasMichel() : PlotUtils::Cut<UNIVERSE, EVENT>("Event Has Michel ") {}

 private:
  using Michel = typename std::remove_reference<decltype(
      std::declval<EVENT>().m_nmichels.front())>::type;

  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    int nmichels = univ.GetNMichels();
    evt.m_bestdist = 9999.;  // setting some default value for best distance
    for (int i = 0; i < nmichels; ++i) {
      Michel current_michel = Michel(univ, i);
      //        if (current_michel.is_fitted != 1) continue;
      //        if (current_michel.vtx_michel_timediff < .400) continue;
      // if (current_michel.pionKE < 80) continue;
      // if (current_michel.true_parentpdg != 211) continue;
      // if (current_michel.true_parentid != 0) continue;
      if (current_michel.true_parentpdg == 211)
        evt.m_ntruepiparents.push_back(current_michel);
      // if (current_michel.overlay_fraction >= 0.5) continue;
      // if (current_michel.true_parentpdg != 211) continue;
      double dist = current_michel.Best3Ddist;  // getting the minimum pion range (vertex to Michel/Clus distance)
      // std::cout << "Pion Energy of Current Michel is " <<
      // current_michel.pionKE << std::endl;
      if (dist <= evt.m_bestdist) {
        evt.m_bestdist = dist;
        evt.m_idx = i;

        evt.m_best_XZ = current_michel.best_XZ;
        evt.m_best_UZ = current_michel.best_UZ;
        evt.m_best_VZ = current_michel.best_VZ;
        evt.m_matchtype = current_michel.BestMatch;
        int bmatch = current_michel.BestMatch;
        if (bmatch == 1 || bmatch == 3) {
          evt.best_x = current_michel.m_x1;
          evt.best_y = current_michel.m_y1;
          evt.best_z = current_michel.m_z1;
        } else if (bmatch == 2 || bmatch == 4) {
          evt.best_x = current_michel.m_x2;
          evt.best_y = current_michel.m_y2;
          evt.best_z = current_michel.m_z2;
        }
        evt.b_truex = current_michel.true_initialx;
        evt.b_truey = current_michel.true_initialy;
        evt.b_truez = current_michel.true_initialz;
      }
      evt.m_nmichels.push_back(current_michel);
    }
    // std::cout << "Now GEtting True Tpi Info from Event" << std::endl;

    double lowtpiinevent = univ.GetTrueTpi();
    evt.lowTpi = lowtpiinevent;
    /*
      std::vector<int> pdgcodes = univ.GetTrueFSPDGCodes();
      int npiplus = 0;
      int npiminus = 0;
      int npi0 = 0;
      int nkaons = 0;
      int nQE = 0;

      for (int j = 0; j < pdgcodes.size(); j++)
      {
         int pdg = pdgcodes[j];
         if (pdg == 211) npiplus++;
         else if (pdg == 111) npi0++;
         else if (pdg == 321 || 311) nkaons++;
      }

      if (npiplus == 1 && npi0 == 0) evt.eventtype = 1;
      else if (npiplus > 0 &&  npi0 > 0) evt.eventtype = 2;
      else if (npiplus == 0 && npi0 > 0) evt.eventtype = 3;
      else if (nkaons > 0) evt.eventtype = 4;
      else {evt.eventtype = 5;}
      //return true;
    */
    return !evt.m_nmichels.empty();
  }
};

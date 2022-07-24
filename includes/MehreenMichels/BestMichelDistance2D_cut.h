// File: BestMichelDistance2D.h
// Brief: Require the the distance to vertex for the best Michel electron is below some value.
// Author: Andrew Olivier aolivier@ur.rochester.edu, Mehreen Sultana msultana@ur.rochester.edu

// PlotUtils includes
#include "PlotUtils/Cut.h"
#include "event/Michel.h"

template <class UNIVERSE, class EVENT>
class BestMichelDistance2D : public PlotUtils::Cut<UNIVERSE, EVENT> {
 public:
  BestMichelDistance2D(const double maxDistance)
      : PlotUtils::Cut<UNIVERSE, EVENT>(
            "Per Michel 2D Distance in at least two views is < " +
            std::to_string(maxDistance) + "mm"),
        m_maxDistance(maxDistance) {}

 private:
  double m_maxDistance;  // Maximum distance from the vertex that the best
                         // Michel can have in mm

  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    // std::cout << "Implementing Michel Cut with 2D distnace of " <<
    // m_maxDistance << " mm" << std::endl;
    std::vector<Michel> nmichelspass;
    for (unsigned int i = 0; i < evt.m_nmichels.size(); i++) {
      // For Vertex Match Check to see if 2D distance cut will
      double upvtxXZ = evt.m_nmichels[i].up_to_vertex_XZ;
      double downvtxXZ = evt.m_nmichels[i].down_to_vertex_XZ;
      double upvtxUZ = evt.m_nmichels[i].up_to_vertex_UZ;
      double downvtxUZ = evt.m_nmichels[i].down_to_vertex_UZ;
      double upvtxVZ = evt.m_nmichels[i].up_to_vertex_VZ;
      double downvtxVZ = evt.m_nmichels[i].down_to_vertex_VZ;

      if (upvtxXZ < m_maxDistance &&
          (upvtxUZ < m_maxDistance || upvtxVZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(0) = 1;
      else if (upvtxUZ < m_maxDistance &&
               (upvtxXZ < m_maxDistance || upvtxVZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(0) = 1;
      else if (upvtxVZ < m_maxDistance &&
               (upvtxXZ < m_maxDistance || upvtxUZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(0) = 1;
      else
        evt.m_nmichels[i].passable_matchtype.at(0) = -1;

      if (downvtxXZ < m_maxDistance &&
          (downvtxUZ < m_maxDistance || downvtxVZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(1) = 2;
      else if (downvtxUZ < m_maxDistance &&
               (downvtxXZ < m_maxDistance || downvtxVZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(1) = 2;
      else if (downvtxVZ < m_maxDistance &&
               (downvtxXZ < m_maxDistance || downvtxUZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(1) = 2;
      else
        evt.m_nmichels[i].passable_matchtype.at(1) = -1;

      double upclusXZ = evt.m_nmichels[i].up_to_clus_XZ;
      double upclusUZ = evt.m_nmichels[i].up_to_clus_UZ;
      double upclusVZ = evt.m_nmichels[i].up_to_clus_VZ;
      double downclusXZ = evt.m_nmichels[i].down_to_clus_XZ;
      double downclusUZ = evt.m_nmichels[i].down_to_clus_UZ;
      double downclusVZ = evt.m_nmichels[i].down_to_clus_VZ;

      if (upclusXZ < m_maxDistance &&
          (upclusUZ < m_maxDistance || upclusVZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(2) = 3;
      else if (upclusUZ < m_maxDistance &&
               (upclusXZ < m_maxDistance || upclusVZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(2) = 3;
      else if (upclusVZ < m_maxDistance &&
               (upclusXZ < m_maxDistance || upclusUZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(2) = 3;
      else
        evt.m_nmichels[i].passable_matchtype.at(2) = -1;

      if (downclusXZ < m_maxDistance &&
          (downclusUZ < m_maxDistance || downclusVZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(3) = 4;
      else if (downclusUZ < m_maxDistance &&
               (downclusXZ < m_maxDistance || downclusVZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(3) = 4;
      else if (downclusVZ < m_maxDistance &&
               (downclusXZ < m_maxDistance || downclusUZ < m_maxDistance))
        evt.m_nmichels[i].passable_matchtype.at(3) = 4;
      else
        evt.m_nmichels[i].passable_matchtype.at(3) = -1;

      std::vector<double> distances3D;

      if (evt.m_nmichels[i].passable_matchtype[0] == 1)
        distances3D.push_back(evt.m_nmichels[i].up_to_vertex_dist3D);
      if (evt.m_nmichels[i].passable_matchtype[1] == 2)
        distances3D.push_back(evt.m_nmichels[i].down_to_vertex_dist3D);
      if (evt.m_nmichels[i].passable_matchtype[2] == 3)
        distances3D.push_back(evt.m_nmichels[i].up_clus_michel_dist3D);
      if (evt.m_nmichels[i].passable_matchtype[3] == 4)
        distances3D.push_back(evt.m_nmichels[i].down_clus_michel_dist3D);
      if (distances3D.empty()) distances3D = {9999., 9999., 9999., 9999.};
      if (distances3D[0] == 9999.)
        continue;  // Comment this line out if you are making efficiency plots
      // std::cout << "Passable Match Types for this Michel are " <<
      // evt.m_nmichels[i].passable_matchtype.at(0) <<
      // evt.m_nmichels[i].passable_matchtype.at(1) <<
      // evt.m_nmichels[i].passable_matchtype.at(2) <<
      // evt.m_nmichels[i].passable_matchtype.at(3) << std::endl;

      std::sort(distances3D.begin(), distances3D.end());

      if (distances3D[0] == evt.m_nmichels[i].up_to_vertex_dist3D) {
        evt.m_nmichels[i].best_XZ = evt.m_nmichels[i].up_to_vertex_XZ;
        evt.m_nmichels[i].best_UZ = evt.m_nmichels[i].up_to_vertex_UZ;
        evt.m_nmichels[i].best_VZ = evt.m_nmichels[i].up_to_vertex_VZ;
        evt.m_nmichels[i].BestMatch = 1;
        evt.m_nmichels[i].Best3Ddist = evt.m_nmichels[i].up_to_vertex_dist3D;

        // std::cout << "This  Michel is UPVTX and has true endpoint " <<
        // evt.m_nmichels[i].trueEndpoint << std::endl;
      } else if (distances3D[0] == evt.m_nmichels[i].down_to_vertex_dist3D) {
        evt.m_nmichels[i].BestMatch = 2;
        evt.m_nmichels[i].best_XZ = evt.m_nmichels[i].down_to_vertex_XZ;
        evt.m_nmichels[i].best_UZ = evt.m_nmichels[i].down_to_vertex_UZ;
        evt.m_nmichels[i].best_VZ = evt.m_nmichels[i].down_to_vertex_VZ;
        evt.m_nmichels[i].Best3Ddist = evt.m_nmichels[i].down_to_vertex_dist3D;
        // std::cout << "This  Michel is DOWNVTX and has true endpoint " <<
        // evt.m_nmichels[i].trueEndpoint << std::endl;

      } else if (distances3D[0] == evt.m_nmichels[i].up_clus_michel_dist3D) {
        evt.m_nmichels[i].BestMatch = 3;
        evt.m_nmichels[i].best_XZ = evt.m_nmichels[i].up_to_clus_XZ;
        evt.m_nmichels[i].best_UZ = evt.m_nmichels[i].up_to_clus_VZ;
        evt.m_nmichels[i].best_VZ = evt.m_nmichels[i].up_to_clus_UZ;
        evt.m_nmichels[i].Best3Ddist = evt.m_nmichels[i].up_clus_michvtx_dist3D;
        // std::cout << "This  Michel is UPCLUS and has true endpoint " <<
        // evt.m_nmichels[i].trueEndpoint << std::endl;

      } else if (distances3D[0] == evt.m_nmichels[i].down_clus_michel_dist3D) {
        evt.m_nmichels[i].BestMatch = 4;
        evt.m_nmichels[i].Best3Ddist =
            evt.m_nmichels[i].down_clus_michvtx_dist3D;
        evt.m_nmichels[i].best_XZ = evt.m_nmichels[i].down_to_clus_XZ;
        evt.m_nmichels[i].best_UZ = evt.m_nmichels[i].down_to_clus_UZ;
        evt.m_nmichels[i].best_VZ = evt.m_nmichels[i].down_to_clus_VZ;
        // std::cout << "This  Michel is DOWNCLUS and has true endpoint " <<
        // evt.m_nmichels[i].trueEndpoint << std::endl;

      } else {
        evt.m_nmichels[i].BestMatch = 0;
        evt.m_nmichels[i].Best3Ddist = 9999.;
        evt.m_nmichels[i].best_XZ = 9999.;
        evt.m_nmichels[i].best_UZ = 9999.;
        evt.m_nmichels[i].best_VZ = 9999.;
        continue;
      }

      nmichelspass.push_back(evt.m_nmichels[i]);
    }

    evt.m_nmichels.clear();         // empty existing vector of Michels
    evt.m_nmichels = nmichelspass;  // replace vector of michels with the vector
                                    // of michels that passed the above cut
    return !evt.m_nmichels.empty();
  }
};

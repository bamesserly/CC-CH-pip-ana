#ifndef CLUSTER_H
#define CLUSTER_H

#include "CVUniverse.h"

// A cluster object will simply be the ith cluster in the event.

// using ClusterMap = std::map<int, Cluster>;

struct Cluster {
  // constructors
  Cluster(const CVUniverse &univ, const int &ci);
  Cluster(){};

  int cluster_idx;

  // Does this michel satisfy our quality
  bool is_quality;

  double energy;

  double time;

  int view;  // 1 = x, 2 = u, 3 = v

  double zpos;  // z position

  double pos;
};

Cluster::Cluster(const CVUniverse &univ, const int &ci) {
  energy = univ.GetVecElem("cluster_energy", ci);  // MeV
  time = univ.GetVecElem("cluster_time", ci) / pow(10, 3);  // microseconds
  pos = univ.GetVecElem("cluster_pos", ci);  // in mm
  zpos = univ.GetVecElem("cluster_z", ci);   // in mm
  view = univ.GetVecElem("cluster_view", ci);  // 1 = X view, 2 = U view, 3 = V view
}

#endif  // CLUSTER_H

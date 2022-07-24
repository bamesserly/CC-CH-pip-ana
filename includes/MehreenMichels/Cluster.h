#ifndef CLUSTER_H
#define CLUSTER_H

// A cluster object will simply be the ith cluster in the event.
class Cluster;

// using ClusterMap = std::map<int, Cluster>;

class Cluster {
 public:
  // constructors
  Cluster(const CVUniverse &univ, int &ci);

  int cluster_idx;
  // Does this michel satisfy our quality
  bool is_quality;

  double energy;

  double time;

  int view;  // 1 = x, 2 = u, 3 = v

  double zpos;  // z position

  double pos;
};

Cluster::Cluster(const CVUniverse &univ, int &ci) {
  energy = univ.GetVecElem("FittedMichel_cluster_energy", ci);  // MeV
  time = univ.GetVecElem("FittedMichel_cluster_time", ci) / pow(10, 3);  // microseconds
  pos = univ.GetVecElem("FittedMichel_cluster_pos", ci);  // in mm
  zpos = univ.GetVecElem("FittedMichel_cluster_z", ci);   // in mm
  view = univ.GetVecElem("FittedMichel_cluster_view", ci);  // 1 = X view, 2 = U view, 3 = V view
}

#endif  // CLUSTER_H

#ifndef CLUSTER_H
#define CLUSTER_H

struct Cluster {
  double energy;  // MeV
  double time;    // microsec
  int view;       // 1 = x, 2 = u, 3 = v
  double zpos;    // mm
  double pos;     // mm
  uint idx;
  Cluster(double energy, double time, int view, double zpos, double pos, uint idx);
  Cluster(){};
  //bool is_quality;
};

Cluster::Cluster(double _energy, double _time, int _view, double _zpos, double _pos, uint _idx)
    : energy(_energy), time(_time), view(_view), zpos(_zpos), pos(_pos), idx(_idx) {}

#endif  // CLUSTER_H

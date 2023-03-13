#ifndef CLUSTER_H
#define CLUSTER_H

struct Cluster {
  Cluster(double energy, double time, double zpos, int view, double pos);
  Cluster(){};
  double energy;  // MeV
  double time;    // microsec
  int view;       // 1 = x, 2 = u, 3 = v
  double zpos;    // mm
  double pos;     // mm
  int cluster_idx;
  bool is_quality;
};

Cluster::Cluster(double _energy, double _time, double _zpos, int _view,
                 double _pos)
    : energy(_energy), time(_time), view(_view), zpos(_zpos), pos(_pos) {}

#endif  // CLUSTER_H

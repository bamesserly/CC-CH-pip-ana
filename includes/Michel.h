#ifndef Michel_h
#define Michel_h

//==============================================================================
// Two classes and helper functions for different michel electron objects
// endpoint: Aaron-style michels matched to track endpoints (written by me)
// trackless: Mehreen-style michels NOT matched to track endpoints (written by
// mehreen)
//==============================================================================

#include "CVUniverse.h"

class Cluster;
class MichelEvent;

//==============================================================================
// Track Endpoint Michel
//==============================================================================
namespace endpoint {
class Michel;  // forward declare
typedef std::map<int, Michel> MichelMap;

class Michel {
 public:
  enum EMatchCategory { kNoMatch, kFit, kNoFit, kOV, kNMatchCategories };

  // constructors
  Michel(const CVUniverse& univ, const int i, const int v);
  Michel()
      : idx(-105),
        vtx(-106),
        had_idx(-107),
        match_category(kNoMatch),
        fit_distance(-1.) {}

  // integer that uniquely identifies cluster of hits
  int idx;

  // integer corresponding to vertex to which cluster was matched.
  // vtx == 0 --> interaction vertex
  // vtx == 1 --> "first" track endpoint, corresponds to hadron index 0.
  //              i.e. hadron index = michel vtx - 1.
  int vtx;

  // hadron index to which this michel is matched (had_idx = vtx - 1)
  int had_idx;

  EMatchCategory match_category;
  double fit_distance;
};

// Michel quality cuts
bool IsQualityMatchedMichel_Fit(const double fit_dist, const double fit_cut);
bool IsQualityMatchedMichel_NoFit(const double nofit_dist,
                                  const double nofit_cut);
bool IsQualityMatchedMichel_OneView(const double ov_dist, const double ov_cut);

// -- Given a single michel cluster matched to two vertices
//    return vertex with the better-matched michel.
Michel CompareMichels(const Michel& r, const Michel& c);

// Add michel to MichelMap. Check if this cluster has already been matched.
// Then use only the best match.
bool AddOrReplaceMichel(MichelMap& mm, const Michel& m);

//  Collect the good michels in this event
MichelMap GetQualityMichels(const CVUniverse& univ);
}  // namespace endpoint

//==============================================================================
// Trackless Michel
//==============================================================================
namespace trackless {
class Michel {
 public:
  // constructors
  Michel(const CVUniverse& univ, const int ci);
  Michel(){};

  // fill in more complicated stuff in "generic info"
  void DoMoreInitializationAndProcessing();

  // fill in "best matching stuff"
  void DoMatching();

  // Gets info for Vtx Match
  void DoesMichelMatchVtx(const CVUniverse& univ);

  // Gets info for ClusterMatch
  void DoesMichelMatchClus(const CVUniverse& univ);

  // get type for best match out of all four saved matches
  void GetBestMatch();

  // Function to calculate pion angle with respect to beam. Simply gets
  // angle between best michel endpoint and vertex.
  void GetPionAngle(const CVUniverse& univ);

  // std::vector<Michel*> CreateMichels(CVUniverse& univ);

  // Data
  std::vector<double> up_location;    // upstream location 0 X 1 U 2 V 3 Z
  std::vector<double> down_location;  // downstream location
  double m_x1 = 9999.;                // Michel Endpoint 1 x
  double m_x2 = 9999.;                // Michel Endpoint 2 x
  double m_y1 = 9999.;                // Michel Endpoint 1 y
  double m_y2 = 9999.;                // Michel Endpoint 2 y
  double m_u1 = 9999.;                // Michel Endpoint 1 u
  double m_u2 = 9999.;                // Michel Endpoint 2 u
  double m_v1 = 9999.;                // Michel Endpoint 1 v
  double m_v2 = 9999.;                // Michel Endpoint 2 v
  double m_z1 = 9999.;                // Mihel Endpoint z 1
  double m_z2 = 9999.;                // Michel Endpoiint z2
  double energy = -999.;              // Michel energy
  double time = -999.;                // Michel Time
  int is_fitted = -1;                 // Is the Michel fitted? 0 no. 1 yes.
  // Following are 2D distances (were stored in vectors but now as explicit data
  // members)
  double up_to_vertex_XZ = 9999.;
  double up_to_vertex_UZ = 9999.;
  double up_to_vertex_VZ = 9999.;
  double down_to_vertex_XZ = 9999.;
  double down_to_vertex_UZ = 9999.;
  double down_to_vertex_VZ = 9999.;
  double down_to_clus_XZ = 9999.;
  double down_to_clus_UZ = 9999.;
  double down_to_clus_VZ = 9999.;
  double up_to_clus_XZ = 9999.;
  double up_to_clus_UZ = 9999.;
  double up_to_clus_VZ = 9999.;
  // Michel End point to Vertex distance  (up = endpoint 1 and down = endpoint
  // 2) TODO: actually find out which end point is upstream or downstream
  double up_to_vertex_dist3D = 9999.;
  double down_to_vertex_dist3D = 9999.;
  // Maybe keep a vector of clusters that matched to each endpoint?
  std::vector<Cluster*> cluster_to_up_match;
  std::vector<Cluster*> cluster_to_down_match;
  // 3D distances between cluster and Michel
  double up_to_cluster_dist3D = 9999.;  // Distance between vertex and cluster
                                        // that was matched to endpoint 1
  double down_to_cluster_dist3D =
      9999.;  // Distance between vertex and cluster matched to endpoint 2
  double up_clus_michel_dist3D =
      9999.;  // Distance between Michel endpoint 1 and clusters
  double down_clus_michel_dist3D =
      9999.;  // Distance between Michel endpoint 2 and clusters
  double up_clus_michvtx_dist3D =
      9999.;  // Distance between the Michel end point 1 that matched to
              // clusters and the vertex - this will be used as pion range
  double down_clus_michvtx_dist3D =
      9999.;  // Distance between the Michel endpoint 2 that matched to clusters
              // and the vertex - this will be used as pion range
  double vtx_michel_timediff = 9999.;

  double overlay_fraction =
      -1.0;  // Overlay fraction of the Michel, Default if Data. 0 if MC 1 if
             // Data... (maybe some events in between?)
  int nclusters = 0;      // number of (non-muon) clusters in the primary event
  int vtx_endpoint = 0;   // 1 or 2 for which Michel end point is closest
  int clus_endpoint = 0;  // 1 or 2 for which Michel endpoint is closest
  // best matching stuff
  // enum *best_cluster_match; // just tells you which of the four matches are
  // the best match enum BestClusterMatch {kUpVtxMatch, kDownVtxMatch,
  // kUpClusMatch, kDownClusMatch, kNClusterMatches}; the following is in place
  // until i figure out how to use the enum.
  int BestMatch = 0;  // 0 = null match, 1= kUpVtxMatch, 2 = kDownVtxMatch, 3 =
                      // kUpClusMatch, 4=  kDownClusMatch,
  int SecondBestMatch = 0;
  int tuple_idx;  // index of the Michel out of all the michels saved in the
                  // tuple
  double Best3Ddist =
      9999.;  // Best 3D distance out of eitehr a vertex or a cluster match
  // best 2D distance for the best type of match
  double best_XZ = 9999.;
  double best_UZ = 9999.;
  double best_VZ = 9999.;
  // Want to save the index of the clusters that the Michel best matched to.
  int xclus_idx;
  int uclus_idx;
  int vclus_idx;

  // True initial position of the michel  TODO: initial the other Michel truth
  // member data here  (energy, time, momentum etc)
  double true_angle = 9999.;
  double true_initialx = 9999.;
  double true_initialy = 9999.;
  double true_initialz = 9999.;
  double true_e = 9999.;
  double true_p = 9999.;
  double true_pdg = -1.0;
  int true_parentid = -1;
  int true_parentpdg = -1;
  double true_parent_energy = -9999.;
  double true_parent_p = -9999.;
  double true_parent_px = -9999.;
  double true_parent_py = -9999.;
  double true_parent_pz = -9999.;
  double true_parent_xi = -9999.;
  double true_parebt_yi = -9999.;
  double true_parent_zi = -9999.;
  double true_parent_xf = -9999.;
  double true_parent_yf = -9999.;
  double true_parent_zf = -9999.;
  // the following member data were created to investigate my weird convoluted
  // way of geting x and y values. Probably dont need them now.
  // TODO: check and remove the following member data
  double best_angle = -9999.;
  double up_clus_x = 9999.;
  double up_clus_y = 9999.;
  double up_clus_z = 9999.;
  double down_clus_x = 9999.;
  double down_clus_y = 9999.;
  double down_clus_z = 9999.;
  double up_vtx_x = 9999.;
  double up_vtx_y = 9999.;
  double up_vtx_z = 9999.;
  double down_vtx_x = 9999.;
  double down_vtx_y = 9999.;
  double down_vtx_z = 9999.;
  double is_overlay = -1;
  double pionKE = -9999.;
  // Adding the following member data to determine the true endpoint position of
  // the Michel
  // 0 = Overlay Michel, 1 = Endpoint 1 is correct intial position of
  // Michel, 2 = Endpoint 2 is correct Initial Position of Michel
  int trueEndpoint = -1;

  // ~DeleteMichel(){delete this;};  // This is going to be the main destructor.

  // -1 is default/NULL like above. 1 = Endpoint 1 is better match, 2 =
  // Endpoint 2 is better match
  int recoEndpoint = -1;

  // This vector will contain a value for each match type -1 or the
  // matchtype 1, 2, 3, 4 for UpVtx, DownVTx, Upclus, DownClus depending of
  // that match type passes our 2D distance cut. (if distance is large, then
  // it'll pass all of them).
  std::vector<int> passable_matchtype{-1, -1, -1, -1};
};

// Create Michel objects for each Michel candidate. Add the good ones to the
// MichelEvent container.
MichelEvent GetQualityMichels(const CVUniverse& univ);
}  // namespace trackless

#endif

#ifndef Michel_H
#define Michel_H

#include <ctime>

#include "CVUniverse.h"
#include "Cluster.h"

namespace endpoint {
bool IsQualityMatchedMichel_Fit(double fit_dist, double fit_cut);
bool IsQualityMatchedMichel_NoFit(double nofit_dist, double nofit_cut);
bool IsQualityMatchedMichel_OneView(double ov_dist, double ov_cut);

class Michel {
 public:
  enum EMatchCategory { kNoMatch, kFit, kNoFit, kOV, kNMatchCategories };
  // constructors
  Michel(const CVUniverse& univ, int i, int v);
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

  // All fit strategies (fitted, 2/3-view no-fit, and one-view no-fit) are
  // calculated and saved in these branches if possible.
  // Get the distance given a strategy and vertex.
  double GetDistMichel(const CVUniverse& univ,
                       const EMatchCategory match_strategy,
                       const unsigned int vertex) const;

  EMatchCategory match_category;
  double fit_distance;
};

typedef std::map<int, Michel> MichelMap;

Michel::Michel(const CVUniverse& univ, int i, int v)
    : idx(i),
      vtx(v),
      had_idx(v - 1),
      match_category(kNoMatch),
      fit_distance(-1.) {
  bool isIntVtx = (vtx == 0);

  // distances for fitted, 2/3-view nofit, and 1-view michels
  double mm_fit_dist = GetDistMichel(univ, kFit, vtx);
  double mm_nofit_dist = GetDistMichel(univ, kNoFit, vtx);
  double mm_ov_dist = GetDistMichel(univ, kOV, vtx);

  // NEW
  const double FIT_CUT = isIntVtx ? 9.0 : 7.5;     // cm
  const double NOFIT_CUT = isIntVtx ? 10.0 : 50.;  // cm
  // OLD
  // const double FIT_CUT   = isIntVtx ? 9.0  :  5.; // cm
  // const double NOFIT_CUT = isIntVtx ? 10.0 : 10.; // cm
  // TODO distinguish between 2/3 view and OV
  // const double FIT_CUT   = isIntVtx ? 9. :  15.0; // cm
  // const double NOFIT_CUT = isIntVtx ? 10. : 15.0; // cm
  if (IsQualityMatchedMichel_Fit(mm_fit_dist, FIT_CUT)) {
    match_category = kFit;
    fit_distance = mm_fit_dist;
  } else if (IsQualityMatchedMichel_NoFit(mm_nofit_dist, NOFIT_CUT)) {
    match_category = kNoFit;
    fit_distance = mm_nofit_dist;
  } else if (IsQualityMatchedMichel_OneView(mm_ov_dist, NOFIT_CUT)) {
    match_category = kOV;
    fit_distance = mm_ov_dist;
  } else {
    match_category = kNoMatch;
    fit_distance = -2.;
  }
}

double Michel::GetDistMichel(const CVUniverse& univ,
                             const EMatchCategory match_strategy,
                             const unsigned int vtx) const {
  std::string branch_name;
  switch (match_strategy) {
    case kFit:
      branch_name = "matched_michel_end_dist";
      break;
    case kNoFit:
      branch_name = "matched_michel_avg_dist";
      break;
    case kOV:
      branch_name = "matched_michel_ov_dist";
      break;
    default:
      return -1.;
  }
  double match_dist = univ.GetVecElem(branch_name.c_str(), vtx);  // mm
  match_dist = match_dist / 10.;                                  // cm

  // IF bogus match distance then throw an error. But, after a bugfix, I no
  // longer have any reason to suspect this will ever occur. Anyway: distances
  // greater than epsilon and less than 5 m.
  bool is_valid_distance = !isnan(match_dist) &&
                           (match_dist == 0. || fabs(match_dist) > 0.0001) &&
                           fabs(match_dist) < 500;
  // assert(is_valid_distance);
  if (!is_valid_distance) {
    std::cerr
        << "WARNING endpoint::michel match distance branch access failure\n";
    match_dist = -1.;  // negative will fail quality cuts.
  }

  return match_dist;
}

//==============================================================================
// MICHEL QUALITY CUTS
//==============================================================================
bool IsQualityMatchedMichel_Fit(const double fit_dist, const double fit_cut) {
  return fit_dist > 0 && fit_dist < fit_cut;
}

bool IsQualityMatchedMichel_NoFit(const double nofit_dist,
                                  const double nofit_cut) {
  return nofit_dist > 0 && nofit_dist < nofit_cut;
}

bool IsQualityMatchedMichel_OneView(const double ov_dist, const double ov_cut) {
  return ov_dist > 0 && ov_dist * (2. / 3.) < ov_cut;
}

// -- Given a single michel cluster matched to two vertices
//    return vertex with the better-matched michel.
Michel CompareMichels(Michel r, Michel c) {
  if (r.match_category > c.match_category)
    return r;
  else if (r.match_category < c.match_category)
    return c;
  else {
    if (r.fit_distance < c.fit_distance)
      return r;
    else if (r.fit_distance > c.fit_distance)
      return c;
    else {
      // This must mean we're comparing the same michels with the same fits.
      // std::cout << "WEIRD COMPAREMICHELS PROBLEM" << std::endl;
      return r;
    }
  }
}

// Add michel to MichelMap. Check if this cluster has already been matched.
// Then use only the best match.
bool AddOrReplaceMichel(MichelMap& mm, Michel m) {
  std::pair<MichelMap::iterator, bool> SC;
  if (mm.count(m.idx) == 0)
    SC = mm.insert(pair<int, Michel>(m.idx, m));
  else {
    Michel reigning_michel = mm.at(m.idx);
    Michel best_michel = CompareMichels(reigning_michel, m);
    mm[m.idx] = best_michel;
  }
  return true;
}

//============================================================================
//  Collect the good michels in this event
//============================================================================
// Get quality michels map<(int michel cluster ID, Michel)>
// * A Michel is uniquely ID-ed by its cluster(of hits) integer index.
//   * The michel tool has already vetted the clusters themselves.
//     Whereas here we evaluate the quality of the fit.
// * At most one interaction vertex michel (which MUST be "fitted")
// * An OV michel will only be kept if no better michels are in the event.
// * In the case of 2+ OV michels, only keep the best one.
// * required: one-to-one michel<->vertex matching:
//    * when one michel cluster is matched to multiple vertices, the vtx
//    with the better match is chosen.
//    * We don't have the problem in the other direction -- michel tool: a
//    vertex cannot be matched to more than one cluster.
// * At the moment, we can return a single michel that is a quality
//   interaction vertex michel. We'd like to call these signal, but we don't
//   know the pion energy...yet. In the meantime, cut them.
MichelMap GetQualityMichels(const CVUniverse& univ) {
  std::map<int, Michel> ret_michels;
  std::vector<int> matched_michel_idxs = univ.GetVec<int>("matched_michel_idx");

  // Loop vertices in the event, i.e. the indices of the michel index vector
  for (uint vtx = 0; vtx < matched_michel_idxs.size(); ++vtx) {
    int mm_idx = matched_michel_idxs[vtx];

    // NO MATCH -- GO TO NEXT VTX. No michel cluster matched to this vertex.
    if (mm_idx < 0) continue;

    // MICHEL CONSTRUCTOR
    // Set match category (e.g. fitted, one-view, etc), match distance.
    Michel current_michel = Michel(univ, mm_idx, vtx);

    // MICHEL MATCH QUALITY CUT -- NEXT VTX
    // A michel cluster was matched to this vertex, but the match doesn't
    // pass match quality cuts.
    if (current_michel.match_category == Michel::kNoMatch) continue;

    // SKIP VERTEX MICHEL -- this is trackless:: namespace territory now. I wash
    // my hands. Skip michels matched this way.
    //
    // For the record, the previous method was to only consider it if the fit
    // category was better than kNoFit.
    if (vtx == 0) {
      continue;
    }

    assert(current_michel.had_idx >= 0 &&
           "endpoint::GetQualityMichels found a vertex michel");

    // ENDPOINT MICHELS

    // ZERO MICHELS FOUND SO FAR
    if (ret_michels.size() == 0)
      ;

    // ONE MICHEL FOUND SO FAR
    // -- If either the michel we have already or this michel is OV, pick
    //    the better of the two. Only one will remain.
    else if (ret_michels.size() == 1) {
      Michel reigning_michel = (ret_michels.begin())->second;
      if (reigning_michel.match_category == Michel::kOV ||
          current_michel.match_category == Michel::kOV)
        ret_michels.clear();
      current_michel = CompareMichels(reigning_michel, current_michel);
    }
    // 2+ MICHELS FOUND SO FAR
    else {
      if (current_michel.match_category == Michel::kOV) continue;
    }

    // ADD THIS MICHEL
    // When a cluster is matched to two vertices, pick the vtx with the
    // better match.
    bool SC = AddOrReplaceMichel(ret_michels, current_michel);
  }  // end vtx loop
  return ret_michels;
}

}  // namespace endpoint

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

// Data container class containing all the michel info needed to do a trackless
// michel analysis.
// The trackless::Michel class initializes these.
struct MichelEvent {
  //~MichelEvent() {
  //  for (auto m : m_nmichelspass) delete m;
  //  for (auto m : m_nmichels) delete m;
  //  for (auto m : m_ntruepiparents) delete m;
  //};
  int m_idx = -1; // Index for Best Michel in nmichels
  double m_bestdist = 9999.;     // in mm
  std::vector<double> m_best2D;  // 0: XZ, 1: UZ, 2:VZ
  double m_best_XZ = 9999.;
  double m_best_UZ = 9999.;
  double m_best_VZ = 9999.;
  int m_matchtype;  // 0 NULL 1 UPVTX 2 DOWNVTX 3 UPCLUS 4 DOWNCLUS
  std::vector<Michel*> m_nmichels;        // nmatched michels
  std::vector<Michel*> m_ntruepiparents;  // michels with true pion parent

  // if some distance cut is applied, we can store the michels that passed for
  // this event in here.
  std::vector<Michel*> m_nmichelspass;

  double best_x = 9999.;
  double best_y = 9999.;
  double best_z = 9999.;
  double b_truex = 9999.;
  double b_truey = 9999.;
  double b_truez = 9999.;
  int bestparentpdg = -1;
  int bestparenttrackid = -1;

  // 0 = null, 1 = only 1 pi+ and no other pion, 2= npi+ and other pion,
  // 3 = npi0 and no other pion, 4 = kaons in event, 5 = other
  int eventtype = 0;

  double lowTpi = 9999.;
};

// Create Michel objects for each Michel candidate. Add the passing ones to
// the MichelEvent container.
// MichelEvent GetQualityMichels(const CVUniverse& u) { return MichelEvent(); }
MichelEvent GetQualityMichels(const CVUniverse& univ) {
  MichelEvent evt{};
  //==========================================================================
  // First: Create a Michel object from each candidate and add them to a
  // container
  //==========================================================================
  std::vector<Michel*> temp_michels;
  for (int i = 0; i < univ.GetNMichels(); ++i) {
    Michel* current_michel = new Michel(univ, i);
    if (current_michel->true_parentpdg == 211)
      evt.m_ntruepiparents.push_back(current_michel);
    double dist =
        current_michel->Best3Ddist;  // getting the minimum pion range (vertex
                                     // to Michel/Clus distance)
    if (dist <= evt.m_bestdist) {
      evt.m_bestdist = dist;
      evt.m_idx = i;
      evt.m_best_XZ = current_michel->best_XZ;
      evt.m_best_UZ = current_michel->best_UZ;
      evt.m_best_VZ = current_michel->best_VZ;
      evt.m_matchtype = current_michel->BestMatch;
      int bmatch = current_michel->BestMatch;
      if (bmatch == 1 || bmatch == 3) {
        evt.best_x = current_michel->m_x1;
        evt.best_y = current_michel->m_y1;
        evt.best_z = current_michel->m_z1;
      } else if (bmatch == 2 || bmatch == 4) {
        evt.best_x = current_michel->m_x2;
        evt.best_y = current_michel->m_y2;
        evt.best_z = current_michel->m_z2;
      }
      evt.b_truex = current_michel->true_initialx;
      evt.b_truey = current_michel->true_initialy;
      evt.b_truez = current_michel->true_initialz;
    }
    temp_michels.push_back(current_michel);
  }

  double lowtpiinevent = univ.GetTrueTpi();
  evt.lowTpi = lowtpiinevent;
  //==========================================================================

  //==========================================================================
  // Second: loop the container, identify passing michels, and set them to the
  // MichelEvent
  //==========================================================================
  const double m_maxDistance = 150;  // Maximum distance from the vertex that
                                     // the best Michel can have in mm
  // loop over Michel* objects
  for (auto michel : temp_michels) {
    // For Vertex Match Check to see if 2D distance cut will
    double upvtxXZ = michel->up_to_vertex_XZ;
    double downvtxXZ = michel->down_to_vertex_XZ;
    double upvtxUZ = michel->up_to_vertex_UZ;
    double downvtxUZ = michel->down_to_vertex_UZ;
    double upvtxVZ = michel->up_to_vertex_VZ;
    double downvtxVZ = michel->down_to_vertex_VZ;

    // Problem happens in here
    if (upvtxXZ < m_maxDistance &&
        (upvtxUZ < m_maxDistance || upvtxVZ < m_maxDistance)) {
      michel->passable_matchtype.at(0) = 1;
    } else if (upvtxUZ < m_maxDistance &&
               (upvtxXZ < m_maxDistance || upvtxVZ < m_maxDistance)) {
      michel->passable_matchtype.at(0) = 1;
    } else if (upvtxVZ < m_maxDistance &&
               (upvtxXZ < m_maxDistance || upvtxUZ < m_maxDistance)) {
      michel->passable_matchtype.at(0) = 1;
    } else {
      michel->passable_matchtype.at(0) = -1;
    }

    if (downvtxXZ < m_maxDistance &&
        (downvtxUZ < m_maxDistance || downvtxVZ < m_maxDistance))
      michel->passable_matchtype.at(1) = 2;
    else if (downvtxUZ < m_maxDistance &&
             (downvtxXZ < m_maxDistance || downvtxVZ < m_maxDistance))
      michel->passable_matchtype.at(1) = 2;
    else if (downvtxVZ < m_maxDistance &&
             (downvtxXZ < m_maxDistance || downvtxUZ < m_maxDistance))
      michel->passable_matchtype.at(1) = 2;
    else
      michel->passable_matchtype.at(1) = -1;

    double upclusXZ = michel->up_to_clus_XZ;
    double upclusUZ = michel->up_to_clus_UZ;
    double upclusVZ = michel->up_to_clus_VZ;
    double downclusXZ = michel->down_to_clus_XZ;
    double downclusUZ = michel->down_to_clus_UZ;
    double downclusVZ = michel->down_to_clus_VZ;

    if (upclusXZ < m_maxDistance &&
        (upclusUZ < m_maxDistance || upclusVZ < m_maxDistance))
      michel->passable_matchtype.at(2) = 3;
    else if (upclusUZ < m_maxDistance &&
             (upclusXZ < m_maxDistance || upclusVZ < m_maxDistance))
      michel->passable_matchtype.at(2) = 3;
    else if (upclusVZ < m_maxDistance &&
             (upclusXZ < m_maxDistance || upclusUZ < m_maxDistance))
      michel->passable_matchtype.at(2) = 3;
    else
      michel->passable_matchtype.at(2) = -1;

    if (downclusXZ < m_maxDistance &&
        (downclusUZ < m_maxDistance || downclusVZ < m_maxDistance))
      michel->passable_matchtype.at(3) = 4;
    else if (downclusUZ < m_maxDistance &&
             (downclusXZ < m_maxDistance || downclusVZ < m_maxDistance))
      michel->passable_matchtype.at(3) = 4;
    else if (downclusVZ < m_maxDistance &&
             (downclusXZ < m_maxDistance || downclusUZ < m_maxDistance))
      michel->passable_matchtype.at(3) = 4;
    else
      michel->passable_matchtype.at(3) = -1;

    std::vector<double> distances3D;

    if (michel->passable_matchtype[0] == 1)
      distances3D.push_back(michel->up_to_vertex_dist3D);
    if (michel->passable_matchtype[1] == 2)
      distances3D.push_back(michel->down_to_vertex_dist3D);
    if (michel->passable_matchtype[2] == 3)
      distances3D.push_back(michel->up_clus_michel_dist3D);
    if (michel->passable_matchtype[3] == 4)
      distances3D.push_back(michel->down_clus_michel_dist3D);
    if (distances3D.empty()) distances3D = {9999., 9999., 9999., 9999.};
    if (distances3D[0] == 9999.)
      continue;  // Comment this line out if you are making efficiency plots

    std::sort(distances3D.begin(), distances3D.end());

    if (distances3D[0] == michel->up_to_vertex_dist3D) {
      michel->best_XZ = michel->up_to_vertex_XZ;
      michel->best_UZ = michel->up_to_vertex_UZ;
      michel->best_VZ = michel->up_to_vertex_VZ;
      michel->BestMatch = 1;
      michel->Best3Ddist = michel->up_to_vertex_dist3D;

    } else if (distances3D[0] == michel->down_to_vertex_dist3D) {
      michel->BestMatch = 2;
      michel->best_XZ = michel->down_to_vertex_XZ;
      michel->best_UZ = michel->down_to_vertex_UZ;
      michel->best_VZ = michel->down_to_vertex_VZ;
      michel->Best3Ddist = michel->down_to_vertex_dist3D;
      // std::cout << "This  Michel is DOWNVTX and has true endpoint " <<
      // michel->trueEndpoint << std::endl;

    } else if (distances3D[0] == michel->up_clus_michel_dist3D) {
      michel->BestMatch = 3;
      michel->best_XZ = michel->up_to_clus_XZ;
      michel->best_UZ = michel->up_to_clus_VZ;
      michel->best_VZ = michel->up_to_clus_UZ;
      michel->Best3Ddist = michel->up_clus_michvtx_dist3D;
      // std::cout << "This  Michel is UPCLUS and has true endpoint " <<
      // michel->trueEndpoint << std::endl;

    } else if (distances3D[0] == michel->down_clus_michel_dist3D) {
      michel->BestMatch = 4;
      michel->Best3Ddist = michel->down_clus_michvtx_dist3D;
      michel->best_XZ = michel->down_to_clus_XZ;
      michel->best_UZ = michel->down_to_clus_UZ;
      michel->best_VZ = michel->down_to_clus_VZ;
      // std::cout << "This  Michel is DOWNCLUS and has true endpoint " <<
      // michel->trueEndpoint << std::endl;

    } else {
      michel->BestMatch = 0;
      michel->Best3Ddist = 9999.;
      michel->best_XZ = 9999.;
      michel->best_UZ = 9999.;
      michel->best_VZ = 9999.;
      continue;
    }

    evt.m_nmichels.push_back(michel);
  }
  //==========================================================================

  return evt;
}
}  // namespace trackless

#endif

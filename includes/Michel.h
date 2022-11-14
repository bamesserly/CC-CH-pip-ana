#ifndef Michel_H
#define Michel_H

#include "CVUniverse.h"
#include "MichelEvent.h"  // trackless::MichelEvent

namespace endpoint {
class Michel;

typedef std::map<int, Michel> MichelMap;

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
  assert(is_valid_distance);

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

    // VERTEX MICHEL -- trackless:: namespace territory now. I wash my hands.
    bool isIntVtx = (vtx == 0);
    if (isIntVtx 
        //&& current_michel.match_category <= Michel::kNoFit
    ) {
      continue;
    }
    // ENDPOINT MICHELS
    else {
      assert(had_idx > 0 && "endpoint::GetQualityMichels found a vertex michel");
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
// Create Michel objects for each Michel candidate. Add the good ones to the
// MichelEvent container.
MichelEvent GetQualityMichels(const CVUniverse& univ) { return MichelEvent(); }
}  // namespace trackless

#endif

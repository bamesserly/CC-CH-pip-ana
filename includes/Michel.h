#ifndef Michel_H
#define Michel_H

#include "CVUniverse.h"

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
  // Get the distance given a stragegy and cluster index.
  double GetDistMichel(const CVUniverse& univ,
                       const EMatchCategory match_strategy,
                       const unsigned int cluster) const;

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
  double mm_fit_dist = GetDistMichel(univ, kFit, idx);
  double mm_nofit_dist = GetDistMichel(univ, kNoFit, idx);
  double mm_ov_dist = GetDistMichel(univ, kOV, idx);

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
                             const unsigned int cluster) const {
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
  double match_dist = univ.GetVecElem(branch_name.c_str(), cluster);  // mm
  match_dist = match_dist / 10.;                                      // cm

  // In MAD p2, elements of these branches are not all necessarily
  // initialized for a given michel cluster, leading to undefined behavior.
  // Must fix this in p3, but in the meantime catch most cases by requiring
  // distances to be greater than epsilon and less than 5 m.
  bool is_valid_distance = !isnan(match_dist) &&
                           (match_dist == 0. || fabs(match_dist) > 0.0001) &&
                           fabs(match_dist) < 500;

  // if (!is_valid_distance) {
  //  for (auto i : univ.GetVec<int>(branch_name.c_str()))
  //    std::cout << i << " ";
  //  std::cout << "\n";
  //}

  match_dist = match_dist ? is_valid_distance : -1.;

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

#endif

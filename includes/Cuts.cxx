#ifndef Cuts_cxx
#define Cuts_cxx

#include "Cuts.h"

#include "TruthCategories/Sidebands.h"  // sidebands::kSidebandCutVal
#include "utilities.h"                  // ContainerEraser

//==============================================================================
// Helpers
//==============================================================================
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

//==============================================================================
// Passes ALL Cuts
//==============================================================================
// cut-by-cut, we fill endpoint_michels and vertex_michels.
// if a track fails a cut, we remove the track's michel from the lists.
// then at the end, return the track indices.
bool PassesCuts(CVUniverse& univ, std::vector<int>& pion_candidate_idxs,
                bool is_mc, SignalDefinition signal_definition,
                std::vector<ECuts> cuts) {
  pion_candidate_idxs.clear();
  static MichelMap endpoint_michels;
  static MichelMap
      vertex_michels;  // Keep track of these, but not used currently
  endpoint_michels.clear();
  vertex_michels.clear();
  bool pass = true;
  for (auto c : cuts) {
    univ.SetPionCandidates(GetHadIdxsFromMichels(
        endpoint_michels));  // Set the pion candidates to the universe
    pass = pass && PassesCut(univ, c, is_mc, signal_definition,
                             endpoint_michels, vertex_michels);
  }

  // Each endpoint michel has an associated hadron track.
  // Our official pion candidates are those tracks.
  pion_candidate_idxs = GetHadIdxsFromMichels(endpoint_michels);

  return pass;
}

// Check all cuts AND whether we are W sideband. Save a lot of time.
// Strategy is:
// 1. check all cuts except for W
// 2. check if W > 1.5 (sideband)
// 3. check if W < 1.4 (signal)
bool PassesCuts(CVUniverse& universe, std::vector<int>& pion_candidate_idxs,
                const bool is_mc, SignalDefinition signal_definition,
                bool& is_w_sideband, std::vector<ECuts> cuts) {
  // is the W cut even in the cuts vector provided?
  bool do_w_cut = std::find(cuts.begin(), cuts.end(), kWexp) != cuts.end();

  // either way, attempt to remove it
  std::vector<ECuts> w_sideband_cuts = kCutsVector;
  w_sideband_cuts.erase(
      std::remove(w_sideband_cuts.begin(), w_sideband_cuts.end(), kWexp),
      w_sideband_cuts.end());

  // check passes all but w cut
  bool passes_all_but_w_cut = PassesCuts(universe, pion_candidate_idxs, is_mc,
                                         signal_definition, w_sideband_cuts);

  // is w sideband = all cuts but W && W > 1.5
  is_w_sideband = passes_all_but_w_cut &&
                  (universe.GetWexp() >= sidebands::kSidebandCutVal);

  // Finally check all cuts == all cuts but W && W
  bool passes_all_cuts = passes_all_but_w_cut;
  if (do_w_cut)
    passes_all_cuts =
        passes_all_but_w_cut && WexpCut(universe, signal_definition);

  return passes_all_cuts;
}

// NEW! Return passes_all_cuts, is_w_sideband, and pion_candidate_indices
std::tuple<bool, bool, std::vector<int>> PassesCuts(
    CVUniverse& universe, const bool is_mc,
    const SignalDefinition signal_definition, std::vector<ECuts> cuts) {
  //============================================================================
  // passes all cuts but w cut
  //============================================================================
  MichelMap endpoint_michels;
  MichelMap vertex_michels;
  bool passes_all_but_w_cut = true;
  for (auto c : GetWSidebandCuts()) {
    // Set the pion candidates to the universe. The values set in early cuts
    // are used for later cuts, which is why we assign them to the CVU.
    universe.SetPionCandidates(GetHadIdxsFromMichels(endpoint_michels));

    passes_all_but_w_cut =
        passes_all_but_w_cut && PassesCut(universe, c, is_mc, signal_definition,
                                          endpoint_michels, vertex_michels);
  }

  //============================================================================
  // The cuts function returns a container of endpoint michels which are
  // matched to hadron tracks that have passed the pion candidate cuts. From
  // here on out, use the pion candidates to calculate pion quantities for this
  // event-universe.
  //============================================================================
  std::vector<int> pion_candidate_idxs =
      GetHadIdxsFromMichels(endpoint_michels);

  //============================================================================
  // is in the w sideband
  //============================================================================
  bool is_w_sideband = passes_all_but_w_cut &&
                       (universe.GetWexp() >= sidebands::kSidebandCutVal);

  //============================================================================
  // finally: check the w cut
  //============================================================================
  // is the W cut in the cuts vector provided?
  bool do_w_cut = std::find(cuts.begin(), cuts.end(), kWexp) != cuts.end();

  bool passes_all_cuts = passes_all_but_w_cut;
  if (do_w_cut)
    passes_all_cuts =
        passes_all_but_w_cut && WexpCut(universe, signal_definition);

  return {passes_all_cuts, is_w_sideband, pion_candidate_idxs};
}

//==============================================================================
// Fuction to count the number of events that pass the cuts
//==============================================================================

EventCount PassedCuts(const CVUniverse& univ,
                      std::vector<int>& pion_candidate_idxs, bool is_mc,
                      SignalDefinition signal_definition,
                      std::vector<ECuts> cuts) {
  pion_candidate_idxs.clear();
  static MichelMap endpoint_michels;
  static MichelMap vertex_michels;
  endpoint_michels.clear();
  vertex_michels.clear();
  EventCount Pass;
  bool pass = true;
  for (auto cu : cuts) Pass[cu] = 0;

  for (auto c : cuts) {
    pass = pass && PassesCut(univ, c, is_mc, signal_definition,
                             endpoint_michels, vertex_michels);
    if (pass) {
      Pass[c] = 1.;
    }
  }

  return Pass;
}

//==============================================================================
// Passes INDIVIDUAL Cut
//==============================================================================
// Updates the michel containers
bool PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
               const SignalDefinition signal_definition,
               MichelMap& endpoint_michels, MichelMap& vertex_michels) {
  const bool useOVMichels = false;
  if (IsPrecut(cut) && !is_mc) return true;

  switch (cut) {
    case kNoCuts:
      return true;

    case kGoodObjects:
      return univ.IsTruth() ? GoodObjectsCut(univ) : true;

    case kGoodVertex:
      return univ.IsTruth() ? GoodVertexCut(univ) : true;

    case kFiducialVolume:
      return univ.IsTruth() ? FiducialVolumeCut(univ) : true;

    case kMinosActivity:
      return univ.IsTruth() ? MinosActivityCut(univ) : true;

    case kPrecuts:
      return univ.IsTruth() ? GoodObjectsCut(univ) && GoodVertexCut(univ) &&
                                  FiducialVolumeCut(univ)
                            : true;
      // MinosActivityCut(univ) : true;

    case kVtx:
      return vtxCut(univ);

    case kMinosMatch:
      return MinosMatchCut(univ);

    case kMinosCharge:
      return MinosChargeCut(univ);

    case kMinosMuon:
      return MinosMatchCut(univ) && MinosChargeCut(univ);

    case kWexp:
      return WexpCut(univ, signal_definition);

    case kIsoProngs:
      return IsoProngCut(univ);

    case kPmu:
      return PmuCut(univ);

    // ==== At Least One Michel ====
    // For now, we need at least one ENDPOINT michel (any # of vtx michels).
    // This cut fills our michel containers, which we use to ID pion tracks
    // and subsequently make track cuts (LLR, node).
    case kAtLeastOneMichel: {
      MichelMap all_michels = GetQualityMichels(univ);
      for (auto m : all_michels) {
        if (m.second.had_idx == -1)
          vertex_michels.insert(m);
        else
          endpoint_michels.insert(m);
      }
      return endpoint_michels.size() > 0;
    }

    case kAtLeastOnePionCandidateTrack:
      return GetQualityPionCandidateIndices(univ).size() > 0;

    // If a michel's pion fails the LLR cut, remove it from the michels
    case kLLR: {
      ContainerEraser::erase_if(endpoint_michels,
                                [&univ](std::pair<int, Michel> mm) {
                                  return !LLRCut(univ, mm.second.had_idx);
                                });
      return endpoint_michels.size() > 0;
    }

    // If a michel's pion fails the node cut, remove it from the michels
    case kNode: {
      ContainerEraser::erase_if(endpoint_michels,
                                [&univ](std::pair<int, Michel> mm) {
                                  return !NodeCut(univ, mm.second.had_idx);
                                });
      return endpoint_michels.size() > 0;
    }

    case kPionMult: {
      if (signal_definition == kOnePi || signal_definition == kOnePiNoW)
        return endpoint_michels.size() == 1 && vertex_michels.size() == 0;
      else
        return endpoint_michels.size() >= 1;
    }

    case kAllCuts:
      return true;

    default:
      std::cout << "PassesCut Error Unknown Cut!" << cut << std::endl;
      return false;
  };
}

//==============================================================================
// Cut functions
//==============================================================================
// Truth precuts
bool GoodObjectsCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_hasGoodObjects");
}
bool GoodVertexCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_isGoodVertex");
}
bool FiducialVolumeCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_isFidVol_smeared");
}
bool MinosActivityCut(const CVUniverse& univ) {
  return univ.GetInt("truth_reco_muon_is_minos_match");
}

// Eventwide reco cuts
bool MinosMatchCut(const CVUniverse& univ) {
  return univ.GetBool("isMinosMatchTrack");
}
// Equivalent to Brandon's, but using standard minos branches
bool MinosChargeCut(const CVUniverse& univ) {
  return univ.GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
}

bool WexpCut(const CVUniverse& univ, SignalDefinition signal_definition) {
  switch (signal_definition) {
    case kOnePi:
    case kNPi:
      return univ.GetWexp() < GetWCutValue(signal_definition);
    case kOnePiNoW:
    case kNPiNoW:
      return true;
    default:
      std::cout << "WexpCut SIGNAL DEF ERROR";
      return false;
  }
}

// cut on max number of iso prongs
// PrimaryBlobProngTool::makeShowerBlobProngs
bool IsoProngCut(const CVUniverse& univ) {
  return univ.GetNIsoProngs() < CCNuPionIncConsts::kIsoProngCutVal;
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

    // VERTEX MICHEL -- must meet the gold standard for its fit
    bool isIntVtx = (vtx == 0);
    if (isIntVtx && current_michel.match_category <= Michel::kNoFit) {
      continue;
    }
    // ENDPOINT MICHELS
    else {
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

bool NodeCut(const CVUniverse& univ, const RecoPionIdx pion_candidate_idx) {
  return 6. < univ.GetEnode01(pion_candidate_idx) &&
         univ.GetEnode01(pion_candidate_idx) < 32. &&
         2. < univ.GetEnode2(pion_candidate_idx) &&
         univ.GetEnode2(pion_candidate_idx) < 22. &&
         0. < univ.GetEnode3(pion_candidate_idx) &&
         univ.GetEnode3(pion_candidate_idx) < 19. &&
         0. < univ.GetEnode4(pion_candidate_idx) &&
         univ.GetEnode4(pion_candidate_idx) < 31. &&
         0. < univ.GetEnode5(pion_candidate_idx) &&
         univ.GetEnode5(pion_candidate_idx) < 60.;
}

bool LLRCut(const CVUniverse& univ, const RecoPionIdx pion_candidate_idx) {
  // if (pion_candidate_idx < 0) return false;
  // else return univ.GetLLRScore(pion_candidate_idx) > 0.;
  return univ.GetLLRScore(pion_candidate_idx) > 0.;
}

// Helper -- get vector of good-enough had tracks, id'ed by unique indices.
std::vector<int> GetQualityPionCandidateIndices(const CVUniverse& univ) {
  std::vector<int> pion_candidate_indices;
  int n_hadrons = univ.GetInt("MasterAnaDev_hadron_number");
  for (int i_hadron = 0; i_hadron != n_hadrons; ++i_hadron)
    if (HadronQualityCuts(univ, i_hadron))
      pion_candidate_indices.push_back(i_hadron);
  return pion_candidate_indices;
}

bool HadronQualityCuts(const CVUniverse& univ,
                       const RecoPionIdx pion_candidate_idx) {
  return univ.GetVecElem("MasterAnaDev_hadron_isForked", pion_candidate_idx) ==
             0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isExiting", pion_candidate_idx) ==
             0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isSideECAL",
                         pion_candidate_idx) == 0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isODMatch", pion_candidate_idx) ==
             0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isTracker", pion_candidate_idx) ==
             1;
};
// Vtx cut for detection volume
bool vtxCut(const CVUniverse& univ) {
  bool pass = true;
  pass = pass && zVertexCut(univ, 8340.0, 5990.0);
  pass = pass && XYVertexCut(univ, 850.0);
  return pass;
}

bool zVertexCut(const CVUniverse& univ, const double upZ, const double downZ) {
  double vtxZ = univ.GetVecElem("vtx", 2);
  if (vtxZ > downZ && vtxZ < upZ)
    return true;
  else
    return false;
}

bool XYVertexCut(const CVUniverse& univ, const double a) {
  const double x = univ.GetVecElem("vtx", 0), y = univ.GetVecElem("vtx", 1);
  if (x < 0) {
    if (x > -a && univ.leftlinesCut(a, x, y))
      return true;
    else
      return false;
  } else {
    if (x < a && univ.rightlinesCut(a, x, y))
      return true;
    else
      return false;
  }
}

bool PmuCut(const CVUniverse& univ) {
  if (univ.GetPmu() / 1000 < 1.5 || 20 < univ.GetPmu() / 1000)
    return false;
  else
    return true;
}

#endif  // Cuts_cxx

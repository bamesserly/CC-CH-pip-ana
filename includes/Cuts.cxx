#ifndef Cuts_cxx
#define Cuts_cxx

#include "Cuts.h"
#include "TruthCategories/Sidebands.h"  // sidebands::kSidebandCutVal
#include "utilities.h"                  // ContainerEraser

//==============================================================================
// Generic Pass Cut(s) Functions
//      * passes, is_sideband, pion_indices = PassesCuts()
//      * bool PassesCut(cut)
//      * PassedCuts <-- just an event counter
//==============================================================================
// Passes All Cuts v3 (latest and greatest)
// return tuple {passes_all_cuts, is_w_sideband, pion_candidate_idxs}
std::tuple<bool, bool, std::vector<int>> PassesCuts(
    CVUniverse& universe, const bool is_mc,
    const SignalDefinition signal_definition, std::vector<ECuts> cuts) {
  //============================================================================
  // passes all cuts but w cut
  //============================================================================
  endpoint::MichelMap endpoint_michels;
  endpoint::MichelMap vertex_michels;
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

// Passes All Cuts v2 (being deprecated)
// fills stuff by reference, instead of returning
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

// Passes All Cuts v1 (being deprecated)
// fills by reference and doesn't check W sideband
bool PassesCuts(CVUniverse& univ, std::vector<int>& pion_candidate_idxs,
                bool is_mc, SignalDefinition signal_definition,
                std::vector<ECuts> cuts) {
  pion_candidate_idxs.clear();
  static endpoint::MichelMap endpoint_michels;
  static endpoint::MichelMap
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

// Count events
EventCount PassedCuts(const CVUniverse& univ,
                      std::vector<int>& pion_candidate_idxs, bool is_mc,
                      SignalDefinition signal_definition,
                      std::vector<ECuts> cuts) {
  pion_candidate_idxs.clear();
  static endpoint::MichelMap endpoint_michels;
  static endpoint::MichelMap vertex_michels;
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

// Pass Single, Given Cut
bool PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
               const SignalDefinition signal_definition,
               endpoint::MichelMap& endpoint_michels, endpoint::MichelMap& vertex_michels) {
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

    case kPmu:
      return PmuCut(univ);

    // ==== At Least One Michel ====
    // For now, we need at least one ENDPOINT michel (any # of vtx michels).
    // This cut fills our michel containers, which we use to ID pion tracks
    // and subsequently make track cuts (LLR, node).
    case kAtLeastOneMichel: {
      endpoint::MichelMap all_michels = endpoint::GetQualityMichels(univ);
      for (auto m : all_michels) {
        if (m.second.had_idx == -1)
          vertex_michels.insert(m);
        else
          endpoint_michels.insert(m);
      }
      vertex::MichelEvent mehreen_michels = vertex::GetQualityMichels(univ);
      return endpoint_michels.size() > 0 /*|| mehreen_michels.size() = 0*/;
    }

    case kAtLeastOnePionCandidateTrack:
      return GetQualityPionCandidateIndices(univ).size() > 0;

    // If a michel's pion fails the LLR cut, remove it from the michels
    case kLLR: {
      ContainerEraser::erase_if(endpoint_michels,
                                [&univ](std::pair<int, endpoint::Michel> mm) {
                                  return !LLRCut(univ, mm.second.had_idx);
                                });
      return endpoint_michels.size() > 0;
    }

    // If a michel's pion fails the node cut, remove it from the michels
    case kNode: {
      ContainerEraser::erase_if(endpoint_michels,
                                [&univ](std::pair<int, endpoint::Michel> mm) {
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
// Cut Definitions
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

// Exclusive - cuts on pion tracks
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


//==============================================================================
// Helper
//==============================================================================
// Get candidate pions that pass the minimal HadronQualityCuts
std::vector<int> GetQualityPionCandidateIndices(const CVUniverse& univ) {
  std::vector<int> pion_candidate_indices;
  int n_hadrons = univ.GetInt("MasterAnaDev_hadron_number");
  for (int i_hadron = 0; i_hadron != n_hadrons; ++i_hadron)
    if (HadronQualityCuts(univ, i_hadron))
      pion_candidate_indices.push_back(i_hadron);
  return pion_candidate_indices;
}

#endif  // Cuts_cxx

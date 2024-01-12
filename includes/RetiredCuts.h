//==============================================================================
// Cleaned out Cuts.* files.
// Don't plan on using these cuts anymore and don't try to compile them.
//==============================================================================
#ifndef RetiredCuts_H
#define RetiredCuts_H

#include "CVUniverse.h"
#include "Constants.h"  // ECuts
#include "SignalDefinition.h"

bool DeadTimeCut(const CVUniverse&);
bool DeadTimeCut(const CVUniverse& univ) { return univ.GetInt("tdead") <= 1; }

bool MinosCoilCut(const CVUniverse&);
bool MinosCoilCut(const CVUniverse& univ) {
  const double MINOS_COIL_RADIUS = 210;  // mm
  const double MAX_MINOS_RADIUS = 2500;  // mm
  const double coilXPos = 1219.0;
  const double coilYPos = 393.0;
  const double minos_x =
      univ.GetDouble("MasterAnaDev_minos_trk_end_x") + coilXPos;
  const double minos_y =
      univ.GetDouble("MasterAnaDev_minos_trk_end_y") + coilYPos;
  double minosR = sqrt(pow(minos_x, 2) + pow(minos_y, 2));
  // if (!((pow(minos_x,2) + pow(minos_y,2) )>= pow(MINOS_COIL_RADIUS, 2)) )
  //  cout << minos_x << " " << minos_y << " " << MINOS_COIL_RADIUS << endl;
  return (minosR > MINOS_COIL_RADIUS && minosR < MAX_MINOS_RADIUS);
}

bool IsoBlobCut(const CVUniverse&);
bool IsoBlobCut(const CVUniverse& univ) {
  return univ.GetInt("n_iso_blob_prongs") <
         1;  // RecoilUtils::createIsoBlobProngs
}

bool IsoProngSepCut(const CVUniverse&);
bool IsoProngSepCut(const CVUniverse& univ) {
  return univ.GetLargestIsoProngSep() < 300;
}

bool ThetamuCut(const CVUniverse&);
bool ThetamuCut(const CVUniverse& univ) { return univ.GetThetamu() < 0.3491; }

bool BrandonMinosChargeCut(const CVUniverse&);
bool BrandonMinosChargeCut(CVUniverse& univ) {
  return univ.GetDouble("MasterAnaDev_muon_qpqpe") < 0.0;
}

bool CCIncMinosChargeCut(const CVUniverse&);
bool CCIncMinosChargeCut(CVUniverse& univ) {
  if (univ.GetBool("MasterAnaDev_minos_used_curvature"))
    return 1. / univ.GetDouble("MasterAnaDev_minos_trk_eqp_qp") < -5.0;
  else if (univ.GetBool("MasterAnaDev_minos_used_range"))
    return univ.GetBool("MasterAnaDev_minos_trk_qp") < 0.0;
  else
    return false;
}

bool ExactlyOneEndpointMichelCut(const CVUniverse&, SignalDefinition);
bool ExactlyOneEndpointMichelCut(const CVUniverse& univ,
                                 SignalDefinition signal_definition) {
  if (signal_definition == kNPi || signal_definition == kNPiNoW) {
    return true;
  } else if (signal_definition == kOnePi || signal_definition == kOnePiNoW) {
    endpoint::MichelMap mm = endpoint::GetQualityMichels(univ);
    if (mm.size() == 1) {  // require only one michel
      endpoint::Michel m = (mm.begin())->second;
      if (m.vtx == 0)
        return false;  // so this has a no-vertex michel cut baked in
      // pion_candidate_idx = m.vtx - 1;   // SELECT OUR PION
      return true;
    } else
      return false;
  } else {
    std::cout << "ExactlyOneEndpointMichelcut SIGNAL DEFINITION ERROR"
              << std::endl;
    return false;
  }
}

bool AtLeastOneBrandonMichelCut(const CVUniverse&);
bool AtLeastOneBrandonMichelCut(const CVUniverse& univ) {
  // Get quality hadron candidates
  std::vector<int> pion_candidate_indices =
      GetQualityPionCandidateIndices(univ);

  // Loop quality hadron candidates to see if they have a good brandon michel
  for (auto pion_candidate_idx : pion_candidate_indices) {
    int michel_views = univ.GetVecElem("MasterAnaDev_hadron_endMichel_category",
                                       pion_candidate_idx);
    int michel_ndigits = univ.GetVecElem(
        "MasterAnaDev_hadron_endMichel_ndigits", pion_candidate_idx);
    double michel_energy =
        univ.GetVecElem("MasterAnaDev_hadron_endMichel_energy",
                        pion_candidate_idx);  /// TODO sys universe function
    double michel_slice_energy = univ.GetVecElem(
        "MasterAnaDev_hadron_endMichel_slice_energy", pion_candidate_idx);

    if (michel_views < 1)
      continue;                    // no michel
    else if (michel_views == 1) {  // 1 view
      if (michel_energy < 55.0 && michel_ndigits < 35 &&
          michel_slice_energy < 100.0)
        return true;
      else
        continue;
    } else if (michel_views > 1) {  // 2+3 views
      if (michel_energy < 55.0 && michel_ndigits < 35 &&
          michel_ndigits >= michel_views)
        return true;
      else
        continue;
    }
  }
  return false;
}

// ntrk are tracks not including the muon prong
bool AtLeastOneAnchoredProngCut(const CVUniverse&);
bool AtLeastOneAnchoredProngCut(const CVUniverse& univ) {
  int ntrk = univ.GetInt("n_anchored_long_trk_prongs") +
             univ.GetInt("n_anchored_short_trk_prongs");
  return ntrk > 0;
}

bool AtLeastOneNodeCandidateCut(const CVUniverse&);
bool AtLeastOneNodeCandidateCut(const CVUniverse& univ) {
  // Get quality hadron candidates
  std::vector<int> pion_candidate_indices =
      GetQualityPionCandidateIndices(univ);
  for (auto pion_candidate_idx : pion_candidate_indices) {
    if (NodeCut(univ, pion_candidate_idx)) return true;
  }
  return false;
}

bool AtLeastOneLLRCandidateCut(const CVUniverse&);
bool AtLeastOneLLRCandidateCut(const CVUniverse& univ) {
  // Get quality hadron candidates
  std::vector<int> pion_candidate_indices =
      GetQualityPionCandidateIndices(univ);
  for (auto pion_candidate_idx : pion_candidate_indices) {
    if (LLRCut(univ, pion_candidate_idx)) return true;
  }
  return false;
}

// Call the cut functions -- fill-in/return the ref to the good pion candidate
// PassesCuts v1 (being deprecated))
bool PassesCuts(CVUniverse&, std::vector<int>& pion_candidate_idxs, bool is_mc,
                SignalDefinition, std::vector<ECuts> cuts = kCutsVector);
// Check all cuts
// this one doesn't keep track of w cut, so we'll be retiring it.
// And I don't think it's used.
bool PassesCuts(CVUniverse& univ, std::vector<int>& pion_candidate_idxs,
                bool is_mc, SignalDefinition signal_definition,
                std::vector<ECuts> cuts) {
  pion_candidate_idxs.clear();
  static endpoint::MichelMap endpoint_michels;
  static endpoint::MichelMap
      vtx_michels;  // Keep track of these, but not used currently
  endpoint_michels.clear();
  vtx_michels.clear();
  bool pass = true;
  for (auto c : cuts) {
    // Set the pion candidates to the universe
    // GetHadIdxsFromMichels takes a trackless::MichelEvent now. To remain
    // backwards compatible, pass it a dummy.
    univ.SetPionCandidates(GetHadIdxsFromMichels(endpoint_michels));
    pass = pass && PassesCut(univ, c, is_mc, signal_definition,
                             endpoint_michels, vtx_michels);
  }

  // Each endpoint michel has an associated hadron track.
  // Our official pion candidates are those tracks.
  pion_candidate_idxs = GetHadIdxsFromMichels(endpoint_michels);

  return pass;
}

// also tell whether we are w sideband
// PassesCuts v2 (being deprecated)
bool PassesCuts(CVUniverse&, std::vector<int>& pion_candidate_idxs,
                const bool is_mc, const SignalDefinition, bool& is_w_sideband,
                std::vector<ECuts> cuts = kCutsVector);
// Check all cuts AND whether we are W sideband. Save a lot of time.
// Strategy is:
// 1. check all cuts except for W
// 2. check if W > 1.5 (sideband)
// 3. check if W < 1.4 (signal)
//
// BUT works with side effects and fills stuff by reference so:
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

// PassesCut v1 (being deprecated))
bool PassesCut(const CVUniverse&, const ECuts cut, const bool is_mc,
               const SignalDefinition, endpoint::MichelMap& endpoint_michels,
               endpoint::MichelMap& vertex_michels);
// Individual cut
// Works with side effects and fills by reference so:
// BEING DEPRECATED
bool PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
               const SignalDefinition signal_definition,
               endpoint::MichelMap& endpoint_michels,
               endpoint::MichelMap& vtx_michels) {
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
      endpoint::MichelMap all_michels = endpoint::GetQualityMichels(univ);
      for (auto m : all_michels) {
        if (m.second.had_idx == -1)
          vtx_michels.insert(m);
        else
          endpoint_michels.insert(m);
      }
      trackless::MichelEvent mehreen_michels =
          trackless::GetQualityMichels(univ);
      return endpoint_michels.size() > 0;  // || mehreen_michels.size() = 0;
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
        return endpoint_michels.size() == 1 && vtx_michels.size() == 0;
      else
        return endpoint_michels.size() >= 1;
    }

    case kAllCuts:
      return true;

    default:
      std::cout << "PassesCut Error Unknown Cut!" << cut << "  "
                << GetCutName(cut) << "\n";
      return false;
  };
}

// Event counter that uses the old functions
//==============================================================================
// Fuction to count the number of events that pass the cuts
// TODO this should be renamed
//==============================================================================
bool PassesCuts(CCPiEvent&, bool& is_w_sideband);
EventCount PassedCuts(const CVUniverse& univ,
                      std::vector<int>& pion_candidate_idxs, bool is_mc,
                      SignalDefinition signal_definition,
                      std::vector<ECuts> cuts) {
  pion_candidate_idxs.clear();
  static endpoint::MichelMap endpoint_michels;
  static endpoint::MichelMap vtx_michels;
  endpoint_michels.clear();
  vtx_michels.clear();
  EventCount Pass;
  bool pass = true;
  for (auto cu : cuts) Pass[cu] = 0;

  for (auto c : cuts) {
    pass = pass && PassesCut(univ, c, is_mc, signal_definition,
                             endpoint_michels, vtx_michels);
    if (pass) {
      Pass[c] = 1.;
    }
  }

  return Pass;
}

// Used in analysis pipeline
// Uses PassesCuts v2. Does check w sideband, but fills by reference instead of
// returning its results. v3 is the future.
bool PassesCuts(CCPiEvent& e, bool& is_w_sideband) {
  return PassesCuts(*e.m_universe, e.m_reco_pion_candidate_idxs, e.m_is_mc,
                    e.m_signal_definition, is_w_sideband);
}

// Uses PassesCuts v1.
bool PassesCuts(CCPiEvent&, std::vector<ECuts> cuts);
// No longer used anywhere. Doesn't check w sideband while looping all cuts.
// Nothing wrong with it per se. Checking the w sideband is just practically
// free. v3 of PassesCuts is the future, anyways.
bool PassesCuts(CCPiEvent& e, std::vector<ECuts> cuts) {
  return PassesCuts(*e.m_universe, e.m_reco_pion_candidate_idxs, e.m_is_mc,
                    e.m_signal_definition, cuts);
}

#endif

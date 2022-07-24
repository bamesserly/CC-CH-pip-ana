//==============================================================================
// Cleaned out Cuts.* files.
// Don't plan on using these cuts anymore and don't try to compile them.
//==============================================================================
#ifndef RetiredCuts_H
#define RetiredCuts_H

#include "Constants.h" // ECuts
#include "CVUniverse.h"
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

#endif

#ifndef PIONCANDIDATE_H
#define PIONCANDIDATE_H

struct PionCandidate {
  bool is_tracked;
  bool is_untracked;
  MichelEvent michel_event;
  MichelMap michel_map;
};

std::vector<int> GetHadIdxsFromMichels(const endpoint::MichelMap& endpoint_michels) {
  std::vector<int> ret;
  for (auto m : endpoint_michels) ret.push_back(m.second.had_idx);
  return ret;
}

std::vector<PionCandidate> GetPionCandidates_Tracked() {
  endpoint::MichelMap endpoint_michels = endpoint::GetQualityMichels(univ);
}





bool NodeCut(const CVUniverse& univ, const RecoPionIdx pidx) {
  return 6. < univ.GetEnode01(pidx) && univ.GetEnode01(pidx) < 32. &&
         2. < univ.GetEnode2(pidx) && univ.GetEnode2(pidx) < 22. &&
         0. < univ.GetEnode3(pidx) && univ.GetEnode3(pidx) < 19. &&
         0. < univ.GetEnode4(pidx) && univ.GetEnode4(pidx) < 31. &&
         0. < univ.GetEnode5(pidx) && univ.GetEnode5(pidx) < 60.;
}

bool LLRCut(const CVUniverse& univ, const RecoPionIdx pion_candidate_idx) {
  // if (pion_candidate_idx < 0) return false;
  // else return univ.GetLLRScore(pion_candidate_idx) > 0.;
  return univ.GetLLRScore(pion_candidate_idx) > 0.;
}

// Get candidate pions that pass the minimal HadronQualityCuts
std::vector<int> GetQualityPionCandidateIndices(const CVUniverse& univ) {
  std::vector<int> pion_candidate_indices;
  int n_hadrons = univ.GetInt("MasterAnaDev_hadron_number");
  for (int i_hadron = 0; i_hadron != n_hadrons; ++i_hadron)
    if (HadronQualityCuts(univ, i_hadron))
      pion_candidate_indices.push_back(i_hadron);
  return pion_candidate_indices;
}

bool HadronQualityCuts(const CVUniverse& univ, const RecoPionIdx pidx) {
  return univ.GetVecElem("MasterAnaDev_hadron_isForked", pidx) == 0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isExiting", pidx) == 0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isSideECAL", pidx) == 0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isODMatch", pidx) == 0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isTracker", pidx) == 1;
};

bool PassesPionMultiplicityCut(const CCPiEvent& event) {
  PassesPionMultiplicityCut(event.m_passed_pion_candidates.size(), signal_definition);
}
bool PassesPionMultiplicityCut(const unsigned int n_passed_pion_candidates, SignalDefinition signal_definition) {
  return signal_definition.m_n_pi_min < n_passed_pion_candidates && n_passed_pion_candidates < signal_definition.m_n_pi_max;
}

#endif

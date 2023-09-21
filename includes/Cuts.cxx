//============================================================================
/*
This file contains the definitions of individual reco (aka event selection)
cuts.

It also contains the PassesCut(s) functions to apply all of them at once to
perform the event selection.

The default cuts used for the analysis are in the kCutsVector, located in
includes/CutUtils.h

A not-brief note on how the "exclusive" pion cuts work:

  The cuts system makes "exclusive" cuts on michels and on their associated
  pion tracks. In the end, PassesCuts spits out vectors of michels that passed
  the exclusive cuts.

  The way it has worked for a long time is that for the rest of the analysis
  (in particular, when calculating the exclusive/pion properties of events),
  we would switch from operating in michels to operating in pion track
  indices.

  But now we've got trackless michels. At the moment I'm restricting to one
  trackless michel -- the best one. When a good trackless michel is found,
  convert/code it to a pion track index of -1. When calculating pion
  quantities with a track index of -1, we'll do something special to calculate
  those quantities.

  Going forward, we could want to use more than one trackless michels.  In
  that case, we may need to move away from calculating pion quantities
  directly from track indices. And instead extract pion quantities from the
  michels themselves.
*/
//============================================================================
#ifndef Cuts_cxx
#define Cuts_cxx

#include "Cuts.h"

#include "CutUtils.h"  // GetHadIdxsFromMichels, IsPrecut, GetWSidebandCuts, kCutsVector
#include "Michel.h"  // endpoint::Michel, endpoint::MichelMap, endpoint::GetQualityMichels
#include "TruthCategories/Sidebands.h"  // sidebands::kSidebandCutVal
#include "utilities.h"                  // ContainerEraser

bool PassesCuts(const CCPiEvent& event,
                const SignalDefinition& signal_definition,
                const std::vector<ECuts>& cuts) {
  bool pass = true;
  for (auto cut : cuts) {
    pass = pass && PassesCut(event, cut, signal_definition);
  }
  return pass;
}

bool PassesCut(const CCPiEvent& event, const ECuts& cut,
               const SignalDefinition& signal_definition) {
  bool pass = false;
  CVUniverse univ = *event.m_universe;

  // throw these away
  endpoint::MichelMap tracked_michels;
  LowRecoilPion::MichelEvent<CVUniverse> vtx_michels;

  if (IsPrecut(cut) && !event.m_is_mc) return true;

  switch (cut) {
    case kNoCuts:
      pass = true;
      break;

    // gaudi cut AKA precut (maybe still used in MAD?)
    case kGoodObjects:
      pass = univ.IsTruth() ? GoodObjectsCut(univ) : true;
      break;

    // gaudi cut AKA precut (maybe still used in MAD?)
    case kGoodVertex:
      pass = univ.IsTruth() ? GoodVertexCut(univ) : true;
      break;

    // gaudi cut AKA precut (probably not used in MAD)
    case kFiducialVolume:
      pass = univ.IsTruth() ? FiducialVolumeCut(univ) : true;
      break;

    // gaudi cut AKA precut (probably not used in MAD)
    case kMinosActivity:
      pass = univ.IsTruth() ? MinosActivityCut(univ) : true;
      break;

    case kPrecuts:
      pass = univ.IsTruth() ? GoodObjectsCut(univ) && GoodVertexCut(univ) &&
                                  FiducialVolumeCut(univ)
                            : true;
      // MinosActivityCut(univ) : true;
      break;

    case kVtx:
      pass = vtxCut(univ);
      break;

    case kMinosMatch:
      pass = MinosMatchCut(univ);
      break;

    case kMinosCharge:
      pass = MinosChargeCut(univ);
      break;

    case kMinosMuon:
      pass = MinosMatchCut(univ) && MinosChargeCut(univ);
      break;

    case kWexp:
      pass = WexpCut(univ, signal_definition);
      break;

    case kIsoProngs:
      pass = IsoProngCut(univ);
      break;

    case kPmu:
      pass = PmuCut(univ);
      break;

    // modify michels
    case kAtLeastOneMichel: {
      tracked_michels = endpoint::GetQualityMichels(univ);
      vtx_michels = LowRecoilPion::MichelEvent<CVUniverse>();  // LowRecoilPion::GetQualityMichels<CVUniverse>(univ);
      pass = tracked_michels.size() > 0;  // || vtx_michels.m_idx != -1;
      break;
    }

    // If a michel's pion fails the LLR cut, remove it from the michels
    case kLLR: {
      ContainerEraser::erase_if(tracked_michels,
                                [&univ](std::pair<int, endpoint::Michel> mm) {
                                  return !LLRCut(univ, mm.second.had_idx);
                                });
      pass = tracked_michels.size() > 0;  // || vtx_michels.m_idx != -1;
      break;
    }

    // modify michels
    // If a michel's pion fails the node cut, remove it from the michels
    case kNode: {
      ContainerEraser::erase_if(tracked_michels,
                                [&univ](std::pair<int, endpoint::Michel> mm) {
                                  return !NodeCut(univ, mm.second.had_idx);
                                });
      pass = tracked_michels.size() > 0;  // || vtx_michels.m_idx != -1;
      break;
    }

    // modify michels
    // If a michel's pion fails track quality, remove it from the michels
    case kTrackQuality: {
      ContainerEraser::erase_if(
          tracked_michels, [&univ](std::pair<int, endpoint::Michel> mm) {
            return !HadronQualityCuts(univ, mm.second.had_idx);
          });
      pass = tracked_michels.size() > 0;  // || vtx_michels.m_idx != -1;
      break;
    }

    // modify michels
    // the quality track, LLR, and node cuts may have removed the michels of
    // failed tracks
    case kAtLeastOnePionCandidate:
      pass = tracked_michels.size() > 0;  // || vtx_michels.m_idx != -1;
      break;

    // TODO check trackless michels and no double-counting
    case kPionMult: {
      pass = signal_definition.m_n_pi_min <= tracked_michels.size() &&
             tracked_michels.size() <= signal_definition.m_n_pi_max;
      break;
    }

    case kAtLeastOnePionCandidateTrack:
      pass = GetQualityPionCandidateIndices(univ).size() >
             0;  // ||
                 // vtx_michels.m_idx != -1;
      break;

    case kAllCuts:
      pass = true;
      break;

    default:
      std::cout << "PassesCut Error Unknown Cut!" << cut << "  "
                << GetCutName(cut) << "\n";
      pass = false;
  };

  return pass;
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

bool WexpCut(const CVUniverse& univ, const SignalDefinition signal_definition) {
  return univ.GetWexp() < signal_definition.m_w_max;
}

// cut on max number of iso prongs
// PrimaryBlobProngTool::makeShowerBlobProngs
bool IsoProngCut(const CVUniverse& univ) {
  return univ.GetNIsoProngs() < CCNuPionIncConsts::kIsoProngCutVal;
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

// Vtx cut for detection volume
bool vtxCut(const CVUniverse& univ) {
  bool pass = true;
  pass = pass && zVertexCut(univ, CCNuPionIncConsts::kZVtxMaxCutVal,
                            CCNuPionIncConsts::kZVtxMinCutVal);
  pass = pass && XYVertexCut(univ, CCNuPionIncConsts::kApothemCutVal);
  return pass;
}

bool zVertexCut(const CVUniverse& univ, const double upZ, const double downZ) {
  double vtxZ = univ.GetVecElem("vtx", 2);
  return vtxZ > downZ && vtxZ < upZ;
}

bool XYVertexCut(const CVUniverse& univ, const double a) {
  const double x = univ.GetVecElem("vtx", 0);
  const double y = univ.GetVecElem("vtx", 1);
  return univ.IsInHexagon(x, y, a);
}

bool PmuCut(const CVUniverse& univ) {
  double pmu = univ.GetPmu();
  return CCNuPionIncConsts::kPmuMinCutVal < pmu &&
         pmu < CCNuPionIncConsts::kPmuMaxCutVal;
}

#endif  // Cuts_cxx

#ifndef Hadron_H
#define Hadron_H

#include "../CVUniverse.h"

enum HadronType {
  kStoppingPip, kInteractPip, kStoppingP, kInteractP, kPim, /*kPi0, kNeutron,
  kKaon, kSigma, kMuon,*/ kOther, kNHadronTypes
};

HadronType GetHadronType(const CVUniverse& universe, RecoPionIdx hadron_id) {
  int pdg = universe.GetVecElem("MasterAnaDev_hadron_tm_PDGCode", hadron_id);
  int destruct_code = universe.GetVecElem("MasterAnaDev_hadron_tm_destructCode", hadron_id);
  //std::cout << hadron_id << " " << pdg << " " <<  destruct_code << std::endl;
    //double frac = universe.GetVecElem("MasterAnaDev_hadron_tm_fraction", hadron_id);
    //if (frac < 0.5) 
    //  std::cout << "frac = " << frac << std::endl;

    //int aaron_pdg = universe.GetVecElem("MasterAnaDev_hadron_PDGCode", hadron_id); // Identical to other method.

    // Scan for specific particles
    //if ( abs(pdg == 321) && universe.ShortName() == "cv") {
    //  int link_size = 200;
    //  char link [link_size];
    //  universe.ArachneLink(link, link_size);
    //  //std::cout.precision(3);
    //  std::cout << link << std::endl;
    //}

  if ( pdg == 211 && destruct_code == 0 )
    return kStoppingPip;

  else if ( pdg == 211 && destruct_code != 0 )
    return kInteractPip;

  else if ( pdg == -211 && destruct_code == 0 )
    return kPim;

  else if ( pdg == -211 && destruct_code != 0 )
    return kPim;

  else if ( pdg == 2212 && destruct_code == 0 )
    return kStoppingP;

  else if ( pdg == 2212 && destruct_code != 0 )
    return kInteractP;

  /*
  else if ( pdg == 111)
    return kPi0;

  else if ( pdg == 2112 )
    return kNeutron;

  else if ( abs(pdg) == 311 || abs(pdg) == 321 )
    return kKaon;

  else if ( abs(pdg) == 3222)
    return kSigma;

  else if ( abs(pdg) == 13) {
    return kMuon;
  }
  */

  else {
    return kOther;
  }
}

std::string GetTruthClassification_LegendLabel(HadronType hadron_category) {
  switch (hadron_category) {
    case kStoppingPip: return "Stopping #pi^{+}";
    case kInteractPip: return "Interacting #pi^{+}";
    case kPim:         return "#pi^{-}";
    case kStoppingP:   return "Stopping Proton";
    case kInteractP:   return "Interacting Proton";
    //case kPi0:         return "#pi^{0}";
    /*
    case kNeutron:     return "n";
    case kKaon:        return "K";
    case kSigma:       return "#Sigma";
    case kMuon:        return "#mu";
    */
    default:           return "Other";
  }
}

std::string GetTruthClassification_Name(HadronType hadron_category) {
  switch(hadron_category) {
    case kStoppingPip: return "stopPip";
    case kInteractPip: return "interactPip";
    case kPim:         return "Pim";
    case kStoppingP:   return "stopP";
    case kInteractP:   return "interactP";
    //case kPi0:         return "Pi0";
    /*
    case kNeutron:     return "n";
    case kKaon:        return "K";
    case kSigma:       return "Sigma";
    case kMuon:        return "mu";
    */
    default:           return "other"; 
  }
}

#endif

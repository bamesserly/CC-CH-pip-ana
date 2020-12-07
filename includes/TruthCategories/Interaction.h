#ifndef Interaction_H
#define Interaction_H

#include <stdlib.h> // exit, EXIT_FAILURE

#include "../CVUniverse.h"

enum FSParticleType
{ 
  // !CCnumuNpi
  kOtherInt, //k0pi, 


  // CCnumuNpi
  kCCnumu_0pip_NpiX, kCCnumu_2pip_NpiX, kCCnumu_1pip_NpiX,  // BGs
  kCCnumu_2pip_0piX, kCCnumu_1pip_0piX,                    // Signals

  //kNC_NpiX, knue_NpiX, knumubar_NpiX,   // almost none of these.
                                          // adding them to kOtherInt

  kNFSParticleTypes
};

enum ChannelType
{ 
  kCCQE, kCOHPI, kNONRES, kRES, kNChannelTypes
};

//==============================================================================

FSParticleType GetFSParticleType(const CVUniverse& universe){
  int N_pip = universe.GetInt("truth_N_pip");
  int N_pim = universe.GetInt("truth_N_pim");
  int N_pi0 = universe.GetInt("truth_N_pi0");

  //Get off the bus in the following order please:
  // 0pi, NC, nue, antinumu, !numu, !CC
  if( N_pip == 0 && N_pim == 0 && N_pi0 == 0 ) return kOtherInt; //k0pi;
  if( universe.GetInt("mc_current") == 2 )     return kOtherInt; //kNC_NpiX;
  if( universe.GetInt("mc_incoming") == 12)    return kOtherInt; //knue_NpiX;
  if( universe.GetInt("mc_incoming") == -14)   return kOtherInt; //knumubar_NpiX;
  if( universe.GetInt("mc_incoming") != 14)    return kOtherInt; 
  if( universe.GetInt("mc_current") != 1)      return kOtherInt;

  //still on the bus? you're CC numu Npi (N>0)
  
  //0pi+, NpiX
  if( N_pip == 0 ) return kCCnumu_0pip_NpiX;

  //1pi+, 0 or N piX
  else if( N_pip == 1 ){
    if( N_pim == 0 && N_pi0 == 0 ) return kCCnumu_1pip_0piX;
    else return kCCnumu_1pip_NpiX;
  }

  //Npi+, 0 or N piX
  else if( N_pip == 2 ){
    if( N_pim == 0 && N_pi0 == 0 ) return kCCnumu_2pip_0piX;
    else return kCCnumu_2pip_NpiX;
  }

  else return kOtherInt;
}

std::string GetTruthClassification_LegendLabel(FSParticleType interaction_category){
  switch (interaction_category){
    //case k0pi:
    //  return "No #pi";
    //case kNC_NpiX:
    //  return "NC N#pi";
    //case knue_NpiX:
    //  return "#nu_{e} N#pi";
    //case knumubar_NpiX: 
    //  return "#bar{#nu_{#mu}} N#pi";
    case kOtherInt:
      return "Other";
    case kCCnumu_0pip_NpiX:
      return "#nu_{#mu}CC0#pi^{+}+N(#pi^{0},#pi^{-})";
    case kCCnumu_1pip_0piX:
      return "#nu_{#mu}CC1#pi^{+}";
    case kCCnumu_1pip_NpiX:
      return "#nu_{#mu}CC1#pi^{+}+N(#pi^{0},#pi^{-})";
    case kCCnumu_2pip_0piX:
      return "#nu_{#mu}CC2#pi^{+}";
    case kCCnumu_2pip_NpiX:
      return "#nu_{#mu}CC2#pi^{+}+N(#pi^{0},#pi^{-})";
    default: 
      return "Other";
  };
}

std::string GetTruthClassification_Name(FSParticleType interaction_category){
  switch (interaction_category){
    //case k0pi:
    //  return "0pi";
    //case kNC_NpiX:
    //  return "NC_NpiX";
    //case knue_NpiX:
    //  return "nue_NpiX";
    //case knumubar_NpiX: 
    //  return "numubar_NpiX";
    case kOtherInt:
      return "other";
    case kCCnumu_0pip_NpiX:
      return "CCnumu_0pip_NpiX";
    case kCCnumu_1pip_0piX:
      return "CCnumu_1pip_0piX"; 
    case kCCnumu_1pip_NpiX:
      return "CCnumu_1pip_NpiX";
    case kCCnumu_2pip_0piX:
      return "CCnumu_2pip_0piX";
    case kCCnumu_2pip_NpiX:
      return "CCnumu_2pip_NpiX";
    default:
      return "other";
  };
}

//==============================================================================

//kNoInteraction=0, kQEL=1, kRES=2, kDIS=3, kCOHPI=4, kAMNUGAMMA=5, kIMD=6, kNUEEL=7, k2P2H=8
ChannelType GetChannelType(const CVUniverse& universe){
  int category = universe.GetInt("mc_intType");
  if(category < 0) exit(EXIT_FAILURE);
  else if (category == 1) return kCCQE;
  else if (category == 2) return kRES;
  else if (category == 4) return kCOHPI;
  else return kNONRES;
}

std::string GetTruthClassification_LegendLabel(ChannelType interaction2_category){
  switch (interaction2_category){
    case kCCQE:
      return "CCQE";
    case kRES:
      return "Resonance";
    case kCOHPI:
      return "Coherent";
    case kNONRES:
      return "Non-Resonance";
    default: 
      return "ERROR";
  };
}

std::string GetTruthClassification_Name(ChannelType interaction2_category){
  switch (interaction2_category){
    case kCCQE:
      return "CCQE";
    case kRES:
      return "Res";
    case kCOHPI:
      return "Coh";
    case kNONRES:
      return "NonRes";
    default: 
      return "ERROR";
  };
}
#endif

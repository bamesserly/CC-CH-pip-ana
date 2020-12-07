#ifndef NPions_H
#define NPions_H

#include <stdlib.h> // exit, EXIT_FAILURE

#include "../CVUniverse.h"

enum NPionsType
{
  kThreeOrMorePions, kTwoPions, kOnePion, kZeroPions, kNNPionsTypes
};

enum NPi0Type
{
  kTwoOrMorePi0, kOnePi0, kZeroPi0, kNNPi0Types
};

enum NPipType
{
  kZeroPip, kOnePip, kTwoPip, kThreeOrMorePip, kNNPipTypes
};

//==============================================================================

NPionsType GetNPionsType(const CVUniverse& universe){
  int n_pions = universe.GetInt("truth_N_pip") + universe.GetInt("truth_N_pim") + universe.GetInt("truth_N_pi0");
  if (n_pions < 0) exit(EXIT_FAILURE);
  else if (n_pions == 0) return kZeroPions;
  else if (n_pions == 1) return kOnePion;
  else if (n_pions == 2) return kTwoPions;
  else return kThreeOrMorePions;

  //else if (n_pions <= 2) return static_cast<NPionsType>(n_pions);
}

std::string GetTruthClassification_LegendLabel(NPionsType category){
  switch(category){
    case kZeroPions:
      return "0 #pi";
    case kOnePion:
      return "1 #pi";
    case kTwoPions:
      return "2 #pi";
    case kThreeOrMorePions:
      return ">2 #pi";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(NPionsType category){
  switch(category){
    case kZeroPions:
      return "zero_pi";
    case kOnePion:
      return "one_pi";
    case kTwoPions:
      return "two_pi";
    case kThreeOrMorePions:
      return "three_or_more_pi";
    default:
      return "ERROR";
  }
}

//==============================================================================

NPi0Type GetNPi0Type(const CVUniverse& universe){
  int n_pions = universe.GetInt("truth_N_pi0");
  if (n_pions < 0) exit(EXIT_FAILURE);
  else if (n_pions == 0 )return kZeroPi0;
  else if (n_pions == 1) return kOnePi0;
  else return kTwoOrMorePi0;

  //else if (n_pions <= 1) return static_cast<NPi0Type>(n_pions);
}

std::string GetTruthClassification_LegendLabel(NPi0Type category){
  switch(category){
    case kZeroPi0:
      return "0 #pi^{0}";
    case kOnePi0:
      return "1 #pi^{0}";
    case kTwoOrMorePi0:
      return ">1 #pi^{0}";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(NPi0Type category){
  switch(category){
    case kZeroPi0:
      return "zero_pi";
    case kOnePi0:
      return "one_pi";
    case kTwoOrMorePi0:
      return "two_or_more_pi";
    default:
      return "ERROR";
  }
}

//==============================================================================

NPipType GetNPipType(const CVUniverse& universe){
  int n_pions = universe.GetInt("truth_N_pip");
  if (n_pions < 0) exit(EXIT_FAILURE);
  else if (n_pions <= 2) return static_cast<NPipType>(n_pions);
  else return kThreeOrMorePip;
}

std::string GetTruthClassification_LegendLabel(NPipType category){
  switch(category){
    case kZeroPip:
      return "0 #pi^{+}";
    case kOnePip:
      return "1 #pi^{+}";
    case kTwoPip:
      return "2 #pi^{+}";
    case kThreeOrMorePip:
      return ">2 #pi^{+}";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(NPipType category){
  switch(category){
    case kZeroPip:
      return "zero_pip";
    case kOnePip:
      return "one_pip";
    case kTwoPip:
      return "two_pip";
    case kThreeOrMorePip:
      return "three_or_more_pip";
    default:
      return "ERROR";
  }
}

#endif

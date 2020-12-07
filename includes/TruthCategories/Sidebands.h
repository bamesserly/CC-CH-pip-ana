#ifndef Sidebands_H
#define Sidebands_H

#include <algorithm>
#include "../CVUniverse.h"
#include "../SignalDefinition.h"
#include "../Cuts.h"

namespace sidebands{
  const int kNWFitCategories = 3;
  const std::string kFitVarString("wexp_fit");
  const double kSidebandCutVal = 1500;
  const int kWSideband_ColorScheme = 3;
}

enum WSidebandType
{
  kWSideband_Signal, kWSideband_Low, kWSideband_Mid, kWSideband_High, kNWSidebandTypes
};


// n_categories: should we chop up the sideband into 2 or 3 categories?
// n_categories == 2: WexpTrue [0, 1.8],[1.8, up]
// n_categories == 3: WexpTrue [0, 1.4],[1.4, 1.8],[1.8, up]
// This is separate from the question of whether we float or fix each category.
// For that question, see the implementation of the fitter.
WSidebandType GetWSidebandType(const CVUniverse& universe, 
                               const SignalDefinition signal_definition, 
                               const int n_categories = 3) {
  if (n_categories != 2 && n_categories != 3)
    throw std::invalid_argument("GetWSidebandType: n_categories is 2 or 3");

  double Wexp_true = universe.GetWexpTrue();

  if (IsSignal(universe, signal_definition))
    return kWSideband_Signal;
  else {
    if (Wexp_true > 1800)
      return kWSideband_High;
    else
      return (n_categories == 3 && Wexp_true > 1400) ? kWSideband_Mid : 
                                                       kWSideband_Low;
  }
}


std::string GetTruthClassification_LegendLabel(const WSidebandType w_category, 
                                               const int n_categories = 3) {
  if (n_categories != 2 && n_categories != 3)
    throw std::invalid_argument("GetWSidebandType: n_categories is 2 or 3");

  switch (w_category) {
    case kWSideband_Low:
    {
      return n_categories == 3 ? "W_{exp} True < 1.4" :
                                 "W_{exp} True < 1.8";
    }
    case kWSideband_Mid:
      return "W_{exp} True < 1.8";
    case kWSideband_High:
      return "W_{exp} True > 1.8";
    case kWSideband_Signal:
      return "Signal";
    default:
      return "ERROR";
  }
}


std::string GetTruthClassification_Name(WSidebandType w_category) {
  switch (w_category) {
    case kWSideband_Low:
      return "WexpTrueLow";
    case kWSideband_Mid:
      return "WexpTrueMid";
    case kWSideband_High:
      return "WexpTrueHigh";
    case kWSideband_Signal:
      return "Signal";
    default:
      return "GetTruthClassification_Name: ERROR, invalid WSidebandType";
  }
}


//==============================================================================

enum Pi0SidebandType
{
  kPi0Sideband_Signal, kPi0Sideband_Pi0, kPi0Sideband_Other, kNPi0SidebandTypes
};

Pi0SidebandType GetPi0SidebandType(const CVUniverse& universe, 
                                   const SignalDefinition signal_definition) {
  int n_pi0 = universe.GetInt("truth_N_pi0");
  if (IsSignal(universe, signal_definition)) return kPi0Sideband_Signal;
  else if (n_pi0 > 0) return kPi0Sideband_Pi0;
  else return kPi0Sideband_Other;
}

std::string GetTruthClassification_LegendLabel(Pi0SidebandType pi0_category) {
  switch (pi0_category) {
    case kPi0Sideband_Pi0:
      return "Has #pi^{0}";
    case kPi0Sideband_Other:
      return "Other Background";
    case kPi0Sideband_Signal:
      return "Signal";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(Pi0SidebandType pi0_category) {
  switch (pi0_category) {
    case kPi0Sideband_Pi0:
      return "HasPi0";
    case kPi0Sideband_Other:
      return "OtherBG";
    case kPi0Sideband_Signal:
      return "Signal";
    default:
      return "ERROR";
  }
}

#endif

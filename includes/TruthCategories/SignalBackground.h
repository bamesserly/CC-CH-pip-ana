#ifndef SignalBackground_H
#define SignalBackground_H

#include "W.h"
#include "../CVUniverse.h"
#include "../SignalDefinition.h"

enum SignalBackgroundType
{
  kB, kS, kNSignalBackgroundTypes
};

enum WBackgroundType
{
  kB_HighW, kB_HighW_Npi, kB_Other, kS_W, kNWBackgroundTypes
};

enum MesonBackgroundType
{
  kB_Meson, kB_NonMeson, kS_Meson, kNMesonBackgroundTypes
};

//==============================================================================
// Simple Signal-Background (no BG breakdown)
//==============================================================================
SignalBackgroundType GetSignalBackgroundType(const CVUniverse& universe,
                                             SignalDefinition signal_definition) {
  if (IsSignal(universe, signal_definition))
    return kS;
  else 
    return kB;
}

std::string GetTruthClassification_LegendLabel(SignalBackgroundType category) {
  switch (category) {
    case kS:
      return "Signal";
    case kB:
      return "Background";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(SignalBackgroundType category) {
  switch (category) {
    case kS:
      return "Signal";
    case kB:
      return "BG";
    default:
      return "ERROR";
  }
}


//==============================================================================
// Signal-Background -- W
//==============================================================================
WBackgroundType GetWBackgroundType(const CVUniverse& universe,
                                   SignalDefinition signal_definition) {
  if (IsSignal(universe, signal_definition))
    return kS_W;
  else if ( GetTruthWType(universe, signal_definition) == kHighW ) {
    if (universe.GetInt("truth_N_pip") > 1) {
      //if (universe.shortName() == "cv") {
      //  int link_size = 200;
      //  char link [link_size];
      //  universe.arachneLink(link, link_size);
      //  std::cout << link << std::endl;
      //}
      return kB_HighW_Npi;
    }
    else
      return kB_HighW;
  }
  else 
    return kB_Other;
}

std::string GetTruthClassification_LegendLabel(WBackgroundType category) {
  switch (category) {
    case kS_W:
      return "Signal";
    case kB_HighW:
      return "High W_{exp} BG, <= 1#pi^{+}";
    case kB_HighW_Npi:
      return "High W_{exp} BG, > 1#pi^{+}";
    case kB_Other:
      return "Other BG";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(WBackgroundType category) {
  switch (category) {
    case kS_W:
      return "Signal";
    case kB_HighW:
      return "HighWBG_1pi";
    case kB_HighW_Npi:
      return "HighWBG_Npi";
    case kB_Other:
      return "OtherBG";
    default:
      return "ERROR";
  }
}
//==============================================================================


//==============================================================================
// Signal-Background -- Meson content
//==============================================================================
MesonBackgroundType GetMesonBackgroundType(const CVUniverse& universe,
                                           SignalDefinition signal_definition) {
  if (IsSignal(universe, signal_definition))
    return kS_Meson;

  std::vector<int> fs_particles = universe.GetVec<int>("mc_FSPartPDG");
  bool has_meson = false;
  // loop fs particles, search for mesons
  for (auto particle : fs_particles) {
    if ( abs(particle) ==  321  || 
         abs(particle) ==  130  ||
         abs(particle) ==  311  ||
         abs(particle) ==  310  ||
         abs(particle) ==  3222 ||
         abs(particle) ==  3122 ||
         abs(particle) ==  111  || 
         particle      == -211 )
    {
      has_meson = true;
    }
  }

  if (has_meson) return kB_Meson;
  else return kB_NonMeson;
}

std::string GetTruthClassification_LegendLabel(MesonBackgroundType category) {
  switch (category) {
    case kS_Meson:
      return "Signal";
    case kB_Meson:
      return "BG with FS meson";
    case kB_NonMeson:
      return "BG with no FS meson";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(MesonBackgroundType category) {
  switch (category) {
    case kS_Meson:
      return "Signal";
    case kB_Meson:
      return "BG_Meson";
    case kB_NonMeson:
      return "BG_NonMeson";
    default:
      return "ERROR";
  }
}



#endif

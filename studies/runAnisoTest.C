//==============================================================================
// This script visualizes diagnostic variables and the variables that we cut on
// before and after the cut is performed.
// Plots are broken down into their truth categories (includes/TruthCategories).
//==============================================================================
#ifndef runAnisoTest_C
#define runAnisoTest_C

#include "includes/CCPiEvent.h"
#include "includes/CCPiMacroUtil.h"
#include "includes/HadronVariable.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"  // GetVar
#include "includes/common_stuff.h"      // EDataMCTruth
#include "plotting_functions.h"
#include "xsec/makeCrossSectionMCInputs.C"  // GetAnalysisVariables

class Variable;
class HadronVariable;

//==============================================================================
// Loop
//==============================================================================
void Loop(const CCPiMacroUtil& util, CVUniverse* universe,
          const EDataMCTruth& type) {
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    std::cout << universe->GetVecElem("truth_genie_wgt_Theta_Delta2Npi", 2)
              << "  "
              << universe->GetVecElem("truth_genie_wgt_Theta_Delta2Npi", 4)
              << "\n";
  }  // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runAnisoTest(int signal_definition_int = 0, const char* plist = "ME1B") {
  // INIT MACRO UTILITY OBJECT
  const std::string macro("runAnisoTest");
  bool do_data = true, do_mc = true, do_truth = false;
  bool do_systematics = false, do_grid = false;
  CCPiMacroUtil util(signal_definition_int, plist, do_data, do_mc, do_truth,
                     do_systematics, do_grid);
  util.PrintMacroConfiguration(macro);
  Loop(util, util.m_error_bands.at("cv").at(0), kMC);
}

#endif  // runAnisoTest.C

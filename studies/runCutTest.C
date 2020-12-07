#ifndef runCutTest_C
#define runCutTest_C

#include "includes/CCPiEvent.h"
#include "includes/CCPiMacroUtil.h"
#include "includes/CVUniverse.h"
#include "includes/Cuts.h"
#include "includes/common_stuff.h"  // typedefs EventCount

//==============================================================================
// Loop and fill
//==============================================================================
void Loop(const CCPiMacroUtil& util, CVUniverse* universe,
          const EDataMCTruth& type,
          std::pair<EventCount*, EventCount*>& counters) {
  double counter = 0;
  double counter_2 = 0;
  double counter_3 = 0;

  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  std::cout << is_mc << "  " << is_truth << "\n";
  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    // New IsWSideband
    std::vector<int> cands;
    bool is_w_sb = false;
    PassesCuts(*universe, cands, is_mc, util.m_signal_definition,
               is_w_sb);  // set is_w_sb
    if (is_w_sb) counter += event.m_weight;

    // Old IsWSideband
    if (IsWSideband(event)) counter_2 += event.m_weight;

    //// Old PassesCuts
    // std::vector<int> cands_2;
    // if (PassesCuts(*universe, cands_2, is_mc, util.m_signal_definition))
    //  counter_2 += event.m_weight;

    //// Loop Cut-by-cut
    // MichelMap dummy1, dummy2;
    // bool pass = true;
    //// Purity and efficiency
    // for (auto i_cut : kCutsVector) {
    //  if (event.m_is_truth != IsPrecut(i_cut)) continue; // truth loop does
    //  precuts pass = pass && PassesCut(*event.m_universe, i_cut,
    //  event.m_is_mc,
    //                           event.m_signal_definition, dummy1, dummy2);
    //  if (pass && i_cut == kNode)
    //    counter_3 += event.m_weight;
    //} // cuts

  }  // events
  std::cout << "*** Done ***\n\n";

  std::cout << is_mc << "  " << is_truth << "  " << counter << "  " << counter_2
            << "\n";
}

//==============================================================================
// Main
//==============================================================================
void runCutTest(int signal_definition_int = 0, const char* plist = "ALL") {
  // INIT MACRO UTILITY OBJECT
  const std::string macro("runCutTest");
  bool do_data = true, do_mc = true, do_truth = true;
  bool do_systematics = false, do_grid = false;
  CCPiMacroUtil util(signal_definition_int, plist, do_data, do_mc, do_truth,
                     do_systematics, do_grid);
  util.PrintMacroConfiguration(macro);

  // EFFICIENCY/PURITY COUNTERS
  // typdef EventCount map<ECut, double>
  EventCount n_remaining_sig, n_remaining_bg, n_remaining_data;
  std::pair<EventCount*, EventCount*> signal_bg_counters(&n_remaining_sig,
                                                         &n_remaining_bg);
  std::pair<EventCount*, EventCount*> data_count(&n_remaining_data, NULL);

  Loop(util, util.m_data_universe, kData, data_count);
  Loop(util, util.m_error_bands.at("cv").at(0), kMC, signal_bg_counters);
  // Loop(util, util.m_error_bands_truth.at("cv").at(0), kTruth,
  // signal_bg_counters);
}

#endif

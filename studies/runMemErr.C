#ifndef runEffPurTable_C
#define runEffPurTable_C

#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/Constants.h"  // typedefs EventCount
#include "includes/Cuts.h"
#include "includes/EventSelectionTable.h"
#include "includes/MacroUtil.h"
#include "includes/Michel.h"

//==============================================================================
// Loop and fill
//==============================================================================
std::tuple<EventCount, EventCount> FillCounters(
    const CCPi::MacroUtil& util, CVUniverse* universe, const EDataMCTruth& type,
    const EventCount& s, const EventCount& b = EventCount()) {
  EventCount signal = s;
  EventCount bg = b;
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  //for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
  for (Long64_t i_event = 0; i_event < 10000; ++i_event) {
    //if(i_event==5 or i_event==4) {
      std::clock_t start;
      if (i_event % 500000 == 0)
        std::cout << (i_event / 1000) << "k " << std::endl;
      universe->SetEntry(i_event);
      trackless::MichelEvent m = trackless::GetQualityMichels(*universe);
    //}
  } // events
  std::cout << "*** Done ***\n\n";
  return {signal, bg};
}

//==============================================================================
// Main
//==============================================================================
void runMemErr(int signal_definition_int = 0, const char* plist = "ALL") {
  // INIT MACRO UTILITY OBJECT
  const std::string macro("runEffPurTable");
  bool do_data = true, do_mc = true, do_truth = true;
  bool do_systematics = false, is_grid = false;
  bool use_xrootd = false;
  std::string data_file_list = GetPlaylistFile(plist, false, use_xrootd);
  std::string mc_file_list = GetPlaylistFile(plist, true, use_xrootd);
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // EFFICIENCY/PURITY COUNTERS
  // typdef EventCount map<ECut, double>
  EventCount n_remaining_sig, n_remaining_bg, n_remaining_data;

  std::tie(n_remaining_data, std::ignore) =
      FillCounters(util, util.m_data_universe, kData, n_remaining_data);
}

#endif

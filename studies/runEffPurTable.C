#ifndef runEffPurTable_C
#define runEffPurTable_C

#include <chrono>

#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/Constants.h"  // typedefs EventCount, t12 globals for timing
#include "includes/Cuts.h"
#include "includes/EventSelectionTable.h"
#include "includes/MacroUtil.h"

//==============================================================================
// Loop and fill
//==============================================================================
std::tuple<EventCount, EventCount> FillCounters(
    const CCPi::MacroUtil& util, CVUniverse* universe, const EDataMCTruth& type,
    const EventCount& s, const EventCount& b = EventCount()) {
  EventCount signal = s;
  EventCount bg = b;
  int passMehreencut = 0;
  int passUntrOnepion = 0;
  int passUntrOnepionSignal = 0;
  int passMehreenCutmore1 = 0;
  int passMehreenCutmore1Signal = 0;
  int passMehreenCutless1 = 0;
  int passMehreenCutless1Signal = 0;

  bool is_mc, is_truth;
  //  if (type == kTruth) is_truth = true;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 100000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    // if (i_event == 100000) break;
    universe->SetEntry(i_event);
    universe->SetTruth(is_truth);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    std::map<ECuts, bool> passMap;
    if (!is_truth) {
      LowRecoilPion::Cluster d;
      LowRecoilPion::Cluster c(*universe, 0);
      LowRecoilPion::Michel<CVUniverse> m(*universe, 0);
      typedef LowRecoilPion::MichelEvent<CVUniverse> MichelEvent;
      typedef LowRecoilPion::hasMichel<CVUniverse, MichelEvent> hasMichel;
      typedef LowRecoilPion::BestMichelDistance2D<CVUniverse, MichelEvent>
          BestMichelDistance2D;
      typedef LowRecoilPion::GetClosestMichel<CVUniverse, MichelEvent>
          GetClosestMichel;
      // Get Quality Michels
      MichelEvent trackless_michels;
      //    bool pass = hasMichel::hasMichelCut(*universe, trackless_michels);
      bool pass = true;
      pass = pass && universe->GetNMichels() == 1;
      passMap.insert(std::make_pair(kOneMichel, pass));
      pass = pass && hasMichel::hasMichelCut(*universe, trackless_michels);
      passMap.insert(std::make_pair(kHasMichel, pass));
      pass = pass && BestMichelDistance2D::BestMichelDistance2DCut(
                         *universe, trackless_michels);
      passMap.insert(std::make_pair(kBestMichelDistance, pass));
      pass = pass && GetClosestMichel::GetClosestMichelCut(*universe,
                                                           trackless_michels);
      passMap.insert(std::make_pair(kClosestMichel, pass));
      universe->SetVtxMichels(trackless_michels);
      if (pass && is_mc) {
        passMehreencut++;
        if (universe->GetNMichels() > 1) {
          passMehreenCutmore1++;
          if (event.m_is_signal) passMehreenCutmore1Signal++;
        }
        if (universe->GetNMichels() < 1) {
          passMehreenCutless1++;
          if (event.m_is_signal) passMehreenCutless1Signal++;
        }
        if (universe->GetNMichels() == 1) {
          passUntrOnepion++;
          if (event.m_is_signal) passUntrOnepionSignal++;
        }
      }
      //    pass = pass && universe->GetNMichels() == 1;
      //    passMap.insert(std::make_pair(kOneMichel, pass));
      pass = pass &&
             universe->GetTpiTrackless() > util.m_signal_definition.m_tpi_min;
      pass = pass &&
             universe->GetTpiTrackless() < util.m_signal_definition.m_tpi_max;
      passMap.insert(std::make_pair(kTpi, pass));
      pass = pass && universe->GetPTmu() < util.m_signal_definition.m_ptmu_max;
      passMap.insert(std::make_pair(kPTmu, pass));
      pass = pass && universe->GetTracklessWexp() > 0.;
      pass = pass && universe->GetTracklessWexp() < 1400;
      passMap.insert(std::make_pair(kUntrackedWexp, pass));

      //      if (pass)
      //        std::cout << "Event = " << i_event << "\n";
      //        for (auto i_cut : kUntrackedCutsVector)
      // 	std::cout << "Cut " << i_cut << " " << passMap[i_cut] << "\n";
    }
    std::tie(signal, bg) = ccpi_event::FillCounters(
        event, signal, bg, passMap);  // Does a lot of work
  }                                   // events

  if (is_mc) {
    std::cout << "Pass Mehreen cuts = " << passMehreencut << "\n";
    std::cout << "Pass Mehreen cuts N michels > 1 = " << passMehreenCutmore1
              << "\n";
    std::cout << "Pass Mehreen cuts N michels > 1 signal = "
              << passMehreenCutmore1Signal << "\n";
    std::cout << "Pass Mehreen cuts N michels < 1 = " << passMehreenCutless1
              << "\n";
    std::cout << "Pass Mehreen cuts N michels < 1 Signal = "
              << passMehreenCutless1Signal << "\n";
    std::cout << "Pass Mehreen cuts N michels = 1 " << passUntrOnepion << "\n";
    std::cout << "Pass Mehreen cuts N michels = 1 Signal = "
              << passUntrOnepionSignal << "\n";
  }

  std::cout << "*** Done ***\n\n";
  return {signal, bg};
}

//==============================================================================
// Main
//==============================================================================
void runEffPurTable(int signal_definition_int = 1, const char* plist = "ALL") {
  auto start = std::chrono::steady_clock::now();
  bool is_mc = true;
  const bool use_xrootd = true;
  const bool do_test_playlist = false;
  std::string mc_file_list = GetPlaylistFile(plist, is_mc, use_xrootd);
  is_mc = false;
  std::string data_file_list = GetPlaylistFile(plist, is_mc, use_xrootd);
  // const int signal_definition_int = SignalDefinition::OnePiTracked().m_id;
  const bool is_grid = false;
  const bool do_truth = true;
  const bool do_systematics = false;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.m_name = "runEffPurTable";
  util.PrintMacroConfiguration();

  // EFFICIENCY/PURITY COUNTERS
  // typdef EventCount map<ECut, double>
  EventCount n_remaining_sig, n_remaining_bg, n_remaining_data;

  std::tie(n_remaining_data, std::ignore) =
      FillCounters(util, util.m_data_universe, kData, n_remaining_data);

  std::tie(n_remaining_sig, n_remaining_bg) =
      FillCounters(util, util.m_error_bands.at("cv").at(0), kMC,
                   n_remaining_sig, n_remaining_bg);

  std::tie(n_remaining_sig, n_remaining_bg) =
      FillCounters(util, util.m_error_bands_truth.at("cv").at(0), kTruth,
                   n_remaining_sig, n_remaining_bg);
  std::cout << "n_remaining_sig.at(kNoCuts) " << n_remaining_sig[kNoCuts]
            << "\n";

  PrintEffPurTable(n_remaining_sig, n_remaining_bg, n_remaining_data,
                   util.m_data_pot, util.m_mc_pot);
  auto stop = std::chrono::steady_clock::now();
  std::cout
      << "running time: "
      << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()
      << "sec\n";
  std::cout << "cluster creation time: " << t01.count() << "sec\n";
  std::cout << "cluster loop 1: " << t12.count() << "sec\n";
  std::cout << "cluster loop 2: " << t23.count() << "sec\n";
}

#endif

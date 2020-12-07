// This script started from a clone of the runEffPurTable script and it still
// produces that table. Additionally, I copy-pasted CVUniverse::GetWeight here
// and plot individual weight components.

// The purpose of this script was to validate the latest PlotUtils, updated
// from my thesis tag.

#ifndef runEffPurTable_C
#define runEffPurTable_C

#include "includes/CCPiMacroUtil.h"
#include "includes/CVUniverse.h"
#include "includes/common_stuff.h" // typedefs EventCount
#include "includes/CCPiEvent.h"
#include "includes/Cuts.h"
#include "../event_selection/EventSelectionTable.h"
#include "plotting_functions.h"


//==============================================================================
// Loop and fill
//==============================================================================
void FillCounters(const CCPiMacroUtil& util, CVUniverse* universe,
                  const EDataMCTruth& type,
                  std::pair<EventCount*, EventCount*>& counters) {

  // STUDY: weights and plotutils versions
  TH1D* h_total_wgts       = new TH1D(Form("total_weights_%i",       type), "total_weights",       50, -1., 5.);
  TH1D* h_total_check_wgts = new TH1D(Form("total_check_weights_%i", type), "total_check_weights", 50, -1., 5.);
  TH1D* h_flux_wgts        = new TH1D(Form("flux_weights_%i",        type), "flux_weights",        50, -1., 5.);
  TH1D* h_genie_wgts       = new TH1D(Form("genie_weights_%i",       type), "genie_weights",       50, -1., 5.);
  TH1D* h_rpa_wgts         = new TH1D(Form("rpa_weights_%i",         type), "rpa_weights",         50, -1., 5.);
  TH1D* h_nrp_wgts         = new TH1D(Form("nrp_weights_%i",         type), "nrp_weights",         50, -1., 5.);
  TH1D* h_2p2h_wgts        = new TH1D(Form("2p2h_weights_%i",        type), "2p2h_weights",        50, -1., 5.);
  TH1D* h_mueff_wgts       = new TH1D(Form("mueff_weights_%i",       type), "mueff_weights",       50, -1., 5.);

  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  for(Long64_t i_event=0; i_event < n_entries; ++i_event){
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    // STUDY: weights and plotutils versions
    if (type == kMC) {

      double wgt_flux=1., wgt_2p2h=1.;
      double wgt_rpa=1.,   wgt_nrp=1.;
      double wgt_genie=1., wgt_mueff=1.;

      double Enu  = universe->GetDouble("mc_incomingE")*1e-3;
      int nu_type = universe->GetInt("mc_incoming");

      wgt_flux = universe->GetFluxAndCVWeight(Enu, nu_type);
      wgt_genie = universe->GetGenieWeight();
      wgt_rpa = universe->GetRPAWeight();
      //wgt_2p2h = universe->Get2p2hWeight(q0, q3);
      //wgt_nrp = universe->GetNonResPiWeight();
      if (!universe->m_is_truth && universe->GetBool("isMinosMatchTrack"))
        wgt_mueff = universe->GetMinosEfficiencyWeight();

      double wgt_total_check = wgt_genie*wgt_flux*wgt_2p2h*wgt_rpa*wgt_nrp*wgt_mueff;
      //if(wgt_total_check != event.m_weight)
      //  std::cout << wgt_total_check/event.m_weight << "\n";

      //std::cout << wgt_flux << "\n";

      h_total_wgts->Fill(event.m_weight);
      h_flux_wgts->Fill(wgt_flux);
      h_genie_wgts->Fill(wgt_genie);
      h_rpa_wgts->Fill(wgt_rpa);
      h_nrp_wgts->Fill(wgt_nrp);
      h_2p2h_wgts->Fill(wgt_2p2h);
      h_mueff_wgts->Fill(wgt_mueff);
      h_total_check_wgts->Fill(wgt_total_check);

    }

    ccpi_event::FillCounters(event, counters); // Does a lot of work
  } // events
  std::cout << "*** Done ***\n\n";

  const bool use_log_scale = true;
  const char* tag;
  if(type == kMC) {
    PlotTH1_1(h_total_wgts,       Form("wgts_total_%s",       tag), -1, use_log_scale);
    PlotTH1_1(h_flux_wgts,        Form("wgts_flux_%s",        tag), -1, use_log_scale);
    PlotTH1_1(h_genie_wgts,       Form("wgts_genie_%s",       tag), -1, use_log_scale);
    PlotTH1_1(h_rpa_wgts,         Form("wgts_rpa_%s",         tag), -1, use_log_scale);
    PlotTH1_1(h_nrp_wgts,         Form("wgts_nrp_%s",         tag), -1, use_log_scale);
    PlotTH1_1(h_2p2h_wgts,        Form("wgts_2p2h_%s",        tag), -1, use_log_scale);
    PlotTH1_1(h_mueff_wgts,       Form("wgts_mueff_%s",       tag), -1, use_log_scale);
    PlotTH1_1(h_total_check_wgts, Form("wgts_total_check_%s", tag), -1, use_log_scale);
  }

}


//==============================================================================
// Main
//==============================================================================
void runWeightCheck(int signal_definition_int = 0, const char* plist = "ALL") {
  // INIT MACRO UTILITY OBJECT
  const std::string macro("runEffPurTable");
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

  FillCounters(util, util.m_data_universe,  kData,  data_count);
  FillCounters(util, util.m_error_bands.at("cv").at(0), kMC,
               signal_bg_counters);
  FillCounters(util, util.m_error_bands_truth.at("cv").at(0), kTruth,
               signal_bg_counters);

  PrintEffPurTable(n_remaining_sig, n_remaining_bg, n_remaining_data,
                   util.m_data_pot, util.m_mc_pot);
}

#endif

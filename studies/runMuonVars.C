//==============================================================================
// Template for a generic loop data/mc and plot script.
// Assume that NO systematics are being analyzed.
// Good for stacked histograms, branch validation, residuals, etc.
//==============================================================================
#include <iostream>
#include <vector>

#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/HadronVariable.h"
#include "includes/MacroUtil.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "xsec/plotting_functions.h"

// Forward declare my variables because we're hiding the header.
class Variable;
class HadronVariable;

// These functions get the muon reco px and py in the DETECTOR coordinates.
// We want them wrt to the BEAM angle (3 deg rotation).
// And further these quantities won't propagate systematics.
// Thus we use the mat official functions.
double GetPXmuMAD(CVUniverse* u) {
  return u->GetDouble("MasterAnaDev_muon_Px");
}
double GetPYmuMAD(CVUniverse* u) {
  return u->GetDouble("MasterAnaDev_muon_Py");
}
double GetPZmuMAD(CVUniverse* u) {
  return u->GetDouble("MasterAnaDev_muon_Pz");
}
double GetPmuMAD(CVUniverse* u) {
  return sqrt(pow(GetPXmuMAD(u), 2.0) + pow(GetPYmuMAD(u), 2.0) +
              pow(GetPZmuMAD(u), 2.0));
}
double GetThetamuMAD(CVUniverse* u) {
  return u->GetDouble("MasterAnaDev_muon_theta");
}

namespace run_muon_vars {
//==============================================================================
// Do some event processing (e.g. make cuts, get best pion) and fill hists
//==============================================================================
void FillVars(CCPiEvent& event, const std::vector<Variable*>& variables) {
  const CVUniverse* universe = event.m_universe;
  const double wgt = event.m_weight;
  const bool is_mc = event.m_is_mc;
  const SignalDefinition sd = event.m_signal_definition;

  if (universe->ShortName() != "cv") return;

  event.m_passes_cuts = PassesCuts(event, event.m_is_w_sideband);
  event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
  if (event.m_passes_cuts) ccpi_event::FillStackedHists(event, variables);
}

//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  typedef Variable Var;
  typedef HadronVariable HVar;
  HVar* thetapi_deg =
      new HVar("thetapi_deg", "#theta_{#pi}", "deg",
               CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);
  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);
  std::vector<Var*> variables = {thetapi_deg, pmu};
  return variables;
}
}  // namespace run_muon_vars

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFill(const CCPi::MacroUtil& util, CVUniverse* universe,
                 const EDataMCTruth& type, std::vector<Variable*>& variables) {
  TH1D* h_pxmu_mad = new TH1D("pxmu_mad", "pxmu_mad", 200, -1000, 1000.);
  TH1D* h_pxmu_new = new TH1D("pxmu_new", "pxmu_new", 200, -1000, 1000.);
  TH1D* h_pxmu_mad_lep =
      new TH1D("pxmu_mad_lep", "pxmu_mad_lep", 200, -1000, 1000.);
  TH1D* h_pxmu_resid = new TH1D("pxmu_resid", "pxmu_resid", 100, -1., 1.);

  TH1D* h_pymu_mad = new TH1D("pymu_mad", "pymu_mad", 200, -1000, 1000.);
  TH1D* h_pymu_new = new TH1D("pymu_new", "pymu_new", 200, -1000, 1000.);
  TH1D* h_pymu_mad_lep =
      new TH1D("pymu_mad_lep", "pymu_mad_lep", 200, -1000, 1000.);
  TH1D* h_pymu_resid = new TH1D("pymu_resid", "pymu_resid", 100, -5., 5.);

  TH1D* h_pzmu_mad = new TH1D("pzmu_mad", "pzmu_mad", 200, 0., 2000.);
  TH1D* h_pzmu_new = new TH1D("pzmu_new", "pzmu_new", 200, 0., 2000.);
  TH1D* h_pzmu_mad_lep =
      new TH1D("pzmu_mad_lep", "pzmu_mad_lep", 200, 0., 2000.);
  TH1D* h_pzmu_resid = new TH1D("pzmu_resid", "pzmu_resid", 100, -1., 1.);

  TH1D* h_pmu_mad = new TH1D("pmu_mad", "pmu_mad", 200, 0., 2000.);
  TH1D* h_pmu_new = new TH1D("pmu_new", "pmu_new", 200, 0., 2000.);
  TH1D* h_pmunom_new = new TH1D("pmunom_new", "pmunom_new", 200, 0., 2000.);
  TH1D* h_pmu_resid = new TH1D("pmu_resid", "pmu_resid", 100, -1., 1.);

  TH1D* h_thmu_mad = new TH1D("thmu_mad", "thmu_mad", 200, 0., 2.);
  TH1D* h_thmu_new = new TH1D("thmu_new", "thmu_new", 200, 0., 2.);
  TH1D* h_thmu_resid = new TH1D("thmu_resid", "thmu_resid", 100, -1., 1.);

  std::cout << "Loop and Fill CutVars\n";
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    // for (Long64_t i_event = 0; i_event < 200; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);

    // For mc, get weight, check signal, and sideband
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    h_pxmu_new->Fill(universe->GetPXmu());
    h_pxmu_mad->Fill(GetPXmuMAD(universe));
    h_pxmu_mad_lep->Fill(universe->GetVecElem("MasterAnaDev_leptonE", 0));
    h_pxmu_resid->Fill(GetPXmuMAD(universe) / universe->GetPXmu() - 1);

    h_pymu_new->Fill(universe->GetPYmu());
    h_pymu_mad->Fill(GetPYmuMAD(universe));
    h_pymu_mad_lep->Fill(universe->GetVecElem("MasterAnaDev_leptonE", 1));
    h_pymu_resid->Fill(GetPYmuMAD(universe) / universe->GetPYmu() - 1);

    h_pzmu_new->Fill(universe->GetPZmu());
    h_pzmu_mad->Fill(GetPZmuMAD(universe));
    h_pzmu_mad_lep->Fill(universe->GetVecElem("MasterAnaDev_leptonE", 2));
    h_pzmu_resid->Fill(GetPZmuMAD(universe) / universe->GetPZmu() - 1);

    h_pmu_mad->Fill(GetPmuMAD(universe));
    h_pmu_new->Fill(universe->GetPmu());
    h_pmunom_new->Fill(universe->GetPmu_nominal());
    h_pmu_resid->Fill(GetPmuMAD(universe) / universe->GetPmu() - 1);

    h_thmu_mad->Fill(GetThetamuMAD(universe));
    h_thmu_new->Fill(universe->GetThetamu());
    h_thmu_resid->Fill(GetThetamuMAD(universe) / universe->GetThetamu() - 1);

    //// WRITE THE FILL FUNCTION
    // run_muon_vars::FillVars(event, variables);

    // std::cout << universe->GetPXmu()    << " ";
    // std::cout << universe->GetPXmuMAD() << " ,";
    // std::cout << universe->GetPYmu()    << " ";
    // std::cout << universe->GetPYmuMAD() << "\n";

  }  // events
  std::cout << "*** Done ***\n\n";

  // PlotTogether(h_pxmu_mad, "pxmu_mad", h_pxmu_new, "pxmu_heidi",
  // Form("pxmu%i", is_mc)); PlotTogether(h_pymu_mad, "pymu_mad", h_pymu_new,
  // "pymu_heidi", Form("pymu%i", is_mc)); PlotTogether(h_pzmu_mad, "pzmu_mad",
  // h_pzmu_new, "pzmu_heidi", Form("pzmu%i", is_mc));

  PlotTH1_1(h_pmu_mad, Form("pmu_mad%i", is_mc));
  PlotTH1_1(h_pmu_new, Form("pmu_new%i", is_mc));
  PlotTH1_1(h_pmunom_new, Form("pmu_mad_lep%i", is_mc));
  PlotTH1_1(h_pmu_resid, Form("pmuresid%i", is_mc), -1, true);

  PlotTH1_1(h_pxmu_mad, Form("pxmu_mad%i", is_mc));
  PlotTH1_1(h_pxmu_new, Form("pxmu_new%i", is_mc));
  PlotTH1_1(h_pxmu_mad_lep, Form("pxmu_mad_lep%i", is_mc));
  PlotTH1_1(h_pxmu_resid, Form("pxmuresid%i", is_mc), -1, true);

  PlotTH1_1(h_pymu_mad, Form("pymu_mad%i", is_mc));
  PlotTH1_1(h_pymu_new, Form("pymu_new%i", is_mc));
  PlotTH1_1(h_pymu_mad_lep, Form("pymu_mad_lep%i", is_mc));
  PlotTH1_1(h_pymu_resid, Form("pymuresid%i", is_mc));

  PlotTH1_1(h_pzmu_mad, Form("pzmu_mad%i", is_mc));
  PlotTH1_1(h_pzmu_new, Form("pzmu_new%i", is_mc));
  PlotTH1_1(h_pzmu_mad_lep, Form("pzmu_mad_lep%i", is_mc));
  PlotTH1_1(h_pzmu_resid, Form("pzmuresid%i", is_mc), -1, true);

  PlotTH1_1(h_thmu_mad, Form("thmu_mad%i", is_mc));
  PlotTH1_1(h_thmu_new, Form("thmu_new%i", is_mc));
  PlotTH1_1(h_thmu_resid, Form("thmuresid%i", is_mc), -1, true);
}

//==============================================================================
// Main
//==============================================================================
void runMuonVars(std::string plist = "ME1L") {
  //=========================================
  // Input tuples
  //=========================================
  bool is_mc = true;
  std::string mc_file_list, data_file_list;
  mc_file_list = GetPlaylistFile(plist, is_mc);
  is_mc = false;
  data_file_list = GetPlaylistFile(plist, is_mc);

  //=========================================
  // Init macro utility
  //=========================================
  const int signal_definition_int = 0;
  const std::string macro("runMuonVars");
  const bool is_grid = false;
  const bool do_truth = false;
  const bool do_systematics = false;

  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  //=========================================
  // Get variables and initialize their hists
  //=========================================
  std::vector<Variable*> variables = run_muon_vars::GetVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  //=========================================
  // Loop and Fill
  //=========================================
  LoopAndFill(util, util.m_data_universe, kData, variables);
  LoopAndFill(util, util.m_error_bands.at("cv").at(0), kMC, variables);

  // for (auto v : variables) {
  //  std::string tag = v->Name();
  //  double ymax = -1;
  //  bool do_bwn = true;
  //  std::cout << "Plotting" << std::endl;
  //  PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kS),
  //             util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
  //             v->Name(), "SSB", ymax, do_bwn);
  //}

  std::cout << "Success" << std::endl;
}

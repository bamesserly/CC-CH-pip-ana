#ifndef makeXsecMCInputs_C
#define makeXsecMCInputs_C

#include <cassert>

#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/Cuts.h"
#include "includes/SignalDefinition.h"
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/common_functions.h"           // GetVar, WritePOT
#include "includes/common_stuff.h"               // typedef EventCount
#include "includes/MacroUtil.h"

#ifndef __CINT__            // CINT doesn't know about std::function
#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/HadronVariable.h"
#include "includes/Variable.h"
#endif

class Variable;
class HadronVariable;

//==============================================================================
// Helper Functions
//==============================================================================
namespace make_xsec_mc_inputs {
typedef Variable Var;
typedef HadronVariable HVar;

std::vector<Variable*> GetOnePiVariables(bool include_truth_vars = true) {
  HVar* tpi = new HVar("tpi", "T_{#pi}", "MeV", GetBinArray("tpi"),
                       &CVUniverse::GetTpi);

  HVar* tpi_mbr = new HVar("tpi_mbr", "T_{#pi} (MBR)", tpi->m_units,
                           GetBinArray("tpi"), &CVUniverse::GetTpiMBR);

  HVar* thetapi_deg =
      new HVar("thetapi_deg", "#theta_{#pi}", "deg", GetBinArray("thetapi_deg"),
               &CVUniverse::GetThetapiDeg);

  Var* pmu =
      new Var("pmu", "p_{#mu}", "MeV", GetBinArray("pmu"), &CVUniverse::GetPmu);

  Var* thetamu_deg =
      new Var("thetamu_deg", "#theta_{#mu}", "deg", GetBinArray("thetamu_deg"),
              &CVUniverse::GetThetamuDeg);

  Var* enu =
      new Var("enu", "E_{#nu}", "MeV", GetBinArray("enu"), &CVUniverse::GetEnu);

  Var* q2 =
      new Var("q2", "Q^{2}", "MeV^{2}", GetBinArray("q2"), &CVUniverse::GetQ2);

  Var* wexp = new Var("wexp", "W_{exp}", "MeV", GetBinArray("wexp"),
                      &CVUniverse::GetWexp);

  Var* wexp_fit = new Var(sidebands::kFitVarString, wexp->m_hists.m_xlabel,
                          wexp->m_units, 32, 0.e3, 3.2e3, &CVUniverse::GetWexp);

  Var* ptmu = new Var("ptmu", "p^{T}_{#mu}", "MeV", GetBinArray("ptmu"),
                      &CVUniverse::GetPTmu);

  Var* pzmu = new Var("pzmu", "p^{z}_{#mu}", "MeV", GetBinArray("pzmu"),
                      &CVUniverse::GetPZmu);

  // True Variables
  bool is_true = true;
  HVar* tpi_true =
      new HVar("tpi_true", "T_{#pi} True", tpi->m_units,
               tpi->m_hists.m_bins_array, &CVUniverse::GetTpiTrue, is_true);

  HVar* thetapi_deg_true =
      new HVar("thetapi_deg_true", "#theta_{#pi} True", thetapi_deg->m_units,
               thetapi_deg->m_hists.m_bins_array,
               &CVUniverse::GetThetapiTrueDeg, is_true);

  Var* pmu_true =
      new Var("pmu_true", "p_{#mu} True", pmu->m_units,
              pmu->m_hists.m_bins_array, &CVUniverse::GetPmuTrue, is_true);

  Var* thetamu_deg_true =
      new Var("thetamu_deg_true", "#theta_{#mu} True", thetamu_deg->m_units,
              thetamu_deg->m_hists.m_bins_array, &CVUniverse::GetThetamuTrueDeg,
              is_true);

  Var* enu_true =
      new Var("enu_true", "E_{#nu} True", enu->m_units,
              enu->m_hists.m_bins_array, &CVUniverse::GetEnuTrue, is_true);

  Var* q2_true =
      new Var("q2_true", "Q^{2} True", q2->m_units, q2->m_hists.m_bins_array,
              &CVUniverse::GetQ2True, is_true);

  Var* wexp_true =
      new Var("wexp_true", "W_{exp} True", wexp->m_units,
              wexp->m_hists.m_bins_array, &CVUniverse::GetWexpTrue, is_true);

  Var* ptmu_true =
      new Var("ptmu_true", "pt_{#mu} True", "MeV", ptmu->m_hists.m_bins_array,
              &CVUniverse::GetPTmuTrue, is_true);

  Var* pzmu_true =
      new Var("pzmu_true", "pz_{#mu} True", "MeV", pzmu->m_hists.m_bins_array,
              &CVUniverse::GetPZmuTrue, is_true);

  std::vector<Var*> variables = {tpi,         tpi_mbr, thetapi_deg, pmu,
                                 thetamu_deg, enu,     q2,          wexp,
                                 wexp_fit,    ptmu,    pzmu};

  if (include_truth_vars) {
    variables.push_back(tpi_true);
    variables.push_back(thetapi_deg_true);
    variables.push_back(pmu_true);
    variables.push_back(thetamu_deg_true);
    variables.push_back(enu_true);
    variables.push_back(q2_true);
    variables.push_back(wexp_true);
    variables.push_back(ptmu_true);
    variables.push_back(pzmu_true);
  }

  return variables;
}

std::map<std::string, Variable*> GetOnePiVariables_Map(
    bool include_truth_vars = true) {
  std::map<std::string, Var*> var_map;
  std::vector<Var*> var_vec = GetOnePiVariables(include_truth_vars);
  for (auto v : var_vec) var_map[v->Name()] = v;
  return var_map;
}

}  // namespace make_xsec_mc_inputs

std::vector<Variable*> GetAnalysisVariables(SignalDefinition signal_definition,
                                            bool include_truth_vars = false) {
  std::vector<Variable*> variables;
  switch (signal_definition) {
    case kOnePi:
      variables = make_xsec_mc_inputs::GetOnePiVariables(include_truth_vars);
      break;
    default:
      std::cerr << "Variables for other SDs not yet implemented.\n";
      std::exit(1);
  }
  return variables;
}

void SyncAllHists(Variable& v) {
  v.m_hists.m_selection_mc.SyncCVHistos();
  v.m_hists.m_bg.SyncCVHistos();
  v.m_hists.m_bg_loW.SyncCVHistos();
  v.m_hists.m_bg_midW.SyncCVHistos();
  v.m_hists.m_bg_hiW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_sig.SyncCVHistos();
  v.m_hists.m_wsidebandfit_loW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_midW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_hiW.SyncCVHistos();
  v.m_hists.m_effnum.SyncCVHistos();
  v.m_hists.m_effden.SyncCVHistos();
}

//==============================================================================
// Loop and Fill
//==============================================================================
#ifndef __CINT__
void LoopAndFillMCXSecInputs(const CCPi::MacroUtil& util,
                             const EDataMCTruth& type,
                             std::vector<Variable*>& variables) {
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  const UniverseMap error_bands =
      is_truth ? util.m_error_bands_truth : util.m_error_bands;
  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;

    // Variables that hold info about whether the CVU passes cuts
    bool checked_cv = false, cv_passes_cuts = false, cv_is_w_sideband = false;
    std::vector<RecoPionIdx> cv_reco_pion_candidate_idxs;

    // Loop universes, make cuts, and fill
    for (auto error_band : error_bands) {
      std::vector<CVUniverse*> universes = error_band.second;
      for (auto universe : universes) {
        universe->SetEntry(i_event);
        if (universe->GetDouble("mc_incoming") == 12 && universe->ShortName() =="cv")
          universe->PrintArachneLink();
        CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

        //===============
        // FILL TRUTH
        //===============
        if (type == kTruth) {
          ccpi_event::FillTruthEvent(event, variables);
          continue;
        }

        //===============
        // CHECK CUTS
        //===============
        if (universe->IsVerticalOnly()) {  // Universe only affects weights
          if (!checked_cv) {  // Only check vertical-only universes once.
            cv_passes_cuts =
                PassesCuts(*universe, cv_reco_pion_candidate_idxs, is_mc,
                           util.m_signal_definition, cv_is_w_sideband);
            checked_cv = true;
          }

          if (checked_cv) {  // Already checked a vertical only universe
            event.m_passes_cuts = cv_passes_cuts;
            event.m_is_w_sideband = cv_is_w_sideband;
            event.m_reco_pion_candidate_idxs = cv_reco_pion_candidate_idxs;
            event.m_highest_energy_pion_idx =
                GetHighestEnergyPionCandidateIndex(event);
          }
        } else {  // Universe shifts something laterally
          event.m_passes_cuts = PassesCuts(event, event.m_is_w_sideband);
          event.m_highest_energy_pion_idx =
              GetHighestEnergyPionCandidateIndex(event);
        }

        //===============
        // FILL RECO
        //===============
        ccpi_event::FillRecoEvent(event, variables);
      }  // universes
    }    // error bands
  }      // events
  std::cout << "*** Done ***\n\n";
}
#endif

//==============================================================================
// Main
//==============================================================================
void makeCrossSectionMCInputs(int signal_definition_int = 0,
                              std::string plist = "ME1A",
                              bool do_systematics = false,
                              bool do_truth = false, bool is_grid = false,
                              std::string input_file = "", int run = 0) {

  // INPUT TUPLES
  const bool is_mc = true;
  std::string mc_file_list;
  assert(!(is_grid && input_file.empty()) && "On the grid, infile must be specified.");
  mc_file_list = input_file.empty() ? GetPlaylistFile(plist, is_mc) : input_file;
    
  // INIT MACRO UTILITY
  const std::string macro("MCXSecInputs");
  //std::string a_file = "root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr/minerva/persistent/users/bmesserl/pions//20200713/merged/mc/ME1A/CCNuPionInc_mc_AnaTuple_run00110000_Playlist.root";
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, plist, do_truth,
                       is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // INIT OUTPUT
  const std::string tag("20200913");
  std::string outfile_name(Form("%s_%d%d%d%d_%s_%d_%s.root", macro.c_str(),
                                signal_definition_int, int(do_systematics), int(do_truth),
                                int(is_grid), plist.c_str(), run, tag.c_str()));
  std::cout << "Saving output to " << outfile_name << "\n\n";
  TFile fout(outfile_name.c_str(), "RECREATE");

  // INIT VARS, HISTOS, AND EVENT COUNTERS
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  // MC RECO and TRUTH
  LoopAndFillMCXSecInputs(util, kMC, variables);
  if (util.m_do_truth) LoopAndFillMCXSecInputs(util, kTruth, variables);

  // WRITE TO FILE
  std::cout << "Synching and Writing\n\n";
  WritePOT(fout, true, util.m_mc_pot);
  fout.cd();
  for (auto v : variables) {
    SyncAllHists(*v);
    v->WriteMCHists(fout);
  }
}

#endif  // makeXsecMCInputs_C

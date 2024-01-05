#ifndef makeXsecMCInputs_C
#define makeXsecMCInputs_C

#include <cassert>
#include <ctime>
#include <functional>

#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/Constants.h"
#include "includes/Cuts.h"
#include "includes/HadronVariable.h"
#include "includes/MacroUtil.h"
#include "includes/SignalDefinition.h"
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/Variable.h"
#include "includes/common_functions.h"  // GetVar, WritePOT

//==============================================================================
// Helper Functions
//==============================================================================
namespace make_xsec_mc_inputs {
typedef Variable Var;
typedef HadronVariable HVar;

std::vector<Variable*> GetOnePiVariables(bool include_truth_vars = true) {
  const int nadphibins = 16;
  const double adphimin = -CCNuPionIncConsts::PI;
  const double adphimax = CCNuPionIncConsts::PI;

  Var* tpi = new HVar("tpi", "T_{#pi}", "MeV", CCPi::GetBinning("tpi"),
                       &CVUniverse::GetTpi);

  //Var* tpi_mbr = new HVar("tpi_mbr", "T_{#pi} (MBR)", tpi->m_units,
  //                         CCPi::GetBinning("tpi"), &CVUniverse::GetTpiMBR);

  Var* thetapi_deg =
      new HVar("thetapi_deg", "#theta_{#pi}", "deg",
               CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);

  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);

  Var* thetamu_deg =
      new Var("thetamu_deg", "#theta_{#mu}", "deg",
              CCPi::GetBinning("thetamu_deg"), &CVUniverse::GetThetamuDeg);

  Var* enu = new Var("enu", "E_{#nu}", "MeV", CCPi::GetBinning("enu"),
                     &CVUniverse::GetEnu);

  Var* q2 = new Var("q2", "Q^{2}", "MeV^{2}", CCPi::GetBinning("q2"),
                    &CVUniverse::GetQ2);

  Var* wexp = new Var("wexp", "W_{exp}", "MeV", CCPi::GetBinning("wexp"));

  Var* wexp_fit =
      new Var(sidebands::kFitVarString, wexp->m_hists.m_xlabel, wexp->m_units,
              CCPi::GetBinning("wexp_fit"), &CVUniverse::GetWexp);

  Var* ptmu = new Var("ptmu", "p^{T}_{#mu}", "MeV", CCPi::GetBinning("ptmu"),
                      &CVUniverse::GetPTmu);

  Var* pzmu = new Var("pzmu", "p^{z}_{#mu}", "MeV", CCPi::GetBinning("pzmu"),
                      &CVUniverse::GetPZmu);

  Var* ehad = new Var("ehad", "ehad", "MeV", CCPi::GetBinning("ehad"),
                      &CVUniverse::GetEhad);

  // True Variables
  bool is_true = true;
  Var* tpi_true =
      new Var("tpi_true", "T_{#pi} True", tpi->m_units,
               tpi->m_hists.m_bins_array);
  tpi_true->m_is_true = true;

  Var* thetapi_deg_true =
      new Var("thetapi_deg_true", "#theta_{#pi} True", thetapi_deg->m_units,
               thetapi_deg->m_hists.m_bins_array);
  thetapi_deg_true->m_is_true = true;

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

  Var* ehad_true =
      new Var("ehad_true", "ehad True", "MeV", ehad->m_hists.m_bins_array,
              &CVUniverse::GetEhadTrue, is_true);

  std::vector<Var*> variables = {tpi,         thetapi_deg, pmu,
                                 thetamu_deg, enu,     q2,          wexp,
                                 wexp_fit,    ptmu,    pzmu,        ehad};

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
    variables.push_back(ehad_true);
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

std::vector<Variable*> GetAnalysisVariables(
    const SignalDefinition& signal_definition,
    const bool include_truth_vars = false) {
  using GetVariablesFn = std::function<std::vector<Variable*>(bool)>;
  std::map<int, GetVariablesFn> get_variables{
      {SignalDefinition::OnePi().m_id, make_xsec_mc_inputs::GetOnePiVariables},
      {SignalDefinition::OnePiTracked().m_id,
       make_xsec_mc_inputs::GetOnePiVariables}};
  return get_variables.at(signal_definition.m_id)(include_truth_vars);
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

// Given a macro, make an output filename with a timestamp
std::string GetOutFilename(const CCPi::MacroUtil& util, const int run = 0) {
  auto time = std::time(nullptr);
  char tchar[100];
  std::strftime(tchar, sizeof(tchar), "%F", std::gmtime(&time));  // YYYY-MM-dd
  const std::string timestamp = tchar;
  const std::string outfile_format = "%s_%d%d%d%d_%s_%d_%s.root";
  return std::string(Form(outfile_format.c_str(), util.m_name.c_str(),
                          util.m_signal_definition.m_id,
                          int(util.m_do_systematics), int(util.m_do_truth),
                          int(util.m_is_grid), util.m_plist_string.c_str(), run,
                          timestamp.c_str()));
}

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillMCXSecInputs(const UniverseMap& error_bands,
                             const Long64_t n_entries, const bool is_truth,
                             const SignalDefinition& signal_definition,
                             std::vector<Variable*>& variables) {
  const bool is_mc = true;

  for (auto band : error_bands) {
    std::vector<CVUniverse*> universes = band.second;
    for (auto universe : universes) universe->SetTruth(is_truth);
  }

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % (n_entries / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;

    // Save vertical-only universe info across universes for optimization --
    // there's no need to recheck vert univ cuts.
    VertUniverseInfo vertical_universe_info;

    // Loop universes, make cuts, and fill
    for (auto error_band : error_bands) {
      std::vector<CVUniverse*> universes = error_band.second;
      for (auto universe : universes) {
        universe->SetEntry(i_event);

        // CCPiEvent keeps track of lots of event properties
        CCPiEvent event(is_mc, is_truth, signal_definition, universe);
        event.m_weight = universe->GetWeight();

        // Fill truth -- internally update the histograms owned by variables 
        if (is_truth) {
          ccpi_event::FillTruthEvent(event, variables);
          continue;
        }

        // Macro-level event reco (computationally intensive).
        // 
        // Construct michels/pions and check cuts.
        //
        // As we fail cuts (and we're not sideband either), don't waste time
        // continuing to process the event.
        LoopStatusCode status;
        std::tie(vertical_universe_info, status) = event.Process(vertical_universe_info);

        if (status == LoopStatusCode::SKIP)
          continue;

        // Re-call GetWeight because the node cut efficiency systematic
        // needs a pion candidate to calculate its weight.
        event.m_weight = universe->GetWeight();

        // Fill reco -- enforce cuts, internally update the histograms owned by variables
        ccpi_event::FillRecoEvent(event, variables);
      }  // universes
    }    // error bands
  } // entries 
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void makeCrossSectionMCInputs(int signal_definition_int = 0,
                              std::string plist = "ME1A",
                              bool do_systematics = false,
                              bool do_truth = false,
                              const bool do_test_playlist = false,
                              bool is_grid = false, std::string input_file = "",
                              int run = 0) {
  // 1. Input Data
  const bool is_mc = true;
  const bool use_xrootd = true;
  assert(!(is_grid && input_file.empty()) &&
         "On the grid, individual infile must be specified.");
  std::string mc_file_list =
      input_file.empty()
          ? CCPi::GetPlaylistFile(plist, is_mc, do_test_playlist, use_xrootd)
          : input_file;

  // 2. Macro Utility/Manager -- keep track of macro-level things, creat input
  // data chain
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, plist, do_truth,
                       is_grid, do_systematics);
  util.m_name = "MCXSecInputs";
  util.PrintMacroConfiguration();

  // 3. Prepare Output
  std::string outfile_name = GetOutFilename(util, run);
  std::cout << "Saving output to " << outfile_name << "\n\n";
  TFile fout(outfile_name.c_str(), "RECREATE");

  // 4. Initialize Variables (and the histograms that they own)
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  // 5. Loop MC Reco -- process events and fill histograms owned by variables
  bool is_truth = false;
  LoopAndFillMCXSecInputs(util.m_error_bands, util.GetMCEntries(), is_truth,
                          util.m_signal_definition, variables);

  // 6. Loop Truth
  if (util.m_do_truth) {
    is_truth = true;
    LoopAndFillMCXSecInputs(util.m_error_bands_truth, util.GetTruthEntries(),
                            is_truth, util.m_signal_definition, variables);
  }

  // 7. Write to file
  std::cout << "Synching and Writing\n\n";
  WritePOT(fout, is_mc, util.m_mc_pot);
  fout.cd();
  for (auto v : variables) {
    SyncAllHists(*v);
    v->WriteMCHists(fout);
  }
}

#endif  // makeXsecMCInputs_C

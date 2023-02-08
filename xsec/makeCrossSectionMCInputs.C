#ifndef makeXsecMCInputs_C
#define makeXsecMCInputs_C

#include <cassert>
#include <ctime>

#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/Constants.h"
#include "includes/Cuts.h"
#include "includes/MacroUtil.h"
#include "includes/SignalDefinition.h"
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/common_functions.h"           // GetVar, WritePOT

#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/HadronVariable.h"
#include "includes/HadronVariable2D.h"
#include "includes/Variable2D.h"
#include "includes/Variable.h"

class Variable;
class Variable2D;
class HadronVariable;
class HadronVariable2D;
//==============================================================================
// Helper Functions
//==============================================================================
namespace make_xsec_mc_inputs {
typedef Variable Var;
typedef Variable2D Var2D;
typedef HadronVariable HVar;
typedef HadronVariable2D HVar2D;

std::vector<Variable*> GetOnePiVariables(bool include_truth_vars = true) {
  const int nadphibins = 16;
  const double adphimin = -CCNuPionIncConsts::PI;
  const double adphimax = CCNuPionIncConsts::PI;

  HVar* tpi = new HVar("tpi", "T_{#pi}", "MeV", CCPi::GetBinning("tpi"),
                       &CVUniverse::GetTpi);

  HVar* tpi_mbr = new HVar("tpi_mbr", "T_{#pi} (MBR)", tpi->m_units,
                           CCPi::GetBinning("tpi"), &CVUniverse::GetTpiMBR);

  HVar* thetapi_deg =
      new HVar("thetapi_deg", "#theta_{#pi}", "deg",
               CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);

  HVar* pimuAngle =
      new HVar("pimuAngle", "p_{#pi}p_{#mu} angle", "deg",
               CCPi::GetBinning("pimuAngle"), &CVUniverse::GetpimuAngle);

  HVar* PT =
      new HVar("PT", "P^{T}", "MeV",
               CCPi::GetBinning("PT"), &CVUniverse::GetPT);
  HVar* ALR = new HVar("ALR", "ALR", "0cp,1L,2R", CCPi::GetBinning("ALR"),
                       &CVUniverse::GetALR);

  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);

  Var* thetamu_deg =
      new Var("thetamu_deg", "#theta_{#mu}", "deg",
              CCPi::GetBinning("thetamu_deg"), &CVUniverse::GetThetamuDeg);

  Var* enu = new Var("enu", "E_{#nu}", "MeV", CCPi::GetBinning("enu"),
                     &CVUniverse::GetEnu);

  Var* q2 = new Var("q2", "Q^{2}", "MeV^{2}", CCPi::GetBinning("q2"),
                    &CVUniverse::GetQ2);

  Var* wexp = new Var("wexp", "W_{exp}", "MeV", CCPi::GetBinning("wexp"),
                      &CVUniverse::GetWexp);

  Var* wexp_fit = new Var(sidebands::kFitVarString, wexp->m_hists.m_xlabel,
                          wexp->m_units, 32, 0.e3, 3.2e3, &CVUniverse::GetWexp);

  Var* ptmu = new Var("ptmu", "p^{T}_{#mu}", "MeV", CCPi::GetBinning("ptmu"),
                      &CVUniverse::GetPTmu);

  Var* pzmu = new Var("pzmu", "p^{z}_{#mu}", "MeV", CCPi::GetBinning("pzmu"),
                      &CVUniverse::GetPZmu);

  HVar* cosadtheta = new HVar("cosadtheta", "cos(#theta_{Adler})", "", CCPi::GetBinning("cosadtheta"),
                      &CVUniverse::GetAdlerCosTheta);

  HVar* adphi = new HVar("adphi", "#phi_{Adler}", "rad", nadphibins, adphimin, adphimax,
                      &CVUniverse::GetAdlerPhi);

  // True Variables
  bool is_true = true;
  HVar* tpi_true =
      new HVar("tpi_true", "T_{#pi} True", tpi->m_units,
               tpi->m_hists.m_bins_array, &CVUniverse::GetTpiTrue, is_true);

  HVar* thetapi_deg_true =
      new HVar("thetapi_deg_true", "#theta_{#pi} True", thetapi_deg->m_units,
               thetapi_deg->m_hists.m_bins_array,
               &CVUniverse::GetThetapiTrueDeg, is_true);

  HVar* ALR_true = new HVar("ALR", "ALR_True", ALR->m_units, ALR->m_hists.m_bins_array,
                       &CVUniverse::GetALRTrue, is_true);

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
      new Var("ptmu_true", "p^{t}_{#mu} True", "MeV", ptmu->m_hists.m_bins_array,
              &CVUniverse::GetPTmuTrue, is_true);

  Var* pzmu_true =
      new Var("pzmu_true", "p^{z}_{#mu} True", "MeV", pzmu->m_hists.m_bins_array,
              &CVUniverse::GetPZmuTrue, is_true);

  HVar* cosadtheta_true = new HVar("cosadtheta_true", "cos(#theta_{Adler}) True", "", cosadtheta->m_hists.m_bins_array,
                      &CVUniverse::GetAdlerCosThetaTrue, is_true);

  HVar* adphi_true = new HVar("adphi_true", "#phi_{Adler} True", "rad",  nadphibins, adphimin, adphimax,
                      &CVUniverse::GetAdlerPhiTrue, is_true);

  HVar* pimuAngle_true =
      new HVar("pimuAngle_true", "p_{#pi}p_{#mu} angle True", "deg",
               pimuAngle->m_hists.m_bins_array, &CVUniverse::GetpimuAngleTrue, is_true);

  HVar* PT_true =
      new HVar("PT_true", "P^{T} True", "MeV",
               PT->m_hists.m_bins_array, &CVUniverse::GetPTTrue, is_true);

  // Ehad variables
  Var* ehad = new Var("ehad", "ehad", "MeV", CCPi::GetBinning("ehad"),
                      &CVUniverse::GetEhad);
  Var* ehad_true = new Var("ehad_true", "ehad True", "MeV", ehad->m_hists.m_bins_array, 
                           &CVUniverse::GetEhadTrue);
  ehad_true->m_is_true = true;

  std::vector<Var*> variables = {tpi,      /* tpi_mbr, thetapi_deg,*/ pmu,
                                 thetamu_deg, /*enu,     q2,          wexp,*/
                                 wexp_fit,    ptmu,    pzmu/*       ehad,*/
				/* cosadtheta,  adphi,   pimuAngle,   PT, ALR*/};

  if (include_truth_vars) {
    variables.push_back(tpi_true);
//  variables.push_back(thetapi_deg_true);
    variables.push_back(pmu_true);
    variables.push_back(thetamu_deg_true);
//  variables.push_back(enu_true);
//  variables.push_back(q2_true);
    variables.push_back(wexp_true);
    variables.push_back(ptmu_true);
    variables.push_back(pzmu_true);
//  variables.push_back(ehad_true);
//  variables.push_back(cosadtheta_true);
//  variables.push_back(adphi_true);
//  variables.push_back(pimuAngle_true);
//  variables.push_back(PT_true);
  //variables.push_back(ALR_true);
  }

  return variables;
}

std::vector<HadronVariable*> GetOnePiHadronVariables(bool include_truth_vars = true){
  // Reco Variables
  HVar* tpi = new HVar("tpi", "T_{#pi}", "MeV", CCPi::GetBinning("tpi"),
                       &CVUniverse::GetTpi);
  // True Variables
  bool is_true = true;
  HVar* tpi_true =
      new HVar("tpi_true", "T_{#pi} True", tpi->m_units,
               tpi->m_hists.m_bins_array, &CVUniverse::GetTpiTrue, is_true);
  std::vector<HVar*> variables = {tpi};
  
  if (include_truth_vars) {
   variables.push_back(tpi_true);
  }
  return variables;
}

std::vector<Variable2D*> GetOnePiVariables2D(bool include_truth_vars = true){

  std::vector<Variable*> var1D = GetOnePiVariables(true);
  std::vector<HadronVariable*> Hvar1D = GetOnePiHadronVariables(true);
  bool is_true = true;
  // Reco 2D Variables 
  Var2D* pzmu_pTmu = new Var2D(var1D[5], var1D[4]);
  Var2D* pmu_thetamu = new Var2D(var1D[1], var1D[2]);

  HVar2D* tpi_thetapi_deg = new HVar2D("tpi", "thetapi_deg", "T_{#pi}",
			       "#theta_{#pi}", "MeV", "deg",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("thetapi_deg"), 
                               &CVUniverse::GetTpi, &CVUniverse::GetThetapiDeg);

  HVar2D* tpi_pmu = new HVar2D("tpi", "pmu", "T_{#pi}", "p_{#mu}", "MeV", "MeV",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("pmu"),
                               &CVUniverse::GetTpi, &CVUniverse::GetPmu);

  HVar2D* pTmu_tpi = new HVar2D("ptmu", "tpi", "p^{t}_{#mu}", "T_{#pi}", "MeV", "MeV",
                               CCPi::GetBinning("ptmu"), CCPi::GetBinning("tpi"),
                               &CVUniverse::GetPTmu, &CVUniverse::GetTpi);

  // True 2d Variables
  Var2D* pzmu_pTmu_true = new Var2D(var1D[11], var1D[10]);
  Var2D* pmu_thetamu_true = new Var2D(var1D[7], var1D[8]);

  HVar2D* tpi_thetapi_deg_true = new HVar2D("tpi_true", "thetapi_deg_true",
                               "T_{#pi} true", "#theta_{#pi} true", "MeV", "deg",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("thetapi_deg"),
                               &CVUniverse::GetTpiTrue, &CVUniverse::GetThetapiTrueDeg,
                               is_true);

  HVar2D* tpi_pmu_true = new HVar2D("tpi_true", "pmu_true", "T_{#pi} True",
                               "p_{#mu} True", "MeV", "MeV",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("pmu"),
                               &CVUniverse::GetTpiTrue, &CVUniverse::GetPmuTrue, is_true, 4);

  HVar2D* pTmu_tpi_true = new HVar2D("ptmu_true", "tpi_true", "p^{t}_{#mu} True",
			       "T_{#pi} True", "MeV", "MeV",
                               CCPi::GetBinning("ptmu"), CCPi::GetBinning("tpi"),
                               &CVUniverse::GetPTmuTrue, &CVUniverse::GetTpiTrue, is_true);

  std::vector<Var2D*> variables2D = {pzmu_pTmu, tpi_thetapi_deg, tpi_pmu, pmu_thetamu, pTmu_tpi};
  if (include_truth_vars){
    variables2D.push_back(pzmu_pTmu_true);
    variables2D.push_back(tpi_thetapi_deg_true);
    variables2D.push_back(tpi_pmu_true);
    variables2D.push_back(pmu_thetamu_true);
    variables2D.push_back(pTmu_tpi_true);
  }
  return variables2D;
}

std::map<std::string, Variable*> GetOnePiVariables_Map(
    bool include_truth_vars = true) {
  std::map<std::string, Var*> var_map;
  std::vector<Var*> var_vec = GetOnePiVariables(include_truth_vars);
  for (auto v : var_vec) var_map[v->Name()] = v;
  return var_map;
}

std::map<std::string, Variable2D*> GetOnePiVariables2D_Map(
    bool include_truth_vars = true) {
  std::map<std::string, Var2D*> var2D_map;
  std::vector<Var2D*> var_vec = GetOnePiVariables2D(include_truth_vars);
  for (auto v : var_vec) var2D_map[v->NameX() + "_vs_" + v->NameY()] = v;
  return var2D_map;
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

std::vector<Variable2D*> GetAnalysisVariables2D(SignalDefinition signal_definition,
                                            bool include_truth_vars = false) {
  std::vector<Variable2D*> variables2D;
  switch (signal_definition) {
    case kOnePi:
      variables2D = make_xsec_mc_inputs::GetOnePiVariables2D(include_truth_vars);
      break;
    default:
      std::cerr << "Variables for other SDs not yet implemented.\n";
      std::exit(1);
  }
  return variables2D;
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

void SyncAllHists2D(Variable2D& v2D){
  v2D.m_hists2D.m_selection_mc.SyncCVHistos();
  v2D.m_hists2D.m_bg.SyncCVHistos();
  v2D.m_hists2D.m_bg_loW.SyncCVHistos();
  v2D.m_hists2D.m_bg_midW.SyncCVHistos();
  v2D.m_hists2D.m_bg_hiW.SyncCVHistos();
//v2D.m_hists2D.m_wsidebandfit_sig.SyncCVHistos();
//v2D.m_hists2D.m_wsidebandfit_loW.SyncCVHistos();
//v2D.m_hists2D.m_wsidebandfit_midW.SyncCVHistos();
//v2D.m_hists2D.m_wsidebandfit_hiW.SyncCVHistos();
  v2D.m_hists2D.m_effnum.SyncCVHistos();
  v2D.m_hists2D.m_effden.SyncCVHistos();
}
//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillMCXSecInputs(const CCPi::MacroUtil& util,
                             const EDataMCTruth& type,
                             std::vector<Variable*>& variables,
			     std::vector<Variable2D*>& variables2D) {
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  const UniverseMap error_bands =
      is_truth ? util.m_error_bands_truth : util.m_error_bands;
  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
  //for (Long64_t i_event = 0; i_event < 5000; ++i_event) {
    if (i_event % (n_entries / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;

//    if (i_event == 5000) break; 
    // Variables that hold info about whether the CVU passes cuts
    PassesCutsInfo cv_cuts_info;
    bool checked_cv = false;

    // Loop universes, make cuts, and fill
    for (auto error_band : error_bands) {
      std::vector<CVUniverse*> universes = error_band.second;
      for (auto universe : universes) {
        universe->SetEntry(i_event);
        //std::cout << universe->ShortName() << "\n";
        //if (universe->GetDouble("mc_incoming") == 12 &&
        //    universe->ShortName() == "cv")
        //  universe->PrintArachneLink();
        CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe); // call GetWeight

        //===============
        // FILL TRUTH
        //===============
        if (type == kTruth) {
          ccpi_event::FillTruthEvent(event, variables);
          ccpi_event::FillTruthEvent2D(event, variables2D);
          continue;
        }

        //===============
        // CHECK CUTS
        //===============
        // Universe only affects weights
        if (universe->IsVerticalOnly()) {
          // Only check vertical-only universes once.
          if (!checked_cv) {
            // Check cuts
            cv_cuts_info = PassesCuts(event);
            checked_cv = true;
          }

          // Already checked a vertical-only universe
          if (checked_cv) {
            std::tie(event.m_passes_cuts, event.m_is_w_sideband, event.m_passes_all_cuts_except_w, event.m_reco_pion_candidate_idxs) = cv_cuts_info.GetAll();
            event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
          }
        // Universe shifts something laterally
        } else {
          PassesCutsInfo cuts_info = PassesCuts(event);
          std::tie(event.m_passes_cuts, event.m_is_w_sideband, event.m_passes_all_cuts_except_w, event.m_reco_pion_candidate_idxs) = cuts_info.GetAll();
          event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
        }

        // The universe needs to know its pion candidates in order to calculate
        // recoil/hadronic energy.
        universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);

        // Need to re-call this because the node cut efficiency systematic
        // needs a pion candidate to calculate its weight.
        event.m_weight = universe->GetWeight();

        //===============
        // FILL RECO
        //===============
        ccpi_event::FillRecoEvent(event, variables);
        ccpi_event::FillRecoEvent2D(event, variables2D);
      }  // universes
    }    // error bands
  }      // events
  std::cout << "*** Done ***\n\n";
}

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
  assert(!(is_grid && input_file.empty()) &&
         "On the grid, infile must be specified.");
  // const bool use_xrootd = false;
  mc_file_list = input_file.empty()
                     ? GetPlaylistFile(plist, is_mc /*, use_xrootd*/)
                     : input_file;

  // INIT MACRO UTILITY
  const std::string macro("MCXSecInputs");
  // std::string a_file =
  // "root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr/minerva/persistent/users/bmesserl/pions//20200713/merged/mc/ME1A/CCNuPionInc_mc_AnaTuple_run00110000_Playlist.root";
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, plist, do_truth,
                       is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // INIT OUTPUT
  auto time = std::time(nullptr);
  char tchar[100];
  std::strftime(tchar, sizeof(tchar), "%F", std::gmtime(&time));  // YYYY-MM-dd
  const std::string tag = tchar;
  std::string outfile_name(Form("%s_%d%d%d%d_%s_%d_%s.root", macro.c_str(),
                                signal_definition_int, int(do_systematics),
                                int(do_truth), int(is_grid), plist.c_str(), run,
                                tag.c_str()));
  std::cout << "Saving output to " << outfile_name << "\n\n";
  TFile fout(outfile_name.c_str(), "RECREATE");

  // INIT VARS, HISTOS, AND EVENT COUNTERS
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  std::vector<Variable2D*> variables2D =
      GetAnalysisVariables2D(util.m_signal_definition, do_truth_vars);

  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  for (auto v2D : variables2D){
    v2D->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);
  }
  // LOOP MC RECO
  for (auto band : util.m_error_bands) {
    std::vector<CVUniverse*> universes = band.second;
    for (auto universe : universes)
      universe->SetTruth(false);
  }
  LoopAndFillMCXSecInputs(util, kMC, variables, variables2D);

  // LOOP TRUTH
  if (util.m_do_truth) {
    // m_is_truth is static, so we turn it on now
    for (auto band : util.m_error_bands_truth) {
      std::vector<CVUniverse*> universes = band.second;
      for (auto universe : universes)
        universe->SetTruth(true);
    }
    LoopAndFillMCXSecInputs(util, kTruth, variables, variables2D);
  }

  // WRITE TO FILE
  std::cout << "Synching and Writing\n\n";
  WritePOT(fout, true, util.m_mc_pot);
  fout.cd();
  for (auto v : variables) {
    SyncAllHists(*v);
    v->WriteMCHists(fout);
  }
  for (auto v2D : variables2D) {
    v2D->m_response.GetMigrationObjects(v2D->m_hists2D.m_response, v2D->m_hists2D.m_responseReco, v2D->m_hists2D.m_responseTrue);
//    v2D->m_hists2D.m_response = v2D->m_response.GetMigrationMatrix();
    SyncAllHists2D(*v2D);
    v2D->WriteMCHists(fout);
  }

}

#endif  // makeXsecMCInputs_C

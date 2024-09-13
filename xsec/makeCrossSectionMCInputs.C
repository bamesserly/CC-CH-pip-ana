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
#include "includes/MacroUtil.h"
#include "includes/SignalDefinition.h"
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/common_functions.h"           // GetVar, WritePOT
#include "includes/utilities.h"

#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/HadronVariable.h"
#include "includes/HadronVariable2D.h"
#include "includes/Variable2D.h"
#include "includes/Variable.h"

//class Variable;
//class Variable2D;
//class HadronVariable;
//class HadronVariable2D;
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

  Var* wexp_fit =
      new Var(sidebands::kFitVarString, wexp->m_hists.m_xlabel, wexp->m_units,
              CCPi::GetBinning("wexp_fit"), &CVUniverse::GetWexp);

  Var* ptmu = new Var("ptmu", "p^{T}_{#mu}", "MeV", CCPi::GetBinning("ptmu"),
                      &CVUniverse::GetPTmu);

  Var* pzmu = new Var("pzmu", "p^{z}_{#mu}", "MeV", CCPi::GetBinning("pzmu"),
                      &CVUniverse::GetPZmu);

  HVar* cosadtheta = new HVar("cosadtheta", "cos(#theta_{Adler})", "", CCPi::GetBinning("cosadtheta"),
                      &CVUniverse::GetAdlerCosTheta);

  HVar* adphi = new HVar("adphi", "#phi_{Adler}", "rad", nadphibins, adphimin, adphimax,
                      &CVUniverse::GetAdlerPhi);

  Var* mehreen_tpi = new Var("mtpi", "Mehreen T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
                     &CVUniverse::GetTpiTrackless);

  Var* mehreen_thetapi_deg = new Var("mthetapi_deg", "Mehreen #theta_{#pi}", "deg",
                     CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapitracklessDeg);

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

  Var* mehreen_tpi_true =
      new Var("mtpi_true", "Mehreen T_{#pi} True", "MeV",
               mehreen_tpi->m_hists.m_bins_array, &CVUniverse::GetTrueTpi,
               is_true);

  Var* mehreen_thetapi_deg_true =
      new Var("mthetapi_deg_true", "Mehreen #theta_{#pi} True", "deg",
               mehreen_thetapi_deg->m_hists.m_bins_array, &CVUniverse::GetThetapitracklessTrueDeg,
               is_true);

  // Ehad variables
  Var* ehad = new Var("ehad", "ehad", "MeV", CCPi::GetBinning("ehad"),
                      &CVUniverse::GetEhad);
  Var* ehad_true = new Var("ehad_true", "ehad True", "MeV", ehad->m_hists.m_bins_array, 
                           &CVUniverse::GetEhadTrue);
  ehad_true->m_is_true = true;

  std::vector<Var*> variables = {tpi,      /* tpi_mbr, thetapi_deg,*/ pmu,
                                 thetamu_deg, enu,/*     q2,          wexp,*/
                                 wexp_fit,    ptmu,    pzmu/*       ehad,*/
				/* cosadtheta,  adphi,   pimuAngle,   PT, ALR*/};

  if (include_truth_vars) {
    variables.push_back(tpi_true);
//  variables.push_back(thetapi_deg_true);
    variables.push_back(pmu_true);
    variables.push_back(thetamu_deg_true);
    variables.push_back(enu_true);
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
  Var2D* pzmu_pTmu = new Var2D(var1D[6], var1D[5]);
  Var2D* pTmu_pzmu = new Var2D(var1D[5], var1D[6]);
//  Var2D* pmu_thetamu = new Var2D(var1D[1], var1D[2]);

  HVar2D* tpi_thetapi_deg = new HVar2D("tpi", "thetapi_deg", "T_{#pi}",
			       "#theta_{#pi}", "MeV", "deg",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("thetapi_deg"), 
                               &CVUniverse::GetTpi, &CVUniverse::GetThetapiDeg);
  HVar2D* thetapi_deg_tpi = new HVar2D("thetapi_deg", "tpi", "#theta_{#pi}","T_{#pi}",
			       "deg", "MeV", 
                               CCPi::GetBinning("thetapi_deg"), CCPi::GetBinning("tpi"), 
                               &CVUniverse::GetThetapiDeg, &CVUniverse::GetTpi);

  HVar2D* tpi_pmu = new HVar2D("tpi", "pmu", "T_{#pi}", "p_{#mu}", "MeV", "MeV",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("pmu_with_tpi"),
                               &CVUniverse::GetTpi, &CVUniverse::GetPmu);
  HVar2D* pmu_tpi = new HVar2D("pmu", "tpi", "p_{#mu}", "T_{#pi}", "MeV", "MeV",
                               CCPi::GetBinning("pmu_with_tpi"), CCPi::GetBinning("tpi"),
                               &CVUniverse::GetPmu, &CVUniverse::GetTpi);


  HVar2D* pTmu_tpi = new HVar2D("ptmu", "tpi", "p^{t}_{#mu}", "T_{#pi}", "MeV", "MeV",
                               CCPi::GetBinning("ptmu_with_tpi"),
                               CCPi::GetBinning("tpi_with_ptmu"),
                               &CVUniverse::GetPTmu, &CVUniverse::GetTpi);
  HVar2D* tpi_pTmu = new HVar2D("tpi", "ptmu", "T_{#pi}", "p^{t}_{#mu}", "MeV", "MeV",
                               CCPi::GetBinning("tpi_with_ptmu"),
                               CCPi::GetBinning("ptmu_with_tpi"),
                               &CVUniverse::GetTpi, &CVUniverse::GetPTmu);

  HVar2D* tpi_enu = new HVar2D("tpi", "enu", "T_{#pi}", "E_{#nu}", "MeV", "MeV",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("enu_with_tpi"),
                               &CVUniverse::GetTpi, &CVUniverse::GetEnu);
  HVar2D* enu_tpi = new HVar2D("enu", "tpi", "E_{#nu}", "T_{#pi}", "MeV", "MeV",
                               CCPi::GetBinning("enu_with_tpi"), CCPi::GetBinning("tpi"),
                               &CVUniverse::GetEnu, &CVUniverse::GetTpi);
//  Var2D* pzmu_thetamu = new Var2D(var1D[5], var1D[2]);
//  Var2D* ptmu_thetamu = new Var2D(var1D[4], var1D[2]);


  // True 2d Variables
  Var2D* pzmu_pTmu_true = new Var2D(var1D[13], var1D[12]);
  Var2D* pTmu_pzmu_true = new Var2D(var1D[12], var1D[13]);
//  Var2D* pmu_thetamu_true = new Var2D(var1D[7], var1D[8]);

  HVar2D* tpi_thetapi_deg_true = new HVar2D("tpi_true", "thetapi_deg_true",
                               "T_{#pi} true", "#theta_{#pi} true", "MeV", "deg",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("thetapi_deg"),
                               &CVUniverse::GetTpiTrue, &CVUniverse::GetThetapiTrueDeg,
                               is_true);
  HVar2D* thetapi_deg_tpi_true = new HVar2D("thetapi_deg_true", "tpi_true", 
                               "#theta_{#pi} true", "T_{#pi} true", "deg",  "MeV",
                               CCPi::GetBinning("thetapi_deg"), CCPi::GetBinning("tpi"),
                               &CVUniverse::GetThetapiTrueDeg, &CVUniverse::GetTpiTrue, 
                               is_true);

  HVar2D* tpi_pmu_true = new HVar2D("tpi_true", "pmu_true", "T_{#pi} True",
                               "p_{#mu} True", "MeV", "MeV",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("pmu_with_tpi"),
                               &CVUniverse::GetTpiTrue, &CVUniverse::GetPmuTrue, is_true);
  HVar2D* pmu_tpi_true = new HVar2D("pmu_true", "tpi_true", "p_{#mu} True", "T_{#pi} True",
                               "MeV", "MeV",
                               CCPi::GetBinning("pmu_with_tpi"), CCPi::GetBinning("tpi"),
                               &CVUniverse::GetPmuTrue, &CVUniverse::GetTpiTrue, is_true);

  HVar2D* pTmu_tpi_true = new HVar2D("ptmu_true", "tpi_true", "p^{t}_{#mu} True",
			       "T_{#pi} True", "MeV", "MeV",
                               CCPi::GetBinning("ptmu_with_tpi"),
                               CCPi::GetBinning("tpi_with_ptmu"),
                               &CVUniverse::GetPTmuTrue, &CVUniverse::GetTpiTrue, is_true);
  HVar2D* tpi_pTmu_true = new HVar2D("tpi_true", "ptmu_true", "T_{#pi} True",
			       "p^{t}_{#mu} True", "MeV", "MeV",
                               CCPi::GetBinning("tpi_with_ptmu"),
                               CCPi::GetBinning("ptmu_with_tpi"),
                               &CVUniverse::GetTpiTrue, &CVUniverse::GetPTmuTrue, is_true);

  HVar2D* tpi_enu_true = new HVar2D("tpi_true", "enu_true", "T_{#pi} True", "E_{#nu} True",
			       "MeV", "MeV",
                               CCPi::GetBinning("tpi"), CCPi::GetBinning("enu_with_tpi"),
                               &CVUniverse::GetTpiTrue, &CVUniverse::GetEnuTrue, is_true);
  HVar2D* enu_tpi_true = new HVar2D("enu_true", "tpi_true", "E_{#nu} True", "T_{#pi} True",
			       "MeV", "MeV",
                               CCPi::GetBinning("enu_with_tpi"), CCPi::GetBinning("tpi"),
                               &CVUniverse::GetEnuTrue, &CVUniverse::GetTpiTrue, is_true);

  pzmu_pTmu_true->m_is_true = true;
//  Var2D* ptmu_thetamu_true = new Var2D(var1D[10], var1D[8]);
  std::vector<Var2D*> variables2D = {pzmu_pTmu, pTmu_pzmu, tpi_thetapi_deg, thetapi_deg_tpi, pmu_tpi, tpi_pmu, tpi_pTmu, pTmu_tpi, tpi_enu, enu_tpi};
  if (include_truth_vars){
    variables2D.push_back(pzmu_pTmu_true);
    variables2D.push_back(pTmu_pzmu_true);
    variables2D.push_back(tpi_thetapi_deg_true);
    variables2D.push_back(thetapi_deg_tpi_true);
    variables2D.push_back(tpi_pmu_true);
    variables2D.push_back(pmu_tpi_true);
    variables2D.push_back(pTmu_tpi_true);
    variables2D.push_back(tpi_pTmu_true);
    variables2D.push_back(tpi_enu_true);
    variables2D.push_back(enu_tpi_true);
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

std::vector<Variable2D*> GetAnalysisVariables2D(const SignalDefinition signal_definition,
                                            bool include_truth_vars = false) {
  std::vector<Variable2D*> variables2D;
  variables2D = make_xsec_mc_inputs::GetOnePiVariables2D(include_truth_vars);
/*  switch (signal_definition) {
    case kOnePi:
      variables2D = make_xsec_mc_inputs::GetOnePiVariables2D(include_truth_vars);
      break;
    default:
      std::cerr << "Variables for other SDs not yet implemented.\n";
      std::exit(1);
  }*/
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
  v2D.m_hists2D.m_migration.SyncCVHistos();
  v2D.m_hists2D.m_migration_reco.SyncCVHistos();
  v2D.m_hists2D.m_migration_true.SyncCVHistos();
  v2D.m_hists2D.m_effden.SyncCVHistos();
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
                             std::vector<Variable*>& variables,
			     std::vector<Variable2D*>& variables2D) {
  const bool is_mc = true;

  for (auto band : error_bands) {
    std::vector<CVUniverse*> universes = band.second;
    for (auto universe : universes) universe->SetTruth(is_truth);
  }
  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
  //for (Long64_t i_event = 0; i_event < 5000; ++i_event) {
    if (i_event % (n_entries / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;

    //if (i_event == 50000) break; 
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

        // CCPiEvent keeps track of lots of event properties
        CCPiEvent event(is_mc, is_truth, signal_definition, universe);

        //===============
        // FILL TRUTH
        //===============
        if (type == kTruth) {
          universe->SetPassesTrakedTracklessCuts(true,true);
          ccpi_event::FillTruthEvent(event, variables);
          ccpi_event::FillTruthEvent2D(event, variables2D);
          continue;
        }

        // Check Cuts -- computationally expensive
        //
        // This looks complicated for optimization reasons.
        // Namely, for all vertical-only universes (meaning only the event
        // weight differs from CV) no need to recheck cuts.
        PassesCutsInfo cuts_info;
        LowRecoilPion::Cluster d;
        LowRecoilPion::Cluster c(*universe,0);
        LowRecoilPion::Michel<CVUniverse> m(*universe,0);
        LowRecoilPion::MichelEvent<CVUniverse> trackless_michels;

        //===============
        // CHECK CUTS
        //===============
        // Universe only affects weights

        bool good_trackless_michels = LowRecoilPion::hasMichel<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::hasMichelCut(*universe, trackless_michels);
        good_trackless_michels = good_trackless_michels && LowRecoilPion::BestMichelDistance2D<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::BestMichelDistance2DCut(*universe, trackless_michels);
        good_trackless_michels = good_trackless_michels && LowRecoilPion::GetClosestMichel<CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::GetClosestMichelCut(*universe, trackless_michels);

//        if (good_trackless_michels) {
//	  trackless_michels.m_allmichels[trackless_michels.m_idx].GetPionAngle(*universe);
//          std::cout << "Angle = " << trackless_michels.m_bestthetaangle << "\n";
          
//        }
        universe->SetVtxMichels(trackless_michels);

        bool pass = true;
        pass = pass && universe->GetNMichels() == 1;
        pass = pass && universe->GetTpiTrackless() < 350.;
        pass = pass && universe->GetPmu() > 1500.;
        pass = pass && universe->GetPmu() < 20000.;
        pass = pass && universe->GetNIsoProngs() < 2;
        pass = pass && universe->IsInHexagon(universe->GetVecElem("vtx", 0), universe->GetVecElem("vtx", 1), 850.);
        pass = pass && universe->GetVecElem("vtx", 2) > 5990.;
        pass = pass && universe->GetVecElem("vtx", 2) < 8340.;
        pass = pass && universe->GetBool("isMinosMatchTrack");
        pass = pass && universe->GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
        pass = pass && universe->GetThetamuDeg() < 20;

        if (universe->IsVerticalOnly()) {
          if (!checked_cv) {
            cv_cuts_info = PassesCuts(event);
            checked_cv = true;
          }
          assert(checked_cv);
          cuts_info = cv_cuts_info;
        } else {
          cuts_info = PassesCuts(event);
        }

        // Save results of cuts to Event and universe
        std::tie(event.m_passes_cuts, event.m_is_w_sideband,
                 event.m_passes_all_cuts_except_w,
                 event.m_reco_pion_candidate_idxs) = cuts_info.GetAll();

        event.m_highest_energy_pion_idx =
            GetHighestEnergyPionCandidateIndex(event);

        universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);

        // Re-call GetWeight because the node cut efficiency systematic
        // Need to re-call this because the node cut efficiency systematic
        // needs a pion candidate to calculate its weight.
        event.m_weight = universe->GetWeight();
        universe->SetVtxMichels(trackless_michels);
        event.m_passes_trackless_cuts_except_w = false;
        event.m_passes_trackless_sideband = false;
        if (pass && universe->GetWexp() > 1400.){
          event.m_passes_trackless_cuts_except_w = true;
          if (universe->GetWexp() > 1500.) event.m_passes_trackless_sideband = true;
          pass = false;
        }
        event.m_passes_trackless_cuts = good_trackless_michels && pass;
        event.m_passes_trackless_sideband = event.m_passes_trackless_sideband && good_trackless_michels;
        event.m_passes_trackless_cuts_except_w = event.m_passes_trackless_cuts_except_w && good_trackless_michels;
        universe->SetPassesTrakedTracklessCuts(event.m_passes_cuts || event.m_passes_all_cuts_except_w, event.m_passes_trackless_cuts || event.m_passes_trackless_cuts_except_w);
        

        //===============
        // FILL RECO
        //===============
/*        if (i_event == 4216){
          std::cout << "Event = " << i_event << "\n";
          std::cout << "Pass the cuts? " << event.m_passes_trackless_cuts << "\n"; 
          std::cout << "Pass the sidebands? " << event.m_passes_trackless_cuts << "\n"; 
          std::cout << "Pass cuts exept W " << event.m_passes_trackless_cuts << "\n"; 
          std::cout << "Pass tracked cuts? " << event.m_passes_cuts << "\n"; 
          std::cout << "Pass tracked cuts? " << event.m_is_w_sideband << "\n"; 
        }*/
/*        if (event.m_is_signal && event.m_passes_trackless_cuts){
          std::cout << "Is Signal and pass cuts, Event " << i_event << "\n";
	  std::cout << "Thetapi reco = " << universe->GetThetapitracklessDeg() << "\n";
  	  std::cout << "Thetapi true = " << universe->GetThetapitracklessTrueDeg() << "\n";
//          std::cout << "thetapi from trackless_michels = " << trackless_michels.best_angle << "\n";
        }*/
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

  std::vector<Variable2D*> variables2D =
      GetAnalysisVariables2D(util.m_signal_definition, do_truth_vars);

  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  for (auto v2D : variables2D){
    v2D->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);
  }
  // 5. Loop MC Reco -- process events and fill histograms owned by variables
  bool is_truth = false;
  LoopAndFillMCXSecInputs(util.m_error_bands, util.GetMCEntries(), is_truth,
                          util.m_signal_definition, variables, variables2D);

  // 6. Loop Truth
  if (util.m_do_truth) {
    is_truth = true;
    LoopAndFillMCXSecInputs(util.m_error_bands_truth, util.GetTruthEntries(),
                            is_truth, util.m_signal_definition, variables, variables2D);
  }

  // 7. Write to file
  std::cout << "Synching and Writing\n\n";
  WritePOT(fout, is_mc, util.m_mc_pot);
  fout.cd();
  for (auto v : variables) {
    SyncAllHists(*v);
    v->WriteMCHists(fout);
  }
  for (auto v2D : variables2D) {
//    v2D->m_response.GetMigrationObjects(v2D->m_hists2D.m_response, v2D->m_hists2D.m_responseReco, v2D->m_hists2D.m_responseTrue);
//    v2D->m_hists2D.m_response = v2D->m_response.GetMigrationMatrix();
    SyncAllHists2D(*v2D);
    v2D->WriteMCHists(fout);
  }

}

#endif  // makeXsecMCInputs_C

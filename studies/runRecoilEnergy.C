//==============================================================================
// Compare new and old methods for calculating Ehad.
// old: cal-corrected sum over non-muon clusters
// new: differently-cal-corrected sum over non-muon and non-pion clusters +
// pion KE.
//==============================================================================
#include <iostream>
#include <vector>

#include "ccpion_common.h" // GetPlaylistFile
#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/Constants.h" // typedef RecoPionIdx 
#include "includes/CVUniverse.h"
#include "includes/MacroUtil.h"
#include "includes/HadronVariable.h"
#include "includes/Variable.h"
#include "plotting_functions.h"

#include "includes/common_functions.h"

// Forward declare my variables because we're hiding the header.
class Variable;
//class HadronVariable;

double mean = 0.;
int counts = 0;

namespace run_recoil_energy {
//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  typedef Variable Var;
  typedef HadronVariable HVar;

  // Ecal
    const int necalbins = 25;
    const double ecalmin = 0.;
    const double ecalmax = 1800.;
    Var* ecalrecoil         = new Var("ecalrecoil",         "ecalrecoil",         "mev", necalbins, ecalmin, ecalmax, &CVUniverse::GetCalRecoilEnergy);
    Var* ecalrecoil_ccinc   = new Var("ecalrecoil_ccinc",   "ecalrecoil_ccinc",   "mev", necalbins, ecalmin, ecalmax, &CVUniverse::GetCalRecoilEnergy_CCIncSpline);
    Var* ecalrecoil_ccpi    = new Var("ecalrecoil_ccpi",    "ecalrecoil_ccpi",    "mev", necalbins, ecalmin, ecalmax, &CVUniverse::GetCalRecoilEnergy_CCPiSpline);
    Var* ecalrecoil_default = new Var("ecalrecoil_default", "ecalrecoil_default", "mev", necalbins, ecalmin, ecalmax, &CVUniverse::GetCalRecoilEnergy_DefaultSpline);

  // Ecal - nopi
    Var* ecalrecoilnopi       = new Var("ecalrecoilnopi",       "ecalrecoilnopi",       "mev", CCPi::GetBinning("ecal_nopi"), &CVUniverse::GetCalRecoilEnergyNoPi_DefaultSpline);
    Var* ecalrecoilnopi_ccinc = new Var("ecalrecoilnopi_ccinc", "ecalrecoilnopi_ccinc", "mev", CCPi::GetBinning("ecal_nopi"), &CVUniverse::GetCalRecoilEnergyNoPi_CCIncSpline);
    Var* ecalrecoilnopi_corr  = new Var("ecalrecoilnopi_corr",  "ecalrecoilnopi_corr",  "mev", CCPi::GetBinning("ecal_nopi"));
    Var* ecalrecoilnopi_true  = new Var("ecalrecoilnopi_true",  "ecalrecoilnopi_true",  "mev", CCPi::GetBinning("ecal_nopi"), &CVUniverse::GetCalRecoilEnergyNoPiTrue);
    ecalrecoilnopi_true->m_is_true = true;

  // Etracks
    const int netracksbins = 25;
    const double etracksmin = 100.;
    const double etracksmax = 900.;
    Var* etrackrecoil = new Var("etrackrecoil", "etrackrecoil", "mev", netracksbins, etracksmin, etracksmax, &CVUniverse::GetTrackRecoilEnergy);
    Var* etracks_true = new Var("etracks_true", "etracks_true", "mev", netracksbins, etracksmin, etracksmax, &CVUniverse::GetAllTrackEnergyTrue);
    etracks_true->m_is_true = true;

  // Ehad
    const int nehadbins = 25;
    const double ehadmin = 200.;
    const double ehadmax = 1800.;
    Var* ehad      = new Var("ehad",      "ehad",      "mev", nehadbins, ehadmin, ehadmax, &CVUniverse::GetEhad);
    Var* ehad_true = new Var("ehad_true", "ehad_true", "mev", nehadbins, ehadmin, ehadmax, &CVUniverse::GetEhadTrue);
    ehad_true->m_is_true = true;

  // pi vars
    const int nepibins = 30;
    const double epimin = 0.;
    const double epimax = 600.;
    HVar* epi_cal  = new HVar("epi_cal",  "epi_cal",  "mev", nepibins, epimin, epimax, &CVUniverse::GetCalEpi);
    HVar* epi_true = new HVar("epi_true", "epi_true", "mev", nepibins, epimin, epimax, &CVUniverse::GetEpiTrueMatched);
    epi_true->m_is_true = true;

  // W vars
    const int nwbins = 26;
    const double wmin = 600;
    const double wmax = 1900;
    Var* wexp      = new Var("wexp",      "wexp",      "mev", nwbins, wmin, wmax, &CVUniverse::GetWexp);
    Var* wexp_true = new Var("wexp_true", "wexp_true", "mev", nwbins, wmin, wmax, &CVUniverse::GetWexpTrue);
    Var* wgenie    = new Var("wgenie",    "wgenie",    "mev", nwbins, wmin, wmax, &CVUniverse::GetWgenie);
    Var* wresid    = new Var("wresid",    "wresid",    "mev", 30,     -0.5,   0.5,    &CVUniverse::GetWexpFResidual);
    wexp_true->m_is_true = true;
    wgenie->m_is_true    = true;
    wresid->m_is_true    = true;

  std::vector<Var*> variables = {ehad, ecalrecoil, etrackrecoil,
      ecalrecoilnopi, ecalrecoil_ccinc,
      ecalrecoil_default, ecalrecoil_ccpi, ecalrecoilnopi_ccinc, ehad_true,
      ecalrecoilnopi_true, epi_cal, epi_true, wexp, wexp_true, wgenie,
      ecalrecoilnopi_corr,etracks_true, 
      wresid};

  return variables;
}
//==============================================================================
// Do some event processing (e.g. make cuts, get best pion) and fill hists.
// What we're looking at in this study:
//==============================================================================
void FillVars(CCPiEvent& event, const std::vector<Variable*>& variables) {
  // Shortening event variables for convenience
  const CVUniverse* universe = event.m_universe;
  const double wgt = event.m_weight;
  const bool is_mc = event.m_is_mc;
  const SignalDefinition sd = event.m_signal_definition;
  RecoPionIdx best_pion = event.m_highest_energy_pion_idx;

  // No systematics considered here
  if (universe->ShortName() != "cv") return;

  /*
  if(is_mc) {
    TruePionIdx true_index = universe->GetVecElem("MasterAnaDev_hadron_tm_trackID", best_pion);
    double tt = true_index >= 0 ? universe->GetTpiTrue(true_index) : -9999.;
    std::cout << universe->GetTpiTrueMatched(best_pion) - tt << "\n";
  }
  */

  double ecalrecoil_nopi = universe->GetCalRecoilEnergyNoPi_DefaultSpline();
  double ecalrecoil_nopi_corr = universe->GetCalRecoilEnergyNoPi_Corrected(ecalrecoil_nopi);

  // Fill stacked histograms
  for (auto v : variables) {
    if(v->Name()=="ecalrecoilnopi_corr") 
      ccpi_event::FillStackedHists(event, v, ecalrecoil_nopi_corr);
    else
      ccpi_event::FillStackedHists(event, v);
  }

  // RETURN -- ONLY MC BEYOND THIS POINT
  if (!is_mc) return;
  

  // Fill migration histograms
    if(HasVar(variables,"ehad"))
      GetVar(variables, "ehad") -> m_hists.m_migration.FillUniverse(
          *universe, GetVar(variables, "ehad")->GetValue(*universe), GetVar(variables, "ehad_true")->GetValue(*universe), event.m_weight);

    if(HasVar(variables,"epi_cal")){
      double r = GetVar(variables, "epi_cal")->GetValue(*universe, best_pion);
      double t = GetVar(variables, "epi_true")->GetValue(*universe, best_pion);
      counts += 1;
      mean += t-r;   
      GetVar(variables, "epi_cal") -> m_hists.m_migration.FillUniverse(*universe, r, t, event.m_weight);
    }

    if(HasVar(variables,"ecalrecoilnopi"))
      GetVar(variables, "ecalrecoilnopi") -> m_hists.m_migration.FillUniverse(
          *universe, GetVar(variables, "ecalrecoilnopi")->GetValue(*universe), GetVar(variables, "ecalrecoilnopi_true")->GetValue(*universe), event.m_weight);

    if(HasVar(variables,"ecalrecoilnopi_ccinc"))
      GetVar(variables, "ecalrecoilnopi_ccinc") -> m_hists.m_migration.FillUniverse(
          *universe, GetVar(variables, "ecalrecoilnopi_ccinc")->GetValue(*universe), GetVar(variables, "ecalrecoilnopi_true")->GetValue(*universe), event.m_weight);

    if(HasVar(variables,"ecalrecoilnopi_corr"))
      GetVar(variables, "ecalrecoilnopi_corr") -> m_hists.m_migration.FillUniverse(
          *universe, ecalrecoil_nopi_corr, GetVar(variables, "ecalrecoilnopi_true")->GetValue(*universe), event.m_weight);

    if(HasVar(variables,"ecalrecoil"))
      GetVar(variables, "ecalrecoil") -> m_hists.m_migration.FillUniverse(
          *universe, GetVar(variables, "ecalrecoil")->GetValue(*universe), GetVar(variables, "ehad_true")->GetValue(*universe), event.m_weight);

    if(HasVar(variables,"ecalrecoil_default"))
      GetVar(variables, "ecalrecoil_default") -> m_hists.m_migration.FillUniverse(
          *universe, GetVar(variables, "ecalrecoil_default")->GetValue(*universe), GetVar(variables, "ehad_true")->GetValue(*universe), event.m_weight);

    if(HasVar(variables,"ecalrecoil_ccpi"))
      GetVar(variables, "ecalrecoil_ccpi") -> m_hists.m_migration.FillUniverse(
          *universe, GetVar(variables, "ecalrecoil_ccpi")->GetValue(*universe), GetVar(variables, "ehad_true")->GetValue(*universe), event.m_weight);

    if(HasVar(variables,"wexp"))
    GetVar(variables, "wexp") -> m_hists.m_migration.FillUniverse(
        *universe, GetVar(variables, "wexp")->GetValue(*universe), GetVar(variables, "wexp_true")->GetValue(*universe), event.m_weight);

    if(HasVar(variables,"etrackrecoil"))
      GetVar(variables, "etrackrecoil") -> m_hists.m_migration.FillUniverse(
          *universe, GetVar(variables, "etrackrecoil")->GetValue(*universe), GetVar(variables, "etracks_true")->GetValue(*universe), event.m_weight);


}
} // namespace run_recoil_energy

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFill(const CCPi::MacroUtil& util, CVUniverse* universe,
                        const EDataMCTruth& type,
                        std::vector<Variable*>& variables) {

  std::cout << "Loop and Fill CutVars\n";
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  std::vector<ECuts> CutsVec = GetCutsVector();
  double EPC = 0.0;
  EventCount counter;
  for (auto c : CutsVec){
    counter[c] = 0.0;
  }

  for(Long64_t i_event=0; i_event < n_entries; ++i_event){
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
    universe->SetEntry(i_event);

    // For mc, get weight, check signal, and sideband
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    // Must pass cuts
    // (Fill event.m_reco_pion_candidate_idxs in the process of checking cuts)
    event.m_passes_cuts = PassesCuts(event, event.m_is_w_sideband);
    if(!event.m_passes_cuts) continue;

    // Get best pion candidate
    event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);

    //// Save pion candidates vector the universe -- needed for Ehad
    //universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);

    if (event.m_passes_cuts){
      EPC = EPC + 1.0;
  //    if (!vtxCut(*universe)){
    //  std::cout << "vtxCut cut is not working" << "\n";
    //  }
    }
    EventCount current_counter = PassedCuts(*universe, event.m_reco_pion_candidate_idxs, is_mc, util.m_signal_definition);
    for (auto c : CutsVec){
      counter[c] = counter[c] + current_counter[c];
    }

    // Fill
    run_recoil_energy::FillVars(event, variables);
  } // events

  for (auto c : CutsVec){
    std::cout << GetCutName(c) << "\t"  << counter[c] << "\n";
  }
  std::cout << "Normal Counter = " << EPC <<"\n";
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runRecoilEnergy(std::string plist = "ME1A") {
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
    const std::string macro("runRecoilEnergy");
    const bool is_grid = false;
    const bool do_truth = false;
    const bool do_systematics = false;

    CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list, plist, do_truth, is_grid, do_systematics);
    util.PrintMacroConfiguration(macro);

  //=========================================
  // Get variables and initialize their hists
  //=========================================
  std::vector<Variable*> variables = run_recoil_energy::GetVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  //=========================================
  // Loop and Fill
  //=========================================
  LoopAndFill(util, util.m_data_universe,              kData, variables);
  LoopAndFill(util, util.m_error_bands.at("cv").at(0), kMC,   variables);

  std::cout << mean << ", " << counts << "\n";
  std::cout << mean/counts << "\n";

  for (auto v : variables) {
    std::string tag = v->Name();
    std::cout << tag << "\n";
    double ymax = -1;
    if(tag == "wexp" || tag == "wexp_true")
      ymax = 1.25;
    if(tag == "etrackrecoil" || tag == "etracks_true")
      ymax = 3.1;
    if(tag == "epi_cal" || tag == "epi_true")
      ymax = 1.7;
    if(tag == "ehad" || tag == "ehad_true")
      ymax = 0.9;
    if(tag == "ecalrecoil_ccinc" || tag == "ecalrecoil_default" || tag == "ecalrecoil_ccpi")
      ymax = 0.9;
    if(tag == "wresid")
      ymax = 1.7e3;
    //if(tag != "epi_aaron" && tag != "epi_ben" && tag != "epi_true")
    //  ymax = 0.725;
    //else
    //  ymax = 2.8;
    bool do_bwn = true;
    std::cout << "Plotting" << std::endl;
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kS),
               util.m_data_pot, util.m_mc_pot, util.m_signal_definition, v->Name(),"SSB", ymax, do_bwn);
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kOtherInt),
               util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag, "FSPart", ymax, do_bwn);
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kRES),
               util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag, "Channel", ymax, do_bwn);
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kPim),
               util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag, "Hadron", ymax, do_bwn);
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kLowW),
               util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag, "Wtrue", ymax, do_bwn);
    if (!v->m_is_true) {
      double zmax = -1;
      if(tag == "ecalrecoilnopi_corr")
        zmax = 50;
      PlotUtils::MnvH2D* mig = (PlotUtils::MnvH2D*)v->m_hists.m_migration.hist->Clone(uniq());
      PlotMigration_AbsoluteBins(mig, v->Name(), zmax);
      PlotMigration_VariableBins(mig, v->Name(), zmax);
    }
  }

  std::cout << "Success" << std::endl;
}

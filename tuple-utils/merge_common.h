#ifndef MERGE_COMMON_H
#define MERGE_COMMON_H

#include "TString.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TList.h"
#include "PlotUtils/MacroUtil.cxx"

#include <iostream>
#include <string>
#include <cassert>

//======================================================================
bool isGoodFile(const char* filename)
{
  TFile f(filename);
  if(f.IsZombie()) return false;
  TTree* meta=(TTree*)f.Get("Meta");
  if(!meta) return false;
  if(!meta->GetBranch("POT_Total")) return false;
  if(!meta->GetBranch("POT_Used")) return false;
  if(!f.Get("Truth")) return false;
  return true;
}

double getTChainPOTXROOTD(TChain& ch){
  double sumPOTUsed=0;
  TObjArray *fileElements=ch.GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    std::string tup(chEl->GetTitle());
    sumPOTUsed += CountPOT(tup);
  }

  return sumPOTUsed;
}

//======================================================================
double getTChainPOT(TChain& ch, const char* branch="POT_Used")
{
  double sumPOTUsed=0;

  TObjArray *fileElements=ch.GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    TFile f(chEl->GetTitle());
    TTree* t=(TTree*)f.Get("Meta");
    if(!t){
      cout << "No Meta tree in file " << chEl->GetTitle() << endl;
      continue;
    }
    assert(t->GetEntries()==1);
    t->GetEntry(0);
    TLeaf* lUsed=t->GetLeaf(branch);
    if(lUsed)         sumPOTUsed+=lUsed->GetValue();
  }

  return sumPOTUsed;
}

//======================================================================

// Get the POT in a Resurrection MC file by parsing the filename,
// since resurrection MC files don't store their POT correctly
double getResurrectionMCFilePOT(const char* filename)
{
  // Parse the filename to work out how many subruns went into this
  // file, and multiply it by the known POT-per-file for MC to get the
  // total POT
  TString basename(gSystem->BaseName(filename));
  // Filenames look like:
  //
  // SIM_minerva_00013200_Subruns_0079-0080-0081-0082-0083_MECAnaTool_Ana_Tuple_v10r8p1.root
  //
  // so first split at underscores, then get the 5th entry and split it at dashes
  TObjArray* split1=basename.Tokenize("_");
  // The 4th item ought to be "Subruns", otherwise something has gone wrong
  TObjString* check=(TObjString*)split1->At(3);
  if(check->GetString()!=TString("Subruns")){
    cout << "Can't parse " << basename << " because item ought to be Subruns but is " << check << endl;
    exit(1);
  }
  TString subrunList=((TObjString*)split1->At(4))->GetString();
  TObjArray* split2=subrunList.Tokenize("-");
  const int nSubruns=split2->GetEntries();
  const double potPerSubrun=1e17;

  delete split1;
  delete split2;

  return nSubruns*potPerSubrun;
}


void setBranchStatuses(TChain& ch)
{
  // In latest ntuples, there aren't any really useless branches, so we'll just keep everything
  return;

  // cout << "Setting branch statuses" << endl;
  // const char* branches[]={
  //   "eventID",
  //   "physEvtNum",
  //   "processType",
  //   "primaryPart",
  //   "n_slices",
  //   "slice_numbers",
  //   "shared_slice",
  //   "vtx",
  //   "vtxErr",
  //   "E",
  //   "found_truth",
  //   "phys_front_activity",
  //   "phys_energy_in_road_upstream_is_rockmuon_consistent",
  //   "rock_muons_removed",
  //   "minos_track_match",
  //   "minos_stub_match",
  //   "unknown_helicity",
  //   "minos_track_inside_partial_plane",
  //   "prim_vtx_has_misassigned_track_direction",
  //   "prim_vtx_has_broken_track",
  //   "pass_MECAna",
  //   "short_track_vtx_used",
  //   "muon_sp_moved",
  //   "vtx_fit_converged",
  //   "muon_is_correct",
  //   "has_int_vtx",
  //   "has_bad_object",
  //   "has_muon",
  //   "muon_has_charge",
  //   "has_good_vtx",
  //   "is_rock_muon",
  //   "is_shared_slice",
  //   "blob_mufuzz_nBlobs",
  //   "blob_mufuzz_nClus",
  //   "blob_mufuzz_nClus_ecal",
  //   "blob_mufuzz_nClus_hcal",
  //   "blob_mufuzz_nClus_nucl",
  //   "blob_mufuzz_nClus_od",
  //   "blob_mufuzz_nClus_tracker",
  //   "blob_recoil_nBlobs",
  //   "blob_recoil_nClus",
  //   "blob_recoil_nClus_ecal",
  //   "blob_recoil_nClus_hcal",
  //   "blob_recoil_nClus_nucl",
  //   "blob_recoil_nClus_od",
  //   "blob_recoil_nClus_tracker",
  //   "blob_vtx_nBlobs",
  //   "blob_vtx_nClus",
  //   "blob_vtx_nClus_ecal",
  //   "blob_vtx_nClus_hcal",
  //   "blob_vtx_nClus_nucl",
  //   "blob_vtx_nClus_od",
  //   "blob_vtx_nClus_tracker",
  //   "broken_track_most_us_plane",
  //   "muon_n_USclusters",
  //   "muon_truthMatch_track_id",
  //   "n_prim_long_tracks",
  //   "n_prim_short_tracks",
  //   "n_start_vertices",
  //   "n_tracks",
  //   "n_tracks_non_prim",
  //   "n_tracks_prim",
  //   "n_tracks_prim_forked",
  //   "n_tracks_prim_kinked",
  //   "n_vertices_startpoint",
  //   "nhits_over_tenmev",
  //   "nhits_over_tenmev_vtx",
  //   "nhits_over_twenty",
  //   "nhits_over_twenty_vtx",
  //   "phys_energy_in_road_downstream_nplanes",
  //   "phys_energy_in_road_upstream_nplanes",
  //   "phys_n_dead_discr_pair",
  //   "phys_n_dead_discr_pair_in_prim_track_region",
  //   "phys_n_dead_discr_pair_two_mod_downstream_prim_track",
  //   "phys_n_dead_discr_pair_two_mod_upstream_prim_vtx",
  //   "phys_n_dead_discr_pair_upstream_prim_track_proj",
  //   "phys_vertex_is_fiducial",
  //   "rock_muon_code",
  //   "truth_has_michel_electron",
  //   "usact_extent",
  //   "usact_n_consecutive",
  //   "usact_n_planes",
  //   "blob_ccqe_recoil_E",
  //   "blob_mufuzz_E",
  //   "blob_mufuzz_E_ecal",
  //   "blob_mufuzz_E_hcal",
  //   "blob_mufuzz_E_nucl",
  //   "blob_mufuzz_E_od",
  //   "blob_mufuzz_E_tracker",
  //   "blob_recoil_E",
  //   "blob_recoil_E_ecal",
  //   "blob_recoil_E_ecal_em",
  //   "blob_recoil_E_ecal_highn",
  //   "blob_recoil_E_ecal_lown",
  //   "blob_recoil_E_ecal_meson",
  //   "blob_recoil_E_ecal_midn",
  //   "blob_recoil_E_ecal_mu",
  //   "blob_recoil_E_ecal_other",
  //   "blob_recoil_E_ecal_p",
  //   "blob_recoil_E_ecal_xtalk",
  //   "blob_recoil_E_em",
  //   "blob_recoil_E_hcal",
  //   "blob_recoil_E_hcal_em",
  //   "blob_recoil_E_hcal_highn",
  //   "blob_recoil_E_hcal_lown",
  //   "blob_recoil_E_hcal_meson",
  //   "blob_recoil_E_hcal_midn",
  //   "blob_recoil_E_hcal_mu",
  //   "blob_recoil_E_hcal_other",
  //   "blob_recoil_E_hcal_p",
  //   "blob_recoil_E_hcal_xtalk",
  //   "blob_recoil_E_highn",
  //   "blob_recoil_E_lown",
  //   "blob_recoil_E_meson",
  //   "blob_recoil_E_midn",
  //   "blob_recoil_E_mu",
  //   "blob_recoil_E_nucl",
  //   "blob_recoil_E_nucl_em",
  //   "blob_recoil_E_nucl_highn",
  //   "blob_recoil_E_nucl_lown",
  //   "blob_recoil_E_nucl_meson",
  //   "blob_recoil_E_nucl_midn",
  //   "blob_recoil_E_nucl_mu",
  //   "blob_recoil_E_nucl_other",
  //   "blob_recoil_E_nucl_p",
  //   "blob_recoil_E_nucl_xtalk",
  //   "blob_recoil_E_od",
  //   "blob_recoil_E_od_em",
  //   "blob_recoil_E_od_highn",
  //   "blob_recoil_E_od_lown",
  //   "blob_recoil_E_od_meson",
  //   "blob_recoil_E_od_midn",
  //   "blob_recoil_E_od_mu",
  //   "blob_recoil_E_od_other",
  //   "blob_recoil_E_od_p",
  //   "blob_recoil_E_od_xtalk",
  //   "blob_recoil_E_other",
  //   "blob_recoil_E_p",
  //   "blob_recoil_E_tracker",
  //   "blob_recoil_E_tracker_em",
  //   "blob_recoil_E_tracker_highn",
  //   "blob_recoil_E_tracker_lown",
  //   "blob_recoil_E_tracker_meson",
  //   "blob_recoil_E_tracker_midn",
  //   "blob_recoil_E_tracker_mu",
  //   "blob_recoil_E_tracker_other",
  //   "blob_recoil_E_tracker_p",
  //   "blob_recoil_E_tracker_xtalk",
  //   "blob_recoil_E_xtalk",
  //   "blob_vtx_E",
  //   "blob_vtx_E_ecal",
  //   "blob_vtx_E_ecal_em",
  //   "blob_vtx_E_ecal_highn",
  //   "blob_vtx_E_ecal_lown",
  //   "blob_vtx_E_ecal_meson",
  //   "blob_vtx_E_ecal_midn",
  //   "blob_vtx_E_ecal_mu",
  //   "blob_vtx_E_ecal_other",
  //   "blob_vtx_E_ecal_p",
  //   "blob_vtx_E_ecal_xtalk",
  //   "blob_vtx_E_em",
  //   "blob_vtx_E_hcal",
  //   "blob_vtx_E_nucl",
  //   "blob_vtx_E_nucl_em",
  //   "blob_vtx_E_nucl_highn",
  //   "blob_vtx_E_nucl_lown",
  //   "blob_vtx_E_nucl_meson",
  //   "blob_vtx_E_nucl_midn",
  //   "blob_vtx_E_nucl_mu",
  //   "blob_vtx_E_nucl_other",
  //   "blob_vtx_E_nucl_p",
  //   "blob_vtx_E_nucl_xtalk",
  //   "blob_vtx_E_od",
  //   "blob_vtx_E_other",
  //   "blob_vtx_E_p",
  //   "blob_vtx_E_tracker",
  //   "blob_vtx_E_tracker_em",
  //   "blob_vtx_E_tracker_highn",
  //   "blob_vtx_E_tracker_lown",
  //   "blob_vtx_E_tracker_meson",
  //   "blob_vtx_E_tracker_midn",
  //   "blob_vtx_E_tracker_mu",
  //   "blob_vtx_E_tracker_other",
  //   "blob_vtx_E_tracker_p",
  //   "blob_vtx_E_tracker_xtalk",
  //   "blob_vtx_E_xtalk",
  //   "energy_at_vtx",
  //   "energy_from_mc",
  //   "energy_from_mc_fraction",
  //   "energy_from_mc_fraction_of_highest",
  //   "lrMuonDistanceFromVertexTracker",
  //   "lrMuonFuzzEnergyEcal",
  //   "lrMuonFuzzEnergyHcal",
  //   "lrMuonFuzzEnergyTracker",
  //   "lrMuonLength",
  //   "lrMuonLengthEcal",
  //   "lrMuonLengthHcal",
  //   "lrMuonLengthTracker",
  //   "mean_time_next_slice",
  //   "mean_time_prev_slice",
  //   "mean_time_this_slice",
  //   "muon_phi",
  //   "muon_theta",
  //   "muon_thetaX",
  //   "muon_thetaY",
  //   "muon_truthMatch_fraction",
  //   "numi_horn_curr",
  //   "numi_pot",
  //   "numi_x",
  //   "numi_x_width",
  //   "numi_y",
  //   "numi_y_width",
  //   "phys_energy_dispersed",
  //   "phys_energy_in_road_downstream",
  //   "phys_energy_in_road_upstream",
  //   "phys_energy_unattached",
  //   "prim_vtx_smallest_opening_angle",
  //   "primary_track_minerva_energy",
  //   "primary_track_minerva_phi",
  //   "primary_track_minerva_theta",
  //   "rik_hcal_recoil",
  //   "rik_hcal_recoil_bad",
  //   "usact_avg_E",
  //   "usact_avg_E_consecutive",
  //   "usact_frac_withE",
  //   "vtx_fit_chi2",
  //   "has_michel_at_vertex",
  //   "has_michel_beginModule",
  //   "has_michel_category_sz",
  //   "has_michel_category",
  //   "has_michel_endModule",
  //   "has_michel_numDigits",
  //   "has_michel_numModules",
  //   "has_michel_numPlanes",
  //   "has_michel_numTracks",
  //   "has_michel_slice_number",
  //   "has_michel_vertex_type",
  //   "has_michel_view_sum",
  //   "has_michel_distance",
  //   "has_michel_energy",
  //   "has_michel_slice_energy",
  //   "has_michel_time_diff",
  //   "muon_corrected_p",
  //   "orig_short_vtx",
  //   "orig_vtx",
  //   "primary_track_minerva_end_position",
  //   "primary_track_minerva_start_position",
  //   "truth_has_michel_from_pion_minus_momentum_sz",
  //   "truth_has_michel_from_pion_minus_momentum",
  //   "truth_has_michel_from_pion_plus_momentum_sz",
  //   "truth_has_michel_from_pion_plus_momentum",
  //   "usact_E_per_plane",
  //   "vtx_blob_radius",
  //   "truth_has_physics_event",
  //   "truth_muon_is_plausible",
  //   "truth_pass_MECAna",
  //   "truth_electronN",
  //   "truth_gammaN",
  //   "truth_in_analyzable_area",
  //   "truth_in_fiducial_area",
  //   "truth_muon_leaving_code",
  //   "truth_muon_track_id",
  //   "truth_neutronN",
  //   "truth_otherN",
  //   "truth_pionN",
  //   "truth_pizeroN",
  //   "truth_protonN",
  //   "truth_reco_has_michel_electron",
  //   "truth_vtx_module",
  //   "truth_vtx_plane",
  //   "truth_electronE",
  //   "truth_gammaE",
  //   "truth_missingE",
  //   "truth_muon_phi",
  //   "truth_muon_theta",
  //   "truth_muon_thetaX",
  //   "truth_muon_thetaY",
  //   "truth_neutronE",
  //   "truth_nu",
  //   "truth_otherE",
  //   "truth_pionE",
  //   "truth_pizeroE",
  //   "truth_protonE",
  //   "truth_q3",
  //   "truth_visibleE",
  //   "truth_muon_end_momentum",
  //   "truth_muon_end_position",
  //   "MECAna_nuFlavor",
  //   "MECAna_nuHelicity",
  //   "MECAna_intCurrent",
  //   "MECAna_intType",
  //   "MECAna_E",
  //   "MECAna_Q2",
  //   "MECAna_x",
  //   "MECAna_y",
  //   "MECAna_W",
  //   "MECAna_score",
  //   "MECAna_leptonE",
  //   "MECAna_vtx",
  //   "MECAna_minos_trk_is_contained",
  //   "MECAna_minos_trk_is_ok",
  //   "MECAna_minos_used_range",
  //   "MECAna_minos_used_curvature",
  //   "MECAna_pass_canonical_cut",
  //   "MECAna_is_cc",
  //   "MECAna_in_analyzable_area",
  //   "MECAna_in_fiducial_area",
  //   "MECAna_minos_trk_end_plane",
  //   "MECAna_minos_trk_quality",
  //   "MECAna_r_minos_trk_vtx_plane",
  //   "MECAna_t_minos_trk_numFSMuons",
  //   "MECAna_t_minos_trk_primFSLeptonPDG",
  //   "MECAna_vtx_module",
  //   "MECAna_vtx_plane",
  //   "MECAna_E_ccqe",
  //   "MECAna_E_wide_window",
  //   "MECAna_Q2_ccqe",
  //   "MECAna_Q2_wide_window",
  //   "MECAna_W_wide_window",
  //   "MECAna_minos_trk_bave",
  //   "MECAna_minos_trk_chi2",
  //   "MECAna_minos_trk_end_u",
  //   "MECAna_minos_trk_end_v",
  //   "MECAna_minos_trk_end_x",
  //   "MECAna_minos_trk_end_y",
  //   "MECAna_minos_trk_end_z",
  //   "MECAna_minos_trk_eqp",
  //   "MECAna_minos_trk_eqp_qp",
  //   "MECAna_minos_trk_fit_pass",
  //   "MECAna_minos_trk_ndf",
  //   "MECAna_minos_trk_p",
  //   "MECAna_minos_trk_p_curvature",
  //   "MECAna_minos_trk_p_range",
  //   "MECAna_minos_trk_qp",
  //   "MECAna_minos_trk_vtx_x",
  //   "MECAna_minos_trk_vtx_y",
  //   "MECAna_minos_trk_vtx_z",
  //   "MECAna_nu_energy_recoil",
  //   "MECAna_q3",
  //   "MECAna_r_minos_trk_bdL",
  //   "MECAna_r_minos_trk_end_dcosx",
  //   "MECAna_r_minos_trk_end_dcosy",
  //   "MECAna_r_minos_trk_end_dcosz",
  //   "MECAna_r_minos_trk_vtx_dcosx",
  //   "MECAna_r_minos_trk_vtx_dcosy",
  //   "MECAna_r_minos_trk_vtx_dcosz",
  //   "MECAna_recoil_E",
  //   "MECAna_recoil_E_wide_window",
  //   "MECAna_vtx_apothem",
  //   "MECAna_x_wide_window",
  //   "MECAna_y_wide_window",
  //   "MECAna_muon_vtx",
  //   "MECAna_sys_muon_curve_energy_shift",
  //   "MECAna_sys_muon_energy_shift",
  //   "MECAna_sys_muon_minerva_energy_shift",
  //   "MECAna_sys_muon_qSquared_shift",
  //   "MECAna_sys_muon_range_energy_shift",
  //   "MECAna_sys_muon_wSquared_shift",
  //   "MECAna_sys_muon_xbj_shift",
  //   "MECAna_sys_muon_y_shift",
  //   "MECAna_sys_nu_energy_shift",
  //   "MECAna_sys_recoil_energy_shift",
  //   "MECAna_sys_recoil_qSquared_shift",
  //   "MECAna_sys_recoil_wSquared_shift",
  //   "MECAna_sys_recoil_xbj_shift",
  //   "MECAna_sys_recoil_y_shift",
  //   "MECAna_sys_total_qSquared_shift",
  //   "MECAna_sys_total_wSquared_shift",
  //   "MECAna_sys_total_xbj_shift",
  //   "MECAna_sys_total_y_shift",
  //   "ev_run",
  //   "ev_subrun",
  //   "ev_detector",
  //   "ev_triggerType",
  //   "ev_gate",
  //   "ev_global_gate",
  //   "ev_gps_time_sec",
  //   "ev_gps_time_usec",
  //   "mc_run",
  //   "mc_subrun",
  //   "mc_nInteractions",
  //   "mc_MIState",
  //   "mc_pot",
  //   "mc_beamConfig",
  //   "mc_processType",
  //   "mc_nthEvtInSpill",
  //   "mc_nthEvtInFile",
  //   "mc_intType",
  //   "mc_current",
  //   "mc_charm",
  //   "mc_weight",
  //   "mc_XSec",
  //   "mc_diffXSec",
  //   "mc_incoming",
  //   "mc_fluxDriverProb",
  //   "mc_targetNucleus",
  //   "mc_targetZ",
  //   "mc_targetA",
  //   "mc_targetNucleon",
  //   "mc_struckQuark",
  //   "mc_seaQuark",
  //   "mc_resID",
  //   "mc_primaryLepton",
  //   "mc_incomingE",
  //   "mc_Bjorkenx",
  //   "mc_Bjorkeny",
  //   "mc_Q2",
  //   "mc_nuT",
  //   "mc_w",
  //   "mc_vtx",
  //   "mc_incomingPartVec",
  //   "mc_initNucVec",
  //   "mc_primFSLepton",
  //   "mc_nFSPart",
  //   "mc_FSPartPx",
  //   "mc_FSPartPy",
  //   "mc_FSPartPz",
  //   "mc_FSPartE",
  //   "mc_FSPartPDG",
  //   "mc_cvweight_total",
  //   "wgt",
  //   "mc_cvweight_totalFlux",
  //   "mc_cvweight_totalXsec",
  //   "mc_cvweight_NA49",
  //   "genie_wgt_n_shifts",
  //   "truth_genie_wgt_AGKYxF1pi",
  //   "truth_genie_wgt_AhtBY",
  //   "truth_genie_wgt_BhtBY",
  //   "truth_genie_wgt_CCQEPauliSupViaKF",
  //   "truth_genie_wgt_CV1uBY",
  //   "truth_genie_wgt_CV2uBY",
  //   "truth_genie_wgt_EtaNCEL",
  //   "truth_genie_wgt_FrAbs_N",
  //   "truth_genie_wgt_FrAbs_pi",
  //   "truth_genie_wgt_FrCEx_N",
  //   "truth_genie_wgt_FrCEx_pi",
  //   "truth_genie_wgt_FrElas_N",
  //   "truth_genie_wgt_FrElas_pi",
  //   "truth_genie_wgt_FrInel_N",
  //   "truth_genie_wgt_FrInel_pi",
  //   "truth_genie_wgt_FrPiProd_N",
  //   "truth_genie_wgt_FrPiProd_pi",
  //   "truth_genie_wgt_MFP_N",
  //   "truth_genie_wgt_MFP_pi",
  //   "truth_genie_wgt_MaCCQE",
  //   "truth_genie_wgt_MaCCQEshape",
  //   "truth_genie_wgt_MaNCEL",
  //   "truth_genie_wgt_MaRES",
  //   "truth_genie_wgt_MvRES",
  //   "truth_genie_wgt_NormCCQE",
  //   "truth_genie_wgt_NormCCRES",
  //   "truth_genie_wgt_NormDISCC",
  //   "truth_genie_wgt_NormNCRES",
  //   "truth_genie_wgt_RDecBR1gamma",
  //   "truth_genie_wgt_Rvn1pi",
  //   "truth_genie_wgt_Rvn2pi",
  //   "truth_genie_wgt_Rvp1pi",
  //   "truth_genie_wgt_Rvp2pi",
  //   "truth_genie_wgt_Theta_Delta2Npi",
  //   "truth_genie_wgt_VecFFCCQEshape",
  //   "truth_genie_wgt_shifts",
  //   0
  // }; 

  // ch.SetBranchStatus("*", 0);
  // int i=0;
  // while(branches[i]){
  //   if(ch.GetBranch(branches[i]) || ch.FindBranch(branches[i])){
  //     ch.SetBranchStatus(branches[i], 1);
  //   }
  //   ++i;
  // }
  // cout << "Done setting branch statuses" << endl;
}


#endif

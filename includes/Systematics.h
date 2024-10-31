#ifndef Systematics_h
#define Systematics_h
//=============================================================================
// Get container of systematics
// Reminder: Universes own a chain. This is so one can get values within a
// a universe.
// If you're not accessing a universe's values, the chain isn't needed. This is
// the case when I load an MnvH1D from a file into a HistWrapper. You can
// instead pass a dummy chain.
// You don't need a real CHW to make the universes and use them with a HW.
//=============================================================================
#include "CVUniverse.h"
#include "CohDiffractiveSystematics.h"
#include "Constants.h"           // typedefs UniverseMap
#include "LateralSystematics.h"  // Detector systematics
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/GeantHadronSystematics.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/MichelSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/TargetMassSystematics.h"

namespace systematics {
const std::vector<std::string> kGenieSystematics_FSI_nucleons = {
    //   nucleon
    // Hadronization model
    "GENIE_FrAbs_N",     // Absorption
    "GENIE_FrCEx_N",     // Charge exchange
    "GENIE_FrElas_N",    // Elastic
    "GENIE_FrInel_N",    // Inelastic
    "GENIE_FrPiProd_N",  // Pion produciton
    "GENIE_MFP_N"        // mean free paths
};

const std::vector<std::string> kGenieSystematics_FSI_pions = {
    "GENIE_AGKYxF1pi",     // Hadronization model
    "GENIE_FrAbs_pi",      // Absorption
    "GENIE_FrCEx_pi",      // Charge exchange
    "GENIE_FrElas_pi",     // Elastic
                           // Inelastic
    "GENIE_FrPiProd_pi",   // Pion produciton
    "GENIE_MFP_pi",        // mean free paths
    "GENIE_RDecBR1gamma",  // Resonant decay photon branching ratio
    // "GENIE_Theta_Delta2Npi"  // anisotropic delta decay (BROKEN)
};

const std::vector<std::string> kGenieSystematics_InteractionModel = {
    "GENIE_AhtBY",
    "GENIE_BhtBY",
    "GENIE_CV1uBY",
    "GENIE_CV2uBY",  // Bodek-yank params
    "GENIE_MaNCEL",
    "GENIE_EtaNCEL",  // masses/form factors
    "GENIE_MaRES",
    "GENIE_MvRES",  // masses
    "GENIE_NormDISCC",
    "GENIE_NormNCRES",          // norm
    "GENIE_VecFFCCQEshape",     // shapes
    "GENIE_CCQEPauliSupViaKF",  // pauli suppression
    "GENIE_D2_MaRES",
    "GENIE_EP_MvRES",
    "GENIE_Theta_Delta2Npi",
    "GENIE_D2_NormCCRES",
    "GENIE_MaCCQE"};

UniverseMap GetSystematicUniversesMap(PlotUtils::ChainWrapper* chain,
                                      bool is_truth = false,
                                      bool do_full_systematics = false) {
  // return map
  UniverseMap error_bands;

  // CV
  error_bands[std::string("cv")].push_back(new CVUniverse(chain));

  // Systematics
  if (do_full_systematics) {
    //========================================================================
    // DETECTOR
    //========================================================================
    std::vector<double> sigmas = {-1., +1.};
    for (auto sigma : sigmas) {
      error_bands[std::string("Birks")].push_back(
          new BirksShiftUniverse(chain, sigma));

      error_bands[std::string("BB")].push_back(
          new BetheBlochShiftCVUniverse(chain, sigma));

      error_bands[std::string("DetMass")].push_back(
          new DetectorMassShiftCVUniverse(chain, sigma));

      error_bands[std::string("TrackAngle")].push_back(
          new TrackAngleShiftCVUniverse(chain, sigma));

      error_bands[std::string("BeamAngle")].push_back(
          new BeamAngleShiftCVUniverse(chain, sigma));

      error_bands[std::string("NodeCutEff")].push_back(
          new NodeCutEffUniverse(chain, sigma));

      error_bands[std::string("Target_Mass_CH")].push_back(
          new PlotUtils::TargetMassScintillatorUniverse<CVUniverse>(chain,
                                                                    sigma));
    }

    UniverseMap geant_bands =
        PlotUtils::GetGeantHadronSystematicsMap<CVUniverse>(chain);
    error_bands.insert(geant_bands.begin(), geant_bands.end());

    //========================================================================
    // FLUX
    //========================================================================
    UniverseMap bands_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(
        chain, CCNuPionIncConsts::kNFluxUniverses);
    error_bands.insert(bands_flux.begin(), bands_flux.end());

    //========================================================================
    // GENIE
    //========================================================================
    UniverseMap genie_error_bands =
        PlotUtils::GetGenieSystematicsMap<CVUniverse>(
            chain, false);  // Not including the new fitted values
    error_bands.insert(genie_error_bands.begin(), genie_error_bands.end());

    // New GENIE MaRES and NormCCRes error bands No Covariance
    UniverseMap new_res_genie_error_bands =
        PlotUtils::GetGenieResPionFitSystematicsMap<CVUniverse>(chain);
    error_bands.insert(new_res_genie_error_bands.begin(),
                       new_res_genie_error_bands.end());

    // New GENIE MvRES
    UniverseMap new_ep_genie_error_bands =
        PlotUtils::GetGenieEPMvResSystematicsMap<CVUniverse>(chain);
    error_bands.insert(new_ep_genie_error_bands.begin(),
                       new_ep_genie_error_bands.end());

    //========================================================================
    // MnvTunes
    //========================================================================
    // RPA
    UniverseMap bands_rpa = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
    error_bands.insert(bands_rpa.begin(), bands_rpa.end());

    // 2P2H
    UniverseMap bands_2p2h =
        PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
    error_bands.insert(bands_2p2h.begin(), bands_2p2h.end());

    //// LowQ2Pi
    // GetLowQ2PiSystematicsMap(typename T::config_t chain);
    // std::vector<CVUniverse*> error_bands_lowq2pi =
    // PlotUtils::GetLowQ2PiSystematics<CVUniverse>(chain);
    // error_bands[std::string("LowQ2Pi")] = error_bands_lowq2pi;
    //    UniverseMap error_bands_lowq2pi =
    //        PlotUtils::GetLowQ2PiSystematicsMap<CVUniverse>(chain);
    //    error_bands.insert(error_bands_lowq2pi.begin(),
    //    error_bands_lowq2pi.end());

    //========================================================================
    // Angle Systematics
    //========================================================================
    ////Angle systematics
    UniverseMap angle_error_bands =
        PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
    error_bands.insert(angle_error_bands.begin(), angle_error_bands.end());

    //========================================================================
    // Muons
    //========================================================================
    ////Muon Angle systematics
    UniverseMap muon_angle_error_bands =
        PlotUtils::GetMuonAngleResolutionSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muon_angle_error_bands.begin(),
                       muon_angle_error_bands.end());

    //// MUON - MINERvA
    UniverseMap bands_muon_minerva =
        PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
    error_bands.insert(bands_muon_minerva.begin(), bands_muon_minerva.end());

    //// MUON - MINOS
    UniverseMap bands_muon_minos =
        PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
    error_bands.insert(bands_muon_minos.begin(), bands_muon_minos.end());

    //// MINOS EFFICIENCY
    UniverseMap bands_minoseff =
        PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
    error_bands.insert(bands_minoseff.begin(), bands_minoseff.end());

    UniverseMap muon_res_error_bands =
        PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muon_res_error_bands.begin(),
                       muon_res_error_bands.end());

    //========================================================================
    // Particle Response
    //========================================================================
    const bool use_neutron = false;
    const bool use_new = true;
    UniverseMap bands_response =
        PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain, use_neutron,
                                                         use_new);
    error_bands.insert(bands_response.begin(), bands_response.end());

    //========================================================================
    // Michel Efficiency Error bands
    //========================================================================
    UniverseMap michel_error_bands =
        PlotUtils::GetMichelEfficiencySystematicsMap<CVUniverse>(chain);
    error_bands.insert(michel_error_bands.begin(), michel_error_bands.end());

    //========================================================================
    // Tpi estimator for untracked pis.
    //========================================================================
    // Due to uncertainty on fit params.
    UniverseMap tpi_estimator_bands =
        PlotUtils::GetTpiMichelRangeEstimatorSystematicsMap<CVUniverse>(chain);
    error_bands.insert(tpi_estimator_bands.begin(), tpi_estimator_bands.end());

    // Due to uncertainty on michel range from scintillator mass
    UniverseMap mich_range_bands =
        PlotUtils::GetMichelRangeDetectorMassSystematicsMap<CVUniverse>(chain);
    error_bands.insert(mich_range_bands.begin(), mich_range_bands.end());

    // Due  to uncertainty on Tpi weight
    // UniverseMap bands_UntrackedPion =
    //	   PlotUtils::GetUntrackedPionSystematicsMap<CVUniverse>(chain);
    // error_bands.insert(bands_UntrackedPion.begin(),
    // bands_UntrackedPion.end());

    /*    std::vector<CVUniverse*> bands_UntrackedPion =
            PlotUtils::GetUntrackedPionSystematics<CVUniverse>(chain);
        error_bands[std::string("UntrackedPi")] = bands_UntrackedPion;
    */
    //========================================================================
    // One systematic to cover both lowq2 weight and mehreen's tpi weight
    //========================================================================
    std::vector<CVUniverse*> bands_pi_tunes =
        PlotUtils::GetChargedPionTuneSystematics<CVUniverse>(chain);
    error_bands[std::string("pitune")] = bands_pi_tunes;

    //========================================================================
    // Diffractive pion production unc
    //========================================================================
    UniverseMap error_bands_cohdiff = GetCohDiffractiveSystematicsMap(chain);
    error_bands.insert(error_bands_cohdiff.begin(), error_bands_cohdiff.end());
  }

  for (auto band : error_bands) {
    std::vector<CVUniverse*> universes = band.second;
    for (auto universe : universes) universe->SetTruth(is_truth);
  }

  return error_bands;
}
}  // namespace systematics

#endif  // Systematics_h

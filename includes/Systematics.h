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
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MuonSystematics.h"
//#include "PlotUtils/ResponseSystematics.h"

#include "Constants.h" // typedefs UniverseMap
#include "CVUniverse.h"
#include "LateralSystematics.h" // Detector systematics

namespace systematics {
  const std::vector<std::string> kGenieSystematics_FSI_nucleons = {
     //   nucleon      
                        // Hadronization model
    "GENIE_FrAbs_N",    // Absorption
    "GENIE_FrCEx_N",    // Charge exchange
    "GENIE_FrElas_N",   // Elastic
    "GENIE_FrInel_N",   // Inelastic
    "GENIE_FrPiProd_N", // Pion produciton
    "GENIE_MFP_N"       // mean free paths
  };

  const std::vector<std::string> kGenieSystematics_FSI_pions = {
    "GENIE_AGKYxF1pi",       // Hadronization model
    "GENIE_FrAbs_pi",        // Absorption
    "GENIE_FrCEx_pi",        // Charge exchange
    "GENIE_FrElas_pi",       // Elastic
                             // Inelastic
    "GENIE_FrPiProd_pi",     // Pion produciton
    "GENIE_MFP_pi",          // mean free paths
    "GENIE_RDecBR1gamma",    // Resonant decay photon branching ratio
    "GENIE_Theta_Delta2Npi"  // anisotropic delta decay (BROKEN)
  };

  const std::vector<std::string> kGenieSystematics_InteractionModel = {
    "GENIE_AhtBY",     "GENIE_BhtBY",     "GENIE_CV1uBY",    "GENIE_CV2uBY", // Bodek-yank params
    "GENIE_MaNCEL",    "GENIE_EtaNCEL",                                      // masses/form factors
    "GENIE_MaRES",     "GENIE_MvRES",                                        // masses
    "GENIE_NormDISCC", "GENIE_NormNCRES",                                    // norm
    "GENIE_VecFFCCQEshape",                                                  // shapes
    "GENIE_CCQEPauliSupViaKF"                                                // pauli suppression
  };

  UniverseMap GetSystematicUniversesMap(PlotUtils::ChainWrapper* chain,
                                        bool is_truth = false,
                                        bool do_full_systematics = false) {
    // return map
    UniverseMap error_bands;

    // CV
    error_bands[std::string("cv")].push_back(new CVUniverse(chain));

    // Systematics
    if(do_full_systematics){
      //========================================================================
      // DETECTOR
      //========================================================================
      std::vector<double> sigmas = {-1., +1.};
      for (auto sigma : sigmas){
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
      }

      //========================================================================
      // FLUX
      //========================================================================
      UniverseMap bands_flux = 
          PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain, CCNuPionIncConsts::kNFluxUniverses);
      error_bands.insert(bands_flux.begin(), bands_flux.end());

      //========================================================================
      // GENIE
      //========================================================================
      // Standard
      UniverseMap bands_genie =
          PlotUtils::GetStandardGenieSystematicsMap<CVUniverse>(chain);
      error_bands.insert(bands_genie.begin(), bands_genie.end());

      // Pion final state normalization
      UniverseMap bands_pi_fs_norm =
          PlotUtils::GetGenieRvx1piSystematicsMap<CVUniverse>(chain);
      error_bands.insert(bands_pi_fs_norm.begin(), bands_pi_fs_norm.end());

      //========================================================================
      // MnvTunes
      //========================================================================
      // RPA
      UniverseMap bands_rpa =
          PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
      error_bands.insert(bands_rpa.begin(), bands_rpa.end());

      // 2P2H
      UniverseMap bands_2p2h =
          PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
      error_bands.insert(bands_2p2h.begin(), bands_2p2h.end());

      //// LowQ2Pi 
      ////GetLowQ2PiSystematicsMap(typename T::config_t chain);
      //std::vector<CVUniverse*> error_bands_lowq2pi = PlotUtils::GetLowQ2PiSystematics<CVUniverse>(chain);
      //error_bands[std::string("LowQ2Pi")] = error_bands_lowq2pi;

      //========================================================================
      // Muons
      //========================================================================
      //// MUON - MINERvA
      //UniverseMap bands_muon_minerva =
      //    PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
      //error_bands.insert(bands_muon_minerva.begin(), bands_muon_minerva.end());

      //// MUON - MINOS
      //UniverseMap bands_muon_minos =
      //    PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
      //error_bands.insert(bands_muon_minos.begin(), bands_muon_minos.end());

      //// MINOS EFFICIENCY
      //UniverseMap bands_minoseff =
      //    PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
      //error_bands.insert(bands_minoseff.begin(), bands_minoseff.end());

      //========================================================================
      // Particle Response
      //========================================================================
      //      const bool use_neutron = false;
      //      const bool use_new = true;
      //    UniverseMap bands_response = 
      //      PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain , use_neutron, use_new/*, proton*/);
      // error_bands.insert(bands_response.begin(), bands_response.end());

    }


    for(auto band : error_bands){
      std::vector<CVUniverse*> universes = band.second;
      for(auto universe : universes) universe->SetTruth(is_truth);
    }
    return error_bands;
  }
}

#endif // Systematics_h

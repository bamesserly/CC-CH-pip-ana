#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include <cstdlib>

#include "GENIEXSecExtract/XSecLooper.h"
#include "PlotUtils/FluxReweighter.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
typedef unsigned int uint;

class MinModDepCCQEXSec : public XSec {
 public:
  MinModDepCCQEXSec(const char* name) : XSec(name){};

  TVector3 GetPmuVector(ChainWrapper& chw, int entry) {  // pmu vector in GeV
    TVector3 pmuVec((double)chw.GetValue("mc_primFSLepton", entry, 0) / 1000.,
                    (double)chw.GetValue("mc_primFSLepton", entry, 1) / 1000.,
                    (double)chw.GetValue("mc_primFSLepton", entry, 2) / 1000.);
    double numi_beam_angle_rad = -0.05887;
    pmuVec.RotateX(numi_beam_angle_rad);
    return pmuVec;
  }

  double thetamudegrees(ChainWrapper& chw,
                        int entry) {  // it returns the value
                                      // of thetamu in degrees
    TVector3 pmuVec = GetPmuVector(chw, entry);
    double thetamu = pmuVec.Theta() * 180. / 3.141592654;
    return thetamu;
  }

  double GetPmu(ChainWrapper& chw, int entry) {  // return pmu in GeV
    TVector3 pmuVec = GetPmuVector(chw, entry);
    return pmuVec.Mag();
  }

  double CalcWexp(double Q2, double Ehad) {  // return Wexp in GeV
    double W = pow(0.9383, 2.0) - Q2 + 2.0 * (0.9383) * Ehad;
    W = W > 0 ? sqrt(W) : 0.0;
    return W;
  }

  double GetEnuTrue(ChainWrapper& chw, int entry) {  // return Enu in GeV
    return (double)chw.GetValue("mc_incomingE", entry) / 1000.;
  }

  double GetQ2True(ChainWrapper& chw, int entry) {  // return Q2 in GeV^2
    return (double)chw.GetValue("mc_Q2", entry) / 1.e6;
  }

  int GetNChargedPionsTrue(ChainWrapper& chw, int entry) {
    return (int)chw.GetValue("truth_N_pip", entry) +
           (int)chw.GetValue("truth_N_pim", entry);
  }

  double GetTpi(ChainWrapper& chw, int entry, int idx) {  // return Tpi in GeV
    double t_pi_E = (double)chw.GetValue("truth_pi_E", entry, idx);
    if (t_pi_E == -1.) {
      std::cerr << "CVU::GetTpi: Default energy.\n"
                   "Tried to access truth pion energy for a nonexistent "
                   "truth pion trajectory.\n";
      throw t_pi_E;
    }
    return (t_pi_E - 139.569) / 1000;
  }

  std::vector<double> GetTpiTrueVec(ChainWrapper& chw, int entry) {
    std::vector<double> ret;
    const int n_true_pions = GetNChargedPionsTrue(chw, entry);
    for (int idx = 0; idx < n_true_pions; ++idx) {
      ret.push_back(GetTpi(chw, entry, idx));
    }
    return ret;
  }

  int GetPiCharge(ChainWrapper& chw, int entry, int idx) const {
    int t_pi_charge = (int)chw.GetValue("truth_pi_charge", entry, idx);
    if (t_pi_charge == 0) {
      std::cerr << "CVU::GetPiCharge: Default charge.\n"
                   "Tried to access truth pion charge for a nonexistent "
                   "truth pion trajectory.\n";
      throw t_pi_charge;
    }
    return t_pi_charge;
  }

  int GetHighestpionEnergyIdx(ChainWrapper& chw,
                              int entry) {  // returns the index for the
    // pion track with the hights
    // energy
    std::vector<double> tpi_vec = GetTpiTrueVec(chw, entry);  // pip and pim
    const int n_true_pions = GetNChargedPionsTrue(chw, entry);

    int reigning_idx = -1;
    double reigning_tpi = 0;
    for (int idx = 0; idx < n_true_pions; ++idx) {
      if (tpi_vec[idx] > reigning_tpi && GetPiCharge(chw, entry, idx) > 0.) {
        reigning_idx = idx;
        reigning_tpi = tpi_vec[idx];
      }
    }
    return reigning_idx;
  }

  bool OnePionEvt(ChainWrapper& chw, int entry) {
    bool pass = true;
    pass = pass && (int)chw.GetValue("truth_N_pi0", entry) == 0;
    pass = pass && (int)chw.GetValue("truth_N_pim", entry) == 0;
    pass = pass && (int)chw.GetValue("truth_N_pip", entry) == 1;
    return pass;
  }

  int NOtherParticles(ChainWrapper& chw, int entry) {
    int n_other_particles = 0;
    n_other_particles = (int)chw.GetValue("truth_N_chargedK", entry) +
                        (int)chw.GetValue("truth_N_K0", entry) +
                        (int)chw.GetValue("truth_N_sigma", entry) +
                        (int)chw.GetValue("truth_N_lambda", entry);
    return n_other_particles;
  }

  double GetElepTrue(ChainWrapper& chw, int entry) {  // return Emu in GeV
    return (double)chw.GetValue("mc_primFSLepton", entry, 3) / 1000;
  }

  double GetEhadTrue(ChainWrapper& chw, int entry) {  // return Ehad in GeV
    return GetEnuTrue(chw, entry) - GetElepTrue(chw, entry);
  }

  double GetWexpTrue(ChainWrapper& chw, int entry) {  // return Wexp in GeV
    return CalcWexp(GetQ2True(chw, entry), GetEhadTrue(chw, entry));
  }

  int NSignalPions(ChainWrapper& chw, int entry) {
    int n_signal_pions = 0;
    int n_true_pions = GetNChargedPionsTrue(chw, entry);
    for (int idx = 0; idx < n_true_pions; ++idx) {
      double t_pi = GetTpi(chw, entry, idx);
      // double theta_pi = GetThetapiTrue(idx);
      if (GetPiCharge(chw, entry, idx) > 0 && t_pi > 0.035 && t_pi < 0.350
          //&& (theta_pi < 1.39626 || 1.74533 < theta_pi))
      )
        ++n_signal_pions;
    }
    return n_signal_pions;
  }

  bool leftlinesCut(const double a, const double x, const double y) {
    double b, yls, yli;
    b = a * (2 * sqrt(3) / 3);
    yls = (sqrt(3) / 3) * x + b;
    yli = -(sqrt(3) / 3) * x - b;
    if (y > yli && y < yls)
      return true;
    else
      return false;
  }

  bool rightlinesCut(const double a, const double x, const double y) {
    double b, yls, yli;
    b = a * (2 * sqrt(3) / 3);
    yls = -(sqrt(3) / 3) * x + b;
    yli = (sqrt(3) / 3) * x - b;
    if (y > yli && y < yls)
      return true;
    else
      return false;
  }

  bool zVertexSig(ChainWrapper& chw, int entry) {
    double vtxZ = (double)chw.GetValue("mc_vtx", entry, 2);
    if (vtxZ > 5990.0 && vtxZ < 8340.0)
      return true;
    else
      return false;
  }
  bool XYVertexSig(ChainWrapper& chw, int entry) {
    const double a = 850.0;
    const double x = (double)chw.GetValue("mc_vtx", entry, 0),
                 y = (double)chw.GetValue("mc_vtx", entry, 1);
    if (x < 0) {
      if (x > -a && leftlinesCut(a, x, y))
        return true;
      else
        return false;
    } else {
      if (x < a && rightlinesCut(a, x, y))
        return true;
      else
        return false;
    }
  }

  bool VtxSignal(ChainWrapper& chw, int entry) {
    bool Pass = true;
    Pass = Pass && zVertexSig(chw, entry);
    Pass = Pass && XYVertexSig(chw, entry);
    return Pass;
  }

  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry) {
    int pionIdx = GetHighestpionEnergyIdx(chw, entry);
    double pmu = GetPmu(chw, entry);
    double Wexp = GetWexpTrue(chw, entry);
    int Npions = NSignalPions(chw, entry);
    if ((int)chw.GetValue("mc_current", entry) != 1) return false;
    if (!chw.GetValue("truth_is_fiducial", entry)) return false;
    if (!VtxSignal(chw, entry)) return false;
    if ((int)chw.GetValue("mc_incoming", entry) != 14) return false;
    if (thetamudegrees(chw, entry) > 20.) return false;
    if (Wexp < 0.) return false;
    if (Wexp > 1.4) return false;
    if (Npions != 1) return false;
    if (NOtherParticles(chw, entry) > 0) return false;
    if (pmu < 1.5 || 20 < pmu) return false;
    if (chw.GetValue("truth_N_pi0", entry) != 0) return false;
    if (chw.GetValue("truth_N_pim", entry) != 0) return false;
    return true;
  }
};

int runXSecLooper() {
  TH1::AddDirectory(kFALSE);  // Needed so that MnvH1D gets to clean up its own
                              // MnvLatErrorBands (which are TH1Ds).

  // const std::string playlistFile =
  //    "/minerva/app/users/granados/cmtuser/MINERvA101/"
  //    "MINERvA-101-Cross-Section/MCME1A.txt";

  // shorter playlist for testing
  const std::string playlistFile = "MCME1A.txt";

  // Create the XSecLooper and tell it the input files
  // Inputs should be the merged ntuples:
  XSecLooper loop(playlistFile.c_str());

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);

  // Setting the number of Universes in the GENIE error band (default 100, put 0
  // if you do not want to include the universes)
  loop.setNumUniv(0);
  loop.setFiducial(5990, 8340);
  // loop.setFiducial(5890, 8467);
  loop.setPlaylist(PlotUtils::FluxReweighter::minervame1A);
  // Add the differential cross section dsigma/ds_dpT
  double pmu_edges[] = {0., 1., 2., 3., 4., 5.5, 7.5, 10., 13., 20., 30.};
  int pmu_nbins = 10;

  std::cout << "pmu \n";

  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpmu = new MinModDepCCQEXSec("pmu");
  ds_dpmu->setBinEdges(pmu_nbins, pmu_edges);
  ds_dpmu->setDimension(1);
  ds_dpmu->setFluxIntLimits(0.0, 100.0);
  ds_dpmu->setIsFluxIntegrated(true);
  ds_dpmu->setNormalizationType(XSec::kPerNucleon);
  ds_dpmu->setUniverses(0);  // default value, put 0 if you do not want
                             // universes to be included.
  ds_dpmu->setVariable(XSec::kPLep);

  loop.addXSec(ds_dpmu);

  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =
      "GENIEXSECEXTRACT_" +
      playlistFile.substr(playlistFile.rfind("/") + 1, playlistFile.find(".")) +
      ".root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  std::cout << "getXSecs vector size " << loop.getXSecs().size() << "\n";

  PlotUtils::MnvH1D* flux = (PlotUtils::MnvH1D*)loop.getFluxHist();

  for (uint i = 0; i < loop.getXSecs().size(); ++i) {
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }

  PlotUtils::FluxReweighter* fluxReweighter = new PlotUtils::FluxReweighter(
      14, true, PlotUtils::FluxReweighter::minervame1A,
      PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6);

  PlotUtils::MnvH1D* Integrated_flux =
      fluxReweighter->GetIntegratedFluxReweighted_FromInputFlux(
          flux, loop.getXSecs()[0]->getXSecHist(), 0.0, 100.0);

  PlotUtils::MnvH1D* integrated_flux_other =
      fluxReweighter->GetIntegratedFluxReweighted(
          14, loop.getXSecs()[0]->getXSecHist(), 0., 100.);

  Integrated_flux->Write();

  PlotUtils::MnvH1D* unfolded =
      (PlotUtils::MnvH1D*)loop.getXSecs()[0]->getXSecHist()->Clone("unfolded");
  unfolded->ClearAllErrorBands();
  unfolded->Reset();
  unfolded->Multiply(loop.getXSecs()[0]->getXSecHist(), Integrated_flux);

  unfolded->Write();
  PlotUtils::MnvH1D* selected = ds_dpmu->getXSecHist()->Clone("selected");
  selected->Write();
  std::cout << "Flag 4\n";

  return 0;
}

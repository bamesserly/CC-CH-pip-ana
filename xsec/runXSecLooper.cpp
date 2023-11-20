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

  bool Is1PiPlus(const std::map<std::string, int>& particles) {
    TString counted_signal("muons piplus piplus_range pions mesons");
    TString any_baryon("protons neutrons nucleons heavy_baryons nuclei");
    int n_bg = 0;
    for (auto p : particles) {
      std::string particle = p.first;
      int count = p.second;
      if (counted_signal.Contains(particle) || any_baryon.Contains(particle))
        continue;
      n_bg += count;
    }
    if (particles.at("muons") == 1 && particles.at("piplus") == 1 &&
        particles.at("pions") == particles.at("piplus") &&
        particles.at("mesons") == particles.at("piplus") && n_bg == 0)
      return true;  // any number of baryons, nothing else

    return false;
  }

  std::map<string, int> GetParticleTopology(
      const std::vector<int>& FS_PDG, const std::vector<double>& FS_energy) {
    std::map<std::string, int> genie_n;

    // Overarching categories: nucleons, mesons
    genie_n["muons"] = 0;     // Muons, photons (do we want electrons...)
    genie_n["photons"] = 0;   // Photons are filled if there are electrons
    genie_n["pi_zeros"] = 0;  // Pions
    genie_n["piplus"] = 0;
    genie_n["piplus_range"] = 0;  // Pi+ that passes the kinematic cuts
    genie_n["piminus"] = 0;
    genie_n["pions"] = 0;
    genie_n["kaons"] = 0;  // Other mesons
    genie_n["charms"] = 0;
    genie_n["mesons"] = 0;
    genie_n["protons"] = 0;  // Baryons
    genie_n["neutrons"] = 0;
    genie_n["nucleons"] = 0;
    genie_n["heavy_baryons"] = 0;
    genie_n["nuclei"] = 0;
    genie_n["others"] = 0;  // others

    double tpi = 0;
    for (uint p = 0; p < FS_PDG.size(); p++) {
      // Get tpi
      tpi = (FS_energy[p] - 139.569);

      // Look at every particle's pdg
      switch (FS_PDG[p]) {
        case 13:
          genie_n["muons"]++;
          break;
        case 22:
          // Check the energy.  If below 10 MeV, it's a nuclear deexcitation,
          // which we ignore
          if (FS_energy[p] > 10) genie_n["photons"]++;
          break;
        case 11:
        case -11:
          genie_n["photons"]++;
          break;
        case 211:
          if (0. < tpi &&
              tpi < 350.)
            genie_n["piplus_range"]++;
          genie_n["piplus"]++;
          genie_n["pions"]++;
          genie_n["mesons"]++;
          break;
        case -211:
          genie_n["piminus"]++;
          genie_n["pions"]++;
          genie_n["mesons"]++;
          break;
        case 111:
          genie_n["pi_zeros"]++;
          // genie_n["photons"]++;
          genie_n["pions"]++;
          genie_n["mesons"]++;
          break;
        case 130:
        case 310:
        case 311:
        case 321:
        case -130:
        case -310:
        case -311:
        case -321:
          genie_n["kaons"]++;
          genie_n["mesons"]++;
          break;
        case 411:
        case 421:
        case 431:
        case -411:
        case -421:
        case -431:
          genie_n["charms"]++;
          genie_n["mesons"]++;
          break;
        case 2212:
          genie_n["protons"]++;
          genie_n["nucleons"]++;
          break;
        case 2112:
          genie_n["neutrons"]++;
          genie_n["nucleons"]++;
          break;
        case 2000000101:  // bindino is binding energy placeholder
          break;
        default:
          if (FS_PDG[p] > 3000 && FS_PDG[p] < 5000)
            genie_n["heavy_baryons"]++;
          else if (FS_PDG[p] > 1000000000 && FS_PDG[p] < 1099999999)
            genie_n["nuclei"]++;
          else
            genie_n["others"]++;
      }
    }

    return genie_n;
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
      if (GetPiCharge(chw, entry, idx) > 0 && t_pi > 0. && t_pi < 0.350
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

  bool IsInHexagon(double x, double y, double apothem) const {
    double lenOfSide = apothem * (2 / sqrt(3));
    double slope = (lenOfSide / 2.0) / apothem;
    double xp = fabs(x);
    double yp = fabs(y);

    if ((xp * xp + yp * yp) < apothem * apothem)
      return true;
    else if (xp <= apothem && yp * yp < lenOfSide / 2.0)
      return true;
    else if (xp <= apothem && yp < lenOfSide - xp * slope)
      return true;

    return false;
  }

  bool XYVertexSig(ChainWrapper& chw, int entry) {
    const double a = 850.0;
    const double x = (double)chw.GetValue("mc_vtx", entry, 0),
                 y = (double)chw.GetValue("mc_vtx", entry, 1);
    return IsInHexagon (x, y, a);
  }

  bool VtxSignal(ChainWrapper& chw, int entry) {
    bool Pass = true;
    Pass = Pass && zVertexSig(chw, entry);
    Pass = Pass && XYVertexSig(chw, entry);
    return Pass;
  }
  std::vector<int> FSPDGvector(ChainWrapper& chw, int entry){
    int nFSPart = chw.GetValue("mc_nFSPart", entry);
    std::vector<int> FSPDG;//(nFSPart); 
    for (int i = 0; i < nFSPart; ++i){
//      FSPDG[i] = chw.GetValue("mc_FSPartPDG", entry, i);
      FSPDG.push_back(chw.GetValue("mc_FSPartPDG", entry, i));
    }
    return FSPDG;
  }
  std::vector<double> FSEvector(ChainWrapper& chw, int entry){
    int nFSPart = chw.GetValue("mc_nFSPart", entry);
    std::vector<double> FSE;//(nFSPart); 
    for (int i = 0; i < nFSPart; ++i){
//      FSE[i] = chw.GetValue("mc_FSPartE", entry, i);
      FSE.push_back(chw.GetValue("mc_FSPartE", entry, i));
    }
    return FSE;
  }
  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry) {
    int pionIdx = GetHighestpionEnergyIdx(chw, entry);
    std::map<std::string, int> particles = GetParticleTopology(
        FSPDGvector(chw, entry),
        FSEvector(chw, entry));
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
//    if (Npions != 1) return false;
//    if (NOtherParticles(chw, entry) > 0) return false;
    if (particles.at("piplus_range") != 1) return false;
    if (!Is1PiPlus(particles)) return false;
    if (pmu < 1.5) return false;
    if (pmu > 20.) return false;
//  if (chw.GetValue("truth_N_pi0", entry) != 0) return false;
//  if (chw.GetValue("truth_N_pim", entry) != 0) return false;
//    std::cout << "Entry = " << entry << "\n";
    return true;
  }
};
//Main
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
//  double var_edges[] = {0., 1., 2., 3., 4., 5.5, 7.5, 10., 13., 20., 30.}; // pmu binning
  double var_edges[] = {0.,  0.02, 0.045, 0.06, 0.075, 0.1, 0.125, 0.166, 0.2, 0.35}; // tpi binning
  int var_nbins = (sizeof(var_edges) / sizeof(*var_edges)) - 1;

  std::cout << "Number of bins = " <<var_nbins << "\n";

  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dvar = new MinModDepCCQEXSec("tpi");
  ds_dvar->setBinEdges(var_nbins, var_edges);
  ds_dvar->setDimension(1);
  ds_dvar->setFluxIntLimits(0.0, 100.0);
  ds_dvar->setIsFluxIntegrated(true);
  ds_dvar->setNormalizationType(XSec::kPerNucleon);
  ds_dvar->setUniverses(0);  // default value, put 0 if you do not want
                             // universes to be included.
  ds_dvar->setVariable(XSec::kTPiPlus);// For tpi
//  ds_dvar->setVariable(XSec::kPLep); // For pmu

  loop.addXSec(ds_dvar);

  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =
      "GENIEXSECEXTRACT_"+
      playlistFile.substr(playlistFile.rfind("/") + 1, playlistFile.find(".")) +
      "_tpi.root";
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
  std::cout << "Ratio integrated Fluxes = " << Integrated_flux->GetBinContent(1) <<"   " << integrated_flux_other->GetBinContent(1) << "\n";

  Integrated_flux->Write();

  PlotUtils::MnvH1D* unfolded =
      (PlotUtils::MnvH1D*)loop.getXSecs()[0]->getXSecHist()->Clone("unfolded");
  unfolded->ClearAllErrorBands();
  unfolded->Reset();
  unfolded->Multiply(loop.getXSecs()[0]->getXSecHist(), Integrated_flux);

  unfolded->Write();
  PlotUtils::MnvH1D* selected = ds_dvar->getXSecHist()->Clone("selected");
  selected->Write();

  return 0;
}

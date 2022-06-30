#include "GENIEXSecExtract/XSecLooper.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include "PlotUtils/FluxReweighter.h"
#include <cstdlib>
typedef unsigned int uint;

class MinModDepCCQEXSec : public XSec
{
public:
  MinModDepCCQEXSec(const char* name)
    :XSec(name)
  {
  };

bool isCCInclusiveSignal( ChainWrapper& chw, int entry )
{
  double theta              = 0.;
  double true_muon_px   = (double)chw.GetValue("mc_primFSLepton",entry,0)/1000;
  double true_muon_py   = (double)chw.GetValue("mc_primFSLepton",entry,1)/1000;
  double true_muon_pz   = (double)chw.GetValue("mc_primFSLepton",entry,2)/1000;
  double numi_beam_angle_rad = -0.05887;
  double pyprime = -1.0*sin(numi_beam_angle_rad)*true_muon_pz + cos(numi_beam_angle_rad)*true_muon_py;
  double pzprime =  1.0*cos(numi_beam_angle_rad)*true_muon_pz + sin(numi_beam_angle_rad)*true_muon_py;
  double pSquare = pow(true_muon_px,2) + pow(pyprime,2) + pow(pzprime,2);
  theta = acos( pzprime / sqrt(pSquare) );
  theta *= 180./3.14159;

  if(!chw.GetValue("truth_is_fiducial",entry)) return false; //Doesn't work for MasterAnaDev tuples.  What does this even mean in the targets anyway? :(
  if( pzprime >= 1.5 && theta <= 20.0 ) return true;
  return false;

}

double CalcWexp(double Q2, double Ehad){
  double W = pow(938.3, 2.0) - Q2 +
             2.0 * (938.3)*Ehad;
  W = W > 0 ? sqrt(W) : 0.0;
  return W;
}

double GetEnuTrue(ChainWrapper& chw, int entry) {
  return (double)chw.GetValue("mc_incomingE",entry);
}

double GetQ2True(ChainWrapper& chw, int entry) { /* MeV^2 */
  return (double)chw.GetValue("mc_Q2",entry);
}

bool OnePionEvt(ChainWrapper& chw, int entry) {
  bool pass = true;
  pass = pass && (int)chw.GetValue("truth_N_pi0",entry) == 0;
  pass = pass && (int)chw.GetValue("truth_N_pim",entry) == 0;
  pass = pass && (int)chw.GetValue("truth_N_pip",entry) == 1;
  return pass;
} 

int NOtherParticles(ChainWrapper& chw, int entry){
  int n_other_particles = 0;
  n_other_particles = (int)chw.GetValue("truth_N_chargedK",entry) +
                      (int)chw.GetValue("truth_N_K0",entry) +
                      (int)chw.GetValue("truth_N_sigma",entry) +
                      (int)chw.GetValue("truth_N_lambda",entry);
  return n_other_particles;
}

double GetElepTrue(ChainWrapper& chw, int entry){
  return (double)chw.GetValue("mc_primFSLepton",entry,3);
}

double GetEhadTrue(ChainWrapper& chw, int entry) { 
  return GetEnuTrue(chw, entry) - GetElepTrue(chw, entry); 
}

double GetWexpTrue(ChainWrapper& chw, int entry){return CalcWexp(GetQ2True(chw, entry), GetEhadTrue(chw, entry));}

bool BenCuts(ChainWrapper& chw, int entry){
  double wexp = GetWexpTrue(chw, entry);
  bool pass = true;
  pass = pass && wexp > 0.;
  pass = pass && wexp < 1400;
  pass = pass && OnePionEvt(chw, entry);
  pass = pass && NOtherParticles(chw, entry) == 0;
  return pass;
}

  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_incoming", entry)!=14) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if(!isCCInclusiveSignal  ( chw, entry ) ) return false;
    if(!BenCuts  ( chw, entry ) ) return false;
    return true;
  }
};

int runXSecLooperMnv101(){
  TH1::AddDirectory(kFALSE); //Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  const std::string playlistFile = "/minerva/app/users/granados/cmtuser/MINERvA101/MINERvA-101-Cross-Section/MCME1A.txt";
  //const std::string playlistFile = "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/MCME1A_short.txt"; // shorter playlist for testing

  // Create the XSecLooper and tell it the input files
  // Inputs should be the merged ntuples:
  XSecLooper loop(playlistFile.c_str());

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);

  // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
  loop.setNumUniv(0); 
  loop.setFiducial(5990, 8340);
  loop.setPlaylist(PlotUtils::FluxReweighter::minervame1A);
  // Add the differential cross section dsigma/ds_dpT
  double pmu_edges[] = { 0., 1., 2., 3., 4., 5.5, 7.5, 10., 13., 20., 30. };
  int pmu_nbins = 10; 
 
  std::cout << "pmu \n";
 
  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpmu = new MinModDepCCQEXSec("pmu");
  ds_dpmu->setBinEdges(pmu_nbins, pmu_edges);
  ds_dpmu->setVariable(XSec::kPLep);
  ds_dpmu->setIsFluxIntegrated(true);
  ds_dpmu->setDimension(1);
  ds_dpmu->setFluxIntLimits(0.0, 100.0);
  ds_dpmu->setNormalizationType(XSec::kPerNucleon);  
  ds_dpmu->setUniverses(100); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpmu);

  loop.runLoop();
  std::cout << "Flag 1 \n";
  // Get the output histograms and save them to file
  string geniefilename =  "GENIEXSECEXTRACT_" + playlistFile.substr(playlistFile.rfind("/")+1, playlistFile.find(".")) + ".root";
  std::cout << "Flag 2\n";
  TFile fout(geniefilename.c_str(), "RECREATE");
  std::cout << "Flag 3\n";
  std::cout << "getXSecs vector size " << loop.getXSecs().size() << "\n";

  PlotUtils::MnvH1D* flux = (PlotUtils::MnvH1D*)loop.getFluxHist();

  for(uint i=0; i<loop.getXSecs().size(); ++i)
  {
    std::cout << "Flag loop first " << i <<"\n"; 
    loop.getXSecs()[i]->getXSecHist()->Write();
    std::cout << "Flag loop second " << i <<"\n";
    loop.getXSecs()[i]->getEvRateHist()->Write();
    std::cout << "Flag loop thrid " << i <<"\n";
  }

  PlotUtils::FluxReweighter* fluxReweighter = new PlotUtils::FluxReweighter( 14, true, PlotUtils::FluxReweighter::minervame1A, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6 );

  PlotUtils::MnvH1D* Integrated_flux = fluxReweighter->GetIntegratedFluxReweighted_FromInputFlux(flux, loop.getXSecs()[0]->getXSecHist(), 0.0, 100.0);

  PlotUtils::MnvH1D* integrated_flux_other = fluxReweighter->GetIntegratedFluxReweighted(14, loop.getXSecs()[0]->getXSecHist(), 0., 100.);
 
  Integrated_flux->Write();

  PlotUtils::MnvH1D* unfolded = (PlotUtils::MnvH1D*)loop.getXSecs()[0]->getXSecHist()->Clone("unfolded");
  unfolded->ClearAllErrorBands();
  unfolded->Reset();
  unfolded->Multiply(loop.getXSecs()[0]->getXSecHist(), Integrated_flux);

  unfolded->Write();
  PlotUtils::MnvH1D* selected = ds_dpmu->getXSecHist()->Clone("selected");
  selected->Write();
  std::cout << "Flag 4\n";
  
  return 0;
}

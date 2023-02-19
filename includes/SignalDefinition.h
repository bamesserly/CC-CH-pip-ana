#ifndef SignalDefinition_H
#define SignalDefinition_H

#include "includes/CVUniverse.h"
#include "includes/Constants.h"

enum SignalDefinition{ 
  kOnePi, kOnePiNoW, kNPi, kNPiNoW, kNSignalDefTypes 
};

double GetWCutValue(SignalDefinition signal_definition) {
  switch (signal_definition) {
    case kOnePi:
      return 1400.;
    case kNPi:
      return 2000.;
    case kOnePiNoW:
    case kNPiNoW:
      return 120000.;
    default:
      std::cout << "ERROR GetWCutValue" << std::endl;
      return -1.;
  }
}

void fillTopology( std::vector<int>& FS_PDG, std::vector<double>& FS_E, std::map< string, int >& genie_n )//Obtained from Aaron's code
{
  //Overarching categories: nucleons, mesons
  genie_n["muons"]         = 0;//Muons, photons (do we want electrons...)
  genie_n["photons"]       = 0;//Photons are filled if there are electrons
  genie_n["pi_zeros"]      = 0;//Pions
  genie_n["piplus"]        = 0;
  genie_n["piplus_range"]  = 0;//Pi+ that passes the kinematic cuts
  genie_n["piminus"]       = 0;
  genie_n["pions"]         = 0;
  genie_n["kaons"]         = 0;//Other mesons
  genie_n["charms"]        = 0;
  genie_n["mesons"]        = 0;
  genie_n["protons"]       = 0;//Baryons
  genie_n["neutrons"]      = 0;
  genie_n["nucleons"]      = 0;
  genie_n["heavy_baryons"] = 0;
  genie_n["nuclei"]        = 0;
  genie_n["others"]        = 0;//others

  double tpi = 0;
  for(uint p = 0; p < FS_PDG.size(); p++) {
    //Get tpi    
    tpi = (FS_E[p]-MinervaUnits::M_pion);

    //Look at every particle's pdg
    switch(FS_PDG[p]){
      case 13:
        genie_n["muons"]++; 
        break;
      case 22:
        //Check the energy.  If below 10 MeV, it's a nuclear deexcitation, which we ignore
        if( FS_E[p] > 10 ) genie_n["photons"]++; 
        break;
      case 11:
      case -11:
        genie_n["photons"]++; 
        break;
      case 211:
        if( 35 < tpi && tpi < 350 ) genie_n["piplus_range"]++;
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
        //genie_n["photons"]++;
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
      case 2000000101: //bindino is binding energy placeholder
        break;
      default:
        if( FS_PDG[p] > 3000 && FS_PDG[p] < 5000 ) genie_n["heavy_baryons"]++;
        else if( FS_PDG[p] > 1000000000 && FS_PDG[p] < 1099999999 ) genie_n["nuclei"]++;
        else genie_n["others"]++;
    }
  }
}

bool is1PiPlus( map< string, int > &particles )
{
//  if( incoming =! 14 || current =! 1 ) return false;

  //which particles are part of your signal (those you count, those you want any)
  //Don't need two variables, just here to keep bookkeeping straight
  TString counted_signal("muons piplus piplus_range pions mesons");
  TString any_baryon("protons neutrons nucleons heavy_baryons nuclei");

  //Count bkgd particles
  int nBkgd=0;
  for(std::map<string,int>::iterator it=particles.begin(); it!=particles.end(); ++it){
    string particle = it->first;
    int num         = it->second;
    if( counted_signal.Contains(particle) || any_baryon.Contains(particle)) continue; 
    nBkgd += num;
  }
                                                    //SIGNAL
  if( particles["muons"]  == 1                   && //1 Muon
      particles["piplus"] == 1                   && //1 Piplus
      particles["pions"]  == particles["piplus"] && //piplus is the only pion/meson
      particles["mesons"] == particles["piplus"] && 
      nBkgd               == 0                      ) return true;   //any number of baryons, nothing else

  return false;
}


// Number of abs(pdg) == 211 true TG4Trajectories which also:
// (1) are pip, (2) satisfy a KE restriction
int NSignalPions(const CVUniverse& univ) {
  int n_signal_pions = 0;
  int n_true_pions = univ.GetNChargedPionsTrue();
  for (TruePionIdx idx = 0; idx<n_true_pions; ++idx) {
    double t_pi = univ.GetTpiTrue(idx);
    double theta_pi = univ.GetThetapiTrue(idx);
    if(univ.GetPiChargeTrue(idx) > 0 
       && t_pi > 35
       && t_pi < 350
       //&& (theta_pi < 1.39626 || 1.74533 < theta_pi))
    )
      ++n_signal_pions;
  }
  return n_signal_pions;
}

int NOtherParticles(const CVUniverse& univ){
  int n_other_particles = 0;
  n_other_particles = univ.GetInt("truth_N_chargedK") +
                      univ.GetInt("truth_N_K0") +
                      univ.GetInt("truth_N_sigma") +
                      univ.GetInt("truth_N_lambda");
  return n_other_particles;
}

bool zVertexSig (const CVUniverse& univ){
        double vtxZ = univ.GetVecElem("mc_vtx",2);
        if (vtxZ > CCNuPionIncConsts::kZVtxMinCutVal && vtxZ < CCNuPionIncConsts::kZVtxMaxCutVal) return true;
        else return false;
}
bool XYVertexSig (const CVUniverse& univ){
  const double a = 850.0;
  const double x = univ.GetVecElem("mc_vtx",0), y = univ.GetVecElem("mc_vtx",1);
  return univ.IsInHexagon( x, y, a);
}

bool VtxSignal(const CVUniverse& univ){
  bool Pass = true;
  Pass = Pass && zVertexSig(univ);
  Pass = Pass && XYVertexSig(univ);
  return Pass;
}

bool PmuSignal(const CVUniverse& univ){
    if( univ.GetPmuTrue()/1000 < 1.5 || 20 < univ.GetPmuTrue()/1000 ) return false;
    else return true;
}

bool IsSignal(const CVUniverse& universe, SignalDefinition signal_definition = kOnePi) {
  int n_signal_pions = NSignalPions(universe);
  std::vector<int> mc_FSPartPDG = universe.GetVec<int>("mc_FSPartPDG");
  std::vector<double> mc_FSPartE = universe.GetVec<double>("mc_FSPartE");
  std::map< string, int > particles;
  fillTopology( mc_FSPartPDG, mc_FSPartE, particles );
  if( universe.GetInt("mc_current")  == 1 
       && universe.GetBool("truth_is_fiducial") 
       && VtxSignal(universe)
       && universe.GetInt("mc_incoming") == 14 
       && universe.GetThetalepTrue() < 0.226892803// 0.3491 // 20 deg
       && universe.GetWexpTrue() > 0
       && universe.GetWexpTrue() < GetWCutValue(signal_definition)
//       && n_signal_pions > 0
//       && NOtherParticles(universe) == 0
       && particles["piplus_range"] == 1
       && PmuSignal(universe)
       && is1PiPlus(particles)
       // && TODO Muon or neutrino energy cut
       //&& 1500. < universe.GetDouble("mc_incomingE") && universe.GetDouble("mc_incomingE") < 10000.
  ) {}
  
  else {
    return false;
  } 


  switch (signal_definition) {
    case kOnePi:/*
      if (n_signal_pions == 1
          && universe.GetInt("truth_N_pi0") == 0
          && universe.GetInt("truth_N_pim") == 0
         )
          return true;
      else return false;*/
    case kOnePiNoW:
      if (n_signal_pions == 1
          && universe.GetInt("truth_N_pi0") == 0                        
          && universe.GetInt("truth_N_pim") == 0
         )
          return true;
      else return false;
    case kNPi:
    case kNPiNoW:
      return true;

    default:
      std::cout << "IsSignal Error Unknown Signal Definition!" << std::endl;
      return false;
  }
}

std::string GetSignalName(SignalDefinition signal_definition){
  switch (signal_definition){
    case kOnePi:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X  (W < 1.4 GeV)";
    case kOnePiNoW:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X";
    case kNPi:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X  (W < 1.8 GeV)";
    case kNPiNoW:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X";
    default:
      return "UNKNOWN SIGNAL";
  }
}

std::string GetSignalFileTag(SignalDefinition signal_definition){
  switch (signal_definition){
    case kOnePi:
      return "1Pi";
    case kOnePiNoW:
      return "1PiNoW";
    case kNPi:
      return "NPi";
    case kNPiNoW:
      return "NPiNoW";
    default:
      return "UNKNOWN SIGNAL";
  }
}

#endif // Signal_H

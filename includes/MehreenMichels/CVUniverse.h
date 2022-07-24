// =============================================================================
// Base class for an un-systematically shifted (i.e. CV) universe. Implement
// "Get" functions for all the quantities that you need for your analysis.
//
// This class inherits from PU::MinervaUniverse, which in turn inherits from
// PU::BaseUniverse. PU::BU defines the interface with anatuples.
// 
// Within the class, "WeightFunctions" and "MuonFunctions" are included to gain
// access to standardized weight and muon variable getters. See:
// https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MinervaUniverse_Structure_
// for a full list of standardized functions you can use. In general, if a
// standard version of a function is available, you should be using it.
// =============================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include <iostream>

#include "PlotUtils/MinervaUniverse.h"

class CVUniverse : public PlotUtils::MinervaUniverse {

  public:
  #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
  #include "PlotUtils/TruthFunctions.h" //Getq3True
  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
      : PlotUtils::MinervaUniverse(chw, nsigma) {}

  virtual ~CVUniverse() {}

  // ========================================================================
  // Quantities defined here as constants for the sake of below. Definition
  // matched to Dan's CCQENuInclusiveME variables from:
  // `/minerva/app/users/drut1186/cmtuser/Minerva_v22r1p1_OrigCCQENuInc/Ana/CCQENu/ana_common/include/CCQENuUtils.h`
  // ========================================================================
  static constexpr double M_n = 939.56536;
  static constexpr double M_p = 938.272013;
  static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

  static constexpr int PDG_n = 2112;
  static constexpr int PDG_p = 2212;

  // ========================================================================
  // Write a "Get" function for all quantities access by your analysis.
  // For composite quantities (e.g. Enu) use a calculator function.
  //
  // In order to properly calculate muon variables and systematics use the
  // various functions defined in MinervaUniverse.
  // E.g. GetPmu, GetEmu, etc.
  // ========================================================================

  // Quantities only needed for cuts
  // Although unlikely, in principle these quanties could be shifted by a
  // systematic. And when they are, they'll only be shifted correctly if we
  // write these accessor functions.
  double GetTrueQ3() const
  {
 	return Getq3True()/pow(10,3);
  } 
  //Muon kinematics
  double GetMuonPT() const //GeV/c
  {
    return GetPmu()/1000. * sin(GetThetamu());
  }

  double GetMuonPz() const //GeV/c
  {
    return GetPmu()/1000. * cos(GetThetamu());
  }

  double GetMuonPTTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * sin(GetThetalepTrue());
  }

  double GetMuonPzTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * cos(GetThetalepTrue());
  }

  double GetEmuGeV() const //GeV
  {
    return GetEmu()/1000.;
  }

  double GetElepTrueGeV() const //GeV
  {
    return GetElepTrue()/1000.;
  }

  int GetInteractionType() const {
    return GetInt("mc_intType");
  }

  int GetTargetNucleon() const {
    return GetInt("mc_targetNucleon");
  }
  
  double GetBjorkenXTrue() const {
    return GetDouble("mc_Bjorkenx");
  }

  double GetBjorkenYTrue() const {
    return GetDouble("mc_Bjorkeny");
  }

  virtual bool IsMinosMatchMuon() const {
    return GetInt("has_interaction_vertex") == 1;
  }
  
  ROOT::Math::XYZTVector GetVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }

  ROOT::Math::XYZTVector GetTrueVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("mc_vtx").data());
    return result;
  }

  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
  }
  
  //TODO: If there was a spline correcting Eavail, it might not really be Eavail.
  //      Our energy correction spline, one of at least 2 I know of, corrects q0
  //      so that we get the right neutrino energy in an inclusive sample.  So,
  //      this function could be correcting for neutron energy which Eavail should
  //      not do.
  virtual double GetEavail() const {
    return GetDouble("recoilE_SplineCorrected");
  }
  
  virtual double GetQ2Reco() const{
    return GetDouble("qsquared_recoil");
  }

  //GetRecoilE is designed to match the NSF validation suite
  virtual double GetRecoilE() const {
    return GetVecElem("recoil_summed_energy", 0);
  }
  
  virtual double Getq3() const{
    double eavail = GetEavail()/pow(10,3);
    double q2 = GetQ2Reco() / pow(10,6);
    double q3mec = sqrt(eavail*eavail + q2);
    return q3mec;
  }
   
  virtual int GetCurrent() const { return GetInt("mc_current"); }

  virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }

  virtual double GetMuonQP() const {
    return GetDouble((GetAnaToolName() + "_minos_trk_qp").c_str());
  }

  //Some functions to match CCQENuInclusive treatment of DIS weighting. Name matches same Dan area as before.
  virtual double GetTrueExperimentersQ2() const {
    double Enu = GetEnuTrue(); //MeV
    double Emu = GetElepTrue(); //MeV
    double thetaMu = GetThetalepTrue();
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double CalcTrueExperimentersQ2(double Enu, double Emu, double thetaMu) const{
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double GetTrueExperimentersW() const {
    double nuclMass = M_nucleon;
    int struckNucl = GetTargetNucleon();
    if (struckNucl == PDG_n){
      nuclMass=M_n;
    }
    else if (struckNucl == PDG_p){
      nuclMass=M_p;
    }
    double Enu = GetEnuTrue();
    double Emu = GetElepTrue();
    double thetaMu = GetThetalepTrue();
    double Q2 = CalcTrueExperimentersQ2(Enu, Emu, thetaMu);
    return TMath::Sqrt(pow(nuclMass,2) + 2.0*(Enu-Emu)*nuclMass - Q2);
  }

  // Adding Branches for Mehreen's Michels
  virtual int GetNMichels() const{
      return GetInt("FittedMichel_michel_fitPass_sz");
  }

  virtual std::vector<int> GetTrueFSPDGCodes() const {
     std::vector<int> pdgcodes;
     int pdgsize = GetInt("mc_nFSPart");
     for (int i = 0; i< pdgsize; i++)
      {
          int pdg = GetVecElem("mc_FSPartPDG", i);
          pdgcodes.push_back(pdg);
      }

      return pdgcodes;
  }       
  
  virtual int GetTrueNMichels() const {
     int ntruemichels = 0;
     int nmichels = GetNMichels();
     for (int i =0; i< nmichels; i++)
     {
       
      double datafrac =  GetVecElem("FittedMichel_michel_datafraction", i);
      int fitted = GetVecElem("FittedMichel_michel_fitPass", i);
      if (fitted == 0) continue;
      if (datafrac > 0.5) continue;
      ntruemichels++;

     }

     return ntruemichels;
  }

 virtual int GetTrueNPionsinEvent() const {
     int npion = 0;
     int pdgsize = GetInt("mc_nFSPart");
     for (int i = 0; i< pdgsize; i++)
      {
          int pdg = GetVecElem("mc_FSPartPDG", i);
          if (pdg == 211) npion++;
      }
      return npion;

 } 
  
  virtual int GetNTruePions() const{
      return GetInt("FittedMichel_all_piontrajectory_trackID_sz");
  }

  virtual int GetPionParentID(int i) const {
     return GetVecElem("FittedMichel_all_piontrajectory_ParentID", i);
  }

  virtual int GetPionPDG(int i) const{
     return GetVecElem("FittedMichel_all_piontrajectory_pdg", i);
  }

  virtual double GetPionE(int i) const{
     return GetVecElem("FittedMichel_all_piontrajectory_energy",i)/pow(10,3);
  }

  virtual double GetPionP(int i) const{
    return GetVecElem("FittedMichel_all_piontrajectory_momentum", i)/pow(10,3);
  }

  virtual double GetPionMass(int i) const{
    double pionmass = sqrt(pow(GetPionE(i), 2) - pow(GetPionP(i), 2));
    return pionmass;
  }

  virtual double GetPionKE(int i) const{
    double pionKE = GetPionE(i) - GetPionMass(i);
  }

  virtual double thetaXWRTBeam(double x, double y, double z) const{

  double pzp = cos( MinervaUnits::numi_beam_angle_rad)*z + sin( MinervaUnits::numi_beam_angle_rad)*y;
  int sign = (x<0.) ? -1 : 1;
  double denom2 = pow(x,2) + pow(pzp,2);
  if( 0. == denom2 ) return -999.;
  else
    return sign * acos(pzp / sqrt(denom2));
  }

  virtual double thetaWRTBeam(double x, double y, double z) const{
      double pyp = -1.0 *sin( MinervaUnits::numi_beam_angle_rad)*z + cos( MinervaUnits::numi_beam_angle_rad )*y;
      double pzp = cos( MinervaUnits::numi_beam_angle_rad )*z + sin( MinervaUnits::numi_beam_angle_rad )*y;
      double denom2 = pow(x,2) + pow(pyp,2) + pow(pzp,2);
      if( 0. == denom2 ) return -9999.;
      else return acos(pzp / sqrt(denom2) );
  } 

  virtual double GetTrueTpi() const {
     int nFSpi = GetNTruePions();
     double pionKE = 9999.;
     for (int i = 0; i < nFSpi; i++){
         int pdg = GetVecElem("FittedMichel_all_piontrajectory_pdg", i);
         int pitrackid = GetVecElem("FittedMichel_all_piontrajectory_ParentID", i);

         double energy = GetVecElem("FittedMichel_all_piontrajectory_energy",i);
         double p = GetVecElem("FittedMichel_all_piontrajectory_momentum", i);
         double mass = sqrt(pow(energy,2) - pow(p, 2));
         double tpi = energy - mass;
         if (tpi <= pionKE) pionKE = tpi;
      }
       
      return pionKE;
   }
 virtual void PrintArachneLink() const {
  int link_size = 200;
  char link[link_size];
  int run = GetInt("mc_run");
  int subrun = GetInt("mc_subrun");
  int gate = GetInt("mc_nthEvtInFile") + 1;
  int slice = GetVecElem("slice_numbers", 0);
  sprintf(link,
          "http://minerva05.fnal.gov/Arachne/"
          "arachne.html\?det=SIM_minerva&recoVer=v21r1p1&run=%d&subrun=%d&gate="
          "%d&slice=%d",
          run, subrun, gate, slice);
       // strncpy(); // FAIL
       // memcpy();  // FAIL
       std::cout << link << std::endl;
       }
 virtual double GetTrueEAvail() const
 {
     double Eavail = 0.0;
     int pdgsize = GetInt("mc_nFSPart");
     for (int i = 0; i< pdgsize; i++)
     {
        int pdg = GetVecElem("mc_FSPartPDG", i);  
        double energy = GetVecElem("mc_FSPartE", i); // hopefully this is in MeV
        if (pdg == 211) Eavail+= energy - 139.57; // subtracting pion mass to get Kinetic energy
        if (pdg == 2212) Eavail += energy - 938.28; // proton
        if (pdg == 111) Eavail += energy; // pi0
        if (pdg == 22) Eavail += energy; // photons
        if (pdg == 311) Eavail += energy - 497.611/2.0; // K0 ???? - Kaon rest mass / 2
        if (pdg == 321) Eavail += energy - 497.611/2.0; // Kaon+ Kinetic Energy  divide by Kmass/2 
        if (pdg == 551) Eavail += energy; //Adding etas
     }

     return Eavail; // in MeV
 }


  //Still needed for some systematics to compile, but shouldn't be used for reweighting anymore.
  protected:
  #include "PlotUtils/WeightFunctions.h" // Get*Weight
};

#endif

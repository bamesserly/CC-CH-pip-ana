#ifndef CohDiffractiveSystematics_CXX
#define CohDiffractiveSystematics_CXX 1

#include "CohDiffractiveSystematics.h"
#include "PlotUtils/TargetUtils.h"
#include <iostream>

std::vector<CVUniverse*> GetCohDiffractiveSystematics( PlotUtils::ChainWrapper* chain ) {
  std::vector<CVUniverse*> ret;
  ret.push_back(new DiffractiveUniverse( chain, -1., fracDiffractiveUnc ));
  ret.push_back(new DiffractiveUniverse( chain, 1. , fracDiffractiveUnc ));

  ret.push_back(new CoherentPiPlasticUniverse( chain, -1., fracCoherentPiUncTracker ));
  ret.push_back(new CoherentPiPlasticUniverse( chain, 1. , fracCoherentPiUncTracker ));

  ret.push_back(new CoherentPiCarbonUniverse( chain, -1., fracCoherentPiUncNuclear ));
  ret.push_back(new CoherentPiCarbonUniverse( chain, 1. , fracCoherentPiUncNuclear ));

  ret.push_back(new CoherentPiWaterUniverse( chain, -1., fracCoherentPiUncNuclear ));
  ret.push_back(new CoherentPiWaterUniverse( chain, 1. , fracCoherentPiUncNuclear ));

  ret.push_back(new CoherentPiIronUniverse( chain, -1., fracCoherentPiUncNuclear ));
  ret.push_back(new CoherentPiIronUniverse( chain, 1. , fracCoherentPiUncNuclear ));

  ret.push_back(new CoherentPiLeadUniverse( chain, -1., fracCoherentPiUncNuclear ));
  ret.push_back(new CoherentPiLeadUniverse( chain, 1. , fracCoherentPiUncNuclear ));

  return ret;
};

std::map< std::string, std::vector<CVUniverse*> > GetCohDiffractiveSystematicsMap(PlotUtils::ChainWrapper* chain) {
  std::map< std::string, std::vector<CVUniverse*> > ret;
  
  ret["DiffractiveModelUnc"].push_back(new DiffractiveUniverse( chain, -1., fracDiffractiveUnc ));
  ret["DiffractiveModelUnc"].push_back(new DiffractiveUniverse( chain, 1. , fracDiffractiveUnc ));

  ret["CoherentPiUncTracker_CH"].push_back(new CoherentPiPlasticUniverse( chain, -1., fracCoherentPiUncTracker ));
  ret["CoherentPiUncTracker_CH"].push_back(new CoherentPiPlasticUniverse( chain, 1. , fracCoherentPiUncTracker ));

  ret["CoherentPiUncNuclear_C"].push_back(new CoherentPiCarbonUniverse( chain, -1., fracCoherentPiUncNuclear ));
  ret["CoherentPiUncNuclear_C"].push_back(new CoherentPiCarbonUniverse( chain, 1. , fracCoherentPiUncNuclear ));

  ret["CoherentPiUncNuclear_H2O"].push_back(new CoherentPiWaterUniverse( chain, -1., fracCoherentPiUncNuclear ));
  ret["CoherentPiUncNuclear_H2O"].push_back(new CoherentPiWaterUniverse( chain, 1. , fracCoherentPiUncNuclear ));

  ret["CoherentPiUncNuclear_Fe"].push_back(new CoherentPiIronUniverse( chain, -1., fracCoherentPiUncNuclear ));
  ret["CoherentPiUncNuclear_Fe"].push_back(new CoherentPiIronUniverse( chain, 1. , fracCoherentPiUncNuclear ));

  ret["CoherentPiUncNuclear_Pb"].push_back(new CoherentPiLeadUniverse( chain, -1., fracCoherentPiUncNuclear ));
  ret["CoherentPiUncNuclear_Pb"].push_back(new CoherentPiLeadUniverse( chain, 1. , fracCoherentPiUncNuclear ));

  return ret;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Class Definitions
// Constructor
DiffractiveUniverse::DiffractiveUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracDiffUnc ) : 
      CVUniverse(chw, nsigma), m_fracDiffUnc(fracDiffUnc) {}

double DiffractiveUniverse::GetDiffractiveWeight() const {
  if( GetInt("mc_intType") != 4) return 1.;
  if( !IsInPlastic() ) return 1.; 

  double cv_wgt = CVUniverse::GetDiffractiveWeight();
  double shift_val = m_nsigma * m_fracDiffUnc*(cv_wgt-1);
  return cv_wgt + shift_val;
}

std::string DiffractiveUniverse::ShortName() const { return "DiffractiveModelUnc"; }

std::string DiffractiveUniverse::LatexName() const { return "Diffractive Pion Uncertainty"; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class Definitions
// Constructor
CoherentPiPlasticUniverse::CoherentPiPlasticUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc ) : 
      CVUniverse(chw, nsigma), m_fracCohPiUnc(fracCohPiUnc) {}

double CoherentPiPlasticUniverse::GetCoherentPiWeight( double thpi_true, double tpi_true ) const {
  if( GetInt("mc_intType") != 4 ) return 1.;

  //In Plastic?
  if( PlotUtils::TargetUtils::Get().InPassiveTargetMC( GetIntVtxXTrue(), 
                                                       GetIntVtxYTrue(), 
                                                       GetIntVtxZTrue(), 
                                                       GetInt("mc_targetZ") ) ) return 1.;

  double shift_val = 1 + m_nsigma * m_fracCohPiUnc;
  return shift_val*CVUniverse::GetCoherentPiWeight(thpi_true,tpi_true);
}

std::string CoherentPiPlasticUniverse::ShortName() const { return "CoherentPiUnc_CH"; }

std::string CoherentPiPlasticUniverse::LatexName() const { return "Coherent Pion Scintillator Uncertainty"; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class Definitions
// Constructor
CoherentPiCarbonUniverse::CoherentPiCarbonUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc ) : 
      CVUniverse(chw, nsigma), m_fracCohPiUnc(fracCohPiUnc) {}

double CoherentPiCarbonUniverse::GetCoherentPiWeight( double thpi_true, double tpi_true ) const {
  if (GetInt("mc_intType") != 4) return 1.;
  
  //In Carbon?
  if( !PlotUtils::TargetUtils::Get().InCarbonTargetVolMC( GetIntVtxXTrue(), 
                                                          GetIntVtxYTrue(), 
                                                          GetIntVtxZTrue() ) ||
                                                        ( GetInt("mc_targetZ") != 6 ) ) return 1.;

  double shift_val = 1 + m_nsigma * m_fracCohPiUnc;
  return shift_val*CVUniverse::GetCoherentPiWeight(thpi_true,tpi_true);
}

std::string CoherentPiCarbonUniverse::ShortName() const { return "CoherentPiUnc_C"; }

std::string CoherentPiCarbonUniverse::LatexName() const { return "Coherent Pion Carbon Uncertainty"; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class Definitions
// Constructor
CoherentPiWaterUniverse::CoherentPiWaterUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc ) : 
      CVUniverse(chw, nsigma), m_fracCohPiUnc(fracCohPiUnc) {}

double CoherentPiWaterUniverse::GetCoherentPiWeight( double thpi_true, double tpi_true ) const {
  if (GetInt("mc_intType") != 4) return 1.;

  //In Water?
  if( !PlotUtils::TargetUtils::Get().InWaterTargetMC( GetIntVtxXTrue(), 
                                                      GetIntVtxYTrue(), 
                                                      GetIntVtxZTrue(), 
                                                      GetInt("mc_targetZ") ) ) return 1.;

  double shift_val = 1 + m_nsigma * m_fracCohPiUnc;
  return shift_val*CVUniverse::GetCoherentPiWeight(thpi_true,tpi_true);
}

std::string CoherentPiWaterUniverse::ShortName() const { return "CoherentPiUnc_H2O"; }

std::string CoherentPiWaterUniverse::LatexName() const { return "Coherent Pion Water Uncertainty"; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class Definitions
// Constructor
CoherentPiIronUniverse::CoherentPiIronUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc ) : 
      CVUniverse(chw, nsigma), m_fracCohPiUnc(fracCohPiUnc) {}

double CoherentPiIronUniverse::GetCoherentPiWeight( double thpi_true, double tpi_true ) const {
  if (GetInt("mc_intType") != 4) return 1.;

  //In Iron?
  if( !PlotUtils::TargetUtils::Get().InIronTargetVolMC( GetIntVtxXTrue(), 
                                                        GetIntVtxYTrue(), 
                                                        GetIntVtxZTrue() ) ||
                                                      ( GetInt("mc_targetZ") != 26 ) ) return 1.;

  double shift_val = 1 + m_nsigma * m_fracCohPiUnc;
  return shift_val*CVUniverse::GetCoherentPiWeight(thpi_true,tpi_true);
}

std::string CoherentPiIronUniverse::ShortName() const { return "CoherentPiUnc_Fe"; }

std::string CoherentPiIronUniverse::LatexName() const { return "Coherent Pion Iron Uncertainty"; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class Definitions
// Constructor
CoherentPiLeadUniverse::CoherentPiLeadUniverse( PlotUtils::ChainWrapper* chw, double nsigma, double fracCohPiUnc ) : 
      CVUniverse(chw, nsigma), m_fracCohPiUnc(fracCohPiUnc) {}

double CoherentPiLeadUniverse::GetCoherentPiWeight( double thpi_true, double tpi_true ) const {
  if (GetInt("mc_intType") != 4) return 1.;

  //In Lead?
  if( !PlotUtils::TargetUtils::Get().InLeadTargetVolMC( GetIntVtxXTrue(), 
                                                        GetIntVtxYTrue(), 
                                                        GetIntVtxZTrue() ) ||
                                                      ( GetInt("mc_targetZ") != 82 ) ) return 1.;

  double shift_val = 1 + m_nsigma * m_fracCohPiUnc;
  return shift_val*CVUniverse::GetCoherentPiWeight(thpi_true,tpi_true);
}

std::string CoherentPiLeadUniverse::ShortName() const { return "CoherentPiUnc_Pb"; }

std::string CoherentPiLeadUniverse::LatexName() const { return "Coherent Pion Lead Uncertainty"; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif // CohDiffractiveSystematics_CXX

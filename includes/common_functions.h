#ifndef common_functions_h
#define common_functions_h

#include <algorithm> // erase, remove_if
#include "TFile.h"
#include "PlotUtils/MnvH1D.h"
#include "Constants.h" // CCNuPionIncConsts::PI
#ifndef __CINT__
#include "Variable.h"
#endif // __CINT__

class Variable;

//==============================================================================
// Misc Utility Functions
// * Write POT number to a root file
// * Manipulate vectors of variables
// * Angles/Unit conversions
//==============================================================================

// Write POT to file as a number
void WritePOT(TFile& fout, const bool is_mc, const float pot) {
  fout.Write(0,TObject::kOverwrite);
  fout.cd();
  const char* name = is_mc ? "mc_pot" : "data_pot";
  PlotUtils::MnvH1D* h_pot = new PlotUtils::MnvH1D(name, name, 1, 0., 1.);
  h_pot->Fill(0.5, pot);
  h_pot->Write();

  //TVectorD mc_pot(1);
  //mc_pot[0] = util.m_mc_pot;
  //mc_pot.Write("mc_pot");
}

// Does a vector of variables contain a certain variable?
bool HasVar(std::vector<Variable*> variables, std::string name) {
#ifndef __CINT__ // CINT doesn't like lambdas
  auto it = find_if (variables.begin(), variables.end(), 
                      [&name](Variable* v) {return v->Name() == name;});
  if (it != variables.end())
    return true;
  else
    return false;
#endif // __CINT__
}

// Get a certain variable from a vector of variables
Variable* GetVar(std::vector<Variable*> variables, std::string name) {
#ifndef __CINT__ // CINT doesn't like lambdas
  auto it = find_if (variables.begin(), variables.end(), 
                      [&name](Variable* v) {return v->Name() == name;});
  if (it != variables.end()) {
    return *it;
  }
  else {
    std::cerr << name << " variable not found!\n";
    return nullptr;
  }
#endif // __CINT__
}

// Angle and Geometry Functions
// Restrict to [0,pi]
double FixAngle(double angle) {
  double ret = angle;
  if (ret < 0.0) ret = -1.0*ret;
  if (ret > CCNuPionIncConsts::PI) ret = 2.0*CCNuPionIncConsts::PI - ret;
  return ret;
}

double ConvertRadToDeg(double rad) { return rad*180./CCNuPionIncConsts::PI; }

#endif // common_functions_h

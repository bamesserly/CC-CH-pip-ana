#ifndef common_functions_h
#define common_functions_h

#include <algorithm> // erase, remove_if
#include "PlotUtils/MnvH1D.h"

#include "TFile.h"
#include "TKey.h"

#include "Constants.h" // CCNuPionIncConsts::PI
#include "MacroUtil.h"
#ifndef __CINT__
#include "Variable.h"
#endif // __CINT__

class Variable;

//==============================================================================
// Misc Utility Functions
// * Write POT number to a root file
// * Make HW from universes and a variable's binning
// * Copy all hists from one file to another
// * Get POT from input file, write it to output file and MacroUtil
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

// Make a HistWrapper from a variable's binning
void InitializeHW(Variable* var, std::string name, std::string label,
                  UniverseMap error_bands, CVHW& hw) {
  MH1D* hist = new MH1D(name.c_str(), label.c_str(), var->NBins(),
                        var->m_hists.m_bins_array.GetArray());

  // make CVHW from MH1D
  const bool clear_bands = true;
  hw = CVHW(hist, error_bands, clear_bands);

  delete hist;
}

// Copy all hists of one file into another file
void CopyHists(TFile& fin, TFile& fout) {
  TIter nextkey(fin.GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)nextkey())) {
    TString name = key->GetName();
    PlotUtils::MnvH1D* h = (PlotUtils::MnvH1D*)key->ReadObj();
    if (!h) {
      std::cout << "Invalid histogram " << name << "\n";
      continue;
    }

    // h->Scale(potScale);
    fout.cd();
    // std::cout << "Writing histogram " << name << "\n";
    h->Write(name);
  }
}

// Get POT from input file branch, save it to our MacroUtil and to output file.
void SetPOT(TFile& fin, TFile& fout, CCPi::MacroUtil& util) {
  util.m_mc_pot = -1;
  if (util.m_data_pot == 0) std::cout << "WARNING DATA POT == 0\n";

  // MC
  // TODO make this safer
  PlotUtils::MnvH1D* h_mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float mc_pot = h_mc_pot->GetBinContent(1);
  util.m_mc_pot = mc_pot;
  util.m_pot_scale = util.m_data_pot / util.m_mc_pot;

  // Data
  bool is_mc = false;
  WritePOT(fout, is_mc, util.m_data_pot);
}

// Loop variables and save specifically the data hists to file
void SaveDataHistsToFile(TFile& fout, std::vector<Variable*> variables) {
  std::cout << "Saving Data Hists\n\n";
  // fout.cd();
  for (auto v : variables) {
    std::string name = v->Name();
    v->m_hists.m_selection_data->GetXaxis()->SetTitle(
        v->m_hists.m_selection_mc.hist->GetXaxis()->GetTitle());
    v->m_hists.m_selection_data->Write(Form("selection_data_%s", name.c_str()));
    if (name == sidebands::kFitVarString) {
      v->m_hists.m_wsidebandfit_data->Write(
          Form("wsidebandfit_data_%s", name.c_str()));
    }
  }
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

#endif // common_functions_h

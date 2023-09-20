#ifndef common_var_functions_h
#define common_var_functions_h

#include "Constants.h" // CCNuPionIncConsts::PI
#include "TFile.h"
#include "Variable.h"

// Make a HistWrapper from a variable's binning
void InitializeHW(VariableBase* var, std::string name, std::string label,
                  UniverseMap error_bands, CVHW& hw) {
  MH1D* hist = new MH1D(name.c_str(), label.c_str(), var->NBins(),
                        var->m_hists.m_bins_array.GetArray());

  // make CVHW from MH1D
  const bool clear_bands = true;
  hw = CVHW(hist, error_bands, clear_bands);

  delete hist;
}

// Loop variables and save specifically the data hists to file
void SaveDataHistsToFile(TFile& fout, std::vector<VariableBase*> variables) {
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
      v->m_hists.m_wsideband_data->Write(
          Form("wsideband_data_%s", name.c_str()));
    }
  }
}

// Does a vector of variables contain a certain variable?
bool HasVar(std::vector<VariableBase*> variables, std::string name) {
  auto it = find_if (variables.begin(), variables.end(), 
                      [&name](VariableBase* v) {return v->Name() == name;});
  if (it != variables.end())
    return true;
  else
    return false;
}

// Get a certain variable from a vector of variables
VariableBase* GetVar(std::vector<VariableBase*> variables, std::string name) {
  auto it = find_if (variables.begin(), variables.end(), 
                      [&name](VariableBase* v) {return v->Name() == name;});
  if (it != variables.end()) {
    return *it;
  }
  else {
    std::cerr << name << " variable not found!\n";
    return nullptr;
  }
}

#endif

#ifndef common_functions_h
#define common_functions_h

#include <algorithm> // erase, remove_if
#include "PlotUtils/MnvH1D.h"

#include "TFile.h"
#include "TKey.h"

#include "Constants.h" // CCNuPionIncConsts::PI
#include "MacroUtil.h"

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
  fout.Flush();

  //TVectorD mc_pot(1);
  //mc_pot[0] = util.m_mc_pot;
  //mc_pot.Write("mc_pot");
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


#endif // common_functions_h

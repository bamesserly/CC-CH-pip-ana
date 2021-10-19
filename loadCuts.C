#include "TString.h"
#include "TSystem.h"
void loadCuts()
{
  TString path(TString::Format("%s/CC-CH-pip-ana/includes/", gSystem->Getenv("TOPDIR")));
  //std::cout << "path " << path << "\n";
  TString oldpath=gSystem->GetIncludePath();
  oldpath+=" -I";
  oldpath+=path;
  gSystem->SetIncludePath(oldpath);
  gSystem->CompileMacro("CVUniverse.cxx",       "kF");
  gSystem->CompileMacro("Cuts.cxx",             "k");
  gSystem->CompileMacro("StackedHistogram.cxx", "k");
  gSystem->CompileMacro("Histograms.cxx",       "k");
  gSystem->CompileMacro("Variable.cxx",         "k");
  gSystem->CompileMacro("HadronVariable.cxx",   "k");
  gSystem->CompileMacro("MacroUtil.cxx",        "kF");
  gSystem->CompileMacro("CCPiEvent.cxx",        "kF");
  gSystem->CompileMacro("WSidebandFitter.cxx",  "k");
}

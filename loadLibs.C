#include <iostream>

#include "TInterpreter.h"
#include "TString.h"
#include "TSystem.h"

void loadIncludes(bool verbose_cvu) {
  const char* cvu_flags = verbose_cvu ? "kf" : "kfg";
  TString path(
      TString::Format("%s/cc-ch-pip-ana/includes/", gSystem->Getenv("TOPDIR")));
  // std::cout << "path " << path << "\n";
  TString oldpath = gSystem->GetIncludePath();
  oldpath += " -I";
  oldpath += path;
  gSystem->SetIncludePath(oldpath);
  gSystem->CompileMacro("CVUniverse.cxx", cvu_flags);
  gSystem->CompileMacro("Cuts.cxx", "k");
  gSystem->CompileMacro("StackedHistogram.cxx", "k");
  gSystem->CompileMacro("Histograms.cxx", "k");
  gSystem->CompileMacro("Variable.cxx", "k");
  gSystem->CompileMacro("HadronVariable.cxx", "k");
  gSystem->CompileMacro("MacroUtil.cxx", "k");
  gSystem->CompileMacro("CCPiEvent.cxx", "k");
  gSystem->CompileMacro("WSidebandFitter.cxx", "k");
  gSystem->CompileMacro("CohDiffractiveSystematics.cxx", "k");
}

void loadLibs(bool verbose_cvu = true) {
  // MnvH1D hides approximately everything, so just turn off the pages
  // of compiler warnings. It would have been easier to do this by
  // using SetFlagsDebug(), but those flags get put before the default
  // settings in the compile line, and so the default settings win
  TString makeSharedLib(gSystem->GetMakeSharedLib());
  makeSharedLib.ReplaceAll("-Woverloaded-virtual", "-Wno-overloaded-virtual");
  gSystem->SetMakeSharedLib(makeSharedLib);

  // compile includes
  loadIncludes(verbose_cvu);

  // gSystem->CompileMacro("util/plot/GridCanvas.cxx", "k");
  // Long complicated reason to do this because of using TExec to set colour
  // palettes
  gSystem->CompileMacro("includes/myPlotStyle.h", "k");
}

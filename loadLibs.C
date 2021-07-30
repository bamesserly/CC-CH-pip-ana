#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"

#include <iostream>

void loadLibs()
{
#ifdef __CINT__
  cout << "Run loadLibs.C+" << endl;
  exit(1);
#else
  gSystem->SetAclicMode(TSystem::kDebug);

  // MnvH1D hides approximately everything, so just turn off the pages
  // of compiler warnings. It would have been easier to do this by
  // using SetFlagsDebug(), but those flags get put before the default
  // settings in the compile line, and so the default settings win
  TString makeSharedLib(gSystem->GetMakeSharedLib());
  makeSharedLib.ReplaceAll("-Woverloaded-virtual", "-Wno-overloaded-virtual");
  gSystem->SetMakeSharedLib(makeSharedLib);

  // Compile local scripts
  gInterpreter->ExecuteMacro("loadCuts.C");

  //gSystem->CompileMacro("util/plot/GridCanvas.cxx", "k");
  // Long complicated reason to do this because of using TExec to set colour palettes
  gSystem->CompileMacro("includes/myPlotStyle.h", "k");
#endif
}



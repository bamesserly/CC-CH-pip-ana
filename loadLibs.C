#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "Cintex/Cintex.h"

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

  // Add to include dirs
  const int kNIncDirs=2;
  TString incDirs[kNIncDirs] = {
    TString::Format("%s", gSystem->Getenv("PLOTUTILSROOT")),
    TString::Format("%s", gSystem->Getenv("UNFOLDUTILSROOT")),
  };

  for(int i=0; i<kNIncDirs; ++i){
    //cout << "Adding " << incDirs[i] << endl;
    gInterpreter->AddIncludePath(incDirs[i]);
  }

  // Link plotutils, mcreweight, and others
  const int kNLibs=3;
  // TODO: Will probably need to add T2KReweight and the neut reweight library to this
  const char* libs[kNLibs] = {
    "libplotutils.so",
    "libCintex.so",
    "libUnfoldUtils.so"
  };

  for(int i=0; i<kNLibs; ++i){
    cout << "Loading " << libs[i] << endl;
    gSystem->Load(libs[i]);
  }

  ROOT::Cintex::Cintex::Enable();

  gSystem->CompileMacro("util/plot/GridCanvas.cxx", "k");
  // Long complicated reason to do this because of using TExec to set colour palettes
  gSystem->CompileMacro("util/plot/myPlotStyle.h", "k");
#endif
}



#define MC_OUT_FILE_NAME \
  "runEventLoopMC_Dec312021_q3bin2_atleast1pi_distcuts.root"
#define DATA_OUT_FILE_NAME \
  "runEventLoopData_Dec312021_q3bin2_atleast1pi_distcuts.root"

#define USAGE                                                                  \
  "\n*** USAGE ***\n"                                                          \
  "runEventLoop <dataPlaylist.txt> <mcPlaylist.txt>\n\n"                       \
  "*** Explanation ***\n"                                                      \
  "Reduce MasterAnaDev AnaTuples to event selection histograms to extract a\n" \
  "single-differential inclusive cross section for the 2021 MINERvA 101 "      \
  "tutorial.\n\n"                                                              \
  "*** The Input Files ***\n"                                                  \
  "Playlist files are plaintext files with 1 file name per line.  Filenames "  \
  "may be\n"                                                                   \
  "xrootd URLs or refer to the local filesystem.  The first playlist file's\n" \
  "entries will be treated like data, and the second playlist's entries "      \
  "must\n"                                                                     \
  "have the \"Truth\" tree to use for calculating the efficiency "             \
  "denominator.\n\n"                                                           \
  "*** Output ***\n"                                                           \
  "Produces a two files, " MC_OUT_FILE_NAME " and " DATA_OUT_FILE_NAME         \
  ", with\n"                                                                   \
  "all histograms needed for the ExtractCrossSection program also built by "   \
  "this\n"                                                                     \
  "package.  You'll need a .rootlogon.C that loads ROOT object definitions "   \
  "from\n"                                                                     \
  "PlotUtils to access systematics information from these files.\n\n"          \
  "*** Environment Variables ***\n"                                            \
  "Setting up this package appends to PATH and LD_LIBRARY_PATH.  "             \
  "PLOTUTILSROOT,\n"                                                           \
  "MPARAMFILESROOT, and MPARAMFILES must be set according to the setup "       \
  "scripts in\n"                                                               \
  "those packages for systematics and flux reweighters to function.\n"         \
  "If MNV101_SKIP_SYST is defined at all, output histograms will have no "     \
  "error bands.\n"                                                             \
  "This is useful for debugging the CV and running warping studies.\n\n"       \
  "*** Return Codes ***\n"                                                     \
  "0 indicates success.  All histograms are valid only in this case.  Any "    \
  "other\n"                                                                    \
  "return code indicates that histograms should not be used.  Error "          \
  "messages\n"                                                                 \
  "about what went wrong will be printed to stderr.  So, they'll end up in "   \
  "your\n"                                                                     \
  "terminal, but you can separate them from everything else with something "   \
  "like:\n"                                                                    \
  "\"runEventLoop data.txt mc.txt 2> errors.txt\"\n"

enum ErrorCodes {
  success = 0,
  badCmdLine = 1,
  badInputFile = 2,
  badFileRead = 3,
  badOutputFile = 4
};

// PlotUtils includes
// No junk from PlotUtils please!  I already
// know that MnvH1D does horrible horrible things.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

// Includes from this package
#include "cuts/BestMichelDistance2D.h"
#include "cuts/MaxPzMu.h"
#include "cuts/SignalDefinition.h"
#include "cuts/VtxMatchFirst.h"
#include "cuts/hasMichel.h"
#include "cuts/hasTruePion.h"
#include "cuts/q3RecoCut.h"
#include "event/CVUniverse.h"
#include "event/Michel.h"
#include "event/MichelEvent.h"
#include "studies/PerMichel2DVar.h"
#include "studies/PerMichelEventVarByGENIELabel.h"
#include "studies/PerMichelVarByGENIELabel.h"
#include "studies/Study.h"
#include "systematics/Systematics.h"
#include "util/GetFluxIntegral.h"
#include "util/GetPlaylist.h"
#include "util/Variable.h"
#include "util/Variable2D.h"
//#include "Binning.h" //TODO: Fix me

// PlotUtils includes
#include "PlotUtils/CCInclusiveCuts.h"
#include "PlotUtils/CCInclusiveSignal.h"
#include "PlotUtils/CrashOnROOTMessage.h"  //Sets up ROOT's debug callbacks by itself
#include "PlotUtils/Cutter.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/makeChainWrapper.h"
#pragma GCC diagnostic pop

// ROOT includes
#include "TParameter.h"

// c++ includes
#include <cstdlib>  //getenv()
#include <iostream>

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*>> error_bands,
    std::vector<Variable*> vars, std::vector<Variable2D*> vars2D,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model) {
  assert(!error_bands["cv"].empty() &&
         "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  const int nEntries = chain->GetEntries();
  for (int i = 0; i < nEntries; ++i) {
    if (i % 1000 == 0)
      std::cout << i << " / " << nEntries << "\r" << std::flush;

    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands) {
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes) {
        MichelEvent
            myevent;  // make sure your event is inside the error band loop.

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);

        // This is where you would Access/create a Michel

        // weight is ignored in isMCSelected() for all but the CV Universe.
        if (!michelcuts.isMCSelected(*universe, myevent, cvWeight).all())
          continue;  // all is another function that will later help me with
                     // sidebands
        const double weight = model.GetWeight(
            *universe, myevent);  // Only calculate the per-universe weight for
                                  // events that will actually use it.
        std::vector<double> q3bin1weights = {
            0.82435, 0.830887, 0.862543, 0.917496, 0.991634, 1.08006,
            1.17502, 1.2697,   1.35885,  1.43734,  1.49575,  1.51875,
            1.47963, 1.34423,  1.13559,  0.918846, 0.788976, 0.735919,
            0.71303, 0.706644, 0.70802,  0.710867, 0.711998};

        std::vector<double> q3bin2weights = {
            0.169059, 0.182445, 0.242976, 0.339379, 0.459126, 0.586182,
            0.708931, 0.82085,  0.924898, 1.03088,  1.14148,  1.24816,
            1.32363,  1.32895,  1.24746,  1.06005,  0.868318, 0.767249,
            0.771477, 0.835023, 0.913111, 0.971778, 0.987021};

        std::vector<double> q3bin3weights = {
            0.406167, 0.400722, 0.381415, 0.359252, 0.345346, 0.355835,
            0.406526, 0.510034, 0.666337, 0.857511, 1.04733,  1.20275,
            1.311,    1.36494,  1.34045,  1.23871,  1.09418,  0.959903,
            0.908177, 0.930722, 0.995452, 1.05769,  1.07703};

        std::vector<double> q3bin4weights = {
            0.465274, 0.475479, 0.521774, 0.59605, 0.685829, 0.781784,
            0.876495, 0.967061, 1.05796,  1.15474, 1.25674,  1.35469,
            1.43084,  1.47119,  1.47522,  1.41,    1.25184,  1.07685,
            0.968015, 0.941743, 0.966135, 1.00764, 1.02239};

        std::vector<double> q3bin5weights = {
            0.549138, 0.562134, 0.624496, 0.72724,  0.859891, 1.00808,
            1.15921,  1.30858,  1.45383,  1.5935,   1.7235,   1.83011,
            1.88988,  1.88183,  1.80408,  1.63456,  1.38423,  1.11548,
            0.902733, 0.778054, 0.732044, 0.731376, 0.738204};

        std::vector<double> tpibin = {.002, .006, .010, .014, .018, .022,
                                      .026, .030, .034, .038, .043, .049,
                                      .061, .075, .090, .125, .175, .225,
                                      .275, .325, .375, .450, .550};
        double weight2 = weight;
        std::vector<double> currentbins;
        double q3_mecAna = universe->Getq3();
        if (q3_mecAna < 0.40)
          currentbins = q3bin1weights;
        else if (q3_mecAna >= 0.40 || q3_mecAna < 0.60)
          currentbins = q3bin2weights;
        else if (q3_mecAna >= 0.60 || q3_mecAna < 0.80)
          currentbins = q3bin3weights;
        else if (q3_mecAna >= 0.80 || q3_mecAna < 1.00)
          currentbins = q3bin4weights;
        else if (q3_mecAna >= 1.00 || q3_mecAna < 1.20)
          currentbins = q3bin5weights;

        int nbins = currentbins.size();

        double tpi = universe->GetTrueTpi();
        /*
          for (int i = 0; i<nbins; i++){
            if (tpi <= tpibin[0]) weight2 = weight*currentbins[0];
            else if (tpi > tpibin[0] && tpi <= 0.500){
              double binedge1 = (tpibin[i]+tpibin[i]);
              double binedge2 = (tpibin[i]+tpibin[i+1]);
              if (tpi > binedge1 && tpi < binedge2 && tpi <
                  0.5) { weight2 = weight*currentbins[i];
              }
            }
            else{
              if (tpi > .500 ) weight2 =
                weight*currentbins[22];
            }
          }
        */
        for (auto& var : vars)
          var->selectedMCReco->FillUniverse(universe,
                                            var->GetRecoValue(*universe),
                                            weight2);  //"Fake data" for closure

        const bool isSignal = michelcuts.isSignal(*universe, weight2);

        if (isSignal) {
          for (auto& study : studies)
            study->SelectedSignal(*universe, myevent, weight2);

          for (auto& var : vars) {
            // Cross section components
            var->efficiencyNumerator->FillUniverse(
                universe, var->GetTrueValue(*universe), weight2);
            var->migration->FillUniverse(universe, var->GetRecoValue(*universe),
                                         var->GetTrueValue(*universe), weight2);
            var->selectedSignalReco->FillUniverse(
                universe, var->GetRecoValue(*universe),
                weight2);  // Efficiency numerator in reco variables.  Useful
                           // for warping studies.
          }

          for (auto& var : vars2D) {
            var->efficiencyNumerator->FillUniverse(
                universe, var->GetTrueValueX(*universe),
                var->GetTrueValueY(*universe), weight2);
          }
        } else {
          int bkgd_ID = -1;
          if (universe->GetCurrent() == 2)
            bkgd_ID = 0;
          else
            bkgd_ID = 1;

          for (auto& var : vars)
            (*var->m_backgroundHists)[bkgd_ID].FillUniverse(
                universe, var->GetRecoValue(*universe), weight2);
          for (auto& var : vars2D)
            (*var->m_backgroundHists)[bkgd_ID].FillUniverse(
                universe, var->GetRecoValueX(*universe),
                var->GetRecoValueY(*universe), weight2);
        }
      }  // End band's universe loop
    }    // End Band loop
  }      // End entries loop
  std::cout << "Finished MC reco loop.\n";
}

void LoopAndFillData(PlotUtils::ChainWrapper* data,
                     std::vector<CVUniverse*> data_band,
                     std::vector<Variable*> vars,
                     std::vector<Variable2D*> vars2D,
                     std::vector<Study*> studies,
                     PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)

{
  std::cout << "Starting data loop...\n";
  const int nEntries = data->GetEntries();
  for (int i = 0; i < data->GetEntries(); ++i) {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if (i % 1000 == 0)
        std::cout << i << " / " << nEntries << "\r" << std::flush;
      // std::cout << "Creating Michel Event" << std::endl;
      MichelEvent myevent;
      // std::cout << "Applying Cuts to Data Event " << std::endl;
      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;
      // std::cout << "Filling Data STudies" << std::endl;
      for (auto& study : studies) study->Selected(*universe, myevent, 1);

      for (auto& var : vars) {
        var->dataHist->FillUniverse(
            universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
      }

      for (auto& var : vars2D) {
        var->dataHist->FillUniverse(universe, var->GetRecoValueX(*universe),
                                    var->GetRecoValueY(*universe), 1);
      }
    }
  }
  std::cout << "Finished data loop.\n";
}

void LoopAndFillEffDenom(
    PlotUtils::ChainWrapper* truth,
    std::map<std::string, std::vector<CVUniverse*>> truth_bands,
    std::vector<Variable*> vars, std::vector<Variable2D*> vars2D,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model) {
  assert(!truth_bands["cv"].empty() &&
         "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();

  std::cout << "Starting efficiency denominator loop...\n";
  const int nEntries = truth->GetEntries();
  for (int i = 0; i < nEntries; ++i) {
    if (i % 1000 == 0)
      std::cout << i << " / " << nEntries << "\r" << std::flush;
    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : truth_bands) {
      std::vector<CVUniverse*> truth_band_universes = band.second;
      for (auto universe : truth_band_universes) {
        MichelEvent myevent;  // Only used to keep the Model happy

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);

        if (!michelcuts.isEfficiencyDenom(*universe, cvWeight))
          continue;  // Weight is ignored for isEfficiencyDenom() in all but the
                     // CV universe
        const double weight = model.GetWeight(
            *universe,
            myevent);  // Only calculate the weight for events that will use it

        // Fill efficiency denominator now:
        for (auto var : vars) {
          var->efficiencyDenominator->FillUniverse(
              universe, var->GetTrueValue(*universe), weight);
        }
        // Fill Studies denominator:
        // for(auto& study: studies) study->SelectedSignal(*universe, cvEvent,
        // weight);

        for (auto var : vars2D) {
          var->efficiencyDenominator->FillUniverse(
              universe, var->GetTrueValueX(*universe),
              var->GetTrueValueY(*universe), weight);
        }
      }
    }
  }
  std::cout << "Finished efficiency denominator loop.\n";
}

// Returns false if recoTreeName could not be inferred
bool inferRecoTreeNameAndCheckTreeNames(const std::string& mcPlaylistName,
                                        const std::string& dataPlaylistName,
                                        std::string& recoTreeName) {
  const std::vector<std::string> knownTreeNames = {"Truth", "Meta", "Header"};
  bool areFilesOK = false;

  std::ifstream playlist(mcPlaylistName);
  std::string firstFile = "";
  playlist >> firstFile;
  auto testFile = TFile::Open(firstFile.c_str());
  if (!testFile) {
    std::cerr << "Failed to open the first MC file at " << firstFile << "\n";
    return false;
  }

  // Does the MC playlist have the Truth tree?  This is needed for the
  // efficiency denominator.
  const auto truthTree = testFile->Get("Truth");
  if (truthTree == nullptr ||
      !truthTree->IsA()->InheritsFrom(TClass::GetClass("TTree"))) {
    std::cerr << "Could not find the \"Truth\" tree in MC file named "
              << firstFile << "\n";
    return false;
  }

  // Figure out what the reco tree name is
  for (auto key : *testFile->GetListOfKeys()) {
    if (static_cast<TKey*>(key)->ReadObj()->IsA()->InheritsFrom(
            TClass::GetClass("TTree")) &&
        std::find(knownTreeNames.begin(), knownTreeNames.end(),
                  key->GetName()) == knownTreeNames.end()) {
      recoTreeName = key->GetName();
      areFilesOK = true;
    }
  }
  delete testFile;
  testFile = nullptr;

  // Make sure the data playlist's first file has the same reco tree
  playlist.open(dataPlaylistName);
  playlist >> firstFile;
  testFile = TFile::Open(firstFile.c_str());
  if (!testFile) {
    std::cerr << "Failed to open the first data file at " << firstFile << "\n";
    return false;
  }

  const auto recoTree = testFile->Get(recoTreeName.c_str());
  if (recoTree == nullptr ||
      !recoTree->IsA()->InheritsFrom(TClass::GetClass("TTree"))) {
    std::cerr << "Could not find the \"" << recoTreeName
              << "\" tree in data file named " << firstFile << "\n";
    return false;
  }

  return areFilesOK;
}

//==============================================================================
// Main
//==============================================================================
int main(const int argc, const char** argv) {
  //============================================================================
  // Setup
  //============================================================================
    TH1::AddDirectory(false);

    // Validate input.
    // I expect a data playlist file name and an MC playlist file name which is
    // exactly 2 arguments.
    const int nArgsExpected = 2;
    // argc is the size of argv.  I check for number of arguments + 1 because
    // argv[0] is always the path to the executable.
    if (argc != nArgsExpected + 1) {
      std::cerr << "Expected " << nArgsExpected << " arguments, but got "
                << argc - 1 << "\n"
                << USAGE << "\n";
      return badCmdLine;
    }

    // One playlist must contain only MC files, and the other must contain only
    // data files. Only checking the first file in each playlist because opening
    // each file an extra time remotely (e.g. through xrootd) can get expensive.
    // TODO: Look in INSTALL_DIR if files not found?
    const std::string mc_file_list = argv[2], data_file_list = argv[1];
    // Check that necessary TTrees exist in the first file of mc_file_list and
    // data_file_list
    std::string reco_tree_name;
    if (!inferRecoTreeNameAndCheckTreeNames(mc_file_list, data_file_list,
                                            reco_tree_name)) {
      std::cerr << "Failed to find required trees in MC playlist " << mc_file_list
                << " and/or data playlist " << data_file_list << ".\n"
                << USAGE << "\n";
      return badInputFile;
    }

    const bool doCCQENuValidation =
        (reco_tree_name == "CCQENu");  // Enables extra histograms and might
                                       // influence which systematics I use.

    const bool is_grid = false;
    PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true, is_grid);
    std::cout << options.m_mc->GetChain()->GetName() << std::endl;
    options.m_plist_string = util::GetPlaylist(*options.m_mc, true);  // TODO: Put GetPlaylist into PlotUtils::MacroUtil

    // You're required to make some decisions
    PlotUtils::MinervaUniverse::SetNuEConstraint(true);
    PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string);  // TODO: Infer this from the files somehow?
    PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
    PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
    PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

    // outfiles
    TFile* mc_MichelStudies = TFile::Open(
        "Dec312021_noreweight_distcuts_q3bin2_atleast1pi_MC.root", "RECREATE");
    TFile* data_MichelStudies = TFile::Open(
        "Dec312021_noreweight_distcuts_q3bin2_atleast1pi_data.root", "RECREATE");

  //============================================================================
  // Cuts
  //============================================================================
    // Now that we've defined what a cross section is, decide which sample and
    // model we're extracting a cross section for.
    PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands, preCuts;
    PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t signalDefinition, phaseSpace;

    const double minZ = 5980, maxZ = 8422, apothem = 850;  // All in mm
    preCuts.emplace_back(
        new reco::ZRange<CVUniverse, MichelEvent>("Tracker", minZ, maxZ));
    preCuts.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
    preCuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(20.));
    preCuts.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>());
    preCuts.emplace_back(
        new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
    preCuts.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>());
    preCuts.emplace_back(new Q3RangeReco<CVUniverse, MichelEvent>(0.4, 0.6));
    // preCuts.emplace_back(new hasMichel<CVUniverse, MichelEvent>());
    // preCuts.emplace_back(new BestMichelDistance2D<CVUniverse,
    // MichelEvent>(100.)); preCuts.emplace_back(new VtxMatchFirst<CVUniverse,
    // MichelEvent>(200., 102.));
    preCuts.emplace_back(new NPiCut<CVUniverse, MichelEvent>(1));
    preCuts.emplace_back(new hasMichel<CVUniverse, MichelEvent>());
    preCuts.emplace_back(new BestMichelDistance2D<CVUniverse, MichelEvent>(150.));

  //============================================================================
  // Signal Definition
  //============================================================================
    signalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
    signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
    signalDefinition.emplace_back(new Q3Limit<CVUniverse>(0.4, 0.6));
    // signalDefinition.emplace_back(new hasTruePion<CVUniverse>());
    phaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Tracker", minZ, maxZ));
    phaseSpace.emplace_back(new truth::Apothem<CVUniverse>(apothem));
    phaseSpace.emplace_back(new truth::MuonAngle<CVUniverse>(20.));
    phaseSpace.emplace_back(new truth::PZMuMin<CVUniverse>(1500.));

  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(preCuts), std::move(sidebands), std::move(signalDefinition), std::move(phaseSpace));

  //============================================================================
  // Weights
  //============================================================================
    std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>> MnvTunev1;
    MnvTunev1.emplace_back(
        new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
    MnvTunev1.emplace_back(
        new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
    MnvTunev1.emplace_back(
        new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
    MnvTunev1.emplace_back(
        new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTunev1.emplace_back(
        new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());

    // TODO: Add my pion reweighter here. - Mehreen S.  Nov 22, 2021
    PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTunev1));

  //============================================================================
  // set up systematics
  //============================================================================
    const bool doSystematics = (getenv("MNV101_SKIP_SYST") == nullptr);
    if (!doSystematics) {
      std::cout << "Skipping systematics (except 1 flux universe) because "
                   "environment variable MNV101_SKIP_SYST is set.\n";
      // Necessary to get Flux integral later...  Doesn't work with just
      // 1 flux universe though because _that_ triggers "spread errors".
      PlotUtils::MinervaUniverse::SetNFluxUniverses(2);  
    }

    std::map<std::string, std::vector<CVUniverse*>> error_bands;
    if (doSystematics) {
      error_bands = GetStandardSystematics(options.m_mc);
    }
    else {
      std::map<std::string, std::vector<CVUniverse*>> band_flux =
          PlotUtils::GetFluxSystematicsMap<CVUniverse>(
              options.m_mc, CVUniverse::GetNFluxUniverses());
      error_bands.insert(
          band_flux.begin(),
          band_flux.end());  // Necessary to get flux integral later...
    }
    error_bands["cv"] = {new CVUniverse(options.m_mc)};
    std::map<std::string, std::vector<CVUniverse*>> truth_bands;
    if (doSystematics) truth_bands = GetStandardSystematics(options.m_truth);
    truth_bands["cv"] = {new CVUniverse(options.m_truth)};

  //============================================================================
  // Set Up a Data Universe
  //============================================================================
    CVUniverse* data_universe = new CVUniverse(options.m_data);
    std::vector<CVUniverse*> data_band = {data_universe};
    std::map<std::string, std::vector<CVUniverse*>> data_error_bands;
    data_error_bands["cv"] = data_band;

  //============================================================================
  // Variables and binning
  //============================================================================
    std::vector<double> dansPTBins = {0,    0.075, 0.15, 0.25, 0.325, 0.4, 0.475,
                                      0.55, 0.7,   0.85, 1,    1.25,  1.5},
                        dansPzBins = {1.5, 2, 2.5, 3,  3.5, 4,  4.5, 5, 6,
                                      7,   8, 9,   10, 15,  20, 40,  60},
                        robsEmuBins = {0,  1,  2,  3,  4,  5,  7,  9,
                                       12, 15, 18, 22, 36, 50, 75, 80},
                        robsRecoilBins;

    std::vector<double> tpibins = {0,    4.,   8.,    12.,  16.,  20.,  24.,
                                   28.,  32.,  36.,   40.,  46.,  52.,  70.,
                                   80.,  100., 150.,  200., 250., 300., 350.,
                                   400., 500., 1000., 2000.};

    const double robsRecoilBinWidth = 50;  // MeV
    for (int whichBin = 0; whichBin < 100 + 1; ++whichBin)
      robsRecoilBins.push_back(robsRecoilBinWidth * whichBin);

    std::vector<Variable*> vars = {
        new Variable("pTmu", "p_{T, #mu} [GeV/c]", dansPTBins,
                     &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue),
        new Variable("q3", "q3 [GeV]", dansPTBins, &CVUniverse::Getq3,
                     &CVUniverse::GetTrueQ3)};

    std::vector<Variable2D*> vars2D;
    if (doCCQENuValidation) {
      std::cerr << "Detected that tree name is CCQENu.  Making validation "
                   "histograms.\n";
      vars.push_back(new Variable("pzmu", "p_{||, #mu} [GeV/c]", dansPzBins,
                                  &CVUniverse::GetMuonPz,
                                  &CVUniverse::GetMuonPzTrue));
      vars.push_back(new Variable("Emu", "E_{#mu} [GeV]", robsEmuBins,
                                  &CVUniverse::GetEmuGeV,
                                  &CVUniverse::GetElepTrueGeV));
      vars.push_back(new Variable(
          "Erecoil", "E_{recoil}", robsRecoilBins, &CVUniverse::GetRecoilE,
          &CVUniverse::GetTrueEAvail));  // TODO: q0 is not the same as recoil
                                         // energy without a spline correction

      vars2D.push_back(new Variable2D(*vars[1], *vars[0]));
    }

  //============================================================================
  // This is where your list of Studies go for PerMichel variables -> Accessed
  // through MichelEvent
  //============================================================================
    std::vector<Study*> studies;

    VarConfig deltat_config{"deltat", "#mus", 30, 0., 9.};
    VarConfig pirange_config{"pirange", "mm", 100, 0.0, 2000.0};
    VarConfig tpi_config{"KE", "meV", 100, 0., 1000.};

    // Variables that we want to plot
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          delta_t = [](const CVUniverse& univ, const MichelEvent& evt,
                       const int whichMichel) {
            int evttype = evt.eventtype;
            double micheltime = evt.m_nmichels[whichMichel].time;
            double vtxtime = univ.GetVertex().t();
            double deltat =
                (micheltime -
                 vtxtime / 1000.);  // hopefully this is in microseconds (mus)
            // if (evttype == 1) return deltat;
            // else return -9999.;
            return deltat;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          permichel_range = [](const CVUniverse& univ, const MichelEvent& evt,
                               const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            // if (evt.eventtype == 1) return micheldist;
            // else return -9999.;
            return micheldist;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          pertruepimichel_range = [](const CVUniverse& univ, const MichelEvent& evt,
                                     const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (evt.m_nmichels[whichMichel].true_parentpdg == 211)
              return micheldist;
            else
              return -9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          permichel_tpi = [](const CVUniverse& univ, const MichelEvent& evt,
                             const int whichMichel) {
            return evt.m_nmichels[whichMichel].pionKE;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          overlay_vtx_tpi = [](const CVUniverse& univ, const MichelEvent& evt,
                               const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            double KE = evt.m_nmichels[whichMichel].pionKE;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            if (overlayfrac > .5 && (matchtype == 1 || matchtype == 2))
              return KE;
            else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          overlay_clus_tpi = [](const CVUniverse& univ, const MichelEvent& evt,
                                const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            double KE = evt.m_nmichels[whichMichel].pionKE;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            if (overlayfrac > .5 && (matchtype == 3 || matchtype == 4))
              return KE;
            else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          truemichel_goodvtx_tpi = [](const CVUniverse& univ,
                                      const MichelEvent& evt,
                                      const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            double KE = evt.m_nmichels[whichMichel].pionKE;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
            int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
            if (overlayfrac < .5 && trueEnd == recoEnd &&
                (matchtype == 1 || matchtype == 2))
              return KE;
            else
              return -9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          truemichel_goodclus_tpi = [](const CVUniverse& univ,
                                       const MichelEvent& evt,
                                       const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            double KE = evt.m_nmichels[whichMichel].pionKE;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
            int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
            if (overlayfrac < .5 && trueEnd == recoEnd &&
                (matchtype == 3 || matchtype == 4))
              return KE;
            else
              return -9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          truemichel_badvtx_tpi = [](const CVUniverse& univ, const MichelEvent& evt,
                                     const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            double KE = evt.m_nmichels[whichMichel].pionKE;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
            int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
            if (overlayfrac < .5 && trueEnd != recoEnd &&
                (matchtype == 1 || matchtype == 2)) {
              // univ.PrintArachneLink();
              // std::cout << "Printing Michel Time for bad VTX match type "  <<
              // evt.m_nmichels[whichMichel].time << std::endl;
              return KE;
            } else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          truemichel_badclus_tpi = [](const CVUniverse& univ,
                                      const MichelEvent& evt,
                                      const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            double KE = evt.m_nmichels[whichMichel].pionKE;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
            int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
            if (overlayfrac < .5 && trueEnd != recoEnd &&
                (matchtype == 3 || matchtype == 4)) {
              // univ.PrintArachneLink();
              // std::cout << "Printing Michel Time for bad CLUS match type "  <<
              // evt.m_nmichels[whichMichel].time << std::endl;
              return KE;
            } else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          overlay_vtx_range = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;

            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            if (overlayfrac > .5 && (matchtype == 1 || matchtype == 2))
              return micheldist;
            else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          overlay_clus_range = [](const CVUniverse& univ, const MichelEvent& evt,
                                  const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;

            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            if (overlayfrac > .5 && (matchtype == 3 || matchtype == 4))
              return micheldist;
            else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          truemichel_goodvtx_range = [](const CVUniverse& univ,
                                        const MichelEvent& evt,
                                        const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;

            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
            int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
            if (overlayfrac < .5 && trueEnd == recoEnd &&
                (matchtype == 1 || matchtype == 2))
              return micheldist;
            else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)> truemichel_goodclus_range 
                                    = [](const CVUniverse& univ,
                                         const MichelEvent& evt,
                                         const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
            int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
            if (overlayfrac < .5 && trueEnd == recoEnd &&
                (matchtype == 3 || matchtype == 4))
              return micheldist;
            else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          truemichel_badvtx_range = [](const CVUniverse& univ,
                                       const MichelEvent& evt,
                                       const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
            int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
            if (overlayfrac < .5 && trueEnd != recoEnd &&
                (matchtype == 1 || matchtype == 2)) {
              //      univ.PrintArachneLink();
              //	std::cout << "Printing Michel Time for bad VTX match type "  <<
              //evt.m_nmichels[whichMichel].time << std::endl;
              return micheldist;
            } else
              return -9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          truemichel_badclus_range = [](const CVUniverse& univ,
                                        const MichelEvent& evt,
                                        const int whichMichel) {
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            int overlayfrac = evt.m_nmichels[whichMichel].overlay_fraction;
            int matchtype = evt.m_nmichels[whichMichel].BestMatch;
            int trueEnd = evt.m_nmichels[whichMichel].trueEndpoint;
            int recoEnd = evt.m_nmichels[whichMichel].recoEndpoint;
            if (overlayfrac < .5 && trueEnd != recoEnd &&
                (matchtype == 3 || matchtype == 4)) {
              //	 univ.PrintArachneLink();
              //       std::cout << "Printing Michel Time for bad CLUS match type "
              //       << evt.m_nmichels[whichMichel].time << std::endl;
              return micheldist;
            } else
              return -9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_XZ = [](const CVUniverse& univ, const MichelEvent& evt,
                         const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].best_XZ;
            return twoDdist;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_UZ = [](const CVUniverse& univ, const MichelEvent& evt,
                         const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].best_UZ;
            return twoDdist;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_VZ = [](const CVUniverse& univ, const MichelEvent& evt,
                         const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].best_VZ;
            return twoDdist;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_XZ_upvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                               const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].up_to_vertex_XZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 1)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_UZ_upvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                               const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].up_to_vertex_UZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 1)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_VZ_upvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                               const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].up_to_vertex_VZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 1)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_XZ_upclus = [](const CVUniverse& univ, const MichelEvent& evt,
                                const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].up_to_clus_XZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 1)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_UZ_upclus = [](const CVUniverse& univ, const MichelEvent& evt,
                                const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].up_to_clus_UZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 1)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_VZ_upclus = [](const CVUniverse& univ, const MichelEvent& evt,
                                const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].up_to_clus_VZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 1)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_XZ_downclus = [](const CVUniverse& univ, const MichelEvent& evt,
                                  const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].down_to_clus_XZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 2)
              return twoDdist;
            else
              return 9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_UZ_downclus = [](const CVUniverse& univ, const MichelEvent& evt,
                                  const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].down_to_clus_UZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 2)
              return twoDdist;
            else
              return 9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_VZ_downclus = [](const CVUniverse& univ, const MichelEvent& evt,
                                  const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].down_to_clus_VZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 2)
              return twoDdist;
            else
              return 9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_XZ_downvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].down_to_vertex_XZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 2)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_UZ_downvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].down_to_vertex_UZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 2)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          michel_VZ_downvtx = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double twoDdist = evt.m_nmichels[whichMichel].down_to_vertex_VZ;
            int trueend = evt.m_nmichels[whichMichel].trueEndpoint;
            if (trueend == 2)
              return twoDdist;
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          pion_angle = [](const CVUniverse& univ, const MichelEvent& evt,
                          const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].best_angle;
            return cos(angle);
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          pion_angle_range1 = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].best_angle;
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (micheldist <= 150.)
              return cos(angle);
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          pion_angle_range2 = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].best_angle;
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (micheldist > 150. && micheldist <= 250.)
              return cos(angle);
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          pion_angle_range3 = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].best_angle;
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (micheldist > 250. && micheldist <= 500.)
              return cos(angle);
            else
              return 9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          pion_angle_range4 = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].best_angle;
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (micheldist > 500.)
              return cos(angle);
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          true_angle = [](const CVUniverse& univ, const MichelEvent& evt,
                          const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].true_angle;
            return cos(angle);
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          true_angle_range1 = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].true_angle;
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (micheldist <= 150.)
              return cos(angle);
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          true_angle_range2 = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].true_angle;
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (micheldist > 150 && micheldist <= 250.)
              return cos(angle);
            else
              return 9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          true_angle_range3 = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].true_angle;
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (micheldist > 250 && micheldist <= 500.)
              return cos(angle);
            else
              return 9999.;
          };

      std::function<double(const CVUniverse&, const MichelEvent&, const int)>
          true_angle_range4 = [](const CVUniverse& univ, const MichelEvent& evt,
                                 const int whichMichel) {
            double angle = evt.m_nmichels[whichMichel].true_angle;
            double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
            if (micheldist > 500)
              return cos(angle);
            else
              return 9999.;
          };
      std::function<double(const CVUniverse&, const MichelEvent&)> best_pionrange =
          [](const CVUniverse& univ, const MichelEvent& evt) {
            int bestidx = evt.m_idx;
            if (bestidx < 0)
              return -9999.;
            else {
              double dist = evt.m_nmichels[bestidx].Best3Ddist;
              return dist;
            }
          };

      std::function<double(const CVUniverse&, const MichelEvent&)>
          best_pionrange_overlay_vtx =
              [](const CVUniverse& univ, const MichelEvent& evt) {
                int bestidx = evt.m_idx;
                if (bestidx < 0)
                  return -9999.;
                else {
                  double micheldist = evt.m_nmichels[bestidx].Best3Ddist;
                  int overlayfrac = evt.m_nmichels[bestidx].overlay_fraction;
                  int matchtype = evt.m_nmichels[bestidx].BestMatch;
                  int trueEnd = evt.m_nmichels[bestidx].trueEndpoint;
                  int recoEnd = evt.m_nmichels[bestidx].recoEndpoint;
                  if (overlayfrac > .5 && (matchtype == 1 || matchtype == 2))
                    return micheldist;
                  else
                    return -9999.;
                }
              };

      std::function<double(const CVUniverse&, const MichelEvent&)>
          best_pionrange_overlay_clus =
              [](const CVUniverse& univ, const MichelEvent& evt) {
                int bestidx = evt.m_idx;
                if (bestidx < 0)
                  return -9999.;
                else {
                  double micheldist = evt.m_nmichels[bestidx].Best3Ddist;
                  int overlayfrac = evt.m_nmichels[bestidx].overlay_fraction;
                  int matchtype = evt.m_nmichels[bestidx].BestMatch;
                  int trueEnd = evt.m_nmichels[bestidx].trueEndpoint;
                  int recoEnd = evt.m_nmichels[bestidx].recoEndpoint;
                  if (overlayfrac > .5 && (matchtype == 3 || matchtype == 4))
                    return micheldist;
                  else
                    return -9999.;
                }
              };

      std::function<double(const CVUniverse&, const MichelEvent&)>
          best_pionrange_truegood_vtx =
              [](const CVUniverse& univ, const MichelEvent& evt) {
                int bestidx = evt.m_idx;
                if (bestidx < 0)
                  return -9999.;
                else {
                  double micheldist = evt.m_nmichels[bestidx].Best3Ddist;
                  int overlayfrac = evt.m_nmichels[bestidx].overlay_fraction;
                  int matchtype = evt.m_nmichels[bestidx].BestMatch;
                  int trueEnd = evt.m_nmichels[bestidx].trueEndpoint;
                  int recoEnd = evt.m_nmichels[bestidx].recoEndpoint;
                  if (overlayfrac < .5 && trueEnd == recoEnd &&
                      (matchtype == 1 || matchtype == 2))
                    return micheldist;
                  else
                    return -9999.;
                }
              };

      std::function<double(const CVUniverse&, const MichelEvent&)>
          best_pionrange_truebad_vtx =
              [](const CVUniverse& univ, const MichelEvent& evt) {
                int bestidx = evt.m_idx;
                if (bestidx < 0)
                  return -9999.;
                else {
                  double micheldist = evt.m_nmichels[bestidx].Best3Ddist;
                  int overlayfrac = evt.m_nmichels[bestidx].overlay_fraction;
                  int matchtype = evt.m_nmichels[bestidx].BestMatch;
                  int trueEnd = evt.m_nmichels[bestidx].trueEndpoint;
                  int recoEnd = evt.m_nmichels[bestidx].recoEndpoint;
                  if (overlayfrac < .5 && trueEnd != recoEnd &&
                      (matchtype == 1 || matchtype == 2))
                    return micheldist;
                  else
                    return -9999.;
                }
              };

      std::function<double(const CVUniverse&, const MichelEvent&)>
          best_pionrange_truegood_clus =
              [](const CVUniverse& univ, const MichelEvent& evt) {
                int bestidx = evt.m_idx;
                if (bestidx < 0)
                  return -9999.;
                else {
                  double micheldist = evt.m_nmichels[bestidx].Best3Ddist;
                  int overlayfrac = evt.m_nmichels[bestidx].overlay_fraction;
                  int matchtype = evt.m_nmichels[bestidx].BestMatch;
                  int trueEnd = evt.m_nmichels[bestidx].trueEndpoint;
                  int recoEnd = evt.m_nmichels[bestidx].recoEndpoint;
                  if (overlayfrac < .5 && trueEnd == recoEnd &&
                      (matchtype == 3 || matchtype == 4))
                    return micheldist;
                  else
                    return -9999.;
                }
              };

      std::function<double(const CVUniverse&, const MichelEvent&)>
          best_pionrange_truebad_clus =
              [](const CVUniverse& univ, const MichelEvent& evt) {
                int bestidx = evt.m_idx;
                if (bestidx < 0)
                  return -9999.;
                else {
                  double micheldist = evt.m_nmichels[bestidx].Best3Ddist;
                  int overlayfrac = evt.m_nmichels[bestidx].overlay_fraction;
                  int matchtype = evt.m_nmichels[bestidx].BestMatch;
                  int trueEnd = evt.m_nmichels[bestidx].trueEndpoint;
                  int recoEnd = evt.m_nmichels[bestidx].recoEndpoint;
                  if (overlayfrac < .5 && trueEnd != recoEnd &&
                      (matchtype == 3 || matchtype == 4))
                    return micheldist;
                  else
                    return -9999.;
                }
              };

      std::function<double(const CVUniverse&, const MichelEvent&)> lowesttpi =
          [](const CVUniverse& univ, const MichelEvent& evt) {
            double lowtpi = evt.lowTpi;
            return lowtpi;
          };

    // Fill studies, use the above-defined functions
      studies.push_back(new PerMichelVarByGENIELabel(
          true_angle, "true_angle", "cos(#theta)", 21, -1.0, 1.0, error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(true_angle_range1, "true_angle_range1",
                                       "cos(#theta)", 21, -1.0, 1., error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(true_angle_range2, "true_angle_range2",
                                       "cos(#theta)", 21, -1.0, 1., error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(true_angle_range3, "true_angle_range3",
                                       "cos(#theta)", 21, -1.0, 1., error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(true_angle_range4, "true_angle_range4",
                                       "cos(#theta)", 21, -1.0, 1., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          pion_angle, "pion_angle", "cos(#theta)", 21, -1.0, 1., error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(pion_angle_range1, "pion_angle_range1",
                                       "cos(#theta)", 21, -1.0, 1., error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(pion_angle_range2, "pion_angle_range2",
                                       "cos(#theta)", 21, -1.0, 1., error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(pion_angle_range3, "pion_angle_range3",
                                       "cos(#theta)", 21, -1.0, 1., error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(pion_angle_range4, "pion_angle_range4",
                                       "cos(#theta)", 21, -1.0, 1., error_bands));

      studies.push_back(new PerMichelVarByGENIELabel(michel_XZ, "best_XZ", "mm",
                                                     100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(michel_UZ, "best_UZ", "mm",
                                                     100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(michel_VZ, "best_VZ", "mm",
                                                     100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_XZ_upvtx, "upvtx_XZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_XZ_upclus, "upclus_XZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_XZ_downclus, "downclus_XZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_XZ_downvtx, "downvtx_XZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_UZ_upvtx, "upvtx_UZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_UZ_upclus, "upclus_UZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_UZ_downclus, "downclus_UZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_UZ_downvtx, "downvtx_UZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_VZ_upvtx, "upvtx_VZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_VZ_upclus, "upclus_VZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_VZ_downclus, "downclus_VZ", "mm", 100, 0.0, 2000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          michel_VZ_downvtx, "downvtx_VZ", "mm", 100, 0.0, 2000., error_bands));

      studies.push_back(new PerMichelVarByGENIELabel(
          delta_t, "michelmuon_deltat", "#mus", 30, 0.0, 9.0, error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(permichel_range,
                                                     "permichel_pirange", "mm", 100,
                                                     0.0, 2000.0, error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(overlay_vtx_range,
                                                     "overlay_vtx_range", "mm", 100,
                                                     0.0, 2000.0, error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(overlay_clus_range, "overlay_clus_range",
                                       "mm", 100, 0.0, 2000.0, error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          truemichel_goodvtx_range, "truemichel_goodvtx_range", "mm", 100, 0.0,
          2000.0, error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          truemichel_goodclus_range, "truemichel_goodclus_range", "mm", 100, 0.0,
          2000.0, error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          truemichel_badvtx_range, "truemichel_badvtx_range", "mm", 100, 0.0,
          2000.0, error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          truemichel_badclus_range, "truemichel_badclus_range", "mm", 100, 0.0,
          2000.0, error_bands));
      studies.push_back(
          new PerMichelVarByGENIELabel(permichel_tpi, "Per_Michel_PrimaryParentKE",
                                       "meV", 100, 0, 1000., error_bands));
      studies.push_back(new PerMichelVarByGENIELabel(
          pertruepimichel_range, "permichel_pirange_truepi", "mm", 100, 0.0, 2000.0,
          error_bands));

      // studies.push_back(new PerMichel2DVar(delta_t, permichel_range,
      // deltat_config, pirange_config, error_bands));
      studies.push_back(new PerMichel2DVar(
          permichel_tpi, permichel_range, tpi_config, pirange_config, error_bands));
      // studies.push_back(new PerMichel2DVar(overlay_vtx_tpi, overlay_vtx_range,
      // tpi_config, pirange_config, error_bands)); studies.push_back(new
      // PerMichel2DVar(overlay_clus_tpi, overlay_clus_range, tpi_config,
      // pirange_config, error_bands)); studies.push_back(new
      // PerMichel2DVar(truemichel_goodvtx_range, truemichel_goodvtx_tpi,
      // tpi_config, pirange_config, error_bands)); studies.push_back(new
      // PerMichel2DVar(truemichel_goodclus_range, truemichel_goodclus_tpi,
      // tpi_config, pirange_config, error_bands)); studies.push_back(new
      // PerMichel2DVar(truemichel_badvtx_range, truemichel_badvtx_tpi, tpi_config,
      // pirange_config, error_bands)); studies.push_back(new
      // PerMichel2DVar(truemichel_badclus_range, truemichel_badclus_tpi,
      // tpi_config, pirange_config, error_bands));
      // Studies for Per Event Variables (Best Michel In Event)

      studies.push_back(new PerMichelEventVarByGENIELabel(
          lowesttpi, "LowestKE_pion", "MeV", 100, 0, 1000., error_bands));

      studies.push_back(new PerMichelEventVarByGENIELabel(
          best_pionrange, "best_pionrange", "mm", 100, 0.0, 2000.0, error_bands));
      studies.push_back(new PerMichelEventVarByGENIELabel(
          best_pionrange_overlay_clus, "best_pionrange_overlay_clus", "mm", 100,
          0.0, 2000.0, error_bands));
      studies.push_back(new PerMichelEventVarByGENIELabel(
          best_pionrange_overlay_vtx, "best_pionrange_overlay_vtx ", "mm", 100, 0.0,
          2000.0, error_bands));
      studies.push_back(new PerMichelEventVarByGENIELabel(
          best_pionrange_truegood_vtx, "best_pionrange_truegood_vtx", "mm", 100,
          0.0, 2000.0, error_bands));
      studies.push_back(new PerMichelEventVarByGENIELabel(
          best_pionrange_truebad_vtx, "best_pionrange_truebad_vtx", "mm", 100, 0.0,
          2000.0, error_bands));
      studies.push_back(new PerMichelEventVarByGENIELabel(
          best_pionrange_truegood_clus, "best_pionrange_truegood_clus", "mm", 100,
          0.0, 2000.0, error_bands));
      studies.push_back(new PerMichelEventVarByGENIELabel(
          best_pionrange_truebad_clus, "best_pionrange_truebad_clus", "mm", 100,
          0.0, 2000.0, error_bands));

    // Data studies
      std::vector<Study*> data_studies;

      data_studies.push_back(
          new PerMichelEventVarByGENIELabel(best_pionrange, "best_pionrange", "mm",
                                            100, 0.0, 2000.0, data_error_bands));
      data_studies.push_back(new PerMichelVarByGENIELabel(
          delta_t, "michelmuon_deltat", "#mus", 30, 0.0, 9.0, data_error_bands));
      data_studies.push_back(
          new PerMichelVarByGENIELabel(permichel_range, "permichel_pirange", "mm",
                                       100, 0.0, 2000.0, data_error_bands));
      data_studies.push_back(new PerMichelVarByGENIELabel(
          pion_angle, "pion_angle", "cos(#theta)", 21, -1.0, 1., data_error_bands));
      data_studies.push_back(new PerMichelVarByGENIELabel(
          pion_angle_range1, "pion_angle_range1", "cos(#theta)", 21, -1.0, 1.,
          data_error_bands));
      data_studies.push_back(new PerMichelVarByGENIELabel(
          pion_angle_range2, "pion_angle_range2", "cos(#theta)", 21, -1.0, 1.,
          data_error_bands));
      data_studies.push_back(new PerMichelVarByGENIELabel(
          pion_angle_range3, "pion_angle_range3", "cos(#theta)", 21, -1.0, 1.,
          data_error_bands));
      data_studies.push_back(new PerMichelVarByGENIELabel(
          pion_angle_range4, "pion_angle_range4", "cos(#theta)", 21, -1.0, 1.,
          data_error_bands));

  //============================================================================
  // Initialize Hists
  //============================================================================
    for (auto& var : vars) var->InitializeMCHists(error_bands, truth_bands);
    for (auto& var : vars) var->InitializeDATAHists(data_band);

    for (auto& var : vars2D) var->InitializeMCHists(error_bands, truth_bands);
    for (auto& var : vars2D) var->InitializeDATAHists(data_band);

  //============================================================================
  // Loop entries and fill
  //============================================================================
    try {
      CVUniverse::SetTruth(false);
      LoopAndFillEventSelection(options.m_mc, error_bands, vars, vars2D, studies,
                                mycuts, model);
      CVUniverse::SetTruth(true);
      LoopAndFillEffDenom(options.m_truth, truth_bands, vars, vars2D, studies,
                          mycuts, model);
      options.PrintMacroConfiguration(argv[0]);
      std::cout << "MC cut summary:\n" << mycuts << "\n";
      mycuts.resetStats();

      // CVUniverse::SetTruth(false);
      // LoopAndFillData(options.m_data, data_band, vars, vars2D, data_studies,
      // mycuts); std::cout << "Data cut summary:\n" << mycuts << "\n";

      // Write MC results
      TFile* mcOutDir = TFile::Open(MC_OUT_FILE_NAME, "RECREATE");
      if (!mcOutDir) {
        std::cerr << "Failed to open a file named " << MC_OUT_FILE_NAME
                  << " in the current directory for writing histograms.\n";
        return badOutputFile;
      }

      // TFile* mc_MichelStudies =
      // TFile::Open("AllMichel_hasFittedMichel_500mm.root", "RECREATE");
      //"ALL2DDistprinted_OnlyPionMichels_tpimorethan80meV_forceendpointmatch_2Ddistcut_mc.root",
      //"RECREATE");
      for (auto& study : studies) study->SaveOrDraw(*mc_MichelStudies);
      for (auto& var : vars) var->WriteMC(*mcOutDir);
      for (auto& var : vars2D) var->Write(*mcOutDir);
      for (auto& study : data_studies) study->SaveOrDraw(*data_MichelStudies);
      // Protons On Target
      auto mcPOT = new TParameter<double>("POTUsed", options.m_mc_pot);
      mcPOT->Write();
      mc_MichelStudies->cd();
      mcPOT->Write();

      PlotUtils::TargetUtils targetInfo;
      assert(error_bands["cv"].size() == 1 &&
             "List of error bands must contain a universe named \"cv\" for the "
             "flux integral.");

      for (const auto& var : vars) {
        // Flux integral only if systematics are being done (temporary solution)
        util::GetFluxIntegral(*error_bands["cv"].front(),
                              var->efficiencyNumerator->hist)
            ->Write((var->GetName() + "_reweightedflux_integrated").c_str());
        // Always use MC number of nucleons for cross section
        auto nNucleons = new TParameter<double>(
            (var->GetName() + "_fiducial_nucleons").c_str(),
            targetInfo.GetTrackerNNucleons(minZ, maxZ, true, apothem));
        nNucleons->Write();
      }

      // Write data results
      TFile* dataOutDir = TFile::Open(DATA_OUT_FILE_NAME, "RECREATE");
      if (!dataOutDir) {
        std::cerr << "Failed to open a file named " << DATA_OUT_FILE_NAME
                  << " in the current directory for writing histograms.\n";
        return badOutputFile;
      }

      for (auto& var : vars) var->WriteData(*dataOutDir);

      // Protons On Target
      auto dataPOT = new TParameter<double>("POTUsed", options.m_data_pot);
      dataPOT->Write();
      data_MichelStudies->cd();
      dataPOT->Write();
      std::cout << "Success" << std::endl;
    } catch (const ROOT::exception& e) {
      std::cerr
          << "Ending on a ROOT error message.  No histograms will be produced.\n"
          << "If the message talks about \"TNetXNGFile\", this could be a "
             "problem with dCache.  The message is:\n"
          << e.what() << "\n"
          << USAGE << "\n";
      return badFileRead;
    }

  return success;
}

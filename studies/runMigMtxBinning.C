
#ifndef runMigMtxBinning_C
#define runMigMtxBinning_C

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TTree.h"
#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif  // MNVROOT6

#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/MacroUtil.h"
#include "includes/CVUniverse.h"
#include "includes/Cuts.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "includes/Variable2D.h"
#include "includes/common_functions.h"
#include "../xsec/makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "../xsec/plotting_functions.h"
#include "../xsec/plotting_functions2D.h"

typedef std::vector<vector<double>> Matrix;


//==============================================================================
// Loop and fill Tree with true and false points
//==============================================================================
void Loop(const CCPi::MacroUtil& util, CVUniverse* universe,
          const EDataMCTruth& type,
          std::vector<Variable2D*>& variables2D,
 	  std::string fName,
          std::string treeName, std::string xvar, std::string yvar){ 
  TFile f(Form("%s.root",fName.c_str()),"recreate");
  typedef struct {Double_t vxt, vxr, vyt, vyr;} POINT;
  POINT point;
  Double_t pntpmuthetamu[4];
  Double_t pntptmutpi[4];
  Double_t pntpzmuptmu[4];
  Double_t pnttpipmu[4];
  Double_t pnttpithetapi[4];
  TTree tree(treeName.c_str(), "Tree with the MigMtx points");
  tree.Branch("pmu_vs_thetamu_deg", &pntpmuthetamu, "pntpmuthetamu[4]/D");
  tree.Branch("ptmu_vs_tpi", &pntptmutpi, "pntptmutpi[4]/D");
  tree.Branch("pzmu_vs_ptmu", &pntpzmuptmu, "pntpzmuptmu[4]/D");
  tree.Branch("tpi_vs_pmu", &pnttpipmu, "pnttpipmu[4]/D");
  tree.Branch("tpi_vs_thetapi_deg", &pnttpithetapi, "pnttpithetapi[4]/D");

  bool is_mc = true, is_truth = false;
  Long64_t n_entries = util.GetMCEntries();
  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % (n_entries / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);
//    if (i_event == 4000) break;
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    
    PassesCutsInfo cuts_info = PassesCuts(event);
    std::tie(event.m_passes_cuts, event.m_is_w_sideband, event.m_passes_all_cuts_except_w, event.m_reco_pion_candidate_idxs) = cuts_info.GetAll();
    event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);

    universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);

    event.m_weight = universe->GetWeight();

    if (event.m_is_signal && event.m_passes_cuts) {
//        std::cout << "simon\n";
    //for (auto var2D : variables2D){
//        if (var2D->m_is_true) continue;      
        Variable2D* reco_pmuthetamu = GetVar2D(variables2D, "pmu", "thetamu_deg");
        Variable2D* true_pmuthetamu = GetVar2D(variables2D, "pmu_true", "thetamu_deg_true");
        if (true_pmuthetamu == 0) continue;
        Variable2D* reco_pzmuptmu = GetVar2D(variables2D, "pzmu", "ptmu");
        Variable2D* true_pzmuptmu = GetVar2D(variables2D, "pzmu_true", "ptmu_true");
        if (true_pzmuptmu == 0) continue;
        Variable2D* reco_tpithetapi = GetVar2D(variables2D, "tpi", "thetapi_deg");
        Variable2D* true_tpithetapi = GetVar2D(variables2D, "tpi_true", "thetapi_deg_true");
        if (true_tpithetapi == 0) continue;
        RecoPionIdx reco_idx = event.m_highest_energy_pion_idx;
        TruePionIdx true_idx = event.m_universe->GetHighestEnergyTruePionIndex();

        double reco_fill_pmu = reco_pmuthetamu->GetValueX(*event.m_universe, reco_idx);
        double reco_fill_thetamu = reco_pmuthetamu->GetValueY(*event.m_universe, reco_idx);
        double reco_fill_ptmu = reco_pzmuptmu->GetValueY(*event.m_universe, reco_idx);
        double reco_fill_pzmu = reco_pzmuptmu->GetValueX(*event.m_universe, reco_idx);
        double reco_fill_tpi = reco_tpithetapi->GetValueY(*event.m_universe, reco_idx);
        double reco_fill_thetapi = reco_tpithetapi->GetValueY(*event.m_universe, reco_idx);

        double true_fill_pmu = true_pmuthetamu->GetValueX(*event.m_universe, true_idx);
        double true_fill_thetamu = true_pmuthetamu->GetValueY(*event.m_universe, true_idx);
        double true_fill_ptmu = true_pzmuptmu->GetValueY(*event.m_universe, true_idx);
        double true_fill_pzmu = true_pzmuptmu->GetValueX(*event.m_universe, true_idx);
        double true_fill_tpi = true_tpithetapi->GetValueY(*event.m_universe, true_idx);
        double true_fill_thetapi = true_tpithetapi->GetValueY(*event.m_universe, true_idx);
        pntpmuthetamu[0] = true_fill_pmu*event.m_weight;
        pntpmuthetamu[1] = reco_fill_pmu*event.m_weight;
        pntpmuthetamu[2] = true_fill_thetamu*event.m_weight;
        pntpmuthetamu[3] = reco_fill_thetamu*event.m_weight;
 
        pntptmutpi[0] = true_fill_ptmu*event.m_weight;
        pntptmutpi[1] = reco_fill_ptmu*event.m_weight;
        pntptmutpi[2] = true_fill_tpi*event.m_weight;
        pntptmutpi[3] = reco_fill_tpi*event.m_weight;

        pntpzmuptmu[0] = true_fill_pzmu*event.m_weight;
        pntpzmuptmu[1] = reco_fill_pzmu*event.m_weight;
        pntpzmuptmu[2] = true_fill_ptmu*event.m_weight;
        pntpzmuptmu[3] = reco_fill_ptmu*event.m_weight;

        pnttpipmu[0] = true_fill_tpi*event.m_weight;
        pnttpipmu[1] = reco_fill_tpi*event.m_weight;
        pnttpipmu[2] = true_fill_pmu*event.m_weight;
        pnttpipmu[3] = reco_fill_pmu*event.m_weight;

        pnttpithetapi[0] = true_fill_tpi*event.m_weight;
        pnttpithetapi[1] = reco_fill_tpi*event.m_weight;
        pnttpithetapi[2] = true_fill_thetapi*event.m_weight;
        pnttpithetapi[3] = reco_fill_thetapi*event.m_weight;
      //  std::cout << pnt[0] << " " << pnt[1] << " " << pnt[2] << " " << pnt[3] << "\n";
        tree.Fill();
    //}
    }

  }  // events
  f.cd();
  tree.Write();
}

void LoopforTrueEvents(std::string MigFileName, std::string TreeName, std::string tag,
                       int NBinsx, int NBinsy, 
                       TArrayD binedgesx, TArrayD binedgesy){

  TFile *MigFile = new TFile(Form("%s.root", MigFileName.c_str()),"update");
  TTree *MigTree = (TTree*)MigFile->Get(TreeName.c_str());
  Int_t cord[3];
  Double_t pnt[4];
  TBranch *bTrueEvbin = MigTree->Branch(Form("TrueEventBinningLocation%s",tag.c_str()),&cord,"cord[3]/I");
  TBranch *MigVal = MigTree->GetBranch("pmu_vs_theta_mu");
  MigVal->SetAddress(&pnt);

  Long64_t nentries = MigTree->GetEntries();
  std::cout << "Entries = " << nentries << "\n";
  for (Long64_t entry = 0;entry < nentries; entry++){
    MigVal->GetEntry(entry);
//    if (entry == 10) break;
    int i, j, x, y;
    for(i = 0; i < NBinsx; ++i){
      if (pnt[0] > binedgesx[i] && pnt[0] < binedgesx[i+1]){
        x = i + 1;
        break;
      }
    } 
    for(j = 0; j < NBinsy; ++j){     
      if (pnt[2] > binedgesy[j] && pnt[2] < binedgesy[j+1]){
        y = j + 1;
        break;
      }
    }
    cord[0] = entry;
    cord[1] = x;
    cord[2] = y;
//    std::cout << cord[0] << " " << cord[1] << " " << cord[2] << "\n";
    bTrueEvbin->Fill();   

  }
  MigFile->cd();
  MigTree->Write();
}

//template <int xbins, int ybins>
std::vector<double> byStaError(std::vector<double> xbinedges, std::vector<double> ybinedges,
                int& xnbins, int& ynbins, Matrix& truemtrx, int x, int y, double nentries){
//  for (int i = 0; i <= ynbins; ++i){
//    std::cout << "que pinche pablo ybinedges = " << ybinedges[i] << "\n";
//  }

  double currval = nentries, val;
  int xhigher, yhigher, z = 0; 
  std::vector <double> NewBinning;
  if (x > 0 && y > 0 && x < xnbins - 1 && y < ynbins - 1){
    for (int j = y - 1; j < y + 2; j++){
      for(int i = x - 1; i < x + 2; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ){ z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }

    if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
    if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}
    if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
    if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 8){
      xbinedges.erase(xbinedges.begin()+x+1); 
      xbinedges.erase(xbinedges.begin()+x); 
      ybinedges.erase(ybinedges.begin()+y);
      ybinedges.erase(ybinedges.begin()+y+1);
      ynbins = ynbins -2;
      xnbins = xnbins -2;
    } 
  }

  else if (x == 0 && y == 0){
    for (int j = y; j < y + 2; j++){
      for(int i = x; i < x + 2; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ) { z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }
    
    if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}

    if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 3){
      xbinedges.erase(xbinedges.begin()+x+1); 
      ybinedges.erase(ybinedges.begin()+y+1);
      ynbins = ynbins -1;
      xnbins = xnbins -1;
    } 
  }

  else if (x > 0 && y == 0 && x < xnbins - 1){
    for (int j = y; j < y + 2; j++){
      for(int i = x - 1; i < x + 2; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ) { z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }
    
    if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
    if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}
    if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 5){
      xbinedges.erase(xbinedges.begin()+x+1); 
      xbinedges.erase(xbinedges.begin()+x); 
      ybinedges.erase(ybinedges.begin()+y+1);
      ynbins = ynbins -1;
      xnbins = xnbins -2;
    }
  } 

  else if (y == 0 && x == xnbins - 1){
    for (int j = y; j < y + 2; j++){
      for(int i = x - 1; i < x + 1; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ) { z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }
    
    if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
    if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 3){
      xbinedges.erase(xbinedges.begin()+x);
      ybinedges.erase(ybinedges.begin()+y+1);
      ynbins = ynbins -1;
      xnbins = xnbins -1;
    }
  }

  else if ( y > 0 && x == xnbins - 1 && y < ynbins - 1){
    for (int j = y - 1; j < y + 2; j++){
      for(int i = x - 1; i < x + 1; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ) { z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }
    if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
    if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
    if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 5){
      xbinedges.erase(xbinedges.begin()+x);
      ybinedges.erase(ybinedges.begin()+y);
      ybinedges.erase(ybinedges.begin()+y+1);
      ynbins = ynbins -2;
      xnbins = xnbins -1;
    }
  }

  else if (x == xnbins - 1 && y == ynbins - 1){
    for (int j = y - 1; j < y + 1; j++){
      for(int i = x - 1; i < x + 1; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ) { z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }
    
    if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
    if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 3){
      xbinedges.erase(xbinedges.begin()+x);
      ybinedges.erase(ybinedges.begin()+y);
      ynbins = ynbins -1;
      xnbins = xnbins -1;
    }
  }

  else if (x > 0 && x < xnbins - 1 && y == ynbins - 1){
    for (int j = y - 1; j < y + 1; j++){
      for(int i = x - 1; i < x + 2; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ) { z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }
    
    if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
    if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}

    if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 5) {
      xbinedges.erase(xbinedges.begin()+x+1); 
      xbinedges.erase(xbinedges.begin()+x); 
      ybinedges.erase(ybinedges.begin()+y);
      ynbins = ynbins -1;
      xnbins = xnbins -2;
    }
  }

  else if (x == 0 && y == ynbins - 1){
    for (int j = y - 1; j < y + 1; j++){
      for(int i = x; i < x + 2; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ) { z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }
    
    if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}

    if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 3){
      xbinedges.erase(xbinedges.begin()+x+1);
      ybinedges.erase(ybinedges.begin()+y);
      ynbins = ynbins -1;
      xnbins = xnbins -1;
    }
  }

  else if (x == 0 && y > 0 && y < ynbins - 1){
    for (int j = y - 1; j < y + 2; j++){
      for(int i = x; i < x + 2; ++i){
        val = truemtrx[i][j];
        if (i == x && j == y) continue;
        if (val == 0 ) { z++; continue;}
        if (val >= currval) continue;
        if (val <= currval){
          currval = val;
          xhigher = i;
          yhigher = j;
        }
      }
    }
//    for (int i = 0; i <= ynbins; ++i){
//      std::cout << "que pedro ybinedges = " << ybinedges[i] << "\n";
//    }    
    if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}
    if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
    if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
    if ((xhigher == x && yhigher == y) || z == 5){
      xbinedges.erase(xbinedges.begin()+x+1);
      ybinedges.erase(ybinedges.begin()+y);
      ybinedges.erase(ybinedges.begin()+y+1);
      xnbins = xnbins -1;
      ynbins = ynbins -2;
    }
//    for (int i = 0; i <= ynbins; ++i){
//      std::cout << "que pex ybinedges = " << ybinedges[i] << "\n";
//    }
  }
  else std::cerr << "Algo esta mal carnal AAAAAAAAAAAAAAAAAAAh!";

//  std::cout << "" << xbinedges.size() << "\n";

  for (int i = 0; i <= xnbins; ++i)
    NewBinning.push_back(xbinedges[i]);

  for (int i = 0; i <= ynbins; ++i){
    NewBinning.push_back(ybinedges[i]);
//    std::cout << "yxbinedges = " << ybinedges[i] << "\n";
  }
  return NewBinning;

}
//==============================================================================
// Main
//==============================================================================
void runMigMtxBinning(int signal_definition_int = 0,
                      std::string plist = "ME1L", bool is_grid = false,
 		      std::string input_file = "", int run = 0) {

  // Macro Utility
  std::string mc_file_list;

  assert(!(is_grid && input_file.empty()) &&
         "On the grid, infile must be specified.");
  // const bool use_xrootd = false;
  mc_file_list = input_file.empty()
                     ? GetPlaylistFile(plist, true /*, use_xrootd*/)
                     : input_file;

  const std::string macro("runMigMtxBinning");
  const bool do_truth = false,
             do_systematics = false;
  CCPi::MacroUtil   util(signal_definition_int, mc_file_list, plist, do_truth,
                         is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  std::string xvar = "ptmu";
  std::string yvar = "tpi";
  std::string MigFileName ("Migration_AllVar_ALL_Weight");
//  std::string MigFileName (Form("Migration_%d",run));
  std::string TreeName = "Migration"; 
  std::string tag = "";

  // Variables and their histograms
  std::vector<Variable2D*> variables2D =
      GetAnalysisVariables2D(util.m_signal_definition, true); 
  // Starting the Loop
  CVUniverse* cvu = util.m_error_bands.at("cv").at(0);
  if (is_grid) Loop(util, cvu, kMC, variables2D, MigFileName, TreeName, xvar, yvar); 
  else {
  Variable2D* var = GetVar2D(variables2D, xvar, yvar);
  TArrayD oldedgesx = CCPi::GetBinning(xvar),
          oldedgesy = CCPi::GetBinning(yvar);
  int NBinsx = var->NBinsX(),
      NBinsy = var->NBinsY();

  double xmin = oldedgesx[0], 
//         xmax = oldedgesx[NBinsx],
         xmax = 5500.,
         ymin = oldedgesy[0],
         ymax = oldedgesy[NBinsy];
         ymin = 9.5;
         ymax = 20.;
//  xmin = 8250;

  const double xstep = (xmax - xmin)/10,
               ystep = (ymax - ymin)/10;

  std::vector<double> xedges, yedges;
  xedges.push_back(xmin);
  xedges.push_back(xstep);
  xedges.push_back(xmax);
  yedges.push_back(ymin);
  yedges.push_back(ystep);
  yedges.push_back(ymax);

  int xbin = 0,
      ybin = 0,
      dxcount = 0,
      dycount = 0;
      
  const double match = 0.5,
               minerr = 0.05;


  std::vector<int> test;
  for (int i = 0; i < 10; ++i)
    test.push_back(i*2);
  
  test.erase(test.end()-2);
  test.erase(test.begin()+3+1);
  for(int i = 0; i < test.size(); i++)
    std::cout << test[i] << " ";
  std::cout << "\n";


  // Taking the tree with the reco and true points
  TFile *f = new TFile(Form("%s.root", MigFileName.c_str()),"READ");
  TTree *T = (TTree*)f->Get(Form("%s%s",TreeName.c_str(),tag.c_str()));

  Double_t pnt[4];

  TBranch *MigVal = T->GetBranch(Form("%s_vs_%s", xvar.c_str(), yvar.c_str()));
  MigVal->SetAddress(&pnt);

  bool Rebinning = true;
  bool allxrange = false;
  bool allyrange = false;
  int count = 0;
  std::cout << "Entries = " << T->GetEntries() << "\n"; 
  Long64_t nentries = T->GetEntries();
  int counter =0;

  // ------------------------------------------------------
  // I'm developing a new algoritm to have a better binning
  // ------------------------------------------------------
  
  int xnbins = 8, ynbins = 4;
             
  const double xbinwidth = (xmax - xmin)/xnbins,
               ybinwidth = (ymax - ymin)/ynbins;

    std::vector<double> xbinedges;
    std::vector<double> ybinedges;
    xbinedges = {0., 250., 350., 400., 500, 600, 8000., 2000., 2500.};
    ybinedges = {35., 100., 140., 180., 350.};
    xbinedges.push_back(xmin);
    ybinedges.push_back(ymin);
/*    for (int i = 1; i <= xnbins;i++){
      xbinedges.push_back(xbinedges[i-1] + xbinwidth);
    }
    for (int j = 1; j <= ynbins;j++){
      ybinedges.push_back(ybinedges[j-1] + ybinwidth);
    }*/

    double currselval, binval, currselvalm, binvalm;
    int x, y, xm, ym, migFactor;
    Matrix truemtrx( xnbins,vector<double>(ynbins,0.));
    Matrix recomtrx(xnbins,vector<double>(ynbins,0.));
    double goodrecomtrx [xnbins][ynbins];
    double badRecoCounter [xnbins][ynbins];
    double badReco [xnbins][ynbins][2];
    double Dx, Dy, medx, medy; 

  while (Rebinning){
    
    for (int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; ++i){
        truemtrx[i][j] = 0.;
        recomtrx[i][j] = 0.;
        goodrecomtrx[i][j] = 0.; 
        badRecoCounter[i][j] = 0.;
      }
    }

    for (int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; ++i){
        badReco [i][j][0] = 0;
        badReco [i][j][1] = 0;
      }
    }

    for (Long64_t entry = 0; entry < nentries; ++entry){
      if (entry % 10000 == 0)
        std::cout << (entry / 10000) << "k" << std::endl;
      MigVal->GetEntry(entry);
      //Getting the reco and true values for each event
      double truex = pnt[0];
      double truey = pnt[2];
      double recox = pnt[1];
      double recoy = pnt[3];
       
      for (int j = 0; j < ynbins; j++){
        for(int i = 0; i < xnbins; ++i){
          if (truex > xbinedges[i] && truex < xbinedges[i+1] &&  
              truey > ybinedges[j] && truey < ybinedges[j+1]){
            truemtrx [i][j] = truemtrx [i][j] + 1.; 
            break;
          }
        }
      } 

      for (int j = 0; j < ynbins; j++){
        for(int i = 0; i < xnbins; ++i){
          if (recox > xbinedges[i] && recox < xbinedges[i+1] &&  
              recoy > ybinedges[j] && recoy < ybinedges[j+1]){
            recomtrx [i][j] = recomtrx [i][j] + 1.; 
            break;
          }
        }
      }
 
      for (int j = 0; j < ynbins; j++){
        for(int i = 0; i < xnbins; ++i){
          if (truex > xbinedges[i] && truex < xbinedges[i+1] &&
              truey > ybinedges[j] && truey < ybinedges[j+1]){
            medx = (xbinedges[i+1] - xbinedges[i])/2 + xbinedges[i];
            medy = (ybinedges[j+1] - ybinedges[j])/2 + ybinedges[j];
            if (recox > xbinedges[i] && recox < xbinedges[i+1] &&
                recoy > ybinedges[j] && recoy < ybinedges[j+1] ){
              goodrecomtrx [i][j] = goodrecomtrx [i][j] + 1; 
              break;
            }
            else {
              Dx = recox - medx;
              Dy = recoy - medy;
              badReco [i][j][0] = badReco [i][j][0] + Dx;
              badReco [i][j][1] = badReco [i][j][1] + Dy;
              badRecoCounter[i][j] = badRecoCounter[i][j] + 1.;
/*              if (j == 0) {
                std::cout << i << ", "<< j << "\n";
                std::cout << "medx = " << medx << " medy = "<< medy << "\n";
                std::cout << "truex = " << truex << " truey = "<< truey << "\n";
                std::cout << "recox = " << recox << " recoy = "<< recoy << "\n";
                std::cout << "Dx = " << Dx << " Dy = "<< Dy << "\n";                
              }
              if (recoy < 0) std::cout << "Negative ptmu >:/ \n";*/
/*            if (recox < xbinedges[i])
     	        badReco [i][j][0] = badReco[i][j][0] - 1.;
              if (recox > xbinedges[i+1])
                badReco [i][j][0] = badReco[i][j][0] + 1.;
              if (recoy < ybinedges[j])
                badReco [i][j][1] = badReco[i][j][1] - 1.;
              if (recoy > ybinedges[j+1])
                badReco [i][j][1] = badReco[i][j][1] + 1.;*/
              break;
            }
          }
        } 
      }
    }// End of loop for the the tree

    //found the smaller bin != to cero
    
    currselval = nentries;
    x = 0;
    y = 0;
    xm = 0;
    ym = 0;
    migFactor = 1;
    for (int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; ++i){
    //    migFactor = goodrecomtrx[i][j]/truemtrx[i][j];
        binval = truemtrx[i][j] * migFactor;
        if (binval == 0 ) continue;
        if (binval >= currselval) continue;
        if (binval <= currselval){  
          currselval = binval;
          x = i;
          y = j;
        }
      }
    } 
    for (int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; ++i){
        binvalm = goodrecomtrx[x][y]/truemtrx[i][j];
        if (binvalm == 0 ) continue;
        if (binvalm >= currselvalm) continue;
        if (binvalm <= currselvalm){  
          currselvalm = binvalm;
          xm = i;
          ym = j;
        }
      }
    } 
/*    
    for(int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; i++){
        std::cout << truemtrx[i][j] << "\t";
      }
      std::cout << "\n";
    }
*/
    bool isgoodreconstructed = false, isgoodstats = false;
    double diagonal = goodrecomtrx[xm][ym]/truemtrx[xm][ym];
    double StatError = 1/sqrt(currselval);
    
    std::cout << "Truth Matrix \n";
    for(int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; i++){
        std::cout << truemtrx[i][j] << "\t";
      }
      std::cout << "\n";
    }
    std::cout << "Reconstructed Matrix \n";
    for(int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; i++){
        std::cout << recomtrx[i][j] << "\t";
      }
      std::cout << "\n";
    }
    std::cout << "Migration Matrix \n";
    for(int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; i++){
        std::cout << goodrecomtrx[i][j]/truemtrx[i][j] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "Truth Error Matrix \n";
    for(int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; i++){
        std::cout << 1/sqrt(truemtrx[i][j]) << " ";
      }
      std::cout << "\n";
    }

    std::cout << "Bad Reco \n";
    for(int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; i++){
//        std::cout << "(" << badReco [i][j][0]/(truemtrx[i][j]-goodrecomtrx[i][j]) << ","<< badReco [i][j][1]/(truemtrx[i][j]-goodrecomtrx[i][j]) << ")" << " ";
        std::cout << "(" << badReco [i][j][0]/badRecoCounter[i][j] << ","<< badReco [i][j][1]/badRecoCounter[i][j] << ")" << " ";
      }
      std::cout << "\n";
    }
    std::cout << "N xbins = "<< xnbins << " N y bins = " << ynbins << "\n"; 
    std::cout << "The smaller bin is " << x << "," << y << "\n";
    std::cout << "The value is " << currselval << "\n";

//    for (int i = 0; i <= ynbins; ++i){
//      std::cout << "que pedro pinche pablo ybinedges = " << ybinedges[i] << "\n";
//    }

    std::cout << "X binning edges ";
    for (int i = 0; i < xnbins + 1; i++)
      std::cout << xbinedges[i] << "\t";
    std::cout << "\n Ybinning edges ";
    for (int i = 0; i < ynbins + 1; i++)
      std::cout << ybinedges[i] << "\t";
    std::cout << "\n";

    if (diagonal >= match) isgoodreconstructed = true;
    if (StatError <= minerr) isgoodstats = true; 
    if (isgoodreconstructed && isgoodstats )Rebinning = false;
    else if (!isgoodstats){ 
//      std::vector <double> newBins = byStaError(xbinedges, ybinedges, xnbins, ynbins, truemtrx,
//                                                x,y, (double)nentries);
//      if (xnbins < xbinedges.size()-1) xbinedges.erase(xbinedges.begin());
//      if (ynbins < ybinedges.size()-1) ybinedges.erase(ybinedges.begin());
     
//      std::cout << "NewBins size = " << newBins.size() << "\n";
//      for (int i = 0; i < newBins.size(); ++i)
//        std::cout << "NewBins = " << newBins[i] << "\n";
/*     
//      std::cout << "xnbins = " << xnbins << "\n";
      for (int i = 0; i <= xnbins ; ++i){
        xbinedges[i] = newBins [i];
//        std::cout << "xbinedges = " << xbinedges[i] << "\n";
      }
      
//      std::cout << "xnbins = " << xnbins << "\n";
      for (int i = 0; i <= ynbins; ++i){ 
        ybinedges[i] = newBins [i+xnbins+1];
//        std::cout << "yxbinedges = " << ybinedges[i] << "\n";
      }*/

     double currval = nentries, val;
      int xhigher, yhigher; 
      if (x > 0 && y > 0 && x < xnbins - 1 && y < ynbins - 1){
        for (int j = y - 1; j < y + 2; j++){
          for(int i = x - 1; i < x + 2; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
 
        if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
	if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}
        if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
	if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y){
 	  xbinedges.erase(xbinedges.begin()+x+1); 
	  xbinedges.erase(xbinedges.begin()+x); 
	  ybinedges.erase(ybinedges.begin()+y);
	  ybinedges.erase(ybinedges.begin()+y+1);
	  ynbins = ynbins -2;
	  xnbins = xnbins -2;
        } 
      }

      else if (x == 0 && y == 0){
        for (int j = y; j < y + 2; j++){
          for(int i = x; i < x + 2; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
        
        if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}

        if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y){
          xbinedges.erase(xbinedges.begin()+x+1); 
          ybinedges.erase(ybinedges.begin()+y+1);
          ynbins = ynbins -1;
          xnbins = xnbins -1;
        } 
      }

      else if (x > 0 && y == 0 && x < xnbins - 1){
        for (int j = y; j < y + 2; j++){
          for(int i = x - 1; i < x + 2; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
        
        if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
        if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}
        if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y){
          xbinedges.erase(xbinedges.begin()+x+1); 
          xbinedges.erase(xbinedges.begin()+x); 
          ybinedges.erase(ybinedges.begin()+y+1);
          ynbins = ynbins -1;
          xnbins = xnbins -2;
        }
      } 

      else if (y == 0 && x == xnbins - 1){
        for (int j = y; j < y + 2; j++){
          for(int i = x - 1; i < x + 1; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
        
        if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
        if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y){
          xbinedges.erase(xbinedges.begin()+x);
          ybinedges.erase(ybinedges.begin()+y+1);
          ynbins = ynbins -1;
          xnbins = xnbins -1;
        }
      }

      else if ( y > 0 && x == xnbins - 1 && y < ynbins - 1){
        for (int j = y - 1; j < y + 2; j++){
          for(int i = x - 1; i < x + 1; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
        std::cout << "Que pex Xh = " << xhigher << " Yh = " << yhigher << " X = " << x << " Y = " << y <<"\n";
        if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
        if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
        if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y){
          xbinedges.erase(xbinedges.begin()+x);
          ybinedges.erase(ybinedges.begin()+y);
          ybinedges.erase(ybinedges.begin()+y+1);
          ynbins = ynbins -2;
          xnbins = xnbins -1;
        }
      }

      else if (x == xnbins - 1 && y == ynbins - 1){
        for (int j = y - 1; j < y + 1; j++){
          for(int i = x - 1; i < x + 1; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
        
        if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
        if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y){
          xbinedges.erase(xbinedges.begin()+x);
          ybinedges.erase(ybinedges.begin()+y);
          ynbins = ynbins -1;
          xnbins = xnbins -1;
        }
      }

      else if (x > 0 && x < xnbins - 1 && y == ynbins - 1){
        for (int j = y - 1; j < y + 1; j++){
          for(int i = x - 1; i < x + 2; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
        
        if (xhigher < x){ xbinedges.erase(xbinedges.begin()+x); xnbins = xnbins -1;}
        if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}

        if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y) {
          xbinedges.erase(xbinedges.begin()+x+1); 
          xbinedges.erase(xbinedges.begin()+x); 
          ybinedges.erase(ybinedges.begin()+y);
          ynbins = ynbins -1;
          xnbins = xnbins -2;
        }
      }

      else if (x == 0 && y == ynbins - 1){
        for (int j = y - 1; j < y + 1; j++){
          for(int i = x; i < x + 2; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
        
        if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}

        if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y){
          xbinedges.erase(xbinedges.begin()+x+1);
          ybinedges.erase(ybinedges.begin()+y);
          ynbins = ynbins -1;
          xnbins = xnbins -1;
        }
      }

      else if (x == 0 && y > 0 && y < ynbins - 1){
        for (int j = y - 1; j < y + 2; j++){
          for(int i = x; i < x + 2; ++i){
            val = truemtrx[i][j];
            if (i == x && j == y) continue;
            if (val == 0 ) continue;
            if (val >= currval) continue;
            if (val <= currval){
              currval = val;
              xhigher = i;
              yhigher = j;
            }
          }
        }
        
        if (xhigher > x){ xbinedges.erase(xbinedges.begin()+x+1); xnbins = xnbins -1;}

        if (yhigher < y){ ybinedges.erase(ybinedges.begin()+y); ynbins = ynbins -1;}
        if (yhigher > y){ ybinedges.erase(ybinedges.begin()+y+1); ynbins = ynbins -1;}
        if (xhigher == x && yhigher == y){
          xbinedges.erase(xbinedges.begin()+x+1);
          ybinedges.erase(ybinedges.begin()+y);
          ybinedges.erase(ybinedges.begin()+y+1);
          xnbins = xnbins -1;
          ynbins = ynbins -2;
        }
      }
      else std::cerr << "Algo esta mal carnal AAAAAAAAAAAAAAAAAAAh!";
    }
    else Rebinning = false;


    Rebinning = false;
  }


    for(int j = 0; j < ynbins; j++){
      for(int i = 0; i < xnbins; i++){
        std::cout << truemtrx[i][j] << "\t";
      }
      std::cout << "\n";
    }


    std::cout << "The smaller bin is " << x << "," << y << "\n";
    std::cout << "The value is " << currselval << "\n";

    std::cout << "X binning edges ";
    for (int i = 0; i < xnbins + 1; i++)
      std::cout << xbinedges[i] << "\t";
    std::cout << "\n Ybinning edges ";
    for (int i = 0; i < ynbins + 1; i++)
      std::cout << ybinedges[i] << "\t";
    std::cout << "\n";
  }//End of !is_grid

//  Starting the Loop to satisfy the requerements 
/*  while (Rebinning){
    int xrecoutL = 0,
        xrecoutR = 0,
        yrecoutD = 0,
        yrecoutU = 0,
        goodReco = 0; 

    dxcount = 0;
    dycount = 0;
    count = 0;

 //   std::vector<int> passedEvents;
//    passedEvents.push_back(-1);
    // Statiting the Loop for all the events

    for (Long64_t entry = 0; entry < nentries; ++entry){
      if (entry % 10000 == 0)
        std::cout << (entry / 10000) << "k" << std::endl;
//      if (passed) continue;
      MigVal->GetEntry(entry);
      //Getting the reco and true values for each event
      double truex = pnt[0];
      double truey = pnt[2];
      double recox = pnt[1];
      double recoy = pnt[3];

      if (truex < xedges[xbin] || truex > xedges[xbin + 1] ||
          truey < yedges[ybin] || truey > yedges[ybin + 1]){
          if (truex > (xedges[xbin + 1] + xstep) || truey > (yedges[ybin + 1] + ystep)){
            continue;}
          else if (truex > xedges[xbin + 1] && truex < (xedges[xbin + 1] + xstep) &&
                   truey > yedges[ybin] && truey < (yedges[ybin + 1] + ystep)){
            dxcount = dxcount + 1;
          }
          else if (truex > xedges[xbin] && truex < (xedges[xbin + 1] + xstep) &&
                   truey > yedges[ybin + 1] && truey < (yedges[ybin + 1] + ystep)){
            dycount = dycount + 1;
          }        
      }
      
//      passedEvents.push_back(entry);
      count++;
      if (recox > xedges[xbin] && recox < xedges[xbin + 1] &&
          recoy > yedges[ybin] && recoy < yedges[ybin + 1]) goodReco++;
      else if (recox < xedges[xbin]) xrecoutL++;
      else if (recox > xedges[xbin + 1]) xrecoutR++;
      else if (recoy < yedges[ybin]) yrecoutD++;
      else if (recoy > yedges[ybin + 1]) yrecoutU++;
      else  std::cerr << "Algo anda mal carnal\n";
      
    }

    double currError = 1/sqrt((double)count);
    double diagonal = (double)goodReco/(double)count;
    bool metmtxcriteria = true;
    bool meterrcriteria = true;
    // we obtain if the migration matrix is goin to Left-Right, Down-Up
    double xLR = (double)xrecoutR/(double)xrecoutL,
           yDU = (double)yrecoutU/(double)yrecoutD;
    if (currError > minerr) meterrcriteria = false;
    if (diagonal < match) metmtxcriteria = false;  

//  std::cout << "Count = " << count << " goodReco = " << goodReco << "\n"; 

    if (metmtxcriteria && meterrcriteria){
      if (xedges[xbin + 1] == xmax && yedges[ybin + 1] == ymax){
        Rebinning = false;
      }
      else if (xedges[xbin + 1] < xmax){
        if (ybin == 0) {
          xbin = xbin + 1;
	  xedges.push_back(xedges[xbin + 1]);
          xedges[xbin + 1] = xedges[xbin] + xstep;
 	}
        else {
          xbin = xbin + 1;
	  allxrange = false;
	}
      }
      else if (xedges[xbin + 1] == xmax){
	xbin = 0;
        ybin = ybin + 1;
      }
      else std::cerr << "Algo raro pasa carnal, se cumplen los dos criterios";
    }
    else {
      if (allxrange && allyrange){
        std::cout << "Ultimo bin no optimizado :(";
        Rebinning = false;
      }
      else if (allxrange){
	yedges[ybin + 1] = yedges[ybin + 1] + ystep;
      }
      else if (allyrange){
	xedges[xbin + 1] = xedges[xbin + 1] + xstep;
      }
      else if ((double)dxcount*xLR > (double)dycount*yDU){
        if (xmax <= xedges[xbin + 1] + xstep){
          xedges.erase(xedges.end()-2);
          allxrange = true; 
        }
        else {
	  xedges[xbin + 1] = xedges[xbin + 1] + xstep;         
        }
      }
      else if ((double)dxcount*xLR < (double)dycount*yDU){
        if (ymax <= yedges[xbin + 1] + ystep){
          yedges.erase(yedges.end()-2);
          allyrange = true;
        }
        else {
          yedges[ybin + 1] = yedges[ybin + 1] + ystep;
        }
      }
      else {
        if (ymax < yedges[xbin + 1] + ystep){
          yedges.erase(yedges.end()-2);
          allyrange = true;
        }
        else if (xmax < xedges[xbin + 1] + xstep){
          xedges.erase(xedges.end()-2);
          allxrange = true;
        }
        else {
          xedges[xbin + 1] = xedges[xbin + 1] + xstep;
          yedges[ybin + 1] = yedges[ybin + 1] + ystep;
        }  
      }

    if(!allxrange || !allyrange){
        if (xmax <= xedges[xbin + 1] + xstep){
          xedges.erase(xedges.end()-2);
          allxrange = true;
        }
        if (ymax < yedges[xbin + 1] + ystep){
          yedges.erase(yedges.end()-2);
          allyrange = true;
        }
      }

    }// end it doesn't met a creteria

    if (counter == 3) Rebinning = false;
    counter ++;

  } // end of while 
*/

  std::cout << "\n It's done \n";
  //============================================================================
}

#endif  // runMigMtxBinning_C


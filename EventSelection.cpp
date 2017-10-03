#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include <TRandom1.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TChain.h>
#include "TThread.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>

#include "AnaTree.h"
#include "Spectrum.hpp"
#include "Spectrum2D.hpp"
#include "PlottingTools.h"

//#include "PlotHandler.hpp"
//#include "SelectionTools.hpp"

const bool _breakdownPlots = true;

double _beamSpillStarts = 3.2;  // us
double _beamSpillEnds   = 4.8;  // us
double _flashShift      = 0.;//4.06; //us

const double targetPOT = 4.95e19;

using namespace std;

//____________________________________________________________________________________________________
void DrawProgressBar(double progress, double barWidth) {
  
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}

//____________________________________________________________________________________________________
void DrawPOT(double pot)
{
  //std::string str = "Simulated POT:" + std::to_string(pot);
  
  std::stringstream sstm;
  sstm << "Simulated POT: " << pot;
  std::string str = sstm.str();
  
  TLatex* pot_latex = new TLatex(.10, .96, str.c_str());
  pot_latex->SetTextColor(kGray+2);
  pot_latex->SetNDC();
  pot_latex->SetTextSize(1/30.);
  pot_latex->SetTextAlign(10); //left adjusted
  pot_latex->Draw();
  
  
  std::stringstream sstm2;
  sstm2 << "Scaled to POT: " << targetPOT;
  str = sstm2.str();
  
  TLatex* pot_latex_2 = new TLatex(.10, .92, str.c_str());
  pot_latex_2->SetTextFont(62);
  pot_latex_2->SetTextColor(kGray+2);
  pot_latex_2->SetNDC();
  pot_latex_2->SetTextSize(1/30.);
  pot_latex_2->SetTextAlign(10);//left adjusted
  pot_latex_2->Draw();
  
}

//____________________________________________________________________________________________________
void ActivateBranches(AnaTree *at) {
  
  //t->fChain->SetBranchStatus("",1);
}

//____________________________________________________________________________________________________
double CalcLength(const double& x_1, const double& y_1, const double& z_1, const double& x_2, const double& y_2, const double& z_2) {
  return sqrt(pow(x_1-x_2, 2) + pow(y_1-y_2, 2) + pow(z_1-z_2, 2));
}




//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
int main(int argc, char* argv[]) {
  
  clock_t begin = clock();
  
  std::string filen     = "ubxsec_output.root";
  bool evalPOT = false;
  int maxEntries = -1;
  
  //*************************
  //* Getting input parameters
  //*************************
  
  int c;
  //int digit_optind = 0;
  
  while (1) {
    //int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"filename",    required_argument, 0,  'f' },
      {"maxentries",  required_argument, 0,  'e' },
      {"flashstart",  required_argument, 0,  'a' },
      {"flashend",    required_argument, 0,  'b' },
      {"flashshift",  required_argument, 0,  's' },
      {0,             0,                 0,   0  }
    };
    
    c = getopt_long(argc, argv, "f:e:a:b:s:p",
                    long_options, &option_index);
    if (c == -1)
      break;
    
    switch (c) {
      case 0:
        std::cout << "Option " << long_options[option_index].name;
        if (optarg)
          std::cout << " with arg " << optarg << std::endl;
        break;
        
      case 'f':
        std::cout << "Input file: " << optarg << std::endl;
        filen = optarg;
        break;
        
      case 'e':
        std::cout << "Run over entries: " << optarg << std::endl;
        maxEntries = atof(optarg);
        break;
        
      case 'a':
        std::cout << "Beam spill start time: " << optarg << std::endl;
        _beamSpillStarts = atof(optarg);
        break;
        
      case 'b':
        std::cout << "Beam spill end time: " << optarg << std::endl;
        _beamSpillEnds = atof(optarg);
        break;
        
      case 's':
        std::cout << "Flash shift: " << optarg << std::endl;
        _flashShift = atof(optarg);
        break;

      case 'p':
        std::cout << "Calculating POT." << std::endl;
        evalPOT = true;
        break;
        
      case '?':
        break;
        
      default:
        printf("?? getopt returned character code 0%o ??\n", c);
    }
  }
  
  if (optind < argc) {
    printf("non-option ARGV-elements: ");
    while (optind < argc)
      printf("%s ", argv[optind++]);
    printf("\n");
  }
  
  std::cout << "_beamSpillStarts is " << _beamSpillStarts << std::endl;
  std::cout << "_beamSpillEnds is   " << _beamSpillEnds << std::endl;

  
  

  //*************************
  //* Starting ROOT application
  //*************************
  
  TApplication* rootapp = new TApplication("ROOT Application",&argc, argv);
  gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS
  gROOT->ProcessLine(".x rootlogon.C");

  gROOT->ProcessLine(".L loader_C.so");
  
  std::cout << "Opening output file ubxsecana_output.root." << std::endl;
  TFile *file_out = new TFile("ubxsecana_output.root","RECREATE");
  if ( file_out->IsOpen() )
    std::cout << "File opened successfully" << std::endl;
  
  string pattern = filen;
  
  
  
  //*************************
  //* Getting POTs
  //*************************
  
  double totalPOT = 0.;
  
  if (maxEntries > 0) evalPOT = false;
  
  if (evalPOT) {
    
    cout << " ----- " << endl;
    cout << "| Calculating POT" << endl;
    cout << "| " << endl;
    
    TChain *cpot;
    cpot = new TChain("UBXSec/pottree");
    cpot->Add(pattern.c_str());
    cout << "| Number of entries in the pot tree: " << cpot->GetEntries() << endl;
    Double_t pot;
    cpot->SetBranchAddress("pot", &pot);
    for (int potEntry = 0; potEntry < cpot->GetEntries(); potEntry++) {
      cpot->GetEntry(potEntry);
      totalPOT += pot;
    } // end loop entries
    cout << "| Total POT: " << totalPOT << endl;
    cout << " ----- " << endl << endl;
  } // end if evalPOT
  else
    totalPOT = -1.;
  
  double pot_scaling = 1.;
  if (evalPOT) pot_scaling = targetPOT/totalPOT;
  
  
  TChain *chain_ubxsec;
  chain_ubxsec = new TChain("UBXSec/tree");
  chain_ubxsec->Add(pattern.c_str());
  
  cout << "Using file: " << pattern << endl;
  
  int Nfiles = chain_ubxsec->GetNtrees();
  cout << "Number of files: " << Nfiles << endl;
  
  int evts = chain_ubxsec -> GetEntries();
  cout << "Number of events used is: " << evts << endl;
  
  AnaTree * t = new AnaTree(chain_ubxsec);
  ActivateBranches(t);
  
  
  /*
   //Create a new file + a clone of old tree in new file
   TFile *oldfile = new TFile(pattern.c_str());
   TTree *oldtree = (TTree*)oldfile->Get("analysistree/at");
   TFile *newfile = new TFile("selectedEntries.root","recreate");
   TTree *newtree = oldtree->CloneTree(0);
   */
  
  //Spectrum* Sflashtime      = new Spectrum("flash_time",      ";Flash Time [#mus];Entries per bin",       300000, -3000, 3000, totalPOT);
  
  int nsignal = 0;
  
  int signal_sel = 0;
  int bkg_anumu_sel = 0;
  int bkg_nue_sel = 0;
  int bkg_nc_sel = 0;
  int bkg_outfv_sel = 0;
  int bkg_cosmic_sel = 0;
  int bkg_cosmic_top_sel = 0;
  
  int nEvtsWFlashInBeamSpill = 0;
  int nNumuCC = 0;
  
  int nue_cc_fv = 0;
  int nue_cc_selected = 0;
  int nue_cc_selected_total = 0;

  int nSignalWMuonReco = 0;
  int nSignalMuonRecoVtxOk = 0;
  
  int nSignalFlashMatched = 0;
  
  int n_slc_nu_origin = 0;
  int n_slc_acpt_tag_nu = 0;
  
  TH1D* h_eff_num = new TH1D("h_eff_num", "h_eff_num", 6, 0, 4);
  TH1D* h_eff_den = new TH1D("h_eff_den", "h_eff_den", 6, 0, 4);
  TEfficiency* pEff = new TEfficiency("eff",";Neutrino Energy (truth) [GeV];Efficiency",6, 0, 4);
  
  
  
  TH1D* h_chi2 = new TH1D("h_chi2", "h_chi2", 50, 0, 50);
  TH1D* h_flsTime = new TH1D("h_flsTime", ";Flash time w.r.t. trigger [#mus];Flashes", 175, 0, 25);
  TH1D* h_flsTime_wcut = new TH1D("h_flsTime_wcut", ";Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 175, 0, 25);
  h_flsTime->Sumw2(); h_flsTime_wcut->Sumw2();
  TH1D* h_nslices = new TH1D("h_nslices", ";Number of slices per event;Entries per bin", 15, 0, 15);
  TH1D* h_vtx_resolution = new TH1D("h_nslh_vtx_resolutionices", ";Vertex resolution (2D) [cm];Entries per bin", 300, 0, 500);
  
  TH2D* h_frac_diff = new TH2D("h_frac_diff", ";PMT ID; Fractional difference", 32, 0, 32, 80, -2, 2);
  TH2D* h_frac_diff_others = new TH2D("h_frac_diff", ";PMT ID; Fractional difference", 32, 0, 32, 80, -2, 2);
  double hypo_spec_x[32], hypo_spec_y[32];
  double meas_spec_x[32], meas_spec_y[32];
  double numc_spec_x[32], numc_spec_y[32];
  
  TH1D* h_xdiff = new TH1D("h_xdiff", "h_xdiff", 1000, -100,100);
  TH1D* h_xdiff_others = new TH1D("h_xdiff_others", "h_xdiff_others", 1000, -100,100);
  TH1D* h_zdiff = new TH1D("h_zdiff", "h_zdiff", 1000, 0,1000);
  TH1D* h_zdiff_others = new TH1D("h_zdiff_others", "h_zdiff_others", 1000, 0,1000);
  std::map<std::string,TH1D*> hmap_xdiff;
  hmap_xdiff["total"] = new TH1D("h_xdiff_total", ";QLL x - TPC x [cm];", 80, -200,200);
  hmap_xdiff["signal"] = new TH1D("h_xdiff_signal", ";QLL x - TPC x [cm];", 80, -200,200);
  hmap_xdiff["background"] = new TH1D("h_xdiff_background", ";QLL x - TPC x [cm];", 80, -200,200);
  std::map<std::string,TH1D*> hmap_zdiff;
  hmap_zdiff["total"] = new TH1D("h_zdiff_total", ";Hypo z - Flash z [cm];", 160, -400,400);
  hmap_zdiff["signal"] = new TH1D("h_zdiff_signal", ";Hypo z - Flash z [cm];", 160, -400,400);
  hmap_zdiff["background"] = new TH1D("h_zdiff_background", ";Hypo z - Flash z [cm];", 160, -400,400);

  
  TH1D* h_vtxcheck_angle_good = new TH1D("h_vtxcheck_angle_good", ";Angle [rad];Entries per bin", 100, 0, 4);
  TH1D* h_vtxcheck_angle_bad  = new TH1D("h_vtxcheck_angle_bad",  ";Angle [rad];Entries per bin",  100, 0, 4);
  
  TH1D* h_muon_track_eff  = new TH1D("h_muon_track_eff",  ";Muon track efficiency;Entries per bin",  100, 0, 1);
  TH1D* h_muon_track_pur  = new TH1D("h_muon_track_pur",  ";Muon track purity;Entries per bin",  100, 0, 1);
  
  TH1D* h_mueff_num = new TH1D("h_mueff_num", "h_mueff_num", 30, 0, 2);
  TH1D* h_mueff_2_num = new TH1D("h_mueff_2_num", "h_mueff_2_num", 30, 0, 2);
  TH1D* h_mueff_den = new TH1D("h_mueff_den", "h_mueff_den", 30, 0, 2);
  
  TH2D* h_mu_eff_mom = new TH2D("h_mu_eff_mom", ";True Muon Momentum [GeV]; Efficiency", 50, 0, 2, 20, 0, 1);
  TH2D* h_mu_pur_mom = new TH2D("h_mu_pur_mom", ";True Muon Momentum [GeV]; Purity", 50, 0, 2, 20, 0, 1);
  
  TH2D* h_mumom_nue = new TH2D("h_mumom_nue", ";True Neutrino Energy [GeV]; True Muon Momentum [GeV]", 50, 0, 2, 50, 0, 4);
  
  TH1D* h_acpt_tagged  = new TH1D("h_acpt_tagged",  ";Tagged TPC Objects;Entries per bin",  10, 0, 10);
  
  TH1D* h_slice_origin = new TH1D("h_slice_origin",  ";;",  3, -0.5, 2.5);
  
  TH1D* h_slice_npfp = new TH1D("h_slice_npfp",  ";npfp;",  10, 0, 10);
  TH1D* h_slice_npfp_others = new TH1D("h_slice_npfp_others",  ";npfp;",  10, 0, 10);
  
  TH1D* h_slice_ntrack = new TH1D("h_slice_ntrack",  ";npfp;",  10, 0, 10);
  TH1D* h_slice_ntrack_others = new TH1D("h_slice_ntrack_others",  ";npfp;",  10, 0, 10);
  
  TH1D* h_fm_score = new TH1D("h_fm_score",  ";fm score;",  500, 0, 10);
  TH1D* h_fm_score_others = new TH1D("h_fm_score_other",  ";fm score;",  500, 0, 10);
  TH2D* h_fm_score_pe = new TH2D("h_fm_score_pe",  ";fm score;Reco PE",  500, 0, 10, 500, 0, 2000);
  
  TH1D* h_n_slc_flsmatch = new TH1D("h_n_slc_flsmatch",  ";n slices flash matched per event;",  10, 0, 10);
  
  std::map<std::string,TH1D*> hmap_trklen;
  hmap_trklen["total"] = new TH1D("h_trklen_total", "; Track length;", 30, 0, 700);
  hmap_trklen["signal"] = new TH1D("h_trklen_signal", "; Track length;", 30, 0, 700);
  hmap_trklen["cosmic"] = new TH1D("h_trklen_cosmic", "; Track length;", 30, 0, 700);
  hmap_trklen["cosmic_stopmu"] = new TH1D("h_trklen_cosmic_stopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["cosmic_nostopmu"] = new TH1D("h_trklen_cosmic_nostopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["outfv"] = new TH1D("h_trklen_outfv", "; Track length;", 30, 0, 700);
  hmap_trklen["outfv_stopmu"] = new TH1D("h_trklen_outfv_stopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["outfv_nostopmu"] = new TH1D("h_trklen_outfv_nostopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["nc"] = new TH1D("h_trklen_nc", "; Track length;", 30, 0, 700);
  hmap_trklen["nc_proton"] = new TH1D("h_trklen_nc_proton", "; Track length;", 30, 0, 700);
  hmap_trklen["nc_pion"] = new TH1D("h_trklen_nc_pion", "; Track length;", 30, 0, 700);
  hmap_trklen["nc_other"] = new TH1D("h_trklen_nc_other", "; Track length;", 30, 0, 700);
  hmap_trklen["anumu"] = new TH1D("h_trklen_anumu", "; Track length;", 30, 0, 700);
  hmap_trklen["nue"] = new TH1D("h_trklen_nue", "; Track length;", 30, 0, 700);
  hmap_trklen["signal_stopmu"] = new TH1D("h_trklen_signal_stopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["signal_nostopmu"] = new TH1D("h_trklen_signal_nostopmu", "; Track length;", 30, 0, 700);

  std::map<std::string,TH1D*> hmap_trkphi;
  hmap_trkphi["total"] = new TH1D("h_trkphi_total", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["signal"] = new TH1D("h_trkphi_signal", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["cosmic"] = new TH1D("h_trkphi_cosmic", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["outfv"] = new TH1D("h_trkphi_outfv", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nc"] = new TH1D("h_trkphi_nc", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["anumu"] = new TH1D("h_trkphi_anumu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nue"] = new TH1D("h_trkphi_nue", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["cosmic_stopmu"] = new TH1D("h_trkphi_cosmic_stopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["cosmic_nostopmu"] = new TH1D("h_trkphi_cosmic_nostopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["outfv_stopmu"] = new TH1D("h_trkphi_outfv_stopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["outfv_nostopmu"] = new TH1D("h_trkphi_outfv_nostopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nc_proton"] = new TH1D("h_trkphi_nc_proton", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nc_pion"] = new TH1D("h_trkphi_nc_pion", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nc_other"] = new TH1D("h_trkphi_nc_other", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["signal_stopmu"] = new TH1D("h_trkphi_signal_stopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["signal_nostopmu"] = new TH1D("h_trkphi_signal_nostopmu", "; Track #phi;", 20, -3.15, 3.15);
  
  std::map<std::string,TH1D*> hmap_trktheta;
  hmap_trktheta["total"] = new TH1D("h_trktheta_total", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["signal"] = new TH1D("h_trktheta_signal", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["cosmic"] = new TH1D("h_trktheta_cosmic", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["outfv"] = new TH1D("h_trktheta_outfv", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["nc"] = new TH1D("h_trktheta_nc", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["anumu"] = new TH1D("h_trktheta_anumu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["nue"] = new TH1D("h_trktheta_nue", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["cosmic_stopmu"] = new TH1D("h_trktheta_cosmic_stopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["cosmic_nostopmu"] = new TH1D("h_trktheta_cosmic_nostopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["outfv_stopmu"] = new TH1D("h_trktheta_outfv_stopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["outfv_nostopmu"] = new TH1D("h_trktheta_outfv_nostopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["nc_proton"] = new TH1D("h_trktheta_nc_proton", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["nc_pion"] = new TH1D("h_trktheta_nc_pion", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["nc_other"] = new TH1D("h_trktheta_nc_other", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["signal_stopmu"] = new TH1D("h_trktheta_signal_stopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta["signal_nostopmu"] = new TH1D("h_trktheta_signal_nostopmu", "; Track cos(#theta);", 30, -1, 1);

  std::map<std::string,TH1D*> hmap_multpfp;
  hmap_multpfp["total"] = new TH1D("h_multpfp_total", "; PFP Multiplicity", 10, 0, 10);
  hmap_multpfp["signal"] = new TH1D("h_multpfp_signal", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["cosmic"] = new TH1D("h_multpfp_cosmic", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["outfv"] = new TH1D("h_multpfp_outfv", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nc"] = new TH1D("h_multpfp_nc", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["anumu"] = new TH1D("h_multpfp_anumu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nue"] = new TH1D("h_multpfp_nue", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["cosmic_stopmu"] = new TH1D("h_multpfp_cosmic_stopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["cosmic_nostopmu"] = new TH1D("h_multpfp_cosmic_nostopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["outfv_stopmu"] = new TH1D("h_multpfp_outfv_stopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["outfv_nostopmu"] = new TH1D("h_multpfp_outfv_nostopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nc_proton"] = new TH1D("h_multpfp_nc_proton", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nc_pion"] = new TH1D("h_multpfp_nc_pion", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nc_other"] = new TH1D("h_multpfp_nc_other", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["signal_stopmu"] = new TH1D("h_multpfp_signal_stopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["signal_nostopmu"] = new TH1D("h_multpfp_signal_nostopmu", "; PFP Multiplicity;", 10, 0, 10);

  std::map<std::string,TH1D*> hmap_multtracktol;
  hmap_multtracktol["total"] = new TH1D("h_multtracktol_total", "; Track Multiplicity (5 cm)", 10, 0, 10);
  hmap_multtracktol["signal"] = new TH1D("h_multtracktol_signal", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["cosmic"] = new TH1D("h_multtracktol_cosmic", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["outfv"] = new TH1D("h_multtracktol_outfv", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nc"] = new TH1D("h_multtracktol_nc", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["anumu"] = new TH1D("h_multtracktol_anumu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nue"] = new TH1D("h_multtracktol_nue", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["cosmic_stopmu"] = new TH1D("h_multtracktol_cosmic_stopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["cosmic_nostopmu"] = new TH1D("h_multtracktol_cosmic_nostopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["outfv_stopmu"] = new TH1D("h_multtracktol_outfv_stopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["outfv_nostopmu"] = new TH1D("h_multtracktol_outfv_nostopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nc_proton"] = new TH1D("h_multtracktol_nc_proton", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nc_pion"] = new TH1D("h_multtracktol_nc_pion", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nc_other"] = new TH1D("h_multtracktol_nc_other", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["signal_stopmu"] = new TH1D("h_multtracktol_signal_stopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["signal_nostopmu"] = new TH1D("h_multtracktol_signal_nostopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);

  TH1D* h_pot = new TH1D("h_pot", "First bin contains number of POT (not valid on data)", 1, 0, 1);
  TH1D* h_nevts = new TH1D("h_nevts", "First bin contains number of events", 1, 0, 1);

  h_pot->SetBinContent(1, totalPOT);
  h_nevts->SetBinContent(1, evts);

  
  int barWidth = 70;
  
  if(maxEntries > 0.) evts = maxEntries;
  
  for(int i = 0; i < evts; i++) {
    
    DrawProgressBar((double)i/(double)evts, barWidth);
    
    chain_ubxsec->GetEntry(i);
    
    //SelectionTools * selection = new SelectionTools(t);
    
    //cout << "***** Event " << i << endl;
    
    // ************************
    //
    // Preliminary distributions
    //
    // ************************
    
    /*
     // Flashes
     for (int fls = 0; fls < t->no_flashes; fls++) {
     Sflashtime      ->Fill(t->flash_time[fls]);
     Sflashpe        ->Fill(t->flash_pe[fls]);
     Sflashycenter   ->Fill(t->flash_ycenter[fls]);
     Sflashzcenter   ->Fill(t->flash_zcenter[fls]);
     Sflashtimewidth ->Fill(t->flash_timewidth[fls]);
     if (t->flash_pe[fls] > 50) {
     Sflashtime50pe->Fill(t->flash_time[fls]);
     int k = 0;
     double distance = sqrt(pow(54.999*100-t->vx_flux[k],2)+pow(74.461*100-t->vy_flux[k],2)+pow(677.611*100-t->vz_flux[k],2));
     Sfls_timeVSnu_distance->Fill(t->flash_time[fls],distance);
     }
     Sfls_timeVSpe->Fill(t->flash_time[fls], t->flash_pe[fls]);
     }
     bool doneForThisEvent = false;
     bool doneForThisEvent_cosmic = false;
     for (int geantpar = 0; geantpar < t->geant_list_size; geantpar++) {
     if (!doneForThisEvent && t->origin[geantpar] == 1 && t->process_primary[geantpar]==1) {
     Sgeanttruetime_neutrino->Fill(t->StartT[geantpar]);
     doneForThisEvent = true;
     }
     if (!doneForThisEvent_cosmic &&t->origin[geantpar] == 2 && t->process_primary[geantpar]==1) {
     Sgeanttruetime_cosmic->Fill(t->StartT[geantpar]);
     doneForThisEvent_cosmic = true;
     }
     }
     */
    
    
    
    if (t->nbeamfls == 0) continue;
    int flashInBeamSpill = -1;
    
    bool goodflash = false;
    for (int fls = 0; fls < t->nbeamfls; fls ++){
      h_flsTime->Fill(t->beamfls_time->at(fls) - _flashShift);
      if(t->beamfls_pe->at(fls) > 50.) {
        h_flsTime_wcut->Fill(t->beamfls_time->at(fls) - _flashShift);
      }
      if (t->beamfls_time->at(fls) > _beamSpillStarts && t->beamfls_time->at(fls) < _beamSpillEnds) {
        
        flashInBeamSpill = fls;
        if (t->beamfls_pe->at(fls) >= 50) {
          goodflash = true;
          nEvtsWFlashInBeamSpill++;
        }
      }
    }
    
    if (flashInBeamSpill == -1) continue;
    
    
    
    
    
    
    bool isSignal = false;
    if (t->nupdg == 14 && t->ccnc == 0 && t->fv == 1){
      nsignal++;
      isSignal = true;
      h_eff_den->Fill(t->nu_e);
      h_mueff_den->Fill(t->true_muon_mom);
      
      if (t->muon_is_reco){
        h_mumom_nue->Fill(t->nu_e, t->true_muon_mom);
        nSignalWMuonReco++;
        h_mueff_num->Fill(t->true_muon_mom);
        for (auto origin : *t->slc_origin){
          if (origin == 0 || origin == 2) {
            h_mueff_2_num->Fill(t->true_muon_mom);
            break;
          }
        }
        if (t->vtx_resolution > -1 && t->vtx_resolution < 10) nSignalMuonRecoVtxOk++;
        
        h_muon_track_eff->Fill(t->muon_reco_eff);
        h_muon_track_pur->Fill(t->muon_reco_pur);
        
        h_mu_eff_mom->Fill(t->true_muon_mom, t->muon_reco_eff);
        h_mu_pur_mom->Fill(t->true_muon_mom, t->muon_reco_pur);
      }
      else{
        //std::cout << "This is a signal event but the muon was not reconstructed. Event: " << event << std::endl;
      }
    }
    if(t->nupdg == 14 && t->ccnc == 0){
      nNumuCC++;
    }
    
    
    bool isNueCC = false;
    if (t->nupdg == 12 && t->ccnc == 0 && t->fv == 1) {
      isNueCC = true;
      nue_cc_fv++;
    }
    //if (isSignal) std::cout << "IS SIGNAL - event " << t->event << std::endl;
    
    
    if (isSignal) h_vtx_resolution->Fill(t->vtx_resolution);
    //if (isSignal && vtx_resolution > 200 && vtx_resolution < 210) std::cout << "vtx_resolution is fucked for event: " << event << std::endl;
    
    
    
    int n_acpt_tagged_per_event = 0;
    for (int slc = 0; slc < t->nslices; slc ++) {
      
      if (t->slc_origin->at(slc) == 0 || t->slc_origin->at(slc) == 2) {
        h_slice_npfp->Fill(t->slc_npfp->at(slc));
        h_slice_ntrack->Fill(t->slc_ntrack->at(slc));
        //if (t->slc_flsmatch_score->at(slc) < 0.0001) std::cout << ">>>> Neutrino has a low score (" << t->slc_flsmatch_score->at(slc) << "), event " << t->event << std::endl;
      } else {
        h_slice_npfp_others->Fill(t->slc_npfp->at(slc));
        h_slice_ntrack_others->Fill(t->slc_ntrack->at(slc));
      }
      
      if (isSignal) {
        if (t->slc_origin->at(slc) == 0 || t->slc_origin->at(slc) == 2) {
          h_fm_score->Fill(t->slc_flsmatch_score->at(slc));
          h_fm_score_pe->Fill(t->slc_flsmatch_score->at(slc), t->beamfls_pe->at(flashInBeamSpill));
        } else {
          h_fm_score_others->Fill(t->slc_flsmatch_score->at(slc));
        }
        
      }
      
      
      if (t->slc_origin->at(slc) == 0) h_slice_origin->Fill(2);
      if (t->slc_origin->at(slc) == 1) h_slice_origin->Fill(0);
      if (t->slc_origin->at(slc) == 2) h_slice_origin->Fill(1);
      
      if (t->slc_acpt_outoftime->at(slc) == 1) n_acpt_tagged_per_event++;
      
      if ((t->slc_origin->at(slc) == 0 || t->slc_origin->at(slc) == 2) && t->fv == 1) {
        n_slc_nu_origin ++;
        if (t->slc_acpt_outoftime->at(slc) == 1) {
          n_slc_acpt_tag_nu ++;
          std::cout << "A neutrino was acpt tagged in event " << t->event << std::endl;
        }
      }
    }
    h_acpt_tagged->Fill(n_acpt_tagged_per_event);
    
    
    
    
    
    
    // ************************
    //
    //  Selection
    //
    // ************************
    
    
    
    
    
    if (!goodflash) continue;
    
    std::vector<bool> isBackground(t->nslices, false);
    
    // Slice loop
    for (int slc = 0; slc < t->nslices; slc ++) {
      
      bool nu_origin = (t->slc_origin->at(slc) == 0 || t->slc_origin->at(slc) == 2);
      
      // PMTs
      /*
      if (flashInBeamSpill > -1 && t->slc_flsmatch_score->at(slc) > -1
          && t->slc_flsmatch_qllx->at(slc)!= -9999 && t->slc_flsmatch_tpcx->at(slc)!=-9999) {
        hmap_xdiff["total"]->Fill(t->slc_flsmatch_qllx->at(slc) - t->slc_flsmatch_tpcx->at(slc));
        hmap_zdiff["total"]->Fill(t->slc_flsmatch_hypoz->at(slc) - t->beamfls_z->at(flashInBeamSpill));
        if ( isSignal && (t->slc_origin->at(slc) == 0 || t->slc_origin->at(slc) == 2)) {
          hmap_xdiff["signal"]->Fill(t->slc_flsmatch_qllx->at(slc) - t->slc_flsmatch_tpcx->at(slc));
          hmap_zdiff["signal"]->Fill(t->slc_flsmatch_hypoz->at(slc) - t->beamfls_z->at(flashInBeamSpill));
        } else {
          hmap_xdiff["background"]->Fill(t->slc_flsmatch_qllx->at(slc) - t->slc_flsmatch_tpcx->at(slc));
          hmap_zdiff["background"]->Fill(t->slc_flsmatch_hypoz->at(slc) - t->beamfls_z->at(flashInBeamSpill));
        }
      }
      */
      //   nu
      if ( isSignal && (t->slc_origin->at(slc) == 0 || t->slc_origin->at(slc) == 2) && flashInBeamSpill > -1 && t->slc_flsmatch_score->at(slc) > -1){
        for (int pmt = 0; pmt < 32; pmt++) {
          if ((t->slc_flshypo_spec->at(slc))[pmt] < 5 || (t->beamfls_spec->at(flashInBeamSpill))[pmt] < 5) continue;
          double mean = ((t->slc_flshypo_spec->at(slc))[pmt] + (t->beamfls_spec->at(flashInBeamSpill))[pmt]) / 2.;
          h_frac_diff->Fill(pmt, ( (t->slc_flshypo_spec->at(slc))[pmt] - (t->beamfls_spec->at(flashInBeamSpill))[pmt] ) / (mean) );
        }
        if (t->slc_flsmatch_qllx->at(slc)!= -9999 && t->slc_flsmatch_tpcx->at(slc)!=-9999){
          h_xdiff->Fill(t->slc_flsmatch_qllx->at(slc) - t->slc_flsmatch_tpcx->at(slc));
          h_zdiff->Fill(t->slc_flsmatch_hypoz->at(slc) - t->beamfls_z->at(flashInBeamSpill));
        }
      }
      //   others
      if (t->slc_origin->at(slc) == 1 && flashInBeamSpill > -1 && t->slc_flsmatch_score->at(slc) > -1){
        for (int pmt = 0; pmt < 32; pmt++) {
          if ((t->slc_flshypo_spec->at(slc))[pmt] < 5 || (t->beamfls_spec->at(flashInBeamSpill))[pmt] < 5) continue;
          double mean = ((t->slc_flshypo_spec->at(slc))[pmt] + (t->beamfls_spec->at(flashInBeamSpill))[pmt]) / 2.;
          h_frac_diff_others->Fill(pmt, ( (t->slc_flshypo_spec->at(slc))[pmt] - (t->beamfls_spec->at(flashInBeamSpill))[pmt] ) / (mean) );
        }
        if (t->slc_flsmatch_qllx->at(slc)!= -9999 && t->slc_flsmatch_tpcx->at(slc)!=-9999){
          h_xdiff_others->Fill(t->slc_flsmatch_qllx->at(slc) - t->slc_flsmatch_tpcx->at(slc));
          h_zdiff_others->Fill(t->slc_flsmatch_hypoz->at(slc) - t->beamfls_z->at(flashInBeamSpill));
        }
      }
      //  spec
      if (t->event == -1/*150801*/) {
        if (flashInBeamSpill > -1 && t->slc_flsmatch_score->at(slc) > -1){
          for (int pmt = 0; pmt < 32; pmt++) {
            std::cout << "tpcx " << t->slc_flsmatch_tpcx->at(0) << std::endl;
            std::cout << "qllx " << t->slc_flsmatch_qllx->at(0) << std::endl;
            hypo_spec_x[pmt] = pmt;
            hypo_spec_y[pmt] = (t->slc_flshypo_spec->at(0))[pmt]; ;//(t->slc_flshypo_spec->at(3))[pmt];
            std::cout << "SCORE IS " << t->slc_flsmatch_score->at(0) << std::endl;
            meas_spec_x[pmt] = pmt;
            meas_spec_y[pmt] = (t->beamfls_spec->at(flashInBeamSpill))[pmt];
            numc_spec_x[pmt] = pmt;
            //numc_spec_y[pmt] = t->numc_flash_spec->at(pmt);
          }
        }
      }
      
      // CheckVertex
      if (isSignal && nu_origin) {
        h_vtxcheck_angle_good->Fill(t->slc_vtxcheck_angle->at(slc));
      } else {
        h_vtxcheck_angle_bad->Fill(t->slc_vtxcheck_angle->at(slc));
      }
      
      // Track chi2
      h_chi2->Fill(t->slc_kalman_chi2->at(slc)/(double)t->slc_kalman_ndof->at(slc));
      
      // Number of TPCObjects per event
      h_nslices->Fill(t->nslices);
      
      
      isBackground.at(slc) = false;
      
      // ACPT
      if (t->slc_acpt_outoftime->at(slc) == 1) {
        //isBackground.at(slc) = true;
        //continue;
      }
      
      // n slices
      if (t->nslices > 9) {
        //isBackground.at(slc) = true;
        //continue;
      }
      
      /* Dead region
       if (slc_nuvtx_closetodeadregion_w->at(slc) == 1 || slc_nuvtx_closetodeadregion_u->at(slc) == 1 || slc_nuvtx_closetodeadregion_v->at(slc) == 1)
       isBackground.at(slc) = true;
       */
      
      // Crosses top
      //if (slc_crosses_top_boundary->at(slc) == 1)
      //isBackground.at(slc) = true;
      
    }
    
    int n_slc_flsmatch = 0;
    
    // Find slice with maximum score that was not tagged as bkg
    double score_max = -1;
    int scl_ll_max = -1;
    for (int slc = 0; slc < t->nslices; slc ++){
      
      if (t->slc_flsmatch_score->at(slc) > 0.00000001) {
        n_slc_flsmatch++;
      }
      
      //if (isBackground.at(slc)) continue;
      //if (t->slc_flsmatch_qllx->at(slc) - t->slc_flsmatch_tpcx->at(slc) > 20) continue;
      //if(t->slc_flsmatch_hypoz->at(slc) - t->beamfls_z->at(flashInBeamSpill) > 100) continue;
      //if(!t->slc_iscontained->at(slc)) continue;
      
      if(t->slc_flsmatch_score->at(slc) > score_max){
        scl_ll_max = slc;
        score_max = t->slc_flsmatch_score->at(slc);
      }
    }
    
    h_n_slc_flsmatch->Fill(n_slc_flsmatch);
    
    //if (n_slc_flsmatch >= 4) continue;
    
    if (scl_ll_max == -1) continue;
    
    if (score_max <= 0.00000001) continue;
    //if (score_max <= 0.001) continue;
    //std::cout << "passed score" << std::endl;
    //if(t->slc_flsmatch_qllx->at(scl_ll_max) - t->slc_flsmatch_tpcx->at(scl_ll_max) > 20) continue;
    //std::cout << "passed x diff" << std::endl;
    //if(t->slc_flsmatch_qllx->at(scl_ll_max) - t->slc_flsmatch_tpcx->at(scl_ll_max) < -80) continue;
    if(abs(t->slc_flsmatch_hypoz->at(scl_ll_max) - t->beamfls_z->at(flashInBeamSpill)) > 100) continue;
    //std::cout << "passed z diff" << std::endl;
    
    if (isSignal && (t->slc_origin->at(scl_ll_max)==0 || t->slc_origin->at(scl_ll_max)==2)) nSignalFlashMatched ++;
    
    //if(isBackground.at(scl_ll_max)) continue;
    
    if(t->slc_nuvtx_fv->at(scl_ll_max) == 0) continue;
    //std::cout << "passed fv" << std::endl;
    
    if(t->slc_vtxcheck_angle->at(scl_ll_max) > 2.9) continue;
    
    if(t->slc_vtxcheck_angle->at(scl_ll_max) < 0.05 && t->slc_vtxcheck_angle->at(scl_ll_max) !=-9999 ) continue;
    
    if(t->slc_ntrack->at(scl_ll_max) == 0) continue;
    
    if(!t->slc_passed_min_track_quality->at(scl_ll_max)) continue;
    
    if(!t->slc_passed_min_vertex_quality->at(scl_ll_max)) continue;
    
    //if(t->slc_ntrack->at(scl_ll_max) == 1 && t->slc_crosses_top_boundary->at(scl_ll_max) == 1) continue;
    
    //if(t->slc_longesttrack_length->at(scl_ll_max) < 25.) continue;
    
    //if(!t->slc_iscontained->at(scl_ll_max)) continue;
    
    //if(t->slc_crosses_top_boundary->at(scl_ll_max) == 1) continue;
    
    
    
    // Event is selected
    
    hmap_trklen["total"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
    hmap_trkphi["total"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
    hmap_trktheta["total"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
    hmap_multpfp["total"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
    hmap_multtracktol["total"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
    
    //if (flashInBeamSpill > -1 && t->slc_flsmatch_score->at(scl_ll_max) > -1
    //    && t->slc_flsmatch_qllx->at(scl_ll_max)!= -9999 && t->slc_flsmatch_tpcx->at(scl_ll_max)!=-9999) {
      hmap_xdiff["total"]->Fill(t->slc_flsmatch_qllx->at(scl_ll_max) - t->slc_flsmatch_tpcx->at(scl_ll_max));
      hmap_zdiff["total"]->Fill(t->slc_flsmatch_hypoz->at(scl_ll_max) - t->beamfls_z->at(flashInBeamSpill));
      if ( isSignal && (t->slc_origin->at(scl_ll_max) == 0 || t->slc_origin->at(scl_ll_max) == 2)) {
        hmap_xdiff["signal"]->Fill(t->slc_flsmatch_qllx->at(scl_ll_max) - t->slc_flsmatch_tpcx->at(scl_ll_max));
        hmap_zdiff["signal"]->Fill(t->slc_flsmatch_hypoz->at(scl_ll_max) - t->beamfls_z->at(flashInBeamSpill));
      } else {
        hmap_xdiff["background"]->Fill(t->slc_flsmatch_qllx->at(scl_ll_max) - t->slc_flsmatch_tpcx->at(scl_ll_max));
        hmap_zdiff["background"]->Fill(t->slc_flsmatch_hypoz->at(scl_ll_max) - t->beamfls_z->at(flashInBeamSpill));
      }
    //}
    
    bool nu_origin = false;
    if ((t->slc_origin->at(scl_ll_max) == 0 || t->slc_origin->at(scl_ll_max) == 2)) nu_origin = true;
    
    if (isNueCC) {
      nue_cc_selected_total++;
    }
    
    // Signal
    if(nu_origin && t->ccnc==0 && t->nupdg==14 && t->fv==1){
      signal_sel ++;
      h_eff_num->Fill(t->nu_e);
      pEff->Fill(true, t->nu_e);
      hmap_trklen["signal"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["signal"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["signal"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["signal"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["signal"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      //std::cout << "Is signal and is selected. event: " << t->event << std::endl;
      //if (t->slc_mult_track_tolerance->at(scl_ll_max) == 0) std::cout << "Is signal with mult_track_tolerance=0. event: " << t->event << std::endl;
      if (t->slc_origin_extra->at(scl_ll_max) == 0) {
        hmap_trklen["signal_stopmu"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["signal_stopmu"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["signal_stopmu"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["signal_stopmu"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["signal_stopmu"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));

      } else {
        hmap_trklen["signal_nostopmu"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["signal_nostopmu"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["signal_nostopmu"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["signal_nostopmu"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["signal_nostopmu"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      }
    }
    // anumu
    else if(nu_origin && t->ccnc==0 && t->nupdg==-14 && t->fv==1){
      bkg_anumu_sel ++;
      pEff->Fill(false, t->nu_e);
      hmap_trklen["anumu"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["anumu"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["anumu"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["anumu"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["anumu"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
    }
    // nue
    else if(nu_origin && t->ccnc==0 && (t->nupdg==-12 || t->nupdg==12) && t->fv==1){
      bkg_nue_sel ++;
      pEff->Fill(false, t->nu_e);
      hmap_trklen["nue"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["nue"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["nue"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["nue"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["nue"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      if (t->nupdg == 12)
        nue_cc_selected++;
    }
    // nc
    else if(nu_origin && t->ccnc==1 && t->fv==1){
      //std::cout << "origin extra is: " << t->slc_origin_extra->at(scl_ll_max) << std::endl;
      bkg_nc_sel ++;
      pEff->Fill(false, t->nu_e);
      hmap_trklen["nc"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["nc"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["nc"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["nc"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["nc"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      // proton
      if (t->slc_origin_extra->at(scl_ll_max) == 3) {
        hmap_trklen["nc_proton"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["nc_proton"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["nc_proton"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["nc_proton"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["nc_proton"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
        //if (t->slc_longesttrack_length->at(scl_ll_max) > 75.)
         // std::cout << ">>>>>>>>> event " << t->event << std::endl;
      }
      //pion
      else if (t->slc_origin_extra->at(scl_ll_max) == 2) {
        hmap_trklen["nc_pion"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["nc_pion"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["nc_pion"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["nc_pion"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["nc_pion"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      }
      // other
      else {
        //std::cout << "Is NC Other, event: " << t->event << std::endl;
        hmap_trklen["nc_other"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["nc_other"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["nc_other"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["nc_other"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["nc_other"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      }
      //std::cout << "Is a nc but is selected. event: " << t->event << std::endl;
    }
    // outfv
    else if(nu_origin && t->fv==0){
      //std::cout << "origin extra is: " << t->slc_origin_extra->at(scl_ll_max) << std::endl;
      bkg_outfv_sel ++;
      //std::cout << "Is OutFV. event: " << event << std::endl;
      pEff->Fill(false, t->nu_e);
      hmap_trklen["outfv"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["outfv"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["outfv"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["outfv"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["outfv"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      
      if (t->slc_origin_extra->at(scl_ll_max) == 0) {
        hmap_trklen["outfv_stopmu"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["outfv_stopmu"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["outfv_stopmu"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["outfv_stopmu"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["outfv_stopmu"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      } else {
        hmap_trklen["outfv_nostopmu"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["outfv_nostopmu"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["outfv_nostopmu"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["outfv_nostopmu"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["outfv_nostopmu"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      }
      //if (t->slc_mult_track_tolerance->at(scl_ll_max) == 2) std::cout << "Is OutFV with mult_track_tolerance=2. event: " << t->event << std::endl;
    }
    // cosmic
    else{
      bkg_cosmic_sel ++;
      if (t->slc_crosses_top_boundary->at(scl_ll_max) == 1 )
        bkg_cosmic_top_sel++;
      pEff->Fill(false, t->nu_e);
      hmap_trklen["cosmic"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["cosmic"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["cosmic"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["cosmic"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["cosmic"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      //std::cout << "Is a cosmic but is selected. event: " << t->event << std::endl;
      
      if (t->slc_origin_extra->at(scl_ll_max) == 0) {
        hmap_trklen["cosmic_stopmu"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["cosmic_stopmu"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["cosmic_stopmu"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["cosmic_stopmu"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["cosmic_stopmu"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      } else {
        hmap_trklen["cosmic_nostopmu"]->Fill(t->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["cosmic_nostopmu"]->Fill(t->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["cosmic_nostopmu"]->Fill(t->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["cosmic_nostopmu"]->Fill(t->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["cosmic_nostopmu"]->Fill(t->slc_mult_track_tolerance->at(scl_ll_max));
      }
    }
    
    
    
    
    
    //newtree->Fill();
    
    //delete selection;
  } // end of event loop
  
  //Sflashtime      ->Save();
  
  
  //newtree->AutoSave();
  //delete oldfile;
  //delete newfile;
  
  cout << endl << endl << "********************************" << endl;
  
  
  
  
  // ************************
  //
  //  Printing
  //
  // ************************
  
  std::cout << "nsignal is " << nsignal << std::endl;
  std::cout << "signal_sel is " << signal_sel << std::endl;
  std::cout << "bkg_anumu_sel is " << bkg_anumu_sel << std::endl;
  std::cout << "bkg_nue_sel is " << bkg_nue_sel << std::endl;
  std::cout << "bkg_nc_sel is " << bkg_nc_sel << std::endl;
  std::cout << "bkg_outfv_sel is " << bkg_outfv_sel << std::endl;
  std::cout << "bkg_cosmic_sel is " << bkg_cosmic_sel << std::endl << std::endl;
  
  std::cout << "Efficiency: " << signal_sel/(double)nsignal << std::endl;
  std::cout << "Purity: " << signal_sel/(double)(signal_sel+bkg_anumu_sel+bkg_nue_sel+bkg_nc_sel+bkg_outfv_sel+bkg_cosmic_sel) << std::endl;
  std::cout << "Cosmic contamination: " << bkg_cosmic_sel/(double)(bkg_anumu_sel+bkg_nue_sel+bkg_nc_sel+bkg_outfv_sel+bkg_cosmic_sel) << std::endl;
  std::cout << "  of which crossing top: " << bkg_cosmic_top_sel/(double)bkg_cosmic_sel << std::endl;
  std::cout << "NC contamination: " << bkg_nc_sel/(double)(bkg_anumu_sel+bkg_nue_sel+bkg_nc_sel+bkg_outfv_sel+bkg_cosmic_sel) << std::endl;
  std::cout << "OUTFV contamination: " << bkg_outfv_sel/(double)(bkg_anumu_sel+bkg_nue_sel+bkg_nc_sel+bkg_outfv_sel+bkg_cosmic_sel) << std::endl << std::endl;
  
  std::cout << "n events with a flash in the beam spill: " << nEvtsWFlashInBeamSpill << std::endl;
  std::cout << "n events numu CC (all voulumes): " << nNumuCC << std::endl;
  std::cout << " Signal events that have a recon muon: " << nSignalWMuonReco << std::endl;
  std::cout << " Signal events that have a recon muon and a recon vertex 10 cm close in YZ plane: " << nSignalMuonRecoVtxOk << std::endl << std::endl;
  
  std::cout << "Number of signal events that were correctly flash-matched: " << nSignalFlashMatched << std::endl << std::endl;
  
  std::cout << "Number of neutrino origin slices in total: " << n_slc_nu_origin << std::endl;
  std::cout << "Number of neutrino origin slices tagged as cosmic by the ACPT algo in total: " << n_slc_acpt_tag_nu << std::endl << std::endl;
  
  
  
  std::cout << "Number of simulated nue CC in FV: " << nue_cc_fv << std::endl;
  std::cout << "Number of selected nue CC in FV (as such):  " << nue_cc_selected << std::endl;
  std::cout << "Number of selected nue CC in FV (total):  " << nue_cc_selected_total << std::endl << std::endl;
  //std::cout << "\t Ratio: " << (double)nue_cc_selected/(double)nue_cc_fv << std::endl;
  
  
  
  // ************************
  //
  //  Plotting
  //
  // ************************
  
  TString temp2;
  
  TCanvas * canvas_efficiency = new TCanvas();
  TEfficiency* pEff2 = new TEfficiency(*h_eff_num,*h_eff_den);
  pEff2->SetTitle(";True Neutrino Energy [GeV];Efficiency");
  pEff2->SetLineColor(kGreen+3);
  pEff2->SetMarkerColor(kGreen+3);
  pEff2->SetMarkerStyle(20);
  pEff2->SetMarkerSize(0.5);
  pEff2->Draw("AP");
  
  temp2 = "./output/efficiency";
  canvas_efficiency->SaveAs(temp2 + ".pdf");
  canvas_efficiency->SaveAs(temp2 + ".C","C");
  
  TCanvas * canvas_muon_reco_efficiency = new TCanvas();
  TEfficiency* pEff3 = new TEfficiency(*h_mueff_num,*h_mueff_den);
  pEff3->SetTitle(";True Muon Momentum [GeV];Reconstruction Efficiency");
  pEff3->Draw("AP");
  
  temp2 = "./output/muon_reco_efficiency";
  canvas_muon_reco_efficiency->SaveAs(temp2 + ".pdf");
  canvas_muon_reco_efficiency->SaveAs(temp2 + ".C","C");
  
  /*TCanvas *c33 = new TCanvas();
  TEfficiency* pEff4 = new TEfficiency(*h_mueff_2_num,*h_mueff_den);
  pEff4->SetTitle(";True Muon Momentum [GeV];Reconstruction Efficiency");
  pEff4->SetLineColor(kRed);
  pEff4->Draw("P same");
   */
  
  
  new TCanvas();
  h_chi2->Draw("histo");
  
  new TCanvas();
  h_flsTime->Draw("histo");
  
  new TCanvas();
  h_nslices->Draw("histo");
  
  new TCanvas();
  h_vtx_resolution->Draw("histo");
  
  new TCanvas();
  h_frac_diff->Draw("colz");
  
  new TCanvas();
  h_frac_diff_others->Draw("colz");
  
  // PE spec
  new TCanvas();
  TGraph* gr = new TGraph(32,hypo_spec_x,hypo_spec_y);
  TGraph* gr2 = new TGraph(32,meas_spec_x,meas_spec_y);
  TGraph* gr3 = new TGraph(32,numc_spec_x,numc_spec_y);
  gr->SetLineColor(kGreen+2);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(kGreen+2);
  gr->SetMarkerSize(1.2);
  gr->SetMarkerStyle(20);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("PMT ID");
  gr->GetYaxis()->SetTitle("PE Count");
  gr->Draw("ALP");
  gr2->SetLineColor(kBlue+2);
  gr2->SetLineWidth(2);
  gr2->SetMarkerColor(kBlue+2);
  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerStyle(20);
  gr2->SetTitle("");
  gr2->GetXaxis()->SetTitle("PMT ID");
  gr2->GetYaxis()->SetTitle("PE Count");
  gr2->Draw("LP");
  gr3->SetLineColor(kRed+2);
  gr3->SetLineWidth(2);
  gr3->SetMarkerColor(kRed+2);
  gr3->SetMarkerSize(1.2);
  gr3->SetMarkerStyle(20);
  gr3->SetTitle("");
  gr3->GetXaxis()->SetTitle("PMT ID");
  gr3->GetYaxis()->SetTitle("PE Count");
  //gr3->Draw("LP");
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(gr,"Hypo flash","l");
  leg->AddEntry(gr2,"Reco flash","l");
  //leg->AddEntry(gr3,"Neutrino MCFlash","l");
  
  leg->Draw();
  
  new TCanvas();
  h_vtxcheck_angle_good->Scale(1./h_vtxcheck_angle_good->Integral());
  h_vtxcheck_angle_bad->Scale(1./h_vtxcheck_angle_bad->Integral());
  h_vtxcheck_angle_good->Draw("histo");
  h_vtxcheck_angle_bad->Draw("histo same");
  h_vtxcheck_angle_bad->SetLineColor(kRed);
  
  
  new TCanvas();
  //h_muon_track_eff->Draw();
  h_mu_eff_mom->Draw("colz");
  
  new TCanvas();
  //h_muon_track_pur->Draw();
  h_mu_pur_mom->Draw("colz");
  
  new TCanvas();
  h_mumom_nue->Draw("colz");
  
  new TCanvas();
  h_acpt_tagged->Draw("histo");
  
  new TCanvas();
  h_xdiff->Draw("histo");
  h_xdiff_others->Draw("histo same");
  h_xdiff_others->SetLineColor(kRed);
  //TCanvas *c16 = new TCanvas();
  //h_xdiff_others->Draw("histo");
  
  new TCanvas();
  h_zdiff->Draw("histo");
  h_zdiff_others->Draw("histo same");
  h_zdiff_others->SetLineColor(kRed);
  
  
  new TCanvas();
  h_slice_origin->Draw("histo");
  
  new TCanvas();
  h_slice_npfp->DrawNormalized("histo");
  h_slice_npfp_others->DrawNormalized("histo same");
  h_slice_npfp_others->SetLineColor(kRed);
  
  new TCanvas();
  h_slice_ntrack->DrawNormalized("histo");
  h_slice_ntrack_others->DrawNormalized("histo same");
  h_slice_ntrack_others->SetLineColor(kRed);
  
  new TCanvas();
  h_fm_score->Draw("histo");
  h_fm_score_others->Draw("histo same");
  h_fm_score_others->SetLineColor(kRed);
  
  new TCanvas();
  h_n_slc_flsmatch->Draw("histo");
  
  new TCanvas();
  h_fm_score_pe->Draw("colz");
  
  
  
  TCanvas * final1 = new TCanvas();
  THStack *hs_trklen = new THStack("hs_trklen",";Candidate Track Length [cm]; Selected Events");
  //DrawTHStack(hs_trklen, pot_scaling, _breakdownPlots, hmap_trklen);
  
  
  // Construct legend
  TLegend* leg2;
  if (_breakdownPlots){
    leg2 = new TLegend(0.56,0.37,0.82,0.82,NULL,"brNDC");
  } else {
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
  }
  std::stringstream sstm;
  // numu
  if (_breakdownPlots) {
  leg2->AddEntry(hmap_trklen["signal_stopmu"],"#nu_{#mu} CC (stopping #mu)","f");
  leg2->AddEntry(hmap_trklen["signal_nostopmu"],"#nu_{#mu} CC (other)","f");
  } else {
    sstm << "#nu_{#mu} CC (signal), " << std::setprecision(2)  << hmap_trklen["signal"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
    leg2->AddEntry(hmap_trklen["signal"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  
  // nue
  sstm << "#nu_{e}, #bar{#nu}_{e} CC, " << std::setprecision(2)  << hmap_trklen["nue"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
  leg2->AddEntry(hmap_trklen["nue"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // anumu
  sstm << "#bar{#nu}_{#mu} CC, " << std::setprecision(2)  << hmap_trklen["anumu"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
  leg2->AddEntry(hmap_trklen["anumu"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // nc, outfv, cosmic
  if (_breakdownPlots) {
  leg2->AddEntry(hmap_trklen["nc_other"],"NC (other)","f");
  leg2->AddEntry(hmap_trklen["nc_pion"],"NC (pion)","f");
  leg2->AddEntry(hmap_trklen["nc_proton"],"NC (proton)","f");
  leg2->AddEntry(hmap_trklen["outfv_stopmu"],"OUTFV (stopping #mu)","f");
  leg2->AddEntry(hmap_trklen["outfv_nostopmu"],"OUTFV (other)","f");
  leg2->AddEntry(hmap_trklen["cosmic_stopmu"],"Cosmic (stopping #mu)","f");
  leg2->AddEntry(hmap_trklen["cosmic_nostopmu"],"Cosmic (other)","f");
  } else {
    sstm << "NC, " << std::setprecision(2)  << hmap_trklen["nc"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
    leg2->AddEntry(hmap_trklen["nc"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "OUTFV, " << std::setprecision(2)  << hmap_trklen["outfv"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
    leg2->AddEntry(hmap_trklen["outfv"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "Cosmic, " << std::setprecision(2)  << hmap_trklen["cosmic"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
    leg2->AddEntry(hmap_trklen["cosmic"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  leg2->AddEntry(hmap_trklen["total"],"MC Stat Unc.","f");
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/trklen";
  final1->SaveAs(temp2 + ".pdf");
  final1->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final2 = new TCanvas();
  THStack *hs_trkphi = new THStack("hs_trkphi",";Candidate Track #phi; Selected Events");
  //DrawTHStack(hs_trkphi, pot_scaling, _breakdownPlots, hmap_trkphi);
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/trkphi";
  final2->SaveAs(temp2 + ".pdf");
  final2->SaveAs(temp2 + ".C","C");
  
  
  
  TCanvas * final3 = new TCanvas();
  THStack *hs_trktheta = new THStack("hs_trktheta",";Candidate Track cos(#theta); Selected Events");
  //DrawTHStack(hs_trktheta, pot_scaling, _breakdownPlots, hmap_trktheta);
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/trktheta";
  final3->SaveAs(temp2 + ".pdf");
  final3->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final4 = new TCanvas();
  THStack *hs_multpfp = new THStack("hs_multpfp",";PFP Multiplicity; Selected Events");
  //DrawTHStack(hs_multpfp, pot_scaling, _breakdownPlots, hmap_multpfp);
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/multpfp";
  final4->SaveAs(temp2 + ".pdf");
  final4->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final5 = new TCanvas();
  THStack *hs_multtracktol = new THStack("hs_multtracktol",";Track Multiplicity (5 cm); Selected Events");
  //DrawTHStack(hs_multtracktol, pot_scaling, _breakdownPlots, hmap_multtracktol);
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/multtracktol";
  final5->SaveAs(temp2 + ".pdf");
  final5->SaveAs(temp2 + ".C","C");
  
  
  
  
  file_out->cd();
  for (auto iter : hmap_trklen) {
    iter.second->Write();
  }
  h_pot->Write();
  h_nevts->Write();

  file_out->WriteObject(&hmap_trklen, "hmap_trklen");
  file_out->WriteObject(&hmap_trktheta, "hmap_trktheta");
  file_out->WriteObject(&hmap_trkphi, "hmap_trkphi");
  file_out->WriteObject(&hmap_multpfp, "hmap_multpfp");
  file_out->WriteObject(&hmap_multtracktol, "hmap_multtracktol");
  
  file_out->WriteObject(&hmap_xdiff, "hmap_xdiff");
  file_out->WriteObject(&hmap_zdiff, "hmap_zdiff");
  
  h_flsTime->Write();
  h_flsTime_wcut->Write();

  file_out->Write();
  
  file_out->Close();
  
  
  
  
  
  
  
  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  //rootapp->Run();
  //rootapp->Terminate(0);
  
  return 0;
}











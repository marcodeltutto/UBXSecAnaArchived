#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>

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
//#include "PlotHandler.hpp"
//#include "SelectionTools.hpp"

const bool _breackdownPlots = true;

const double _beamSpillStarts = 3.2; // us
const double _beamSpillEnds   = 4.8; // us

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
  
  //at->fChain->SetBranchStatus("",1);
}


double CalcLength(const double& x_1, const double& y_1, const double& z_1, const double& x_2, const double& y_2, const double& z_2) {
  return sqrt(pow(x_1-x_2, 2) + pow(y_1-y_2, 2) + pow(z_1-z_2, 2));
}

//____________________________________________________________________________________________________
void DrawTHStack(THStack *hs_trklen,
                   double pot_scaling,
                   std::map<std::string,TH1D*> themap);

//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
int main(int argc, char* argv[]) {
  
  clock_t begin = clock();
  
  if (argc != 4) {
    cout << "Provide 3 arguments!" << endl;
    cout << "./NuMICCInclusive bnb pot maxentries." << endl;
    exit(0);
  }
  
  string filen      = argv[1];
  string dopot      = argv[2];
  int    maxEntries = atoi(argv[3]);
  cout << endl << "File name " << filen << " beam." << endl;
  
  TApplication* rootapp = new TApplication("ROOT Application",&argc, argv);
  //gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS
  gROOT->ProcessLine(".x rootlogon.C");
  
  std::cout << "Opening output file ubxsecana_output.root." << std::endl;
  TFile *file_out = new TFile("ubxsecana_output.root","NEW");
  if ( file_out->IsOpen() )
    std::cout << "File opened successfully" << std::endl;
  
  string pattern = filen;
  //if(beam == "numi") pattern = "/data/t2k/lar/uboone/prodgenie_numi_nu_uboone_MCC7/prodgenie_numi_nu_cosmic_uboone_merged_gen_g4_detsim_reco1_reco2_ana.root";
  //if(beam == "bnb")  pattern = "/pnfs/uboone/persistent/users/aschu/MC_BNB_Cosmic/prodgenie_bnb_nu_cosmic_uboone_v05_08_00_at.root";
  //if(beam == "bnb")  pattern = "files/output9850.root";
  
  
  bool evalPOT = false;
  double totalPOT = 0.;
  
  if (dopot == "pot") evalPOT = true;
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
  
  
  TChain *chain_ubxsec;
  chain_ubxsec = new TChain("UBXSec/tree");
  chain_ubxsec->Add(pattern.c_str());
  
  cout << "Using file: " << pattern << endl;
  
  int Nfiles = chain_ubxsec->GetNtrees();
  cout << "Number of files: " << Nfiles << endl;
  
  int evts = chain_ubxsec -> GetEntries();
  cout << "Number of events used is: " << evts << endl;
  
  AnaTree * at = new AnaTree(chain_ubxsec);
  ActivateBranches(at);
  
  std::string TrackProdName  = "pandoraNu";
  std::string VertexProdName = "pandoraNu";
  
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
  
  /*
   TH1D* h_bkg_cosmic_sel = new TH1D("h_bkg_cosmic_sel", "h_bkg_cosmic_sel", 6, 0, 4);
   TH1D* h_bkg_cosmic_all = new TH1D("h_bkg_cosmic_all", "h_bkg_cosmic_all", 6, 0, 4);
   
   TH1D* h_bkg_outfv_sel = new TH1D("h_bkg_outfv_sel", "h_bkg_outfv_sel", 6, 0, 4);
   TH1D* h_bkg_outfv_all = new TH1D("h_bkg_outfv_all", "h_bkg_outfv_all", 6, 0, 4);
   
   TH1D* h_bkg_nc_sel = new TH1D("h_bkg_nc_sel", "h_bkg_nc_sel", 6, 0, 4);
   TH1D* h_bkg_nc_all = new TH1D("h_bkg_nc_all", "h_bkg_nc_all", 6, 0, 4);
   
   TH1D* h_bkg_nue_sel = new TH1D("h_bkg_nue_sel", "h_bkg_nue_sel", 6, 0, 4);
   TH1D* h_bkg_nue_all = new TH1D("h_bkg_nue_all", "h_bkg_nue_all", 6, 0, 4);
   
   TH1D* h_bkg_anumu_sel = new TH1D("h_bkg_anumu_sel", "h_bkg_anumu_sel", 6, 0, 4);
   TH1D* h_bkg_anumu_all = new TH1D("h_bkg_anumu_all", "h_bkg_anumu_all", 6, 0, 4);
   */
  
  
  TH1D* h_chi2 = new TH1D("h_chi2", "h_chi2", 50, 0, 50);
  TH1D* h_flsTime = new TH1D("h_flsTime", ";Flash time w.r.t. trigger [#mus];Events", 100, 0, 25);
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

  
  int barWidth = 70;
  
  if(maxEntries > 0.) evts = maxEntries;
  
  for(int i = 0; i < evts; i++) {
    
    DrawProgressBar((double)i/(double)evts, barWidth);
    
    chain_ubxsec->GetEntry(i);
    
    //SelectionTools * selection = new SelectionTools(at);
    
    //cout << "***** Event " << i << endl;
    
    // ************************
    //
    // Preliminary distributions
    //
    // ************************
    
    /*
     // Flashes
     for (int fls = 0; fls < at->no_flashes; fls++) {
     Sflashtime      ->Fill(at->flash_time[fls]);
     Sflashpe        ->Fill(at->flash_pe[fls]);
     Sflashycenter   ->Fill(at->flash_ycenter[fls]);
     Sflashzcenter   ->Fill(at->flash_zcenter[fls]);
     Sflashtimewidth ->Fill(at->flash_timewidth[fls]);
     if (at->flash_pe[fls] > 50) {
     Sflashtime50pe->Fill(at->flash_time[fls]);
     int k = 0;
     double distance = sqrt(pow(54.999*100-at->vx_flux[k],2)+pow(74.461*100-at->vy_flux[k],2)+pow(677.611*100-at->vz_flux[k],2));
     Sfls_timeVSnu_distance->Fill(at->flash_time[fls],distance);
     }
     Sfls_timeVSpe->Fill(at->flash_time[fls], at->flash_pe[fls]);
     }
     bool doneForThisEvent = false;
     bool doneForThisEvent_cosmic = false;
     for (int geantpar = 0; geantpar < at->geant_list_size; geantpar++) {
     if (!doneForThisEvent && at->origin[geantpar] == 1 && at->process_primary[geantpar]==1) {
     Sgeanttruetime_neutrino->Fill(at->StartT[geantpar]);
     doneForThisEvent = true;
     }
     if (!doneForThisEvent_cosmic &&at->origin[geantpar] == 2 && at->process_primary[geantpar]==1) {
     Sgeanttruetime_cosmic->Fill(at->StartT[geantpar]);
     doneForThisEvent_cosmic = true;
     }
     }
     */
    
    
    
    if (at->nbeamfls == 0) continue;
    int flashInBeamSpill = -1;
    
    bool goodflash = false;
    for (int fls = 0; fls < at->nbeamfls; fls ++){
      h_flsTime->Fill(at->beamfls_time->at(fls));
      if (at->beamfls_time->at(fls) > _beamSpillStarts && at->beamfls_time->at(fls) < _beamSpillEnds) {
        
        flashInBeamSpill = fls;
        if (at->beamfls_pe->at(fls) >= 50) {
          goodflash = true;
          nEvtsWFlashInBeamSpill++;
        }
      }
    }
    
    if (flashInBeamSpill == -1) continue;
    
    
    
    
    
    
    bool isSignal = false;
    if (at->nupdg == 14 && at->ccnc == 0 && at->fv == 1){
      nsignal++;
      isSignal = true;
      h_eff_den->Fill(at->nu_e);
      h_mueff_den->Fill(at->true_muon_mom);
      
      if (at->muon_is_reco){
        h_mumom_nue->Fill(at->nu_e, at->true_muon_mom);
        nSignalWMuonReco++;
        h_mueff_num->Fill(at->true_muon_mom);
        for (auto origin : *at->slc_origin){
          if (origin == 0 || origin == 2) {
            h_mueff_2_num->Fill(at->true_muon_mom);
            break;
          }
        }
        if (at->vtx_resolution > -1 && at->vtx_resolution < 10) nSignalMuonRecoVtxOk++;
        
        h_muon_track_eff->Fill(at->muon_reco_eff);
        h_muon_track_pur->Fill(at->muon_reco_pur);
        
        h_mu_eff_mom->Fill(at->true_muon_mom, at->muon_reco_eff);
        h_mu_pur_mom->Fill(at->true_muon_mom, at->muon_reco_pur);
      }
      else{
        //std::cout << "This is a signal event but the muon was not reconstructed. Event: " << event << std::endl;
      }
    }
    if(at->nupdg == 14 && at->ccnc == 0){
      nNumuCC++;
    }
    
    
    bool isNueCC = false;
    if (at->nupdg == 12 && at->ccnc == 0 && at->fv == 1) {
      isNueCC = true;
      nue_cc_fv++;
    }
    //if (isSignal) std::cout << "IS SIGNAL - event " << at->event << std::endl;
    
    
    if (isSignal) h_vtx_resolution->Fill(at->vtx_resolution);
    //if (isSignal && vtx_resolution > 200 && vtx_resolution < 210) std::cout << "vtx_resolution is fucked for event: " << event << std::endl;
    
    
    
    int n_acpt_tagged_per_event = 0;
    for (int slc = 0; slc < at->nslices; slc ++) {
      
      if (at->slc_origin->at(slc) == 0 || at->slc_origin->at(slc) == 2) {
        h_slice_npfp->Fill(at->slc_npfp->at(slc));
        h_slice_ntrack->Fill(at->slc_ntrack->at(slc));
        //if (at->slc_flsmatch_score->at(slc) < 0.0001) std::cout << ">>>> Neutrino has a low score (" << at->slc_flsmatch_score->at(slc) << "), event " << at->event << std::endl;
      } else {
        h_slice_npfp_others->Fill(at->slc_npfp->at(slc));
        h_slice_ntrack_others->Fill(at->slc_ntrack->at(slc));
      }
      
      if (isSignal) {
        if (at->slc_origin->at(slc) == 0 || at->slc_origin->at(slc) == 2) {
          h_fm_score->Fill(at->slc_flsmatch_score->at(slc));
          h_fm_score_pe->Fill(at->slc_flsmatch_score->at(slc), at->beamfls_pe->at(flashInBeamSpill));
        } else {
          h_fm_score_others->Fill(at->slc_flsmatch_score->at(slc));
        }
        
      }
      
      
      if (at->slc_origin->at(slc) == 0) h_slice_origin->Fill(2);
      if (at->slc_origin->at(slc) == 1) h_slice_origin->Fill(0);
      if (at->slc_origin->at(slc) == 2) h_slice_origin->Fill(1);
      
      if (at->slc_acpt_outoftime->at(slc) == 1) n_acpt_tagged_per_event++;
      
      if ((at->slc_origin->at(slc) == 0 || at->slc_origin->at(slc) == 2) && at->fv == 1) {
        n_slc_nu_origin ++;
        if (at->slc_acpt_outoftime->at(slc) == 1) {
          n_slc_acpt_tag_nu ++;
          std::cout << "A neutrino was acpt tagged in event " << at->event << std::endl;
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
    
    std::vector<bool> isBackground(at->nslices, false);
    
    // Slice loop
    for (int slc = 0; slc < at->nslices; slc ++) {
      
      bool nu_origin = (at->slc_origin->at(slc) == 0 || at->slc_origin->at(slc) == 2);
      
      // PMTs
      //   nu
      if ( isSignal && (at->slc_origin->at(slc) == 0 || at->slc_origin->at(slc) == 2) && flashInBeamSpill > -1 && at->slc_flsmatch_score->at(slc) > -1){
        for (int pmt = 0; pmt < 32; pmt++) {
          if ((at->slc_flshypo_spec->at(slc))[pmt] < 5 || (at->beamfls_spec->at(flashInBeamSpill))[pmt] < 5) continue;
          double mean = ((at->slc_flshypo_spec->at(slc))[pmt] + (at->beamfls_spec->at(flashInBeamSpill))[pmt]) / 2.;
          h_frac_diff->Fill(pmt, ( (at->slc_flshypo_spec->at(slc))[pmt] - (at->beamfls_spec->at(flashInBeamSpill))[pmt] ) / (mean) );
        }
        if (at->slc_flsmatch_qllx->at(slc)!= -9999 && at->slc_flsmatch_tpcx->at(slc)!=-9999){
          h_xdiff->Fill(at->slc_flsmatch_qllx->at(slc) - at->slc_flsmatch_tpcx->at(slc));
          h_zdiff->Fill(at->slc_flsmatch_hypoz->at(slc) - at->beamfls_z->at(flashInBeamSpill));
        }
      }
      //   others
      if (at->slc_origin->at(slc) == 1 && flashInBeamSpill > -1 && at->slc_flsmatch_score->at(slc) > -1){
        for (int pmt = 0; pmt < 32; pmt++) {
          if ((at->slc_flshypo_spec->at(slc))[pmt] < 5 || (at->beamfls_spec->at(flashInBeamSpill))[pmt] < 5) continue;
          double mean = ((at->slc_flshypo_spec->at(slc))[pmt] + (at->beamfls_spec->at(flashInBeamSpill))[pmt]) / 2.;
          h_frac_diff_others->Fill(pmt, ( (at->slc_flshypo_spec->at(slc))[pmt] - (at->beamfls_spec->at(flashInBeamSpill))[pmt] ) / (mean) );
        }
        if (at->slc_flsmatch_qllx->at(slc)!= -9999 && at->slc_flsmatch_tpcx->at(slc)!=-9999){
          h_xdiff_others->Fill(at->slc_flsmatch_qllx->at(slc) - at->slc_flsmatch_tpcx->at(slc));
          h_zdiff_others->Fill(at->slc_flsmatch_hypoz->at(slc) - at->beamfls_z->at(flashInBeamSpill));
        }
      }
      //  spec
      if (at->event == -1/*150801*/) {
        if (flashInBeamSpill > -1 && at->slc_flsmatch_score->at(slc) > -1){
          for (int pmt = 0; pmt < 32; pmt++) {
            std::cout << "tpcx " << at->slc_flsmatch_tpcx->at(0) << std::endl;
            std::cout << "qllx " << at->slc_flsmatch_qllx->at(0) << std::endl;
            hypo_spec_x[pmt] = pmt;
            hypo_spec_y[pmt] = (at->slc_flshypo_spec->at(0))[pmt]; ;//(at->slc_flshypo_spec->at(3))[pmt];
            std::cout << "SCORE IS " << at->slc_flsmatch_score->at(0) << std::endl;
            meas_spec_x[pmt] = pmt;
            meas_spec_y[pmt] = (at->beamfls_spec->at(flashInBeamSpill))[pmt];
            numc_spec_x[pmt] = pmt;
            //numc_spec_y[pmt] = at->numc_flash_spec->at(pmt);
          }
        }
      }
      
      /* CheckVertex
       if (nu_origin && at->fv == 1 && at->ccnc == 0) {
       if (at->vtx_resolution <= 10.){
       h_vtxcheck_angle_good->Fill(at->slc_vtxcheck_angle->at(slc));
       //if (at->slc_vtxcheck_angle->at(slc) > 3)
       //std::cout << "Angle is about 180 for a good vertex for event: " << at->event << std::endl;
       }
       if (at->vtx_resolution > 10.){
       h_vtxcheck_angle_bad->Fill(at->slc_vtxcheck_angle->at(slc));
       //if (at->slc_vtxcheck_angle->at(slc) < 2)
       //std::cout << "Angle is not 180 for a bad vertex for event: " << at->event << std::endl;
       }
       }*/
      // CheckVertex
      if (isSignal && nu_origin) {
        h_vtxcheck_angle_good->Fill(at->slc_vtxcheck_angle->at(slc));
      } else {
        h_vtxcheck_angle_bad->Fill(at->slc_vtxcheck_angle->at(slc));
      }
      
      
      h_chi2->Fill(at->slc_kalman_chi2->at(slc)/(double)at->slc_kalman_ndof->at(slc));
      //if (at->slc_kalman_chi2->at(slc)/(double)at->slc_kalman_ndof->at(slc) > 15)
      //std::cout << "This event (" << event << ") has a track with chi2/ndof = " << slc_kalman_chi2->at(slc)/(double)slc_kalman_ndof->at(slc) << std::endl;
      h_nslices->Fill(at->nslices);
      
      isBackground.at(slc) = false;
      
      
      // ACPT
      if (at->slc_acpt_outoftime->at(slc) == 1) {
        //isBackground.at(slc) = true;
        //continue;
      }
      
      // n slices
      if (at->nslices > 9) {
        //isBackground.at(slc) = true;
        //continue;
      }
      
      /* Cosmic flash match
       if (slc_flsmatch_cosmic_t0->at(slc) < 3 || slc_flsmatch_cosmic_t0->at(slc) > 5) {
       
       if (slc_flsmatch_cosmic_score->at(slc) > 0.01) {
       
       if (slc_flsmatch_cosmic_score->at(slc) > slc_flsmatch_score->at(slc))
       isBackground.at(slc) = true;
       }
       }*/
      
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
    for (int slc = 0; slc < at->nslices; slc ++){
      
      if (at->slc_flsmatch_score->at(slc) > 0.00000001) {
        n_slc_flsmatch++;
      }
      
      //if (isBackground.at(slc)) continue;
      //if (at->slc_flsmatch_qllx->at(slc) - at->slc_flsmatch_tpcx->at(slc) > 20) continue;
      //if(at->slc_flsmatch_hypoz->at(slc) - at->beamfls_z->at(flashInBeamSpill) > 100) continue;
      //if(!at->slc_iscontained->at(slc)) continue;
      
      if(at->slc_flsmatch_score->at(slc) > score_max){
        scl_ll_max = slc;
        score_max = at->slc_flsmatch_score->at(slc);
      }
    }
    
    h_n_slc_flsmatch->Fill(n_slc_flsmatch);
    
    //if (n_slc_flsmatch >= 4) continue;
    
    if (scl_ll_max == -1) continue;
    
    if (score_max <= 0.00000001) continue;
    //if (score_max <= 0.001) continue;
    //std::cout << "passed score" << std::endl;
    if(at->slc_flsmatch_qllx->at(scl_ll_max) - at->slc_flsmatch_tpcx->at(scl_ll_max) > 20) continue;
    //std::cout << "passed x diff" << std::endl;
    //if(at->slc_flsmatch_qllx->at(scl_ll_max) - at->slc_flsmatch_tpcx->at(scl_ll_max) < -80) continue;
    if(at->slc_flsmatch_hypoz->at(scl_ll_max) - at->beamfls_z->at(flashInBeamSpill) > 100) continue;
    //std::cout << "passed z diff" << std::endl;
    
    if (isSignal && (at->slc_origin->at(scl_ll_max)==0 || at->slc_origin->at(scl_ll_max)==2)) nSignalFlashMatched ++;
    
    //if(isBackground.at(scl_ll_max)) continue;
    
    if(at->slc_nuvtx_fv->at(scl_ll_max) == 0) continue;
    //std::cout << "passed fv" << std::endl;
    
    if(at->slc_vtxcheck_angle->at(scl_ll_max) > 2.9) continue;
    
    if(at->slc_vtxcheck_angle->at(scl_ll_max) < 0.05 && at->slc_vtxcheck_angle->at(scl_ll_max) !=-9999 ) continue;
    
    if(at->slc_ntrack->at(scl_ll_max) == 0) continue;
    
    if(!at->slc_passed_min_track_quality->at(scl_ll_max)) continue;
    
    if(!at->slc_passed_min_vertex_quality->at(scl_ll_max)) continue;
    
    //if(at->slc_ntrack->at(scl_ll_max) == 1 && at->slc_crosses_top_boundary->at(scl_ll_max) == 1) continue;
    
    //if(at->slc_longesttrack_length->at(scl_ll_max) < 25.) continue;
    
    //if(!at->slc_iscontained->at(scl_ll_max)) continue;
    
    //if(at->slc_crosses_top_boundary->at(scl_ll_max) == 1) continue;
    
    
    
    // Event is selected
    
    hmap_trklen["total"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
    hmap_trkphi["total"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
    hmap_trktheta["total"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
    hmap_multpfp["total"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
    hmap_multtracktol["total"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
    
    bool nu_origin = false;
    if ((at->slc_origin->at(scl_ll_max) == 0 || at->slc_origin->at(scl_ll_max) == 2)) nu_origin = true;
    
    if (isNueCC) {
      nue_cc_selected_total++;
    }
    
    // Signal
    if(nu_origin && at->ccnc==0 && at->nupdg==14 && at->fv==1){
      signal_sel ++;
      h_eff_num->Fill(at->nu_e);
      pEff->Fill(true, at->nu_e);
      hmap_trklen["signal"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["signal"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["signal"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["signal"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["signal"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      //std::cout << "Is signal and is selected. event: " << at->event << std::endl;
      //if (at->slc_mult_track_tolerance->at(scl_ll_max) == 0) std::cout << "Is signal with mult_track_tolerance=0. event: " << at->event << std::endl;
      if (at->slc_origin_extra->at(scl_ll_max) == 0) {
        hmap_trklen["signal_stopmu"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["signal_stopmu"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["signal_stopmu"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["signal_stopmu"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["signal_stopmu"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));

      } else {
        hmap_trklen["signal_nostopmu"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["signal_nostopmu"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["signal_nostopmu"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["signal_nostopmu"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["signal_nostopmu"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      }
    }
    // anumu
    else if(nu_origin && at->ccnc==0 && at->nupdg==-14 && at->fv==1){
      bkg_anumu_sel ++;
      pEff->Fill(false, at->nu_e);
      hmap_trklen["anumu"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["anumu"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["anumu"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["anumu"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["anumu"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
    }
    // nue
    else if(nu_origin && at->ccnc==0 && (at->nupdg==-12 || at->nupdg==12) && at->fv==1){
      bkg_nue_sel ++;
      pEff->Fill(false, at->nu_e);
      hmap_trklen["nue"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["nue"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["nue"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["nue"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["nue"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      if (at->nupdg == 12)
        nue_cc_selected++;
    }
    // nc
    else if(nu_origin && at->ccnc==1 && at->fv==1){
      //std::cout << "origin extra is: " << at->slc_origin_extra->at(scl_ll_max) << std::endl;
      bkg_nc_sel ++;
      pEff->Fill(false, at->nu_e);
      hmap_trklen["nc"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["nc"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["nc"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["nc"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["nc"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      // proton
      if (at->slc_origin_extra->at(scl_ll_max) == 3) {
        hmap_trklen["nc_proton"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["nc_proton"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["nc_proton"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["nc_proton"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["nc_proton"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
        //if (at->slc_longesttrack_length->at(scl_ll_max) > 75.)
         // std::cout << ">>>>>>>>> event " << at->event << std::endl;
      }
      //pion
      else if (at->slc_origin_extra->at(scl_ll_max) == 2) {
        hmap_trklen["nc_pion"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["nc_pion"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["nc_pion"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["nc_pion"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["nc_pion"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      }
      // other
      else {
        //std::cout << "Is NC Other, event: " << at->event << std::endl;
        hmap_trklen["nc_other"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["nc_other"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["nc_other"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["nc_other"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["nc_other"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      }
      //std::cout << "Is a nc but is selected. event: " << at->event << std::endl;
    }
    // outfv
    else if(nu_origin && at->fv==0){
      //std::cout << "origin extra is: " << at->slc_origin_extra->at(scl_ll_max) << std::endl;
      bkg_outfv_sel ++;
      //std::cout << "Is OutFV. event: " << event << std::endl;
      pEff->Fill(false, at->nu_e);
      hmap_trklen["outfv"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["outfv"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["outfv"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["outfv"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["outfv"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      
      if (at->slc_origin_extra->at(scl_ll_max) == 0) {
        hmap_trklen["outfv_stopmu"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["outfv_stopmu"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["outfv_stopmu"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["outfv_stopmu"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["outfv_stopmu"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      } else {
        hmap_trklen["outfv_nostopmu"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["outfv_nostopmu"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["outfv_nostopmu"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["outfv_nostopmu"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["outfv_nostopmu"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      }
      //if (at->slc_mult_track_tolerance->at(scl_ll_max) == 2) std::cout << "Is OutFV with mult_track_tolerance=2. event: " << at->event << std::endl;
    }
    // cosmic
    else{
      bkg_cosmic_sel ++;
      if (at->slc_crosses_top_boundary->at(scl_ll_max) == 1 )
        bkg_cosmic_top_sel++;
      pEff->Fill(false, at->nu_e);
      hmap_trklen["cosmic"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      hmap_trkphi["cosmic"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      hmap_trktheta["cosmic"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      hmap_multpfp["cosmic"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
      hmap_multtracktol["cosmic"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      //std::cout << "Is a cosmic but is selected. event: " << at->event << std::endl;
      
      if (at->slc_origin_extra->at(scl_ll_max) == 0) {
        hmap_trklen["cosmic_stopmu"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["cosmic_stopmu"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["cosmic_stopmu"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["cosmic_stopmu"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["cosmic_stopmu"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      } else {
        hmap_trklen["cosmic_nostopmu"]->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        hmap_trkphi["cosmic_nostopmu"]->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        hmap_trktheta["cosmic_nostopmu"]->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        hmap_multpfp["cosmic_nostopmu"]->Fill(at->slc_mult_pfp->at(scl_ll_max));
        hmap_multtracktol["cosmic_nostopmu"]->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
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
  std::cout << "Number of selected nue CC in FV (total):  " << nue_cc_selected_total << std::endl;
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
  DrawTHStack(hs_trklen,
                targetPOT/totalPOT,
                hmap_trklen);
  
  // Construct legend
  TLegend* leg2;
  if (_breackdownPlots){
    leg2 = new TLegend(0.56,0.37,0.82,0.82,NULL,"brNDC");
  } else {
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
  }
  std::stringstream sstm;
  // numu
  if (_breackdownPlots) {
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
  if (_breackdownPlots) {
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
  DrawTHStack(hs_trkphi,
              targetPOT/totalPOT,
              hmap_trkphi);
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/trkphi";
  final2->SaveAs(temp2 + ".pdf");
  final2->SaveAs(temp2 + ".C","C");
  
  
  
  TCanvas * final3 = new TCanvas();
  THStack *hs_trktheta = new THStack("hs_trktheta",";Candidate Track cos(#theta); Selected Events");
  DrawTHStack(hs_trktheta,
              targetPOT/totalPOT,
              hmap_trktheta);
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/trktheta";
  final3->SaveAs(temp2 + ".pdf");
  final3->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final4 = new TCanvas();
  THStack *hs_multpfp = new THStack("hs_multpfp",";PFP Multiplicity; Selected Events");
  DrawTHStack(hs_multpfp,
              targetPOT/totalPOT,
              hmap_multpfp);
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/multpfp";
  final4->SaveAs(temp2 + ".pdf");
  final4->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final5 = new TCanvas();
  THStack *hs_multtracktol = new THStack("hs_multtracktol",";Track Multiplicity (5 cm); Selected Events");
  DrawTHStack(hs_multtracktol,
              targetPOT/totalPOT,
              hmap_multtracktol);
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/multtracktol";
  final5->SaveAs(temp2 + ".pdf");
  final5->SaveAs(temp2 + ".C","C");
  
  
  
  
  file_out->cd();
  //h_trklen_total->Write();
  file_out->Close();
  
  
  
  
  
  
  
  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  rootapp->Run();
  rootapp->Terminate(0);
  
  return 0;
}








//**********************************************

void DrawTHStack(THStack *hs_trklen,
                   double pot_scaling,
                   std::map<std::string,TH1D*> themap){
  
  
  for (auto iter : themap) {
    iter.second->Scale(pot_scaling);
  }
  
  
  if (_breackdownPlots) {
    themap["cosmic_nostopmu"]->SetLineColor(kBlue+2);
    themap["cosmic_nostopmu"]->SetFillColor(kBlue+2);
    hs_trklen->Add(themap["cosmic_nostopmu"]);
    themap["cosmic_stopmu"]->SetLineColor(kBlue);
    themap["cosmic_stopmu"]->SetFillColor(kBlue);
    hs_trklen->Add(themap["cosmic_stopmu"]);
    themap["outfv_nostopmu"]->SetLineColor(kGreen+3);
    themap["outfv_nostopmu"]->SetFillColor(kGreen+3);
    hs_trklen->Add(themap["outfv_nostopmu"]);
    themap["outfv_stopmu"]->SetLineColor(kGreen+2);
    themap["outfv_stopmu"]->SetFillColor(kGreen+2);
    hs_trklen->Add(themap["outfv_stopmu"]);
    themap["nc_proton"]->SetLineColor(kGray+2);
    themap["nc_proton"]->SetFillColor(kGray+2);
    hs_trklen->Add(themap["nc_proton"]);
    themap["nc_pion"]->SetLineColor(kGray+1);
    themap["nc_pion"]->SetFillColor(kGray+1);
    hs_trklen->Add(themap["nc_pion"]);
    themap["nc_other"]->SetLineColor(kGray);
    themap["nc_other"]->SetFillColor(kGray);
    hs_trklen->Add(themap["nc_other"]);
  }
  else {
    themap["cosmic"]->SetLineColor(kBlue+2);
    themap["cosmic"]->SetFillColor(kBlue+2);
    hs_trklen->Add(themap["cosmic"]);
    themap["outfv"]->SetLineColor(kGreen+2);
    themap["outfv"]->SetFillColor(kGreen+2);
    hs_trklen->Add(themap["outfv"]);
    themap["nc"]->SetLineColor(kGray);
    themap["nc"]->SetFillColor(kGray);
    hs_trklen->Add(themap["nc"]);
  }
  
  themap["anumu"]->SetLineColor(kOrange-3);
  themap["anumu"]->SetFillColor(kOrange-3);
  hs_trklen->Add(themap["anumu"]);
  themap["nue"]->SetLineColor(kMagenta+1);
  themap["nue"]->SetFillColor(kMagenta+1);
  hs_trklen->Add(themap["nue"]);
  
  if (_breackdownPlots) {
    themap["signal_nostopmu"]->SetLineColor(kRed+2);
    themap["signal_nostopmu"]->SetFillColor(kRed+2);
    hs_trklen->Add(themap["signal_nostopmu"]);
    themap["signal_stopmu"]->SetLineColor(kRed);
    themap["signal_stopmu"]->SetFillColor(kRed);
    hs_trklen->Add(themap["signal_stopmu"]);
  }
  else {
    themap["signal"]->SetLineColor(kRed);
    themap["signal"]->SetFillColor(kRed);
    hs_trklen->Add(themap["signal"]);
  }
  hs_trklen->Draw();
  
  //h_trklen_total->DrawCopy("hist");
  
  //gStyle->SetHatchesLineWidth(1);
  themap["total"]->SetFillColor(kBlack);
  themap["total"]->SetFillStyle(3005);
  themap["total"]->Draw("E2 same");

}










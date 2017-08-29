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

#include "AnaTree.h"
#include "Spectrum.hpp"
#include "Spectrum2D.hpp"
//#include "PlotHandler.hpp"
//#include "SelectionTools.hpp"

const double _beamSpillStarts = 3.2; // us
const double _beamSpillEnds   = 4.8; // us

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
  
  TLatex* prelim = new TLatex(.35, .95, str.c_str());
  prelim->SetTextColor(kGray+2);
  prelim->SetNDC();
  prelim->SetTextSize(1/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

//____________________________________________________________________________________________________
void ActivateBranches(AnaTree *at) {
  
  //at->fChain->SetBranchStatus("",1);
}


double CalcLength(const double& x_1, const double& y_1, const double& z_1, const double& x_2, const double& y_2, const double& z_2) {
  return sqrt(pow(x_1-x_2, 2) + pow(y_1-y_2, 2) + pow(z_1-z_2, 2));
}

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
  
  
  TH1D* h_trklen_total = new TH1D("h_trklen_total", "; Track length;", 30, 0, 700);
  TH1D* h_trklen_signal = new TH1D("h_trklen_signal", "; Track length;", 30, 0, 700);
  TH1D* h_trklen_cosmic = new TH1D("h_trklen_cosmic", "; Track length;", 30, 0, 700);
  TH1D* h_trklen_outfv = new TH1D("h_trklen_outfv", "; Track length;", 30, 0, 700);
  TH1D* h_trklen_nc = new TH1D("h_trklen_nc", "; Track length;", 30, 0, 700);
  TH1D* h_trklen_anumu = new TH1D("h_trklen_anumu", "; Track length;", 30, 0, 700);
  TH1D* h_trklen_nue = new TH1D("h_trklen_nue", "; Track length;", 30, 0, 700);
  TH1D* h_trklen_cosmic_stopmu = new TH1D("h_trklen_cosmic_stopmu", "; Track length;", 30, 0, 700);
  TH1D* h_trklen_cosmic_nostopmu = new TH1D("h_trklen_cosmic_nostopmu", "; Track length;", 30, 0, 700);
  
  TH1D* h_trkphi_total = new TH1D("h_trkphi_total", "; Track #phi;", 20, -3.15, 3.15);
  TH1D* h_trkphi_signal = new TH1D("h_trkphi_signal", "; Track #phi;", 20, -3.15, 3.15);
  TH1D* h_trkphi_cosmic = new TH1D("h_trkphi_cosmic", "; Track #phi;", 20, -3.15, 3.15);
  TH1D* h_trkphi_outfv = new TH1D("h_trkphi_outfv", "; Track #phi;", 20, -3.15, 3.15);
  TH1D* h_trkphi_nc = new TH1D("h_trkphi_nc", "; Track #phi;", 20, -3.15, 3.15);
  TH1D* h_trkphi_anumu = new TH1D("h_trkphi_anumu", "; Track #phi;", 20, -3.15, 3.15);
  TH1D* h_trkphi_nue = new TH1D("h_trkphi_nue", "; Track #phi;", 20, -3.15, 3.15);
  TH1D* h_trkphi_cosmic_stopmu = new TH1D("h_trkphi_cosmic_stopmu", "; Track #phi;", 20, -3.15, 3.15);
  TH1D* h_trkphi_cosmic_nostopmu = new TH1D("h_trkphi_cosmic_nostopmu", "; Track #phi;", 20, -3.15, 3.15);
  
  TH1D* h_trktheta_total = new TH1D("h_trktheta_total", "; Track cos(#theta);", 30, -1, 1);
  TH1D* h_trktheta_signal = new TH1D("h_trktheta_signal", "; Track cos(#theta);", 30, -1, 1);
  TH1D* h_trktheta_cosmic = new TH1D("h_trktheta_cosmic", "; Track cos(#theta);", 30, -1, 1);
  TH1D* h_trktheta_outfv = new TH1D("h_trktheta_outfv", "; Track cos(#theta);", 30, -1, 1);
  TH1D* h_trktheta_nc = new TH1D("h_trktheta_nc", "; Track cos(#theta);", 30, -1, 1);
  TH1D* h_trktheta_anumu = new TH1D("h_trktheta_anumu", "; Track cos(#theta);", 30, -1, 1);
  TH1D* h_trktheta_nue = new TH1D("h_trktheta_nue", "; Track cos(#theta);", 30, -1, 1);
  TH1D* h_trktheta_cosmic_stopmu = new TH1D("h_trktheta_cosmic_stopmu", "; Track cos(#theta);", 30, -1, 1);
  TH1D* h_trktheta_cosmic_nostopmu = new TH1D("h_trktheta_cosmic_nostopmu", "; Track cos(#theta);", 30, -1, 1);
  
  TH1D* h_multpfp_total = new TH1D("h_multpfp_total", "; PFP Multiplicity", 10, 0, 10);
  TH1D* h_multpfp_signal = new TH1D("h_multpfp_signal", "; PFP Multiplicity;", 10, 0, 10);
  TH1D* h_multpfp_cosmic = new TH1D("h_multpfp_cosmic", "; PFP Multiplicity;", 10, 0, 10);
  TH1D* h_multpfp_outfv = new TH1D("h_multpfp_outfv", "; PFP Multiplicity;", 10, 0, 10);
  TH1D* h_multpfp_nc = new TH1D("h_multpfp_nc", "; PFP Multiplicity;", 10, 0, 10);
  TH1D* h_multpfp_anumu = new TH1D("h_multpfp_anumu", "; PFP Multiplicity;", 10, 0, 10);
  TH1D* h_multpfp_nue = new TH1D("h_multpfp_nue", "; PFP Multiplicity;", 10, 0, 10);
  TH1D* h_multpfp_cosmic_stopmu = new TH1D("h_multpfp_cosmic_stopmu", "; PFP Multiplicity;", 10, 0, 10);
  TH1D* h_multpfp_cosmic_nostopmu = new TH1D("h_multpfp_cosmic_nostopmu", "; PFP Multiplicity;", 10, 0, 10);
  
  TH1D* h_multtracktol_total = new TH1D("h_multtracktol_total", "; Track Multiplicity (5 cm)", 10, 0, 10);
  TH1D* h_multtracktol_signal = new TH1D("h_multtracktol_signal", "; Track Multiplicity (5 cm);", 10, 0, 10);
  TH1D* h_multtracktol_cosmic = new TH1D("h_multtracktol_cosmic", "; Track Multiplicity (5 cm);", 10, 0, 10);
  TH1D* h_multtracktol_outfv = new TH1D("h_multtracktol_outfv", "; Track Multiplicity (5 cm);", 10, 0, 10);
  TH1D* h_multtracktol_nc = new TH1D("h_multtracktol_nc", "; Track Multiplicity (5 cm);", 10, 0, 10);
  TH1D* h_multtracktol_anumu = new TH1D("h_multtracktol_anumu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  TH1D* h_multtracktol_nue = new TH1D("h_multtracktol_nue", "; Track Multiplicity (5 cm);", 10, 0, 10);
  TH1D* h_multtracktol_cosmic_stopmu = new TH1D("h_multtracktol_cosmic_stopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  TH1D* h_multtracktol_cosmic_nostopmu = new TH1D("h_multtracktol_cosmic_nostopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  
  
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
      if (at->event == 150801) {
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
      
      // CheckVertex
      if (at->slc_origin->at(slc) == 0 && at->fv == 1 && at->slc_vtxcheck_angle->at(slc) > -9999) {
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
    
    if(at->slc_ntrack->at(scl_ll_max) == 0) continue;
    
    if(!at->slc_passed_min_track_quality->at(scl_ll_max)) continue;
    
    //if(at->slc_ntrack->at(scl_ll_max) == 1 && at->slc_crosses_top_boundary->at(scl_ll_max) == 1) continue;
    
    //if(at->slc_longesttrack_length->at(scl_ll_max) < 25.) continue;
    
    //if(!at->slc_iscontained->at(scl_ll_max)) continue;
    
    //if(at->slc_crosses_top_boundary->at(scl_ll_max) == 1) continue;
    
    
    
    // Event is selected
    
    
    h_trklen_total->Fill(at->slc_longesttrack_length->at(scl_ll_max));
    h_trkphi_total->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
    h_trktheta_total->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
    h_multpfp_total->Fill(at->slc_mult_pfp->at(scl_ll_max));
    h_multtracktol_total->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
    
    bool nu_origin = false;
    if ((at->slc_origin->at(scl_ll_max) == 0 || at->slc_origin->at(scl_ll_max) == 2)) nu_origin = true;
    
    // Signal
    if(nu_origin && at->ccnc==0 && at->nupdg==14 && at->fv==1){
      signal_sel ++;
      h_eff_num->Fill(at->nu_e);
      pEff->Fill(true, at->nu_e);
      h_trklen_signal->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      h_trkphi_signal->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      h_trktheta_signal->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      h_multpfp_signal->Fill(at->slc_mult_pfp->at(scl_ll_max));
      h_multtracktol_signal->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      //std::cout << "Is signal and is selected. event: " << at->event << std::endl;
      //if (at->slc_mult_track_tolerance->at(scl_ll_max) == 0) std::cout << "Is signal with mult_track_tolerance=0. event: " << at->event << std::endl;
    }
    // anumu
    else if(nu_origin && at->ccnc==0 && at->nupdg==-14 && at->fv==1){
      bkg_anumu_sel ++;
      pEff->Fill(false, at->nu_e);
      h_trklen_anumu->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      h_trkphi_anumu->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      h_trktheta_anumu->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      h_multpfp_anumu->Fill(at->slc_mult_pfp->at(scl_ll_max));
      h_multtracktol_anumu->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
    }
    // nue
    else if(nu_origin && at->ccnc==0 && (at->nupdg==-12 || at->nupdg==12) && at->fv==1){
      bkg_nue_sel ++;
      pEff->Fill(false, at->nu_e);
      h_trklen_nue->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      h_trkphi_nue->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      h_trktheta_nue->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      h_multpfp_nue->Fill(at->slc_mult_pfp->at(scl_ll_max));
      h_multtracktol_nue->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
    }
    // nc
    else if(nu_origin && at->ccnc==1 && at->fv==1){
      bkg_nc_sel ++;
      pEff->Fill(false, at->nu_e);
      h_trklen_nc->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      h_trkphi_nc->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      h_trktheta_nc->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      h_multpfp_nc->Fill(at->slc_mult_pfp->at(scl_ll_max));
      h_multtracktol_nc->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      //std::cout << "Is a nc but is selected. event: " << at->event << std::endl;
    }
    // outfv
    else if(nu_origin && at->fv==0){
      bkg_outfv_sel ++;
      //std::cout << "Is OutFV. event: " << event << std::endl;
      pEff->Fill(false, at->nu_e);
      h_trklen_outfv->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      h_trkphi_outfv->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      h_trktheta_outfv->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      h_multpfp_outfv->Fill(at->slc_mult_pfp->at(scl_ll_max));
      h_multtracktol_outfv->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      //if (at->slc_mult_track_tolerance->at(scl_ll_max) == 2) std::cout << "Is OutFV with mult_track_tolerance=2. event: " << at->event << std::endl;
    }
    // cosmic
    else{
      bkg_cosmic_sel ++;
      if (at->slc_crosses_top_boundary->at(scl_ll_max) == 1 )
        bkg_cosmic_top_sel++;
      pEff->Fill(false, at->nu_e);
      //if (at->muon_is_reco == 1)
      h_trklen_cosmic->Fill(at->slc_longesttrack_length->at(scl_ll_max));
      h_trkphi_cosmic->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
      h_trktheta_cosmic->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
      h_multpfp_cosmic->Fill(at->slc_mult_pfp->at(scl_ll_max));
      h_multtracktol_cosmic->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      //std::cout << "Is a cosmic but is selected. event: " << at->event << std::endl;
      
      if (at->slc_origin_extra->at(scl_ll_max) == 0) {
        h_trklen_cosmic_stopmu->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        h_trkphi_cosmic_stopmu->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        h_trktheta_cosmic_stopmu->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        h_multpfp_cosmic_stopmu->Fill(at->slc_mult_pfp->at(scl_ll_max));
        h_multtracktol_cosmic_stopmu->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
      } else {
        h_trklen_cosmic_nostopmu->Fill(at->slc_longesttrack_length->at(scl_ll_max));
        h_trkphi_cosmic_nostopmu->Fill(at->slc_longesttrack_phi->at(scl_ll_max));
        h_trktheta_cosmic_nostopmu->Fill(at->slc_longesttrack_theta->at(scl_ll_max));
        h_multpfp_cosmic_nostopmu->Fill(at->slc_mult_pfp->at(scl_ll_max));
        h_multtracktol_cosmic_nostopmu->Fill(at->slc_mult_track_tolerance->at(scl_ll_max));
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
  std::cout << "Number of neutrino origin slices tagged as cosmic by the ACPT algo in total: " << n_slc_acpt_tag_nu << std::endl;
  
  
  
  
  
  
  
  // ************************
  //
  //  Plotting
  //
  // ************************
  
  new TCanvas();
  TEfficiency* pEff2 = new TEfficiency(*h_eff_num,*h_eff_den);
  pEff2->SetTitle(";True Neutrino Energy [GeV];Efficiency");
  pEff2->Draw("AP");
  
  new TCanvas();
  TEfficiency* pEff3 = new TEfficiency(*h_mueff_num,*h_mueff_den);
  pEff3->SetTitle(";True Muon Momentum [GeV];Reconstruction Efficiency");
  pEff3->Draw("AP");
  
  //TCanvas *c33 = new TCanvas();
  TEfficiency* pEff4 = new TEfficiency(*h_mueff_2_num,*h_mueff_den);
  pEff4->SetTitle(";True Muon Momentum [GeV];Reconstruction Efficiency");
  pEff4->SetLineColor(kRed);
  pEff4->Draw("P same");
  
  
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
  h_trklen_cosmic_nostopmu->SetLineColor(kBlue+2);
  h_trklen_cosmic_nostopmu->SetFillColor(kBlue+2);
  hs_trklen->Add(h_trklen_cosmic_nostopmu);
  h_trklen_cosmic_stopmu->SetLineColor(kBlue);
  h_trklen_cosmic_stopmu->SetFillColor(kBlue);
  hs_trklen->Add(h_trklen_cosmic_stopmu);
  h_trklen_outfv->SetLineColor(kOrange+3);
  h_trklen_outfv->SetFillColor(kOrange+3);
  hs_trklen->Add(h_trklen_outfv);
  h_trklen_nc->SetLineColor(kGray);
  h_trklen_nc->SetFillColor(kGray);
  hs_trklen->Add(h_trklen_nc);
  h_trklen_anumu->SetLineColor(kOrange-3);
  h_trklen_anumu->SetFillColor(kOrange-3);
  hs_trklen->Add(h_trklen_anumu);
  h_trklen_nue->SetLineColor(kGreen+3);
  h_trklen_nue->SetFillColor(kGreen+3);
  hs_trklen->Add(h_trklen_nue);
  h_trklen_signal->SetLineColor(kRed);
  h_trklen_signal->SetFillColor(kRed);
  hs_trklen->Add(h_trklen_signal);
  hs_trklen->Draw();
  h_trklen_total->Draw("E1 X0 same");
  
  TLegend* leg2 = new TLegend(0.58,0.51,0.82,0.82,NULL,"brNDC");
  leg2->AddEntry(h_trklen_signal,"#nu_{#mu} CC (signal)","f");
  leg2->AddEntry(h_trklen_nue,"#nu_{e}, #bar{#nu}_{e} CC","f");
  leg2->AddEntry(h_trklen_anumu,"#bar{#nu}_{#mu} CC","f");
  leg2->AddEntry(h_trklen_nc,"NC","f");
  leg2->AddEntry(h_trklen_outfv,"OUTFV","f");
  leg2->AddEntry(h_trklen_cosmic_stopmu,"Cosmic (stopping #mu)","f");
  leg2->AddEntry(h_trklen_cosmic_nostopmu,"Cosmic (other)","f");
  leg2->AddEntry(h_trklen_total,"MC Stat Unc.","e");
  leg2->Draw();
  DrawPOT(totalPOT);
  
  TString temp2 = "./output/trklen";
  final1->SaveAs(temp2 + ".pdf");
  final1->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final2 = new TCanvas();
  THStack *hs_trkphi = new THStack("hs_trkphi",";Candidate Track #phi; Selected Events");
  h_trkphi_cosmic_nostopmu->SetLineColor(kBlue+2);
  h_trkphi_cosmic_nostopmu->SetFillColor(kBlue+2);
  hs_trkphi->Add(h_trkphi_cosmic_nostopmu);
  h_trkphi_cosmic_stopmu->SetLineColor(kBlue);
  h_trkphi_cosmic_stopmu->SetFillColor(kBlue);
  hs_trkphi->Add(h_trkphi_cosmic_stopmu);
  h_trkphi_outfv->SetLineColor(kOrange+3);
  h_trkphi_outfv->SetFillColor(kOrange+3);
  hs_trkphi->Add(h_trkphi_outfv);
  h_trkphi_nc->SetLineColor(kGray);
  h_trkphi_nc->SetFillColor(kGray);
  hs_trkphi->Add(h_trkphi_nc);
  h_trkphi_anumu->SetLineColor(kOrange-3);
  h_trkphi_anumu->SetFillColor(kOrange-3);
  hs_trkphi->Add(h_trkphi_anumu);
  h_trkphi_nue->SetLineColor(kGreen+3);
  h_trkphi_nue->SetFillColor(kGreen+3);
  hs_trkphi->Add(h_trkphi_nue);
  h_trkphi_signal->SetLineColor(kRed);
  h_trkphi_signal->SetFillColor(kRed);
  hs_trkphi->Add(h_trkphi_signal);
  hs_trkphi->Draw();
  h_trkphi_total->Draw("E1 X0 same");
  
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/trkphi";
  final2->SaveAs(temp2 + ".pdf");
  final2->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final3 = new TCanvas();
  THStack *hs_trktheta = new THStack("hs_trktheta",";Candidate Track cos(#theta); Selected Events");
  h_trktheta_cosmic_nostopmu->SetLineColor(kBlue+2);
  h_trktheta_cosmic_nostopmu->SetFillColor(kBlue+2);
  hs_trktheta->Add(h_trktheta_cosmic_nostopmu);
  h_trktheta_cosmic_stopmu->SetLineColor(kBlue);
  h_trktheta_cosmic_stopmu->SetFillColor(kBlue);
  hs_trktheta->Add(h_trktheta_cosmic_stopmu);
  h_trktheta_outfv->SetLineColor(kOrange+3);
  h_trktheta_outfv->SetFillColor(kOrange+3);
  hs_trktheta->Add(h_trktheta_outfv);
  h_trktheta_nc->SetLineColor(kGray);
  h_trktheta_nc->SetFillColor(kGray);
  hs_trktheta->Add(h_trktheta_nc);
  h_trktheta_anumu->SetLineColor(kOrange-3);
  h_trktheta_anumu->SetFillColor(kOrange-3);
  hs_trktheta->Add(h_trktheta_anumu);
  h_trktheta_nue->SetLineColor(kGreen+3);
  h_trktheta_nue->SetFillColor(kGreen+3);
  hs_trktheta->Add(h_trktheta_nue);
  h_trktheta_signal->SetLineColor(kRed);
  h_trktheta_signal->SetFillColor(kRed);
  hs_trktheta->Add(h_trktheta_signal);
  hs_trktheta->Draw();
  h_trktheta_total->Draw("E1 X0 same");
  
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/trktheta";
  final3->SaveAs(temp2 + ".pdf");
  final3->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final4 = new TCanvas();
  THStack *hs_multpfp = new THStack("hs_multpfp",";PFP Multiplicity; Selected Events");
  h_multpfp_cosmic_nostopmu->SetLineColor(kBlue+2);
  h_multpfp_cosmic_nostopmu->SetFillColor(kBlue+2);
  hs_multpfp->Add(h_multpfp_cosmic_nostopmu);
  h_multpfp_cosmic_stopmu->SetLineColor(kBlue);
  h_multpfp_cosmic_stopmu->SetFillColor(kBlue);
  hs_multpfp->Add(h_multpfp_cosmic_stopmu);
  h_multpfp_outfv->SetLineColor(kOrange+3);
  h_multpfp_outfv->SetFillColor(kOrange+3);
  hs_multpfp->Add(h_multpfp_outfv);
  h_multpfp_nc->SetLineColor(kGray);
  h_multpfp_nc->SetFillColor(kGray);
  hs_multpfp->Add(h_multpfp_nc);
  h_multpfp_anumu->SetLineColor(kOrange-3);
  h_multpfp_anumu->SetFillColor(kOrange-3);
  hs_multpfp->Add(h_multpfp_anumu);
  h_multpfp_nue->SetLineColor(kGreen+3);
  h_multpfp_nue->SetFillColor(kGreen+3);
  hs_multpfp->Add(h_multpfp_nue);
  h_multpfp_signal->SetLineColor(kRed);
  h_multpfp_signal->SetFillColor(kRed);
  hs_multpfp->Add(h_multpfp_signal);
  hs_multpfp->Draw();
  h_multpfp_total->Draw("E1 X0 same");
  
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/multpfp";
  final4->SaveAs(temp2 + ".pdf");
  final4->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final5 = new TCanvas();
  THStack *hs_multtracktol = new THStack("hs_multtracktol",";Track Multiplicity (5 cm); Selected Events");
  h_multtracktol_cosmic_nostopmu->SetLineColor(kBlue+2);
  h_multtracktol_cosmic_nostopmu->SetFillColor(kBlue+2);
  hs_multtracktol->Add(h_multtracktol_cosmic_nostopmu);
  h_multtracktol_cosmic_stopmu->SetLineColor(kBlue);
  h_multtracktol_cosmic_stopmu->SetFillColor(kBlue);
  hs_multtracktol->Add(h_multtracktol_cosmic_stopmu);
  h_multtracktol_outfv->SetLineColor(kOrange+3);
  h_multtracktol_outfv->SetFillColor(kOrange+3);
  hs_multtracktol->Add(h_multtracktol_outfv);
  h_multtracktol_nc->SetLineColor(kGray);
  h_multtracktol_nc->SetFillColor(kGray);
  hs_multtracktol->Add(h_multtracktol_nc);
  h_multtracktol_anumu->SetLineColor(kOrange-3);
  h_multtracktol_anumu->SetFillColor(kOrange-3);
  hs_multtracktol->Add(h_multtracktol_anumu);
  h_multtracktol_nue->SetLineColor(kGreen+3);
  h_multtracktol_nue->SetFillColor(kGreen+3);
  hs_multtracktol->Add(h_multtracktol_nue);
  h_multtracktol_signal->SetLineColor(kRed);
  h_multtracktol_signal->SetFillColor(kRed);
  hs_multtracktol->Add(h_multtracktol_signal);
  hs_multtracktol->Draw();
  h_multtracktol_total->Draw("E1 X0 same");
  
  leg2->Draw();
  DrawPOT(totalPOT);
  
  temp2 = "./output/multtracktol";
  final5->SaveAs(temp2 + ".pdf");
  final5->SaveAs(temp2 + ".C","C");
  
  
  
  
  
  
  
  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  rootapp->Run();
  rootapp->Terminate(0);
  
  return 0;
}

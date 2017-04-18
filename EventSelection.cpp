#include <iostream>
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

#include "AnaTree.h"
#include "Spectrum.hpp"
#include "Spectrum2D.hpp"
//#include "PlotHandler.hpp"
//#include "SelectionTools.hpp"

using namespace std;


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
  
  string beam       = argv[1];
  string dopot      = argv[2];
  int    maxEntries = atoi(argv[3]);
  cout << endl << "Working with the " << beam << " beam." << endl;
  
  TApplication* rootapp = new TApplication("ROOT Application",&argc, argv);
  //gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS
  
  string pattern;
  if(beam == "numi") pattern = "/data/t2k/lar/uboone/prodgenie_numi_nu_uboone_MCC7/prodgenie_numi_nu_cosmic_uboone_merged_gen_g4_detsim_reco1_reco2_ana.root";
  //if(beam == "bnb")  pattern = "/pnfs/uboone/persistent/users/aschu/MC_BNB_Cosmic/prodgenie_bnb_nu_cosmic_uboone_v05_08_00_at.root";
  if(beam == "bnb")  pattern = "files/output9850.root";
  
  
  bool evalPOT = false;
  double totalPOT = 0.;
  /*
   if (dopot == "pot") evalPOT = true;
   if (maxEntries > 0) evalPOT = false;
   
   if (evalPOT) {
   
   cout << " ----- " << endl;
   cout << "| Calculating simulated POT." << endl;
   cout << "| " << endl;
   
   TChain *cpot;
   cpot = new TChain("analysistree/pottree");
   cpot->Add(pattern.c_str());
   cout << "| Number of entries in the pot tree: " << cpot->GetEntries() << endl;
   Double_t pot;
   cpot->SetBranchAddress("pot", &pot);
   for (int potEntry = 0; potEntry < cpot->GetEntries(); potEntry++) {
   cpot->GetEntry(potEntry);
   totalPOT += pot;
   } // end loop entries
   cout << "| Simulated POT: " << totalPOT << endl;
   cout << " ----- " << endl << endl;
   } // end if evalPOT
   else totalPOT = -1.;
   */
  
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
  
  int selectedEvents = 0.;
  int counter = 0.;
  int EventsWithFlash = 0, EventsVtxInFV = 0, EventsFlashMatched = 0, EventsTracksInFV = 0, EventsTrackLong = 0;
  int noHitsOnUplane = 0, noHitsOnVplane = 0, noHitsOnWplane = 0;
  
  Spectrum* Sflashtime      = new Spectrum("flash_time",      ";Flash Time [#mus];Entries per bin",       300000, -3000, 3000, totalPOT);
  
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
  
  TH1D* h_eff_num = new TH1D("h_eff_num", "h_eff_num", 6, 0, 4);
  TH1D* h_eff_den = new TH1D("h_eff_den", "h_eff_den", 6, 0, 4);
  TEfficiency* pEff = new TEfficiency("eff",";Neutrino Energy (truth) [GeV];Efficiency",6, 0, 4);
  
  TH1D* h_chi2 = new TH1D("h_chi2", "h_chi2", 50, 0, 50);
  TH1D* h_flsTime = new TH1D("h_flsTime", ";Flash time w.r.t. trigger [#mus];Events", 100, 0, 25);
  TH1D* h_nslices = new TH1D("h_nslices", ";Number of slices per event;Entries per bin", 15, 0, 15);
  TH1D* h_vtx_resolution = new TH1D("h_nslh_vtx_resolutionices", ";Vertex resolution (2D) [cm];Entries per bin", 300, 0, 500);
  
  TH2D* h_frac_diff = new TH2D("h_frac_diff", ";PMT ID; Fractional difference", 32, 0, 32, 80, -2, 2);
  TH2D* h_frac_diff_others = new TH2D("h_frac_diff", ";PMT ID; Fractional difference", 32, 0, 32, 80, -2, 2);
  double hypo_spec_x[32], hypo_spec_y[32];
  double meas_spec_x[32], meas_spec_y[32];
  
  TH1D* h_vtxcheck_angle_good = new TH1D("h_vtxcheck_angle_good", ";Angle [rad];Entries per bin", 100, 0, 4);
  TH1D* h_vtxcheck_angle_bad  = new TH1D("h_vtxcheck_angle_bad",  ";Angle [rad];Entries per bin",  100, 0, 4);
  
  
  
  
  
  
  if(maxEntries > 0.) evts = maxEntries;
  
  for(int i = 0; i < evts; i++) {
    
    if(i%100 == 0) cout << "\t... " << i << endl;
    
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
    
    
    // ************************
    //
    //  Selection
    //
    // ************************
    
    
    
    
    
    bool isSignal = false;
    if (at->nupdg == 14 && at->ccnc == 0 && at->fv == 1){
      nsignal++;
      isSignal = true;
      h_eff_den->Fill(at->nu_e);
      
      if (at->muon_is_reco){
        nSignalWMuonReco++;
        if (at->vtx_resolution > -1 && at->vtx_resolution < 10)
          nSignalMuonRecoVtxOk++;
      }
      else{
        //std::cout << "This is a signal event but the muon was not reconstructed. Event: " << event << std::endl;
      }
    }
    if(at->nupdg == 14 && at->ccnc == 0){
      nNumuCC++;
    }
    
    
    if (isSignal) h_vtx_resolution->Fill(at->vtx_resolution);
    //if (isSignal && vtx_resolution > 200 && vtx_resolution < 210) std::cout << "vtx_resolution is fucked for event: " << event << std::endl;
    
    
    
    // Selection
    
    
    if (at->nbeamfls == 0) continue;
    int flashInBeamSpill = -1;
    
    bool goodflash = false;
    for (int fls = 0; fls < at->nbeamfls; fls ++){
      h_flsTime->Fill(at->beamfls_time->at(fls));
      if (at->beamfls_time->at(fls) > 3 && at->beamfls_time->at(fls) < 5) {
        
        flashInBeamSpill = fls;
        if (at->beamfls_pe->at(fls) >= 50) {
          goodflash = true;
          nEvtsWFlashInBeamSpill++;
        }
      }
    }
    
    
    if (!goodflash) continue;
    
    std::vector<bool> isBackground(at->nslices, false);
    
    // Slice loop
    for (int slc = 0; slc < at->nslices; slc ++) {
      
      // PMTs
      //   nu
      if (at->slc_origin->at(slc) == 0 && flashInBeamSpill > -1 && at->slc_flsmatch_score->at(slc) > -1){
        for (int pmt = 0; pmt < 32; pmt++) {
          if ((at->slc_flshypo_spec->at(slc))[pmt] < 5 || (at->beamfls_spec->at(flashInBeamSpill))[pmt] < 5) continue;
          double mean = ((at->slc_flshypo_spec->at(slc))[pmt] + (at->beamfls_spec->at(flashInBeamSpill))[pmt]) / 2.;
          h_frac_diff->Fill(pmt, ( (at->slc_flshypo_spec->at(slc))[pmt] - (at->beamfls_spec->at(flashInBeamSpill))[pmt] ) / (mean) );
        }
      }
      //   others
      if (at->slc_origin->at(slc) == 1 && flashInBeamSpill > -1 && at->slc_flsmatch_score->at(slc) > -1){
        for (int pmt = 0; pmt < 32; pmt++) {
          if ((at->slc_flshypo_spec->at(slc))[pmt] < 5 || (at->beamfls_spec->at(flashInBeamSpill))[pmt] < 5) continue;
          double mean = ((at->slc_flshypo_spec->at(slc))[pmt] + (at->beamfls_spec->at(flashInBeamSpill))[pmt]) / 2.;
          h_frac_diff_others->Fill(pmt, ( (at->slc_flshypo_spec->at(slc))[pmt] - (at->beamfls_spec->at(flashInBeamSpill))[pmt] ) / (mean) );
        }
      }
      //  spec
      if (i == 0) {
        if (flashInBeamSpill > -1 && at->slc_flsmatch_score->at(slc) > -1){
          for (int pmt = 0; pmt < 32; pmt++) {
            hypo_spec_x[pmt] = pmt;
            hypo_spec_y[pmt] = (at->slc_flshypo_spec->at(3))[pmt];
            meas_spec_x[pmt] = pmt;
            meas_spec_y[pmt] = (at->beamfls_spec->at(flashInBeamSpill))[pmt];
          }
        }
      }
      
      // CheckVertex
      if (at->slc_origin->at(slc) == 0 && at->fv == 1 && at->slc_vtxcheck_angle->at(slc) > -9999) {
        if (at->vtx_resolution <= 10.){
          h_vtxcheck_angle_good->Fill(at->slc_vtxcheck_angle->at(slc));
          if (at->slc_vtxcheck_angle->at(slc) > 3)
            std::cout << "Angle is about 180 for a good vertex for event: " << at->event << std::endl;
        }
        if (at->vtx_resolution > 10.){
          h_vtxcheck_angle_bad->Fill(at->slc_vtxcheck_angle->at(slc));
          if (at->slc_vtxcheck_angle->at(slc) < 2)
            std::cout << "Angle is not 180 for a bad vertex for event: " << at->event << std::endl;
        }
      }
      
      h_chi2->Fill(at->slc_kalman_chi2->at(slc)/(double)at->slc_kalman_ndof->at(slc));
      if (at->slc_kalman_chi2->at(slc)/(double)at->slc_kalman_ndof->at(slc) > 15)
        //std::cout << "This event (" << event << ") has a track with chi2/ndof = " << slc_kalman_chi2->at(slc)/(double)slc_kalman_ndof->at(slc) << std::endl;
        h_nslices->Fill(at->nslices);
      
      isBackground.at(slc) = false;
      
      
      // ACPT
      if (at->slc_acpt_outoftime->at(slc) == 1) {
        //isBackground.at(slc) = true;
        continue;
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
    
    // Find slice with maximum score that was not tagged as bkg
    double score_max = -1;
    int scl_ll_max = -1;
    for (int slc = 0; slc < at->nslices; slc ++){
      
      //if (isBackground.at(slc)) continue;
      
      if(at->slc_flsmatch_score->at(slc) > score_max){
        scl_ll_max = slc;
        score_max = at->slc_flsmatch_score->at(slc);
      }
    }
    
    
    if (scl_ll_max == -1) continue;
    
    if (score_max <= 0.00000001) continue;
    
    if (isSignal && at->slc_origin->at(scl_ll_max)==0) nSignalFlashMatched ++;
    //else if (isSignal) std::cout << "***This event " << event << " was incorrecly flash-matched." << std::endl;
    
    if(isBackground.at(scl_ll_max)) continue;
    
    if(at->slc_nuvtx_fv->at(scl_ll_max) == 0) continue;
    
    if(at->slc_origin->at(scl_ll_max) == 0 && at->ccnc==0 && at->nupdg==14 && at->fv==1){
      signal_sel ++;
      h_eff_num->Fill(at->nu_e);
      pEff->Fill(true, at->nu_e);
    }
    else if(at->slc_origin->at(scl_ll_max) == 0 && at->ccnc==0 && at->nupdg==-14 && at->fv==1){
      bkg_anumu_sel ++;
      pEff->Fill(false, at->nu_e);
    }
    else if(at->slc_origin->at(scl_ll_max) == 0 && at->ccnc==0 && (at->nupdg==-12 || at->nupdg==12) && at->fv==1){
      bkg_nue_sel ++;
      pEff->Fill(false, at->nu_e);
    }
    else if(at->slc_origin->at(scl_ll_max) == 0 && at->ccnc==1 && at->fv==1){
      bkg_nc_sel ++;
      pEff->Fill(false, at->nu_e);
    }
    else if(at->slc_origin->at(scl_ll_max) == 0 && at->fv==0){
      bkg_outfv_sel ++;
      //std::cout << "Is OutFV. event: " << event << std::endl;
      pEff->Fill(false, at->nu_e);
    }
    else{
      bkg_cosmic_sel ++;
      if (at->slc_crosses_top_boundary->at(scl_ll_max) == 1 )
        bkg_cosmic_top_sel++;
      pEff->Fill(false, at->nu_e);
    }
    
    
    
    

    
    
    
    
    
    
    
    
    //newtree->Fill();
    
    //delete selection;
  } // end of event loop
  
  //Sflashtime      ->Save();
  
  
  //newtree->AutoSave();
  //delete oldfile;
  //delete newfile;
  
  cout << endl << endl << "********************************" << endl;
  
  
  
  
  
  
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
  
  std::cout << "Number of signal events that were correctly flash-matched: " << nSignalFlashMatched << std::endl;
  
  
  
  
  
  //h_eff_num->Divide(h_eff_den);
  //h_eff_num->Draw();
  
  
  
  TCanvas *c3 = new TCanvas();
  TEfficiency* pEff2 = new TEfficiency(*h_eff_num,*h_eff_den);
  pEff2->SetTitle(";True Neutrino Energy [GeV];Efficiency");
  pEff2->Draw("AP");
  
  TCanvas *c2 = new TCanvas();
  h_chi2->Draw("histo");
  
  TCanvas *c4 = new TCanvas();
  h_flsTime->Draw("histo");
  
  TCanvas *c5 = new TCanvas();
  h_nslices->Draw("histo");
  
  TCanvas *c6 = new TCanvas();
  h_vtx_resolution->Draw("histo");
  
  TCanvas *c7 = new TCanvas();
  h_frac_diff->Draw("colz");
  
  TCanvas *c8 = new TCanvas();
  h_frac_diff_others->Draw("colz");
  
  // PE spec
  TCanvas *c9 = new TCanvas();
  TGraph* gr = new TGraph(32,hypo_spec_x,hypo_spec_y);
  TGraph* gr2 = new TGraph(32,meas_spec_x,meas_spec_y);
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
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(gr,"Hypo flash","l");
  leg->AddEntry(gr2,"Reco flash","l");
  leg->Draw();
  
  TCanvas *c10 = new TCanvas();
  h_vtxcheck_angle_good->Draw("histo");
  h_vtxcheck_angle_bad->Draw("histo same");
  h_vtxcheck_angle_bad->SetLineColor(kRed);

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  rootapp->Run();
  rootapp->Terminate(0);
  
  return 0;
}

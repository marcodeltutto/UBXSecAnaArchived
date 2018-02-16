#define ShowerTree_cxx
#include "ShowerTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ShowerTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ShowerTree.C
//      root> ShowerTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TH2D * h_sh_theta_length_numu = new TH2D("h_sh_theta_length_numu", "#nu_{#mu};Shower #theta;Shower Length [cm]", 30, -1, 1, 30, 0, 400);
   TH2D * h_sh_theta_length_nue = new TH2D("h_sh_theta_length_nue", "#nu_{e};Shower #theta;Shower Length [cm]", 30, -1, 1, 30, 0, 400);

   TH2D * h_sh_theta_openangle_numu = new TH2D("h_sh_theta_openangle_numu", "#nu_{#mu};Shower #theta;Shower Opening Angle;", 30, -1, 1, 30, 0, 45);
   TH2D * h_sh_theta_openangle_nue = new TH2D("h_sh_theta_openangle_nue", "#nu_{e};Shower #theta;Shower Opening Angle;", 30, -1, 1, 30, 0, 45);

   TH2D * h_sh_theta_phi_numu = new TH2D("h_sh_theta_phi_numu", "#nu_{#mu};Shower #theta;Shower #phi", 30, -1, 1, 30, -3.1415, 3.1415);
   TH2D * h_sh_theta_phi_nue = new TH2D("h_sh_theta_phi_nue", "#nu_{e};Shower #theta;Shower #phi", 30, -1, 1, 30, -3.1415, 3.1415);

   TH2D * h_sh_theta_deltaz_numu = new TH2D("h_sh_theta_deltaz_numu", "#nu_{#mu};Shower #theta;|Shower Z Start - Flash Z Centre| [cm]", 30, -1, 1, 30, 0, 600);
   TH2D * h_sh_theta_deltaz_nue = new TH2D("h_sh_theta_deltaz_nue", "#nu_{e};Shower #theta;|Shower Z Start - Flash Z Centre| [cm]", 30, -1, 1, 30, 0, 600);

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "Entries: " << nentries << std::endl;

   int counter_n_muon_sel = 0;
   int counter_n_shower_if_muon = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (!(s_track_pdg==13||s_track_pdg==-13)) continue;

      counter_n_muon_sel++;

      if (s_shower_length==-1 || s_shower_theta==-1) continue;

      counter_n_shower_if_muon++;
      

      if (s_nupdg==12||s_nupdg==-12) {
      	h_sh_theta_length_nue->Fill(s_shower_theta, s_shower_length);
      	h_sh_theta_openangle_nue->Fill(s_shower_theta, s_shower_openangle/3.1415*180);
      	h_sh_theta_phi_nue->Fill(s_shower_theta, s_shower_phi);
      	h_sh_theta_deltaz_nue->Fill(s_shower_theta, std::abs(s_shower_startz - s_flash_z));
      }

      if (s_nupdg==14||s_nupdg==-14) {
      	h_sh_theta_length_numu->Fill(s_shower_theta, s_shower_length);
      	h_sh_theta_openangle_numu->Fill(s_shower_theta, s_shower_openangle/3.1415*180);
      	h_sh_theta_phi_numu->Fill(s_shower_theta, s_shower_phi);
      	h_sh_theta_deltaz_numu->Fill(s_shower_theta, std::abs(s_shower_startz - s_flash_z));
      }

   }

   std::cout << "counter_n_muon_sel: " << counter_n_muon_sel << std::endl;
   std::cout << "counter_n_shower_if_muon: " << counter_n_shower_if_muon << std::endl;


   TCanvas * c_nue_1 = new TCanvas();
   h_sh_theta_length_nue->Draw("colz");
   c_nue_1->SaveAs("output/sh_theta_length_nue.pdf");
   TCanvas * c_nue_2 = new TCanvas();
   h_sh_theta_openangle_nue->Draw("colz");
   c_nue_2->SaveAs("output/sh_theta_openangle_nue.pdf");
   TCanvas * c_nue_3 = new TCanvas();
   h_sh_theta_phi_nue->Draw("colz");
   c_nue_3->SaveAs("output/sh_theta_phi_nue.pdf");
   TCanvas * c_nue_4 = new TCanvas();
   h_sh_theta_deltaz_nue->Draw("colz");
   c_nue_4->SaveAs("output/sh_theta_deltaz_nue.pdf");

   TCanvas * c_numu_1 = new TCanvas();
   h_sh_theta_length_numu->Draw("colz");
   c_numu_1->SaveAs("output/sh_theta_length_numu.pdf");
   TCanvas * c_numu_2 = new TCanvas();
   h_sh_theta_openangle_numu->Draw("colz");
   c_numu_2->SaveAs("output/sh_theta_openangle_numu.pdf");
   TCanvas * c_numu_3 = new TCanvas();
   h_sh_theta_phi_numu->Draw("colz");
   c_numu_3->SaveAs("output/sh_theta_phi_numu.pdf");
   TCanvas * c_numu_4 = new TCanvas();
   h_sh_theta_deltaz_numu->Draw("colz");
   c_numu_4->SaveAs("output/sh_theta_deltaz_numu.pdf");
}

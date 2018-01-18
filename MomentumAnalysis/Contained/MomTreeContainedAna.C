#define MomTreeContainedAna_cxx
#include "MomTreeContainedAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

std::ofstream _txtfile;

void MomTreeContainedAna::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L MomTreeContainedAna.C
  //      root> MomTreeContainedAna t
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
  
  _txtfile.open ("output.txt");
  _txtfile << "name, mean, std" << std::endl;
  
  std::map<std::string, TH1D*> mom_map;
  TString title = ";MCS Muon Momentum [GeV];;";
  mom_map["0.0-0.1"] = new TH1D("mom_0.0-0.2", title, 20, 0, 0.5);
  mom_map["0.1-0.2"] = new TH1D("mom_0.1-0.2", title, 20, 0, 0.5);
  mom_map["0.2-0.3"] = new TH1D("mom_0.2-0.3", title, 20, 0, 0.5);
  mom_map["0.3-0.4"] = new TH1D("mom_0.3-0.4", title, 20, 0, 0.7);
  mom_map["0.4-0.5"] = new TH1D("mom_0.4-0.5", title, 20, 0, 0.9);
  mom_map["0.5-0.6"] = new TH1D("mom_0.5-0.6", title, 20, 0, 1.1);
  mom_map["0.6-0.7"] = new TH1D("mom_0.6-0.7", title, 20, 0, 1.3);
  mom_map["0.7-0.8"] = new TH1D("mom_0.7-0.8", title, 20, 0, 1.5);
  mom_map["0.8-0.9"] = new TH1D("mom_0.8-0.9", title, 20, 0, 1.7);
  mom_map["0.9-1.0"] = new TH1D("mom_0.9-1.0", title, 20, 0, 1.9);
  mom_map["1.0-1.1"] = new TH1D("mom_1.0-1.1", title, 20, 0, 2.1);
  mom_map["1.1-1.2"] = new TH1D("mom_1.1-1.2", title, 20, 0, 2.3);
  mom_map["1.2-1.3"] = new TH1D("mom_1.2-1.3", title, 20, 0, 2.5);
  mom_map["1.3-1.4"] = new TH1D("mom_1.3-1.4", title, 20, 0, 2.7);
  mom_map["1.4-1.5"] = new TH1D("mom_1.4-1.5", title, 20, 0, 2.9);

  TH2D* mom_true_reco = new TH2D("mom_true_reco", "mom_true_reco", 70, 0, 1, 70, 0, 1);
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    mom_true_reco->Fill(mom_true_contained, mom_mcs_contained);

    if (mom_true_contained > 0.0 && mom_true_contained < 0.1)
      mom_map["0.0-0.1"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.1 && mom_true_contained < 0.2)
      mom_map["0.1-0.2"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.2 && mom_true_contained < 0.3)
      mom_map["0.2-0.3"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.3 && mom_true_contained < 0.4)
      mom_map["0.3-0.4"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.4 && mom_true_contained < 0.5)
      mom_map["0.4-0.5"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.5 && mom_true_contained < 0.6)
      mom_map["0.5-0.6"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.6 && mom_true_contained < 0.7)
      mom_map["0.6-0.7"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.7 && mom_true_contained < 0.8)
      mom_map["0.7-0.8"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.8 && mom_true_contained < 0.9)
      mom_map["0.8-0.9"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 0.9 && mom_true_contained < 1.0)
      mom_map["0.9-1.0"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 1.0 && mom_true_contained < 1.1)
      mom_map["1.0-1.1"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 1.1 && mom_true_contained < 1.2)
      mom_map["1.1-1.2"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 1.2 && mom_true_contained < 1.3)
      mom_map["1.2-1.3"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 1.3 && mom_true_contained < 1.4)
      mom_map["1.3-1.4"]->Fill(mom_mcs_contained);
    
    if (mom_true_contained > 1.4 && mom_true_contained < 1.5)
      mom_map["1.4-1.5"]->Fill(mom_mcs_contained);
    

    // if (Cut(ientry) < 0) continue;
  }
  
  
  double points[17] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45};
  double error_x[17] = {0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05};
  double mean_array[17], sigma_array[17];
  double array[2];

  int counter = 0;
  for (auto m : mom_map) {
    DoFit(m.second, m.first, array);
    mean_array[counter]  = array[0];
    sigma_array[counter] = array[1];
    counter++;
  }
  

  TCanvas * c = new TCanvas;
  
  TGraphErrors* gr = new TGraphErrors(17,points,mean_array,error_x,sigma_array);
  gr->SetMarkerStyle(kFullCircle);
  gr->SetMarkerSize(1.0);
  gr->SetLineWidth(2.0);

  gStyle->SetPalette(kBird);
  mom_true_reco->Draw("colz");
  gr->Draw("P");
  
  gr->Fit("pol1", "", "", 0.1, 0.9);
  TF1 *myfunc = gr->GetFunction("pol1");
  
  TString n = "final";
  
  c->SaveAs("output/" + n + ".pdf");
  c->SaveAs("output/" + n + ".C");


}


void MomTreeContainedAna::DoFit(TH1D * the_histo, std::string name, double *array) {
  
  TCanvas * c = new TCanvas;
  
  the_histo->Fit("gaus");
  TF1 *myfunc = the_histo->GetFunction("gaus");
  the_histo->Draw("histo");
  myfunc->Draw("same");
  
  std::cout << "mean: " << myfunc->GetParameter(1) << std::endl;
  std::cout << "sigma: " << myfunc->GetParameter(2) << std::endl;
  array[0] = myfunc->GetParameter(1);
  array[1] = myfunc->GetParameter(2);
  
  std::string j = name + " GeV";
  
  TLatex* range = new TLatex(.8, .75, j.c_str());
  //range->SetTextColor(kBlue);
  range->SetNDC();
  //range->SetTextSize(2/30.);
  range->SetTextAlign(32);
  range->Draw();
  
  TString n = name;
  n = n + "_gaus";
  
  c->SaveAs("output/" + n + ".pdf");

  _txtfile << name << " GeV |" << array[0] << "|" << array[1] << std::endl;

}


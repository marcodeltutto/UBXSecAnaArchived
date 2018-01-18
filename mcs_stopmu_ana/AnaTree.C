#define AnaTree_cxx
#include "AnaTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnaTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L AnaTree.C
//      root> AnaTree t
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

	TH1D * delta_ll_cosmic_stopmu = new TH1D("delta_ll_cosmic_stopmu", "delta_ll_cosmic_stopmu", 100, -50, 50);
	TH1D * delta_ll_cosmic_nostopmu = new TH1D("delta_ll_cosmic_nostopmu", "delta_ll_cosmic_nostopmu", 100, -50, 50);
	TH1D * delta_ll_neutrino = new TH1D("delta_ll_neutrino", "delta_ll_neutrino", 100, -50, 50);

   //TH1D * score_h = new TH1D("score_h", "score_h", 20, -10, 10);
   //TH1D * signal_perc_h = new TH1D("signal_perc_h", "signal_perc_h", 20, -10, 10);

	Long64_t nentries = fChain->GetEntriesFast();

	double B = 0, S = 0, B_stopmu = 0;
	double total_S = 0, total_B = 0, total_B_stopmu = 0;
	double score = -1;
	std::vector<double> cut_v;
	std::vector<double> score_v;
	std::vector<double> signal_perc_v;
	std::vector<double> bkg_rej_v;
	std::vector<double> bkg_stopmu_rej_v;


	for (double cut_value = 5; cut_value >= -8; cut_value -=0.01) {

		B = 0;
		S = 0;
		B_stopmu = 0;
		total_S = 0;
		total_B = 0;
		total_B_stopmu = 0;

		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			if (delta_ll == -9999) continue;

			if (origin == 1)
				total_B++;

			if (origin == 1 && origin_extra == 0) {
				total_B_stopmu++;
				delta_ll_cosmic_stopmu->Fill(delta_ll);
			}
			if (origin == 1 && origin_extra != 0) {
				delta_ll_cosmic_nostopmu->Fill(delta_ll);
			}
			if (origin == 0) {
				total_S ++;
				delta_ll_neutrino->Fill(delta_ll);
			}


			if (delta_ll > cut_value) {
				if (origin == 0)
					S++;
				if (origin == 1)
					B++;
				if (origin == 1 && origin_extra == 0)
					B_stopmu++;
			}
      // if (Cut(ientry) < 0) continue;
		}

		score = S / sqrt(S + B);
		std::cout << "Cut at " << cut_value << " ==> Score is " << score  
		          << ",  signal perc is " << S/total_S  
		          << ",  bkg rej is " << 1-B/total_B 
		          << ",  bkg rej (stopmu) is " << 1-B_stopmu/total_B_stopmu << std::endl;

		cut_v.push_back(cut_value);
		score_v.push_back(score);
		signal_perc_v.push_back(S/total_S);
		bkg_rej_v.push_back(1-B/total_B);
		bkg_stopmu_rej_v.push_back(1-B_stopmu/total_B_stopmu);

   	//score_h-SetBinContent();
	}

	c2 = new TCanvas("c2","gerrors2",200,10,700,500);

	delta_ll_neutrino->SetFillColor(kRed+1);
	delta_ll_neutrino->SetLineColor(kRed+1);
	delta_ll_neutrino->SetFillStyle(3354);
	delta_ll_neutrino->DrawNormalized();
	delta_ll_cosmic_stopmu->SetLineColor(kBlue+1);
	delta_ll_cosmic_stopmu->SetLineWidth(3);
	delta_ll_cosmic_nostopmu->SetLineColor(kGreen+3);
	delta_ll_cosmic_nostopmu->SetLineWidth(3);
	delta_ll_cosmic_stopmu->DrawNormalized("same");
	delta_ll_cosmic_nostopmu->DrawNormalized("same");



  // Legend
  l = new TLegend(0.1,0.7,0.48,0.9);
  l->AddEntry(delta_ll_neutrino,"Neutrino Origin (vertex in FV)","f");
  l->AddEntry(delta_ll_cosmic_stopmu,"Cosmic Origin (stopping in TPC)","l");
  l->AddEntry(delta_ll_cosmic_nostopmu,"Cosmic Origin (through-going)","l");
  l->Draw();

  TLatex* prelim = new TLatex(.9, .95, "Area Normalised");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();


	double* x = &cut_v[0];
	double* y = &score_v[0];
	double* y_2 = &signal_perc_v[0];
	double* y_3 = &bkg_rej_v[0];
	double* y_4 = &bkg_stopmu_rej_v[0];

	int n = (int) score_v.size();



	// Signal efficiency
	TGraph* gr2 = new TGraph(n,x,y_2);
	gr2->SetLineColor(kRed+2);
	gr2->SetMarkerColor(kRed+2);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerSize(0.5);

	// Background rejection (all)
	TGraph* gr0 = new TGraph(n,x,y_3);
	gr0->SetLineColor(kGreen+3);
	gr0->SetMarkerColor(kGreen+3);
	gr0->SetMarkerStyle(20);
	gr0->SetMarkerSize(0.5);

	// Background rejection (stop mu)
	TGraph* gr5 = new TGraph(n,x,y_4);
	gr5->SetLineColor(kBlue+1);
	gr5->SetMarkerColor(kBlue+1);
	gr5->SetMarkerStyle(20);
	gr5->SetMarkerSize(0.5);

	// Score function
	TGraph* gr = new TGraph(n,x,y);
	gr->SetLineColor(kGreen+3);
	gr->SetMarkerColor(kGreen+3);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.5);



	// Plot efficinecy and rejection all together
	c3= new TCanvas("c3","c3",200,10,700,500);
  gr0->Draw("APL");
	gr0->GetXaxis()->SetRangeUser(-8,6);
  gr2->Draw("LP");
  gr5->Draw("LP");
  ll = new TLegend(0.1,0.7,0.48,0.9);
  ll->AddEntry(gr2,"Efficiency of Retaining Signal","lp");
  ll->AddEntry(gr0,"Background Rejection (all cosmics)","lp");
  ll->AddEntry(gr5,"Background Rejection (stopping muons)","lp");
  ll->Draw();





  // Plot with two Y axis

  //gROOT->ForceStyle(0);
	c1 = new TCanvas("c1","gerrors2",200,10,700,500);
	TPad *pad = new TPad("pad","",0,0,1,1);
  //pad->SetFillColor(42);
  //pad->SetGrid();
	pad->Draw();
	pad->cd();

  // draw a frame to define the range
	TH1F *hr = c1->DrawFrame(-0.4,0,1.2,12);
	hr->SetXTitle("X title");
	hr->SetYTitle("Y title");
	pad->GetFrame()->SetFillColor(21);
	pad->GetFrame()->SetBorderSize(12);

  // First graph
  /*
	TGraph* gr = new TGraph(n,x,y);
	gr->SetLineColor(kGreen+3);
	gr->SetMarkerColor(kGreen+3);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.5);
	gr->Draw("APL");
	gr->GetXaxis()->SetRangeUser(-8,6);


  // Second graph
	TGraph* gr2 = new TGraph(n,x,y_2);
	gr2->SetLineColor(kRed+2);
	gr2->SetMarkerColor(kRed+2);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerSize(0.5);
  //gr2->Draw("PL");
*/

  gr->GetXaxis()->SetRangeUser(-9,6);
	gr->Draw("APL");

  pad->Update();

  //create a transparent pad drawn on top of the main pad
	c1->cd();
	TPad *overlay = new TPad("overlay","",0,0,1,1);
	overlay->SetFillStyle(0);
	overlay->SetFillColor(0);
	overlay->SetFrameFillStyle(0);
	overlay->Draw("FA");
	overlay->cd();
  Double_t xmin = pad->GetUxmin();
  Double_t ymin = 0;
  Double_t xmax = pad->GetUxmax();
  Double_t ymax = 1;
  TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
  hframe->GetXaxis()->SetLabelOffset(99);
  hframe->GetYaxis()->SetLabelOffset(99);
  gr2->Draw("LP");


  //Draw an axis on the right side
  TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,506,"+L");
  axis->SetLineColor(kRed+1);
  axis->SetLabelColor(kRed+1);
  axis->SetTextFont(42);
  axis->SetLabelSize(.04);
  axis->SetLabelOffset(.005);
  axis->Draw();

  hframe->GetYaxis()->SetTickLength(0);
  hframe->GetXaxis()->SetTickLength(0);

  hframe->GetXaxis()->SetTitle("Cut Value on #DeltaLL");
  hframe->GetYaxis()->SetTitle("");
  hframe->SetTitle("");


  // Legend
  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(gr,"Score Function S/#sqrt{S+B}","lp");
  leg->AddEntry(gr2,"Efficiency of retaining signal","lp");
  leg->Draw();






/*
TCanvas *c1 = new TCanvas("c1","hists with different scales",600,400);


  TGraph* gr = new TGraph(n,x,y);
  gr->Draw("APL");
  gr->SetLineColor(kGreen+3);
  gr->SetMarkerColor(kGreen+3);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.5);

  TGraph* gr2 = new TGraph(n,x,y_2);
  gr2->Draw("PL");
  gr2->SetLineColor(kRed+2);
  gr2->SetMarkerColor(kRed+2);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.5);

  Float_t rightmax = 80;
  std::cout << "null->GetUxmax()" << null->GetUxmax()<< std::endl;
    std::cout << "null->GetUymin()" << null->GetUymin()<< std::endl;
  std::cout << "null->GetUxmin()" << null->GetUxmin()<< std::endl;
  std::cout << "null->GetUymax()" << null->GetUymax()<< std::endl;

  TGaxis *axis = new TGaxis(null->GetUxmax(),null->GetUymin(),
         null->GetUxmax(), null->GetUymax(),0,rightmax,510,"+L");
   axis->SetLineColor(kRed);
   axis->SetLabelColor(kRed);
   axis->Draw();
*/


 }

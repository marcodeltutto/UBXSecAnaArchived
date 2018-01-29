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

	TH1D * h_pe_bkg = new TH1D("h_pe_bkg", ";Reconstructed Flash PE;", 750, 0, 10000);
	TH1D * h_pe_sig = new TH1D("h_pe_sig", ";Reconstructed Flash PE;", 750, 0, 10000);

	TH1D * h_n_fls_beamspill = new TH1D("h_n_fls_beamspill", ";Number of Flash in the Beam Spill;", 10, 0, 10);

	Long64_t nentries = fChain->GetEntriesFast();

	double _beamSpillStarts = 3.2;
	double _beamSpillEnds = 4.8;

  std::vector<double> pe_cut_v;
	std::vector<double> signal_evts_v, bkg_events_v;
	std::vector<double> signal_evts_sel_v, bkg_events_sel_v;
	std::vector<double> eff_v, score_v;



	for (double pe_cut = 0; pe_cut < 30; pe_cut+=30) {

		int signal_evts = 0, bkg_evts = 0;
		int signal_evts_sel = 0, bkg_evts_sel = 0;

		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      //
      // Find flash in beam spill
      //

			int flashInBeamSpill = -1;
			double old_pe = -1;

			int flash_counter = 0;

			for (int fls = 0; fls < nbeamfls; fls ++) {
				if (beamfls_time.at(fls) > _beamSpillStarts && beamfls_time.at(fls) < _beamSpillEnds) {
					flash_counter++;
					if (beamfls_pe.at(fls) > old_pe) {
						flashInBeamSpill = fls;
						old_pe = beamfls_pe.at(fls);
					}
				}
			}


      //
      // Check if signal
      //

			bool isSignal = false;
			if (nupdg == 14 && ccnc == 0 && fv == 1)
				isSignal = true;


      //
      // Make PE plot and flash number
      //

      if (pe_cut == 0)
      	h_n_fls_beamspill->Fill(flash_counter);

			if (pe_cut == 0 && flashInBeamSpill != -1) {
				if (isSignal)
					h_pe_sig->Fill(beamfls_pe.at(flashInBeamSpill));
				else 
					h_pe_bkg->Fill(beamfls_pe.at(flashInBeamSpill));
			}



      //
      // Check truth signal and bkg events
      //

			if (isSignal)
				signal_evts++;
			else 
				bkg_evts++;



			//
			// Keep only events with flash in beam spill
			//

			if (flashInBeamSpill == -1) continue;


      //
      // Check PE cut
      //

			if (pe_cut == 0 && beamfls_pe.at(flashInBeamSpill) < 100 && isSignal)
				std::cout << "low pe for signal, pe = " << beamfls_pe.at(flashInBeamSpill) << std::endl;

			if (beamfls_pe.at(flashInBeamSpill) > pe_cut) {
				if (isSignal)
				  signal_evts_sel++;
				else 
					bkg_evts_sel++;
			}



		}

		std::cout << "pe_cut = " << pe_cut << " ==> signal_evts = " << signal_evts 
		          << ", signal_evts_sel = " << signal_evts_sel 
		          << ", eff = " << (double)signal_evts_sel/(double)signal_evts
		          << ", score = " << (double)signal_evts_sel/std::sqrt((double)signal_evts_sel + (double)bkg_evts_sel) << std::endl;
    
    pe_cut_v.emplace_back(pe_cut);
		signal_evts_v.emplace_back(signal_evts);
		signal_evts_sel_v.emplace_back(signal_evts_sel);
		signal_evts_v.emplace_back(bkg_evts);
		signal_evts_sel_v.emplace_back(bkg_evts_sel);
		eff_v.emplace_back((double)signal_evts_sel/(double)signal_evts);
		score_v.emplace_back((double)signal_evts_sel/std::sqrt((double)signal_evts_sel + (double)bkg_evts_sel));
	}

  //
  // TGraph of eff
  //

  double* x = &pe_cut_v[0];
  double* y = &eff_v[0];
  double* y2 = &score_v[0];
  int n = (int) eff_v.size();

	TGraph* g = new TGraph(n,x,y);
	g->SetLineColor(kRed+2);
  g->SetMarkerColor(kRed+2);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetTitle(";Flash PE Cut;Efficiency");
  TCanvas * c1 = new TCanvas();
  g->Draw("APL");

  TGraph* g2 = new TGraph(n,x,y2);
	g2->SetLineColor(kGreen+2);
  g2->SetMarkerColor(kGreen+2);
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(0.5);
  g2->SetTitle(";Flash PE Cut;S/#sqrt{S+B}");
  TCanvas * c2 = new TCanvas();
  g2->Draw("APL");


  //
  // Plot of PE
  //

  TCanvas * c3 = new TCanvas();
	h_pe_bkg->Draw();
	h_pe_sig->SetLineColor(kRed+2);
	h_pe_sig->Draw("same");
	TLegend* l = new TLegend(0.1,0.7,0.48,0.9);
  l->AddEntry(h_pe_sig,"#nu_{#mu} CC in FV","l");
  l->AddEntry(h_pe_bkg,"Background","l");
  l->Draw();
  /*TPad *subpad = new TPad("subpad","",0.6,0.3,0.95,0.65); 
  subpad->Draw(); 
  subpad->cd(); 
  h_pe_bkg->Draw();
	h_pe_sig->SetLineColor(kRed+2);
	h_pe_sig->Draw("same");*/

  //
  // Plot n flashes
  //

  TCanvas * c4 = new TCanvas();
	h_n_fls_beamspill->Draw();



  // Plot with two Y axis
  
  //gROOT->ForceStyle(0);
  c0 = new TCanvas("c0","canvas",200,10,700,500);
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
  
  g2->GetXaxis()->SetRangeUser(0,2000);
  g2->Draw("APL");

  g2->GetYaxis()->SetTickLength(0);
  g2->GetYaxis()->SetLabelSize(0);
  g2->GetXaxis()->SetTitle("");
  g2->GetYaxis()->SetTitle("");
  g2->SetTitle("");
  
  pad->Update();
  
  //create a transparent pad drawn on top of the main pad
  c0->cd();
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
  g->Draw("LP");
  
  
  //Draw an axis on the right side
  TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,506,"+L");
  axis->SetLineColor(kRed+2);
  axis->SetLabelColor(kRed+2);
  axis->SetTitleColor(kRed+2);  
  axis->SetTextFont(42);
  axis->SetLabelSize(.04);
  axis->SetLabelOffset(.005);
  axis->SetTitle("Efficiency");
  axis->SetTitleSize(.055);
  axis->SetTitleOffset(.8);
  axis->Draw();

  //Draw an axis on the left side
  ymin = pad->GetUymin();
  ymax = pad->GetUymax();
  std::cout << pad->GetUymin() << " " << pad->GetUymax() << std::endl;
  TGaxis *axisL = new TGaxis(xmin,0,xmin, 1,ymin,ymax,506,"");
  axisL->SetLineColor(kGreen+2);
  axisL->SetLabelColor(kGreen+2);
  axisL->SetTitleColor(kGreen+2);
  axisL->SetTextFont(42);
  axisL->SetLabelSize(.04);
  axisL->SetLabelOffset(.005);
  axisL->SetTitle("S/#sqrt{S+B}");
  axisL->SetTitleSize(.055);
  axisL->SetTitleOffset(.8);
  axisL->Draw();
  
  hframe->GetYaxis()->SetTickLength(0);
  hframe->GetXaxis()->SetTickLength(0);
  
  hframe->GetXaxis()->SetTitle("Flash PE Cut");
  hframe->GetYaxis()->SetTitle("");
  hframe->SetTitle("");
  
  
  // Legend
  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(g2, "Score Function S/#sqrt{S+B}","lp");
  leg->AddEntry(g, "Efficiency","lp");
  leg->Draw();
  
}

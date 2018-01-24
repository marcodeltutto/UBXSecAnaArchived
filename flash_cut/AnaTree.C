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

	TH1D * h_pe_bkg = new TH1D("h_pe_bkg", "h_pe_bkg", 750, 0, 10000);
	TH1D * h_pe_sig = new TH1D("h_pe_sig", "h_pe_sig", 750, 0, 10000);

	Long64_t nentries = fChain->GetEntriesFast();

	double _beamSpillStarts = 3.2;
	double _beamSpillEnds = 4.8;

  std::vector<double> pe_cut_v;
	std::vector<double> signal_evts_v;
	std::vector<double> signal_evts_sel_v;
	std::vector<double> eff_v;



	for (double pe_cut = 0; pe_cut < 300; pe_cut+=5) {

		int signal_evts = 0;
		int signal_evts_sel = 0;

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

			for (int fls = 0; fls < nbeamfls; fls ++) {
				if (beamfls_time.at(fls) > _beamSpillStarts && beamfls_time.at(fls) < _beamSpillEnds) {
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
      // Make PE plot
      //

			if (pe_cut == 0 && flashInBeamSpill != -1) {
				if (isSignal)
					h_pe_sig->Fill(beamfls_pe.at(flashInBeamSpill));
				else 
					h_pe_bkg->Fill(beamfls_pe.at(flashInBeamSpill));
			}



      //
      // Keep only signal events
      //

			if (!isSignal) continue;

			signal_evts++;



			//
			// Keep only events with flash in beam spill
			//

			if (flashInBeamSpill == -1) continue;


      //
      // Check PE cut
      //

			if (pe_cut == 0 && beamfls_pe.at(flashInBeamSpill) < 100)
				std::cout << "low pe for signal, pe = " << beamfls_pe.at(flashInBeamSpill) << std::endl;

			if (beamfls_pe.at(flashInBeamSpill) > pe_cut)
				signal_evts_sel ++;



		}

		std::cout << "pe_cut = " << pe_cut << " ==> signal_evts = " << signal_evts << ", signal_evts_sel = " << signal_evts_sel << ", eff = " << (double)signal_evts_sel/(double)signal_evts<< std::endl;
    pe_cut_v.emplace_back(pe_cut);
		signal_evts_v.emplace_back(signal_evts);
		signal_evts_sel_v.emplace_back(signal_evts_sel);
		eff_v.emplace_back((double)signal_evts_sel/(double)signal_evts);
	}

  //
  // TGraph of eff
  //

  double* x = &pe_cut_v[0];
  double* y = &eff_v[0];
  int n = (int) eff_v.size();

	TGraph* g = new TGraph(n,x,y);
	g->SetLineColor(kRed+2);
  g->SetMarkerColor(kRed+2);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  TCanvas * c1 = new TCanvas();
  g->Draw("APL");


  //
  // Plot of PE
  //

  TCanvas * c2 = new TCanvas();
	h_pe_bkg->Draw();
	h_pe_sig->SetLineColor(kRed+2);
	h_pe_sig->Draw("same");
	TLegend* l = new TLegend(0.1,0.7,0.48,0.9);
  l->AddEntry(h_pe_sig,"#nu_{#mu} CC in FV","l");
  l->AddEntry(h_pe_bkg,"Background","l");
  l->Draw();
}

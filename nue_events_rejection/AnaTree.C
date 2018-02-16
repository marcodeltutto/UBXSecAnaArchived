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

	Long64_t nentries = fChain->GetEntriesFast();

	TH1D * h_residuals_std_mu = new TH1D("h_residuals_std_mu", "Muons;#sigma_{r_{i}};", 100, 0, 10);
	TH1D * h_residuals_mean_mu = new TH1D("h_residuals_mean_mu", "Muons;<r_{i}>;", 100, -10, 10);
	TH2D * h_residuals_mean_std_mu = new TH2D("h_residuals_mean_std_mu", "Muons;<r_{i}>;#sigma_{r_{i}}", 100, -5, 5, 100, 0, 10);
	TH1D * h_residuals_n_used_hits_mu = new TH1D("h_residuals_n_used_hits_mu", "Muons;Fraction of used hits in cluster;", 100, 0, 1);
	TH1D * h_residuals_r_mu = new TH1D("h_residuals_r_mu", "Muons;Linear correlation coefficient;", 100, -1, 1);
	//TH2D * h_gap_mu = new TH2D("h_gap_mu", "Muons;Number of wires in biggest gap;Number of dead wires in biggest gap", 30, 0, 100, 30, 0, 50);
	TH1D * h_gap_mu = new TH1D("h_gap_mu", "Muons;Number of wires in biggest gap;", 30, 0, 100);
	TH1D * h_maxscatteringangle_mu = new TH1D("h_maxscatteringangle_mu", "Muons;Max Scattering Angle;", 200, 0, 200);

	TH1D * h_residuals_std_el = new TH1D("h_residuals_std_el", "Electrons;#sigma_{r_{i}};", 100, 0, 10);
	TH1D * h_residuals_mean_el = new TH1D("h_residuals_mean_el", "Electrons;<r_{i}>;", 100, -10, 10);
	TH2D * h_residuals_mean_std_el = new TH2D("h_residuals_mean_std_el", "Electrons;<r_{i}>;#sigma_{r_{i}}", 100, -5, 5, 100, 0, 10);
	TH1D * h_residuals_n_used_hits_el = new TH1D("h_residuals_n_used_hits_el", "Electrons;Fraction of used hits in cluster;", 100, 0, 1);
	TH1D * h_residuals_r_el = new TH1D("h_residuals_r_el", "Electrons;Linear correlation coefficient;", 100, -1, 1);
	//TH2D * h_gap_el = new TH2D("h_gap_el", "Electrons;Number of wires in biggest gap;Number of dead wires in biggest gap", 30, 0, 100, 30, 0, 50);
  TH1D * h_gap_el = new TH1D("h_gap_el", "Electrons;Number of wires in biggest gap;", 30, 0, 100);
  TH1D * h_maxscatteringangle_el = new TH1D("h_maxscatteringangle_el", "Electrons;Max Scattering Angle;", 200, 0, 200);

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		for (size_t slice = 0; slice < nslices; slice++) {
			if (std::abs(slc_muoncandidate_truepdg.at(slice)) == 13 && slc_muoncandidate_trueorigin.at(slice) == 0) {
				h_residuals_std_mu->Fill(slc_muoncandidate_residuals_std.at(slice));
				h_residuals_mean_mu->Fill(slc_muoncandidate_residuals_mean.at(slice));
				h_residuals_mean_std_mu->Fill(slc_muoncandidate_residuals_mean.at(slice), slc_muoncandidate_residuals_std.at(slice));
				h_residuals_n_used_hits_mu->Fill(slc_muoncandidate_perc_used_hits_in_cluster.at(slice));
				h_residuals_r_mu->Fill(slc_muoncandidate_linearity.at(slice));
				//h_gap_mu->Fill(slc_muoncandidate_wiregap.at(slice), slc_muoncandidate_wiregap_dead.at(slice));
				h_gap_mu->Fill(slc_muoncandidate_wiregap.at(slice));
				h_maxscatteringangle_mu->Fill(slc_muoncandidate_maxscatteringangle.at(slice));
			}
			if (std::abs(slc_muoncandidate_truepdg.at(slice)) == 11) {
				h_residuals_std_el->Fill(slc_muoncandidate_residuals_std.at(slice));
				h_residuals_mean_el->Fill(slc_muoncandidate_residuals_mean.at(slice));
				h_residuals_mean_std_el->Fill(slc_muoncandidate_residuals_mean.at(slice), slc_muoncandidate_residuals_std.at(slice));
				h_residuals_n_used_hits_el->Fill(slc_muoncandidate_perc_used_hits_in_cluster.at(slice));
				h_residuals_r_el->Fill(slc_muoncandidate_linearity.at(slice));
				//h_gap_el->Fill(slc_muoncandidate_wiregap.at(slice), slc_muoncandidate_wiregap_dead.at(slice));
				h_gap_el->Fill(slc_muoncandidate_wiregap.at(slice));
				h_maxscatteringangle_el->Fill(slc_muoncandidate_maxscatteringangle.at(slice));
			}
		}
	}

	new TCanvas;
	h_residuals_std_mu->DrawNormalized("histo");
	h_residuals_std_el->SetLineColor(kRed+1);
	h_residuals_std_el->DrawNormalized("histo same");

	new TCanvas;
	h_residuals_mean_mu->DrawNormalized("histo");
	h_residuals_mean_el->SetLineColor(kRed+1);
	h_residuals_mean_el->DrawNormalized("histo same");

	new TCanvas;
	h_residuals_mean_std_mu->Draw("colz");
	new TCanvas;
	//h_residuals_mean_std_el->SetLineColor(kRed+1);
	h_residuals_mean_std_el->Draw("colz same");

	new TCanvas;
	h_residuals_n_used_hits_mu->DrawNormalized("histo");
	h_residuals_n_used_hits_el->SetLineColor(kRed+1);
	h_residuals_n_used_hits_el->DrawNormalized("histo same");

	new TCanvas;
	h_residuals_r_mu->DrawNormalized("histo");
	h_residuals_r_el->SetLineColor(kRed+1);
	h_residuals_r_el->DrawNormalized("histo same");

	new TCanvas;
	h_gap_mu->Draw("histo");
	h_gap_el->SetLineColor(kRed+1);
	h_gap_el->Draw("histo same");

	new TCanvas;
	h_maxscatteringangle_mu->Draw("histo");
	h_maxscatteringangle_el->SetLineColor(kRed+1);
	h_maxscatteringangle_el->Draw("histo same");

}

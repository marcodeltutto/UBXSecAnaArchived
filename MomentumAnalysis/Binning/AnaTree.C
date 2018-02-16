#define AnaTree_cxx
#include "AnaTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnaTree::Loop(double down_limit, std::vector<double> &cut_value_v_out, std::vector<double> &truth_fraction_v_out, double & best_up_limit)
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

	TH2D * h_true_reco_mom = new TH2D("h_true_reco_mom", "h_true_reco_mom", 80, 0, 2, 80, 0, 2);
	TH1D * h_mom_resol = new TH1D("h_mom_resol", "h_mom_resol", 80, -2, 2);


	std::vector<double> cut_value_v;
	std::vector<double> truth_fraction_v;


	for (double mom_reco_cut = 5; mom_reco_cut>=down_limit; mom_reco_cut-=0.01){

		int total_truth = 0;
		int selected_truth = 0;

		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

			double true_mom = mom_tree_true;
			double mcs_mom = mom_tree_mcs;

			if (mom_reco_cut == 5) {
				h_true_reco_mom->Fill(true_mom, mcs_mom);
				h_mom_resol->Fill(true_mom - mcs_mom);


			}


			if (mcs_mom < down_limit || mcs_mom > mom_reco_cut) continue;

			total_truth++;

			if (true_mom < down_limit || true_mom > mom_reco_cut) continue;

			selected_truth++;

		}
		cut_value_v.push_back(mom_reco_cut);
		truth_fraction_v.push_back((double)selected_truth/(double)total_truth);

		if ((double)selected_truth/(double)total_truth < 0.6827) {
			best_up_limit = mom_reco_cut;
			break;
		}
	}



	//
  // TGraph
  //

	cut_value_v_out = cut_value_v;
	truth_fraction_v_out = truth_fraction_v;



	double* x = &cut_value_v[0];
	double* y = &truth_fraction_v[0];
	int n = (int) cut_value_v.size();

/*
	TGraph *g = new TGraph(n,x,y);
	g->SetLineColor(kRed+2);
	g->SetMarkerColor(kRed+2);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.5);
	g->SetTitle(";Bla");
	TCanvas * c1 = new TCanvas();
	g->Draw("APL");
	*/



	new TCanvas;
	TLine *line = new TLine(0,0,2,2);
	line->SetLineColor(kRed);
	if (down_limit == 0) {
		h_true_reco_mom->Draw("colz");
		line->Draw();
	}

	new TCanvas;
	if (down_limit == 0) {
		h_mom_resol->Draw("histo");
	}


	std::cout << "Upper limit found: " << best_up_limit << std::endl;

}





void AnaTree::Run()
{

	std::vector<TGraph*> tgraph_v;
	std::vector<int> color_v = {kRed+2, kBlue+2, kGreen+2, kYellow+2, kRed+2, kBlue+2, kGreen+2, kYellow+2};

	double down_limit = 0.;
	double best_up_limit;
	int counter = 0;

	bool exit_flag = false;
	while(!exit_flag) {

		std::vector<double> cut_value_v;
		std::vector<double> truth_fraction_v;  
		this->Loop(down_limit, cut_value_v, truth_fraction_v, best_up_limit);

		std::cout << "after Upper limit found: " << best_up_limit << std::endl;
		down_limit = best_up_limit;


		double* x = &cut_value_v[0];
		double* y = &truth_fraction_v[0];
		int n = (int) cut_value_v.size();


		TGraph *g = new TGraph(n,x,y);
		g->SetLineColor(color_v.at(counter));
		g->SetMarkerColor(color_v.at(counter));
		g->SetMarkerStyle(20);
		g->SetMarkerSize(0.5);
		g->SetTitle(";Muon Momentum [GeV];Percentage of True Values in Reco Bin");
		//TCanvas * c1 = new TCanvas();
		//g->Draw("APL");
		tgraph_v.push_back(g);
		counter++;

		if (best_up_limit == 5) exit_flag = true;
	}

  TCanvas * c1 = new TCanvas();
  tgraph_v.at(0)->Draw("APL");
  for (size_t i = 0; i < tgraph_v.size(); i++) {
  	tgraph_v.at(i)->Draw("PL");
  }
  TLine *line = new TLine(0,0.6827,5,0.6827);
	line->SetLineColor(kRed+1);
	line->SetLineWidth(2);
	line->Draw();

}

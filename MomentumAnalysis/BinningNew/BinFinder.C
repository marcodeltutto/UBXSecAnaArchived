#define BinFinder_cxx
#include "BinFinder.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void BinFinder::Loop()
{
//   In a ROOT session, you can do:
//      root> .L BinFinder.C
//      root> BinFinder t
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

	gStyle->SetPalette(kBird);
	gROOT->SetBatch(kTRUE);

	TH2D * h_true_reco_mom = new TH2D("h_true_reco_mom", ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", 120, 0, 2, 120, 0, 2);
	double bins_mumom[7] = {0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50};
	TH2D * h_true_reco_mom_rightbin = new TH2D("h_true_reco_mom_rightbin", ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", 6, bins_mumom, 6, bins_mumom);

	TH1D * h_reco_mom_pre_truth = new TH1D("h_reco_mom_pre_truth", ";MCS Muon Momentum [GeV];;", 30, 0, 2);

	Long64_t nentries = fChain->GetEntriesFast();

	double lower_bin = 0.0;
	double upper_bin = 2;

	std::vector<double> bin_v;
	bin_v.push_back(lower_bin);

	std::vector<double> points_v, error_x_v, mean_v, std_v;

	bool exit_flag = false;

	while (!exit_flag) {

		std::cout << ">>> Using lower_bin = " << lower_bin << ", and upper_bin = " << upper_bin << std::endl;

		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			double true_mom = mom_tree_true;
			double mcs_mom = mom_tree_mcs;

			h_true_reco_mom->Fill(true_mom, mcs_mom);
			h_true_reco_mom_rightbin->Fill(true_mom, mcs_mom);

			if (true_mom > lower_bin && true_mom < upper_bin) {
				h_reco_mom_pre_truth->Fill(mcs_mom);
			}

		}

		double array[2];
		DoFit(h_reco_mom_pre_truth, "test", array, false);
		std::cout << "mean " << array[0] << ", std " << array[1] << std::endl;
		std::cout << "truth bin: " << upper_bin - lower_bin << std::endl;
		std::cout << "reco bin: " << 2*array[1] << std::endl;
		double truth_bin = upper_bin - lower_bin;
		double reco_bin = 2*array[1];

		if (std::abs(truth_bin - reco_bin) < 0.005) {
			bin_v.push_back(upper_bin);
			points_v.push_back((upper_bin + lower_bin)/2);
			error_x_v.push_back((upper_bin - lower_bin)/2);
			mean_v.push_back(array[0]);
			std_v.push_back(array[1]);
			// Rerun fitting just to save the plot
			ostringstream oss;
			oss <<  lower_bin << "-" << upper_bin;
			std::string label = oss.str();
			DoFit(h_reco_mom_pre_truth, label, array, true);

			lower_bin = upper_bin;
			upper_bin = lower_bin + 2;
		} else {
			upper_bin-=0.005;
		}
		if (lower_bin > 1.7) exit_flag = true;
		h_reco_mom_pre_truth->Reset();
		//exit_flag = true;
	}

	std::cout << "Bins: " << std::endl;
	for (auto b : bin_v)
		std::cout << "\t" << b << std::endl;



	TCanvas * c0 = new TCanvas;
	TLine *line = new TLine(0,0,2,2);
	line->SetLineColor(kRed);
	h_true_reco_mom->Draw("colz");
	line->Draw();

	TString n = "final_nopoints";

	c0->SaveAs("output/" + n + ".pdf");
	c0->SaveAs("output/" + n + ".C");


	TCanvas * c00 = new TCanvas;
	TLine *line2 = new TLine(0,0,2.5,2.5);
	line2->SetLineColor(kRed);
	h_true_reco_mom_rightbin->Draw("colz");
	line2->Draw();

	n = "final_nopoints_rightbin";

	c00->SaveAs("output/" + n + ".pdf");
	c00->SaveAs("output/" + n + ".C");



	double* points = &points_v[0];
	double* mean_array = &mean_v[0];
	double* error_x = &error_x_v[0];
	double* sigma_array = &std_v[0];
	int pts = points_v.size();

	TCanvas * c = new TCanvas;
	TGraphErrors* gr = new TGraphErrors(pts,points,mean_array,error_x,sigma_array);
	gr->SetMarkerStyle(kFullCircle);
	gr->SetMarkerSize(1.0);
	gr->SetLineWidth(2.0);

	gStyle->SetPalette(kBird);
	h_true_reco_mom->Draw("colz");
	gr->Draw("P");

	line->Draw();

  //gr->Fit("pol1", "", "", 0.1, 0.9);
  //TF1 *myfunc = gr->GetFunction("pol1");

	n = "final";

	c->SaveAs("output/" + n + ".pdf");
	c->SaveAs("output/" + n + ".C");
}


void BinFinder::DoFit(TH1D * the_histo, std::string name, double *array, bool save_plot) {

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

	TLatex* range = new TLatex(.9, .95, j.c_str());
	range->SetTextColor(kGray+2);
	range->SetNDC();
	range->SetNDC();
  //range->SetTextSize(2/30.);
	range->SetTextAlign(32);
	range->Draw();

	TString n = name;
	n = n + "_gaus";

	if (save_plot)
		c->SaveAs("output/" + n + ".pdf");

  //_txtfile << name << " GeV |" << array[0] << "|" << array[1] << std::endl;


}


void BinFinder::SmearingMatrix() {


	//int matrix_size = 6 * 7;



	Long64_t nentries = fChain->GetEntriesFast();

  //                           1     2     3     4     5     6       7(overflow)
	double bins_mumom[7] = {0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.5};
	TH2D * h_true_reco_mom_rightbin = new TH2D("h_true_reco_mom_rightbin", ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", 6, bins_mumom, 6, bins_mumom);


	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;


		double true_mom = mom_tree_true;
		double mcs_mom = mom_tree_mcs;

		h_true_reco_mom_rightbin->Fill(true_mom, mcs_mom);


	}

	std::vector<std::vector<double>> S_temp;

	double data[49] = {0.};


	for (int i = 1; i < 8; i++) {      // Reco bin


		std::cout << "This is reco bin " << i << std::endl;

    std::vector<double> p_v;
    p_v.resize(8);

    double sum = 0;

		for (int j = 1; j < 8; j++) {    // True bin

			std::cout << "\tThis is true bin " << j << std::endl;

      p_v.at(j) = h_true_reco_mom_rightbin->GetBinContent(j, i);
      sum += p_v.at(j);

      std::cout << "\tValue is " << p_v.at(j) << std::endl;

		}

		std::cout << "\t>>> Sum is " << sum << std::endl;

    double tot_prob = 0;

		for (int j = 1; j < 8; j++) {
			p_v.at(j) /= sum;

			std::cout << "\t\tProbability at " << j << " is " << p_v.at(j) << std::endl;
			tot_prob += p_v.at(j);

      int row_offset = (i-1)*7;
			data[row_offset + j-1] = p_v.at(j);
		}
		std::cout << "\t\t> Total Probability is " << tot_prob << std::endl;

    S_temp.emplace_back(p_v);

	} // reco bin

	//TMatrixD a(6,6);
	//a.SetMatrixArray(data);
	//a.Print();

  ROOT::Math::SMatrix<double,7> S(data, 49);

  std::cout << S << std::endl;

  ROOT::Math::SVector<double,7> eff_true (0.263441, 0.497797, 0.587106, 0.589119, 0.59103, 0.569231, 1);

  std::cout << S * eff_true << std::endl;

  //double det;
  //S.Det(det);
  //std::cout << "determinant is " << det << std::endl;


  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../ubxsecana_output_bnbcosmic_mcc8.6_stopmu_test15.root");
  if (!f || !f->IsOpen()) {
    f = new TFile("../../ubxsecana_output_bnbcosmic_mcc8.6_stopmu_test15.root");
  }

  TH1D* h_eff_mumom_num = (TH1D*)f->Get("h_eff_mumom_num");
  TH1D* h_eff_mumom_den = (TH1D*)f->Get("h_eff_mumom_den");

  // 
  // Plot without smearing
  //

  TCanvas * c = new TCanvas;
  h_eff_mumom_den->SetTitle("");
  h_eff_mumom_den->GetXaxis()->SetTitle("p_{#mu}^{truth} [GeV]");
  h_eff_mumom_den->GetYaxis()->SetTitle("Events");
  h_eff_mumom_den->SetFillColorAlpha(30, 0.35);
  h_eff_mumom_den->SetLineColor(30);
  h_eff_mumom_den->SetLineWidth(3);

	h_eff_mumom_den->Draw("histo");

  h_eff_mumom_num->SetFillColorAlpha(9, 0.35);
  h_eff_mumom_num->SetLineColor(9);
  h_eff_mumom_num->SetLineWidth(3);

	h_eff_mumom_num->Draw("histo same");

	TLegend * ll = new TLegend(0.5315186,0.7515789,0.8696275,0.8821053,NULL,"brNDC");
  ll->AddEntry(h_eff_mumom_den,"Generated #nu_{#mu} CC in FV","f");
  ll->AddEntry(h_eff_mumom_num,"Selected #nu_{#mu} CC in FV","f");
  ll->Draw();

	TString n = "all_selected";
	n = n + "_true";
	c->SaveAs("output/" + n + ".pdf");

	std::cout << "h_eff_mumom_num->Integral(): " << h_eff_mumom_num->Integral() << std::endl;
	std::cout << "h_eff_mumom_den->Integral(): " << h_eff_mumom_den->Integral() << std::endl;

  //
  // Efficiency (true)
  //

	TEfficiency* teff_true = new TEfficiency(*h_eff_mumom_num,*h_eff_mumom_den);

  TCanvas * c_eff_true = new TCanvas;
  teff_true->SetTitle(";True Muon Momentum [GeV];Efficiency");
  teff_true->SetLineColor(kGreen+3);
  teff_true->SetMarkerColor(kGreen+3);
  teff_true->SetMarkerStyle(20);
  teff_true->SetMarkerSize(0.5);
  teff_true->Draw("AP");

  n = "efficiecy_mumon_true";
	c->SaveAs("output/" + n + ".pdf");

	// 
  // Do the smearing
  //

  ROOT::Math::SVector<double,7> eff_num_true;
  ROOT::Math::SVector<double,7> eff_den_true;

  for (int bin = 1; bin < 8; bin++) {
  	eff_num_true[bin-1] = h_eff_mumom_num->GetBinContent(bin);
  	eff_den_true[bin-1] = h_eff_mumom_den->GetBinContent(bin);
  }

  ROOT::Math::SVector<double,7> eff_num_smear = S * eff_num_true;
  ROOT::Math::SVector<double,7> eff_den_smear = S * eff_den_true;

  TH1D* h_eff_mumom_num_smear = (TH1D*) h_eff_mumom_num->Clone("h_eff_mumom_num_smear");
  TH1D* h_eff_mumom_den_smear = (TH1D*) h_eff_mumom_den->Clone("h_eff_mumom_den_smear");
   
  for (int bin = 1; bin < 8; bin++) {
  	h_eff_mumom_num_smear->SetBinContent(bin, eff_num_smear(bin-1));
  	h_eff_mumom_den_smear->SetBinContent(bin, eff_den_smear(bin-1));
  }

  TCanvas * c_smear = new TCanvas;
  h_eff_mumom_den_smear->SetTitle("");
  h_eff_mumom_den_smear->GetXaxis()->SetTitle("p_{#mu}^{reco} [GeV]");
  h_eff_mumom_den_smear->GetYaxis()->SetTitle("Events");
	h_eff_mumom_den_smear->Draw("histo");
	h_eff_mumom_num_smear->Draw("histo same");

	ll->Draw();

	n = "all_selected_smear";
	c_smear->SaveAs("output/" + n + ".pdf");

  std::cout << "h_eff_mumom_num_smear->Integral(): " << h_eff_mumom_num_smear->Integral() << std::endl;
	std::cout << "h_eff_mumom_den_smear->Integral(): " << h_eff_mumom_den_smear->Integral() << std::endl;


	//
  // Efficiency (reco)
  //

	TEfficiency* teff_reco = new TEfficiency(*h_eff_mumom_num_smear,*h_eff_mumom_den_smear);

  TCanvas * c_eff_reco = new TCanvas;
  teff_reco->SetTitle(";Reco Muon Momentum [GeV];Efficiency");
  teff_reco->SetLineColor(kGreen+3);
  teff_reco->SetMarkerColor(kGreen+3);
  teff_reco->SetMarkerStyle(20);
  teff_reco->SetMarkerSize(0.5);
  teff_reco->Draw("AP");

  n = "efficiecy_mumon_reco";
	c->SaveAs("output/" + n + ".pdf");

}








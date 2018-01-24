void DrawPOT(double pot) {
  
  std::stringstream sstm2;
  sstm2 << "Accumulated POT: " << pot;
  std::string str = sstm2.str();
  
  TLatex* pot_latex_2 = new TLatex(.10, .92, str.c_str());
  pot_latex_2->SetTextFont(62);
  pot_latex_2->SetTextColor(kGray+2);
  pot_latex_2->SetNDC();
  pot_latex_2->SetTextSize(1/30.);
  pot_latex_2->SetTextAlign(10);//left adjusted
  pot_latex_2->Draw();
  
}

void data_mc() {
  
  TFile *f_mc = new TFile("mcs_ana_output_mc_bnbcosmic_finer2.root");
  TFile *f_data_bnbon = new TFile("mcs_ana_output_data_bnbon_finer2.root");
  TFile *f_data_extbnb = new TFile("mcs_ana_output_data_extbnb_finer2.root");
  
  TH1D* h_mc = (TH1D*)f_mc->Get("delta_ll_all");
  TH1D* h_mc_nu = (TH1D*)f_mc->Get("delta_ll_neutrino");
  TH1D* h_mc_cosmic_stopmu = (TH1D*)f_mc->Get("delta_ll_cosmic_stopmu");
  TH1D* h_mc_cosmic_nostopmu = (TH1D*)f_mc->Get("delta_ll_cosmic_nostopmu");

  TH1D* h_data_bnbon = (TH1D*)f_data_bnbon->Get("delta_ll_all");
  TH1D* h_data_extbnb = (TH1D*)f_data_extbnb->Get("delta_ll_all");
  
  double scale_mc = 0.237946;
  double scale_ext = 0.676202;
  
  h_mc->Scale(scale_mc);
  h_data_extbnb->Scale(scale_ext);
  
  TH1D* h_data = (TH1D*) h_data_bnbon->Clone("h_data");
  h_data->Add(h_data_extbnb, -1);
  
  h_mc->SetLineColor(kGreen+2);
  h_mc->SetFillColor(kGreen+2);
  
  
  h_mc_cosmic_stopmu->SetLineColor(kBlue+2);
  h_mc_cosmic_stopmu->SetFillColor(kBlue+2);
  h_mc_cosmic_stopmu->Scale(scale_mc);
  h_mc_cosmic_nostopmu->SetLineColor(kGreen+2);
  h_mc_cosmic_nostopmu->SetFillColor(kGreen+2);
  h_mc_cosmic_nostopmu->Scale(scale_mc);
  h_mc_nu->SetLineColor(kRed+2);
  h_mc_nu->SetFillColor(kRed+2);
  h_mc_nu->Scale(scale_mc);
  
  
  
  /*
  THStack *hs = new THStack("hs",";MCS #DeltaLL;Events");
  h_data_extbnb->SetLineColor(kBlue+2);
  h_data_extbnb->SetFillColor(kBlue+2);
  hs->Add(h_data_extbnb);
  hs->Add(h_mc);
  
  hs->Draw();

  
  h_data_bnbon->SetMarkerStyle(kFullCircle);
  h_data_bnbon->SetMarkerSize(0.6);
  
  h_data_bnbon->Draw("E1 same");
  */
  
  
  THStack *hs2 = new THStack("hs2",";MCS #DeltaLL;Events");
  hs2->Add(h_mc_nu);
  hs2->Add(h_mc_cosmic_nostopmu);
  hs2->Add(h_mc_cosmic_stopmu);
  
  hs2->Draw();
  
  h_data->SetMarkerStyle(kFullCircle);
  h_data->SetMarkerSize(0.6);
  
  h_data->Draw("E1 same");
  
  ll = new TLegend(0.1,0.7,0.48,0.9);
  ll->AddEntry(h_mc_nu,"MC (Neutrino Origin)","f");
  ll->AddEntry(h_mc_cosmic_nostopmu,"MC (Cosmic Origin Through-Going)","f");
  ll->AddEntry(h_mc_cosmic_stopmu,"MC (Cosmic Origin Stopping in TPC)","f");
  ll->AddEntry(h_data,"Data (Beam-on - Beam-off)","lep");
  ll->Draw();
  
  DrawPOT(1.048e+19);

  
  std::cout << "mc integral till -0.5: " << h_mc->Integral(1, h_mc->GetXaxis()->FindBin(-0.5)) << std::endl;
  std::cout << "data integral till -0.5: " << h_data->Integral(1, h_data->GetXaxis()->FindBin(-0.5)) << std::endl;


}




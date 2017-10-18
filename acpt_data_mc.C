void acpt_data_mc() {

  TFile * bnbon_f  = TFile::Open("files/mcc8.3_v1.0.2/ubxsec_output_data_bnbon_mcc8.3_v1.0.2.root", "READ");
  TFile * extbnb_f = TFile::Open("files/mcc8.3_v1.0.2/ubxsec_output_data_extbnb_mcc8.3_v1.0.2.root", "READ");

  bnbon_f->cd("pandoraCosmicACPTTagger");
  extbnb_f->cd("pandoraCosmicACPTTagger");
  
  TH1D * h_bnbon = (TH1D*) bnbon_f->Get("pandoraCosmicACPTTagger/h_diff");
  TH1D * h_extbnb = (TH1D*) extbnb_f->Get("pandoraCosmicACPTTagger/h_diff");

  std::cout << "here1" << std::endl;

  //h_bnbon->Draw();
  h_bnbon->Scale(1);
  h_extbnb->Scale(1.23377);
  std::cout << "here2" << std::endl;

  TH1D* h_data = (TH1D*)h_bnbon->Clone("h_data");
  std::cout << "here3" << std::endl;

  h_data->Sumw2();
  h_data->Add(h_extbnb, -1.);
  std::cout << "here4" << std::endl;


  h_data->Draw();
}

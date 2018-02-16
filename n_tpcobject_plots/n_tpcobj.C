void n_tpcobj() {

  TH1D * h_mcc7 = new TH1D("h_mcc7", "", 15, 0, 15);
  
  h_mcc7->SetBinContent(1, 17); 
  h_mcc7->SetBinContent(2, 26);
  h_mcc7->SetBinContent(3, 78);
  h_mcc7->SetBinContent(4, 139);
  h_mcc7->SetBinContent(5, 212);
  h_mcc7->SetBinContent(6, 299);
  h_mcc7->SetBinContent(7, 455);
  h_mcc7->SetBinContent(8, 376);
  h_mcc7->SetBinContent(9, 271);
  h_mcc7->SetBinContent(10, 189);
  h_mcc7->SetBinContent(11, 102);
  h_mcc7->SetBinContent(12, 77);
  h_mcc7->SetBinContent(13, 59);
  h_mcc7->SetBinContent(14, 25);
  h_mcc7->SetBinContent(15, 20);
  h_mcc7->SetBinContent(16, 5);


  TH1D * h_mcc8 = new TH1D("h_mcc8", "", 15, 0, 15);

  h_mcc8->SetBinContent(1, 1042);
  h_mcc8->SetBinContent(2, 9057);
  h_mcc8->SetBinContent(3, 11498);
  h_mcc8->SetBinContent(4, 9467);
  h_mcc8->SetBinContent(5, 6012);
  h_mcc8->SetBinContent(6, 3093);
  h_mcc8->SetBinContent(7, 1067);
  h_mcc8->SetBinContent(8, 546);
  h_mcc8->SetBinContent(9, 112);
  h_mcc8->SetBinContent(10, 27);
  h_mcc8->SetBinContent(11, 2);
  h_mcc8->SetBinContent(12, 0);
  h_mcc8->SetBinContent(13, 0);
  h_mcc8->SetBinContent(14, 0);
  h_mcc8->SetBinContent(15, 0);
  h_mcc8->SetBinContent(16, 0);


  TH1D * h_mcc8_old = new TH1D("h_mcc8", "", 15, 0, 15);

  h_mcc8_old->SetBinContent(1, 27);
  h_mcc8_old->SetBinContent(2, 1109);
  h_mcc8_old->SetBinContent(3, 3310);
  h_mcc8_old->SetBinContent(4, 5001);
  h_mcc8_old->SetBinContent(5, 4567);
  h_mcc8_old->SetBinContent(6, 3278);
  h_mcc8_old->SetBinContent(7, 1956);
  h_mcc8_old->SetBinContent(8, 992);
  h_mcc8_old->SetBinContent(9, 476);
  h_mcc8_old->SetBinContent(10, 325);
  h_mcc8_old->SetBinContent(11, 17);
  h_mcc8_old->SetBinContent(12, 3);
  h_mcc8_old->SetBinContent(13, 0);
  h_mcc8_old->SetBinContent(14, 0);
  h_mcc8_old->SetBinContent(15, 0);
  h_mcc8_old->SetBinContent(16, 0);


  h_mcc8->GetXaxis()->SetTitle("Number of TPCObjects per Event");
  h_mcc8->SetFillColorAlpha(30, 0.35);
  h_mcc8->SetLineColor(30);
  h_mcc8->SetLineWidth(2);
  h_mcc8->DrawNormalized();

  h_mcc8_old->GetXaxis()->SetTitle("Number of TPCObjects per Event");
  h_mcc8_old->SetFillColorAlpha(50, 0.35);
  h_mcc8_old->SetLineColor(50);
  h_mcc8_old->SetLineWidth(2);
  h_mcc8_old->DrawNormalized();

  h_mcc7->GetXaxis()->SetTitle("Number of TPCObjects per Event");
  h_mcc7->SetFillColorAlpha(9, 0.35);
  h_mcc7->SetLineColor(9);
  h_mcc7->SetLineWidth(2);
   
  h_mcc7->DrawNormalized();
  //h_mcc8_old->DrawNormalized("same");
  h_mcc8->DrawNormalized("same");

  TLegend * ll = new TLegend(0.5315186,0.7515789,0.8696275,0.8821053,NULL,"brNDC");
  ll->AddEntry(h_mcc7,"MCC7","f");
  ll->AddEntry(h_mcc8,"This Selection","f");
  ll->Draw();
}


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


//**********************************************************

TLegend* DrawTHStack(THStack *hs_trklen,
                   double pot_scaling,
                   bool _breakdownPlots,
                   std::map<std::string,TH1D*> themap){
  
  
  for (auto iter : themap) {
    iter.second->Scale(pot_scaling);
  }
  
  
  if (_breakdownPlots) {
    themap["cosmic_nostopmu"]->SetLineColor(kBlue+2);
    themap["cosmic_nostopmu"]->SetFillColor(kBlue+2);
    hs_trklen->Add(themap["cosmic_nostopmu"]);
    themap["cosmic_stopmu"]->SetLineColor(kBlue);
    themap["cosmic_stopmu"]->SetFillColor(kBlue);
    hs_trklen->Add(themap["cosmic_stopmu"]);
    themap["outfv_nostopmu"]->SetLineColor(kGreen+3);
    themap["outfv_nostopmu"]->SetFillColor(kGreen+3);
    hs_trklen->Add(themap["outfv_nostopmu"]);
    themap["outfv_stopmu"]->SetLineColor(kGreen+2);
    themap["outfv_stopmu"]->SetFillColor(kGreen+2);
    hs_trklen->Add(themap["outfv_stopmu"]);
    themap["nc_proton"]->SetLineColor(kGray+2);
    themap["nc_proton"]->SetFillColor(kGray+2);
    hs_trklen->Add(themap["nc_proton"]);
    themap["nc_pion"]->SetLineColor(kGray+1);
    themap["nc_pion"]->SetFillColor(kGray+1);
    hs_trklen->Add(themap["nc_pion"]);
    themap["nc_other"]->SetLineColor(kGray);
    themap["nc_other"]->SetFillColor(kGray);
    hs_trklen->Add(themap["nc_other"]);
  }
  else {
    themap["cosmic"]->SetLineColor(kBlue+2);
    themap["cosmic"]->SetFillColor(kBlue+2);
    hs_trklen->Add(themap["cosmic"]);
    themap["outfv"]->SetLineColor(kGreen+2);
    themap["outfv"]->SetFillColor(kGreen+2);
    hs_trklen->Add(themap["outfv"]);
    themap["nc"]->SetLineColor(kGray);
    themap["nc"]->SetFillColor(kGray);
    hs_trklen->Add(themap["nc"]);
  }
  
  themap["anumu"]->SetLineColor(kOrange-3);
  themap["anumu"]->SetFillColor(kOrange-3);
  hs_trklen->Add(themap["anumu"]);
  themap["nue"]->SetLineColor(kMagenta+1);
  themap["nue"]->SetFillColor(kMagenta+1);
  hs_trklen->Add(themap["nue"]);
  
  if (_breakdownPlots) {
    themap["signal_nostopmu"]->SetLineColor(kRed+2);
    themap["signal_nostopmu"]->SetFillColor(kRed+2);
    hs_trklen->Add(themap["signal_nostopmu"]);
    themap["signal_stopmu"]->SetLineColor(kRed);
    themap["signal_stopmu"]->SetFillColor(kRed);
    hs_trklen->Add(themap["signal_stopmu"]);
  }
  else {
    themap["signal"]->SetLineColor(kRed);
    themap["signal"]->SetFillColor(kRed);
    hs_trklen->Add(themap["signal"]);
  }
  hs_trklen->Draw();
  
  //h_trklen_total->DrawCopy("hist");
  
  //gStyle->SetHatchesLineWidth(1);
  themap["total"]->SetFillColor(kBlack);
  themap["total"]->SetFillStyle(3005);
  themap["total"]->Draw("E2 same");
  
  
  
  
  
  
  
  TLegend* leg2;
  if (_breakdownPlots){
    leg2 = new TLegend(0.56,0.37,0.82,0.82,NULL,"brNDC");
  } else {
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
  }
  std::stringstream sstm;
  // numu
  if (_breakdownPlots) {
    leg2->AddEntry(themap["signal_stopmu"],"#nu_{#mu} CC (stopping #mu)","f");
    leg2->AddEntry(themap["signal_nostopmu"],"#nu_{#mu} CC (other)","f");
  } else {
    sstm << "#nu_{#mu} CC (signal), " << std::setprecision(2)  << themap["signal"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["signal"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  
  // nue
  sstm << "#nu_{e}, #bar{#nu}_{e} CC, " << std::setprecision(2)  << themap["nue"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["nue"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // anumu
  sstm << "#bar{#nu}_{#mu} CC, " << std::setprecision(2)  << themap["anumu"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["anumu"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // nc, outfv, cosmic
  if (_breakdownPlots) {
    leg2->AddEntry(themap["nc_other"],"NC (other)","f");
    leg2->AddEntry(themap["nc_pion"],"NC (pion)","f");
    leg2->AddEntry(themap["nc_proton"],"NC (proton)","f");
    leg2->AddEntry(themap["outfv_stopmu"],"OUTFV (stopping #mu)","f");
    leg2->AddEntry(themap["outfv_nostopmu"],"OUTFV (other)","f");
    leg2->AddEntry(themap["cosmic_stopmu"],"Cosmic (stopping #mu)","f");
    leg2->AddEntry(themap["cosmic_nostopmu"],"Cosmic (other)","f");
  } else {
    sstm << "NC, " << std::setprecision(2)  << themap["nc"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["nc"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "OUTFV, " << std::setprecision(2)  << themap["outfv"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["outfv"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "Cosmic, " << std::setprecision(2)  << themap["cosmic"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["cosmic"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  leg2->AddEntry(themap["total"],"MC Stat Unc.","f");
  leg2->Draw();
  
  return leg2;

}





//**********************************************************
TLegend* DrawTHStack2(THStack *hs_trklen,
                  double pot_scaling,
                  bool _breakdownPlots,
                  std::map<std::string,TH1D*> themap){
  
  
  for (auto iter : themap) {
    iter.second->Scale(pot_scaling);
  }
  
  
  
  themap["background"]->SetLineColor(kBlue+2);
  themap["background"]->SetFillColor(kBlue+2);
  hs_trklen->Add(themap["background"]);
  
  themap["signal"]->SetLineColor(kRed+2);
  themap["signal"]->SetFillColor(kRed+2);
  hs_trklen->Add(themap["signal"]);
  
  hs_trklen->Draw();
  
  themap["total"]->SetFillColor(kBlack);
  themap["total"]->SetFillStyle(3005);
  themap["total"]->Draw("E2 same");
  
  
  
  
  TLegend* leg2;
  leg2 = new TLegend(0.13,0.69,0.45,0.87,NULL,"brNDC");

  std::stringstream sstm;
  
  sstm << "Signal, " << std::setprecision(2)  << themap["signal"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["signal"],sstm.str().c_str(),"f");
  sstm.str("");
  
  sstm << "Background, " << std::setprecision(2)  << themap["background"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["background"],sstm.str().c_str(),"f");
  sstm.str("");
  
  leg2->AddEntry(themap["total"],"MC Stat Unc.","f");
  leg2->Draw();
  
  return leg2;

}

void DrawDataHisto(TH1D* histo) {

  histo->SetMarkerStyle(kFullCircle);
  histo->SetMarkerSize(0.6);

  histo->Draw("E1 same");
  
}




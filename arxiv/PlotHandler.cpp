#define PlotHandler_cxx

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

using namespace std;

#include "PlotHandler.hpp"
#include "AnaTree.h"
#include "Spectrum.hpp"
#include "SelectionTools.hpp"

#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TFile.h"


//____________________________________________________________________
void PlotHandler::InstantiateIntermidiatePlots(double totalPOT) {

  std::vector<std::string> cutname;
  cutname.resize(6);
  cutname.at(0) = "flashtag";
  cutname.at(1) = "vertexcontained";
  cutname.at(2) = "selectbesttrack";
  cutname.at(3) = "flashmatch";
  cutname.at(4) = "trackcontained";
  cutname.at(5) = "longtrack";

  for (unsigned int i = 0; i < cutname.size(); i++) {
    CutPlots cutplots;
    cutplots.Snuenergy_numu    = new Spectrum("nuenergy_numu_"+cutname.at(i),  "#nu_{#mu};Neutrino Energy [GeV];Selected Events up to "+cutname.at(i),80,0,10,totalPOT);
    cutplots.Snuenergy_anumu   = new Spectrum("nuenergy_anumu_"+cutname.at(i), "#bar{#nu}_{#mu};Neutrino Energy [GeV];Selected Events up to "+cutname.at(i),80,0,10,totalPOT);
    cutplots.Snuenergy_nue     = new Spectrum("nuenergy_nue_"+cutname.at(i),   "#nu_{e}/#bar{#nu}_{e};Neutrino Energy [GeV];Selected Events up to "+cutname.at(i),80,0,10, totalPOT);
    cutplots.Snuenergy_nc      = new Spectrum("nuenergy_nc_"+cutname.at(i),    "NC;Neutrino Energy [GeV];Selected Events up to "+cutname.at(i),80,0,10,totalPOT);
    cutplots.Snuenergy_cosmics = new Spectrum("nuenergy_cosmics_"+cutname.at(i),"Cosmics;Neutrino Energy [GeV];Selected Events up to "+cutname.at(i),80,0,10,totalPOT);
    cutplots.Snuenergy_outfv   = new Spectrum("nuenergy_outfv_"+cutname.at(i),"Outside FV;Neutrino Energy [GeV];Selected Events up to "+cutname.at(i),80,0,10,totalPOT);


    cutplots.Strklen_numu      = new Spectrum("trklen_numu_"+cutname.at(i),  "#nu_{#mu};Track Length [cm];Selected Events up to "+cutname.at(i),20,0,1036.8,totalPOT);
    cutplots.Strklen_anumu     = new Spectrum("trklen_anumu_"+cutname.at(i), "#bar{#nu}_{#mu};Track Length [cm];Selected Events up to "+cutname.at(i),20,0,1036.8,totalPOT);
    cutplots.Strklen_nue       = new Spectrum("trklen_nue_"+cutname.at(i),   "#nu_{e};Track Length [cm];Selected Events up to "+cutname.at(i),20,0,1036.8,totalPOT);
    cutplots.Strklen_nc        = new Spectrum("trklen_nc_"+cutname.at(i),    "NC;Track Length [cm];Selected Events up to "+cutname.at(i),20,0,1036.8,totalPOT);
    cutplots.Strklen_cosmics   = new Spectrum("trklen_cosmics_"+cutname.at(i),"Cosmics;Track Length [cm];Selected Events up to "+cutname.at(i),20,0,1036.8,totalPOT);
    cutplots.Strklen_outfv     = new Spectrum("trklen_outfv_"+cutname.at(i),"Outside FV;Track Length [cm];Selected Events up to "+cutname.at(i),20,0,1036.8,totalPOT);

    cutplots.Scostheta_numu      = new Spectrum("costheta_numu_"+cutname.at(i),  "#nu_{#mu};cos#theta;Selected Events up to "+cutname.at(i),20,-1,1,totalPOT);
    cutplots.Scostheta_anumu     = new Spectrum("costheta_anumu_"+cutname.at(i), "#bar{#nu}_{#mu};cos#theta;Selected Events up to "+cutname.at(i),20,-1,1,totalPOT);
    cutplots.Scostheta_nue       = new Spectrum("costheta_nue_"+cutname.at(i),   "#nu_{e};cos#theta;Selected Events up to "+cutname.at(i),20,-1,1,totalPOT);
    cutplots.Scostheta_nc        = new Spectrum("costheta_nc_"+cutname.at(i),    "NC;cos#theta;Selected Events up to "+cutname.at(i),20,-1,1,totalPOT);
    cutplots.Scostheta_cosmics   = new Spectrum("costheta_cosmics_"+cutname.at(i),"Cosmics;cos#theta;Selected Events up to "+cutname.at(i),20,-1,1,totalPOT);
    cutplots.Scostheta_outfv     = new Spectrum("costheta_outfv_"+cutname.at(i),"Outside FV;cos#theta;Selected Events up to "+cutname.at(i),20,-1,1,totalPOT);

    cutToPlotsMap.emplace(cutname[i],cutplots);

  }


}


//____________________________________________________________________
std::string PlotHandler::MakeIntermidiatePlots(std::string cutString, AnaTree * anatree, int bestTrackID){

  CutPlots cutplots = cutToPlotsMap.find(cutString)->second;
  std::string status;

  // *******************************
  // We don't have a candidate track
  // *******************************
  if(bestTrackID == -1) {
    if(anatree->ccnc_truth[0] == 0) {

      if(!SelectionTools::InFV(anatree->nuvtxx_truth[0],
                               anatree->nuvtxy_truth[0],
                               anatree->nuvtxz_truth[0])
         && anatree->nuPDG_truth[0]==14) {

        cutplots.Snuenergy_outfv -> Fill(anatree->enu_truth[0]);
      }
      else if(anatree->nuPDG_truth[0]==14) cutplots.Snuenergy_numu   -> Fill(anatree->enu_truth[0]);
      else if(anatree->nuPDG_truth[0]==-14) cutplots.Snuenergy_anumu -> Fill(anatree->enu_truth[0]);
      else if(anatree->nuPDG_truth[0]==12 || anatree->nuPDG_truth[0]==-12) cutplots.Snuenergy_nue -> Fill(anatree->enu_truth[0]);

    } // end if cc
    else if(anatree->ccnc_truth[0] == 1) cutplots.Snuenergy_nc -> Fill(anatree->enu_truth[0]);
  }
  // *******************************
  // We have a candidate track
  // *******************************
  else if(bestTrackID > -1) {

    double trklen = SelectionTools::GetTrackLength(anatree,bestTrackID);
    double cosTheta = anatree->trktheta_pandoraNu[bestTrackID];

    // Is CC and from neutrino
    if(anatree->ccnc_truth[0] == 0
       && anatree->trkorigin_pandoraNu[bestTrackID][anatree->trkpidbestplane_pandoraNu[bestTrackID]] == 1) {

      if(!SelectionTools::InFV(anatree->nuvtxx_truth[0],
                               anatree->nuvtxy_truth[0],
                               anatree->nuvtxz_truth[0])
         && anatree->nuPDG_truth[0]==14){
        cutplots.Snuenergy_outfv -> Fill(anatree->enu_truth[0]);
        cutplots.Strklen_outfv -> Fill(trklen);
        cutplots.Scostheta_outfv -> Fill(cosTheta);
      }
      else if(anatree->nuPDG_truth[0]==14) {
        cutplots.Snuenergy_numu   -> Fill(anatree->enu_truth[0]);
        cutplots.Strklen_numu -> Fill(trklen);
        cutplots.Scostheta_numu -> Fill(cosTheta);
        //cout << "Is signal." << endl;
        status = "isSignal";
      }
      else if(anatree->nuPDG_truth[0]==-14) {
        cutplots.Snuenergy_anumu -> Fill(anatree->enu_truth[0]);
        cutplots.Strklen_anumu -> Fill(trklen);
        cutplots.Scostheta_anumu -> Fill(cosTheta);
      }
      else if(anatree->nuPDG_truth[0]==12 || anatree->nuPDG_truth[0]==-12) {
        cutplots.Snuenergy_nue -> Fill(anatree->enu_truth[0]);
        cutplots.Strklen_nue -> Fill(trklen);
        cutplots.Scostheta_nue -> Fill(cosTheta);
      }
    }
    // Is NC and from neutrino
    else if(anatree->ccnc_truth[0] == 1
            && anatree->trkorigin_pandoraNu[bestTrackID][anatree->trkpidbestplane_pandoraNu[bestTrackID]] == 1) {
      cutplots.Snuenergy_nc -> Fill(anatree->enu_truth[0]);
      cutplots.Strklen_nc -> Fill(trklen);
      cutplots.Scostheta_nc -> Fill(cosTheta);
    }
    // Is from cosmic
    else if(anatree->trkorigin_pandoraNu[bestTrackID][anatree->trkpidbestplane_pandoraNu[bestTrackID]] != 1){
      cutplots.Snuenergy_cosmics -> Fill(anatree->enu_truth[0]);
      cutplots.Strklen_cosmics -> Fill(trklen);
      cutplots.Scostheta_cosmics -> Fill(cosTheta);
      //cout << "Is cosmic. trkorigin is " << anatree->trkorigin_pandoraNu[bestTrackID][anatree->trkpidbestplane_pandoraNu[bestTrackID]] << endl;
      status = "isCosmic";
    }
  }
  return status;
}

//____________________________________________________________________
THStack* PlotHandler::MakeStackHisto(std::string name, std::string label, TLegend* leg, Spectrum *s1, Spectrum *s2, Spectrum *s3, Spectrum *s4, Spectrum *s5, Spectrum *s6) {

    THStack *hs = new THStack(("hs_"+label).c_str(),(name).c_str());
    TH1D* h1 = s1->ToTH1D();
    h1->SetLineColor(kBlack);
    h1->SetFillColor(kGray);
    hs->Add(h1);
    TH1D* h2 = s2->ToTH1D();
    h2->SetLineColor(kBlack);
    h2->SetFillColor(kOrange+3);
    hs->Add(h2);
    TH1D* h3 = s3 ->ToTH1D();
    h3->SetLineColor(kBlack);
    h3->SetFillColor(kGreen+2);
    hs->Add(h3);
    TH1D* h4 = s4 ->ToTH1D();
    h4->SetLineColor(kBlack);
    h4->SetFillColor(kOrange);
    hs->Add(h4);
    TH1D* h5 = s5 ->ToTH1D();
    h5->SetLineColor(kBlack);
    h5->SetFillColor(kBlue);
    hs->Add(h5);
    TH1D* h6 = s6 ->ToTH1D();
    h6->SetLineColor(kBlack);
    h6->SetFillColor(kRed+1);
    hs->Add(h6);

    leg->AddEntry(h6,"#nu_{#mu}","f");
    leg->AddEntry(h3,"NC","f");
    leg->AddEntry(h2,"#nu_{e}","f");
    leg->AddEntry(h1,"#bar{#nu}_{#mu}","f");
    leg->AddEntry(h4,"#nu_{#mu} CC Out of FV","f");

    leg->SetBorderSize(1);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);

    if(label != "flashtag" || label != "vertexcontained") leg->AddEntry(h5,"Cosmic bgr events","f");

    return hs;

}




//____________________________________________________________________
void PlotHandler::SaveIntermidiatePlots() {

  TFile *f = new TFile("thstacktest.root","RECREATE");
  f->cd();

  for (std::map<std::string,CutPlots>::iterator it=cutToPlotsMap.begin(); it!=cutToPlotsMap.end(); ++it) {
    std::string label = it->first;
    std::cout << "Saving plot " << label << endl;
    CutPlots cutplots = it->second;

    TLegend *leg = new TLegend(0.5716332,0.5242105,0.8567335,0.8484211,NULL,"brNDC");
    THStack *hs_nuenergy = this->MakeStackHisto(";True Neutrino Energy [GeV];Selected events up to " + label,"nuenergy_"+label,leg,
                                                                                                    cutplots.Snuenergy_anumu,
                                                                                                    cutplots.Snuenergy_nue,
                                                                                                    cutplots.Snuenergy_nc,
                                                                                                    cutplots.Snuenergy_outfv,
                                                                                                    cutplots.Snuenergy_cosmics,
                                                                                                    cutplots.Snuenergy_numu);
    TLegend *leg2 = new TLegend(0.5716332,0.5242105,0.8567335,0.8484211,NULL,"brNDC");
    THStack *hs_trklen = this->MakeStackHisto(";Track Length [cm];Selected events up to " + label,"trklen_"+label,leg2,
                                                                                                    cutplots.Strklen_anumu,
                                                                                                    cutplots.Strklen_nue,
                                                                                                    cutplots.Strklen_nc,
                                                                                                    cutplots.Strklen_outfv,
                                                                                                    cutplots.Strklen_cosmics,
                                                                                                    cutplots.Strklen_numu);

    TLegend *leg3 = new TLegend(0.5716332,0.5242105,0.8567335,0.8484211,NULL,"brNDC");
    THStack *hs_costheta = this->MakeStackHisto(";cos#theta;Selected events up to " + label,"costheta_"+label,leg2,
                                                                                                    cutplots.Scostheta_anumu,
                                                                                                    cutplots.Scostheta_nue,
                                                                                                    cutplots.Scostheta_nc,
                                                                                                    cutplots.Scostheta_outfv,
                                                                                                    cutplots.Scostheta_cosmics,
                                                                                                    cutplots.Scostheta_numu);

    hs_nuenergy->Write();
    hs_trklen->Write();
    hs_costheta->Write();
    leg->Write("legend");

    cutplots.Snuenergy_numu      ->Save();
    cutplots.Snuenergy_anumu     ->Save();
    cutplots.Snuenergy_nue       ->Save();
    cutplots.Snuenergy_nc        ->Save();
    cutplots.Snuenergy_outfv     ->Save();
    cutplots.Snuenergy_cosmics   ->Save();


  }
  f->Close();
}




//____________________________________________________________________
void PlotHandler::MakeRatioPlots(std::string stage1, std::string stage2) {


  std::map<std::string,CutPlots>::iterator it = cutToPlotsMap.find(stage1);
  if (it == cutToPlotsMap.end()) {
    std::cout << "Error in PlotHandler::MakeRatioPlots. Cannot find key " << stage1 << endl;
    exit(0);
  }
  std::map<std::string,CutPlots>::iterator it2 = cutToPlotsMap.find(stage2);
  if (it2 == cutToPlotsMap.end()) {
    std::cout << "Error in PlotHandler::MakeRatioPlots. Cannot find key " << stage2 << endl;
    exit(0);
  }

  CutPlots cutplots1 = cutToPlotsMap.find(stage1)->second;
  CutPlots cutplots2 = cutToPlotsMap.find(stage2)->second;

  // NuMu
  Spectrum *Snuenergy_numu_ratio = Spectrum::MakeRatio(cutplots1.Snuenergy_numu,cutplots2.Snuenergy_numu, "nuenergy_numu_ratio", stage1+"_"+stage2+"numu_ratio;True Neutrino Energy [GeV];Ratio");
  Spectrum *Strklen_numu_ratio = Spectrum::MakeRatio(cutplots1.Strklen_numu,cutplots2.Strklen_numu, "trklen_numu_ratio", stage1+"_"+stage2+"numu_ratio;Track Range [cm];Ratio");
  Spectrum *Scostheta_numu_ratio = Spectrum::MakeRatio(cutplots1.Scostheta_numu,cutplots2.Scostheta_numu, "costheta_numu_ratio", stage1+"_"+stage2+"numu_ratio;cos#theta;Ratio");

  // NuE
  Spectrum *Snuenergy_nue_ratio = Spectrum::MakeRatio(cutplots1.Snuenergy_nue,cutplots2.Snuenergy_nue, "nuenergy_nue_ratio", stage1+"_"+stage2+"nue_ratio;True Neutrino Energy [GeV];Ratio");
  Spectrum *Strklen_nue_ratio = Spectrum::MakeRatio(cutplots1.Strklen_nue,cutplots2.Strklen_nue, "trklen_nue_ratio", stage1+"_"+stage2+"nue_ratio;True Neutrino Energy[GeV];Ratio");
  Spectrum *Scostheta_nue_ratio = Spectrum::MakeRatio(cutplots1.Scostheta_nue,cutplots2.Scostheta_nue, "costheta_nue_ratio", stage1+"_"+stage2+"nue_ratio;cos#theta;Ratio");

  // ANuMu
  Spectrum *Snuenergy_anumu_ratio = Spectrum::MakeRatio(cutplots1.Snuenergy_anumu,cutplots2.Snuenergy_anumu, "nuenergy_anumu_ratio", stage1+"_"+stage2+"anumu_ratio;True Neutrino Energy [GeV];Ratio");
  Spectrum *Strklen_anumu_ratio = Spectrum::MakeRatio(cutplots1.Strklen_anumu,cutplots2.Strklen_anumu, "trklen_anumu_ratio", stage1+"_"+stage2+"anumu_ratio;True Neutrino Energy[GeV];Ratio");
  Spectrum *Scostheta_anumu_ratio = Spectrum::MakeRatio(cutplots1.Scostheta_anumu,cutplots2.Scostheta_anumu, "costheta_anumu_ratio", stage1+"_"+stage2+"anumu_ratio;cos#theta;Ratio");

  // NC
  Spectrum *Snuenergy_nc_ratio = Spectrum::MakeRatio(cutplots1.Snuenergy_nc,cutplots2.Snuenergy_nc, "nuenergy_nc_ratio", stage1+"_"+stage2+"nc_ratio;True Neutrino Energy [GeV];Ratio");
  Spectrum *Strklen_nc_ratio = Spectrum::MakeRatio(cutplots1.Strklen_nc,cutplots2.Strklen_nc, "trklen_nc_ratio", stage1+"_"+stage2+"nc_ratio;True Neutrino Energy[GeV];Ratio");
  Spectrum *Scostheta_nc_ratio = Spectrum::MakeRatio(cutplots1.Scostheta_nc,cutplots2.Scostheta_nc, "costheta_nc_ratio", stage1+"_"+stage2+"nc_ratio;cos#theta;Ratio");

  // Cosmic
  Spectrum *Snuenergy_cosmics_ratio = Spectrum::MakeRatio(cutplots1.Snuenergy_cosmics,cutplots2.Snuenergy_cosmics, "nuenergy_cosmics_ratio", stage1+"_"+stage2+"cosmics_ratio;True Neutrino Energy [GeV];Ratio");
  Spectrum *Strklen_cosmics_ratio = Spectrum::MakeRatio(cutplots1.Strklen_cosmics,cutplots2.Strklen_cosmics, "trklen_cosmics_ratio", stage1+"_"+stage2+"cosmics_ratio;True Neutrino Energy[GeV];Ratio");
  Spectrum *Scostheta_cosmics_ratio = Spectrum::MakeRatio(cutplots1.Scostheta_cosmics,cutplots2.Scostheta_cosmics, "costheta_cosmics_ratio", stage1+"_"+stage2+"cosmics_ratio;cos#theta;Ratio");

  // OutFV
  Spectrum *Snuenergy_outfv_ratio = Spectrum::MakeRatio(cutplots1.Snuenergy_outfv,cutplots2.Snuenergy_cosmics, "nuenergy_outfv_ratio", stage1+"_"+stage2+"outfv_ratio;True Neutrino Energy [GeV];Ratio");
  Spectrum *Strklen_outfv_ratio = Spectrum::MakeRatio(cutplots1.Strklen_outfv,cutplots2.Strklen_outfv, "trklen_outfv_ratio", stage1+"_"+stage2+"outfv_ratio;True Neutrino Energy[GeV];Ratio");
  Spectrum *Scostheta_outfv_ratio = Spectrum::MakeRatio(cutplots1.Scostheta_outfv,cutplots2.Scostheta_outfv, "costheta_outfv_ratio", stage1+"_"+stage2+"outfv_ratio;cos#theta;Ratio");


  TFile *f = new TFile((stage1+"_"+stage2+"_ratioplots.root").c_str(),"RECREATE");
  f->cd();

  TH1D* h1 = Strklen_numu_ratio->ToTH1D(); h1->Write();
  h1 = Scostheta_numu_ratio->ToTH1D(); h1->Write();
  h1 = Snuenergy_numu_ratio->ToTH1D(); h1->Write();

  TH1D* h2 = Strklen_anumu_ratio->ToTH1D(); h2->Write();
  h2 = Scostheta_anumu_ratio->ToTH1D(); h2->Write();
  h2 = Scostheta_anumu_ratio->ToTH1D(); h2->Write();

  TH1D* h3 = Strklen_nue_ratio->ToTH1D(); h3->Write();
  h3 = Scostheta_nue_ratio->ToTH1D(); h3->Write();
  h3 = Snuenergy_nue_ratio->ToTH1D(); h3->Write();

  TH1D* h4 = Strklen_nc_ratio->ToTH1D(); h4->Write();
  h4 = Scostheta_nc_ratio->ToTH1D(); h4->Write();
  h4 = Snuenergy_nc_ratio->ToTH1D(); h4->Write();

  TH1D* h5 = Strklen_cosmics_ratio->ToTH1D(); h5->Write();
  h5 = Scostheta_cosmics_ratio->ToTH1D(); h5->Write();
  h5 = Snuenergy_cosmics_ratio->ToTH1D(); h5->Write();

  TH1D* h6 = Strklen_outfv_ratio->ToTH1D(); h6->Write();
  h6 = Scostheta_outfv_ratio->ToTH1D(); h6->Write();
  h6 = Snuenergy_outfv_ratio->ToTH1D(); h6->Write();

  f->Close();

  Strklen_numu_ratio->Save();
  Strklen_nue_ratio->Save();
  Strklen_anumu_ratio->Save();
  Strklen_nc_ratio->Save();
  Strklen_cosmics_ratio->Save();
  Strklen_outfv_ratio->Save();




/*
  TFile *f = new TFile("thstacktest.root","READ");

  std::string stage1 = "longtrack";
  THStack *hs_1 = (THStack*) f->Get(("hs_nuenergy_" + stage1).c_str());
  TList *tl_1 = hs_1->GetHists();
  TIter next(tl_1);
  TObject *obj;
  while ((obj = next())) {
    if (obj->GetTitle() == "#bar{#nu}_{#mu}")
  std::cout << "TITLE is " << obj->GetTitle();
  }
*/
}












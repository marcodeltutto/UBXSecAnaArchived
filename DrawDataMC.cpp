#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include <TRandom1.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TChain.h>
#include "TThread.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>

#include "AnaTree.h"
#include "Spectrum.hpp"
#include "Spectrum2D.hpp"
#include "PlottingTools.h"


const bool _breakdownPlots = true;
const double targetPOT = 4.95e19;


using namespace std;


//____________________________________________________________________________________________________
void DrawPOT(double pot)
{
  //std::string str = "Simulated POT:" + std::to_string(pot);
  
  std::stringstream sstm;
  sstm << "Simulated POT: " << pot;
  std::string str = sstm.str();
  
  TLatex* pot_latex = new TLatex(.10, .96, str.c_str());
  pot_latex->SetTextColor(kGray+2);
  pot_latex->SetNDC();
  pot_latex->SetTextSize(1/30.);
  pot_latex->SetTextAlign(10); //left adjusted
  pot_latex->Draw();
  
  
  std::stringstream sstm2;
  sstm2 << "Scaled to POT: " << targetPOT;
  str = sstm2.str();
  
  TLatex* pot_latex_2 = new TLatex(.10, .92, str.c_str());
  pot_latex_2->SetTextFont(62);
  pot_latex_2->SetTextColor(kGray+2);
  pot_latex_2->SetNDC();
  pot_latex_2->SetTextSize(1/30.);
  pot_latex_2->SetTextAlign(10);//left adjusted
  pot_latex_2->Draw();
  
}


//____________________________________________________________________________________________________
void DrawTHStack(THStack *hs_trklen,
                 double pot_scaling,
                 std::map<std::string,TH1D*> themap);




//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
int main(int argc, char* argv[]) {
  
  clock_t begin = clock();
  
  std::string mc_file_name     = "ubxsecana_output.root";
  std::string bnbon_file_name  = "ubxsecana_output.root";
  std::string extbnb_file_name = "ubxsecana_output.root";
  double bnbon_pot_meas        = -1;
  
  
  int opt;
  while ((opt = getopt (argc, argv, "m:b:e:p:")) != -1)
  {
    switch (opt)
    {
      case 'm':
        std::cout << "MC file: " << optarg << std::endl;
        mc_file_name = optarg;
        break;
      case 'b':
        std::cout << "Data BNBON file: " << optarg << std::endl;
        bnbon_file_name = optarg;
        break;
      case 'e':
        std::cout << "Data EXTBNB file: " << optarg << std::endl;
        extbnb_file_name = optarg;
        break;
      case 'p':
        std::cout << "BNBON POT: " << optarg << std::endl;
        bnbon_pot_meas = atof(optarg);
        break;
      default:
        std::cerr << "Usage: " << argv[0] << " [-m mc file] [-b bnbon data file] [-e extbnb data file]" << std::endl;
        exit(EXIT_FAILURE);
    }
    
  }
  
  if (bnbon_pot_meas < 0) {
    std::cerr << "Provide POT for the BNBON sample." << std::endl;
    exit(EXIT_FAILURE);
  }
  

  // Reset input in order to use TApplication
  char* temp[1] = {argv[0]};
  argc = 1;
  argv = temp;

  
  TApplication* rootapp = new TApplication("ROOT Application",&argc, argv);
  //gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS
  gROOT->ProcessLine(".x rootlogon.C");

  gROOT->ProcessLine(".L loader_C.so");
  
  int bnbon_total_events = 1000;
  int extbnb_total_events = 1000;
  double mc_pot_sim = 6.e19;
  
  // *************************************
  // Opening files
  // *************************************
  TFile* mc_file = TFile::Open(mc_file_name.c_str(), "READ");
  TFile* bnbon_file = TFile::Open(bnbon_file_name.c_str(), "READ");
  TFile* extbnb_file = TFile::Open(extbnb_file_name.c_str(), "READ");
  
  
  // *************************************
  // Getting number of events for bnbon
  // *************************************
  TH1D* h_nevts_bnbon = (TH1D*)bnbon_file->Get("h_nevts");
  bnbon_total_events = h_nevts_bnbon->GetBinContent(1);
  std::cout << "Number of events (BNBON):  " << bnbon_total_events << std::endl;
  
  
  // *************************************
  // Getting number of events for extbnb
  // *************************************
  TH1D* h_nevts_extbnb = (TH1D*)extbnb_file->Get("h_nevts");
  extbnb_total_events = h_nevts_extbnb->GetBinContent(1);
  std::cout << "Number of events (EXTBNB): " << extbnb_total_events << std::endl << std::endl;
  
  
  // *************************************
  // Getting simulated POT
  // *************************************
  TH1D* h_simpot = (TH1D*)mc_file->Get("h_pot");
  mc_pot_sim = h_simpot->GetBinContent(1);
  std::cout << "Simulated POT:      " << mc_pot_sim << std::endl;
  std::cout << "BNBON Measured POT: " << bnbon_pot_meas << std::endl;
  
  // *************************************
  // Calculating scale factors
  // *************************************
  double scale_factor_extbnb = 1.23 * (382718./(double)extbnb_total_events) * ((double)bnbon_total_events/547616.);
  double scale_factor_bnbon = 1.; //bnbon_pot_meas * bnbon_pot_target;
  double scale_factor_mc = bnbon_pot_meas / mc_pot_sim;
  std::cout << "Data Scale Factors:" << std::endl;
  std::cout << "\t BNBON: " << scale_factor_bnbon << std::endl;
  std::cout << "\t EXTBNB: " << scale_factor_extbnb << std::endl;
  std::cout << "MC Scale Factors:" << std::endl;
  std::cout << "\t BNBCOSMIC: " << scale_factor_mc << std::endl;

  
  // *************************************
  // Getting the relevant histograms
  // *************************************
  std::map<std::string,TH1D*>* temp_map;
  mc_file->GetObject("hmap_trklen", temp_map);
  std::map<std::string,TH1D*> hmap_trklen_mc = *temp_map;
  //TH1D* h_trklen_total_mc = (TH1D*)mc_file->Get("h_trklen_total");
  TH1D* h_trklen_total_bnbon = (TH1D*)bnbon_file->Get("h_trklen_total");
  TH1D* h_trklen_total_extbnb = (TH1D*)extbnb_file->Get("h_trklen_total");
  
  
  // *************************************
  // Doing beam-on minus beam-off subctraction
  // *************************************
  h_trklen_total_extbnb->Scale(scale_factor_extbnb);
  h_trklen_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_trklen_data = (TH1D*)h_trklen_total_bnbon->Clone("h_trklen_data");
  h_trklen_total_bnbon->Add(h_trklen_total_extbnb, -1.);
  
  
  // *************************************
  // Plotting data and MC distribution
  // *************************************
  TCanvas* canvas_trklen = new TCanvas();
  THStack *hs_trklen_mc = new THStack("hs_trklen",";Candidate Track Length [cm]; Selected Events");
  DrawTHStack(hs_trklen_mc, scale_factor_mc, true, hmap_trklen_mc);
  h_trklen_data->Draw("E1 same");
  
  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  rootapp->Run();
  //rootapp->Terminate(0);
  
  return 0;
}








//**********************************************

void DrawTHStack(THStack *hs_trklen,
                 double pot_scaling,
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
  
}










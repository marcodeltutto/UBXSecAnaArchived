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
  std::cout << "BNBON Measured POT: " << bnbon_pot_meas << std::endl << std::endl;
  
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
  mc_file->GetObject("hmap_trktheta", temp_map);
  std::map<std::string,TH1D*> hmap_trktheta_mc = *temp_map;
  mc_file->GetObject("hmap_trkphi", temp_map);
  std::map<std::string,TH1D*> hmap_trkphi_mc = *temp_map;
  mc_file->GetObject("hmap_multpfp", temp_map);
  std::map<std::string,TH1D*> hmap_multpfp_mc = *temp_map;
  mc_file->GetObject("hmap_multtracktol", temp_map);
  std::map<std::string,TH1D*> hmap_multtracktol_mc = *temp_map;
  mc_file->GetObject("hmap_xdiff", temp_map);
  std::map<std::string,TH1D*> hmap_xdiff_mc = *temp_map;
  mc_file->GetObject("hmap_zdiff", temp_map);
  std::map<std::string,TH1D*> hmap_zdiff_mc = *temp_map;
  //TH1D* h_trklen_total_mc = (TH1D*)mc_file->Get("h_trklen_total");
  TH1D* h_trklen_total_bnbon = (TH1D*)bnbon_file->Get("h_trklen_total");
  TH1D* h_trklen_total_extbnb = (TH1D*)extbnb_file->Get("h_trklen_total");
  TH1D* h_trktheta_total_bnbon = (TH1D*)bnbon_file->Get("h_trktheta_total");
  TH1D* h_trktheta_total_extbnb = (TH1D*)extbnb_file->Get("h_trktheta_total");
  TH1D* h_trkphi_total_bnbon = (TH1D*)bnbon_file->Get("h_trkphi_total");
  TH1D* h_trkphi_total_extbnb = (TH1D*)extbnb_file->Get("h_trkphi_total");
  TH1D* h_multpfp_total_bnbon = (TH1D*)bnbon_file->Get("h_multpfp_total");
  TH1D* h_multpfp_total_extbnb = (TH1D*)extbnb_file->Get("h_multpfp_total");
  TH1D* h_multtracktol_total_bnbon = (TH1D*)bnbon_file->Get("h_multtracktol_total");
  TH1D* h_multtracktol_total_extbnb = (TH1D*)extbnb_file->Get("h_multtracktol_total");
  TH1D* h_xdiff_total_bnbon = (TH1D*)bnbon_file->Get("h_xdiff_total");
  TH1D* h_xdiff_total_extbnb = (TH1D*)extbnb_file->Get("h_xdiff_total");
  TH1D* h_zdiff_total_bnbon = (TH1D*)bnbon_file->Get("h_zdiff_total");
  TH1D* h_zdiff_total_extbnb = (TH1D*)extbnb_file->Get("h_zdiff_total");
  
  
  // *************************************
  // Doing beam-on minus beam-off subctraction
  // *************************************
  h_trklen_total_extbnb->Scale(scale_factor_extbnb);
  h_trklen_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_trklen_data = (TH1D*)h_trklen_total_bnbon->Clone("h_trklen_data");
  h_trklen_data->Sumw2();
  h_trklen_data->Add(h_trklen_total_extbnb, -1.);
  
  h_trktheta_total_extbnb->Scale(scale_factor_extbnb);
  h_trktheta_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_trktheta_data = (TH1D*)h_trktheta_total_bnbon->Clone("h_trktheta_data");
  h_trktheta_data->Sumw2();
  h_trktheta_data->Add(h_trktheta_total_extbnb, -1.);
  
  h_trkphi_total_extbnb->Scale(scale_factor_extbnb);
  h_trkphi_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_trkphi_data = (TH1D*)h_trkphi_total_bnbon->Clone("h_trklen_data");
  h_trkphi_data->Sumw2();
  h_trkphi_data->Add(h_trkphi_total_extbnb, -1.);
  
  h_multpfp_total_extbnb->Scale(scale_factor_extbnb);
  h_multpfp_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_multpfp_data = (TH1D*)h_multpfp_total_bnbon->Clone("h_multpfp_data");
  h_multpfp_data->Sumw2();
  h_multpfp_data->Add(h_multpfp_total_extbnb, -1.);
  
  h_multtracktol_total_extbnb->Scale(scale_factor_extbnb);
  h_multtracktol_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_multtracktol_data = (TH1D*)h_multtracktol_total_bnbon->Clone("h_multtracktol_data");
  h_multtracktol_data->Sumw2();
  h_multtracktol_data->Add(h_multtracktol_total_extbnb, -1.);
  
  h_xdiff_total_extbnb->Scale(scale_factor_extbnb);
  h_xdiff_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_xdiff_data = (TH1D*)h_xdiff_total_bnbon->Clone("h_xdiff_data");
  h_xdiff_data->Sumw2();
  h_xdiff_data->Add(h_xdiff_total_extbnb, -1.);
  
  h_zdiff_total_extbnb->Scale(scale_factor_extbnb);
  h_zdiff_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_zdiff_data = (TH1D*)h_zdiff_total_bnbon->Clone("h_zdiff_data");
  h_zdiff_data->Sumw2();
  h_zdiff_data->Add(h_zdiff_total_extbnb, -1.);
  
  
  // *************************************
  // Plotting data and MC distribution
  // *************************************
  TCanvas* canvas_trklen = new TCanvas();
  THStack *hs_trklen_mc = new THStack("hs_trklen",";Candidate Track Length [cm]; Selected Events");
  std::cout << "before hmap_trklen_mc[\"total\"]->Integral()" << hmap_trklen_mc["total"]->Integral() << std::endl;
  DrawTHStack(hs_trklen_mc, scale_factor_mc, true, hmap_trklen_mc);
  h_trklen_data->Draw("E1 same");
  std::cout << "h_trklen_data->Integral()" << h_trklen_data->Integral() << std::endl;
  std::cout << "after hmap_trklen_mc[\"total\"]->Integral()" << hmap_trklen_mc["total"]->Integral() << std::endl;

  
  TCanvas* canvas_trktheta = new TCanvas();
  THStack *hs_trktheta_mc = new THStack("hs_trktheta",";Candidate Track cos(#theta); Selected Events");
  DrawTHStack(hs_trktheta_mc, scale_factor_mc, true, hmap_trktheta_mc);
  h_trktheta_data->Draw("E1 same");
  
  TCanvas* canvas_trkphi = new TCanvas();
  THStack *hs_trkphi_mc = new THStack("hs_trkphi",";Candidate Track #phi; Selected Events");
  DrawTHStack(hs_trkphi_mc, scale_factor_mc, true, hmap_trkphi_mc);
  h_trkphi_data->Draw("E1 same");
  
  TCanvas* canvas_multpfp = new TCanvas();
  THStack *hs_multpfp_mc = new THStack("hs_multpfp",";PFP Multiplicity; Selected Events");
  DrawTHStack(hs_multpfp_mc, scale_factor_mc, true, hmap_multpfp_mc);
  h_multpfp_data->Draw("E1 same");
  
  TCanvas* canvas_multtracktol = new TCanvas();
  THStack *hs_multtracktol_mc = new THStack("hs_multtracktol",";Track Multiplicity (5 cm); Selected Events");
  DrawTHStack(hs_multtracktol_mc, scale_factor_mc, true, hmap_multtracktol_mc);
  h_multtracktol_data->Draw("E1 same");
  
  TCanvas* canvas_xdiff = new TCanvas();
  THStack *hs_xdiff_mc = new THStack("hs_xdiff",";xdiff; TPCObjects");
  DrawTHStack2(hs_xdiff_mc, scale_factor_mc, true, hmap_xdiff_mc);
  h_xdiff_data->Draw("E1 same");
  std::cout << "h_xdiff_data->Integral()" << h_xdiff_data->Integral() << std::endl;
  std::cout << "hmap_xdiff_mc[\"total\"]->Integral()" << hmap_xdiff_mc["total"]->Integral() << std::endl;

  TCanvas* canvas_zdiff = new TCanvas();
  THStack *hs_zdiff_mc = new THStack("hs_zdiff",";zdiff; TPCObjects");
  DrawTHStack2(hs_zdiff_mc, scale_factor_mc, true, hmap_zdiff_mc);
  h_zdiff_data->Draw("E1 same");
  std::cout << "h_zdiff_data->Integral()" << h_zdiff_data->Integral() << std::endl;
  std::cout << "hmap_zdiff_mc[\"total\"]->Integral()" << hmap_zdiff_mc["total"]->Integral() << std::endl;

  
  // *************************************
  // Other data/MC distributions
  // *************************************
  TH1D* h_flsTime_mc = (TH1D*)mc_file->Get("h_flsTime_wcut");
  TH1D* h_flsTime_bnbon = (TH1D*)bnbon_file->Get("h_flsTime_wcut");
  TH1D* h_flsTime_extbnb = (TH1D*)extbnb_file->Get("h_flsTime_wcut");
  h_flsTime_extbnb->Scale(scale_factor_extbnb);
  h_flsTime_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_flsTime_data = (TH1D*)h_flsTime_bnbon->Clone("h_flsTime_data");
  h_flsTime_data->Sumw2();
  h_flsTime_data->Add(h_flsTime_extbnb, -1.);
  TCanvas* canvas_flsTime = new TCanvas();
  h_flsTime_mc->Draw("histo");
  h_flsTime_data->SetLineColor(kRed);
  h_flsTime_data->Draw("E1 same");
  new TCanvas();
  h_flsTime_mc->SetLineColor(kBlack);
  h_flsTime_bnbon->SetLineColor(kBlue);
  h_flsTime_extbnb->SetLineColor(kRed);
  h_flsTime_mc->Draw("histo");
  h_flsTime_bnbon->Draw("E1 same");
  h_flsTime_extbnb->Draw("E1 same");


  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  rootapp->Run();
  //rootapp->Terminate(0);
  
  return 0;
}









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
  mc_file->GetObject("hmap_xdiff_b", temp_map);
  std::map<std::string,TH1D*> hmap_xdiff_b_mc = *temp_map;
  mc_file->GetObject("hmap_zdiff_b", temp_map);
  std::map<std::string,TH1D*> hmap_zdiff_b_mc = *temp_map;
  mc_file->GetObject("hmap_xdiff", temp_map);
  std::map<std::string,TH1D*> hmap_xdiff_mc = *temp_map;
  mc_file->GetObject("hmap_zdiff", temp_map);
  std::map<std::string,TH1D*> hmap_zdiff_mc = *temp_map;
  mc_file->GetObject("hmap_vtxx", temp_map);
  std::map<std::string,TH1D*> hmap_vtxx_mc = *temp_map;
  mc_file->GetObject("hmap_vtxy", temp_map);
  std::map<std::string,TH1D*> hmap_vtxy_mc = *temp_map;
  mc_file->GetObject("hmap_vtxz", temp_map);
  std::map<std::string,TH1D*> hmap_vtxz_mc = *temp_map;
  mc_file->GetObject("hmap_flsmatch_score", temp_map);
  std::map<std::string,TH1D*> hmap_flsmatch_score_mc = *temp_map;
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
  TH1D* h_xdiff_b_total_bnbon = (TH1D*)bnbon_file->Get("h_xdiff_total_b");
  TH1D* h_xdiff_b_total_extbnb = (TH1D*)extbnb_file->Get("h_xdiff_total_b");
  TH1D* h_zdiff_b_total_bnbon = (TH1D*)bnbon_file->Get("h_zdiff_total_b");
  TH1D* h_zdiff_b_total_extbnb = (TH1D*)extbnb_file->Get("h_zdiff_total_b");
  TH1D* h_xdiff_total_bnbon = (TH1D*)bnbon_file->Get("h_xdiff_total");
  TH1D* h_xdiff_total_extbnb = (TH1D*)extbnb_file->Get("h_xdiff_total");
  TH1D* h_zdiff_total_bnbon = (TH1D*)bnbon_file->Get("h_zdiff_total");
  TH1D* h_zdiff_total_extbnb = (TH1D*)extbnb_file->Get("h_zdiff_total");
  TH1D* h_vtxx_total_bnbon = (TH1D*)bnbon_file->Get("h_vtxx_total");
  TH1D* h_vtxx_total_extbnb = (TH1D*)extbnb_file->Get("h_vtxx_total");
  TH1D* h_vtxy_total_bnbon = (TH1D*)bnbon_file->Get("h_vtxy_total");
  TH1D* h_vtxy_total_extbnb = (TH1D*)extbnb_file->Get("h_vtxy_total");
  TH1D* h_vtxz_total_bnbon = (TH1D*)bnbon_file->Get("h_vtxz_total");
  TH1D* h_vtxz_total_extbnb = (TH1D*)extbnb_file->Get("h_vtxz_total");
  TH1D* h_flsmatch_score_total_bnbon = (TH1D*)bnbon_file->Get("h_flsmatch_score_total");
  TH1D* h_flsmatch_score_total_extbnb = (TH1D*)extbnb_file->Get("h_flsmatch_score_total");
  
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
  
  h_xdiff_b_total_extbnb->Scale(scale_factor_extbnb);
  h_xdiff_b_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_xdiff_b_data = (TH1D*)h_xdiff_b_total_bnbon->Clone("h_xdiff_b_data");
  h_xdiff_b_data->Sumw2();
  h_xdiff_b_data->Add(h_xdiff_b_total_extbnb, -1.);
  
  h_zdiff_b_total_extbnb->Scale(scale_factor_extbnb);
  h_zdiff_b_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_zdiff_b_data = (TH1D*)h_zdiff_b_total_bnbon->Clone("h_zdiff_b_data");
  h_zdiff_b_data->Sumw2();
  h_zdiff_b_data->Add(h_zdiff_b_total_extbnb, -1.);
  
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
  
  h_vtxx_total_extbnb->Scale(scale_factor_extbnb);
  h_vtxx_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_vtxx_data = (TH1D*)h_vtxx_total_bnbon->Clone("h_vtxx_data");
  h_vtxx_data->Sumw2();
  h_vtxx_data->Add(h_vtxx_total_extbnb, -1.);
  
  h_vtxy_total_extbnb->Scale(scale_factor_extbnb);
  h_vtxy_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_vtxy_data = (TH1D*)h_vtxy_total_bnbon->Clone("h_vtxx_data");
  h_vtxy_data->Sumw2();
  h_vtxy_data->Add(h_vtxy_total_extbnb, -1.);
  
  h_vtxz_total_extbnb->Scale(scale_factor_extbnb);
  h_vtxz_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_vtxz_data = (TH1D*)h_vtxz_total_bnbon->Clone("h_vtxx_data");
  h_vtxz_data->Sumw2();
  h_vtxz_data->Add(h_vtxz_total_extbnb, -1.);
  
  h_flsmatch_score_total_extbnb->Scale(scale_factor_extbnb);
  h_flsmatch_score_total_bnbon->Scale(scale_factor_bnbon);
  TH1D* h_flsmatch_score_data = (TH1D*)h_flsmatch_score_total_bnbon->Clone("h_flsmatch_score_data");
  h_flsmatch_score_data->Sumw2();
  h_flsmatch_score_data->Add(h_flsmatch_score_total_extbnb, -1.);
  
  // *************************************
  // Plotting data and MC distribution
  // *************************************
  TLegend* leg;
  TCanvas* canvas_trklen = new TCanvas();
  THStack *hs_trklen_mc = new THStack("hs_trklen",";Candidate Track Length [cm]; Selected Events");
  std::cout << "before hmap_trklen_mc[\"total\"]->Integral()" << hmap_trklen_mc["total"]->Integral() << std::endl;
  leg = DrawTHStack(hs_trklen_mc, scale_factor_mc, true, hmap_trklen_mc);
  DrawDataHisto(h_trklen_data);
  leg->AddEntry(h_trklen_data,"Data (Beam-on - Beam-off)","lep");
  std::cout << "h_trklen_data->Integral()" << h_trklen_data->Integral() << std::endl;
  std::cout << "after hmap_trklen_mc[\"total\"]->Integral()" << hmap_trklen_mc["total"]->Integral() << std::endl;
  DrawPOT(bnbon_pot_meas);

  TCanvas* canvas_trktheta = new TCanvas();
  THStack *hs_trktheta_mc = new THStack("hs_trktheta",";Candidate Track cos(#theta); Selected Events");
  leg = DrawTHStack(hs_trktheta_mc, scale_factor_mc, true, hmap_trktheta_mc);
  DrawDataHisto(h_trktheta_data);
  leg->AddEntry(h_trktheta_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);

  TCanvas* canvas_trkphi = new TCanvas();
  THStack *hs_trkphi_mc = new THStack("hs_trkphi",";Candidate Track #phi; Selected Events");
  leg = DrawTHStack(hs_trkphi_mc, scale_factor_mc, true, hmap_trkphi_mc);
  DrawDataHisto(h_trkphi_data);
  leg->AddEntry(h_trkphi_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  TCanvas* canvas_multpfp = new TCanvas();
  THStack *hs_multpfp_mc = new THStack("hs_multpfp",";PFP Multiplicity; Selected Events");
  leg = DrawTHStack(hs_multpfp_mc, scale_factor_mc, true, hmap_multpfp_mc);
  DrawDataHisto(h_multpfp_data);
  leg->AddEntry(h_multpfp_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  TCanvas* canvas_multtracktol = new TCanvas();
  THStack *hs_multtracktol_mc = new THStack("hs_multtracktol",";Track Multiplicity (5 cm); Selected Events");
  leg = DrawTHStack(hs_multtracktol_mc, scale_factor_mc, true, hmap_multtracktol_mc);
  DrawDataHisto(h_multtracktol_data);
  leg->AddEntry(h_multtracktol_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  TCanvas* canvas_xdiff_b = new TCanvas();
  THStack *hs_xdiff_b_mc = new THStack("hs_xdiff_b",";QLL x - TPC x [cm]; TPCObjects (Before Selection)");
  leg = DrawTHStack2(hs_xdiff_b_mc, scale_factor_mc, true, hmap_xdiff_b_mc);
  DrawDataHisto(h_xdiff_b_data);
  leg->AddEntry(h_xdiff_b_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  TCanvas* canvas_b_zdiff = new TCanvas();
  THStack *hs_zdiff_b_mc = new THStack("hs_zdiff_b",";Hypo z - Flash z [cm]; TPCObjects (Before Selection)");
  leg = DrawTHStack2(hs_zdiff_b_mc, scale_factor_mc, true, hmap_zdiff_b_mc);
  DrawDataHisto(h_zdiff_b_data);
  leg->AddEntry(h_zdiff_b_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  TCanvas* canvas_xdiff = new TCanvas();
  THStack *hs_xdiff_mc = new THStack("hs_xdiff",";QLL x - TPC x [cm]; Selected Events");
  leg = DrawTHStack2(hs_xdiff_mc, scale_factor_mc, true, hmap_xdiff_mc);
  DrawDataHisto(h_xdiff_data);
  leg->AddEntry(h_xdiff_data,"Data (Beam-on - Beam-off)","lep");
  std::cout << "h_xdiff_data->Integral()" << h_xdiff_data->Integral() << std::endl;
  std::cout << "hmap_xdiff_mc[\"total\"]->Integral()" << hmap_xdiff_mc["total"]->Integral() << std::endl;
  DrawPOT(bnbon_pot_meas);

  TCanvas* canvas_zdiff = new TCanvas();
  THStack *hs_zdiff_mc = new THStack("hs_zdiff",";Hypo z - Flash z [cm]; Selected Events");
  leg = DrawTHStack2(hs_zdiff_mc, scale_factor_mc, true, hmap_zdiff_mc);
  DrawDataHisto(h_zdiff_data);
  leg->AddEntry(h_zdiff_data,"Data (Beam-on - Beam-off)","lep");
  std::cout << "h_zdiff_data->Integral()" << h_zdiff_data->Integral() << std::endl;
  std::cout << "hmap_zdiff_mc[\"total\"]->Integral()" << hmap_zdiff_mc["total"]->Integral() << std::endl;
  DrawPOT(bnbon_pot_meas);

  TCanvas* canvas_vtxx = new TCanvas();
  THStack *hs_vtxx_mc = new THStack("hs_vtxx",";Candidate Neutrino Vertex X [cm]; Selected Events");
  leg = DrawTHStack2(hs_vtxx_mc, scale_factor_mc, true, hmap_vtxx_mc);
  DrawDataHisto(h_vtxx_data);
  hs_vtxx_mc->SetMaximum(600);
  leg->AddEntry(h_vtxx_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  TCanvas* canvas_vtxy = new TCanvas();
  THStack *hs_vtxy_mc = new THStack("hs_vtxy",";Candidate Neutrino Vertex Y [cm]; Selected Events");
  leg = DrawTHStack2(hs_vtxy_mc, scale_factor_mc, true, hmap_vtxy_mc);
  hs_vtxy_mc->SetMaximum(650);
  DrawDataHisto(h_vtxy_data);
  leg->AddEntry(h_vtxy_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  TCanvas* canvas_vtxz = new TCanvas();
  THStack *hs_vtxz_mc = new THStack("hs_vtxz",";Candidate Neutrino Vertex Z [cm]; Selected Events");
  leg = DrawTHStack2(hs_vtxz_mc, scale_factor_mc, true, hmap_vtxz_mc);
  //h_vtxz_data->Draw("E1 same");
  DrawDataHisto(h_vtxz_data);
  hs_vtxz_mc->SetMaximum(900);
  leg->AddEntry(h_vtxz_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  TCanvas* canvas_flsmatch_score = new TCanvas();
  THStack *hs_flsmatch_score_mc = new THStack("hs_flsmatch_score",";1/(-log(L)); Selected Events");
  leg = DrawTHStack2(hs_flsmatch_score_mc, scale_factor_mc, true, hmap_flsmatch_score_mc);
  //h_vtxz_data->Draw("E1 same");
  DrawDataHisto(h_flsmatch_score_data);
  //hs_flsmatch_score_mc->SetMaximum(900);
  leg->AddEntry(h_flsmatch_score_data,"Data (Beam-on - Beam-off)","lep");
  DrawPOT(bnbon_pot_meas);
  
  
  
  // *************************************
  // Other data/MC distributions
  // *************************************
  TH1D* h_flsTime_mc = (TH1D*)mc_file->Get("h_flsTime_wcut");
  TH1D* h_flsTime_bnbon = (TH1D*)bnbon_file->Get("h_flsTime_wcut");
  TH1D* h_flsTime_extbnb = (TH1D*)extbnb_file->Get("h_flsTime_wcut");
  h_flsTime_extbnb->Scale(scale_factor_extbnb);
  h_flsTime_bnbon->Scale(scale_factor_bnbon);
  h_flsTime_mc->Scale(scale_factor_mc);
  TH1D* h_flsTime_data = (TH1D*)h_flsTime_bnbon->Clone("h_flsTime_data");
  h_flsTime_data->Sumw2();
  h_flsTime_data->Add(h_flsTime_extbnb, -1.);
  
  new TCanvas();
  h_flsTime_mc->Draw("histo");
  h_flsTime_data->SetLineColor(kRed);
  h_flsTime_data->SetMarkerColor(kRed);
  DrawDataHisto(h_flsTime_data);
  TLegend* leg2;
  leg2 = new TLegend(0.13,0.69,0.45,0.87,NULL,"brNDC");
  leg2->AddEntry(h_flsTime_mc,"MC BNB+COSMIC","l");
  leg2->AddEntry(h_flsTime_data,"Data (Beam-on - Beam-off)","lep");
  leg2->Draw();
  DrawPOT(bnbon_pot_meas);

  new TCanvas();
  h_flsTime_mc->SetLineColor(kBlack);
  h_flsTime_bnbon->SetLineColor(kRed);
  h_flsTime_extbnb->SetLineColor(kBlue);
  h_flsTime_mc->Draw("histo");
  h_flsTime_bnbon->Draw("histo same");
  h_flsTime_extbnb->Draw("histo same");
  TLegend* leg3;
  leg3 = new TLegend(0.13,0.69,0.45,0.87,NULL,"brNDC");
  leg3->AddEntry(h_flsTime_mc,"MC BNB+COSMIC","l");
  leg3->AddEntry(h_flsTime_bnbon,"DATA BNBON","l");
  leg3->AddEntry(h_flsTime_extbnb,"DATA EXTBNB","l");
  leg3->Draw();
  DrawPOT(bnbon_pot_meas);

  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << endl << endl << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  rootapp->Run();
  //rootapp->Terminate(0);
  
  return 0;
}









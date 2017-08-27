//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Aug 26 15:56:11 2017 by ROOT version 6.06/06
// from TTree tree/
// found on file: ../files/ubxsec_output_mc_7.root
//////////////////////////////////////////////////////////

#ifndef AnaTree_h
#define AnaTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using namespace std;

class AnaTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Int_t           muon_is_reco;
   Double_t        muon_reco_pur;
   Double_t        muon_reco_eff;
   Double_t        true_muon_mom;
   Double_t        true_muon_mom_matched;
   Int_t           nPFPtagged;
   Int_t           muon_is_flash_tagged;
   Double_t        muon_tag_score;
   Double_t        fm_score;
   Int_t           fv;
   Int_t           ccnc;
   Int_t           nupdg;
   Bool_t          is_signal;
   Double_t        nu_e;
   Double_t        recon_muon_start_x;
   Double_t        recon_muon_start_y;
   Double_t        recon_muon_start_z;
   Double_t        recon_muon_end_x;
   Double_t        recon_muon_end_y;
   Double_t        recon_muon_end_z;
   Double_t        mc_muon_start_x;
   Double_t        mc_muon_start_y;
   Double_t        mc_muon_start_z;
   Double_t        mc_muon_end_x;
   Double_t        mc_muon_end_y;
   Double_t        mc_muon_end_z;
   Int_t           mc_muon_contained;
   Int_t           is_swtriggered;
   Double_t        vtx_resolution;
   Int_t           nslices;
   vector<double>  *slc_flsmatch_score;
   vector<double>  *slc_flsmatch_qllx;
   vector<double>  *slc_flsmatch_tpcx;
   vector<double>  *slc_flsmatch_t0;
   vector<double>  *slc_flsmatch_hypoz;
   vector<double>  *slc_flsmatch_xfixed_chi2;
   vector<double>  *slc_flsmatch_xfixed_ll;
   vector<double>  *slc_flsmatch_cosmic_score;
   vector<double>  *slc_flsmatch_cosmic_t0;
   vector<double>  *slc_nuvtx_x;
   vector<double>  *slc_nuvtx_y;
   vector<double>  *slc_nuvtx_z;
   vector<int>     *slc_nuvtx_fv;
   vector<double>  *slc_vtxcheck_angle;
   vector<int>     *slc_origin;
   vector<int>     *slc_nhits_u;
   vector<int>     *slc_nhits_v;
   vector<int>     *slc_nhits_w;
   vector<double>  *slc_longesttrack_length;
   vector<double>  *slc_longesttrack_phi;
   vector<double>  *slc_longesttrack_theta;
   vector<bool>    *slc_longesttrack_iscontained;
   vector<int>     *slc_acpt_outoftime;
   vector<int>     *slc_crosses_top_boundary;
   vector<int>     *slc_nuvtx_closetodeadregion_u;
   vector<int>     *slc_nuvtx_closetodeadregion_v;
   vector<int>     *slc_nuvtx_closetodeadregion_w;
   vector<double>  *slc_kalman_chi2;
   vector<int>     *slc_kalman_ndof;
   vector<bool>    *slc_passed_min_track_quality;
   vector<double>  *slc_n_intime_pe_closestpmt;
   vector<double>  *slc_maxdistance_vtxtrack;
   vector<int>     *slc_npfp;
   vector<int>     *slc_ntrack;
   vector<int>     *slc_nshower;
   vector<bool>    *slc_iscontained;
   vector<int>     *slc_mult_pfp;
   vector<int>     *slc_mult_track;
   vector<int>     *slc_mult_shower;
   vector<int>     *slc_mult_track_tolerance;
   Int_t           nbeamfls;
   vector<double>  *beamfls_time;
   vector<double>  *beamfls_pe;
   vector<double>  *beamfls_z;
   Bool_t          no_mcflash_but_op_activity;
   vector<vector<double> > *beamfls_spec;
   vector<double>  *numc_flash_spec;
   vector<vector<double> > *slc_flshypo_xfixed_spec;
   vector<vector<double> > *slc_flshypo_spec;
   Int_t           nsignal;
   vector<double>  *mctrk_start_x;
   vector<double>  *mctrk_start_y;
   vector<double>  *mctrk_start_z;
   vector<double>  *trk_start_x;
   vector<double>  *trk_start_y;
   vector<double>  *trk_start_z;
   vector<double>  *vtx_x;
   vector<double>  *vtx_y;
   vector<double>  *vtx_z;
   vector<double>  *tvtx_x;
   vector<double>  *tvtx_y;
   vector<double>  *tvtx_z;
   Double_t        pot;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_muon_is_reco;   //!
   TBranch        *b_muon_reco_pur;   //!
   TBranch        *b_muon_reco_eff;   //!
   TBranch        *b_true_muon_mom;   //!
   TBranch        *b_true_muon_mom_matched;   //!
   TBranch        *b_nPFPtagged;   //!
   TBranch        *b_muon_is_flash_tagged;   //!
   TBranch        *b_muon_tag_score;   //!
   TBranch        *b_fm_score;   //!
   TBranch        *b_fv;   //!
   TBranch        *b_ccnc;   //!
   TBranch        *b_nupdg;   //!
   TBranch        *b_is_signal;   //!
   TBranch        *b_nu_e;   //!
   TBranch        *b_recon_muon_start_x;   //!
   TBranch        *b_recon_muon_start_y;   //!
   TBranch        *b_recon_muon_start_z;   //!
   TBranch        *b_recon_muon_end_x;   //!
   TBranch        *b_recon_muon_end_y;   //!
   TBranch        *b_recon_muon_end_z;   //!
   TBranch        *b_mc_muon_start_x;   //!
   TBranch        *b_mc_muon_start_y;   //!
   TBranch        *b_mc_muon_start_z;   //!
   TBranch        *b_mc_muon_end_x;   //!
   TBranch        *b_mc_muon_end_y;   //!
   TBranch        *b_mc_muon_end_z;   //!
   TBranch        *b_mc_muon_contained;   //!
   TBranch        *b_is_swtriggered;   //!
   TBranch        *b_vtx_resolution;   //!
   TBranch        *b_nslices;   //!
   TBranch        *b_slc_flsmatch_score;   //!
   TBranch        *b_slc_flsmatch_qllx;   //!
   TBranch        *b_slc_flsmatch_tpcx;   //!
   TBranch        *b_slc_flsmatch_t0;   //!
   TBranch        *b_slc_flsmatch_hypoz;   //!
   TBranch        *b_slc_flsmatch_xfixed_chi2;   //!
   TBranch        *b_slc_flsmatch_xfixed_ll;   //!
   TBranch        *b_slc_flsmatch_cosmic_score;   //!
   TBranch        *b_slc_flsmatch_cosmic_t0;   //!
   TBranch        *b_slc_nuvtx_x;   //!
   TBranch        *b_slc_nuvtx_y;   //!
   TBranch        *b_slc_nuvtx_z;   //!
   TBranch        *b_slc_nuvtx_fv;   //!
   TBranch        *b_slc_vtxcheck_angle;   //!
   TBranch        *b_slc_origin;   //!
   TBranch        *b_slc_nhits_u;   //!
   TBranch        *b_slc_nhits_v;   //!
   TBranch        *b_slc_nhits_w;   //!
   TBranch        *b_slc_longesttrack_length;   //!
   TBranch        *b_slc_longesttrack_phi;   //!
   TBranch        *b_slc_longesttrack_theta;   //!
   TBranch        *b_slc_longesttrack_iscontained;   //!
   TBranch        *b_slc_acpt_outoftime;   //!
   TBranch        *b_slc_crosses_top_boundary;   //!
   TBranch        *b_slc_nuvtx_closetodeadregion_u;   //!
   TBranch        *b_slc_nuvtx_closetodeadregion_v;   //!
   TBranch        *b_slc_nuvtx_closetodeadregion_w;   //!
   TBranch        *b_slc_kalman_chi2;   //!
   TBranch        *b_slc_kalman_ndof;   //!
   TBranch        *b_slc_passed_min_track_quality;   //!
   TBranch        *b_slc_n_intime_pe_closestpmt;   //!
   TBranch        *b_slc_maxdistance_vtxtrack;   //!
   TBranch        *b_slc_npfp;   //!
   TBranch        *b_slc_ntrack;   //!
   TBranch        *b_slc_nshower;   //!
   TBranch        *b_slc_iscontained;   //!
   TBranch        *b_slc_mult_pfp;   //!
   TBranch        *b_slc_mult_track;   //!
   TBranch        *b_slc_mult_shower;   //!
   TBranch        *b_slc_mult_track_tolerance;   //!
   TBranch        *b_nbeamfls;   //!
   TBranch        *b_beamfls_time;   //!
   TBranch        *b_beamfls_pe;   //!
   TBranch        *b_beamfls_z;   //!
   TBranch        *b_no_mcflash_but_op_activity;   //!
   TBranch        *b_beamfls_spec;   //!
   TBranch        *b_numc_flash_spec;   //!
   TBranch        *b_slc_flshypo_xfixed_spec;   //!
   TBranch        *b_slc_flshypo_spec;   //!
   TBranch        *b_nsignal;   //!
   TBranch        *b_mctrk_start_x;   //!
   TBranch        *b_mctrk_start_y;   //!
   TBranch        *b_mctrk_start_z;   //!
   TBranch        *b_trk_start_x;   //!
   TBranch        *b_trk_start_y;   //!
   TBranch        *b_trk_start_z;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_tvtx_x;   //!
   TBranch        *b_tvtx_y;   //!
   TBranch        *b_tvtx_z;   //!
   TBranch        *b_pot;   //!

   AnaTree(TTree *tree=0);
   virtual ~AnaTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnaTree_cxx
AnaTree::AnaTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../files/ubxsec_output_mc_7.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../files/ubxsec_output_mc_7.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../files/ubxsec_output_mc_7.root:/UBXSec");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

AnaTree::~AnaTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnaTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   slc_flsmatch_score = 0;
   slc_flsmatch_qllx = 0;
   slc_flsmatch_tpcx = 0;
   slc_flsmatch_t0 = 0;
   slc_flsmatch_hypoz = 0;
   slc_flsmatch_xfixed_chi2 = 0;
   slc_flsmatch_xfixed_ll = 0;
   slc_flsmatch_cosmic_score = 0;
   slc_flsmatch_cosmic_t0 = 0;
   slc_nuvtx_x = 0;
   slc_nuvtx_y = 0;
   slc_nuvtx_z = 0;
   slc_nuvtx_fv = 0;
   slc_vtxcheck_angle = 0;
   slc_origin = 0;
   slc_nhits_u = 0;
   slc_nhits_v = 0;
   slc_nhits_w = 0;
   slc_longesttrack_length = 0;
   slc_longesttrack_phi = 0;
   slc_longesttrack_theta = 0;
   slc_longesttrack_iscontained = 0;
   slc_acpt_outoftime = 0;
   slc_crosses_top_boundary = 0;
   slc_nuvtx_closetodeadregion_u = 0;
   slc_nuvtx_closetodeadregion_v = 0;
   slc_nuvtx_closetodeadregion_w = 0;
   slc_kalman_chi2 = 0;
   slc_kalman_ndof = 0;
   slc_passed_min_track_quality = 0;
   slc_n_intime_pe_closestpmt = 0;
   slc_maxdistance_vtxtrack = 0;
   slc_npfp = 0;
   slc_ntrack = 0;
   slc_nshower = 0;
   slc_iscontained = 0;
   slc_mult_pfp = 0;
   slc_mult_track = 0;
   slc_mult_shower = 0;
   slc_mult_track_tolerance = 0;
   beamfls_time = 0;
   beamfls_pe = 0;
   beamfls_z = 0;
   beamfls_spec = 0;
   numc_flash_spec = 0;
   slc_flshypo_xfixed_spec = 0;
   slc_flshypo_spec = 0;
   mctrk_start_x = 0;
   mctrk_start_y = 0;
   mctrk_start_z = 0;
   trk_start_x = 0;
   trk_start_y = 0;
   trk_start_z = 0;
   vtx_x = 0;
   vtx_y = 0;
   vtx_z = 0;
   tvtx_x = 0;
   tvtx_y = 0;
   tvtx_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("muon_is_reco", &muon_is_reco, &b_muon_is_reco);
   fChain->SetBranchAddress("muon_reco_pur", &muon_reco_pur, &b_muon_reco_pur);
   fChain->SetBranchAddress("muon_reco_eff", &muon_reco_eff, &b_muon_reco_eff);
   fChain->SetBranchAddress("true_muon_mom", &true_muon_mom, &b_true_muon_mom);
   fChain->SetBranchAddress("true_muon_mom_matched", &true_muon_mom_matched, &b_true_muon_mom_matched);
   fChain->SetBranchAddress("nPFPtagged", &nPFPtagged, &b_nPFPtagged);
   fChain->SetBranchAddress("muon_is_flash_tagged", &muon_is_flash_tagged, &b_muon_is_flash_tagged);
   fChain->SetBranchAddress("muon_tag_score", &muon_tag_score, &b_muon_tag_score);
   fChain->SetBranchAddress("fm_score", &fm_score, &b_fm_score);
   fChain->SetBranchAddress("fv", &fv, &b_fv);
   fChain->SetBranchAddress("ccnc", &ccnc, &b_ccnc);
   fChain->SetBranchAddress("nupdg", &nupdg, &b_nupdg);
   fChain->SetBranchAddress("is_signal", &is_signal, &b_is_signal);
   fChain->SetBranchAddress("nu_e", &nu_e, &b_nu_e);
   fChain->SetBranchAddress("recon_muon_start_x", &recon_muon_start_x, &b_recon_muon_start_x);
   fChain->SetBranchAddress("recon_muon_start_y", &recon_muon_start_y, &b_recon_muon_start_y);
   fChain->SetBranchAddress("recon_muon_start_z", &recon_muon_start_z, &b_recon_muon_start_z);
   fChain->SetBranchAddress("recon_muon_end_x", &recon_muon_end_x, &b_recon_muon_end_x);
   fChain->SetBranchAddress("recon_muon_end_y", &recon_muon_end_y, &b_recon_muon_end_y);
   fChain->SetBranchAddress("recon_muon_end_z", &recon_muon_end_z, &b_recon_muon_end_z);
   fChain->SetBranchAddress("mc_muon_start_x", &mc_muon_start_x, &b_mc_muon_start_x);
   fChain->SetBranchAddress("mc_muon_start_y", &mc_muon_start_y, &b_mc_muon_start_y);
   fChain->SetBranchAddress("mc_muon_start_z", &mc_muon_start_z, &b_mc_muon_start_z);
   fChain->SetBranchAddress("mc_muon_end_x", &mc_muon_end_x, &b_mc_muon_end_x);
   fChain->SetBranchAddress("mc_muon_end_y", &mc_muon_end_y, &b_mc_muon_end_y);
   fChain->SetBranchAddress("mc_muon_end_z", &mc_muon_end_z, &b_mc_muon_end_z);
   fChain->SetBranchAddress("mc_muon_contained", &mc_muon_contained, &b_mc_muon_contained);
   fChain->SetBranchAddress("is_swtriggered", &is_swtriggered, &b_is_swtriggered);
   fChain->SetBranchAddress("vtx_resolution", &vtx_resolution, &b_vtx_resolution);
   fChain->SetBranchAddress("nslices", &nslices, &b_nslices);
   fChain->SetBranchAddress("slc_flsmatch_score", &slc_flsmatch_score, &b_slc_flsmatch_score);
   fChain->SetBranchAddress("slc_flsmatch_qllx", &slc_flsmatch_qllx, &b_slc_flsmatch_qllx);
   fChain->SetBranchAddress("slc_flsmatch_tpcx", &slc_flsmatch_tpcx, &b_slc_flsmatch_tpcx);
   fChain->SetBranchAddress("slc_flsmatch_t0", &slc_flsmatch_t0, &b_slc_flsmatch_t0);
   fChain->SetBranchAddress("slc_flsmatch_hypoz", &slc_flsmatch_hypoz, &b_slc_flsmatch_hypoz);
   fChain->SetBranchAddress("slc_flsmatch_xfixed_chi2", &slc_flsmatch_xfixed_chi2, &b_slc_flsmatch_xfixed_chi2);
   fChain->SetBranchAddress("slc_flsmatch_xfixed_ll", &slc_flsmatch_xfixed_ll, &b_slc_flsmatch_xfixed_ll);
   fChain->SetBranchAddress("slc_flsmatch_cosmic_score", &slc_flsmatch_cosmic_score, &b_slc_flsmatch_cosmic_score);
   fChain->SetBranchAddress("slc_flsmatch_cosmic_t0", &slc_flsmatch_cosmic_t0, &b_slc_flsmatch_cosmic_t0);
   fChain->SetBranchAddress("slc_nuvtx_x", &slc_nuvtx_x, &b_slc_nuvtx_x);
   fChain->SetBranchAddress("slc_nuvtx_y", &slc_nuvtx_y, &b_slc_nuvtx_y);
   fChain->SetBranchAddress("slc_nuvtx_z", &slc_nuvtx_z, &b_slc_nuvtx_z);
   fChain->SetBranchAddress("slc_nuvtx_fv", &slc_nuvtx_fv, &b_slc_nuvtx_fv);
   fChain->SetBranchAddress("slc_vtxcheck_angle", &slc_vtxcheck_angle, &b_slc_vtxcheck_angle);
   fChain->SetBranchAddress("slc_origin", &slc_origin, &b_slc_origin);
   fChain->SetBranchAddress("slc_nhits_u", &slc_nhits_u, &b_slc_nhits_u);
   fChain->SetBranchAddress("slc_nhits_v", &slc_nhits_v, &b_slc_nhits_v);
   fChain->SetBranchAddress("slc_nhits_w", &slc_nhits_w, &b_slc_nhits_w);
   fChain->SetBranchAddress("slc_longesttrack_length", &slc_longesttrack_length, &b_slc_longesttrack_length);
   fChain->SetBranchAddress("slc_longesttrack_phi", &slc_longesttrack_phi, &b_slc_longesttrack_phi);
   fChain->SetBranchAddress("slc_longesttrack_theta", &slc_longesttrack_theta, &b_slc_longesttrack_theta);
   fChain->SetBranchAddress("slc_longesttrack_iscontained", &slc_longesttrack_iscontained, &b_slc_longesttrack_iscontained);
   fChain->SetBranchAddress("slc_acpt_outoftime", &slc_acpt_outoftime, &b_slc_acpt_outoftime);
   fChain->SetBranchAddress("slc_crosses_top_boundary", &slc_crosses_top_boundary, &b_slc_crosses_top_boundary);
   fChain->SetBranchAddress("slc_nuvtx_closetodeadregion_u", &slc_nuvtx_closetodeadregion_u, &b_slc_nuvtx_closetodeadregion_u);
   fChain->SetBranchAddress("slc_nuvtx_closetodeadregion_v", &slc_nuvtx_closetodeadregion_v, &b_slc_nuvtx_closetodeadregion_v);
   fChain->SetBranchAddress("slc_nuvtx_closetodeadregion_w", &slc_nuvtx_closetodeadregion_w, &b_slc_nuvtx_closetodeadregion_w);
   fChain->SetBranchAddress("slc_kalman_chi2", &slc_kalman_chi2, &b_slc_kalman_chi2);
   fChain->SetBranchAddress("slc_kalman_ndof", &slc_kalman_ndof, &b_slc_kalman_ndof);
   fChain->SetBranchAddress("slc_passed_min_track_quality", &slc_passed_min_track_quality, &b_slc_passed_min_track_quality);
   fChain->SetBranchAddress("slc_n_intime_pe_closestpmt", &slc_n_intime_pe_closestpmt, &b_slc_n_intime_pe_closestpmt);
   fChain->SetBranchAddress("slc_maxdistance_vtxtrack", &slc_maxdistance_vtxtrack, &b_slc_maxdistance_vtxtrack);
   fChain->SetBranchAddress("slc_npfp", &slc_npfp, &b_slc_npfp);
   fChain->SetBranchAddress("slc_ntrack", &slc_ntrack, &b_slc_ntrack);
   fChain->SetBranchAddress("slc_nshower", &slc_nshower, &b_slc_nshower);
   fChain->SetBranchAddress("slc_iscontained", &slc_iscontained, &b_slc_iscontained);
   fChain->SetBranchAddress("slc_mult_pfp", &slc_mult_pfp, &b_slc_mult_pfp);
   fChain->SetBranchAddress("slc_mult_track", &slc_mult_track, &b_slc_mult_track);
   fChain->SetBranchAddress("slc_mult_shower", &slc_mult_shower, &b_slc_mult_shower);
   fChain->SetBranchAddress("slc_mult_track_tolerance", &slc_mult_track_tolerance, &b_slc_mult_track_tolerance);
   fChain->SetBranchAddress("nbeamfls", &nbeamfls, &b_nbeamfls);
   fChain->SetBranchAddress("beamfls_time", &beamfls_time, &b_beamfls_time);
   fChain->SetBranchAddress("beamfls_pe", &beamfls_pe, &b_beamfls_pe);
   fChain->SetBranchAddress("beamfls_z", &beamfls_z, &b_beamfls_z);
   fChain->SetBranchAddress("no_mcflash_but_op_activity", &no_mcflash_but_op_activity, &b_no_mcflash_but_op_activity);
   fChain->SetBranchAddress("beamfls_spec", &beamfls_spec, &b_beamfls_spec);
   fChain->SetBranchAddress("numc_flash_spec", &numc_flash_spec, &b_numc_flash_spec);
   fChain->SetBranchAddress("slc_flshypo_xfixed_spec", &slc_flshypo_xfixed_spec, &b_slc_flshypo_xfixed_spec);
   fChain->SetBranchAddress("slc_flshypo_spec", &slc_flshypo_spec, &b_slc_flshypo_spec);
   fChain->SetBranchAddress("nsignal", &nsignal, &b_nsignal);
   fChain->SetBranchAddress("mctrk_start_x", &mctrk_start_x, &b_mctrk_start_x);
   fChain->SetBranchAddress("mctrk_start_y", &mctrk_start_y, &b_mctrk_start_y);
   fChain->SetBranchAddress("mctrk_start_z", &mctrk_start_z, &b_mctrk_start_z);
   fChain->SetBranchAddress("trk_start_x", &trk_start_x, &b_trk_start_x);
   fChain->SetBranchAddress("trk_start_y", &trk_start_y, &b_trk_start_y);
   fChain->SetBranchAddress("trk_start_z", &trk_start_z, &b_trk_start_z);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("tvtx_x", &tvtx_x, &b_tvtx_x);
   fChain->SetBranchAddress("tvtx_y", &tvtx_y, &b_tvtx_y);
   fChain->SetBranchAddress("tvtx_z", &tvtx_z, &b_tvtx_z);
   fChain->SetBranchAddress("pot", &pot, &b_pot);
   Notify();
}

Bool_t AnaTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnaTree_cxx

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 13 10:53:40 2018 by ROOT version 6.06/06
// from TTree shower_tree/shower_tree
// found on file: ubxsecana_output_merged_mcc8.6.root
//////////////////////////////////////////////////////////

#ifndef ShowerTree_h
#define ShowerTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ShowerTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        s_nupdg;
   Double_t        s_track_pdg;
   Double_t        s_tpcobj_origin;
   Double_t        s_shower_length;
   Double_t        s_shower_phi;
   Double_t        s_shower_theta;
   Double_t        s_shower_openangle;
   Double_t        s_shower_startx;
   Double_t        s_shower_starty;
   Double_t        s_shower_startz;
   Double_t        s_flash_z;

   // List of branches
   TBranch        *b_s_nupdg;   //!
   TBranch        *b_s_track_pdg;   //!
   TBranch        *b_s_tpcobj_origin;   //!
   TBranch        *b_s_shower_length;   //!
   TBranch        *b_s_shower_phi;   //!
   TBranch        *b_s_shower_theta;   //!
   TBranch        *b_s_shower_openangle;   //!
   TBranch        *b_s_shower_startx;   //!
   TBranch        *b_s_shower_starty;   //!
   TBranch        *b_s_shower_startz;   //!
   TBranch        *b_s_flash_z;   //!

   ShowerTree(TTree *tree=0);
   virtual ~ShowerTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ShowerTree_cxx
ShowerTree::ShowerTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ubxsecana_output_merged_mcc8.6_dist.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ubxsecana_output_merged_mcc8.6_dist.root");
      }
      f->GetObject("shower_tree",tree);

   }
   Init(tree);
}

ShowerTree::~ShowerTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ShowerTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ShowerTree::LoadTree(Long64_t entry)
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

void ShowerTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("s_nupdg", &s_nupdg, &b_s_nupdg);
   fChain->SetBranchAddress("s_track_pdg", &s_track_pdg, &b_s_track_pdg);
   fChain->SetBranchAddress("s_tpcobj_origin", &s_tpcobj_origin, &b_s_tpcobj_origin);
   fChain->SetBranchAddress("s_shower_length", &s_shower_length, &b_s_shower_length);
   fChain->SetBranchAddress("s_shower_phi", &s_shower_phi, &b_s_shower_phi);
   fChain->SetBranchAddress("s_shower_theta", &s_shower_theta, &b_s_shower_theta);
   fChain->SetBranchAddress("s_shower_openangle", &s_shower_openangle, &b_s_shower_openangle);
   fChain->SetBranchAddress("s_shower_startx", &s_shower_startx, &b_s_shower_startx);
   fChain->SetBranchAddress("s_shower_starty", &s_shower_starty, &b_s_shower_starty);
   fChain->SetBranchAddress("s_shower_startz", &s_shower_startz, &b_s_shower_startz);
   fChain->SetBranchAddress("s_flash_z", &s_flash_z, &b_s_flash_z);
   Notify();
}

Bool_t ShowerTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ShowerTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ShowerTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ShowerTree_cxx

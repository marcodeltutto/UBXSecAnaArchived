//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 30 11:34:20 2017 by ROOT version 6.06/06
// from TTree mom_tree_contained/
// found on file: ../files/mcc8.3_v2.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v2.0.2.root
//////////////////////////////////////////////////////////

#ifndef MomTreeContainedAna_h
#define MomTreeContainedAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MomTreeContainedAna {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  Int_t           run;
  Int_t           subrun;
  Int_t           event;
  Double_t        mom_true_contained;
  Double_t        mom_mcs_contained;
  Double_t        mom_range_contained;
  
  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_subrun;   //!
  TBranch        *b_event;   //!
  TBranch        *b_mom_true_contained;   //!
  TBranch        *b_mom_mcs_contained;   //!
  TBranch        *b_mom_range_contained;   //!
  
  MomTreeContainedAna(TTree *tree=0);
  virtual ~MomTreeContainedAna();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  
  void DoFit(TH1D*, std::string, double *);
};

#endif

#ifdef MomTreeContainedAna_cxx
MomTreeContainedAna::MomTreeContainedAna(TTree *tree) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../files/mcc8.3_v2.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v2.0.2.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("../../files/mcc8.3_v2.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v2.0.2.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("../../files/mcc8.3_v2.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v2.0.2.root:/UBXSec");
    dir->GetObject("mom_tree_contained",tree);
    
  }
  Init(tree);
}

MomTreeContainedAna::~MomTreeContainedAna()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t MomTreeContainedAna::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t MomTreeContainedAna::LoadTree(Long64_t entry)
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

void MomTreeContainedAna::Init(TTree *tree)
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
  
  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("mom_true_contained", &mom_true_contained, &b_mom_true_contained);
  fChain->SetBranchAddress("mom_mcs_contained", &mom_mcs_contained, &b_mom_mcs_contained);
  fChain->SetBranchAddress("mom_range_contained", &mom_range_contained, &b_mom_range_contained);
  Notify();
}

Bool_t MomTreeContainedAna::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void MomTreeContainedAna::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t MomTreeContainedAna::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef MomTreeContainedAna_cxx

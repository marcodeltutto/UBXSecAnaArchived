//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 30 15:03:21 2017 by ROOT version 6.06/06
// from TTree mom_tree_uncontained/
// found on file: ../../files/mcc8.3_v2.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v2.0.2.root
//////////////////////////////////////////////////////////

#ifndef MomTreeUnContainedAna_h
#define MomTreeUnContainedAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <fstream>


// Header file for the classes stored in the TTree if any.

class MomTreeUnContainedAna {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  Int_t           run;
  Int_t           subrun;
  Int_t           event;
  Double_t        mom_true_uncontained;
  Double_t        mom_mcs_uncontained;
  
  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_subrun;   //!
  TBranch        *b_event;   //!
  TBranch        *b_mom_true_uncontained;   //!
  TBranch        *b_mom_mcs_uncontained;   //!
  
  MomTreeUnContainedAna(TTree *tree=0);
  virtual ~MomTreeUnContainedAna();
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

#ifdef MomTreeUnContainedAna_cxx
MomTreeUnContainedAna::MomTreeUnContainedAna(TTree *tree) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../files/mcc8.3_v2.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v2.0.2.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("../../files/mcc8.3_v2.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v2.0.2.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("../../files/mcc8.3_v2.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v2.0.2.root:/UBXSec");
    dir->GetObject("mom_tree_uncontained",tree);
    
  }
  Init(tree);
}

MomTreeUnContainedAna::~MomTreeUnContainedAna()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t MomTreeUnContainedAna::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t MomTreeUnContainedAna::LoadTree(Long64_t entry)
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

void MomTreeUnContainedAna::Init(TTree *tree)
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
  fChain->SetBranchAddress("mom_true_uncontained", &mom_true_uncontained, &b_mom_true_uncontained);
  fChain->SetBranchAddress("mom_mcs_uncontained", &mom_mcs_uncontained, &b_mom_mcs_uncontained);
  Notify();
}

Bool_t MomTreeUnContainedAna::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void MomTreeUnContainedAna::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t MomTreeUnContainedAna::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef MomTreeUnContainedAna_cxx

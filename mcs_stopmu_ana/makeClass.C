void makeClass() {

  TFile *f = new TFile("../files/mcc8.3_stopmu_test/ubxsec_output_mc_bnbcosmic_mcc8.3_stopmu_test10.root");
  TTree *v = (TTree*)f->Get("pandoraCosmicStoppingMu/tree");
  v->MakeClass("AnaTree");
}

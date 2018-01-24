void makeClass() {

  TFile *f = new TFile("../files/mcc8.3_stopmu_test/ubxsec_output_mc_bnbcosmic_mcc8.3_stopmu_test11.root");
  TTree *v = (TTree*)f->Get("UBXSec/tree");
  v->MakeClass("AnaTree");
}

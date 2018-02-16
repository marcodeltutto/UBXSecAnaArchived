void makeClass() {

  TFile *f = new TFile("../files/mcc8.6_stopmu_test/ubxsec_output_mc_bnbnuecosmic_mcc8.6_stopmu_test17.root");
  TTree *v = (TTree*)f->Get("UBXSec/tree");
  v->MakeClass("AnaTree");
}

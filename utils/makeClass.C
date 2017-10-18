void makeClass() {

  TFile *f = new TFile("../files/mcc8.3_v1.0.2/ubxsec_output_mc_bnbcosmic_mcc8.3_v1.0.2.root");
  TTree *v = (TTree*)f->Get("UBXSec/tree");
  v->MakeClass("AnaTree");
}

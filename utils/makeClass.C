void makeClass() {

  TFile *f = new TFile("../files/ubxsec_output_mc_7.root");
  TTree *v = (TTree*)f->Get("UBXSec/tree");
  v->MakeClass("AnaTree");
}

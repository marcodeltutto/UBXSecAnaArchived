void makeClass() {

  TFile *f = new TFile("../files/output_apr24.root");
  TTree *v = (TTree*)f->Get("UBXSec/tree");
  v->MakeClass("AnaTree");
}

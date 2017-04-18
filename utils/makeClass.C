void makeClass() {

  TFile *f = new TFile("../files/output9850.root");
  TTree *v = (TTree*)f->Get("UBXSec/tree");
  v->MakeClass("AnaTree");
}

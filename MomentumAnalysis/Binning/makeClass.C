void makeClass() {

  TFile *f = new TFile("./ubxsecana_output.root");
  TTree *v = (TTree*)f->Get("mom_tree");
  v->MakeClass("AnaTree");
}

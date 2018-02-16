void makeClass() {

  TFile *f = new TFile("ubxsecana_output_merged_mcc8.6.root");
  TTree *v = (TTree*)f->Get("shower_tree");
  v->MakeClass("ShowerTree");
}

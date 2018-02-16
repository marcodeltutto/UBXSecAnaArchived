void makeClass() {

  TFile *f = new TFile("./ubxsecana_output_bnbcosmic_mcc8.6.root ");
  TTree *v = (TTree*)f->Get("mom_tree");
  v->MakeClass("BinFinder");
}

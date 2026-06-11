void kaon_decay_data(){
  const char* file = "t110_graph_344_k.root";
  TFile* f = TFile::Open(file, "READ");
  
  TH2D* hist_bac_btof = (TH2D*)f->Get("hist_bac_btof");
  TH2D* hist_bac_btof_pass = (TH2D*)f->Get("hist_bac_btof_pass");

  double total = hist_bac_btof->Integral();
  double pass = hist_bac_btof_pass->Integral();

  cout<<pass / total*100<<endl;

  
}

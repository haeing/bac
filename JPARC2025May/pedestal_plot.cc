void pedestal_plot(){
  TFile *f = TFile::Open("t110_graph_354.root","READ");
  TH1D *hist = (TH1D*)f->Get("hist_bac_npe_s");

  hist->GetXaxis()->SetTitle("N_{p.e.}");
  hist->GetYaxis()->SetTitle("Counts");
  hist->SetLineColor(kBlack);
  
  auto c1 = new TCanvas("c1","c1");
  hist->Draw();
  

}

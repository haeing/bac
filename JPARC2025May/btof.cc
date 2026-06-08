void btof(){
  gROOT->SetBatch(false);
  
  int runnumber = 344;
  TString dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025May_root";
  TFile *file = new TFile(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber));
  TTree *data = (TTree*)file->Get("hodo");

  double btof0;
  data->SetBranchAddress("btof0",&btof0);

  TBox *box_pi = new TBox(-0.9,0,1.3,90000);
  box_pi->SetFillColorAlpha(kBlue,0.12);
  box_pi->SetLineColor(0);

  TBox *box_k = new TBox(4.5,0,5.6,90000);
  box_k->SetFillColorAlpha(kRed,0.12);
  box_k->SetLineColor(0);

  
  TH1D *h_btof = new TH1D("h_btof","h_btof",100,-2,7);
  for(int n=0;n<data->GetEntries();n++){
    data->GetEntry(n);
    h_btof->Fill(btof0*-1);
  }
  h_btof->SetTitle(";Time of Flight [ns];Counts");
  auto c1 = new TCanvas("c1","c1");
  h_btof->SetStats(0);
  h_btof->GetYaxis()->SetRangeUser(0,90000);
  h_btof->Draw();
  h_btof->SetLineColor(kBlack);
  box_pi->Draw("same");
  box_k->Draw("same");
  h_btof->Draw("same");
  //c1->SetLogy();
  c1->SaveAs("btof.pdf");
}

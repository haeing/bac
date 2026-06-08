void btof(){
  gROOT->SetBatch(false);
  
  int runnumber = 344;
  TString dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025May_root";
  TFile *file = new TFile(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber));
  TTree *data = (TTree*)file->Get("hodo");

  double btof0;
  data->SetBranchAddress("btof0",&btof0);

  TH1D *h_btof = new TH1D("h_btof","h_btof",100,-2,7);
  for(int n=0;n<data->GetEntries();n++){
    data->GetEntry(n);
    h_btof->Fill(btof0*-1);
  }
  h_btof->SetTitle(";Time of Flight [ns];Counts");
  auto c1 = new TCanvas("c1","c1");
  h_btof->SetStats(0);
  h_btof->Draw();
  h_btof->SetLineColor(kBlack);
  //c1->SetLogy();
  c1->SaveAs("btof.pdf");
}

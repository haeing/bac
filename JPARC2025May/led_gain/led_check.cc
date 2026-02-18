void led_check(int runnumber){

  gROOT->SetBatch(kTRUE);
  
  TString dir = "/gpfs/group/had/sks/Users/haein/JPARC2025May_root";
  TFile *file = new TFile(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber));
  TTree *data = (TTree*)file->Get("hodo");

  vector<double> *bac_adc_u = nullptr;
  TBranch *b_bac_adc_u = data->GetBranch("bac_adc_u");
  b_bac_adc_u->SetAddress(&bac_adc_u);

  TH1D* hist[4];
  for(int i=0;i<4;i++){
    hist[i] = new TH1D(Form("hist%d",i),Form("hist%d",i),500,0,500);
  }
  
  for(int n=0;n<data->GetEntries();n++){
    data->GetEntry(n);
    for(int i=0;i<4;i++){
      hist[i]->Fill((*bac_adc_u)[i]);
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(2,2);
  for(int i=0;i<4;i++){
    c1->cd(i+1);
    hist[i]->Draw();
      
  }

  c1->SaveAs(Form("fig/run00%d_bac_led.pdf",runnumber));
}

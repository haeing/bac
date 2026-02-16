void npe_picture(){

  gStyle->SetOptStat(0);

  TFile *file1 = new TFile("build/two_layer_E72.root");
  TFile *file2 = new TFile("build/three_layer_E72.root");

  TTree *tree1 = (TTree*)file1->Get("tree");
  TTree *tree2 = (TTree*)file2->Get("tree");

  int nhit1;
  int nhit2;
  double evtposx1;
  double evtposy1;
  double evtposx2;
  double evtposy2;
  tree1->SetBranchAddress("nhMppc",&nhit1);
  tree1->SetBranchAddress("evtposx",&evtposx1);
  tree1->SetBranchAddress("evtposy",&evtposy1);
  tree2->SetBranchAddress("nhMppc",&nhit2);
  tree2->SetBranchAddress("evtposx",&evtposx2);
  tree2->SetBranchAddress("evtposy",&evtposy2);

  TH1D *hist1 = new TH1D("hist1","hist1",100,0,100);
  TH1D *hist2 = new TH1D("hist2","hist2",100,0,100);

  double to1=0;
  double count1=0;

  double to2=0;
  double count2=0;
  for(int n=0;n<10000;n++){
    tree1->GetEntry(n);
    if(evtposx1>-57.5&&evtposx1<57.5&&evtposy1>-57.5&&evtposy1<57.5){
      hist1->Fill(nhit1);
      to1++;
      if(nhit1>15)count1++;
     }
    tree2->GetEntry(n);
    if(evtposx2>-57.5&&evtposx2<57.5&&evtposy2>-57.5&&evtposy2<57.5){
      hist2->Fill(nhit2);
      to2++;
      if(nhit2>15)count2++;
      }
  }


  cout<<"two : "<<count1/to1*100<<endl;
  cout<<"three : "<<count2/to2*100<<endl;
  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);
  le->AddEntry(hist1,"Aerogel d = 24.6 mm");
  le->AddEntry(hist2,"Aerogel d = 37.0 mm");

  hist2->SetLineColor(kBlue);
  hist1->SetLineWidth(2);
  hist2->SetLineWidth(2);
  hist1->SetTitle(";N_{p.e.};n");
  TCanvas *c1 = new TCanvas ("c1","c1");
  hist1->Draw();
  hist2->Draw("same");
  le->Draw();

  
}

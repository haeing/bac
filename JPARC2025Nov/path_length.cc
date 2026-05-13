void path_length(){
  TString dir = "/gpfs/group/had/sks/Users/haein/simul-data/E72_Simul";
  TFile *file = new TFile(Form("%s/beam_735_beam.root",dir.Data()));
  TTree *tree = (TTree*)file->Get("g4hyptpc");

  vector<TParticle> *BAC = nullptr;
  TBranch *b_BAC = tree->GetBranch("BAC");
  b_BAC->SetAddress(&BAC);

  TH2D *hist_pos = new TH2D("hist_pos","hist_pos",100,-520,-475,100,-80,50);
  TH1D *hist_path = new TH1D("hist_path","hist_path",50,0,50);
  TH1D *hist_dist = new TH1D("hist_dist","hist_dist",50,0,50);
  TCanvas *c = new TCanvas("c","c");
  
  auto mg = new TMultiGraph();
  for(int n=0;n<tree->GetEntries();n++){
  //for(int n=0;n<100;n++){
    tree->GetEntry(n);
    if(BAC->size()<2)continue;
    TGraph *g_track = new TGraph();
    TParticle &p_first = BAC->front();
    TParticle &p_last  = BAC->back();
    double dx = p_last.Vx() - p_first.Vx();
    double dy = p_last.Vy() - p_first.Vy();
    double dz = p_last.Vz() - p_first.Vz();
    
    double distance = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
    hist_dist->Fill(distance);
    double track_length = 0.;
    for(size_t j=0;j<BAC->size();j++){
      TParticle &p_bac = (*BAC)[j];
      g_track->SetPoint(g_track->GetN(),p_bac.Vz(),p_bac.Vx());
      g_track->SetMarkerStyle(20);
      mg->Add(g_track,"APL");
      hist_pos->Fill(p_bac.Vz(),p_bac.Vx());
      if(j>0){
	TParticle &p1 = (*BAC)[j-1];
	TParticle &p2 = (*BAC)[j];
      
	double dx = p2.Vx() - p1.Vx();
	double dy = p2.Vy() - p1.Vy();
	double dz = p2.Vz() - p1.Vz();
	double dl = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
	track_length += dl;
      }
    }
    hist_path->Fill(track_length);
  }
  c->cd();
  mg->Draw("AL");
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  c1->cd(1);
  hist_pos->Draw("colz");
  c1->cd(2);
  hist_path->Draw();
  c1->cd(3);
  hist_dist->Draw();

}

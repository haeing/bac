void draw_graph(){
  TFile* f = TFile::Open("t110_graph_321.root");
  TGraphErrors* g_bac_npe_mean = (TGraphErrors*)f->Get("g_bac_npe_mean");
  TGraphAsymmErrors *g_eff = (TGraphAsymmErrors*)f->Get("g_eff");

  TFile *f_beam = new TFile("t110_beam_735.root");
  TTree *t_beam = (TTree*)f_beam->Get("tree");
  int seg_bh2;
  t_beam->SetBranchAddress("seg_bh2",&seg_bh2);
  
  TFile *f_simul = TFile::Open("../../../data/JPARC2025May/g4_root/T110_simul.root");
  TTree* t_simul = (TTree*)f_simul->Get("tree");
  int nhMppc;
  double evtposx,evtposy;
  t_simul->SetBranchAddress("nhMppc",&nhMppc);
  t_simul->SetBranchAddress("evtposx",&evtposx);
  t_simul->SetBranchAddress("evtposy",&evtposy);
  

  TH1D *h_npe[12];
  TH2D *h_xy[12];
  TF1 *f_npe[12];
  for(int i=0;i<12;i++){
    h_npe[i] = new TH1D(Form("h_npe%d",i),Form("h_npe%d",i),100,0,100);
    h_xy[i] = new TH2D(Form("h_xy%d",i),Form("h_xy%d",i),200,-200,200,200,-200,200);
    f_npe[i] = new TF1(Form("f_npe%d",i),"gaus",0,100);
  }
			    
  for(int n=0;n<t_beam->GetEntries();n++){
    t_beam->GetEntry(n);
    t_simul->GetEntry(n);
    h_npe[seg_bh2]->Fill(nhMppc);
    h_xy[seg_bh2]->Fill(evtposx,evtposy);
  }

  TCanvas *c0 = new TCanvas("c0","c0");
  c0->Divide(4,3);
  TCanvas *c10 = new TCanvas("c10","c10");
  c10->Divide(4,3);
  TGraphErrors *g_simul_npe = new TGraphErrors();
  double* x_npe  = g_bac_npe_mean->GetX();
  double* ex_npe = g_bac_npe_mean->GetEX();
  double m_npe[7];
  double em_npe[7];
  int pid = 0;
  for(int i=0;i<12;i++){
    c0->cd(i+1);
    h_npe[i]->Fit(f_npe[i],"","",30,50);
    if(i>2 && i < 10){
      m_npe[pid] = f_npe[i]->GetParameter(1);
      em_npe[pid] = f_npe[i]->GetParError(1);
      cout<<x_npe[pid]<<endl;
      g_simul_npe->SetPoint(pid,x_npe[pid],m_npe[pid]);
      g_simul_npe->SetPointError(pid,ex_npe[pid],em_npe[pid]);
      pid++;
    }
    c10->cd(i+1);
    h_xy[i]->Draw("colz");
  }
  
  
  TCanvas *c1 = new TCanvas("c1","c1");
  
  g_simul_npe->SetMarkerStyle(21);
  //g_bac_npe_mean->SetMarkerColor(kRed);
  //g_simul_npe->SetMarkerSize(2);
  int n = g_bac_npe_mean->GetN();
  g_bac_npe_mean->Set(n-1);
  TLegend* leg = new TLegend(0.45, 0.55, 0.88, 0.82);
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);
  
  leg->AddEntry(g_bac_npe_mean,  "Experiment", "p");
  leg->AddEntry(g_simul_npe,  "Simulation", "p");
 
  auto mg = new TMultiGraph();
  mg->Add(g_bac_npe_mean);
  //mg->Add(g_simul_npe);
  mg->GetYaxis()->SetRangeUser(10,50);
  mg->GetXaxis()->SetTitle("X [mm]");
  mg->GetYaxis()->SetTitle("Np.e.");
			   
  mg->Draw("AP");
  //leg->Draw("same");
  
  TCanvas *c2 = new TCanvas("c2","c2");
  g_eff->GetXaxis()->SetTitle("X [mm]");
  g_eff->GetYaxis()->SetTitle("Pion Efficiency");
  g_eff->Set(n-1);
  g_eff->Draw("AP");
  
  
}

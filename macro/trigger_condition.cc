static const double mK = 493.677;
static const double mp = 938.27208816;
static const double mn = 939.5654205;
static const double mpi = 139.57039;
static const double mpi0 = 134.9768;
static const double me = 0.51099895;
static const double mmu = 105.66;

static const double meta =547.862;
static const double mLambda = 1115.683;
static const double BAC_thre = 1/1.115;
static const double KVC_thre = 1/1.45;

double GetBeta(int pdg, double mom){
  double mass = -9999;
  if(TMath::Abs(pdg) == 321)mass = mK;
  else if(TMath::Abs(pdg) == 2212)mass = mp;
  else if(TMath::Abs(pdg) == 2112)mass = mn;
  else if(TMath::Abs(pdg) == 211)mass = mpi;
  else if(TMath::Abs(pdg) == 111)mass = mpi0;
  else if(TMath::Abs(pdg) == 11)mass = me;
  else if(TMath::Abs(pdg) == 13)mass = mmu;

  if(mass == -9999)return -9999;
  else{ return mom/TMath::Sqrt(mass*mass + mom*mom);}
  
}


void trigger_condition(){

  string file_directory = "/gpfs/group/had/sks/Users/haein/E72_BAC_KVC_beta/";
  
  gStyle->SetOptStat(0);


  double bac_min = 0.7;
  double bac_max = 1.05;

  double kvc_min = 0.1;
  double kvc_max = 1.05;
  double bin = 0.005;

  //Histogram
  TH1D *bac_K = new TH1D("bac_K","bac_K",int((bac_max - bac_min)/bin),bac_min,bac_max);
  TH1D *bac_pi = new TH1D("bac_pi","bac_pi",int((bac_max - bac_min)/bin),bac_min,bac_max);
  TH1D *bac_mu = new TH1D("bac_mu","bac_mu",int((bac_max - bac_min)/bin),bac_min,bac_max);
  TH1D *bac_R_p = new TH1D("bac_R_p","bac_R_p",int((bac_max - bac_min)/bin),bac_min,bac_max);
  TH1D *bac_R_pim = new TH1D("bac_R_pim","bac_R_pim",int((bac_max - bac_min)/bin),bac_min,bac_max);
  TH1D *bac_R_pip = new TH1D("bac_R_pip","bac_R_pip",int((bac_max - bac_min)/bin),bac_min,bac_max);
 

  TH1D *kvc_K = new TH1D("kvc_K","kvc_K",int((kvc_max - kvc_min)/bin),kvc_min,kvc_max);
  TH1D *kvc_pi = new TH1D("kvc_pi","kvc_pi",int((kvc_max - kvc_min)/bin),kvc_min,kvc_max);
  TH1D *kvc_mu = new TH1D("kvc_mu","kvc_mu",int((kvc_max - kvc_min)/bin),kvc_min,kvc_max);
  TH1D *kvc_R_p = new TH1D("kvc_R_p","kvc_R_p",int((kvc_max - kvc_min)/bin),kvc_min,kvc_max);
  TH1D *kvc_R_pim = new TH1D("kvc_R_pim","kvc_R_pim",int((kvc_max - kvc_min)/bin),kvc_min,kvc_max);
  TH1D *kvc_R_pip = new TH1D("kvc_R_pip","kvc_R_pip",int((kvc_max - kvc_min)/bin),kvc_min,kvc_max);

  //kaon beam
  TFile *file_K = new TFile((file_directory+"kaon_beam.root").c_str());
  TTree *tree_K = (TTree*)file_K->Get("g4hyptpc");

  TFile *file_r = new TFile((file_directory+"Reaction.root").c_str());
  TTree *tree_r = (TTree*)file_r->Get("g4hyptpc");

  vector<TParticle>* BAC_K = nullptr;
  vector<TParticle>* KVC_K = nullptr;
  vector<TParticle>* TGT_K = nullptr;
  

  TBranch *b_BAC_K = tree_K->GetBranch("BAC");
  b_BAC_K->SetAddress(&BAC_K);

  TBranch *b_KVC_K = tree_K->GetBranch("KVC");
  b_KVC_K->SetAddress(&KVC_K);

  TBranch *b_TGT_K = tree_K->GetBranch("TGT");
  b_TGT_K->SetAddress(&TGT_K);

  int count_tgt = 0;
  for(int n=0;n<tree_K->GetEntries();n++){
    tree_K->GetEntry(n);

    if(TGT_K){
      const auto &p_TGT_K = (*TGT_K)[0];
      double x = p_TGT_K.Vx();
      double z = p_TGT_K.Vz()+143.0;
      if(TMath::Sqrt(p_TGT_K.Vx()*p_TGT_K.Vx()+z*z)<=40 && TMath::Abs(p_TGT_K.Vy())<=50)
	count_tgt++;
    }
    
    if(BAC_K){
      //for(const auto &p_BAC_K : *BAC_K){
      const auto &p_BAC_K = (*BAC_K)[0];
      int BAC_pdg_K = p_BAC_K.GetPdgCode();
      double BAC_mom_K = TMath::Sqrt(p_BAC_K.Px()*p_BAC_K.Px()+p_BAC_K.Py()*p_BAC_K.Py()+p_BAC_K.Pz()*p_BAC_K.Pz());

      if(BAC_pdg_K == -321){
	bac_K->Fill(GetBeta(BAC_pdg_K,BAC_mom_K));
      }
      else if(BAC_pdg_K == -211){
	bac_pi->Fill(GetBeta(BAC_pdg_K,BAC_mom_K));
      }
      else if(BAC_pdg_K == 13){
	bac_mu->Fill(GetBeta(BAC_pdg_K,BAC_mom_K));
      }
    }


    //no-reaction
    if(n>tree_r->GetEntries()){
      if(KVC_K){
	const auto &p_KVC_K = (*KVC_K)[0];
	int KVC_pdg_K = p_KVC_K.GetPdgCode();
	double KVC_mom_K = TMath::Sqrt(p_KVC_K.Px()*p_KVC_K.Px()+p_KVC_K.Py()*p_KVC_K.Py()+p_KVC_K.Pz()*p_KVC_K.Pz());

	if(KVC_pdg_K == -321){
	  kvc_K->Fill(GetBeta(KVC_pdg_K,KVC_mom_K));
	}
	else if(KVC_pdg_K == -211){
	  kvc_pi->Fill(GetBeta(KVC_pdg_K,KVC_mom_K));
	}
	else if(KVC_pdg_K == 13){
	  kvc_mu->Fill(GetBeta(KVC_pdg_K,KVC_mom_K));
	}
      }
    }
  }

  cout<<"Rbeam = "<<(double)count_tgt/(double)tree_K->GetEntries()<<endl;
  //pion beam
  TFile *file_pi = new TFile((file_directory+"pion_beam.root").c_str());
  TTree *tree_pi = (TTree*)file_pi->Get("g4hyptpc");

  vector<TParticle>* BAC_pi = nullptr;
  vector<TParticle>* KVC_pi = nullptr;

  TBranch *b_BAC_pi = tree_pi->GetBranch("BAC");
  b_BAC_pi->SetAddress(&BAC_pi);

  TBranch *b_piVC_pi = tree_pi->GetBranch("KVC");
  b_piVC_pi->SetAddress(&KVC_pi);
  
  for(int n=0;n<tree_pi->GetEntries();n++){
    tree_pi->GetEntry(n);

    if(BAC_pi){
      const auto &p_BAC_pi = (*BAC_pi)[0];
      int BAC_pdg_pi = p_BAC_pi.GetPdgCode();
      double BAC_mom_pi = TMath::Sqrt(p_BAC_pi.Px()*p_BAC_pi.Px()+p_BAC_pi.Py()*p_BAC_pi.Py()+p_BAC_pi.Pz()*p_BAC_pi.Pz());

      if(BAC_pdg_pi == -321){
	bac_K->Fill(GetBeta(BAC_pdg_pi,BAC_mom_pi));
      }
      else if(BAC_pdg_pi == -211){
	bac_pi->Fill(GetBeta(BAC_pdg_pi,BAC_mom_pi));
      }
      else if(BAC_pdg_pi == 13){
	bac_mu->Fill(GetBeta(BAC_pdg_pi,BAC_mom_pi));
      }
    }

    if(KVC_pi){
      const auto &p_KVC_pi = (*KVC_pi)[0];
      int KVC_pdg_pi = p_KVC_pi.GetPdgCode();
      double KVC_mom_pi = TMath::Sqrt(p_KVC_pi.Px()*p_KVC_pi.Px()+p_KVC_pi.Py()*p_KVC_pi.Py()+p_KVC_pi.Pz()*p_KVC_pi.Pz());

      if(KVC_pdg_pi == -321){
	kvc_K->Fill(GetBeta(KVC_pdg_pi,KVC_mom_pi));
      }
      else if(KVC_pdg_pi == -211){
	kvc_pi->Fill(GetBeta(KVC_pdg_pi,KVC_mom_pi));
      }
      else if(KVC_pdg_pi == 13){
	kvc_mu->Fill(GetBeta(KVC_pdg_pi,KVC_mom_pi));
      }
    }
  }

  
  //etalambda reaction
  

  vector<TParticle>* BAC_r = nullptr;
  vector<TParticle>* KVC_r = nullptr;

  TBranch *b_BAC_r = tree_r->GetBranch("BAC");
  b_BAC_r->SetAddress(&BAC_r);

  TBranch *b_rVC_r = tree_r->GetBranch("KVC");
  b_rVC_r->SetAddress(&KVC_r);
  
  for(int n=0;n<tree_r->GetEntries();n++){
    tree_r->GetEntry(n);

    if(BAC_r && !BAC_r->empty()){
      const auto &p_BAC_r = (*BAC_r)[0];
      int BAC_pdg_r = p_BAC_r.GetPdgCode();
      double BAC_mom_r = TMath::Sqrt(p_BAC_r.Px()*p_BAC_r.Px()+p_BAC_r.Py()*p_BAC_r.Py()+p_BAC_r.Pz()*p_BAC_r.Pz());
      
      if(BAC_pdg_r == -211){
	bac_R_pim->Fill(GetBeta(BAC_pdg_r,BAC_mom_r));
      }
      else if(BAC_pdg_r == 211){
	bac_R_pip->Fill(GetBeta(BAC_pdg_r,BAC_mom_r));
      }
      else if(BAC_pdg_r == 2212){
	bac_R_p->Fill(GetBeta(BAC_pdg_r,BAC_mom_r));
      }

    }

    if(KVC_r && !KVC_r->empty()){
      const auto &p_KVC_r = (*KVC_r)[0];
      int KVC_pdg_r = p_KVC_r.GetPdgCode();
      double KVC_mom_r = TMath::Sqrt(p_KVC_r.Px()*p_KVC_r.Px()+p_KVC_r.Py()*p_KVC_r.Py()+p_KVC_r.Pz()*p_KVC_r.Pz());

      if(KVC_pdg_r == -211){
	kvc_R_pim->Fill(GetBeta(KVC_pdg_r,KVC_mom_r));
      }
      else if(KVC_pdg_r == 211){
	kvc_R_pip->Fill(GetBeta(KVC_pdg_r,KVC_mom_r));
      }
      else if(KVC_pdg_r == 2212){
	kvc_R_p->Fill(GetBeta(KVC_pdg_r,KVC_mom_r));
      }
      
    }
  }


  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd();
  gPad->SetLogy();
  bac_K->SetLineColor(kBlack);
  bac_pi->SetLineColor(kRed);
  bac_mu->SetLineColor(kCyan);
  bac_R_p->SetLineColor(kMagenta);
  bac_R_pip->SetLineColor(kGreen);
  bac_R_pim->SetLineColor(kBlue);
  bac_K->SetLineWidth(2);
  bac_pi->SetLineWidth(2);
  bac_mu->SetLineWidth(2);
  bac_R_p->SetLineWidth(2);
  bac_R_pip->SetLineWidth(2);
  bac_R_pim->SetLineWidth(2);
  bac_pi->Draw();
  
  bac_mu->Draw("same");
  bac_K->Draw("same");
  //bac_R_p->Draw("same");
  //bac_R_pip->Draw("same");
  //bac_R_pim->Draw("same");
  
  
  bac_pi->SetTitle(";#beta;Counts");

  TLine *line1 = new TLine(BAC_thre, 0.5, BAC_thre, 5*1e4);
  line1->SetLineColor(kRed); 
  line1->SetLineStyle(2); 
  line1->SetLineWidth(2); 
  line1->Draw();

  TLatex latex_bac;
  latex_bac.SetTextSize(0.06); 
  latex_bac.SetTextFont(42);
  latex_bac.SetNDC(); 

  latex_bac.SetTextColor(kRed); 
  latex_bac.DrawLatex(0.45, 0.80, "BAC threshold");
  
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->cd();
  gPad->SetLogy();
  kvc_K->SetLineColor(kBlack);
  kvc_pi->SetLineColor(kRed);
  kvc_mu->SetLineColor(kCyan);
  kvc_R_p->SetLineColor(kMagenta);
  kvc_R_pip->SetLineColor(kGreen);
  kvc_R_pim->SetLineColor(kBlue);
  kvc_K->SetLineWidth(2);
  kvc_pi->SetLineWidth(2);
  kvc_mu->SetLineWidth(2);
  kvc_R_p->SetLineWidth(2);
  kvc_R_pip->SetLineWidth(2);
  kvc_R_pim->SetLineWidth(2);
  kvc_pi->Draw();
  
  kvc_mu->Draw("same");
  kvc_R_p->Draw("same");
  kvc_R_pip->Draw("same");
  kvc_R_pim->Draw("same");
  kvc_K->Draw("same");
	       
  kvc_pi->SetTitle(";#beta;Counts");

  TLine *line2 = new TLine(KVC_thre, 0.5, KVC_thre, 5*1e4);
  line2->SetLineColor(kBlue); 
  line2->SetLineStyle(2); 
  line2->SetLineWidth(2); 
  line2->Draw();

  TLatex latex;
  latex.SetTextSize(0.06); 
  latex.SetTextFont(42);
  latex.SetNDC(); 

  // 텍스트를 추가 (좌측 상단)
  latex.SetTextColor(kBlack); 
  latex.DrawLatex(0.15, 0.80, "K^{-} from beam");
  latex.SetTextColor(kRed); 
  latex.DrawLatex(0.15, 0.72, "#pi^{-} from beam");
  latex.SetTextColor(kCyan); 
  latex.DrawLatex(0.15, 0.64, "#mu^{-} from beam");
  latex.SetTextColor(kMagenta); 
  latex.DrawLatex(0.15, 0.56, "p from #Lambda");
  latex.SetTextColor(kBlue); 
  latex.DrawLatex(0.15, 0.48, "#pi^{-} from #Lambda");
  latex.SetTextColor(kGreen); 
  latex.DrawLatex(0.15, 0.40, "#pi^{+} from #Lambda");


  TLatex latex_kvc;
  latex_kvc.SetTextSize(0.06); 
  latex_kvc.SetTextFont(42);
  latex_kvc.SetNDC(); 

  latex_kvc.SetTextColor(kBlue); 
  latex_kvc.DrawLatex(0.5, 0.80, "KVC threshold");
  //latex->Draw("same");
  
}

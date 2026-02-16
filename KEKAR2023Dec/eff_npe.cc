void eff_npe(){

  double x_err = 4.5;
  double y_err = 13/2;

  TFile *file = new TFile("kek_result.root");
  TTree *tree = (TTree*)file->Get("tree");

  int run_out,ped_out,x_out,y_out,bac_hv_out,bac_thre_out,kvc_hv_out,kvc_thre_out,bac_thick_out,kvc_thick_out;
  double bac_npe_out,bac_npe_err_out,bac_eff_out,bac_eff_low_out,bac_eff_hi_out;

  tree->SetBranchAddress("runnum",&run_out);
  tree->SetBranchAddress("pednum",&ped_out);
  tree->SetBranchAddress("x_pos",&x_out);
  tree->SetBranchAddress("y_pos",&y_out);
  tree->SetBranchAddress("bac_HV",&bac_hv_out);
  tree->SetBranchAddress("bac_thre",&bac_thre_out);
  tree->SetBranchAddress("bac_thick",&bac_thick_out);

  tree->SetBranchAddress("bac_npe",&bac_npe_out);
  tree->SetBranchAddress("bac_npe_err",&bac_npe_err_out);
  tree->SetBranchAddress("bac_eff",&bac_eff_out);
  tree->SetBranchAddress("bac_eff_err_low",&bac_eff_low_out);
  tree->SetBranchAddress("bac_eff_err_hi",&bac_eff_hi_out);
    
    
  tree->SetBranchAddress("kvc_HV",&kvc_hv_out);
  tree->SetBranchAddress("kvc_thre",&kvc_thre_out);
  tree->SetBranchAddress("kvc_thick",&kvc_thick_out);

  double total = tree->GetEntries();

  TGraphErrors* bac_npe[2][11];
  TGraphAsymmErrors* bac_eff[2][11];
  int y_pos[2][11]={{-54,-51,-45,-36,-18,0,18,36,45,51,54}
		    ,{-54,-36,-18,0,18,36,51,54,-9999,-9999,-9999}}; //layer2 , layer3


  for(int i=0;i<2;i++){
    for(int j=0;j<11;j++){
      bac_npe[i][j] = new TGraphErrors();
      bac_eff[i][j] = new TGraphAsymmErrors();
    }
  }


  int count = 0;
  double npe = 0;

  for(int n=0;n<total;n++){
    tree->GetEntry(n);
    if(x_out==-64||x_out==64||y_out==-54||y_out==54)continue;
    else if(bac_thick_out==2&&bac_hv_out==58&&bac_thre_out==60){
      if(x_out>0&&x_out<25){
	count++;
	npe+=bac_npe_out;
      }
      for(int i=0;i<11;i++){
	if(y_out==y_pos[0][i]){
	  
	  bac_npe[0][i]->SetPoint(bac_npe[0][i]->GetN(),x_out,bac_npe_out);
	  bac_npe[0][i]->SetPointError(bac_npe[0][i]->GetN()-1,x_err,bac_npe_err_out);

	  bac_eff[0][i]->SetPoint(bac_eff[0][i]->GetN(),x_out,bac_eff_out);
	  bac_eff[0][i]->SetPointError(bac_eff[0][i]->GetN()-1,x_err,x_err,bac_eff_low_out,bac_eff_hi_out);
	  
	  break;
	  
	}
      
      
      }
    }
    else if(bac_thick_out==3&&bac_hv_out==58&&bac_thre_out==75){
      for(int i=0;i<8;i++){
	if(y_out==y_pos[1][i]){
	  bac_npe[1][i]->SetPoint(bac_npe[1][i]->GetN(),x_out,bac_npe_out);
	  bac_npe[1][i]->SetPointError(bac_npe[1][i]->GetN()-1,x_err,bac_npe_err_out);

	  bac_eff[1][i]->SetPoint(bac_eff[1][i]->GetN(),x_out,bac_eff_out);
	  bac_eff[1][i]->SetPointError(bac_eff[1][i]->GetN()-1,x_err,x_err,bac_eff_low_out,bac_eff_hi_out);
	  break;
	  
	}
	
      }
      
    }
  }

  cout<<npe/count<<endl;

  TLegend *le[2][2];
  
  TMultiGraph *mg_bac_npe[2];
  TMultiGraph *mg_bac_eff[2];
  for(int i=0;i<2;i++){
    mg_bac_npe[i] = new TMultiGraph();
    mg_bac_eff[i] = new TMultiGraph();
    le[0][i] = new TLegend(0.8,0.5,0.48,0.6);
    le[1][i] = new TLegend(0.8,0.5,0.48,0.6);
    mg_bac_npe[i]->SetTitle(Form("Layer %d;x [mm];N_{p.e.}",i+2));
    mg_bac_eff[i]->SetTitle(Form("Layer %d;x [mm];Efficiency",i+2));
    for(int j=0;j<11;j++){
      //if(i==1&&j>7)continue;
     
      bac_npe[i][j]->SetMarkerColor(j+1);
      bac_npe[i][j]->SetMarkerStyle(3);

      bac_npe[i][j]->SetLineColor(j+1);

      bac_eff[i][j]->SetMarkerColor(j+1);
      bac_eff[i][j]->SetMarkerStyle(3);

      bac_eff[i][j]->SetLineColor(j+1);

      if(j==9){
	bac_npe[i][j]->SetMarkerColor(38);
	bac_npe[i][j]->SetLineColor(38);
	bac_eff[i][j]->SetMarkerColor(38);
	bac_eff[i][j]->SetLineColor(38);
      }
      mg_bac_npe[i]->Add(bac_npe[i][j]);
      mg_bac_eff[i]->Add(bac_eff[i][j]);
      le[0][i]->AddEntry(bac_npe[i][j],Form("y = %d mm",y_pos[i][j]));
      le[1][i]->AddEntry(bac_eff[i][j],Form("y = %d mm",y_pos[i][j]));
    }
  }


  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(2);
  c1->cd(1);
  c1->SetMargin(0.13,0.05,0.15,0.15);
  mg_bac_npe[0]->Draw("ap");
  le[0][0]->Draw();

  c1->cd(2);
  
  mg_bac_npe[1]->Draw("ap");
  le[0][1]->Draw();
  

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(2);
  c2->cd(1);
  c2->SetMargin(0.13,0.05,0.15,0.15);
  mg_bac_eff[0]->Draw("ap");
  le[1][0]->Draw();

  c2->cd(2);
  mg_bac_eff[1]->Draw("ap");
  le[1][1]->Draw();
  
  
    
      
    
  
}
  

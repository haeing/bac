#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TGraph.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>



using namespace std;



int main(int argc, char* argv[]){

  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <specific_number>" << endl;
    return 1;
  }

  // Get the specific number from the command-line argument
  int targetline = stoi(argv[1]);
  //order of ID
  //KVC : S1,S2,S3,S4,a,b,c,d
  //BAC : S,1,2,3,4
  //Get run information
 int pedrun, runnumber, x_pos, y_pos, bac_hv, bac_thre, BAC_thick, KVC_thick, kvc_hv, kvc_thre;




  ifstream inputFile("run_info.txt");

  if(!inputFile){
    cerr<< "Error opening file." <<endl;
    return 1;
  }

  string line;
  int currentline = 0;
  bool linefind = false;

  while(getline(inputFile,line)){
    if(line.empty() || line[0] == '#'){
      continue;
    }
    
    istringstream iss(line);
    int firstNumber;
    if(iss >> firstNumber){
      if(firstNumber == targetline){
	int param1, param2, param3, param4, param5, param6, param7, param8, param9;
	if (iss >> param1 >> param2 >> param3 >> param4 >> param5 >> param6 >> param7 >> param8 >> param9) {
	   runnumber = targetline;
	   pedrun = param1;
	   x_pos = param2;
	   y_pos = param3;
	   bac_hv = param4;
	   bac_thre = param5;
	   BAC_thick = param6;
	   kvc_hv = param7;
	   kvc_thre = param8;
	   KVC_thick = param9;

	   // Use the parameters in your code
	   cout << "Parameters extracted from line " << currentline << ": "
		<< param1 << " " << param2 << " " << param3 << " "
		<< param4 << " " << param5  <<" " << param6 <<" " << param7<<" "<< param8 << " " << param9 << endl;
	   linefind = true;
	   inputFile.close();
	   break;
	   //return 0;
	 }
	 else {
	   cerr << "Error reading parameters from line " << currentline << endl;
	   // Close the file and exit the program with an error code
	   //inputFile.close();
	   //return 1;
	 }
      }
    }
  }

  // If the specified number is not found
  if(!linefind){
    cerr << "Specified number not found." << endl;
    inputFile.close();
    return 1;
  }


  int bac_ch = 5;
  int kvc_ch = 8;
  double gain_bac[bac_ch];

  gain_bac[1] = 14.7686;
  gain_bac[2] = 14.6856;
  gain_bac[3] = 15.1184;
  gain_bac[4] = 13.1112;

  double pxt_bac[5];
  pxt_bac[0] = 0.0;
  pxt_bac[1] = 0.452525;
  pxt_bac[2] = 0.506778;
  pxt_bac[3] = 0.509973;
  pxt_bac[4] = 0.350078;
  
  double gain_kvc[kvc_ch];
  gain_kvc[0] = 20;
  gain_kvc[1] = 20;
  gain_kvc[2] = 20;
  gain_kvc[3] = 20;
  gain_kvc[4] = 14;
  gain_kvc[5] = 14;
  gain_kvc[6] = 14;
  gain_kvc[7] = 14;


  int range_baca_pe[2];
  range_baca_pe[0] = 100;
  range_baca_pe[1] = 500;

  int range_kvca_pe[2];
  range_kvca_pe[0] = 50;
  range_kvca_pe[1] = 300;

  double range_Tt[2];
  range_Tt[0] = 120000;
  range_Tt[1] = 135000;


  double range_Ta[2];
  range_Ta[0] = 200;
  range_Ta[1] = 1400;

  double range_baca[2];
  range_baca[0] = 100;
  range_baca[1] = 4000;

  double range_bacn[2];
  range_bacn[0] = 0;
  range_bacn[1] = 300;

  double range_bact[2];
  range_bact[0] = 120000;
  range_bact[1] = 140000;

  double range_kvca[2];
  range_kvca[0] = 100;
  range_kvca[1] = 4000;

  double range_kvcn[2];
  range_kvcn[0] = 0;
  range_kvcn[1] = 50;

  double range_kvct[2];
  range_kvct[0] = 300;
  range_kvct[1] = 500;

  TString dir = "../../../data/KEKAR2023Dec/KEKAR2023Dec_root";
  TFile *file_pe = new TFile(Form("%s/kekar_run00%d.root",dir.Data(),pedrun));

  TTree  *data_pe = (TTree*)file_pe->Get("tree");
  double bacsuma_pe_raw[1];
  double baca_pe_raw[bac_ch-1][1];
  double kvcsuma_pe_raw[kvc_ch-4][1];
  double kvca_pe_raw[kvc_ch-4][1];

  double baca_pe[bac_ch][1];
  double kvca_pe[kvc_ch][1];
  
  data_pe->SetBranchAddress("bacsuma",bacsuma_pe_raw);
  data_pe->SetBranchAddress("baca",baca_pe_raw);
  data_pe->SetBranchAddress("kvcsuma",kvcsuma_pe_raw);
  data_pe->SetBranchAddress("kvca",kvca_pe_raw);
		    
  TH1D* hist_baca_pe[bac_ch];
  TF1* fit_baca_pe[bac_ch];
  TH1D* hist_kvca_pe[kvc_ch];
  TF1* fit_kvca_pe[kvc_ch];
  for(int i=0;i<bac_ch;i++){
    hist_baca_pe[i]= new TH1D(Form("hist_baca_pe%d",i),Form("hist_baca_pe%d",i),range_baca_pe[1]-range_baca_pe[0],range_baca_pe[0],range_baca_pe[1]);
    fit_baca_pe[i] = new TF1(Form("fit_baca_pe%d",i),"gaus(0)",range_baca_pe[0],range_baca_pe[1]);
  }
  for(int i=0;i<kvc_ch;i++){
    hist_kvca_pe[i] = new TH1D(Form("hist_kvca_pe%d",i),Form("hist_kvca_pe%d",i),range_kvca_pe[1]-range_kvca_pe[0],range_kvca_pe[0],range_kvca_pe[1]);
    fit_kvca_pe[i]  = new TF1(Form("fit_kvca_pe%d",i),"gaus(0)",range_kvca_pe[0],range_kvca_pe[1]);
  }
  
  double total_pe = data_pe->GetEntries();

  //Pedestal run
  for(int n=0;n<total_pe;n++){
    data_pe->GetEntry(n);
    
    for(int i=0;i<bac_ch;i++){
      if(i==0)baca_pe[0][0] = bacsuma_pe_raw[0];
      else if(i>0)baca_pe[i][0] = baca_pe_raw[i-1][0];
      
      hist_baca_pe[i]->Fill(baca_pe[i][0]);
    }
    for(int i=0;i<kvc_ch;i++){
      if(i<4)kvca_pe[i][0] = kvcsuma_pe_raw[i][0];
      else if(i>=4)kvca_pe[i][0] = kvca_pe_raw[i-4][0];
      
      hist_kvca_pe[i]->Fill(kvca_pe[i][0]);
    }
  }


  TCanvas *c_bac_pe = new TCanvas("c_bac_pe","c_bac_pe");
  c_bac_pe->Divide(3,2);
  
  double pa_baca_pe[bac_ch][3];
  double pa_kvca_pe[kvc_ch][3];
  for(int i=0;i<bac_ch;i++){
    c_bac_pe->cd(i+1);
    hist_baca_pe[i]->Draw();
    hist_baca_pe[i]->Fit(fit_baca_pe[i],"Q","",range_baca_pe[0],range_baca_pe[1]);
    fit_baca_pe[i]->GetParameters(pa_baca_pe[i]);
    
  }
  c_bac_pe->SaveAs(Form("picture/run00%d_BAC_pedestal.pdf",runnumber));

  
  TCanvas *c_kvc_pe = new TCanvas("c_kvc_pe","c_kvc_pe");
  c_kvc_pe->Divide(4,2);
  for(int i=0;i<kvc_ch;i++){
    c_kvc_pe->cd(i+1);
    hist_kvca_pe[i]->Draw();
    hist_kvca_pe[i]->Fit(fit_kvca_pe[i],"Q","",range_kvca_pe[0],range_kvca_pe[1]);
    fit_kvca_pe[i]->GetParameters(pa_kvca_pe[i]);
    
  }
  
  c_kvc_pe->SaveAs(Form("picture/run00%d_KVC_pedestal.pdf",runnumber));



  TFile *file = new TFile(Form("%s/kekar_run00%d.root",dir.Data(),runnumber));
  TTree *data = (TTree*)file->Get("tree");

  double Ta[4][1];
  
  double Tt[5][16];

  double bacsuma_raw[2][1];
  double baca_raw[bac_ch-1];
  double bacsumt_raw[16];
  double bact_raw[bac_ch-1][16];
  
  double baca[bac_ch];
  double bact[16];

  double kvcsuma_raw[kvc_ch/2];
  double kvca_raw[kvc_ch/2];
  double kvcsumt_raw[kvc_ch/2][16];
  double kvct_raw[kvc_ch/2][16];

  double kvca[kvc_ch][1];
  double kvct[kvc_ch][1];
  
  
  data->SetBranchAddress("t1a",Ta[0]);
  data->SetBranchAddress("t2a",Ta[1]);
  data->SetBranchAddress("t3a",Ta[2]);
  data->SetBranchAddress("t4a",Ta[3]);
  
  data->SetBranchAddress("t1t",Tt[1]);
  data->SetBranchAddress("t2t",Tt[2]);
  data->SetBranchAddress("t3t",Tt[3]);
  data->SetBranchAddress("t4t",Tt[4]);
  
  data->SetBranchAddress("bacsuma",bacsuma_raw[1]);
  data->SetBranchAddress("baca",baca_raw);
  data->SetBranchAddress("bacsumt",bact);

  data->SetBranchAddress("kvcsuma",kvcsuma_raw);
  data->SetBranchAddress("kvca",kvca_raw);
  data->SetBranchAddress("kvcsumt",kvcsumt_raw);
  data->SetBranchAddress("kvct",kvct_raw);


  TH1D *hist_Tt[4];
  TF1 *fit_Tt[4];
  TH1D *hist_Ta[4];

  TH1D *hist_baca_ins = new TH1D("hist_baca_ins","hist_baca_ins",(range_baca[1]-range_baca[0]),range_baca[0],range_baca[1]);
  TH1D *hist_baca[bac_ch];
  TH1D *hist_baca_full = new TH1D("hist_baca_full", "hist_baca_full",(range_bacn[1]-range_bacn[0])/2,range_bacn[0],range_bacn[1]);
  TH1D *hist_baca_cut = new TH1D("hist_baca_cut", "hist_baca_cut",(range_bacn[1]-range_bacn[0])/2,range_bacn[0],range_bacn[1]);
  TF1 *fit_baca = new TF1("fit_baca","gaus(0)",range_bacn[0],range_bacn[1]);

  TF1 *fit_baca_ind = new TF1("fit_baca_ind","gaus(0)",range_baca[0],range_baca[1]);
  
  for(int i=0;i<bac_ch;i++){
    hist_baca[i] = new TH1D(Form("hist_baca%d",i),Form("hist_baca%d",i),(range_bacn[1]-range_bacn[0])/2,range_bacn[0],range_bacn[1]);
    
  }
  
  TH1D *hist_bact = new TH1D("hist_bact","hist_bact",(range_bact[1]-range_bact[0])/5,range_bact[0],range_bact[1]);

  
  TH1D *hist_bact_cut = new TH1D("hist_bact_cut","hist_bact_cut",(range_bact[1]-range_bact[0])/2,range_bact[0],range_bact[1]);

  TH1D *hist_bact_thre = new TH1D("hist_bact_thre","hist_bact_thre",(range_bact[1]-range_bact[0])/2,range_bact[0],range_bact[1]);

  TH2D *hist_bac_sum_ind = new TH2D("hist_bac_sum_ind","hist_bac_sum_ind",(range_baca[1]-range_baca[0])/2,range_baca[0],range_baca[1],(range_baca[1]-range_baca[0])/2,range_baca[0],range_baca[1]);

  TF1 *fit_sum_ind= new TF1("fit_sum_ind","[0]*x+[1]",range_baca[0],range_baca[1]);

  TF1 *fit_bact = new TF1("fit_bact","gaus(0)",range_bact[0],range_bact[1]);
  
  TH1D *hist_kvca[kvc_ch];
  TH1D *hist_kvct[kvc_ch/2];
  TH1D *hist_kvca_cut[kvc_ch];
  TH1D *hist_kvct_cut[kvc_ch/2];
  TH1D *hist_kvct_thre[kvc_ch/2];
  TH2D *hist_kvc_sum_ind = new TH2D("hist_kvc_sum_ind%d","hist_kvc_sum_ind",(range_kvca[1]-range_kvca[0])/2,range_kvca[0],range_kvca[1],(range_kvcn[1]-range_kvcn[0])/2,range_kvcn[0],range_kvcn[1]);

  TF1 *fit_kvca[kvc_ch];
  TF1 *fit_kvct[kvc_ch/2];
  
  for(int i=0;i<kvc_ch;i++){
    
    hist_kvca[i] = new TH1D(Form("hist_kvca%d",i),Form("hist_kvca%d",i),(range_kvca[1]-range_kvca[0])/2,range_kvca[0],range_kvca[1]);
    hist_kvca_cut[i] = new TH1D(Form("hist_kvca_cut%d",i),Form("hist_kvca_cut%d",i),(range_kvca[1]-range_kvca[0])/2,range_kvca[0],range_kvca[1]);
    fit_kvca[i] = new TF1(Form("fit_kvca%d",i),"gaus(0)",range_kvca[0],range_kvca[1]);
  }
  
  for(int i=0;i<kvc_ch/2;i++){
    hist_Tt[i] = new TH1D(Form("hist_Tt%d",i),Form("hist_Tt%d",i),(range_Tt[1]-range_Tt[0])/50,range_Tt[0],range_Tt[1]);
    hist_Ta[i] = new TH1D(Form("hist_Ta%d",i),Form("hist_Ta%d",i),(range_Ta[1]-range_Ta[0])/5,range_Ta[0],range_Ta[1]);
    fit_Tt[i] = new TF1(Form("fit_Tt%d",i),"gaus(0)",range_Tt[0],range_Tt[1]);
    hist_kvct[i] = new TH1D(Form("hist_kvct%d",i),Form("hist_kvct%d",i),(range_kvct[1]-range_kvct[0])/2,range_kvct[0],range_kvct[1]);
    hist_kvct_cut[i] = new TH1D(Form("hist_kvct_cut%d",i),Form("hist_kvct_cut%d",i),(range_kvct[1]-range_kvct[0])/2,range_kvct[0],range_kvct[1]);
    hist_kvct_thre[i] = new TH1D(Form("hist_kvct_thre%d",i),Form("hist_kvct_thre%d",i),(range_kvct[1]-range_kvct[0])/2,range_kvct[0],range_kvct[1]);
    fit_kvct[i] = new TF1(Form("fit_kvct%d",i),"gaus(0)",range_kvct[0],range_kvct[1]);
  }

  double total = data->GetEntries();
  
  //Total raw
  double bac_adc_summing;
  double kvc_adc_summing;
  
  for(int n=0;n<total;n++){
    data->GetEntry(n);
    bac_adc_summing = 0;
    for(int i=0;i<4;i++){
      hist_Tt[i]->Fill(Tt[i+1][0]);
      hist_Ta[i]->Fill(Ta[i][0]);
    }
    for(int j=0;j<4;j++){
      bac_adc_summing +=baca_raw[j]-pa_baca_pe[j+1][1];
    }
    hist_baca_ins->Fill(bac_adc_summing);
    hist_bac_sum_ind->Fill(bac_adc_summing,(bacsuma_raw[1][0]-pa_baca_pe[0][1]));
  }
  
  TCanvas *c_trig = new TCanvas("c_trig","c_trig");
  c_trig->Divide(4,2);
  


  double pa_Tt[4][3];
  for(int i=0;i<4;i++){
    c_trig->cd(i+1);
    hist_Tt[i]->Draw();
    
    
    hist_Tt[i]->Fit(fit_Tt[i],"Q","",range_Tt[0],range_Tt[1]);
    fit_Tt[i]->GetParameters(pa_Tt[i]);
    c_trig->cd(i+5);
    hist_Ta[i]->Draw();

  }

  c_trig->SaveAs(Form("picture/run00%d_Trigger.pdf",runnumber));

  double pa_baca_ind[3];
  hist_baca_ins->Fit(fit_baca_ind,"Q","",range_baca[0],range_baca[1]);
  fit_baca_ind->GetParameters(pa_baca_ind);

  TCanvas *c_sum_ind = new TCanvas("c_sum_ind","c_sum_ind");
  double pa_bac_sum_ind[2];
  c_sum_ind->Divide(2);
  c_sum_ind->cd(1);
  hist_bac_sum_ind->Draw();
  hist_bac_sum_ind->Fit(fit_sum_ind,"","colz",range_baca[0],range_baca[1]);
  fit_sum_ind->GetParameters(pa_bac_sum_ind);

  c_sum_ind->cd(2);
  hist_baca_ins->Draw();

  c_sum_ind->SaveAs(Form("picture/run00%d_SUM_Ind.pdf",runnumber));

  gain_bac[0] = (gain_bac[1]+gain_bac[2]+gain_bac[3]+gain_bac[4])/4*(pa_bac_sum_ind[0]+pa_bac_sum_ind[1]/pa_baca_ind[2]);
  for(int i=0;i<4;i++)pxt_bac[0] += pxt_bac[i+1];
  pxt_bac[0] /= 4.;
  std::cout<<"pxt_bac : "<<pxt_bac[0]<<std::endl;
  


  

  
  int trig_pass;


  
  //T1&T2&T3&T4 cut
  
  for(int n=0;n<total;n++){
    data->GetEntry(n);
    trig_pass=0;
    for(int i=0;i<4;i++){
      if(Tt[i+1][0]>(pa_Tt[i][1]-5*pa_Tt[i][2])&&Tt[i+1][0]<(pa_Tt[i][1]+5*pa_Tt[i][2])){
	//if(Ta[i][0]>400&&Ta[i][0]<600)trig_pass++; tight cut
	if(Ta[i][0]>400)trig_pass++;
      }
    }
    if(trig_pass==4){
      bac_adc_summing = 0;
      kvc_adc_summing = 0;
      for(int bac=0;bac<bac_ch;bac++){
	if(bac==0){
	  baca[0] = bacsuma_raw[1][0];
	  
	}
	else if(bac>0){
	  baca[bac] = baca_raw[bac-1];
	}
	
	
	hist_baca[bac]->Fill((baca[bac]-pa_baca_pe[bac][1])/gain_bac[bac]*pxt_bac[bac]);

      }
      hist_bact->Fill(bact[0]);

      

      for(int kvc=0;kvc<kvc_ch;kvc++){
	if(kvc<4){
	  kvca[kvc][0] = kvcsuma_raw[kvc];
	  for(int j=0;j<16;j++)kvct[kvc][j] = kvcsumt_raw[kvc][j];
	}
	else if(kvc>=4){
	  kvca[kvc][0] = kvca_raw[kvc-4];
	  for(int j=0;j<16;j++)kvct[kvc][j] = kvct_raw[kvc-4][j];
	}
	hist_kvca[kvc]->Fill((kvca[kvc][0]-pa_kvca_pe[kvc][1]));
	hist_kvct[kvc]->Fill(kvct[kvc][0]);
      }
      if(KVC_thick==1){
	hist_kvc_sum_ind->Fill((kvca[1][0]-pa_kvca_pe[1][1]),(kvca[4][0]-pa_kvca_pe[4][1])/gain_kvc[4]+(kvca[5][0]-pa_kvca_pe[5][1])/gain_kvc[5]);
      }
      else if(KVC_thick==2){
	for(int i=4;i<8;i++)kvc_adc_summing+=(kvca[i][0]-pa_kvca_pe[i][1])/gain_kvc[i];
	hist_kvc_sum_ind->Fill((kvca[1][0]-pa_kvca_pe[1][1]),kvc_adc_summing);
      }
    }
  }
			    


  double pa_baca[3];
  double pa_bact[3];
  
  double pa_kvca[kvc_ch][3];
  double pa_kvct[kvc_ch/2][3];
  
  for(int i=0;i<kvc_ch/2;i++){
    hist_kvct[i]->Fit(fit_kvct[i],"Q0","",range_kvct[0],range_kvct[1]);
    fit_kvct[i]->GetParameters(pa_kvct[i]);
  }
  hist_bact->Fit(fit_bact,"Q","",range_bact[0],range_bact[1]);
  fit_bact->GetParameters(pa_bact);

  
  

  TCanvas *c_bact = new TCanvas("c_bact","c_bact");

  TLine *bac_cutt1 = new TLine(pa_bact[1]-5*pa_bact[2],0,pa_bact[1]-5*pa_bact[2],hist_bact->GetMaximum());
  bac_cutt1->SetLineColor(kRed);
  bac_cutt1->SetLineWidth(2);
  bac_cutt1->SetLineStyle(9);

  TLine *bac_cutt2 = new TLine(pa_bact[1]+5*pa_bact[2],0,pa_bact[1]+5*pa_bact[2],hist_bact->GetMaximum());
  bac_cutt2->SetLineColor(kRed);
  bac_cutt2->SetLineWidth(2);
  bac_cutt2->SetLineStyle(9);
  
  hist_bact->SetTitle("BAC TDC cut condition;TDC [ch];Counts");
  hist_bact->Draw();
  bac_cutt1->Draw();
  bac_cutt2->Draw();
  
  c_bact->SaveAs(Form("picture/run00%d_bact.pdf",runnumber));
  
  
  
  
  TEfficiency *eff_bac = new TEfficiency("eff_bac","",1,-100,100);
  Bool_t bac_pass = 0;

  //T1&T2&T3&T4&detector cut
  for(int n=0;n<total;n++){
    data->GetEntry(n);
    trig_pass=0;
    bac_pass = 0;
    for(int i=0;i<4;i++){
      if(Tt[i+1][0]>(pa_Tt[i][1]-5*pa_Tt[i][2])&&Tt[i+1][0]<(pa_Tt[i][1]+5*pa_Tt[i][2])){
	//if(Ta[i][0]>400&&Ta[i][0]<600)trig_pass++; //tight cut
	if(Ta[i][0]>400)trig_pass++;
      }
    }
    if(trig_pass==4){
      hist_baca_full->Fill((baca[0]-pa_baca_pe[0][1])/gain_bac[0]*pxt_bac[0]);
      if(trig_pass==4&&bact[0]>(pa_bact[1]-5*pa_bact[2])&&bact[0]<(pa_bact[1]+5*pa_bact[2])){
	bac_pass = 1;
	for(int bac=0;bac<bac_ch;bac++){
	  if(bac==0){
	    baca[0] = bacsuma_raw[1][0];
	  }
	  else if(bac>0){
	    baca[bac] = baca_raw[bac-1];
	  }
	}
	hist_baca_cut->Fill((baca[0]-pa_baca_pe[0][1])/gain_bac[0]*pxt_bac[0]);
	hist_bact_cut->Fill(bact[0]);
      }
      eff_bac->Fill(bac_pass,x_pos);
      /*
	if(trig_pass==4&&kvct[0]>pa_kvct[1]-5*pa_kvct[2]&&kvct[0]<pa_kvct[1]+5*pa_kvct[2]){
	for(int kvc=0;kvc<kvc_ch;kvc++)hist_kvca_cut[kvc]->Fill((kvca[kvc][0]-pa_kvca_pe[kvc][1])/gain_kvc[kvc]);
	for(int kvc=0;kvc<kvc_ch/2;kvc++)hist_kvct_cut[kvc]->Fill(kvct[kvc]);
	}
      */
      
    }
  }

  for(int i=0;i<kvc_ch;i++){
    hist_kvca_cut[i]->Fit(fit_kvca[i],"Q0","",range_kvca[0],range_kvca[1]);
    fit_kvca[i]->GetParameters(pa_kvca[i]);
  }

  hist_baca_cut->Fit(fit_baca,"Q","",range_bacn[0],range_bacn[1]);
  fit_baca->GetParameters(pa_baca);
  double pa_baca_err;
  pa_baca_err = fit_baca->GetParError(1);


  TCanvas *c_bac_npe = new TCanvas("c_bac_npe","c_bac_npe");
  hist_baca[0]->Draw();
  hist_baca_cut->SetFillColorAlpha(kBlue, 0.35);
  hist_baca_cut->Draw("same");x
  c_bac_npe->SaveAs(Form("picture/run00%d_bac_npe.pdf",runnumber));

  
  TFile *output = new TFile("kek_result.root","UPDATE");
  if (!output || output->IsZombie()){
    //if (output->IsZombie()){
    cerr << "Error opening or creating output.root file." << endl;
    return 1;
  }

  TTree *tree_out = dynamic_cast<TTree*>(output->Get("tree"));

  
  
  if(!tree_out){

    tree_out = new TTree("tree","tree");
  

    int runnum, pednum, x_pos, y_pos, bac_HV, bac_thre, bac_thick, kvc_HV, kvc_thre, kvc_thick;
    double bac_npe, bac_npe_err, bac_eff, bac_eff_err_low, bac_eff_err_hi;
    
    tree_out->Branch("runnum",&runnum,"runnum/I");
    tree_out->Branch("pednum",&pednum,"pednum/I");
    tree_out->Branch("x_pos",&x_pos,"x_pos/I");
    tree_out->Branch("y_pos",&x_pos,"y_pos/I");
    tree_out->Branch("bac_HV",&bac_HV,"bac_HV/I");
    tree_out->Branch("bac_thre",&bac_thre,"bac_thre/I");
    tree_out->Branch("bac_thick",&bac_thick,"bac_thick/I");

    tree_out->Branch("bac_npe",&bac_npe,"bac_npe/D");
    tree_out->Branch("bac_npe_err",&bac_npe_err,"bac_npe_err/D");
    tree_out->Branch("bac_eff",&bac_eff,"bac_eff/D");
    tree_out->Branch("bac_eff_err_low",&bac_eff_err_low,"bac_eff_err_low/D");
    tree_out->Branch("bac_eff_err_hi",&bac_eff_err_hi,"bac_eff_err_hi/D");


    tree_out->Branch("kvc_HV",&kvc_HV,"kvc_HV/I");
    tree_out->Branch("kvc_thre",&kvc_thre,"kvc_thre/I");
    tree_out->Branch("kvc_thick",&kvc_thick,"kvc_thick/I");
    /*
      tree_out->Branch("kvc_npe",&kvc_npe_out,"kvc_npe/D");
      tree_out->Branch("kvc_npe_err",&kvc_npe_errout,"kvc_npe_err/D");
      tree_out->Branch("kvc_eff",&kvc_eff_out,"kvc_eff/D");
      tree_out->Branch("kvc_eff_err_low",&kvc_eff_err_low_out,"kvc_eff_err_low/D");
      tree_out->Branch("kvc_eff_err_hi",&kvc_eff_err_hi_out,"kvc_eff_err_hi/D");
    */
  }

  int run_out,ped_out,x_out,y_out,bac_hv_out,bac_thre_out,kvc_hv_out,kvc_thre_out,bac_thick_out,kvc_thick_out;
  double bac_npe_out,bac_npe_err_out,bac_eff_out,bac_eff_low_out,bac_eff_hi_out;

  run_out = runnumber;
  ped_out = pedrun;
  x_out = x_pos;
  y_out = y_pos;
  bac_hv_out = bac_hv;
  bac_thre_out = bac_thre;
  bac_thick_out = BAC_thick;
  bac_npe_out = pa_baca[1];
  bac_npe_err_out = pa_baca_err;
  bac_eff_out = eff_bac->GetEfficiency(1);
  bac_eff_low_out = eff_bac->GetEfficiencyErrorLow(1);
  bac_eff_hi_out = eff_bac->GetEfficiencyErrorUp(1);

  kvc_hv_out = kvc_hv;
  kvc_thre_out = kvc_thre;
  kvc_thick_out = KVC_thick;
  
  tree_out->SetBranchAddress("runnum",&run_out);
  tree_out->SetBranchAddress("pednum",&ped_out);
  tree_out->SetBranchAddress("x_pos",&x_out);
  tree_out->SetBranchAddress("y_pos",&y_out);
  tree_out->SetBranchAddress("bac_HV",&bac_hv_out);
  tree_out->SetBranchAddress("bac_thre",&bac_thre_out);
  tree_out->SetBranchAddress("bac_thick",&bac_thick_out);

  tree_out->SetBranchAddress("bac_npe",&bac_npe_out);
  tree_out->SetBranchAddress("bac_npe_err",&bac_npe_err_out);
  tree_out->SetBranchAddress("bac_eff",&bac_eff_out);
  tree_out->SetBranchAddress("bac_eff_err_low",&bac_eff_low_out);
  tree_out->SetBranchAddress("bac_eff_err_hi",&bac_eff_hi_out);
    
    
  tree_out->SetBranchAddress("kvc_HV",&kvc_hv_out);
  tree_out->SetBranchAddress("kvc_thre",&kvc_thre_out);
  tree_out->SetBranchAddress("kvc_thick",&kvc_thick_out);
  

    

  tree_out->Fill();
  tree_out->Write("",TObject::kOverwrite);
  output->Close();

	   
  
  
  

  



  
  
  

  
  
}

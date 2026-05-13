#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>

bool kaon = false;

int npe_threshold = 5;
const double T0_z = -1100.0;
const double BH2_z = -595.7;
const double BAC_x = -7.; //until 11/10
//const double BAC_x = -17.25; //from 11/11
const double BAC_z = -498.3;

const int NumOfSegT0 = 5;
const int NumOfSegBH2 = 15;
const int NumOfSegBAC = 4;                                                                    

const double t0_tdc_min = 692000.;
const double t0_tdc_max = 697000.;

const double bh2_tdc_min = 705000.;
const double bh2_tdc_max = 720000.;

const double bh2_x = 14.;
const double bh2_y = 100.;

const double bac_tdc_min = 700000.;
const double bac_tdc_max = 720000.;

const double tdc_step = 100.;

struct WStat {
  double sumw  = 0.0;
  double sumwx = 0.0;

  void Fill(double x, double err){
    if (!std::isfinite(x) || !std::isfinite(err)) return;
    if (err <= 0) return;
    const double w = 1.0/(err*err);
    sumw  += w;
    sumwx += w*x;
  }

  bool Has() const { return sumw > 0; }
  double Mean() const { return sumwx/sumw; }
  double Err()  const { return std::sqrt(1.0/sumw); }
};

struct Key {
  int board;
  int ch;
  bool operator<(const Key& o) const {
    if (board != o.board) return board < o.board;
    return ch < o.ch;
  }
};

static std::vector<TString> SplitCSV(const TString& line){
  std::vector<TString> out;
  TObjArray* arr = line.Tokenize(",");
  for (int i=0;i<arr->GetEntriesFast();++i)
    out.push_back(((TObjString*)arr->At(i))->GetString().Strip(TString::kBoth));
  arr->Delete(); delete arr;
  return out;
}

static double FindLandauRange(TH1* h)
{
  TSpectrum spec(3);
  int nfound = spec.Search(h, 2, "", 0.05);
  double *xpeaks = spec.GetPositionX();
  vector<double> peaks;
  for(int j=0;j<nfound;j++)
    peaks.push_back(xpeaks[j]);
  sort(peaks.begin(), peaks.end());
  double peak = peaks[0];

  return peak+30;
}

static double FindGausPeak(TH1 *h)
{
  TSpectrum spec(3);
  int nfound = spec.Search(h, 2, "", 0.1);

  double *xpeaks = spec.GetPositionX();

  double maxPeakX = -1;
  double maxHeight = -1;

  for(int i=0;i<nfound;i++)
    {
      double x = xpeaks[i];
      int bin = h->GetXaxis()->FindBin(x);
      double height = h->GetBinContent(bin);

      if(height > maxHeight)
	{
	  maxHeight = height;
	  maxPeakX = x;
	}
    }
  return maxPeakX;
}

void DrawLine(double min, TH1* h)
{
  TLine *line = new TLine(min, 0, min, h->GetMaximum());
  line->SetLineColor(kBlue);
  line->SetLineWidth(2);
  line->Draw("same");
}

void analysis_e72(int runnumber, int runnumber_ped)
{
  double run_mom;
  double bac_tdc_cut[2]={700000.,720000.};
  int n_beam = 0;
  
  if(runnumber <2100){
    run_mom = 735;
    bac_tdc_cut[0] = 704000;
    bac_tdc_cut[1] = 718000;
    
  }
  else if(runnumber == 2489){
    run_mom = 1000;
    bac_tdc_cut[0] = 704000;
    bac_tdc_cut[1] = 718000;
    
  }
  else if(runnumber == 2502){
    run_mom = 814;
    bac_tdc_cut[0] = 704000;
    bac_tdc_cut[1] = 718000;
  }
  else if(runnumber == 2509){
    run_mom = 645;
    bac_tdc_cut[0] = 704000;
    bac_tdc_cut[1] = 718000;
    
  }
  else if(runnumber == 2512){
    run_mom = 400;
    bac_tdc_cut[0] = 695000;
    bac_tdc_cut[1] = 718000;
  }
  else if(runnumber == 2514){
    run_mom = 300;
    bac_tdc_cut[0] = 695000;
    bac_tdc_cut[1] = 718000;
  }

  
  gROOT->SetBatch(kTRUE);
  double bac_ped_mean[NumOfSegBAC]={0.};
  double bh2_adc_cut[2][NumOfSegBH2]={0.};
  double bh2_tdc_cut[NumOfSegBH2][2]={0.};
  double bac_gain[NumOfSegBAC]={0.};
  double bac_pxt[NumOfSegBAC]={0.};
  
  bac_gain[0] = 12.;
  bac_gain[1] = 11.;
  bac_gain[2] = 9.8;
  bac_gain[3] = 11;
  bac_pxt[0] = 0.48;
  bac_pxt[1] = 0.49;
  bac_pxt[2] = 0.44;
  bac_pxt[3] = 0.40;
  

  
  TString dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/bac_study";
  TFile *file = new TFile(Form("%s/run0%d_Hodoscope.root",dir.Data(),runnumber));
  TTree *data = (TTree*)file->Get("hodo");

  TFile *file_ped = new TFile(Form("%s/run0%d_Hodoscope.root",dir.Data(),runnumber_ped));
  TTree *data_ped = (TTree*)file_ped->Get("hodo");

  TFile *file_bcout = new TFile(Form("%s/run0%d_BcOutTracking.root",dir.Data(),runnumber));
  TTree *bcout = (TTree*)file_bcout->Get("bcout");

  vector<double>* bac_adc_u_ped = nullptr;
  vector<vector<double>>* bac_tdc_u_ped = nullptr;

  data_ped->SetBranchAddress("bac_adc_u",&bac_adc_u_ped);
  data_ped->SetBranchAddress("bac_tdc_u",&bac_tdc_u_ped);

  double btof0;
  
  
  vector<double> *bh2_adc_u = nullptr;
  vector<double> *bh2_adc_d = nullptr;

  vector<vector<double>>* bh2_tdc_s = nullptr;
  
  
  vector<double>* bac_adc_u = nullptr;
  vector<vector<double>>* bac_tdc_u = nullptr;

  data->SetBranchAddress("btof0",&btof0);
  

  data->SetBranchAddress("bh2_adc_u",&bh2_adc_u);
  data->SetBranchAddress("bh2_adc_d",&bh2_adc_d);
  data->SetBranchAddress("bh2_tdc_s",&bh2_tdc_s);

  data->SetBranchAddress("bac_adc_u",&bac_adc_u);
  data->SetBranchAddress("bac_tdc_u",&bac_tdc_u);

  vector<double>* x0 = nullptr;
  vector<double>* y0 = nullptr;
  vector<double>* u0 = nullptr;
  vector<double>* v0 = nullptr;

  int ntrack;

  bcout->SetBranchAddress("ntrack",&ntrack);
  bcout->SetBranchAddress("x0",&x0);
  bcout->SetBranchAddress("y0",&y0);
  bcout->SetBranchAddress("u0",&u0);
  bcout->SetBranchAddress("v0",&v0);

  //Pedestal Study
  TH1D *hist_bac_adc_s_ped = new TH1D("hist_bac_adc_s_ped","hist_bac_adc_s_ped",200,0,1000);
  TF1 *f_bac_adc_s_ped = new TF1("f_bac_adc_s_ped","gaus",100,1000);
  

  for(int n=0;n<100000;n++){
    data_ped->GetEntry(n);
    hist_bac_adc_s_ped->Fill((*bac_adc_u_ped)[4]);
  }

  TString out_pdf = Form("result/run%d.pdf",runnumber);
  TCanvas* c1 = new TCanvas("c1","c1");
  TDatime now;
  TString datetime = Form("%04d-%02d-%02d  %02d:%02d:%02d",
                        now.GetYear(),
                        now.GetMonth(),
                        now.GetDay(),
                        now.GetHour(),
                        now.GetMinute(),
			now.GetSecond());
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.04);

  text.DrawLatex(0.3,0.7,"Run summary");
  text.DrawLatex(0.3,0.6,Form("Run number : %d", runnumber));
  text.DrawLatex(0.3,0.5,Form("Date : %s", datetime.Data()));
  
  c1->Print(out_pdf +"(");

  c1->Clear();

  
  double peak = FindGausPeak(hist_bac_adc_s_ped);
  f_bac_adc_s_ped->SetRange(peak-50,peak+50);
  hist_bac_adc_s_ped->Fit(f_bac_adc_s_ped,"RQ");
  hist_bac_adc_s_ped->Draw();
  double bac_ped_mean_s = f_bac_adc_s_ped->GetParameter(1);

  c1->Print(out_pdf);

  
  

  

  TH1D *hist_bh2_adc_u[NumOfSegBH2];
  TH1D *hist_bh2_adc_d[NumOfSegBH2];
  TH1D *hist_bh2_tdc_s[NumOfSegBH2];

  TF1 *f_bh2_adc_u[NumOfSegBH2];
  TF1 *f_bh2_adc_d[NumOfSegBH2];
  TF1 *f_bh2_tdc_s[NumOfSegBH2];

  TH1D *hist_bac_npe[NumOfSegBAC];
  TH1D *hist_bac_npe_s = new TH1D("hist_bac_npe_s","hist_bac_npe_s",110,-100,1000);
  TH1D *hist_bac_npe_s_total = new TH1D("hist_bac_npe_s_total","hist_bac_npe_s_total",110,-100,1000);
  TH1D *hist_bac_npe_s_bh2[NumOfSegBH2];
  TH1D *hist_bac_npe_s_bh2_pass[NumOfSegBH2];
  TH1D *hist_bac_npe_s_pass = new TH1D("hist_bac_npe_s_pass","hist_bac_npe_s_pass",110,-100,1000);
  TH1D *hist_bac_tdc_s = new TH1D("hist_bac_tdc_s","hist_bac_tdc_s",(bac_tdc_max - bac_tdc_min)/tdc_step,bac_tdc_min,bac_tdc_max);

  TF1 *f_bac_npe_s_bh2[NumOfSegBH2];

  TH2D *hist_bcout_bh2_2d[NumOfSegBH2];
  TH2D *hist_bcout_bac_2d[NumOfSegBH2];
  TH1D *hist_bcout_bac[NumOfSegBH2];
  
  
  for(int i=0;i<NumOfSegBH2;i++){
    hist_bh2_adc_u[i] = new TH1D(Form("hist_bh2_adc_u%d",i),Form("hist_bh2_adc_u%d",i),900,100,1000);
    hist_bh2_adc_d[i] = new TH1D(Form("hist_bh2_adc_d%d",i),Form("hist_bh2_adc_d%d",i),900,100,1000);
    hist_bh2_tdc_s[i] = new TH1D(Form("hist_bh2_tdc_s%d",i),Form("hist_bh2_tdc_s%d",i),(bh2_tdc_max - bh2_tdc_min)/tdc_step,bh2_tdc_min,bh2_tdc_max);
    hist_bac_npe_s_bh2[i] = new TH1D(Form("hist_bac_npe_s_bh2%d",i),Form("hist_bac_npe_s_bh2%d",i),210,-100,2000);
    hist_bac_npe_s_bh2_pass[i] = new TH1D(Form("hist_bac_npe_s_bh2_pass%d",i),Form("hist_bac_npe_s_bh2_pass%d",i),210,-100,2000);
    hist_bcout_bac[i] = new TH1D(Form("hist_bcout_bac%d",i),Form("hist_bcout_bac%d",i),100,-150,150);
    hist_bcout_bh2_2d[i] = new TH2D(Form("hist_bcout_bh2_2d%d",i),Form("hist_bcout_bh2_2d%d",i),100,-150,150,100,-150,150);
    hist_bcout_bac_2d[i] = new TH2D(Form("hist_bcout_bac_2d%d",i),Form("hist_bcout_bac_2d%d",i),100,-150,150,100,-150,150);
    
    f_bh2_adc_u[i] = new TF1(Form("f_bh2_adc_u%d",i),"landau",100,1000);
    f_bh2_adc_d[i] = new TF1(Form("f_bh2_adc_d%d",i),"landau",100,1000);
    f_bh2_tdc_s[i] = new TF1(Form("f_bh2_tdc_s%d",i),"gaus",bh2_tdc_min,bh2_tdc_max);
    f_bac_npe_s_bh2[i] = new TF1(Form("f_bac_npe_s_bh2%d",i),"gaus",0,2000);
  }

  for(int i=0;i<NumOfSegBAC;i++){
    hist_bac_npe[i] = new TH1D(Form("hist_bac_npe%d",i),Form("hist_bac_npe%d",i),70,-20,50);
  }

  for(int n=0;n<data->GetEntries();n++){
    data->GetEntry(n);
    if(kaon){
      if(btof0 >-3)continue;
    }
    else if(!kaon){
      if(btof0 <-3)continue;
    }
    if(n%10000 == 0)cout<<"Entry "<<n<<std::endl;
    
    for(int i=0;i<NumOfSegBH2;i++){
      hist_bh2_adc_u[i]->Fill((*bh2_adc_u)[i]);
      hist_bh2_adc_d[i]->Fill((*bh2_adc_d)[i]);
      for(int j=0;j<(*bh2_tdc_s)[i].size();j++){
	hist_bh2_tdc_s[i]->Fill((*bh2_tdc_s)[i][j]);
      }
    }
    for(int j=0;j<(*bac_tdc_u)[4].size();j++)hist_bac_tdc_s->Fill((*bac_tdc_u)[4][j]);
    hist_bac_npe_s->Fill((*bac_adc_u)[4]);
  }
  
  

  double mpv = 0;
  double sigma = 0;
  
  
  
  c1->Clear();
  
  c1->Divide(4,4);
  for(int i=0;i<NumOfSegBH2;i++){
    c1->cd(i+1);
    gPad->SetLogy();
    hist_bh2_adc_u[i]->Draw();
    double fit_min = FindLandauRange(hist_bh2_adc_u[i]);
    f_bh2_adc_u[i]->SetRange(fit_min,fit_min+1000);
    hist_bh2_adc_u[i]->Fit(f_bh2_adc_u[i],"RQ");
    mpv   = f_bh2_adc_u[i]->GetParameter(1);
    sigma = f_bh2_adc_u[i]->GetParameter(2);
    DrawLine(mpv-sigma*2,hist_bh2_adc_u[i]);
    bh2_adc_cut[0][i] = mpv-sigma*2;
  }
  c1->Print(out_pdf);
  c1->Clear();
  
  c1->Divide(4,4);
  for(int i=0;i<NumOfSegBH2;i++){
    c1->cd(i+1);
    gPad->SetLogy();
    hist_bh2_adc_d[i]->Draw();
    double fit_min = FindLandauRange(hist_bh2_adc_d[i]);
    f_bh2_adc_d[i]->SetRange(fit_min,fit_min+1000);
    hist_bh2_adc_d[i]->Fit(f_bh2_adc_d[i],"RQ");
    mpv   = f_bh2_adc_d[i]->GetParameter(1);
    sigma = f_bh2_adc_d[i]->GetParameter(2);
    DrawLine(mpv-sigma*2,hist_bh2_adc_d[i]);
    bh2_adc_cut[1][i] = mpv-sigma*2;
  }
  
  c1->Print(out_pdf);
  c1->Clear();
  c1->Divide(4,4);
  double mean;
  for(int i=0;i<NumOfSegBH2;i++){
    c1->cd(i+1);
    hist_bh2_tdc_s[i]->Draw();
    double peak = FindGausPeak(hist_bh2_tdc_s[i]);

    f_bh2_tdc_s[i]->SetRange(peak-300,peak+300);
    hist_bh2_tdc_s[i]->Fit(f_bh2_tdc_s[i],"RQ");
    mean = f_bh2_tdc_s[i]->GetParameter(1);
    sigma = f_bh2_tdc_s[i]->GetParameter(2);
    DrawLine(mean-3*sigma,hist_bh2_tdc_s[i]);
    DrawLine(mean+3*sigma,hist_bh2_tdc_s[i]);
    bh2_tdc_cut[i][0] = mean-3*sigma;
    bh2_tdc_cut[i][1] = mean+3*sigma;
  }
  c1->Print(out_pdf);

  c1->Clear();
  c1->Divide(3,2);
  for(int i=0;i<NumOfSegBAC;i++){
    c1->cd(i+1);
    gPad->SetLogy();
    hist_bac_npe[i]->Draw();
  }
  c1->cd(5);
  hist_bac_npe_s->Draw();
  c1->cd(6);
  hist_bac_tdc_s->Draw();
  c1->Print(out_pdf);

  //Efficiency, Npe check
  
  bool bh2_pass[NumOfSegBH2] = {false};
  bool bac_pass[NumOfSegBAC] = {false};
  int eff_total[NumOfSegBH2] = {0};
  int eff_pass[NumOfSegBH2] = {0};
  for(int n=0;n<data->GetEntries();n++){
    if(n%10000 == 0)cout<<"Entry "<<n<<std::endl;
    data->GetEntry(n);
    bcout->GetEntry(n);
    if(kaon){
      if(btof0 >-3)continue;
    }

    //BH2 cut start w/ BcOut
    for(int i=0;i<NumOfSegBH2;i++){
      bh2_pass[i] = false;
      if((*bh2_adc_u)[i]>bh2_adc_cut[0][i] && (*bh2_adc_d)[i]>bh2_adc_cut[1][i]){
	for(int j=0;j<(*bh2_tdc_s).size();j++){
	  if((*bh2_tdc_s)[i][j]>bh2_tdc_cut[i][0] && (*bh2_tdc_s)[i][j]<bh2_tdc_cut[i][1]){
	    double bh2_seg_x_min = -1*NumOfSegBH2*bh2_x/2.+i*bh2_x;
	    double bh2_seg_x_max = -1*NumOfSegBH2*bh2_x/2.+(i+1)*bh2_x;
	    if(ntrack>0){
	      if((*x0)[0]+BH2_z*(*u0)[0] > bh2_seg_x_min && (*x0)[0]+BH2_z*(*u0)[0] < bh2_seg_x_max){
		if((*y0)[0]+BH2_z*(*v0)[0] > -1*bh2_y/2. && (*y0)[0] < bh2_y/2.){
		bh2_pass[i] = true;
		break;
		}
	      }
	    }
	  }
	}
      }
    } //BH2 cut end
    
  
    
    
    for(int i=0;i<NumOfSegBH2;i++){
      if(bh2_pass[i]){
	if(x0->size()>0){
	  hist_bcout_bh2_2d[i]->Fill((*x0)[0]+BH2_z*(*u0)[0],(*y0)[0]+BH2_z*(*v0)[0]);
	  hist_bcout_bac[i]->Fill((*x0)[0]+BAC_z*(*u0)[0]);
	  hist_bcout_bac_2d[i]->Fill((*x0)[0]+BAC_z*(*u0)[0],(*y0)[0]+BAC_z*(*v0)[0]);
	}
      }
    } //Bh2 seg
  }

  c1->Clear();
  TCanvas *c2 = new TCanvas("c2","c2");
  TCanvas *c3 = new TCanvas("c3","c3");
  TCanvas *c4 = new TCanvas("c4","c4");
  c1->Divide(4,4);
  c2->Divide(4,4);
  c3->Divide(4,4);
  c4->Divide(4,4);

  
  double bac_x_cut[NumOfSegBH2][2];
  
  for(int i=0;i<NumOfSegBH2;i++){
    c2->cd(i+1);
    hist_bcout_bh2_2d[i]->Draw("colz");
    
    c4->cd(i+1);
    gPad->SetLogy();
    hist_bcout_bac[i]->Draw();
    double probs[2] = {0.05, 0.95};
    hist_bcout_bac[i]->GetQuantiles(2, bac_x_cut[i], probs);
    DrawLine(bac_x_cut[i][0],hist_bcout_bac[i]);
    DrawLine(bac_x_cut[i][1],hist_bcout_bac[i]);
    c3->cd(i+1);
    hist_bcout_bac_2d[i]->Draw("colz");
  }
  c2->Print(out_pdf);
  c3->Print(out_pdf);
  c4->Print(out_pdf);

  //Make beam file for simulation
  
  TFile *file_beam_old = new TFile("/home/had/haein/Work/E72/Simul/k18geant4/hyptpc-11.0.2/param/BEAM/beam.k.run69_0130.root");
  TTree *tree_beam_old = (TTree*)file_beam_old->Get("tr");
  double pInx,pIny,pInz;
  tree_beam_old->SetBranchAddress("pInx",&pInx);
  tree_beam_old->SetBranchAddress("pIny",&pIny);
  tree_beam_old->SetBranchAddress("pInz",&pInz);

  TFile* file_beam = new TFile(Form("e72_beam_%d.root",runnumber), "RECREATE");
  TTree* tree_beam = new TTree("tree", "tree");
  double x_beam,y_beam,z_beam,u0_beam,v0_beam,p_beam;
  int seg_bh2;
  tree_beam->Branch("x",&x_beam,"x/D");
  tree_beam->Branch("y",&y_beam,"y/D");
  tree_beam->Branch("z",&z_beam,"z/D");
  tree_beam->Branch("u0",&u0_beam,"u0/D");
  tree_beam->Branch("v0",&v0_beam,"v0/D");
  tree_beam->Branch("p",&p_beam,"p/D");
  tree_beam->Branch("seg_bh2",&seg_bh2,"seg_bh2/I");


  
			     
  //again bac cut
  
  for(int n=0;n<data->GetEntries();n++){
    if(n%10000 == 0)cout<<"Entry "<<n<<std::endl;
    data->GetEntry(n);
    bcout->GetEntry(n);
    if(kaon){
      if(btof0 >-4)continue;
    }
    

    //BH2, BAC cut start w/ BcOut
    for(int i=0;i<NumOfSegBH2;i++){
      bh2_pass[i] = false;
      bac_pass[i] = false;
      if((*bh2_adc_u)[i]>bh2_adc_cut[0][i] && (*bh2_adc_d)[i]>bh2_adc_cut[1][i]){
	for(int j=0;j<(*bh2_tdc_s).size();j++){
	  if((*bh2_tdc_s)[i][j]>bh2_tdc_cut[i][0] && (*bh2_tdc_s)[i][j]<bh2_tdc_cut[i][1]){
	    double bh2_seg_x_min = -1*NumOfSegBH2*bh2_x/2.+i*bh2_x;
	    double bh2_seg_x_max = -1*NumOfSegBH2*bh2_x/2.+(i+1)*bh2_x;
	    if((*x0)[0]+BH2_z*(*u0)[0] > bh2_seg_x_min && (*x0)[0]+BH2_z*(*u0)[0] < bh2_seg_x_max){
	      if((*y0)[0]+BH2_z*(*v0)[0] > -1*bh2_y/2. && (*y0)[0] < bh2_y/2.){
		bh2_pass[i] = true;
		break;
	      }
	    }
	  }
	}
      }
      if((*x0)[0]+BAC_z*(*u0)[0] > bac_x_cut[i][0] && (*x0)[0]+BAC_z*((*u0)[0] < bac_x_cut[i][1])){
	bac_pass[i] = true;
      }
    } //BH2 cut end
    
      
    
    //Pure one track cut
    if(ntrack != 1)continue;
    
    for(int i=0;i<NumOfSegBH2;i++){
      if(bh2_pass[i] && bac_pass[i]){
	//Make beam file start
	z_beam = -50.; //mm
	double bcout_z = BAC_z-50;
	tree_beam_old->GetEntry(n_beam);
	n_beam++;
	u0_beam = (*u0)[0];
	v0_beam = (*v0)[0];
	p_beam = TMath::Sqrt(pInx*pInx + pIny*pIny + pInz*pInz)*run_mom/0.906;
	
	if(runnumber<2000)
	  x_beam = (*x0)[0]+bcout_z*(*u0)[0]+7.0;
	else if(runnumber>2000)
	  x_beam = (*x0)[0]+bcout_z*(*u0)[0]+17.25;
	
	y_beam = (*y0)[0]+bcout_z*(*v0)[0];
	seg_bh2 = i;
	
	tree_beam->Fill();
	if(n_beam >= tree_beam_old->GetEntries())n_beam = 0;
	//Make beam file end
	eff_total[i]++;
	hist_bac_npe_s_bh2[i]->Fill((*bac_adc_u)[4] - bac_ped_mean_s);
	if(runnumber<2000){
	  if(i >=4 && i <=10)
	    hist_bac_npe_s_total->Fill((*bac_adc_u)[4] - bac_ped_mean_s);
	}
	else if(runnumber>2000){
	  if(i>=3 && i<=9)
	    hist_bac_npe_s_total->Fill((*bac_adc_u)[4] - bac_ped_mean_s);
	}
	
	//BAC npe cut offline
	//if(bac_npe <npe_threshold)continue;
	for(int j=0;j<(*bac_tdc_u).size();j++){
	  if((*bac_tdc_u)[4][j]>bac_tdc_cut[0] && (*bac_tdc_u)[4][j]<bac_tdc_cut[1]){
	    hist_bac_npe_s_bh2_pass[i]->Fill((*bac_adc_u)[4] - bac_ped_mean_s);
	    if(runnumber<2000){
	      if(i >=4 && i <=10)
		hist_bac_npe_s_pass->Fill((*bac_adc_u)[4] - bac_ped_mean_s);
	    }
	    if(runnumber>2000){
	      if(i>=3 && i<=9)
		hist_bac_npe_s_pass->Fill((*bac_adc_u)[4] - bac_ped_mean_s);
	    }
	    eff_pass[i]++;
	    break;
	  }
	}
      }
    } //Bh2 seg
    
  }
  
  tree_beam->Write();
  file_beam->Close();
  

  //BAC x cut start

  c1->Clear();
  
  c1->Divide(4,4);
  //c2->Divide(3,4);
  //c3->Divide(3,4);
  //c4->Divide(3,4);

  TGraphErrors *g_bac_npe_mean = new TGraphErrors();
  double bac_x_mean[NumOfSegBH2];
  double bac_x_rms[NumOfSegBH2];
  int pid = 0;
  for(int i=0;i<NumOfSegBH2;i++){
    c1->cd(i+1);
    gPad->SetLogy();
    hist_bac_npe_s_bh2[i]->Draw();
    f_bac_npe_s_bh2[i]->SetRange(0,1000);
    hist_bac_npe_s_bh2[i]->Fit(f_bac_npe_s_bh2[i],"RQ");
    hist_bac_npe_s_bh2_pass[i]->SetFillColor(kRed);
    hist_bac_npe_s_bh2_pass[i]->Draw("same");
    const double mean     = f_bac_npe_s_bh2[i]->GetParameter(1);
    const double mean_err = f_bac_npe_s_bh2[i]->GetParError(1);
    
    const double x = (-1*(double)NumOfSegBH2/2.+0.5+i)*bh2_x;
    const double ex = bh2_x/2.;
    
    if(i < 4 || i > 10)continue;
    g_bac_npe_mean->SetPoint(pid, (bac_x_cut[i][0]+bac_x_cut[i][1])/2., mean);
    g_bac_npe_mean->SetPointError(pid, (bac_x_cut[i][1] - bac_x_cut[i][0])/2., mean_err);
    pid++;
  }
  //c2->Print(out_pdf);
  //c3->Print(out_pdf);
  //c4->Print(out_pdf);
  c1->Print(out_pdf);

  //BAC x cut end


  
  
  c1->Clear();
  auto c5 = new TCanvas("c5","c5",1000,500);
  c5->Divide(2);
  auto *g_eff = new TGraphAsymmErrors();
  pid = 0;
  for(int i=0;i<NumOfSegBH2;i++){
    if(i < 4 || i >10)continue;
    TEfficiency tmp("tmp","tmp",1,0,1);
    tmp.SetStatisticOption(TEfficiency::kFCP);

    double eff = (double)eff_pass[i] / (double)eff_total[i];
    tmp.SetTotalEvents(1, eff_total[i]);
    tmp.SetPassedEvents(1, eff_pass[i]);
    double elo = tmp.GetEfficiencyErrorLow(1);
    double ehi = tmp.GetEfficiencyErrorUp(1);
    g_eff->SetPoint(pid, (bac_x_cut[i][0]+bac_x_cut[i][1])/2., eff);
    g_eff->SetPointError(pid, (bac_x_cut[i][1]-bac_x_cut[i][0])/2., (bac_x_cut[i][1]-bac_x_cut[i][0])/2., elo, ehi);
    pid++;
  }
  c5->cd(1);
  c5->SetLeftMargin(0.15);
  c5->SetBottomMargin(0.15);
  g_bac_npe_mean->SetMarkerStyle(20);
  g_bac_npe_mean->GetXaxis()->SetTitle("X [mm]");
  g_bac_npe_mean->GetYaxis()->SetTitle("Np.e.");
  g_bac_npe_mean->GetXaxis()->SetTitleSize(0.1);
  g_bac_npe_mean->GetYaxis()->SetTitleSize(0.1);
  g_bac_npe_mean->Draw("AP");
  c5->cd(2);
  c5->SetLeftMargin(0.15);
  c5->SetBottomMargin(0.15);
  g_eff->SetMarkerStyle(20);
  g_eff->GetXaxis()->SetTitle("X [mm]");
  g_eff->GetXaxis()->SetTitle("Pion Efficiency");
  g_eff->GetXaxis()->SetTitleSize(0.1);
  g_eff->GetYaxis()->SetTitleSize(0.1);
  g_eff->Draw("AP");
  c5->Print(out_pdf + ")");
  //Save graphs
  TFile* f_hist = new TFile(Form("e72_hist_%d.root",runnumber),"RECREATE");
  for(int i=0;i<NumOfSegBH2;i++){
    hist_bac_npe_s_bh2[i]->Write(Form("hist_bac_npe_s%d",i));
    hist_bac_npe_s_bh2_pass[i]->Write(Form("hist_bac_npe_s_bh2_pass%d",i));
  }
  hist_bac_npe_s_total->Write("hist_bac_npe_s_total");
  hist_bac_npe_s_pass->Write("hist_bac_npe_s_pass");
  
  f_hist->Close();


  
  
}


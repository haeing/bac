#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>

bool kaon = false;

int npe_threshold = -15;
const double T0_z = -1100.0;
const double BH2_z = -554.62;
const double BAC_z = -1041.63;

const int NumOfSegT0 = 5;
const int NumOfSegBH2 = 11;
const int NumOfSegBAC = 4;                                                                    

const double t0_tdc_min = 692000.;
const double t0_tdc_max = 697000.;

const double bh2_tdc_min = 705000.;
const double bh2_tdc_max = 720000.;

const double bh2_x = 14.;
const double bh2_y = 100.;

const double bac_tdc_min = 730000.;
const double bac_tdc_max = 740000.;

const double kvc_tdc_min = 710000.;
const double kvc_tdc_max = 725000.;

const double sac_tdc_min = 680000.;
const double sac_tdc_max = 700000.;

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

void analysis_t110(int runnumber, int runnumber_ped)
{
  gROOT->SetBatch(kTRUE);
  double bac_ped_mean[NumOfSegBAC]={0.};
  double t0_adc_cut[2][NumOfSegT0]={0.};
  double t0_tdc_cut[NumOfSegT0][2]={0.};
  double bh2_adc_cut[2][NumOfSegBH2]={0.};
  double bh2_tdc_cut[NumOfSegBH2][2]={0.};
  double bac_gain[NumOfSegBAC]={0.};
  double bac_pxt[NumOfSegBAC]={0.};
  double bac_tdc_cut[2]={725000.,740000.};
  
  //LED Gain

  const char* csv = "led_gain/spe_result_update/spe_summary.csv";
  double target_mppc_hv = 58.0;
  double hv_tol = 1e-6;
  std::ifstream fin(csv);
  if(!fin.is_open()){
    std::cerr << "[ERR] cannot open " << csv << std::endl;
    return;
  }


  std::map<Key, WStat> q1stat;
  std::map<Key, WStat> pxtstat;

  std::string line;

  std::getline(fin,line);
  while (std::getline(fin, line)){
    TString t(line.c_str());
    if (t.Strip(TString::kBoth).IsNull()) continue;

    auto tok = SplitCSV(t);
    //if ((int)tok.size() <= 13)continue;
    
    int board = tok[1].Atoi();
    int ch = tok[2].Atoi();
    double mhv = tok[5].Atoi();

    double q1 = tok[9].Atof();
    double q1e = tok[10].Atof();
    double pxt = tok[11].Atof();
    double pxte = tok[12].Atof();

    if(ch == -1)continue;
    Key k{board, ch};
    q1stat[k].Fill(q1, q1e);
    pxtstat[k].Fill(pxt, pxte);
  }
  fin.close();


  std::vector<int>    v_board, v_ch;
  std::vector<double> v_q1, v_q1err, v_pxt, v_pxterr;

  for (const auto& it : q1stat){
    const Key& k = it.first;
    if (!q1stat[k].Has() || !pxtstat[k].Has()) continue;

    v_board.push_back(k.board);
    v_ch.push_back(k.ch);

    v_q1.push_back(q1stat[k].Mean());
    v_q1err.push_back(q1stat[k].Err());

    v_pxt.push_back(pxtstat[k].Mean());
    v_pxterr.push_back(pxtstat[k].Err());
  }

  for (size_t i=0;i<v_board.size();++i){
    bac_gain[v_board[i]]+=v_q1[i];
    bac_pxt[v_board[i]]+=v_pxt[i];
  }
  for(int i=0;i<NumOfSegBAC;i++){
    bac_gain[i] /= 16.;
    bac_gain[i] *= 1.02; //Temperature Effect(Gain increased)
    bac_pxt[i] /=16.;
    bac_pxt[i] *=(0.45 / 0.65);
    std::cout<<"gain : "<<bac_gain[i]<<", pxt : "<<bac_pxt[i]<<std::endl;
  }
  
  

  
  TString dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025May_root";
  TFile *file = new TFile(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber));
  TTree *data = (TTree*)file->Get("hodo");

  TFile *file_ped = new TFile(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber_ped));
  TTree *data_ped = (TTree*)file_ped->Get("hodo");

  TFile *file_bcout = new TFile(Form("%s/run00%d_BcOutTracking.root",dir.Data(),runnumber));
  TTree *bcout = (TTree*)file_bcout->Get("bcout");

  vector<double>* bac_adc_u_ped = nullptr;
  vector<vector<double>>* bac_tdc_u_ped = nullptr;

  data_ped->SetBranchAddress("bac_adc_u",&bac_adc_u_ped);
  data_ped->SetBranchAddress("bac_tdc_u",&bac_tdc_u_ped);

  double btof0;
  
  vector<double>* t0_adc_u = nullptr;
  vector<double>* t0_adc_d = nullptr;

  vector<vector<double>>* t0_tdc_s = nullptr;
  
  vector<double> *bh2_adc_u = nullptr;
  vector<double> *bh2_adc_d = nullptr;

  vector<vector<double>>* bh2_tdc_s = nullptr;
  
  
  vector<double>* bac_adc_u = nullptr;
  vector<vector<double>>* bac_tdc_u = nullptr;

  vector<double>* sac_adc_u = nullptr;
  vector<vector<double>>* sac_tdc_u = nullptr;

  vector<double>* kvc2_adc_s = nullptr;
  //vector<vector<double>>* kvc_tdc_s = nullptr;

  

  data->SetBranchAddress("btof0",&btof0);
  
  data->SetBranchAddress("t0_adc_u",&t0_adc_u);
  data->SetBranchAddress("t0_adc_d",&t0_adc_d);
  data->SetBranchAddress("t0_tdc_s",&t0_tdc_s);

  data->SetBranchAddress("bh2_adc_u",&bh2_adc_u);
  data->SetBranchAddress("bh2_adc_d",&bh2_adc_d);
  data->SetBranchAddress("bh2_tdc_s",&bh2_tdc_s);

  data->SetBranchAddress("bac_adc_u",&bac_adc_u);
  data->SetBranchAddress("bac_tdc_u",&bac_tdc_u);

  data->SetBranchAddress("sac_adc_u",&sac_adc_u);
  data->SetBranchAddress("sac_tdc_u",&sac_tdc_u);
  
  data->SetBranchAddress("kvc2_adc_s",&kvc2_adc_s);
  //data->SetBranchAddress("kvc2_tdc_s",&kvc2_tdc_s);
  
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
  TH1D *hist_bac_adc_ped[NumOfSegBAC];
  TF1 *f_bac_adc_ped[NumOfSegBAC];
  
  for(int i=0;i<NumOfSegBAC;i++){
    hist_bac_adc_ped[i] = new TH1D(Form("hist_bac_adc_ped%d",i),Form("hist_bac_adc_ped%d",i),200,0,1000);
    f_bac_adc_ped[i] = new TF1(Form("f_bac_adc_ped%d",i),"gaus",100,1000);
  }

  for(int n=0;n<data_ped->GetEntries();n++){
    data_ped->GetEntry(n);
    
    for(int i=0;i<NumOfSegBAC;i++){
      hist_bac_adc_ped[i]->Fill((*bac_adc_u_ped)[i]);
    }
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
  c1->Divide(2,2);
  
  for(int i=0;i<NumOfSegBAC;i++){
    c1->cd(i+1);
    hist_bac_adc_ped[i]->Draw();
    double peak = FindGausPeak(hist_bac_adc_ped[i]);
    f_bac_adc_ped[i]->SetRange(peak-50,peak+50);
    hist_bac_adc_ped[i]->Fit(f_bac_adc_ped[i],"RQ");
    double mean = f_bac_adc_ped[i]->GetParameter(1);
    DrawLine(mean,hist_bac_adc_ped[i]);
    bac_ped_mean[i] = mean;
  }
  c1->Print(out_pdf);
  
  
  TH1D *hist_t0_adc_u[NumOfSegT0];
  TH1D *hist_t0_adc_d[NumOfSegT0];
  TH1D *hist_t0_tdc_s[NumOfSegT0];

  TF1 *f_t0_adc_u[NumOfSegT0];
  TF1 *f_t0_adc_d[NumOfSegT0];
  TF1 *f_t0_tdc_s[NumOfSegT0];
  

  TH1D *hist_bh2_adc_u[NumOfSegBH2];
  TH1D *hist_bh2_adc_d[NumOfSegBH2];
  TH1D *hist_bh2_tdc_s[NumOfSegBH2];

  TF1 *f_bh2_adc_u[NumOfSegBH2];
  TF1 *f_bh2_adc_d[NumOfSegBH2];
  TF1 *f_bh2_tdc_s[NumOfSegBH2];
  
  TH1D *hist_bac_npe[NumOfSegBAC];
  TH1D *hist_bac_npe_s = new TH1D("hist_bac_npe_s","hist_bac_npe_s",100,-10,90);
  TH1D *hist_bac_npe_s_pass = new TH1D("hist_bac_npe_s_pass","hist_bac_npe_s_pass",100,-10,90);
  TH1D *hist_bac_npe_s_kaon = new TH1D("hist_bac_npe_s_kaon","hist_bac_npe_s_kaon",100,-10,90);
  TH1D *hist_bac_npe_s_bh2[NumOfSegBH2];
  TH1D *hist_bac_npe_s_bh2_pass[NumOfSegBH2];
  TH1D *hist_bac_tdc_s = new TH1D("hist_bac_tdc_s","hist_bac_tdc_s",(bac_tdc_max - bac_tdc_min)/tdc_step,bac_tdc_min,bac_tdc_max);
  TH1D *hist_btof = new TH1D("hist_btof","hist_btof",100, -2,7);
  TH1D *hist_btof_pass = new TH1D("hist_btof_pass","hist_btof_pass",100,-2,7);

  TF1 *f_bac_npe_s_bh2[NumOfSegBH2];

  TH2D *hist_bcout_bh2_2d[NumOfSegBH2];
  TH2D *hist_bcout_bac_2d[NumOfSegBH2];
  TH1D *hist_bcout_bac[NumOfSegBH2];

  TH1D *hist_sac_adc = new TH1D("hist_sac_adc","hist_sac_adc",100,0,1000);
  TH1D *hist_sac_tdc = new TH1D("hist_sac_tdc","hist_sac_tdc",100,sac_tdc_min,sac_tdc_max);
  TH2D *hist_bac_sac = new TH2D("hist_bac_sac","hist_bac_sac",100,0,1400,100,0,300);
  
  for(int i=0;i<4;i++){
    //hist_kvc_adc[i] = new TH1D(Form("hist_kvc_adc%d",i),Form("hist_kvc_adc%d",i),
  }
  
  for(int i=0;i<NumOfSegT0;i++){
    hist_t0_adc_u[i] = new TH1D(Form("hist_t0_adc_u%d",i),Form("hist_t0_adc_u%d",i),250,100,350);
    hist_t0_adc_d[i] = new TH1D(Form("hist_t0_adc_d%d",i),Form("hist_t0_adc_d%d",i),250,100,350);
    hist_t0_tdc_s[i] = new TH1D(Form("hist_t0_tdc_s%d",i),Form("hist_t0_tdc_s%d",i),(t0_tdc_max - t0_tdc_min)/tdc_step,t0_tdc_min,t0_tdc_max);
    f_t0_adc_u[i] = new TF1(Form("f_t0_adc_u%d",i),"landau",100,400);
    f_t0_adc_d[i] = new TF1(Form("f_t0_adc_d%d",i),"landau",100,400);
    f_t0_tdc_s[i] = new TF1(Form("f_t0_tdc_s%d",i),"gaus",t0_tdc_min,t0_tdc_max);

  }
  for(int i=0;i<NumOfSegBH2;i++){
    hist_bh2_adc_u[i] = new TH1D(Form("hist_bh2_adc_u%d",i),Form("hist_bh2_adc_u%d",i),900,100,1000);
    hist_bh2_adc_d[i] = new TH1D(Form("hist_bh2_adc_d%d",i),Form("hist_bh2_adc_d%d",i),900,100,1000);
    hist_bh2_tdc_s[i] = new TH1D(Form("hist_bh2_tdc_s%d",i),Form("hist_bh2_tdc_s%d",i),(bh2_tdc_max - bh2_tdc_min)/tdc_step,bh2_tdc_min,bh2_tdc_max);
    hist_bac_npe_s_bh2[i] = new TH1D(Form("hist_bac_npe_s_bh2%d",i),Form("hist_bac_npe_s_bh2%d",i),70,-20,50);
    hist_bac_npe_s_bh2_pass[i] = new TH1D(Form("hist_bac_npe_s_bh2_pass%d",i),Form("hist_bac_npe_s_bh2_pass%d",i),70,-20,50);
    hist_bcout_bac[i] = new TH1D(Form("hist_bcout_bac%d",i),Form("hist_bcout_bac%d",i),100,-150,150);
    hist_bcout_bh2_2d[i] = new TH2D(Form("hist_bcout_bh2_2d%d",i),Form("hist_bcout_bh2_2d%d",i),100,-150,150,100,-150,150);
    hist_bcout_bac_2d[i] = new TH2D(Form("hist_bcout_bac_2d%d",i),Form("hist_bcout_bac_2d%d",i),100,-150,150,100,-150,150);
    
    f_bh2_adc_u[i] = new TF1(Form("f_bh2_adc_u%d",i),"landau",100,1000);
    f_bh2_adc_d[i] = new TF1(Form("f_bh2_adc_d%d",i),"landau",100,1000);
    f_bh2_tdc_s[i] = new TF1(Form("f_bh2_tdc_s%d",i),"gaus",bh2_tdc_min,bh2_tdc_max);
    f_bac_npe_s_bh2[i] = new TF1(Form("f_bac_npe_s_bh2%d",i),"gaus",0,80);
  }

  for(int i=0;i<NumOfSegBAC;i++){
    hist_bac_npe[i] = new TH1D(Form("hist_bac_npe%d",i),Form("hist_bac_npe%d",i),70,-20,50);
  }

  for(int n=0;n<data->GetEntries();n++){
    data->GetEntry(n);
    if(kaon){
      if(btof0 >-4)continue;
      for(int j=0;j<(*sac_tdc_u)[8].size();j++){
	if((*sac_tdc_u)[8][j]>684000 && (*sac_tdc_u)[8][j] <688000)continue;
      }
      if((*sac_adc_u)[8]>120)continue;
    }
    if(n%10000 == 0)cout<<"Entry "<<n<<std::endl;
    for(int i=0;i<NumOfSegT0;i++){
      hist_t0_adc_u[i]->Fill((*t0_adc_u)[i]);
      hist_t0_adc_d[i]->Fill((*t0_adc_d)[i]);
      for(int j=0;j<(*t0_tdc_s)[i].size();j++){
	hist_t0_tdc_s[i]->Fill((*t0_tdc_s)[i][j]);
      }
    }
    for(int i=0;i<NumOfSegBH2;i++){
      hist_bh2_adc_u[i]->Fill((*bh2_adc_u)[i]);
      hist_bh2_adc_d[i]->Fill((*bh2_adc_d)[i]);
      for(int j=0;j<(*bh2_tdc_s)[i].size();j++){
	hist_bh2_tdc_s[i]->Fill((*bh2_tdc_s)[i][j]);
      }
    }
    double bac_npe = 0;
    double bac_adc = 0;
    for(int i=0;i<NumOfSegBAC;i++){
      hist_bac_npe[i]->Fill(((*bac_adc_u)[i]-bac_ped_mean[i])/bac_gain[i]*(1-bac_pxt[i]));
      bac_npe+=((*bac_adc_u)[i]-bac_ped_mean[i])/bac_gain[i]*(1-bac_pxt[i]);
      bac_adc+=((*bac_adc_u)[i]-bac_ped_mean[i]);
    }
    for(int j=0;j<(*bac_tdc_u)[4].size();j++){
      hist_bac_tdc_s->Fill((*bac_tdc_u)[4][j]);
    }
    hist_bac_npe_s->Fill(bac_npe);

    //check sac bac correlation
    
    hist_bac_sac->Fill(bac_adc,(*sac_adc_u)[8]);
    hist_sac_adc->Fill((*sac_adc_u)[8]);
    for(int j=0;j<(*sac_tdc_u)[8].size();j++)
      hist_sac_tdc->Fill((*sac_tdc_u)[8][j]);
    //hist_btof->Fill(btof0*-1);
    for(int j=0;j<(*bac_tdc_u)[4].size();j++){
      if((*bac_tdc_u)[4][j]>bac_tdc_cut[0] && (*bac_tdc_u)[4][j]<bac_tdc_cut[1]){
	hist_bac_npe_s_pass->Fill(bac_npe);
	//hist_btof_pass->Fill(btof0*-1);
      }
    }
    
  }
  
  

  c1->Clear();
  c1->Divide(3,4);
  double mpv = 0;
  double sigma = 0;
  
  for(int i=0;i<NumOfSegT0;i++){
    c1->cd(i+1);
    gPad->SetLogy();
    hist_t0_adc_u[i]->Draw();
    double fit_min = FindLandauRange(hist_t0_adc_u[i]);
    f_t0_adc_u[i]->SetRange(fit_min,fit_min+300);
    hist_t0_adc_u[i]->Fit(f_t0_adc_u[i],"RQ");
    mpv   = f_t0_adc_u[i]->GetParameter(1);
    sigma = f_t0_adc_u[i]->GetParameter(2);
    DrawLine(mpv-sigma*4,hist_t0_adc_u[i]);
    t0_adc_cut[0][i] = mpv-sigma*4;
    
    c1->cd(i+7);
    gPad->SetLogy();
    hist_t0_adc_d[i]->Draw();
    fit_min = FindLandauRange(hist_t0_adc_d[i]);
    f_t0_adc_d[i]->SetRange(fit_min,fit_min+300);
    hist_t0_adc_d[i]->Fit(f_t0_adc_d[i],"RQ");
    mpv   = f_t0_adc_d[i]->GetParameter(1);
    sigma = f_t0_adc_d[i]->GetParameter(2);
    DrawLine(mpv-sigma*4,hist_t0_adc_d[i]);
    t0_adc_cut[1][i] = mpv-sigma*4;
  }
  c1->Print(out_pdf);
  c1->Clear();
  c1->Divide(3,2);
  double mean = 0;
  for(int i=0;i<NumOfSegT0;i++){
    c1->cd(i+1);
    hist_t0_tdc_s[i]->Draw();
    double peak = FindGausPeak(hist_t0_tdc_s[i]);

    f_t0_tdc_s[i]->SetRange(peak-1000,peak+1000);
    hist_t0_tdc_s[i]->Fit(f_t0_tdc_s[i],"RQ");
    mean = f_t0_tdc_s[i]->GetParameter(1);
    sigma = f_t0_tdc_s[i]->GetParameter(2);
    DrawLine(mean-3*sigma,hist_t0_tdc_s[i]);
    DrawLine(mean+3*sigma,hist_t0_tdc_s[i]);
    t0_tdc_cut[i][0] = mean-3*sigma;
    t0_tdc_cut[i][1] = mean+3*sigma;
  }
  c1->Print(out_pdf);
  
  c1->Clear();
  
  c1->Divide(3,4);
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
  
  c1->Divide(3,4);
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
  c1->Divide(3,4);
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
  //hist_bac_npe_s_kaon->Draw();
  //hist_bac_npe_s_pass->SetLineColor(kRed);
  //hist_bac_npe_s_pass->SetFillColor(kRed);
  auto g_bac_npe_s_pass = new TF1("g_bac_npe_s_pass","gaus(0)",20,50);
  hist_bac_npe_s_pass->Fit(g_bac_npe_s_pass,"R");
  //hist_bac_npe_s_pass->Draw();
  c1->cd(6);
  hist_bac_tdc_s->Draw();
  c1->Print(out_pdf);
  c1->Clear();
  auto eff = new TEfficiency(*hist_bac_npe_s_pass, *hist_bac_npe_s);
  TF1* fturn = new TF1("fturn",
		       "[0] + [1]/(1.0 + exp(-(x-[2])/[3]))",
		       -10, 50);
  // [0]: noise floor, [1]: amplitude, [2]: turn-on center, [3]: width
  fturn->SetParameters(0.02, 0.95, 5.0, 1.0);
  auto gr = eff->CreateGraph();
  gr->Fit(fturn, "R");
  c1->Divide(2);
  c1->cd(1);
  eff->Draw("AP");
  c1->cd(2);
  gr->Draw("AP");
  c1->Print(out_pdf);

  //Efficiency, Npe check

  bool t0_pass[NumOfSegT0] = {false};
  bool bh2_pass[NumOfSegBH2] = {false};
  bool bac_pass[NumOfSegBAC] = {false};
  int t0_pass_count = 0;
  int eff_total[NumOfSegBH2] = {0};
  int eff_pass[NumOfSegBH2] = {0};
  for(int n=0;n<data->GetEntries();n++){
    if(n%10000 == 0)cout<<"Entry "<<n<<std::endl;
    data->GetEntry(n);
    bcout->GetEntry(n);
    if(kaon){
      if(btof0 >-4)continue;
      for(int j=0;j<(*sac_tdc_u)[8].size();j++){
	if((*sac_tdc_u)[8][j]>684000 && (*sac_tdc_u)[8][j] <688000)continue;
      }
      if((*sac_adc_u)[8]>120)continue;
      
    }
    //T0 cut start
    t0_pass_count = 0;
    for(int i=0;i<NumOfSegT0;i++){
      t0_pass[i] = false;
      if((*t0_adc_u)[i]>t0_adc_cut[0][i] && (*t0_adc_d)[i]>t0_adc_cut[1][i]){
	for(int j=0;j<(*t0_tdc_s).size();j++){
	  if((*t0_tdc_s)[i][j]>t0_tdc_cut[i][0] && (*t0_tdc_s)[i][j]<t0_tdc_cut[i][1]){
	    t0_pass[i] = true;
	    t0_pass_count++;
	    break;
	  }
	}
      }
    } //T0 cut end

    //BH2 cut start w/ BcOut
    for(int i=0;i<NumOfSegBH2;i++){
      bh2_pass[i] = false;
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
    } //BH2 cut end
    
  
    
    if(t0_pass_count==0)continue;
    
    for(int i=0;i<NumOfSegBH2;i++){
      if(bh2_pass[i]){
	//eff_total[i]++;
	hist_btof->Fill(btof0*-1);
	double bac_npe = 0;
	for(int j=0;j<NumOfSegBAC;j++){
	  bac_npe+=((*bac_adc_u)[j]-bac_ped_mean[j])/bac_gain[j]*(1-bac_pxt[j]);
	  
	}
	//hist_bac_npe_s_bh2[i]->Fill(bac_npe);
	for(int j=0;j<(*bac_tdc_u).size();j++){
	  if((*bac_tdc_u)[4][j]>bac_tdc_cut[0] && (*bac_tdc_u)[4][j]<bac_tdc_cut[1]){
	    //hist_bac_npe_s_bh2_pass[i]->Fill(bac_npe);
	    //eff_pass[i]++;
	    hist_btof_pass->Fill(btof0*-1);
	    break;
	  }
	}
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
  c1->Divide(3,4);
  c2->Divide(3,4);
  c3->Divide(3,4);
  c4->Divide(3,4);

  //TGraphErrors *g_bac_npe_mean = new TGraphErrors();
  //double bac_x_mean[NumOfSegBH2];
  //double bac_x_rms[NumOfSegBH2];
  //int pid = 0;
  double bac_x_cut[NumOfSegBH2][2];
  
  for(int i=0;i<NumOfSegBH2;i++){
    /*
    c1->cd(i+1);
    gPad->SetLogy();
    hist_bac_npe_s_bh2[i]->Draw();
    f_bac_npe_s_bh2[i]->SetRange(10,60);
    hist_bac_npe_s_bh2[i]->Fit(f_bac_npe_s_bh2[i],"RQ");
    hist_bac_npe_s_bh2_pass[i]->SetFillColor(kRed);
    hist_bac_npe_s_bh2_pass[i]->Draw("same");
    const double mean     = f_bac_npe_s_bh2[i]->GetParameter(1);
    const double mean_err = f_bac_npe_s_bh2[i]->GetParError(1);
    
    const double x = (-1*(double)NumOfSegBH2/2.+0.5+i)*bh2_x;
    const double ex = bh2_x/2.;
    */
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
    //bac_x_mean[i] = hist_bcout_bac_2d[i]->GetMean(1);
    //bac_x_rms[i] = hist_bcout_bac_2d[i]->GetRMS(1);
    //if(i <2 || i == NumOfSegBH2-1)continue;
    //g_bac_npe_mean->SetPoint(pid, bac_x_mean[i], mean);
    //g_bac_npe_mean->SetPointError(pid, bac_x_rms[i], mean_err);
    //pid++;
  }
  c2->Print(out_pdf);
  c3->Print(out_pdf);
  c4->Print(out_pdf);
  //c1->Print(out_pdf);



  //Make beam file for simulation
  
  //TFile *file_beam_old = new TFile("/home/had/haein/Work/E72/Simul/k18geant4/hyptpc-11.0.2/param/BEAM/beam.k.run69_0130.root");
  TFile *file_beam_old = new TFile("/hsm/had/sks/E72/JPARC2025Nov/beam_simul/beam.k.run69_0130.root");
  TTree *tree_beam_old = (TTree*)file_beam_old->Get("tr");
  double pInx,pIny,pInz;
  tree_beam_old->SetBranchAddress("pInx",&pInx);
  tree_beam_old->SetBranchAddress("pIny",&pIny);
  tree_beam_old->SetBranchAddress("pInz",&pInz);

  TFile* file_beam = new TFile("t110_beam_735.root", "RECREATE");
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
  int n_old_beam = 0;
  for(int n=0;n<data->GetEntries();n++){
    if(n%10000 == 0)cout<<"Entry "<<n<<std::endl;
    data->GetEntry(n);
    bcout->GetEntry(n);
    if(kaon){
      if(btof0 >-4)continue;
      for(int j=0;j<(*sac_tdc_u)[8].size();j++){
	if((*sac_tdc_u)[8][j]>684000 && (*sac_tdc_u)[8][j] <688000)continue;
      }
      //if((*sac_adc_u)[8]>120)continue;

    }
    
    //T0 cut start
    t0_pass_count = 0;
    for(int i=0;i<NumOfSegT0;i++){
      t0_pass[i] = false;
      if((*t0_adc_u)[i]>t0_adc_cut[0][i] && (*t0_adc_d)[i]>t0_adc_cut[1][i]){
	for(int j=0;j<(*t0_tdc_s).size();j++){
	  if((*t0_tdc_s)[i][j]>t0_tdc_cut[i][0] && (*t0_tdc_s)[i][j]<t0_tdc_cut[i][1]){
	    t0_pass[i] = true;
	    t0_pass_count++;
	    break;
	  }
	}
      }
    } //T0 cut end
    

    //BH2, BAC cut start w/ BcOut
    int seg_bh2 = -9999;
    for(int i=0;i<NumOfSegBH2;i++){
      cout<<"i"<<endl;
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
		seg_bh2 = i;
		break;
	      }
	    }
	  }
	}
      }
      if((*x0)[0]+BAC_z*(*u0)[0] > bac_x_cut[i][0] && (*x0)[0]+BAC_z*((*u0)[0] < bac_x_cut[i][1])){
	bac_pass[i] = true;
      }
      else{seg_bh2 = -9999;}
    } //BH2 cut end

    
    
    
    if(t0_pass_count==0)continue;
    //Pure one track cut
    if(ntrack != 1)continue;

    if(seg_bh2 == -9999)continue;
    
    for(int i=0;i<NumOfSegBH2;i++){
      if(bh2_pass[i] && bac_pass[i]){
	//Make beam file start
	z_beam = -50.; //mm
	double bcout_z = BAC_z-50;
	tree_beam_old->GetEntry(n_old_beam);
	u0_beam = (*u0)[0];
	v0_beam = (*v0)[0];
	p_beam = TMath::Sqrt(pInx*pInx + pIny*pIny + pInz*pInz);
	x_beam = (*x0)[0]+bcout_z*(*u0)[0];
	y_beam = (*y0)[0]+bcout_z*(*v0)[0];
	seg_bh2 = i;
	tree_beam->Fill();
	n_old_beam++;
	if(n_old_beam >= tree_beam_old->GetEntries())n_old_beam = 0;
	//Make beam file end
	eff_total[i]++;
	double bac_npe = 0;
	for(int j=0;j<NumOfSegBAC;j++){
	  bac_npe+=((*bac_adc_u)[j]-bac_ped_mean[j])/bac_gain[j]*(1-bac_pxt[j]);
	}
	hist_bac_npe_s_bh2[i]->Fill(bac_npe);
	hist_bac_npe_s_kaon->Fill(bac_npe);
	
	//BAC npe cut offline
	if(bac_npe <npe_threshold)continue;
	for(int j=0;j<(*bac_tdc_u).size();j++){
	  if((*bac_tdc_u)[4][j]>bac_tdc_cut[0] && (*bac_tdc_u)[4][j]<bac_tdc_cut[1]){
	    hist_bac_npe_s_bh2_pass[i]->Fill(bac_npe);
	    
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
  
  c1->Divide(3,4);
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
    f_bac_npe_s_bh2[i]->SetRange(10,60);
    hist_bac_npe_s_bh2[i]->Fit(f_bac_npe_s_bh2[i],"RQ");
    hist_bac_npe_s_bh2_pass[i]->SetFillColor(kRed);
    hist_bac_npe_s_bh2_pass[i]->Draw("same");
    const double mean     = f_bac_npe_s_bh2[i]->GetParameter(1);
    const double mean_err = f_bac_npe_s_bh2[i]->GetParError(1);
    
    const double x = (-1*(double)NumOfSegBH2/2.+0.5+i)*bh2_x;
    const double ex = bh2_x/2.;
    
    if(i < 3 || i == NumOfSegBH2-1)continue;
    g_bac_npe_mean->SetPoint(pid, (bac_x_cut[i][0]+bac_x_cut[i][1])/2., mean);
    g_bac_npe_mean->SetPointError(pid, (bac_x_cut[i][1] - bac_x_cut[i][0])/2., mean_err);
    pid++;
  }
  //c2->Print(out_pdf);
  //c3->Print(out_pdf);
  //c4->Print(out_pdf);
  c1->Print(out_pdf);


  //Threshold check
  
  //BAC x cut end


  
  
  c1->Clear();
  auto c5 = new TCanvas("c5","c5",1000,500);
  c5->Divide(2);
  auto *g_eff = new TGraphAsymmErrors();
  pid = 0;
  for(int i=0;i<NumOfSegBH2;i++){
    if(i < 3 || i == NumOfSegBH2-1)continue;
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
  c5->Print(out_pdf);
  c5->Clear();
  c5->SetLogy();
  hist_btof->Draw();
  hist_btof_pass->SetFillColor(kRed);
  hist_btof_pass->Draw("same");
  c5->Print(out_pdf);
  c5->Clear();
  c5->SetLogy(false);
  hist_bac_sac->Draw("colz");
  c5->Print(out_pdf);
  c5->Clear();
  c5->Divide(2);
  c5->cd(1);
  hist_sac_adc->Draw("colz");
  c5->cd(2);
  hist_sac_tdc->Draw("colz");
  c5->Print(out_pdf);
  c5->Clear();
  hist_bac_npe_s_kaon->Draw();
  c5->Print(out_pdf + ")");
  //Save graphs
  TFile* f_graph = new TFile(Form("t110_graph_%d.root",runnumber),"RECREATE");
  g_bac_npe_mean->Write("g_bac_npe_mean");
  gr->Write("g_thr");
  g_eff->Write("g_eff");
  hist_btof->Write();
  hist_bac_npe_s_pass->Write();
  f_graph->Close();
  
  
}


#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>

bool kaon = false;
bool pion = false;

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

void analysis_pedestal(int runnumber, int runnumber_ped)
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
  //double bac_tdc_cut[2]={733000.,738000.};
  
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

  data->SetBranchAddress("btof0",&btof0);
  
  data->SetBranchAddress("t0_adc_u",&t0_adc_u);
  data->SetBranchAddress("t0_adc_d",&t0_adc_d);
  data->SetBranchAddress("t0_tdc_s",&t0_tdc_s);

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
  TH1D *hist_bac_npe_s = new TH1D("hist_bac_npe_s","hist_bac_npe_s",60,-20,40);
  TH1D *hist_bac_npe_s_pass = new TH1D("hist_bac_npe_s_pass","hist_bac_npe_s_pass",100,-10,90);
  TH1D *hist_bac_npe_s_bh2[NumOfSegBH2];
  TH1D *hist_bac_npe_s_bh2_pass[NumOfSegBH2];
  TH1D *hist_bac_tdc_s = new TH1D("hist_bac_tdc_s","hist_bac_tdc_s",(bac_tdc_max - bac_tdc_min)/tdc_step,bac_tdc_min,bac_tdc_max);
  TH1D *hist_bac_npe_s_particle = new TH1D("hist_bac_npe_s_particle","hist_bac_npe_s_particle",100,-10,90);
  TH1D *hist_bac_npe_s_particle_pass = new TH1D("hist_bac_npe_s_particle_pass","hist_bac_npe_s_particle_pass",100,-10,90);
  TH2D *hist_bac_btof = new TH2D("hist_bac_btof","hist_bac_btof",130,-10,120,100,-2,7);
  TH1D *hist_btof = new TH1D("hist_btof","hist_btof",100,-2,7);
  TH1D *hist_btof_pass = new TH1D("hist_btof_pass","hist_btof_pass",100,-2,7);
  TH2D *hist_bac_btof_pass = new TH2D("hist_bac_btof_pass","hist_bac_btof_pass",130,-10,120,100,-2,7);
  TF1 *f_bac_npe_s_bh2[NumOfSegBH2];

  TH2D *hist_bcout_bh2_2d[NumOfSegBH2];
  TH2D *hist_bcout_bac_2d[NumOfSegBH2];
  TH1D *hist_bcout_bac[NumOfSegBH2];
  
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
    
    double bac_npe = 0;
    for(int i=0;i<NumOfSegBAC;i++){
      hist_bac_npe[i]->Fill(((*bac_adc_u)[i]-bac_ped_mean[i])/bac_gain[i]*(1-bac_pxt[i]));
      bac_npe+=((*bac_adc_u)[i]-bac_ped_mean[i])/bac_gain[i]*(1-bac_pxt[i]);
    }
    for(int j=0;j<(*bac_tdc_u)[4].size();j++){
      hist_bac_tdc_s->Fill((*bac_tdc_u)[4][j]);
    }
    hist_bac_npe_s->Fill(bac_npe);
    
  }

  TFile* f_graph = new TFile(Form("t110_graph_%d.root",runnumber),"RECREATE");
  hist_bac_npe_s->Write();
  f_graph->Close();
  
}
  
    
  


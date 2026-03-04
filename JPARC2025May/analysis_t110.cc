#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>

const int NumOfSegT0 = 5;
const int NumOfSegBH2 = 12;
const int NumOfSegBAC = 4;

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

void analysis_t110()
{

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

  std::cout << "board  ch   <Q1> ± err     <PXT> ± err\n";
  for (size_t i=0;i<v_board.size();++i){
    std::cout
      << v_board[i] << "  "
      << v_ch[i] << "   "
      << v_q1[i] << " ± " << v_q1err[i] << "   "
      << v_pxt[i] << " ± " << v_pxterr[i] << "\n";
  }

  int runnumber = 303;
  int runnumber_ped = 306;
  TString dir = "/gpfs/group/had/sks/Users/haein/JPARC2025May_root";
  TFile *file = new TFile(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber));
  TTree *data = (TTree*)file->Get("hodo");

  TFile *file_ped = new TFile(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber_ped));
  TTree *data_ped = (TTree*)file_ped->Get("hodo");

  vector<double>* t0_adc_u = nullptr;
  vector<double>* t0_adc_d = nullptr;

  vector<vector<double>>* t0_tdc_s = nullptr;

  vector<double> *bh2_adc_u = nullptr;
  vector<double> *bh2_adc_d = nullptr;

  vector<vector<double>>* bh2_tdc_s = nullptr;
  

  vector<double>* bac_adc_u = nullptr;
  vector<vector<double>>*bac_tdc_u = nullptr;

  data->SetBranchAddress("t0_adc_u",&t0_adc_u);
  data->SetBranchAddress("t0_adc_d",&t0_adc_d);
  data->SetBranchAddress("t0_tdc_s",&t0_tdc_s);

  data->SetBranchAddress("bh2_adc_u",&bh2_adc_u);
  data->SetBranchAddress("bh2_adc_d",&bh2_adc_d);
  data->SetBranchAddress("bh2_tdc_s",&bh2_tdc_s);

  TH1D *hist_t0_adc_u[NumOfSegT0];
  TH1D *hist_t0_adc_d[NumOfSegT0];
  TH1D *hist_t0_tdc_s[NumOfSegT0];

  TH1D *hist_bh2_adc_u[NumOfSegBH2];
  TH1D *hist_bh2_adc_d[NumOfSegBH2];
  TH1D *hist_bh2_tdc_s[NumOfSegBH2];
  
  for(int i=0;i<NumOfSegT0;i++){
    hist_t0_adc_u[i] = new TH1D(Form("hist_t0_adc_u%d",i),Form("hist_t0_adc_u%d",i),400,0,400);
    hist_t0_adc_d[i] = new TH1D(Form("hist_t0_adc_d%d",i),Form("hist_t0_adc_d%d",i),400,0,400);
    hist_t0_tdc_s[i] = new TH1D(Form("hist_t0_tdc_s%d",i),Form("hist_t0_tdc_s%d",i),1000,680000,700000);

  }
  for(int i=0;i<NumOfSegBH2;i++){
    hist_bh2_adc_u[i] = new TH1D(Form("hist_bh2_adc_u%d",i),Form("hist_bh2_adc_u%d",i),400,0,400);
    hist_bh2_adc_d[i] = new TH1D(Form("hist_bh2_adc_d%d",i),Form("hist_bh2_adc_d%d",i),400,0,400);
    hist_bh2_tdc_s[i] = new TH1D(Form("hist_bh2_tdc_s%d",i),Form("hist_bh2_tdc_s%d",i),1000,680000,700000);
  }

  for(int n=0;n<data->GetEntries();n++){
    data->GetEntry(n);
    for(int i=0;i<NumOfSegT0;i++){
      hist_t0_adc_u[i]->Fill((*t0_adc_u)[i]);
      hist_t0_adc_d[i]->Fill((*t0_adc_d)[i]);
    }
    for(int i=0;i<NumOfSegBH2;i++){
    }
  }
  TString out_pdf = "result.pdf";
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Print(out_pdf +"[");

  c1->Divide(3,4);
  for(int i=0;i<NumOfSegT0;i++){
    c1->cd(i+1);
    hist_t0_adc_u[i]->Draw();
    c1->cd(i+7);
    hist_t0_adc_d[i]->Draw();
  }
  c1->Print(out_pdf);
  c1->Print(out_pdf + "]");
  
  
}


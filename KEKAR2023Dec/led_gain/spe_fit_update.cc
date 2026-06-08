#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TString.h"

static int gNmax  = 8;
static int gCompN = 0;

static const bool kek = true;

static inline double Ped_pdf(double x,
                             double m0, double s0,
                             double At, double tau, double x0, double w)
{
  // 1) unnormalized pedestal shape (hOff와 동일한 정의)
  auto f_un = [&](double xx){
    const double gaus = TMath::Gaus(xx, m0, s0, false); // <-- hOff와 동일

    // sigmoid turn-on
    const double S = 1.0 / (1.0 + TMath::Exp(-(xx - m0 - x0)/w));

    // tail: (xx < m0+x0 에서 폭주 방지)
    const double dxp = std::max(0.0, xx - (m0 + x0));
    const double tail = At * S * TMath::Exp(-dxp / tau);

    return gaus + tail;
  };

  // 2) numeric normalization over fixed window (hOff와 동일 window 권장)
  const double xmin = m0 - 100;   // 또는 hOff에서 쓰는 xmin/xmax를 그대로 넣어도 됨
  const double xmax = m0 + 300;
  const int N = 800;
  const double dx = (xmax - xmin)/N;

  double norm = 0.0;
  for(int i=0;i<N;i++){
    const double xx = xmin + (i+0.5)*dx;
    norm += f_un(xx);
  }
  norm *= dx;
  if(norm<=0) norm = 1.0;

  return f_un(x) / norm;  // <-- PDF (적분=1)
}



static inline double g_cluster(int k, double p){
  if(k < 1) return 0.0;
  return std::exp(-p) * std::pow(p, k-1) / TMath::Factorial(k-1);
}

static bool FitPedestalSmooth(TH1* hOff,
                              double& m0, double& s0,
                              double& At, double& tau,
                              double& x0, double& w)
{
  if(!hOff) return false;

  // rough estimates
  const double xPeak = hOff->GetBinCenter(hOff->GetMaximumBin());
  const double rms   = hOff->GetRMS();

  // fit range: pedestal + tail이 충분히 들어가도록
  const double xmin = xPeak - 6.0*rms;
  const double xmax = xPeak + 12.0*rms;

  TF1* fOff = new TF1("fOff_pedSmooth",
    "[0]*(TMath::Gaus(x,[1],[2],false) + [3]*(1.0/(1.0+TMath::Exp(-(x-[1]-[5])/[6])))*TMath::Exp(-(x-[1]-[5])/[4]))",
    xmin, xmax);

  // init
  const double x0_init = 1.5*rms;
  fOff->SetParameters(
    hOff->GetMaximum(),   // [0] Ag
    xPeak,                // [1] m0
    0.6*rms,              // [2] s0
    0.1, // [3] At(scale)
    2.0*rms,              // [4] tau
    x0_init,              // [5] x0
    0.5*rms               // [6] w
  );

  // limits (중요)
  fOff->SetParLimits(2, 0.1, 50.0);                         // s0
  fOff->SetParLimits(3, 0.0, 0.5);       // At (scale)
  fOff->SetParLimits(4, 0.5*rms, 30.0*rms);                 // tau
  fOff->SetParLimits(5, 0.1*rms, 2.0*rms);                  // x0
  fOff->SetParLimits(6, 0.05*rms, 2.0*rms);                 // w

  TFitResultPtr r = hOff->Fit(fOff, "RS"); // Q: quiet

  if(int(r) != 0) {
    // fit 실패 시에도 파라미터를 읽을 수는 있지만, 일단 실패 처리
    return false;
  }

  // output
  m0  = fOff->GetParameter(1);
  s0  = std::fabs(fOff->GetParameter(2));

  // 여기 주의: fOff의 [3]은 "tail amplitude(스케일)"이야.
  // 네 Qmodel의 pedestal PDF 정규화 구조에선 At를 "weight"로 쓰게 될 텐데,
  // 지금 단계에선 일단 fOff 결과를 그대로 At_init로 넣고,
  // LED-on에서는 Fix해서 shape만 유지하는 용도로 쓰면 충분히 잘 동작함.
  At  = fOff->GetParameter(3);
  tau = std::fabs(fOff->GetParameter(4));
  x0  = fOff->GetParameter(5);
  w   = std::fabs(fOff->GetParameter(6));

  return true;
}

// ----- compound Poisson probabilities P(n) up to Nmax
// Primaries N~Poisson(mu0), each primary contributes Y>=1.
// Panjer recursion (compound Poisson):
// P0 = exp(-mu0)
// Pn = (mu0/n) * sum_{k=1..n} k*g_k * P_{n-k}
static void CompoundPoissonPn(std::vector<double>& Pn, int Nmax, double mu0, double pxt){
  Pn.assign(Nmax+1, 0.0);
  if(pxt < 0) pxt = 0;
  if(pxt > 0.999) pxt = 0.999;

  Pn[0] = std::exp(-mu0);
  for(int n=1;n<=Nmax;++n){
    double s = 0.0;
    for(int k=1;k<=n;++k){
      s += k * g_cluster(k, pxt) * Pn[n-k];
    }
    Pn[n] = (mu0 / n) * s;
  }
  // normalize (numerical)
  double sum = 0.0;
  for(int n=0;n<=Nmax;++n) sum += Pn[n];
  if(sum>0) for(int n=0;n<=Nmax;++n) Pn[n] /= sum;
}

// ----- Total model (Gaussian pedestal only)
// par: [0]=A0 [1]=mu0 [2]=m0 [3]=s0 [4]=Q1 [5]=s1 [6]=pxt
static double Qmodel(double *xx, double *par){
  const double x   = xx[0];
  const double A0  = par[0];
  const double mu0 = par[1];
  const double m0  = par[2];
  const double s0  = std::fabs(par[3]);
  const double Q1  = par[4];
  const double s1  = std::fabs(par[5]);
  const double pxt = par[6];
  const double At  = par[7];
  const double tau = std::fabs(par[8]);
  const double x0  = par[9];
  const double w   = std::fabs(par[10]);
  
  std::vector<double> Pn;
  CompoundPoissonPn(Pn, gNmax, mu0, pxt);

  auto G = [](double x, double mean, double sig){
    return (1.0/(std::sqrt(2.0*TMath::Pi())*sig)) *
           std::exp(-0.5*(x-mean)*(x-mean)/(sig*sig));
  };

  double sum = 0.0;

  // n=0 pedestal
  
  //sum += Pn[0] * G(x, m0, s0);
  sum += Pn[0] * Ped_pdf(x, m0, s0, At, tau, x0, w);

  // n>=1
  for(int n=1;n<=gNmax;++n){
    const double var  = s0*s0 + n*s1*s1;
    const double sig  = std::sqrt(std::max(var, 1e-12));
    const double mean = m0 + n*Q1;
    sum += Pn[n] * G(x, mean, sig);
  }

  return A0 * sum;
}

// ----- Component for drawing (single n)
static double Qcomp(double *xx, double *par){
  const double x   = xx[0];
  const double A0  = par[0];
  const double mu0 = par[1];
  const double m0  = par[2];
  const double s0  = std::fabs(par[3]);
  const double Q1  = par[4];
  const double s1  = std::fabs(par[5]);
  const double pxt = par[6];
  const double At  = par[7];
  const double tau = std::fabs(par[8]);
  const double x0  = par[9];
  const double w   = std::fabs(par[10]);
  const int comp = par[11];
  

  std::vector<double> Pn;
  CompoundPoissonPn(Pn, gNmax, mu0, pxt);

  auto G = [](double x, double mean, double sig){
    return (1.0/(std::sqrt(2.0*TMath::Pi())*sig)) *
           std::exp(-0.5*(x-mean)*(x-mean)/(sig*sig));
  };

  //int n = gCompN;
  int n = comp;
  if(n<0 || n>gNmax) return 0.0;

  if(n==0){
    return A0 * Pn[0] * Ped_pdf(x, m0, s0, At, tau, x0, w);
  }else{
    const double var  = s0*s0 + n*s1*s1;
    const double sig  = std::sqrt(std::max(var, 1e-12));
    const double mean = m0 + n*Q1;
    return A0 * Pn[n] * G(x, mean, sig);
  }
}

// ----- Fit pedestal (LED-off) with Gaussian only
static bool FitPedestalGaussian(TH1* hOff, double &m0, double &s0){
  if(!hOff) return false;

  const int bmax = hOff->GetMaximumBin();
  const double xpk = hOff->GetXaxis()->GetBinCenter(bmax);
  const double rms = std::max(hOff->GetRMS(), 1.0);

  TF1 g("g_ped","gaus", xpk-3*rms, xpk+3*rms);
  g.SetParameters(hOff->GetMaximum(), xpk, rms);
  TFitResultPtr r = hOff->Fit(&g, "QSNR");
  if(int(r)!=0) return false;

  m0 = g.GetParameter(1);
  s0 = std::fabs(g.GetParameter(2));
  return true;
}

// ----- Estimate Q1 from peak finding
static bool EstimateQ1(TH1* hOn, double &Q1_init, double &m0_init){
  TSpectrum spec(10);
  int nfound = spec.Search(hOn, 2.0, "nobackground", 0.05);
  
  if(nfound < 1){
    Q1_init = 12.0;
    m0_init = 200;
    return true;
  }
  double* px = spec.GetPositionX();
  std::vector<double> peaks(px, px+nfound);
  std::sort(peaks.begin(), peaks.end());
  if(nfound == 1){
    Q1_init = 12.0;
    m0_init = peaks[0];
    return true;
  }
  m0_init  = peaks[0];
  Q1_init = peaks[1] - peaks[0];
  Q1_init= std::max(Q1_init, 10.0);
  return true;
}


// ----- Fit full spectrum
TF1* FitSPE(TH1* hOn, TH1* hOff, TDirectory* outDir, const char* key){
  double m0=0, s0=0;
  double At=0, tau=0, x0=0, w=0;

  
  if(!FitPedestalSmooth(hOff, m0, s0, At, tau, x0, w)){
    std::cerr << "Pedestal smooth fit failed.\n";
    return nullptr;
  }

  double Q1_init=0;
  EstimateQ1(hOn, Q1_init, m0);
  std::cout<<"q0 init : "<<Q1_init<<", m0 init : "<<m0<<std::endl;
  if(Q1_init<=0) Q1_init = std::max(5.0*s0, 12.0);

  double mu_init = (hOn->GetMean() - hOff->GetMean())/Q1_init;
  if(mu_init < 0.05) mu_init = 0.2;

  int Nguess = (int)std::ceil(mu_init + 10.0*std::sqrt(std::max(mu_init,0.1)));
  if(Nguess < 4)  Nguess = 4;
  if(Nguess > 12) Nguess = 12;
  gNmax = Nguess+1;

  const double xmin = hOn->GetXaxis()->GetXmin();
  const double xmax = hOn->GetXaxis()->GetXmax();

  //TF1* f = new TF1("fSPE_gausPed_xtalk", Qmodel, xmin, xmax, 7);
  TF1* f = new TF1("fSPE_gausPed_xtalk", Qmodel, xmin, xmax, 11);
  f->SetParNames("A0","mu0","m0","sigma0","Q1","sigma1","pxt","At","tau","x0","w");

  const double binw = hOn->GetXaxis()->GetBinWidth(1);
  f->SetParameter(0, hOn->GetEntries()*binw);
  f->SetParLimits(0, hOn->GetEntries()*binw-10000,hOn->GetEntries()*binw+10000);

  f->SetParameter(1, mu_init);
  f->SetParLimits(1, mu_init*0.3,mu_init*7.0);
  //f->SetParLimits(1, mu_init*0.3,mu_init*1.5);
  
  //m0 -=Q1_init*1.;
  f->SetParameter(2, m0);
  f->SetParLimits(2, m0-1,m0+1);
  //f->SetParLimits(2, m0-10,m0+10);

  f->SetParameter(3, s0);
  //f->SetParLimits(3, 0.8*s0, 1.2*s0);
  f->SetParLimits(3, 0.99*s0, 1.01*s0);

  f->SetParameter(4, Q1_init);
  f->SetParLimits(4, 0.95*Q1_init, 1.05*Q1_init);
  //f->SetParLimits(4, 0.7*Q1_init, 1.1*Q1_init);

  f->SetParameter(5, std::max(0.3*Q1_init, s0));
  f->SetParLimits(5, 0.01*s0, 0.5*s0);

  // crosstalk: float near 0.05 (seed), but free
  f->SetParameter(6, 0.1);
  f->SetParLimits(6, 0.0, 0.9);



  f->SetParLimits(7, At*0.99, At*1.01);
  f->SetParLimits(8, tau*0.99,tau*1.01);
  f->SetParLimits(9, x0*0.99,x0*1.01);
  f->SetParLimits(10, w*0.99,w*1.01);



  // fit range
  const double fitL = std::max(xmin, m0 - 10.0*s0);
  const double fitR = std::min(xmax, m0 + (gNmax + 5.0)*Q1_init);
  f->SetRange(fitL, fitR);

  TFitResultPtr rr = hOn->Fit(f, "SNR");
  std::cout << "Fit status = " << int(rr) << "\n";
  std::cout << "m0=" << f->GetParameter(2) << " s0=" << f->GetParameter(3)
            << " Q1=" << f->GetParameter(4) << " s1=" << f->GetParameter(5)
            << " mu0=" << f->GetParameter(1) << " pxt=" << f->GetParameter(6) << "\n";

  if (outDir) {
    outDir->cd();
    rr->Write(key);
    hOn->Write(Form("%s_hist", key));
    f->Write(Form("%s_func", key));
  }
  
  
  return f;

  
}
void spe_fit_update(int board, int channel){
  TFile *fin;
  if(!kek)fin = TFile::Open("BAC_LED_57V.root", "READ");
  if(kek)fin = TFile::Open("BAC_LED_KEK.root", "READ");
  
  TTree *tree = nullptr;
  fin->GetObject("tree", tree);

  
  double ADC_ped[4][3];
  double ADC_led[4][3];
  
  tree->SetBranchAddress("ADC_ped", ADC_ped);
  tree->SetBranchAddress("ADC_led", ADC_led);

  TH1D *hOff = new TH1D("hOff","hOff",500,0,500);
  TH1D *hOn = new TH1D("hOn","hOn",500,0,500);
  
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEntry(n);
    hOff->Fill(ADC_ped[board][channel]);
    hOn->Fill(ADC_led[board][channel]);
  }
  
  if(!hOff || !hOn){
    std::cerr << "Need TH1 named h_ledoff and h_ledon.\n";
    return;
  }

  TFile *fout;
  if(!kek)fout = new TFile(Form("spe_result_update/spe_result%d_%d.root",board,channel), "RECREATE");
  if(kek)fout = new TFile(Form("spe_result_update/spe_result%d_hv%d.root",board,channel+56), "RECREATE");
  
  
  const char* key;
  if(!kek)key = Form("mppc%d_%d", board, channel);
  if(kek)key = Form("mppc%d_hv%d", board, channel+56);
  TF1* f = FitSPE(hOn, hOff, fout, key);

  fout->Close();
  delete fout;
  
  if(!f) return;

  TCanvas* c_ped = new TCanvas("c_ped","Pedestal fitting");
  c_ped->cd(1);
  hOff->SetTitle("LED off (pedestal): Gaussian");
  hOff->Draw();
  hOff->GetFunction("g_ped");
  TCanvas* c = new TCanvas("c_spe","SPE fit: Gaussian pedestal + crosstalk");
  
  c->cd(1);
  //hOn->SetTitle("LED on: total fit + components");
  hOn->SetFillStyle(0);
  hOn->SetLineColor(kBlack);
  hOn->SetLineWidth(2);
  hOn->Draw("hist");
  hOn->SetTitle(";ADC [ch];Counts");

  //int Ndraw = std::min(gNmax, 8);
  int Ndraw = gNmax;
  std::vector<int> cols = {kGray+2,kBlue+1,kGreen+2,kMagenta+1,kOrange+7,kCyan+2,kViolet+1,kTeal+3};

  double fitL, fitR;
  f->GetRange(fitL, fitR);
  
  f->SetLineColor(kRed);
  f->SetLineWidth(3);
  f->SetNpx(800);
  f->Draw("same");

  TLegend* leg = new TLegend(0.60,0.55,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(f,   "Total fit", "l");

  TF1* fn[Ndraw];

  //for(int n=0;n<Ndraw;++n){
  for(int n=0;n<7;++n){
    gCompN = n;
    TString nm;
    nm.Form("fn%d", n);
    //fn[n] = new TF1(nm.Data(), Qcomp, fitL, fitR, 8);
    fn[n] = new TF1(nm.Data(), Qcomp, fitL, fitR, 12);
    for(int ip=0;ip<11;++ip) fn[n]->SetParameter(ip, f->GetParameter(ip));
    fn[n]->SetParameter(11,n);
    fn[n]->SetLineStyle(2);
    fn[n]->SetLineWidth(3);
    fn[n]->SetLineColor(cols[n % cols.size()]);
    fn[n]->SetNpx(800);
    fn[n]->Draw("same");
    TString lab = (n==0) ? "Pedestal" : Form("%d p.e.", n);
    leg->AddEntry(fn[n], lab.Data(), "l");
  }

  
  leg->Draw("same");

  
}

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


// ----- cluster size: Y = 1 + Poisson(pxt)  (k>=1)
static inline double g_cluster(int k, double p){
  if(k < 1) return 0.0;
  return std::exp(-p) * std::pow(p, k-1) / TMath::Factorial(k-1);
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

  std::vector<double> Pn;
  CompoundPoissonPn(Pn, gNmax, mu0, pxt);

  auto G = [](double x, double mean, double sig){
    return (1.0/(std::sqrt(2.0*TMath::Pi())*sig)) *
           std::exp(-0.5*(x-mean)*(x-mean)/(sig*sig));
  };

  double sum = 0.0;

  // n=0 pedestal
  sum += Pn[0] * G(x, m0, s0);

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
  const int comp = par[7];

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
    return A0 * Pn[0] * G(x, m0, s0);
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
static bool EstimateQ1(TH1* hOn, double m0, double &Q1_init){
  TSpectrum spec(10);
  int nfound = spec.Search(hOn, 2.0, "nobackground", 0.08);
  if(nfound<=0){
    Q1_init = std::max(hOn->GetRMS()/5.0, 5.0);
    return true;
  }
  double* px = spec.GetPositionX();
  std::vector<double> peaks(px, px+nfound);
  std::sort(peaks.begin(), peaks.end());
  for(double p: peaks){
    if(p > m0 + 1.0){
      Q1_init = p - m0;
      if(Q1_init>0) return true;
    }
  }
  Q1_init = std::max(hOn->GetRMS()/5.0, 5.0);
  return true;
}

// ----- Fit full spectrum
TF1* FitSPE(TH1* hOn, TH1* hOff){
  double m0=0, s0=0;
  if(!FitPedestalGaussian(hOff, m0, s0)){
    std::cerr << "Pedestal Gaussian fit failed.\n";
    return nullptr;
  }

  double Q1_init=0;
  EstimateQ1(hOn, m0, Q1_init);
  if(Q1_init<=0) Q1_init = std::max(5.0*s0, 5.0);

  double mu_init = (hOn->GetMean() - hOff->GetMean())/Q1_init;
  if(mu_init < 0.05) mu_init = 0.2;

  int Nguess = (int)std::ceil(mu_init + 10.0*std::sqrt(std::max(mu_init,0.1)));
  if(Nguess < 4)  Nguess = 4;
  if(Nguess > 12) Nguess = 12;
  gNmax = Nguess+1;

  const double xmin = hOn->GetXaxis()->GetXmin();
  const double xmax = hOn->GetXaxis()->GetXmax();

  TF1* f = new TF1("fSPE_gausPed_xtalk", Qmodel, xmin, xmax, 7);
  f->SetParNames("A0","mu0","m0","sigma0","Q1","sigma1","pxt");

  const double binw = hOn->GetXaxis()->GetBinWidth(1);
  f->SetParameter(0, hOn->GetEntries()*binw);
  f->SetParLimits(0, hOn->GetEntries()*binw-10000,hOn->GetEntries()*binw+10000);

  f->SetParameter(1, mu_init);
  f->SetParLimits(1, 0.5,1.5);
  //f->SetParLimits(1, mu_init*0.5,mu_init*1.5);

  f->SetParameter(2, m0-20.);
  f->SetParLimits(2, m0-5*s0,m0);

  f->SetParameter(3, s0);
  f->SetParLimits(3, 0.8*s0, 1.2*s0);

  f->SetParameter(4, Q1_init);
  f->SetParLimits(4, 0.8*Q1_init, 1.2*Q1_init);

  f->SetParameter(5, std::max(0.3*Q1_init, s0));
  f->SetParLimits(5, 0.01*s0, 0.5*s0);

  // crosstalk: float near 0.05 (seed), but free
  f->SetParameter(6, 0.1);
  f->SetParLimits(6, 0.0, 0.3);

  f->SetParameter(0,111797);
  f->SetParameter(1,0.877);
  f->SetParameter(2,189.941);
  f->SetParameter(3,4.51);
  f->SetParameter(4,13.4252);
  f->SetParameter(5,2.2381);
  f->SetParameter(6,0.2463);
  

  // fit range
  const double fitL = std::max(xmin, m0 - 10.0*s0);
  const double fitR = std::min(xmax, m0 + (gNmax + 4.0)*Q1_init);
  f->SetRange(fitL, fitR);

  TFitResultPtr rr = hOn->Fit(f, "SNR");
  std::cout << "Fit status = " << int(rr) << "\n";
  std::cout << "m0=" << f->GetParameter(2) << " s0=" << f->GetParameter(3)
            << " Q1=" << f->GetParameter(4) << " s1=" << f->GetParameter(5)
            << " mu0=" << f->GetParameter(1) << " pxt=" << f->GetParameter(6) << "\n";

  return f;
}

// ----- draw components
/*
static void DrawComponents(TH1* hOn, TF1* fTotal, int Ndraw){
  std::vector<int> cols = {kGray+2,kRed+1,kBlue+1,kGreen+2,kMagenta+1,kOrange+7,kCyan+2,kViolet+1,kTeal+3};

  double fitL, fitR;
  fTotal->GetRange(fitL, fitR);
  
  fTotal->SetLineColor(kBlack);
  fTotal->SetLineWidth(3);
  fTotal->Draw("same");

  TLegend* leg = new TLegend(0.60,0.55,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  //leg->AddEntry(hOn, "LED on", "l");
  leg->AddEntry(fTotal,   "Total fit", "l");

  TF1 * fn[Ndraw];
  
  for(int n=0;n<=Ndraw;++n){
    gCompN = n;
    TString nm;
    nm.Form("fn%d", n);
    //TF1* fn = new TF1(nm.Data(), Qcomp, fitL, fitR, 7);
    fn[n] = new TF1(nm.Data(), Qcomp, fitL, fitR, 7);
    for(int ip=0;ip<7;++ip) fn[n]->SetParameter(ip, fTotal->GetParameter(ip));
    fn[n]->SetLineStyle(2);
    fn[n]->SetLineWidth(2);
    //fn[n]->SetLineColor(cols[n % cols.size()]);
    fn[n]->SetNpx(800);
    //fn[n]->Draw("same");
    TString lab = (n==0) ? "n=0 (pedestal)" : Form("n=%d p.e.", n);
    leg->AddEntry(fn[n], lab.Data(), "l");
  }

  
  //leg->AddEntry((TObject*)0, Form("fit p_{xt}=%.3f", f->GetParameter(6)), "");
  leg->Draw("same");
}
*/

// ----- entry point
void fit_gain(){
  TFile *fin = TFile::Open("led_data.root");
  TH1 *hOff = (TH1*)fin->Get("hPedestal_LEDoff");
  TH1 *hOn  = (TH1*)fin->Get("hADC_LEDon");


  if(!hOff || !hOn){
    std::cerr << "Need TH1 named h_ledoff and h_ledon.\n";
    return;
  }

  TF1* f = FitSPE(hOn, hOff);
  if(!f) return;

  TCanvas* c = new TCanvas("c_spe","SPE fit: Gaussian pedestal + crosstalk",1100,450);
  /*
  c->Divide(2,1);

  c->cd(1);
  hOff->SetTitle("LED off (pedestal): Gaussian");
  hOff->Draw("hist");
  hOff->GetFunction("g_ped"); // just to keep
  // (if you want to draw the pedestal fit, repeat it here or store it separately)

  c->cd(2);
  */
  c->cd(1);
  //hOn->SetTitle("LED on: total fit + components");
  hOn->SetFillStyle(0);
  hOn->SetLineColor(kBlack);
  hOn->SetLineWidth(2);
  hOn->Draw("hist");

  int Ndraw = std::min(gNmax, 8);
  //DrawComponents(hOn, f, Ndraw);
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
  //leg->AddEntry(hOn, "LED on", "l");
  leg->AddEntry(f,   "Total fit", "l");

  TF1* fn[Ndraw];

  for(int n=0;n<Ndraw;++n){
    gCompN = n;
    TString nm;
    nm.Form("fn%d", n);
    //TF1* fn = new TF1(nm.Data(), Qcomp, fitL, fitR, 7);
    fn[n] = new TF1(nm.Data(), Qcomp, fitL, fitR, 8);
    for(int ip=0;ip<7;++ip) fn[n]->SetParameter(ip, f->GetParameter(ip));
    fn[n]->SetParameter(7,n);
    fn[n]->SetLineStyle(2);
    fn[n]->SetLineWidth(2);
    fn[n]->SetLineColor(cols[n % cols.size()]);
    fn[n]->SetNpx(800);
    fn[n]->Draw("same");
    TString lab = (n==0) ? "n=0 (pedestal)" : Form("n=%d p.e.", n);
    leg->AddEntry(fn[n], lab.Data(), "l");
  }

  
  //leg->AddEntry((TObject*)0, Form("fit p_{xt}=%.3f", f->GetParameter(6)), "");
  leg->Draw("same");
  //c->Update();


}

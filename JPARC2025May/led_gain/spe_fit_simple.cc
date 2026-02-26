// spe_fit_pois_gausped.C
// Poisson * (sum of Gaussians) model, pedestal is pure Gaussian (LED-off).
// Keeps your hOn/hOff usage and ROOT file saving structure.

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TString.h"

static int gNmax  = 8;

// ---------------------- helpers
static inline double PoisPn(int n, double mu){
  if(n < 0) return 0.0;
  return std::exp(-mu) * std::pow(mu, n) / TMath::Factorial(n);
}

// normalized Gaussian PDF
static inline double Gpdf(double x, double mean, double sig){
  sig = std::fabs(sig);
  if(sig < 1e-12) sig = 1e-12;
  return (1.0/(std::sqrt(2.0*TMath::Pi())*sig)) *
         std::exp(-0.5*(x-mean)*(x-mean)/(sig*sig));
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

// ----- Estimate Q1 from peak finding (LED-on)
static bool EstimateQ1(TH1* hOn, double &Q1_init, double &m0_init){
  if(!hOn) return false;

  TSpectrum spec(10);
  int nfound = spec.Search(hOn, 2.0, "nobackground", 0.05);

  if(nfound < 1){
    Q1_init = 12.0;
    m0_init = hOn->GetBinCenter(hOn->GetMaximumBin());
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
  Q1_init  = peaks[1] - peaks[0];
  Q1_init  = std::max(Q1_init, 10.0);
  return true;
}

// ---------------------- Model: Poisson mixture of Gaussians
// par: [0]=A0 [1]=mu [2]=m0 [3]=s0 [4]=Q1 [5]=s1
static double Qmodel_pois_gausPed(double *xx, double *par){
  const double x   = xx[0];
  const double A0  = par[0];
  const double mu  = par[1];
  const double m0  = par[2];
  const double s0  = std::fabs(par[3]);
  const double Q1  = par[4];
  const double s1  = std::fabs(par[5]);

  double sum = 0.0;

  // n=0 pedestal (Gaussian)
  sum += PoisPn(0, mu) * Gpdf(x, m0, s0);

  // n>=1
  for(int n=1; n<=gNmax; ++n){
    const double var  = s0*s0 + n*s1*s1;
    const double sig  = std::sqrt(std::max(var, 1e-12));
    const double mean = m0 + n*Q1;
    sum += PoisPn(n, mu) * Gpdf(x, mean, sig);
  }

  return A0 * sum;
}

// Component for drawing single n
// par: [0..5] same as above, [6]=comp index
static double Qcomp_pois_gausPed(double *xx, double *par){
  const double x   = xx[0];
  const double A0  = par[0];
  const double mu  = par[1];
  const double m0  = par[2];
  const double s0  = std::fabs(par[3]);
  const double Q1  = par[4];
  const double s1  = std::fabs(par[5]);
  const int comp   = (int)par[6];

  if(comp < 0 || comp > gNmax) return 0.0;

  if(comp == 0){
    return A0 * PoisPn(0, mu) * Gpdf(x, m0, s0);
  }else{
    const int n = comp;
    const double var  = s0*s0 + n*s1*s1;
    const double sig  = std::sqrt(std::max(var, 1e-12));
    const double mean = m0 + n*Q1;
    return A0 * PoisPn(n, mu) * Gpdf(x, mean, sig);
  }
}

// ---------------------- Full fit
TF1* FitSPE_Poisson_GausPed(TH1* hOn, TH1* hOff, TDirectory* outDir, const char* key){
  if(!hOn || !hOff) return nullptr;

  // 1) pedestal from LED-off
  double m0=0, s0=0;
  if(!FitPedestalGaussian(hOff, m0, s0)){
    std::cerr << "Pedestal Gaussian fit failed.\n";
    return nullptr;
  }

  // 2) Q1 seed from LED-on
  double Q1_init=0, m0_on_seed=m0;
  EstimateQ1(hOn, Q1_init, m0_on_seed);

  if(Q1_init <= 0) Q1_init = std::max(5.0*s0, 12.0);

  // 3) mu seed
  double mu_init = (hOn->GetMean() - hOff->GetMean())/Q1_init;
  if(mu_init < 0.05) mu_init = 0.2;

  // 4) choose Nmax from mu
  int Nguess = (int)std::ceil(mu_init + 10.0*std::sqrt(std::max(mu_init,0.1)));
  if(Nguess < 4)  Nguess = 4;
  if(Nguess > 12) Nguess = 12;
  gNmax = Nguess + 1;

  const double xmin = hOn->GetXaxis()->GetXmin();
  const double xmax = hOn->GetXaxis()->GetXmax();

  TF1* f = new TF1("fSPE_pois_gausPed", Qmodel_pois_gausPed, xmin, xmax, 6);
  f->SetParNames("A0","mu","m0","sigma0","Q1","sigma1");

  // A0 ~ entries * binwidth (to convert PDF -> counts)
  const double binw = hOn->GetXaxis()->GetBinWidth(1);
  f->SetParameter(0, hOn->GetEntries()*binw);
  f->SetParLimits(0, hOn->GetEntries()*binw - 10000, hOn->GetEntries()*binw + 10000);

  f->SetParameter(1, mu_init);
  f->SetParLimits(1, mu_init*0.3, mu_init*7.0);

  // Fix/tighten pedestal mean/sigma using LED-off result
  f->SetParameter(2, m0);
  f->SetParLimits(2, m0-1, m0+1);

  f->SetParameter(3, s0);
  f->SetParLimits(3, 0.99*s0, 1.01*s0);

  f->SetParameter(4, Q1_init);
  f->SetParLimits(4, 0.8*Q1_init, 1.2*Q1_init);

  // sigma1: keep it flexible enough (you can tighten later)
  f->SetParameter(5, std::max(0.5*s0, 0.25*Q1_init));
  f->SetParLimits(5, 0.2*s0, 5.0*s0);

  // fit range
  const double fitL = std::max(xmin, m0 - 10.0*s0);
  const double fitR = std::min(xmax, m0 + (gNmax + 5.0)*Q1_init);
  f->SetRange(fitL, fitR);

  TFitResultPtr rr = hOn->Fit(f, "SNR");
  std::cout << "Fit status = " << int(rr) << "\n";
  std::cout << "m0=" << f->GetParameter(2) << " s0=" << f->GetParameter(3)
            << " Q1=" << f->GetParameter(4) << " s1=" << f->GetParameter(5)
            << " mu=" << f->GetParameter(1) << "\n";

  if(outDir) {
    outDir->cd();
    rr->Write(key);
    hOn->Write(Form("%s_hist", key));
    f->Write(Form("%s_func", key));
  }

  return f;
}

// ---------------------- Your wrapper (kept structure)
void spe_fit_simple(int pednumber, int runnumber, int board){
  TString dir = "/gpfs/group/had/sks/Users/haein/JPARC2025May_root";
  TFile *fin_ped = TFile::Open(Form("%s/run00%d_Hodoscope.root",dir.Data(),pednumber),"READ");
  TFile *fin     = TFile::Open(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber),"READ");

  if(!fin_ped || fin_ped->IsZombie() || !fin || fin->IsZombie()){
    std::cerr << "File open failed.\n";
    return;
  }
  std::cout<<"file open"<<std::endl;

  TTree *data_ped = nullptr;
  TTree *data     = nullptr;
  fin_ped->GetObject("hodo", data_ped);
  fin->GetObject("hodo", data);
  if(!data_ped || !data){
    std::cerr << "Tree (hodo) not found.\n";
    return;
  }
  std::cout<<"tree open"<<std::endl;

  std::vector<double> *bac_adc_u_ped = nullptr;
  TBranch *b_bac_adc_u_ped = data_ped->GetBranch("bac_adc_u");
  if(!b_bac_adc_u_ped){
    std::cerr << "Branch bac_adc_u not found in ped tree.\n";
    return;
  }
  b_bac_adc_u_ped->SetAddress(&bac_adc_u_ped);

  std::vector<double> *bac_adc_u = nullptr;
  TBranch *b_bac_adc_u = data->GetBranch("bac_adc_u");
  if(!b_bac_adc_u){
    std::cerr << "Branch bac_adc_u not found in led tree.\n";
    return;
  }
  b_bac_adc_u->SetAddress(&bac_adc_u);

  TH1D *hOff = new TH1D("hOff","hOff",500,0,500);
  TH1D *hOn  = new TH1D("hOn","hOn",500,0,500);

  for(Long64_t n=0; n<data_ped->GetEntries(); n++){
    data_ped->GetEntry(n);
    if(!bac_adc_u_ped) continue;
    if((int)bac_adc_u_ped->size() <= board) continue;
    hOff->Fill((*bac_adc_u_ped)[board]);
  }

  for(Long64_t n=0; n<data->GetEntries(); n++){
    data->GetEntry(n);
    if(!bac_adc_u) continue;
    if((int)bac_adc_u->size() <= board) continue;
    hOn->Fill((*bac_adc_u)[board]);
  }

  if(!hOff || !hOn){
    std::cerr << "Need hOff and hOn.\n";
    return;
  }

  TFile *fout = new TFile(Form("spe_result_simple/spe_result_ped%d_led%d_board%d.root",
                               pednumber,runnumber,board), "RECREATE");
  const char* key = Form("ped%d_led%d_board%d", pednumber,runnumber,board);

  TF1* f = FitSPE_Poisson_GausPed(hOn, hOff, fout, key);

  fout->Close();
  delete fout;

  if(!f) return;

  // ---- Draw pedestal fit (Gaussian)
  TCanvas* c_ped = new TCanvas("c_ped","Pedestal fitting");
  c_ped->cd(1);
  hOff->SetTitle("LED off (pedestal): Gaussian");
  hOff->Draw("hist");
  // g_ped is the name used inside FitPedestalGaussian (TF1 local), so it won't be attached.
  // If you want to draw it, you can refit here, but leaving as-is per request.

  // ---- Draw SPE fit + components
  TCanvas* c = new TCanvas("c_spe","SPE fit: Poisson * Gaussians (Gaussian pedestal)");
  c->cd(1);

  hOn->SetFillStyle(0);
  hOn->SetLineColor(kBlack);
  hOn->SetLineWidth(2);
  hOn->Draw("hist");

  double fitL, fitR;
  f->GetRange(fitL, fitR);

  f->SetLineColor(kRed);
  f->SetLineWidth(3);
  f->SetNpx(800);
  f->Draw("same");

  int Ndraw = std::min(gNmax, 8);
  std::vector<int> cols = {kGray+2,kBlue+1,kGreen+2,kMagenta+1,kOrange+7,kCyan+2,kViolet+1,kTeal+3};

  TLegend* leg = new TLegend(0.60,0.55,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(f, "Total fit", "l");

  TF1* fn[16] = {nullptr};

  for(int n=0; n<Ndraw; ++n){
    TString nm; nm.Form("fn%d", n);

    fn[n] = new TF1(nm.Data(), Qcomp_pois_gausPed, fitL, fitR, 7); // 6 params + comp
    for(int ip=0; ip<6; ++ip) fn[n]->SetParameter(ip, f->GetParameter(ip));
    fn[n]->SetParameter(6, n);

    fn[n]->SetLineStyle(2);
    fn[n]->SetLineWidth(3);
    fn[n]->SetLineColor(cols[n % (int)cols.size()]);
    fn[n]->SetNpx(800);
    fn[n]->Draw("same");

    TString lab = (n==0) ? "Pedestal" : Form("%d p.e.", n);
    leg->AddEntry(fn[n], lab.Data(), "l");
  }

  leg->Draw("same");
}

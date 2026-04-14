// draw_npe_vs_invbeta2.C
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TMath.h"

struct RunInfo {
  int run;
  double mom; // MeV/c
};

double CalcBeta(double p_mev, double mass_mev = 139.57039) {
  // beta = p / sqrt(p^2 + m^2)
  return p_mev / std::sqrt(p_mev * p_mev + mass_mev * mass_mev);
}

void beta_compare(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  // 조건:
  // particle = pi, trig = beam, threshold = 30 인 run만 사용
  std::vector<RunInfo> runs = {
    {1928, 735},
    {2489, 1000},
    {2502, 814},
    {2509, 645},
    {2512, 400},
    {2514, 300}
  };

  std::vector<double> x, y, ex, ey;

  for (const auto& r : runs) {
    TString fname = Form("e72_hist_%d.root", r.run);
    TFile* fin = TFile::Open(fname, "READ");

    if (!fin || fin->IsZombie()) {
      std::cerr << "Cannot open file: " << fname << std::endl;
      continue;
    }

    TH1* h = dynamic_cast<TH1*>(fin->Get("hist_bac_npe_s_pass"));
    if (!h) {
      std::cerr << "Histogram hist_bac_npe_s_pass not found in " << fname << std::endl;
      fin->Close();
      continue;
    }

    // fit을 위해 clone
    TH1* hfit = dynamic_cast<TH1*>(h->Clone(Form("h_run%d", r.run)));
    hfit->SetDirectory(0);
    fin->Close();

    // peak 자동 탐색
    int maxBin = hfit->GetMaximumBin();
    double peakX = hfit->GetBinCenter(maxBin);
    double rms   = hfit->GetRMS();

    // 너무 넓거나 좁지 않게 fit 범위 설정
    double fitMin = peakX - 1.0 * rms;
    double fitMax = peakX + 1.0 * rms;

    // 보호장치
    if (rms <= 0) {
      fitMin = peakX - 5.0;
      fitMax = peakX + 5.0;
    }

    if (fitMin < hfit->GetXaxis()->GetXmin()) fitMin = hfit->GetXaxis()->GetXmin();
    if (fitMax > hfit->GetXaxis()->GetXmax()) fitMax = hfit->GetXaxis()->GetXmax();

    TF1* fgaus = new TF1(Form("fgaus_%d", r.run), "gaus", fitMin, fitMax);
    fgaus->SetParameters(hfit->GetMaximum(), peakX, (rms > 0 ? rms/2.0 : 2.0));

    // 먼저 조용히 fit
    hfit->Fit(fgaus, "RQ0");

    // 첫 fit 결과 기반으로 fit 범위 재조정
    double mean1  = fgaus->GetParameter(1);
    double sigma1 = std::fabs(fgaus->GetParameter(2));

    double refitMin = mean1 - 1.5 * sigma1;
    double refitMax = mean1 + 1.5 * sigma1;

    if (sigma1 > 0 && refitMin < refitMax) {
      fgaus->SetRange(refitMin, refitMax);
      hfit->Fit(fgaus, "RQ0");
    }

    double mean    = fgaus->GetParameter(1);
    double meanErr = fgaus->GetParError(1);

    double beta = CalcBeta(r.mom);
    double invbeta2 = 1.0 / (beta * beta);

    std::cout << "Run " << r.run
              << "  p = " << r.mom << " MeV/c"
              << "  beta = " << beta
              << "  1/beta^2 = " << invbeta2
              << "  mean = " << mean << " +/- " << meanErr
              << std::endl;

    x.push_back(invbeta2);
    y.push_back(mean);
    ex.push_back(0.0);
    ey.push_back(meanErr);
  }

  if (x.empty()) {
    std::cerr << "No valid points found." << std::endl;
    return;
  }

  TGraphErrors* gr = new TGraphErrors(
    (int)x.size(),
    x.data(), y.data(),
    ex.data(), ey.data()
  );

  gr->SetName("gr_npe_mean_vs_invbeta2");
  gr->SetTitle("BAC NPE mean vs 1/#beta^{2};1/#beta^{2};Gaussian mean of hist_bac_npe_s_pass");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.2);
  gr->SetLineWidth(2);

  TCanvas* c1 = new TCanvas("c1", "npe vs 1/beta^2", 800, 600);
  gr->Draw("AP");

  // 원하면 선형 fit도 가능
  TF1* flin = new TF1("flin", "[0] + [1]*x",
                      *std::min_element(x.begin(), x.end()) - 0.02,
                      *std::max_element(x.begin(), x.end()) + 0.02);
  gr->Fit(flin, "R");
  gr->Draw("AP");

  c1->SaveAs("npe_mean_vs_invbeta2.pdf");
  c1->SaveAs("npe_mean_vs_invbeta2.png");

  TFile* fout = new TFile("npe_mean_vs_invbeta2.root", "RECREATE");
  gr->Write();
  if (flin) flin->Write("fit_linear");
  c1->Write();
  fout->Close();
}

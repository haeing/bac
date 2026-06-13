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


const double npe_factor = 37.7331 / 398.085;
const double p_ref = 906.70;
const double sig_p_ref = 11.03;
double m_pi = 139.57039; 
struct RunInfo {
  int run;
  double mom; // MeV/c
  bool bfield;
};

double CalcBeta(double p_mev, double mass_mev = 139.57039) {
  // beta = p / sqrt(p^2 + m^2)
  return p_mev / std::sqrt(p_mev * p_mev + mass_mev * mass_mev);
}

void beta_compare(){
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1111);
  //gStyle->SetOptFit(0);

  // 조건:
  // particle = pi, trig = beam, threshold = 30 인 run만 사용
  std::vector<RunInfo> runs = {
    //HS off

    //{2596, 1000},
    {2585, 933,0},
    //{2502, 814},
    {2587, 755,0},
    {2580, 735,0},
    {2589, 715,0},
    //{2592, 645},
    


    //HS On
    //{2940,842},
    {2815,790,1},
    {2610,715,1},
    {2847,600,1},
    {2884,685,1},
    {3005,665,1},
    {2977,645,1},
    {2512, 400,1},

  };

  std::vector<double> x, y, ex, ey;
  std::vector<double> x1, y1, ex1, ey1;
  std::vector<double> mom,emom;
  std::vector<double> eff,eeff;

  TCanvas *c1 = new TCanvas("c1","c1");
  TString out_pdf = "beta_compare.pdf";
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.04);
  TDatime now;
  TString datetime = Form("%04d-%02d-%02d  %02d:%02d:%02d",
                        now.GetYear(),
                        now.GetMonth(),
                        now.GetDay(),
                        now.GetHour(),
                        now.GetMinute(),
			now.GetSecond());
  text.DrawLatex(0.3,0.5,Form("Date : %s", datetime.Data()));

  c1->Print(out_pdf+"(");
  
  
  
  for (const auto& r : runs) {
    TString fname = Form("e72_hist_%d.root", r.run);
    TFile* fin = TFile::Open(fname, "READ");

    if (!fin || fin->IsZombie()) {
      std::cerr << "Cannot open file: " << fname << std::endl;
      continue;
    }

    //TH1* h = dynamic_cast<TH1*>(fin->Get("hist_bac_npe_s_bh2_pass7"));
    TH1* h;
    if(r.run==2940)h= dynamic_cast<TH1*>(fin->Get("hist_bac_npe_s"));
    else{h= dynamic_cast<TH1*>(fin->Get("hist_bac_npe_s_pass"));}
    TH1* h_total = dynamic_cast<TH1*>(fin->Get("hist_bac_npe_s_total"));

    double Ntotal = h_total->GetEntries();
    double Npass  = h->GetEntries();
      
    
    //TH1* h = dynamic_cast<TH1*>(fin->Get("hist_bac_npe_s7"));
    if (!h) {
      std::cerr << "Histogram hist_bac_npe_s_pass not found in " << fname << std::endl;
      fin->Close();
      continue;
    }

    TH1* hfit = dynamic_cast<TH1*>(h->Clone(Form("h_run%d", r.run)));
    hfit->SetDirectory(0);
    fin->Close();


    int maxBin = hfit->GetMaximumBin();
    double peakX = hfit->GetBinCenter(maxBin);
    double rms   = hfit->GetRMS();

    double fitMin = peakX - 10.0 * rms;
    double fitMax = peakX + 10.0 * rms;

    if (rms <= 0) {
      fitMin = peakX - 5.0;
      fitMax = peakX + 5.0;
    }

    if (fitMin < hfit->GetXaxis()->GetXmin()) fitMin = hfit->GetXaxis()->GetXmin();
    if (fitMax > hfit->GetXaxis()->GetXmax()) fitMax = hfit->GetXaxis()->GetXmax();

    if(r.run == 1928){
      fitMin = 300;
      fitMax = 550;
    }
    else if(r.run == 2489){
      fitMin = 300;
      fitMax = 600;
    }
    else if(r.run == 2502){
      fitMin = 300;
      fitMax = 550;
    }
    else if(r.run == 2509){
      fitMin = 350;
      fitMax = 600;
    }
    else if(r.run ==2512){
      fitMin = 100;
      fitMax = 350;
    }
    else if(r.run ==2596){
      fitMin = 200;
      fitMax = 700;
    }
    else if(r.run ==2592){
      fitMin = 300;
      fitMax = 700;
    }
    
    
    TF1 * fgaus;
    fgaus = new TF1(Form("fgaus_%d", r.run), "gaus", fitMin, fitMax);
    //fgaus->SetParameters(hfit->GetMaximum(), peakX, (rms > 0 ? rms/2.0 : 2.0));
    
    hfit->Fit(fgaus, "RQ");
    c1->Clear();
    hfit->Draw();
    c1->Print(out_pdf);

    double mean1  = fgaus->GetParameter(1);
    double sigma1 = std::fabs(fgaus->GetParameter(2));
    
    double refitMin = mean1 - 1.5 * sigma1;
    double refitMax = mean1 + 1.5 * sigma1;

    /*
    if (sigma1 > 0 && refitMin < refitMax) {
      fgaus->SetRange(refitMin, refitMax);
      hfit->Fit(fgaus, "RQ0");
    }
    */
    
    double mean    = fgaus->GetParameter(1);
    double meanErr = fgaus->GetParError(1);
    double sig_npe = fgaus->GetParameter(2);

    double beta = CalcBeta(r.mom);
    double invbeta2 = 1.0 / (beta * beta);

    double sig_p = r.mom * sig_p_ref / p_ref;
    double sig_invbeta2 = 2.0 * m_pi*m_pi / (r.mom*r.mom*r.mom) * sig_p;

    std::cout << "Run " << r.run
              << "  p = " << r.mom << " MeV/c"
              << "  beta = " << beta
              << "  1/beta^2 = " << invbeta2
              << "  mean = " << mean << " +/- " << meanErr
              << std::endl;
    if(r.bfield==0){
      x.push_back(invbeta2);
      y.push_back(mean*npe_factor);
      ex.push_back(sig_invbeta2);
      ey.push_back(meanErr*npe_factor);
      mom.push_back(r.mom);
      emom.push_back(sig_p);
      double eff_p = Npass / Ntotal;
      eff.push_back(eff_p);
      eeff.push_back(sqrt(eff_p*(1.0-eff_p)/Ntotal));
    }
    else if(r.bfield==1){
      x1.push_back(invbeta2);
      y1.push_back(mean*npe_factor);
      ex1.push_back(sig_invbeta2);
      ey1.push_back(meanErr*npe_factor);
    }
  }


  //Add 600 MeV/c point in T110
  /*
  double beta = CalcBeta(600);
  double invbeta2 = 1.0 / (beta * beta);

  double sig_p = 600 * sig_p_ref / p_ref;
  double sig_invbeta2 = 2.0 * m_pi*m_pi / (600*600*600) * sig_p;
    
  x.push_back(invbeta2);
  y.push_back(35.8075);
  ex.push_back(sig_invbeta2);
  ey.push_back(0.096);

  TFile* fin = TFile::Open("../JPARC2025May/t110_graph_351.root", "READ");
  TH1* h= dynamic_cast<TH1*>(fin->Get("hist_bac_npe_s_particle_pass"));
  TH1* h_total = dynamic_cast<TH1*>(fin->Get("hist_bac_npe_s_particle"));
  double Ntotal = h_total->GetEntries();
  double Npass  = h->GetEntries();
  mom.push_back(600.);
  emom.push_back(sig_p);
  double eff_p = Npass / Ntotal;
  eff.push_back(eff_p);
  eeff.push_back(sqrt(eff_p*(1.0-eff_p)/Ntotal));
  */

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
  gr->SetTitle(";1/#beta^{2};Np.e.");
  gr->GetXaxis()->SetLimits(1.01, 1.15);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);
  gr->SetMarkerSize(1.2);
  gr->SetLineWidth(2);

  TGraphErrors* gr1 = new TGraphErrors(
    (int)x1.size(),
    x1.data(), y1.data(),
    ex1.data(), ey1.data()
  );

  //gr->SetName("gr_npe_mean_vs_invbeta2");
  //gr->SetTitle(";1/#beta^{2};Np.e.");
  //gr->GetXaxis()->SetLimits(1.01, 1.15);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(kBlue);
  gr1->SetLineColor(kBlue);
  gr1->SetMarkerSize(1.2);
  gr1->SetLineWidth(2);

  auto mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr1);
  mg->SetName("gr_npe_mean_vs_invbeta2");
  mg->SetTitle(";1/#beta^{2};Np.e.");
  
   

  
  //TCanvas* c1 = new TCanvas("c1", "npe vs 1/beta^2", 800, 600);
  //gr->Draw("AP");

  // 원하면 선형 fit도 가능
  TF1* flin = new TF1("flin", "[0] + [1]*x",
                      *std::min_element(x.begin(), x.end()) - 0.02,
                      //*std::max_element(x.begin(), x.end()) + 0.02);
		      1.30);
  //gr->Fit(flin, "R","",1.00,1.15);
  flin->SetParameters(207.,-164.);
  gr->Fit(flin, "R0","",1.00,1.30);
  //gr->Draw("AP");
  
  //c1->SaveAs("npe_mean_vs_invbeta2.pdf");
  //c1->SaveAs("npe_mean_vs_invbeta2.png");

  c1->Clear();
  gr->GetXaxis()->SetLimits(1.00,1.30);
  gr->GetYaxis()->SetRangeUser(0, 42);
  //gr->Draw("AP");

  mg->GetXaxis()->SetLimits(1.00,1.30);
  mg->GetYaxis()->SetRangeUser(0, 42);
  //mg->Draw("AP");
  TBox *box = new TBox(1.029, 0, 1.047, 42);
  box->SetFillColorAlpha(kRed, 0.15); 
  box->SetLineColor(0);               

  //box->Draw("same");
  //gr->Draw("P same");
  //  mg->Draw("P same");
  

  TMultiGraph *gr_zoom =(TMultiGraph*)mg->Clone("gr_zoom");
  //gr_zoom->Draw("AP");
  gr_zoom->GetXaxis()->SetLimits(1.03, 1.047);
  gr_zoom->GetYaxis()->SetRangeUser(35, 40);
  

  
  TF1* fit_zoom = new TF1("fit_zoom", "[0] + [1]*x",
                      *std::min_element(x.begin(), x.end()) - 0.02,
                      *std::max_element(x.begin(), x.end()) + 0.02);
  fit_zoom->SetParameter(0, flin->GetParameter(0));
  fit_zoom->SetParameter(1, flin->GetParameter(1));
  fit_zoom->SetLineColor(kRed);
  fit_zoom->SetLineWidth(2);
  fit_zoom->SetNpx(1000);
  gr_zoom->Draw("AP");      // 먼저 graph
  fit_zoom->Draw("same");
  
  //gr_zoom->SetMarkerSize(1.2);
  
  /*
  gr_zoom->GetXaxis()->SetNdivisions(505); // major 5개 정도
  gr_zoom->GetYaxis()->SetNdivisions(505); // major 5개 정도
  gr_zoom->GetXaxis()->SetTickLength(0.06);
  gr_zoom->GetYaxis()->SetTickLength(0.06);
  gr_zoom->GetXaxis()->SetLabelSize(0.08);
  gr_zoom->GetYaxis()->SetLabelSize(0.08);
  gr_zoom->GetXaxis()->SetTitleSize(0.08);
  gr_zoom->GetYaxis()->SetTitleSize(0.08);
  gr_zoom->GetXaxis()->SetTitleOffset(0.9);
  gr_zoom->GetXaxis()->SetLabelOffset(0.02);
  gr_zoom->GetYaxis()->SetTitleOffset(0.75);
  */
  //fit_zoom->SetRange(1.03, 1.047);
  
  //fit_zoom->Draw("same");
  //fit_zoom->Draw();
  //gr_zoom->Draw("P same");
  
  TLegend *leg = new TLegend(0.15, 0.15, 0.35, 0.30);
  leg->AddEntry(gr, "B = 0 T");
  leg->AddEntry(gr1, "B = 0.8 T (SHS field)");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("same");

  
  
  //TPad *pad = new TPad("pad","pad",0.5,0.5,0.9,0.9);
  /*
  TPad *pad = new TPad("pad","pad",0.1,0.1,0.5,0.5);
  pad->SetFillStyle(4000); // transparent
  pad->SetFrameFillStyle(0);
  pad->SetLeftMargin(0.12);
  pad->SetRightMargin(0.04);
  pad->SetBottomMargin(0.25);
  pad->SetTopMargin(0.04);
  pad->Draw();
  pad->cd();
  gr->SetMarkerSize(1.2);
  gr1->SetMarkerSize(1.2);
  mg->Draw("AP");
  box->Draw("same");
  gPad->Modified();
  gPad->Update();
  */
  //TGraphErrors *gr_zoom =(TGraphErrors*)gr->Clone("gr_zoom");
  
  TCanvas *c2 = new TCanvas("c2","c2");
  
  mg->Draw("AP");
  //gr->SetLineColor(0);
  //gr1->SetLineColor(0);
  
  box->Draw("same");
  
  flin->Draw("same");
  // after drawing graph
  gPad->Update();

  //gPad->Update();
  
  
  

  

  c1->cd();

 
  c1->Print(out_pdf+")");
  c1->SaveAs("beta_npe.pdf");

  TCanvas *c3 = new TCanvas("c3","c3");
  TGraphErrors* gr_eff = new TGraphErrors(
    (int)mom.size(),
    mom.data(), eff.data(),
    emom.data(), eeff.data()
  );
  gr_eff->SetMarkerStyle(20);
  gr_eff->Draw("AP");
  /*
    TFile* fout = new TFile("npe_mean_vs_invbeta2.root", "RECREATE");
    gr->Write();
    if (flin) flin->Write("fit_linear");
    c1->Write();
    fout->Close();
  */
}

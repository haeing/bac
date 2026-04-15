#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

void transmittance()
{
  // =========================
  // Input
  // =========================
  std::string fname = "ag22k001.txt";
  double t_mm = 12.4;        // sample thickness [mm]
  double lam_eval_nm = 400.; // wavelength for A_T

  // thickness for transmission length
  double d_cm = t_mm / 10.0; // 10 mm = 1 cm

  // =========================
  // Load data
  // =========================
  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cerr << "Error: cannot open file " << fname << std::endl;
    return;
  }

  std::vector<double> lam_nm;
  std::vector<double> T;

  double lam, trans_percent;
  while (fin >> lam >> trans_percent) {
    double trans = trans_percent / 100.0; // convert % -> fraction
    if (trans > 0.0 && trans < 1.0) {
      lam_nm.push_back(lam);
      T.push_back(trans);
    }
  }
  fin.close();

  if (lam_nm.empty()) {
    std::cerr << "Error: no valid data loaded from " << fname << std::endl;
    return;
  }

  // =========================
  // TGraph
  // =========================
  TGraph* g = new TGraph(lam_nm.size());
  for (size_t i = 0; i < lam_nm.size(); ++i) {
    g->SetPoint(i, lam_nm[i], T[i]);
  }

  // =========================
  // Rayleigh model
  // T = A * exp(-B / lambda^4)
  // =========================
  double xmin = *std::min_element(lam_nm.begin(), lam_nm.end())-10;
  double xmax = *std::max_element(lam_nm.begin(), lam_nm.end());

  TF1* f = new TF1("f", "[0]*exp(-[1]/pow(x,4))", xmin, xmax);
  f->SetParName(0, "A");
  f->SetParName(1, "B");
  f->SetParameters(0.95, 1e10);

  TFitResultPtr fitres = g->Fit(f, "S R");

  double A_fit = f->GetParameter(0);
  double B_fit = f->GetParameter(1);
  double A_err = f->GetParError(0);
  double B_err = f->GetParError(1);

  // =========================
  // Rayleigh coefficient C
  // =========================
  double C_nm4_per_mm = B_fit / t_mm;
  double C_nm4_per_mm_err = B_err / t_mm;

  // unit conversion: nm^4/mm -> um^4/cm
  double C_um4_per_cm = C_nm4_per_mm * 1e-11;
  double C_um4_per_cm_err = C_nm4_per_mm_err * 1e-11;

  // =========================
  // Transmittance at 400 nm
  // =========================
  double T400 = f->Eval(lam_eval_nm);

  // =========================
  // Transmission length A_T(400)
  // A_T = -d / ln T
  // =========================
  double AT400_cm = std::numeric_limits<double>::quiet_NaN();
  if (T400 > 0.0 && T400 < 1.0) {
    AT400_cm = -d_cm / std::log(T400);
  }

  // =========================
  // Print results
  // =========================
  std::cout << "==== Rayleigh scattering fit (paper-ready) ====" << std::endl;
  std::cout << "A = " << A_fit << " ± " << A_err << std::endl;
  std::cout << "C = (" << C_um4_per_cm << " ± " << C_um4_per_cm_err
            << ") um^4 / cm" << std::endl;
  std::cout << std::endl;
  std::cout << "T(" << lam_eval_nm << " nm) = " << T400
            << "  (" << T400 * 100.0 << " %)" << std::endl;
  std::cout << "A_T(" << lam_eval_nm << " nm) = " << AT400_cm << " cm" << std::endl;
  std::cout << "(using thickness d = " << d_cm << " cm = " << t_mm << " mm)" << std::endl;

  // =========================
  // Style
  // =========================
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetTextFont(132);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  gStyle->SetLegendFont(132);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // =========================
  // Plot
  // =========================
  TCanvas* c = new TCanvas("c", "Rayleigh fit", 800, 500);
  c->SetMargin(0.12, 0.05, 0.12, 0.05);

  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);

  g->Draw("AP");
  g->GetXaxis()->SetTitle("Wavelength [nm]");
  g->GetYaxis()->SetTitle("Transmittance");
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetYaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetLabelSize(0.045);
  g->GetYaxis()->SetLabelSize(0.045);

  g->GetYaxis()->SetRangeUser(0.0, 1.05);

  f->SetLineColor(kRed + 1);
  f->SetLineWidth(3);
  f->Draw("SAME");

  c->SaveAs("ag22k001_rayleigh_fit.pdf");
  //c->SaveAs("ag22k001_rayleigh_fit_root.png");
}

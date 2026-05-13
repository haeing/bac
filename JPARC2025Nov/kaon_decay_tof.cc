#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TMath.h"

void kaon_decay_tof()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // =========================
  // Constants
  // =========================
  const double c_m_per_ns = 0.299792458; // m/ns

  const double L = 7.0;          // total flight path [m]
  const double pK_lab = 735.0;   // kaon momentum [MeV/c]

  const double mK   = 493.677;   // MeV/c2
  const double mpiC = 139.57039; // charged pion
  const double mpi0 = 134.9768;  // neutral pion
  const double mmu  = 105.65837; // muon
  const double mnu  = 0.0;       // neutrino

  const int NevtPerStep = 30000;

  // =========================
  // Parent kaon
  // =========================
  const double EK = std::sqrt(pK_lab*pK_lab + mK*mK);
  TLorentzVector K_lab(0.0, 0.0, pK_lab, EK);

  const double betaK = pK_lab / EK;

  // =========================
  // Reference TOF: 735 MeV/c pion over 7 m
  // =========================
  const double Epi_ref = std::sqrt(pK_lab*pK_lab + mpiC*mpiC);
  const double betaPi_ref = pK_lab / Epi_ref;
  const double tofPiRef = L / (betaPi_ref * c_m_per_ns);

  std::cout << "K momentum        = " << pK_lab << " MeV/c" << std::endl;
  std::cout << "beta_K            = " << betaK << std::endl;
  std::cout << "beta_pi(735)      = " << betaPi_ref << std::endl;
  std::cout << "TOF pi(735), 7 m  = " << tofPiRef << " ns" << std::endl;

  // =========================
  // Graphs per decay position
  // =========================
  std::vector<TGraph*> gPi;
  std::vector<TGraph*> gMu;

  std::vector<int> colors = {
    kBlack, kRed+1, kBlue+1, kGreen+2,
    kMagenta+1, kOrange+7, kCyan+2, kViolet+1
  };

  // =========================
  // Loop over decay position: 0,1,...,7 m
  // =========================
  for (int istep = 0; istep <= 7; ++istep) {

    const double xdecay = (double)istep; // m
    const double Lrem   = L - xdecay;    // m

    const double tofK_before_decay =
      xdecay / (betaK * c_m_per_ns); // ns

    TGraph* g_pi_step = new TGraph();
    TGraph* g_mu_step = new TGraph();

    g_pi_step->SetName(Form("g_pi_xdecay_%d", istep));
    g_mu_step->SetName(Form("g_mu_xdecay_%d", istep));

    int col = colors[istep % colors.size()];

    g_pi_step->SetMarkerStyle(20);
    g_pi_step->SetMarkerSize(0.35);
    g_pi_step->SetMarkerColor(col);
    g_pi_step->SetLineColor(col);

    g_mu_step->SetMarkerStyle(20);
    g_mu_step->SetMarkerSize(0.35);
    g_mu_step->SetMarkerColor(col);
    g_mu_step->SetLineColor(col);

    // -------------------------------------------------
    // K- -> pi- pi0
    // -------------------------------------------------
    {
      TGenPhaseSpace event;
      double masses[2] = {mpiC, mpi0};
      event.SetDecay(K_lab, 2, masses);

      for (int iev = 0; iev < NevtPerStep; ++iev) {
        event.Generate();

        TLorentzVector* pi_lab = event.GetDecay(0); // charged pi-

        double p_pi = pi_lab->P(); // MeV/c
        double E_pi = pi_lab->E();
        double beta_pi = p_pi / E_pi;

        double tof =
          tofK_before_decay
          + Lrem / (beta_pi * c_m_per_ns);

        double dtof = tof - tofPiRef;

        int n = g_pi_step->GetN();
        g_pi_step->SetPoint(n, p_pi, dtof);
      }
    }

    // -------------------------------------------------
    // K- -> mu- anti-nu_mu
    // -------------------------------------------------
    {
      TGenPhaseSpace event;
      double masses[2] = {mmu, mnu};
      event.SetDecay(K_lab, 2, masses);

      for (int iev = 0; iev < NevtPerStep; ++iev) {
        event.Generate();

        TLorentzVector* mu_lab = event.GetDecay(0); // mu-

        double p_mu = mu_lab->P(); // MeV/c
        double E_mu = mu_lab->E();
        double beta_mu = p_mu / E_mu;

        double tof =
          tofK_before_decay
          + Lrem / (beta_mu * c_m_per_ns);

        double dtof = tof - tofPiRef;

        int n = g_mu_step->GetN();
        g_mu_step->SetPoint(n, p_mu, dtof);
      }
    }

    gPi.push_back(g_pi_step);
    gMu.push_back(g_mu_step);
  }

  // =========================
  // Plot 1: K -> pi pi0
  // =========================
  TCanvas* c1 = new TCanvas("c1", "Delta TOF vs pion momentum", 1000, 700);
  c1->SetMargin(0.13, 0.05, 0.12, 0.05);
  c1->SetGrid();

  TH2D* hframe1 = new TH2D("hframe1", "",
                           100, 0, 1000,
                           100, -5, 15);
  hframe1->GetXaxis()->SetTitle("decay #pi^{-} momentum [MeV/c]");
  hframe1->GetYaxis()->SetTitle("#DeltaTOF from 735 MeV/c #pi over 7 m [ns]");
  hframe1->GetXaxis()->CenterTitle();
  hframe1->GetYaxis()->CenterTitle();
  hframe1->Draw();

  for (auto g : gPi) g->Draw("P SAME");

  TLine* zero1 = new TLine(0, 0, 1000, 0);
  zero1->SetLineColor(kGray+2);
  zero1->SetLineStyle(2);
  zero1->SetLineWidth(2);
  zero1->Draw("SAME");

  TLegend* leg1 = new TLegend(0.62, 0.55, 0.88, 0.88);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  for (int i = 0; i <= 7; ++i) {
    leg1->AddEntry(gPi[i], Form("x_{decay} = %d m", i), "p");
  }
  leg1->Draw();

  TLatex lat;
  lat.SetTextFont(132);
  lat.SetTextSize(0.045);
  lat.DrawLatex(60, 12.5, "K^{-} #rightarrow #pi^{-} #pi^{0}");

  //c1->SaveAs("dtof_vs_pion_momentum_Kpipi.pdf");
  //c1->SaveAs("dtof_vs_pion_momentum_Kpipi.png");

  // =========================
  // Plot 2: K -> mu nu
  // =========================
  TCanvas* c2 = new TCanvas("c2", "Delta TOF vs muon momentum", 1000, 700);
  c2->SetMargin(0.13, 0.05, 0.12, 0.05);
  c2->SetGrid();

  TH2D* hframe2 = new TH2D("hframe2", "",
                           100, 0, 1100,
                           100, -5, 15);
  hframe2->GetXaxis()->SetTitle("decay #mu^{-} momentum [MeV/c]");
  hframe2->GetYaxis()->SetTitle("#DeltaTOF from 735 MeV/c #pi over 7 m [ns]");
  hframe2->GetXaxis()->CenterTitle();
  hframe2->GetYaxis()->CenterTitle();
  hframe2->Draw();

  for (auto g : gMu) g->Draw("P SAME");

  TLine* zero2 = new TLine(0, 0, 1100, 0);
  zero2->SetLineColor(kGray+2);
  zero2->SetLineStyle(2);
  zero2->SetLineWidth(2);
  zero2->Draw("SAME");

  TLegend* leg2 = new TLegend(0.62, 0.55, 0.88, 0.88);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  for (int i = 0; i <= 7; ++i) {
    leg2->AddEntry(gMu[i], Form("x_{decay} = %d m", i), "p");
  }
  leg2->Draw();

  lat.DrawLatex(60, 12.5, "K^{-} #rightarrow #mu^{-} #bar{#nu}_{#mu}");

  //c2->SaveAs("dtof_vs_muon_momentum_Kmunu.pdf");
  //c2->SaveAs("dtof_vs_muon_momentum_Kmunu.png");

  // =========================
  // Combined canvas
  // =========================
  TCanvas* c3 = new TCanvas("c3", "Delta TOF combined", 1600, 650);
  c3->Divide(2,1);

  c3->cd(1);
  gPad->SetMargin(0.13, 0.05, 0.12, 0.05);
  gPad->SetGrid();
  hframe1->Draw();
  for (auto g : gPi) g->Draw("P SAME");
  zero1->Draw("SAME");
  leg1->Draw();
  lat.DrawLatex(60, 12.5, "K^{-} #rightarrow #pi^{-} #pi^{0}");

  c3->cd(2);
  gPad->SetMargin(0.13, 0.05, 0.12, 0.05);
  gPad->SetGrid();
  hframe2->Draw();
  for (auto g : gMu) g->Draw("P SAME");
  zero2->Draw("SAME");
  leg2->Draw();
  lat.DrawLatex(60, 12.5, "K^{-} #rightarrow #mu^{-} #bar{#nu}_{#mu}");

  //c3->SaveAs("dtof_vs_decay_momentum_combined.pdf");
  //c3->SaveAs("dtof_vs_decay_momentum_combined.png");
}

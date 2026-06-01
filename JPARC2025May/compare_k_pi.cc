#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "iostream"

void compare_k_pi()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const char* file_pass = "t110_graph_321.root";
  const char* file_kaon = "t110_graph_344.root";

  const char* hist_pass_name = "hist_bac_npe_s_pass";
  const char* hist_kaon_name = "hist_bac_npe_s_kaon";

  TFile* f_pass = TFile::Open(file_pass, "READ");
  TFile* f_kaon = TFile::Open(file_kaon, "READ");

  if (!f_pass || f_pass->IsZombie()) {
    std::cerr << "Cannot open " << file_pass << std::endl;
    return;
  }

  if (!f_kaon || f_kaon->IsZombie()) {
    std::cerr << "Cannot open " << file_kaon << std::endl;
    return;
  }

  TH1* h_pass = (TH1*)f_pass->Get(hist_pass_name);
  TH1* h_kaon = (TH1*)f_kaon->Get(hist_kaon_name);

  if (!h_pass) {
    std::cerr << "Cannot find " << hist_pass_name << " in " << file_pass << std::endl;
    return;
  }

  if (!h_kaon) {
    std::cerr << "Cannot find " << hist_kaon_name << " in " << file_kaon << std::endl;
    return;
  }

  // Clone so that histograms remain valid after closing files
  TH1* h_pass_n = (TH1*)h_pass->Clone("h_pass_n");
  TH1* h_kaon_n = (TH1*)h_kaon->Clone("h_kaon_n");

  h_pass_n->SetDirectory(nullptr);
  h_kaon_n->SetDirectory(nullptr);

  // Normalize to unit area
  double int_pass = h_pass_n->Integral();
  double int_kaon = h_kaon_n->Integral();

  if (int_pass > 0) h_pass_n->Scale(1.0 / int_pass);
  if (int_kaon > 0) h_kaon_n->Scale(1.0 / int_kaon);

  // Style
  h_pass_n->SetLineColor(kRed + 1);
  h_pass_n->SetLineWidth(3);
  h_pass_n->SetFillStyle(0);

  h_kaon_n->SetLineColor(kBlack);
  h_kaon_n->SetLineWidth(3);
  h_kaon_n->SetFillStyle(0);

  h_pass_n->GetXaxis()->SetTitle("N_{p.e.}");
  h_pass_n->GetYaxis()->SetTitle("Normalized counts");

  h_pass_n->GetXaxis()->SetRangeUser(-10, 90);

  double ymax = std::max(h_pass_n->GetMaximum(), h_kaon_n->GetMaximum());
  h_pass_n->SetMaximum(ymax * 1.25);
  h_pass_n->SetMinimum(0);
  h_pass_n->GetXaxis()->SetRangeUser(0,90);

  TCanvas* c1 = new TCanvas("c1", "Normalized Npe", 900, 700);
  c1->SetMargin(0.13, 0.05, 0.12, 0.05);

  h_pass_n->Draw("hist");
  h_kaon_n->Draw("hist same");

  TLegend* leg = new TLegend(0.58, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_pass_n, "Pion", "l");
  leg->AddEntry(h_kaon_n, "Kaon", "l");
  leg->Draw();

  c1->SaveAs("normalized_bac_npe_pion_kaon.pdf");
  c1->SaveAs("normalized_bac_npe_pion_kaon.png");

  f_pass->Close();
  f_kaon->Close();
}

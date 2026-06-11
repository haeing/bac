#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TPaletteAxis.h"
#include "TLine.h"
#include "iostream"

void compare_k_pi()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const char* file = "t110_graph_344_k_pi.root";
  const char* hist_name = "hist_bac_btof";
  
  TFile* f = TFile::Open(file, "READ");
  if (!f || f->IsZombie()) {
    std::cout << "Cannot open file: " << file << std::endl;
    return;
  }

  TH2D* hist_bac_btof = (TH2D*)f->Get(hist_name);
  if (!hist_bac_btof) {
    std::cout << "Cannot find histogram: " << hist_name << std::endl;
    return;
  }

  TH1D* hpx = hist_bac_btof->ProjectionX("hpx_bac_btof");
  TH1D* hpy = hist_bac_btof->ProjectionY("hpy_bac_btof");

  TCanvas* c = new TCanvas("c", "c", 1000, 850);

  // canvas layout
  const double x_y0   = 0.08;
  const double x_y1   = 0.28;

  const double x_2d0  = 0.28;
  const double x_2d1  = 0.90;

  const double x_cb0  = 0.90;
  const double x_cb1  = 0.94;

  const double y_2d0  = 0.08;
  const double y_2d1  = 0.75;

  const double y_x0   = 0.75;
  const double y_x1   = 0.95;

  TPad* pad_y  = new TPad("pad_y",  "Y projection", x_y0,  y_2d0, x_y1,  y_2d1);
  TPad* pad_2d = new TPad("pad_2d", "2D",           x_2d0, y_2d0, x_2d1, y_2d1);
  TPad* pad_x  = new TPad("pad_x",  "X projection", x_2d0, y_x0,  x_2d1, y_x1);

  pad_y->Draw();
  pad_2d->Draw();
  pad_x->Draw();

  // =====================
  // Y projection, left
  // =====================
  pad_y->cd();
  pad_y->SetLeftMargin(0.30);
  pad_y->SetRightMargin(0.00);
  pad_y->SetTopMargin(0.00);
  pad_y->SetBottomMargin(0.12);
  pad_y->SetTicks(1, 1);

  hpy->SetTitle("");
  hpy->SetLineColor(kBlue+1);
  hpy->SetLineWidth(2);

  hpy->GetYaxis()->SetRangeUser(hist_bac_btof->GetYaxis()->GetXmin(),
                                hist_bac_btof->GetYaxis()->GetXmax());

  // reverse count axis: max on left, 0 on right
  double ymax_count = hpy->GetMaximum() * 1.10;
  hpy->GetXaxis()->SetLimits(ymax_count, 0);

  hpy->GetXaxis()->SetTitle("Counts");
  hpy->GetXaxis()->SetTitleSize(0.06);
  hpy->GetXaxis()->SetTitleOffset(0.9);
  hpy->GetXaxis()->SetLabelSize(0.05);

  hpy->GetYaxis()->SetTitle(hist_bac_btof->GetYaxis()->GetTitle());
  hpy->GetYaxis()->SetTitleSize(0.08);
  hpy->GetYaxis()->SetTitleOffset(0.75);
  hpy->GetYaxis()->SetLabelSize(0.065);

  hpy->Draw("hist hbar");

  // =====================
  // Main 2D
  // =====================
  pad_2d->cd();
  pad_2d->SetLeftMargin(0.00);
  pad_2d->SetRightMargin(0.08);
  pad_2d->SetTopMargin(0.00);
  pad_2d->SetBottomMargin(0.12);
  pad_2d->SetTicks(1, 1);

  hist_bac_btof->SetTitle("");

  // remove y labels/title from main 2D
  hist_bac_btof->GetYaxis()->SetLabelSize(0);
  hist_bac_btof->GetYaxis()->SetTitleSize(0);
  hist_bac_btof->GetYaxis()->SetTickLength(0);

  hist_bac_btof->GetXaxis()->SetTitle(hist_bac_btof->GetXaxis()->GetTitle());
  hist_bac_btof->GetXaxis()->SetTitleSize(0.055);
  hist_bac_btof->GetXaxis()->SetLabelSize(0.055);

  hist_bac_btof->Draw("colz");

  c->Update();

  TPaletteAxis* palette =
    (TPaletteAxis*)hist_bac_btof->GetListOfFunctions()->FindObject("palette");

  if (palette) {
    palette->SetX1NDC(0.915);
    palette->SetX2NDC(0.955);
    palette->SetY1NDC(y_2d0);
    palette->SetY2NDC(y_2d1);
  }

  // =====================
  // X projection, top
  // =====================
  pad_x->cd();
  pad_x->SetLeftMargin(0.00);
  pad_x->SetRightMargin(0.08);
  pad_x->SetTopMargin(0.05);
  pad_x->SetBottomMargin(0.00);
  pad_x->SetTicks(1, 1);

  hpx->SetTitle("");
  hpx->SetLineColor(kBlue+1);
  hpx->SetLineWidth(2);

  hpx->GetXaxis()->SetRangeUser(hist_bac_btof->GetXaxis()->GetXmin(),
                                hist_bac_btof->GetXaxis()->GetXmax());

  // share x-axis with main: hide top x labels/title
  hpx->GetXaxis()->SetLabelSize(0);
  hpx->GetXaxis()->SetTitleSize(0);
  hpx->GetXaxis()->SetTickLength(0.04);

  hpx->GetYaxis()->SetTitle("");
  hpx->GetYaxis()->SetLabelSize(0.045);
  hpx->GetYaxis()->SetTitleSize(0);

  hpx->Draw("hist");

  c->cd();
  c->Update();

  c->SaveAs("compare_k_pi.pdf");
}

// draw_thr.C
void draw_thr()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFont(132, "XYZ");
  gStyle->SetLabelFont(132, "XYZ");

  TFile *fin = TFile::Open("t110_graph_344.root");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open file" << std::endl;
    return;
  }

  auto *g_thr = (TGraphAsymmErrors*)fin->Get("g_thr");
  if (!g_thr) {
    std::cerr << "Cannot find g_thr" << std::endl;
    return;
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
  //c1->SetMargin(0.13, 0.05, 0.12, 0.06);
  //c1->SetTicks(1, 1);

  TH1F *frame = c1->DrawFrame(-10, -0.05, 60, 1.05);
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("N_{p.e.}");
  frame->GetYaxis()->SetTitle("Efficiency");
  /*
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleOffset(1.10);
  frame->GetXaxis()->SetNdivisions(510);
  frame->GetYaxis()->SetNdivisions(510);
  */

  g_thr->SetMarkerStyle(20);
  g_thr->SetMarkerSize(1.2);
  g_thr->SetMarkerColor(kBlack);
  g_thr->SetLineColor(kBlack);
  //g_thr->SetLineWidth(1);
  g_thr->GetYaxis()->SetLimits(0.0, 1);
  //g_thr->Draw("P SAME");

  TF1 *fit = g_thr->GetFunction("turn");
  if (!fit) {
    // 그림에 보이는 sigmoid 형태. 필요하면 파라미터/fit range 조정.
    fit = new TF1("turn",
                  "[0]+[1]/(1.0+exp(-(x-[2])/[3]))",
                  -10, 60);
    fit->SetParameters(0.0, 1.0, 15.0, 2.0);
    g_thr->Fit(fit, "R0");
  }

  auto* threshold_line = new TLine(14.59, -0.05, 14.59, 1.05);
  threshold_line->SetLineStyle(2);
  threshold_line->SetLineWidth(2);
  threshold_line->Draw();

  fit->SetLineColor(kRed);
  fit->SetLineWidth(2);
  fit->SetRange(-10, 60);
  fit->Draw("SAME");

  g_thr->Draw("P SAME"); // 점을 fit 위에 다시 올림

  c1->Modified();
  c1->Update();

  c1->SaveAs("draw_thr.pdf");
  //c1->SaveAs("g_thr_efficiency.png");
}

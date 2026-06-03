void n_index_mom()
{
  
  const double n_aero = 1.115;
  const double pmin = 530.0;
  const double pmax = 1030.0;
  const double pcenter = 735.0;

  const double sep_p1 = 645.0;
  const double sep_p2 = 825.0;

  const double m_pi = 139.57039;
  const double m_K  = 493.677;

  auto n_threshold = [](double p, double m) {
    return std::sqrt(1.0 + (m/p)*(m/p));
  };

  const int N = 3000;
  TGraph* g_pi = new TGraph();
  TGraph* g_K  = new TGraph();

  for (int i = 0; i < N; ++i) {
    double p = 20.0 + (pmax - 20.0) * i / (N - 1.0);
    g_pi->SetPoint(i, p, n_threshold(p, m_pi));
    g_K ->SetPoint(i, p, n_threshold(p, m_K));
  }

  TCanvas* c = new TCanvas("c", "pi K threshold", 850, 520);
  c->SetMargin(0.12, 0.05, 0.12, 0.05);
  //c->SetGrid(1, 1);

  TH1D* frame = new TH1D("frame", "", 100, pmin, pmax);
  frame->SetMinimum(1.0);
  frame->SetMaximum(1.35);
  frame->GetXaxis()->SetTitle("Momentum [MeV/c]");
  frame->GetYaxis()->SetTitle("Refractive index");
  frame->Draw();

  // Threshold curves
  g_pi->SetLineColor(kBlue + 1);
  g_pi->SetLineWidth(3);
  g_pi->Draw("L SAME");

  g_K->SetLineColor(kOrange + 7);
  g_K->SetLineWidth(3);
  g_K->Draw("L SAME");

  // Aerogel refractive index line
  TLine* line_n = new TLine(pmin, n_aero, pmax, n_aero);
  line_n->SetLineColor(kBlack);
  line_n->SetLineStyle(2);
  line_n->SetLineWidth(3);
  line_n->Draw("SAME");

  

  // Red separation region
  double rect_height = 0.01;
  TBox* box = new TBox(sep_p1, n_aero - rect_height/2.0,
                       sep_p2, n_aero + rect_height/2.0);
  box->SetFillColorAlpha(kRed, 0.25);
  box->SetLineColor(0);
  box->Draw("SAME");

  // Vertical guide lines
  TLine* v1 = new TLine(sep_p1, 1.0, sep_p1, 1.35);
  v1->SetLineStyle(9);
  v1->SetLineWidth(2);
  v1->Draw("SAME");

  TLine* v2 = new TLine(sep_p2, 1.0, sep_p2, 1.35);
  v2->SetLineStyle(9);
  v2->SetLineWidth(2);
  v2->Draw("SAME");

  //Main mom line
  TLine* v3 = new TLine(pcenter, 1.0, pcenter, 1.35);
  //v2->SetLineStyle(9);
  v3->SetLineWidth(4);
  v3->SetLineColor(kRed);
  v3->Draw("SAME");

  TLatex latex;
  latex.SetTextSize(0.05);
  latex.SetTextColor(kRed);
  latex.SetTextFont(132);
  latex.DrawLatex(742, 1.305, "735 MeV/#it{c}");

  // Legend
  TLegend* leg = new TLegend(0.63, 0.78, 0.90, 0.92);
  leg->AddEntry(g_pi, "#pi threshold", "l");
  leg->AddEntry(g_K,  "K threshold", "l");
  leg->AddEntry(line_n, "Silica aerogel (n = 1.115)", "l");
  leg->Draw();

  c->SaveAs("pi_k_threshold.pdf");
}

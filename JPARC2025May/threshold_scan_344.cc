// Threshold scan for run 0344.
//
// This macro uses hist_bac_btof made by analysis_t110.cc with kaon = pion = true.
// BTOF defines the samples:
//   pion: -0.9 < BTOF < 1.3 ns
//   kaon:  4.5 < BTOF < 5.6 ns
// A BAC veto (NPE < threshold) is treated as the kaon selection.
// Hence, at every threshold:
//   kaon purity    = N_K(veto) / [N_K(veto) + N_pi(veto)]
//   pion efficiency = N_pi(BAC ON) / N_pi(all), where BAC ON means NPE >= threshold.
//
// Usage in ROOT:
//   root -l -b -q 'threshold_scan_344.cc+("t110_graph_344_k_pi.root", "threshold_scan_344.root")'

#include <algorithm>
#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TTree.h"
#include "TEfficiency.h"

namespace {
void SetBinomialErrors(long long pass, long long total, double& low, double& high)
{
  low = 0.0;
  high = 0.0;
  if (total <= 0 || pass < 0 || pass > total) return;

  TEfficiency interval("interval", "interval", 1, 0., 1.);
  interval.SetStatisticOption(TEfficiency::kFCP);
  interval.SetTotalEvents(1, total);
  interval.SetPassedEvents(1, pass);
  low = interval.GetEfficiencyErrorLow(1);
  high = interval.GetEfficiencyErrorUp(1);
}
}

void threshold_scan_344(const char* input_file = "t110_graph_344_k_pi.root",
                        const char* output_file = "threshold_scan_344.root",
                        double threshold_min = -10.,
                        double threshold_max = 60.)
{
  TFile input(input_file, "READ");
  if (input.IsZombie()) {
    std::cerr << "[ERROR] Cannot open " << input_file << std::endl;
    return;
  }

  TH2D* h_npe_btof = dynamic_cast<TH2D*>(input.Get("hist_bac_btof"));
  if (!h_npe_btof) {
    std::cerr << "[ERROR] hist_bac_btof is not in " << input_file << std::endl;
    return;
  }

  const TAxis* xaxis = h_npe_btof->GetXaxis();
  const TAxis* yaxis = h_npe_btof->GetYaxis();
  const int pion_y_first = yaxis->FindFixBin(-0.9 + 1.e-9);
  const int pion_y_last  = yaxis->FindFixBin( 1.3 - 1.e-9);
  const int kaon_y_first = yaxis->FindFixBin( 4.5 + 1.e-9);
  const int kaon_y_last  = yaxis->FindFixBin( 5.6 - 1.e-9);

  const int first_x = std::max(1, xaxis->FindFixBin(threshold_min + 1.e-9));
  const int last_x  = std::min(xaxis->GetNbins(), xaxis->FindFixBin(threshold_max - 1.e-9));
  if (first_x > last_x) {
    std::cerr << "[ERROR] Requested threshold range is outside the NPE histogram." << std::endl;
    return;
  }

  const long long pion_all = std::llround(h_npe_btof->Integral(1, xaxis->GetNbins(), pion_y_first, pion_y_last));
  const long long kaon_all = std::llround(h_npe_btof->Integral(1, xaxis->GetNbins(), kaon_y_first, kaon_y_last));
  if (pion_all == 0 || kaon_all == 0) {
    std::cerr << "[ERROR] Empty pion or kaon BTOF sample (Npi=" << pion_all
              << ", NK=" << kaon_all << ")." << std::endl;
    return;
  }

  TFile output(output_file, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "[ERROR] Cannot create " << output_file << std::endl;
    return;
  }

  TGraphAsymmErrors g_kaon_purity;
  g_kaon_purity.SetName("g_kaon_purity_vs_npe_threshold");
  g_kaon_purity.SetTitle("Kaon purity of BAC-veto selection;BAC threshold [N_{p.e.}];Kaon purity");
  TGraphAsymmErrors g_pion_efficiency;
  g_pion_efficiency.SetName("g_pion_efficiency_vs_npe_threshold");
  g_pion_efficiency.SetTitle("Pion efficiency of BAC-ON selection;BAC threshold [N_{p.e.}];Pion efficiency");
  TGraphAsymmErrors g_kaon_misid;
  g_kaon_misid.SetName("g_kaon_misid_as_pion_vs_npe_threshold");
  g_kaon_misid.SetTitle("Kaon mis-identification as pion;Threshold [Np.e.];Probability");

  TTree scan("threshold_scan", "Run 0344 BAC threshold scan");
  double threshold = 0., kaon_purity = 0., kaon_purity_err_low = 0., kaon_purity_err_high = 0.;
  double pion_efficiency = 0., pion_efficiency_err_low = 0., pion_efficiency_err_high = 0.;
  double kaon_misid_as_pion = 0., kaon_misid_err_low = 0., kaon_misid_err_high = 0.;
  long long n_kaon_veto = 0, n_pion_veto = 0, n_kaon_total = kaon_all, n_pion_total = pion_all;
  scan.Branch("threshold_npe", &threshold);
  scan.Branch("kaon_purity", &kaon_purity);
  scan.Branch("kaon_purity_err_low", &kaon_purity_err_low);
  scan.Branch("kaon_purity_err_high", &kaon_purity_err_high);
  scan.Branch("pion_efficiency", &pion_efficiency);
  scan.Branch("pion_efficiency_err_low", &pion_efficiency_err_low);
  scan.Branch("pion_efficiency_err_high", &pion_efficiency_err_high);
  scan.Branch("kaon_misid_as_pion", &kaon_misid_as_pion);
  scan.Branch("kaon_misid_err_low", &kaon_misid_err_low);
  scan.Branch("kaon_misid_err_high", &kaon_misid_err_high);
  scan.Branch("n_kaon_veto", &n_kaon_veto);
  scan.Branch("n_pion_veto", &n_pion_veto);
  scan.Branch("n_kaon_total", &n_kaon_total);
  scan.Branch("n_pion_total", &n_pion_total);

  int point = 0;
  for (int xbin = first_x; xbin <= last_x; ++xbin) {
    // The bin upper edge is the threshold: bins through xbin are the BAC-veto sample.
    threshold = xaxis->GetBinUpEdge(xbin);
    n_kaon_veto = std::llround(h_npe_btof->Integral(1, xbin, kaon_y_first, kaon_y_last));
    n_pion_veto = std::llround(h_npe_btof->Integral(1, xbin, pion_y_first, pion_y_last));
    const long long veto_total = n_kaon_veto + n_pion_veto;
    const long long pion_on = pion_all - n_pion_veto;
    const long long kaon_on = kaon_all - n_kaon_veto;

    kaon_purity = veto_total > 0 ? static_cast<double>(n_kaon_veto) / veto_total : 0.;
    pion_efficiency = static_cast<double>(pion_on) / pion_all;
    kaon_misid_as_pion = static_cast<double>(kaon_on) / kaon_all;
    SetBinomialErrors(n_kaon_veto, veto_total, kaon_purity_err_low, kaon_purity_err_high);
    SetBinomialErrors(pion_on, pion_all, pion_efficiency_err_low, pion_efficiency_err_high);
    SetBinomialErrors(kaon_on, kaon_all, kaon_misid_err_low, kaon_misid_err_high);

    g_kaon_purity.SetPoint(point, threshold, kaon_purity);
    g_kaon_purity.SetPointError(point, 0., 0., kaon_purity_err_low, kaon_purity_err_high);
    g_pion_efficiency.SetPoint(point, threshold, pion_efficiency);
    g_pion_efficiency.SetPointError(point, 0., 0., pion_efficiency_err_low, pion_efficiency_err_high);
    g_kaon_misid.SetPoint(point, threshold, kaon_misid_as_pion);
    g_kaon_misid.SetPointError(point, 0., 0., kaon_misid_err_low, kaon_misid_err_high);
    scan.Fill();
    ++point;
  }

  TNamed definition("selection_definition",
    "pi: -0.9 < BTOF < 1.3 ns; K: 4.5 < BTOF < 5.6 ns; "
    "K selection=BAC veto (NPE < threshold); pi selection=BAC ON (NPE >= threshold). "
    "Errors are 68.3% Clopper-Pearson binomial intervals.");
  TNamed source("source_histogram", "hist_bac_btof from analysis_t110.cc output");
  TCanvas c_tradeoff("c_pion_efficiency_vs_kaon_misid", "BAC threshold trade-off", 900, 700);
  TH1D frame("h_tradeoff_frame", ";Threshold [Np.e.];Probability", 1, threshold_min, threshold_max);
  frame.SetMinimum(-0.05);
  frame.SetMaximum(1.08);
  frame.Draw();
  g_pion_efficiency.SetMarkerStyle(20);
  g_pion_efficiency.SetMarkerColor(kBlack);
  g_pion_efficiency.SetLineColor(kBlack);
  g_kaon_misid.SetMarkerStyle(24);
  g_kaon_misid.SetMarkerColor(kBlue + 1);
  g_kaon_misid.SetLineColor(kBlue + 1);
  g_pion_efficiency.Draw("PLE SAME");
  g_kaon_misid.Draw("PLE SAME");
  auto* threshold_line = new TLine(14.59, -0.05, 14.59, 1.08);
  threshold_line->SetLineStyle(2);
  threshold_line->SetLineWidth(2);
  threshold_line->Draw();
  TLegend legend(0.48, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(&g_pion_efficiency, "pion efficiency", "pl");
  legend.AddEntry(&g_kaon_misid, "kaon misidentification", "pl");
  legend.Draw();
  h_npe_btof->Write("hist_bac_btof_source");
  g_kaon_purity.Write();
  g_pion_efficiency.Write();
  g_kaon_misid.Write();
  c_tradeoff.Write();
  scan.Write();
  definition.Write();
  source.Write();
  output.Close();

  std::cout << "Saved " << point << " threshold points to " << output_file << std::endl;
  std::cout << "BTOF samples: Npi=" << pion_all << ", NK=" << kaon_all << std::endl;
}

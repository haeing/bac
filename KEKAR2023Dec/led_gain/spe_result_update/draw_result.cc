#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <string>

const bool kek = false;

// -------------------------------------------
static int FindParIndexByName(const TF1* f, const char* parName)
{
  if (!f) return -1;
  for (int i = 0; i < f->GetNpar(); ++i) {
    const char* name = f->GetParName(i);
    if (name && std::string(name) == parName) return i;
  }
  return -1;
}
// -------------------------------------------

void draw_result()
{
  // ====== FIXED SETTINGS ======
  const int N_BOARD = 4;
  int N_CH    = 16;
  if(kek)N_CH = 3;
  

  const char* FILE_PREFIX = "spe_result";
  char* OUT_PDF = "all_mppc_fit.pdf";
  char* OUT_CSV = "Q1_pxt_all.csv";
  if(kek){
    OUT_PDF = "all_mppc_fit_kek.pdf";
    OUT_CSV = "Q1_pxt_all_kek.csv";
  }
  
  // ============================

  std::ofstream csv(OUT_CSV);
  if(!kek)csv << "board,channel,Q1,Q1_err,pxt,pxt_err,chi2,ndf,status\n";
  else if(kek)csv << "board,channel,HV,Q1,Q1_err,pxt,pxt_err,chi2,ndf,status\n";

  TCanvas* c = new TCanvas("c","c",900,700);
  c->SetGrid(0,0);

  // multipage pdf start
  c->Print(Form("%s[", OUT_PDF));

  std::cout << "board ch    Q1        dQ1        pxt       dpxt     chi2/ndf  status\n";
  std::cout << "-------------------------------------------------------------------\n";

  for (int b = 0; b < N_BOARD; ++b) {
    for (int ch = 0; ch < N_CH; ++ch) {

      std::string fname = Form("%s%d_%d.root", FILE_PREFIX, b, ch);
      if(kek) fname = Form("%s%d_hv%d.root", FILE_PREFIX, b, ch+56);
      if (gSystem->AccessPathName(fname.c_str())) {
        std::cerr << "[SKIP] missing file: " << fname << "\n";
        continue;
      }

      TFile* fin = TFile::Open(fname.c_str(), "READ");
      if (!fin || fin->IsZombie()) {
        std::cerr << "[SKIP] cannot open: " << fname << "\n";
        delete fin;
        continue;
      }

      std::string kHist = Form("mppc%d_%d_hist", b, ch);
      std::string kFunc = Form("mppc%d_%d_func", b, ch);
      std::string kRes  = Form("mppc%d_%d", b, ch);
      if(kek){
	kHist = Form("mppc%d_hv%d_hist", b, ch+56);
	kFunc = Form("mppc%d_hv%d_func", b, ch+56);
	kRes  = Form("mppc%d_hv%d", b, ch+56);
      }

      TH1* h = (TH1*)fin->Get(kHist.c_str());
      TF1* f = (TF1*)fin->Get(kFunc.c_str());
      TFitResult* r = (TFitResult*)fin->Get(kRes.c_str());

      if (!h || !f) {
        std::cerr << "[SKIP] missing objects in " << fname << "\n";
        fin->Close();
        delete fin;
        continue;
      }

      // ---- extract parameters ----
      int iQ1  = FindParIndexByName(f, "Q1");
      int iPxt = FindParIndexByName(f, "pxt");

      double Q1   = (iQ1  >= 0) ? f->GetParameter(iQ1) : 0;
      double eQ1  = (iQ1  >= 0) ? f->GetParError(iQ1)  : 0;
      double pxt  = (iPxt >= 0) ? f->GetParameter(iPxt): 0;
      double epxt = (iPxt >= 0) ? f->GetParError(iPxt) : 0;

      double chi2 = f->GetChisquare();
      int    ndf  = f->GetNDF();
      int    status = (r ? r->Status() : -1);

      // ---- print & save ----
      std::cout << Form("%5d %2d  %9.4f  %9.4f  %9.5f  %9.5f  %7.2f/%d   %d",
                        b, ch, Q1, eQ1, pxt, epxt, chi2, ndf, status)
                << "\n";

      if(!kek){
      csv << b << "," << ch << ","
          << Q1 << "," << eQ1 << ","
          << pxt << "," << epxt << ","
          << chi2 << "," << ndf << ","
          << status << "\n";
      }
      if(kek){
	csv << b << "," << 1 << ","<< ch+56 << ","
          << Q1 << "," << eQ1 << ","
          << pxt << "," << epxt << ","
          << chi2 << "," << ndf << ","
          << status << "\n";
      }

      // ---- draw ----
      c->Clear();
      h->SetStats(0);
      h->Draw("hist");
      f->SetLineWidth(2);
      f->SetNpx(5000);
      f->Draw("same");

      TLegend* leg = new TLegend(0.60,0.72,0.89,0.89);
      leg->SetBorderSize(0);
      leg->AddEntry(h, kHist.c_str(), "l");
      leg->AddEntry(f, kFunc.c_str(), "l");
      leg->Draw();

      TLatex lat;
      lat.SetNDC(true);
      lat.SetTextSize(0.035);
      lat.DrawLatex(0.12,0.92,Form("board %d  channel %d", b, ch));
      if(kek)lat.DrawLatex(0.12,0.92,Form("board %d  channel 1 HV %d", b, ch+56));
      lat.DrawLatex(0.12,0.87,Form("Q1  = %.4f #pm %.4f", Q1, eQ1));
      lat.DrawLatex(0.12,0.82,Form("pxt = %.5f #pm %.5f", pxt, epxt));
      lat.DrawLatex(0.12,0.77,Form("#chi^{2}/NDF = %.2f / %d  (status=%d)",
                                    chi2, ndf, status));

      c->Print(OUT_PDF);

      delete leg;
      fin->Close();
      delete fin;
    }
  }

  // multipage pdf end
  c->Print(Form("%s]", OUT_PDF));
  csv.close();

  std::cout << "\nDONE.\n";
  std::cout << "  PDF : " << OUT_PDF << "\n";
  std::cout << "  CSV : " << OUT_CSV << "\n";
}

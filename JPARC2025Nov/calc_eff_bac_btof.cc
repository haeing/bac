void calc_eff_bac_btof()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::vector<std::tuple<int,int,int>> runs = {
    {2815,790,1},
    {2610,715,1},
    //{2847,600,1},
    {2884,685,1},
    {3005,665,1},
    {2977,645,1}
  };

  const double nsigma = 3.0;

  std::vector<double> v_mom;
  std::vector<double> v_pi_eff, v_pi_err;
  std::vector<double> v_k_eff,  v_k_err;
  std::vector<double> v_ex;

  for(auto &[run,mom,use] : runs){
    if(!use) continue;

    TString fname = Form("e72_hist_%d_pion_kaon.root",run);
    TFile *f = TFile::Open(fname,"READ");

    if(!f || f->IsZombie()){
      std::cout << "Cannot open " << fname << std::endl;
      continue;
    }

    TH2D *h_tot  = (TH2D*)f->Get("hist_bac_btof");
    TH2D *h_pass = (TH2D*)f->Get("hist_bac_btof_pass");

    if(!h_tot || !h_pass){
      std::cout << "Missing hist in " << fname << std::endl;
      f->Close();
      continue;
    }

    TH1D *hy = h_tot->ProjectionY(Form("hy_run%d",run));

    TF1 *fpi = new TF1(Form("fpi_run%d",run),"gaus",-2.0,1.5);
    hy->Fit(fpi,"RQ0");

    double pi_mean  = fpi->GetParameter(1);
    double pi_sigma = fabs(fpi->GetParameter(2));

    TF1 *fk = new TF1(Form("fk_run%d",run),"gaus",3.0,6.5);
    hy->Fit(fk,"RQ0+");

    double k_mean  = fk->GetParameter(1);
    double k_sigma = fabs(fk->GetParameter(2));

    double pi_ymin = pi_mean - nsigma*pi_sigma;
    double pi_ymax = pi_mean + nsigma*pi_sigma;

    double k_ymin = k_mean - nsigma*k_sigma;
    double k_ymax = k_mean + nsigma*k_sigma;

    int xbin_min = 1;
    int xbin_max = h_tot->GetNbinsX();

    int pi_bin_min = h_tot->GetYaxis()->FindBin(pi_ymin);
    int pi_bin_max = h_tot->GetYaxis()->FindBin(pi_ymax);

    int k_bin_min = h_tot->GetYaxis()->FindBin(k_ymin);
    int k_bin_max = h_tot->GetYaxis()->FindBin(k_ymax);

    double pi_total = h_tot ->Integral(xbin_min,xbin_max,pi_bin_min,pi_bin_max);
    double pi_pass  = h_pass->Integral(xbin_min,xbin_max,pi_bin_min,pi_bin_max);

    double k_total = h_tot ->Integral(xbin_min,xbin_max,k_bin_min,k_bin_max);
    double k_pass  = h_pass->Integral(xbin_min,xbin_max,k_bin_min,k_bin_max);

    double pi_eff = (pi_total>0) ? pi_pass/pi_total : 0;
    double k_eff  = (k_total >0) ? k_pass /k_total  : 0;

    double pi_err = (pi_total>0) ? sqrt(pi_eff*(1.0-pi_eff)/pi_total) : 0;
    double k_err  = (k_total >0) ? sqrt(k_eff *(1.0-k_eff )/k_total ) : 0;

    v_mom.push_back(mom);
    v_ex.push_back(0.0);

    v_pi_eff.push_back(pi_eff);
    v_pi_err.push_back(pi_err);

    v_k_eff.push_back(k_eff);
    v_k_err.push_back(k_err);

    std::cout << run << "  "
              << mom << " MeV/c  "
              << "pi eff = " << pi_eff << " +/- " << pi_err
              << "   K eff = " << k_eff << " +/- " << k_err
              << std::endl;

    f->Close();
  }

  int n = v_mom.size();

  TGraphErrors *g_pi = new TGraphErrors(n);
  TGraphErrors *g_k  = new TGraphErrors(n);

  for(int i=0; i<n; i++){
    g_pi->SetPoint(i, v_mom[i], v_pi_eff[i]);
    g_pi->SetPointError(i, v_ex[i], v_pi_err[i]);

    g_k->SetPoint(i, v_mom[i], v_k_eff[i]);
    g_k->SetPointError(i, v_ex[i], v_k_err[i]);
  }

  TCanvas *c = new TCanvas("c_eff","c_eff",900,700);
  c->SetTicks(1,1);
  c->SetGrid(0,0);

  TH1F *frame = c->DrawFrame(580,0.0,810,1.08);
  frame->SetTitle(";Beam momentum [MeV/c];Efficiency");

  frame->GetXaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);

  g_pi->SetMarkerStyle(20);
  g_pi->SetMarkerSize(1.2);
  g_pi->SetMarkerColor(kBlue+1);
  g_pi->SetLineColor(kBlue+1);
  g_pi->SetLineWidth(2);

  g_k->SetMarkerStyle(21);
  g_k->SetMarkerSize(1.2);
  g_k->SetMarkerColor(kRed+1);
  g_k->SetLineColor(kRed+1);
  g_k->SetLineWidth(2);

  g_pi->Draw("P SAME");
  g_k->Draw("P SAME");

  TLegend *leg = new TLegend(0.18,0.18,0.50,0.32);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->AddEntry(g_pi,"#pi efficiency","pl");
  leg->AddEntry(g_k,"K efficiency","pl");
  leg->Draw();


}

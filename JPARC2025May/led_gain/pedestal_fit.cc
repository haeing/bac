void FitPedestalOff_Smooth(TH1* hOff)
{
  if(!hOff){ printf("hOff is null\n"); return; }

  double xPeak = hOff->GetBinCenter(hOff->GetMaximumBin());
  double rms   = hOff->GetRMS();

  double xmin = xPeak - 6.0*rms;
  double xmax = xPeak + 10.0*rms;

  // params:
  // [0] Ag, [1] mu, [2] sigma
  // [3] At, [4] tau
  // [5] x0 (tail turn-on center)
  // [6] w  (turn-on smoothness)

  TF1* fPed = new TF1("fPed_smooth",
		      "[0]*TMath::Gaus(x,[1],[2],false) + "
		      "[3]*(1.0/(1.0+TMath::Exp(-(x-[5])/[6])))*TMath::Exp(-(x-[5])/[4])",
		      xmin, xmax);



  double x0_init = xPeak + 1.5*rms;

  fPed->SetParameters(
    hOff->GetMaximum(),   // Ag
    xPeak,                // mu
    0.6*rms,              // sigma
    0.2*hOff->GetMaximum(), // At
    2.0*rms,              // tau
    x0_init,              // x0
    0.5*rms               // w
  );

  // limits for stability
  fPed->SetParLimits(2, 0.1, 50.0);                    // sigma
  fPed->SetParLimits(4, 0.1, 300.0);                   // tau
  fPed->SetParLimits(5, xPeak + 0.2*rms, xPeak + 5*rms); // x0
  fPed->SetParLimits(6, 0.05*rms, 2.0*rms);            // w

  TFitResultPtr r = hOff->Fit(fPed, "RS");

  TCanvas* c = new TCanvas("c_pedoff","Pedestal fit (LED off)",900,700);
  hOff->Draw("HIST");
  fPed->SetLineColor(kRed);
  fPed->SetLineWidth(2);
  fPed->Draw("SAME");

  printf("\n==== Smooth pedestal+tail fit ====\n");
  printf("mu    = %.6f\n", fPed->GetParameter(1));
  printf("sigma = %.6f\n", fPed->GetParameter(2));
  printf("tau   = %.6f\n", fPed->GetParameter(4));
  printf("x0    = %.6f\n", fPed->GetParameter(5));
  printf("w     = %.6f\n", fPed->GetParameter(6));
  printf("chi2/ndf = %.2f / %d\n", fPed->GetChisquare(), fPed->GetNDF());

  c->SaveAs("pedestal_off_fit_smooth.pdf");
}


void pedestal_fit(int runnumber, int board)
{
  TString dir = "/gpfs/group/had/sks/Users/haein/JPARC2025May_root";
  TFile *fin_ped = TFile::Open(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber),"READ");
  TFile *fin = TFile::Open(Form("%s/run00%d_Hodoscope.root",dir.Data(),runnumber),"READ");

  TTree *data_ped = nullptr;
  TTree *data = nullptr;
  fin_ped->GetObject("hodo", data_ped);
  fin->GetObject("hodo", data);

  vector<double> *bac_adc_u_ped = nullptr;
  TBranch *b_bac_adc_u_ped = data_ped->GetBranch("bac_adc_u");
  b_bac_adc_u_ped->SetAddress(&bac_adc_u_ped);

  vector<double> *bac_adc_u = nullptr;
  TBranch *b_bac_adc_u = data->GetBranch("bac_adc_u");
  b_bac_adc_u->SetAddress(&bac_adc_u);
  

  TH1D *hOff = new TH1D("hOff","hOff",500,0,500);
  TH1D *hOn = new TH1D("hOn","hOn",500,0,500);
  
  for(int n=0;n<data_ped->GetEntries();n++){
    data_ped->GetEntry(n);
    hOff->Fill((*bac_adc_u_ped)[board]);
  }

  FitPedestalOff_Smooth(hOff);

  
}

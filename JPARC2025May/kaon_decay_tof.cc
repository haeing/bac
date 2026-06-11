#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TStyle.h"

//////////////////////////////////////////////////////////////
// USER CONFIGURATION
//////////////////////////////////////////////////////////////

// histogram binning
const int    NBin_TOF = 300;
const double TOF_Min  = -2.0;
const double TOF_Max  = 10.0;

// detector distance
const double Ldet = 8000.0; // mm

// rear detector size

const double rearX = 115.0; // mm
const double rearY = 100.0; // mm
/*
const double rearX = 1000.0; // mm
const double rearY = 1000.0; // mm
*/

// beam profile sigma
const double beamSigmaX = 20.0; // mm
const double beamSigmaY = 20.0; // mm

// rear yields
const double N_pi_rear_target = 2.5e5;
const double N_K_rear_target  = 3.0e4;

// aerogel
const double nAerogel = 1.115;

//////////////////////////////////////////////////////////////

const double c_mm_ns = 299.792458;

// masses [MeV/c2]
const double mK   = 493.677;
const double mpi  = 139.57039;
const double mpi0 = 134.9768;
const double mmu  = 105.65837;
const double me   = 0.51099895;
const double mnu  = 0.0;

// K± ctau
const double cTauK = 3712.0; // mm

// branching ratios
const double BR_munu     = 0.6356;
const double BR_pipi0    = 0.2067;
const double BR_pipipi   = 0.0558;
const double BR_pipi0pi0 = 0.0176;
const double BR_pi0enu   = 0.0507;
const double BR_pi0munu  = 0.0335;

TRandom3 rnd(0);

std::vector<double> p_center;
std::vector<double> p_weight;

enum DecayMode {
    kMuNu,
    kPiPi0,
    kPiPiPi,
    kPiPi0Pi0,
    kPi0ENu,
    kPi0MuNu,
    kOther
};

//////////////////////////////////////////////////////////////

bool InRect(double x, double y, double sx, double sy)
{
    return (std::abs(x) < sx/2.0 &&
            std::abs(y) < sy/2.0);
}

double Beta(double p, double m)
{
    return p / std::sqrt(p*p + m*m);
}

double TOF(double path_mm, double p, double m)
{
    return path_mm / (Beta(p,m) * c_mm_ns);
}

bool PassAerogel(double p, double m)
{
  return Beta(p,m) > 1.0 / nAerogel;
  //return Beta(p,m) > 0.95;
}

//////////////////////////////////////////////////////////////

void LoadMomentumCSV(const char* csv)
{
    std::ifstream fin(csv);

    if (!fin.is_open()) {
        std::cerr << "Cannot open " << csv << std::endl;
        return;
    }

    std::string line;
    std::getline(fin, line);

    int bin;
    double center, value, error;
    char comma;

    while (fin >> bin >> comma
               >> center >> comma
               >> value >> comma
               >> error)
    {
        if (value > 0) {
            p_center.push_back(center);
            p_weight.push_back(value);
        }
    }

    std::cout << "Loaded momentum bins : "
              << p_center.size() << std::endl;
}

//////////////////////////////////////////////////////////////

double SampleMomentum()
{
    double sum = 0;

    for (auto w : p_weight)
        sum += w;

    double r = rnd.Uniform(sum);

    double acc = 0;

    for (size_t i = 0; i < p_weight.size(); i++) {

        acc += p_weight[i];

        if (r < acc)
            return p_center[i];
    }

    return p_center.back();
}

//////////////////////////////////////////////////////////////

double AverageKaonSurvival(int Ntest=200000)
{
    double sum = 0;

    for (int i = 0; i < Ntest; i++) {

        double pK = SampleMomentum();

        double meanDecay =
            (pK / mK) * cTauK;

        sum += std::exp(-Ldet / meanDecay);
    }

    return sum / Ntest;
}

//////////////////////////////////////////////////////////////

DecayMode SelectDecayMode()
{
    double r = rnd.Uniform();

    double acc = 0;

    acc += BR_munu;
    if (r < acc) return kMuNu;

    acc += BR_pipi0;
    if (r < acc) return kPiPi0;

    acc += BR_pipipi;
    if (r < acc) return kPiPiPi;

    acc += BR_pipi0pi0;
    if (r < acc) return kPiPi0Pi0;

    acc += BR_pi0enu;
    if (r < acc) return kPi0ENu;

    acc += BR_pi0munu;
    if (r < acc) return kPi0MuNu;

    return kOther;
}

//////////////////////////////////////////////////////////////

int PropagateDecayProducts(
    double x0,
    double y0,
    double zdecay,
    double pK,
    DecayMode mode,
    TH1D* hMode,
    TH1D* hModeBAC,
    TH1D* hBACAll
)
{
    double EK =
        std::sqrt(pK*pK + mK*mK);

    TLorentzVector kaon(
        0, 0, pK, EK
    );

    std::vector<double> masses;
    std::vector<bool> charged;

    //////////////////////////////////////////////////////////

    if (mode == kMuNu) {

        masses  = {mmu, mnu};
        charged = {true, false};

    } else if (mode == kPiPi0) {

        masses  = {mpi, mpi0};
        charged = {true, false};

    } else if (mode == kPiPiPi) {

        masses  = {mpi, mpi, mpi};
        charged = {true, true, true};

    } else if (mode == kPiPi0Pi0) {

        masses  = {mpi, mpi0, mpi0};
        charged = {true, false, false};

    } else if (mode == kPi0ENu) {

        masses  = {mpi0, me, mnu};
        charged = {false, true, false};

    } else if (mode == kPi0MuNu) {

        masses  = {mpi0, mmu, mnu};
        charged = {false, true, false};

    } else {

        return 0;
    }

    //////////////////////////////////////////////////////////

    int nD = masses.size();

    double massArray[4];

    for (int i = 0; i < nD; i++)
        massArray[i] = masses[i];

    TGenPhaseSpace decay;
    
    decay.SetDecay(
		   kaon,
		   nD,
		   massArray
		   );
    

    decay.Generate();

    //////////////////////////////////////////////////////////

    double tofRef = TOF(Ldet, 735.0, mpi);

    int nAccepted = 0;

    for (int i = 0; i < nD; i++) {

        if (!charged[i])
            continue;

        TLorentzVector* d =
            decay.GetDecay(i);

        double px = d->Px();
        double py = d->Py();
        double pz = d->Pz();
        double p  = d->P();

        if (pz <= 0)
            continue;

        //////////////////////////////////////////////////////

        double remainZ =
            Ldet - zdecay;

        double scale =
            remainZ / pz;

        double xRear =
            x0 + px * scale;

        double yRear =
            y0 + py * scale;

	
        if (!InRect(
                xRear,
                yRear,
                rearX,
                rearY))
            continue;
	

        //////////////////////////////////////////////////////

        double pathDaughter =
            remainZ * p / pz;

        double tof =
            TOF(zdecay, pK, mK)
          + TOF(pathDaughter,
                p,
                masses[i]);

        double dtof =
            tof - tofRef;

        //////////////////////////////////////////////////////

        hMode->Fill(dtof);

        if (PassAerogel(
                p,
                masses[i]))
        {
            hModeBAC->Fill(dtof);
            hBACAll->Fill(dtof);
        }

        nAccepted++;
    }

    return nAccepted;
}

//////////////////////////////////////////////////////////////

void kaon_decay_tof(
    const char* csv="E72_735.csv")
{
    gStyle->SetOptStat(0);

    //////////////////////////////////////////////////////////

    LoadMomentumCSV(csv);

    if (p_center.empty())
        return;

    //////////////////////////////////////////////////////////

    double avgSurvival =
        AverageKaonSurvival();

    int Npi =
        (int)N_pi_rear_target;

    int NK_initial =
        (int)std::round(
            N_K_rear_target
            / avgSurvival
        );

    //////////////////////////////////////////////////////////

    std::cout << "\n";
    std::cout << "Average K survival = "
              << avgSurvival
              << std::endl;

    std::cout << "Estimated initial K = "
              << NK_initial
              << std::endl;

    //////////////////////////////////////////////////////////

    TH1D* hPi = new TH1D(
        "hPi",
        "TOF at rear detector;"
        "TOF - TOF_{#pi(735)} [ns];"
        "Counts / spill",
        NBin_TOF,
        TOF_Min,
        TOF_Max
    );

    TH1D* hK = new TH1D(
        "hK",
        "Pure K",
        NBin_TOF,
        TOF_Min,
        TOF_Max
    );

    //////////////////////////////////////////////////////////

    TH1D* hMuNu =
        new TH1D(
            "hMuNu",
            "K #rightarrow #mu#nu",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPiPi0 =
        new TH1D(
            "hPiPi0",
            "K #rightarrow #pi#pi^{0}",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPiPiPi =
        new TH1D(
            "hPiPiPi",
            "K #rightarrow #pi#pi#pi",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPiPi0Pi0 =
        new TH1D(
            "hPiPi0Pi0",
            "K #rightarrow #pi#pi^{0}#pi^{0}",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPi0ENu =
        new TH1D(
            "hPi0ENu",
            "K #rightarrow #pi^{0}e#nu",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPi0MuNu =
        new TH1D(
            "hPi0MuNu",
            "K #rightarrow #pi^{0}#mu#nu",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    //////////////////////////////////////////////////////////
    // BAC accepted
    //////////////////////////////////////////////////////////

    TH1D* hMuNuBAC =
        new TH1D(
            "hMuNuBAC",
            "BAC K #rightarrow #mu#nu",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPiPi0BAC =
        new TH1D(
            "hPiPi0BAC",
            "BAC K #rightarrow #pi#pi^{0}",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPiPiPiBAC =
        new TH1D(
            "hPiPiPiBAC",
            "BAC K #rightarrow #pi#pi#pi",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPiPi0Pi0BAC =
        new TH1D(
            "hPiPi0Pi0BAC",
            "BAC K #rightarrow #pi#pi^{0}#pi^{0}",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPi0ENuBAC =
        new TH1D(
            "hPi0ENuBAC",
            "BAC K #rightarrow #pi^{0}e#nu",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    TH1D* hPi0MuNuBAC =
        new TH1D(
            "hPi0MuNuBAC",
            "BAC K #rightarrow #pi^{0}#mu#nu",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    //////////////////////////////////////////////////////////

    TH1D* hBACAll =
        new TH1D(
            "hBACAll",
            "All decay daughters with BAC light",
            NBin_TOF,
            TOF_Min,
            TOF_Max
        );

    //////////////////////////////////////////////////////////

    double tofRef =
        TOF(Ldet, 735.0, mpi);

    //////////////////////////////////////////////////////////
    // beam pion
    //////////////////////////////////////////////////////////

    for (int i = 0; i < Npi; i++) {

        double pPi =
            SampleMomentum();

        double dtof =
            TOF(Ldet, pPi, mpi)
            - tofRef;

        hPi->Fill(dtof);
    }

    //////////////////////////////////////////////////////////

    int nPureK = 0;
    int nDecay = 0;
    int nOther = 0;

    //////////////////////////////////////////////////////////

    for (int i = 0;
         i < NK_initial;
         i++)
    {
        double pK =
            SampleMomentum();

        double x0 =
            rnd.Gaus(
                0,
                beamSigmaX
            );

        double y0 =
            rnd.Gaus(
                0,
                beamSigmaY
            );

        if (!InRect(
                x0,
                y0,
                rearX,
                rearY))
            continue;

        //////////////////////////////////////////////////////

        double meanDecay =
            (pK / mK)
            * cTauK;

        double zdecay =
            rnd.Exp(meanDecay);

        //////////////////////////////////////////////////////

        if (zdecay > Ldet) {

            nPureK++;

            double dtof =
                TOF(
                    Ldet,
                    pK,
                    mK
                ) - tofRef;

            hK->Fill(dtof);

            continue;
        }

        //////////////////////////////////////////////////////

        nDecay++;

        DecayMode mode =
            SelectDecayMode();

        if (mode == kOther) {
            nOther++;
            continue;
        }

        //////////////////////////////////////////////////////

        TH1D* h = nullptr;
        TH1D* hBAC = nullptr;

        if (mode == kMuNu) {
            h = hMuNu;
            hBAC = hMuNuBAC;
        }
        else if (mode == kPiPi0) {
            h = hPiPi0;
            hBAC = hPiPi0BAC;
        }
        else if (mode == kPiPiPi) {
            h = hPiPiPi;
            hBAC = hPiPiPiBAC;
        }
        else if (mode == kPiPi0Pi0) {
            h = hPiPi0Pi0;
            hBAC = hPiPi0Pi0BAC;
        }
        else if (mode == kPi0ENu) {
            h = hPi0ENu;
            hBAC = hPi0ENuBAC;
        }
        else if (mode == kPi0MuNu) {
            h = hPi0MuNu;
            hBAC = hPi0MuNuBAC;
        }

        //////////////////////////////////////////////////////

        PropagateDecayProducts(
            x0,
            y0,
            zdecay,
            pK,
            mode,
            h,
            hBAC,
            hBACAll
        );
    }

    //////////////////////////////////////////////////////////
    // style
    //////////////////////////////////////////////////////////

    hPi->SetLineColor(kMagenta+2);
    hK->SetLineColor(kBlack);

    hMuNu->SetLineColor(kGreen+2);
    hPiPi0->SetLineColor(kBlue+1);
    hPiPiPi->SetLineColor(kOrange+7);
    hPiPi0Pi0->SetLineColor(kCyan+2);
    hPi0ENu->SetLineColor(kRed+1);
    hPi0MuNu->SetLineColor(kViolet+1);

    hBACAll->SetLineColor(kRed);
    hBACAll->SetLineWidth(3);
    hBACAll->SetLineStyle(2);

    //////////////////////////////////////////////////////////

    std::vector<TH1D*> hs = {
        hPi,
        hK,
        hMuNu,
        hPiPi0,
        hPiPiPi,
        hPiPi0Pi0,
        hPi0ENu,
        hPi0MuNu
    };

    for (auto h : hs)
        h->SetLineWidth(2);

    //////////////////////////////////////////////////////////
    // canvas
    //////////////////////////////////////////////////////////

    TCanvas* c1 =
        new TCanvas(
            "c1",
            "TOF full decay composition",
            1200,
            800
        );

    c1->SetLogy();

    hPi->GetYaxis()->SetRangeUser(1,1000000);
    hPi->Draw("hist");
    hK->Draw("hist same");

    hMuNu->Draw("hist same");
    hPiPi0->Draw("hist same");
    hPiPiPi->Draw("hist same");
    //hPiPi0Pi0->Draw("hist same");
    //hPi0ENu->Draw("hist same");
    //hPi0MuNu->Draw("hist same");

    hBACAll->Draw("hist same");

    //////////////////////////////////////////////////////////

    TLegend* leg =
        new TLegend(
            0.50,
            0.48,
            0.88,
            0.88
        );

    leg->AddEntry(hPi,
        "#pi", "l");

    leg->AddEntry(hK,
        "#it{K} surviving at BH2", "l");

    leg->AddEntry(hMuNu,
        "#it{K #rightarrow #mu#nu}", "l");

    leg->AddEntry(hPiPi0,
        "#it{K #rightarrow #pi#pi^{0}}", "l");

    leg->AddEntry(hPiPiPi,
		  "#it{K #rightarrow #pi#pi#pi}", "l");
    /*
    leg->AddEntry(hPiPi0Pi0,
        "K #rightarrow #pi#pi^{0}#pi^{0}", "l");

    leg->AddEntry(hPi0ENu,
        "K #rightarrow #pi^{0}e#nu", "l");

    leg->AddEntry(hPi0MuNu,
        "K #rightarrow #pi^{0}#mu#nu", "l");
    */
    leg->AddEntry(hBACAll,
		  "Decay dauthers above BAC Cherenkov threshold",
        "l");

    leg->Draw();

    //////////////////////////////////////////////////////////

    std::cout << "\n";
    std::cout << "========== Summary =========="
              << std::endl;

    std::cout << "Input rear pion yield : "
              << Npi
              << std::endl;

    std::cout << "Estimated initial K : "
              << NK_initial
              << std::endl;

    std::cout << "Pure K survived : "
              << nPureK
              << std::endl;

    std::cout << "K decayed before rear : "
              << nDecay
              << std::endl;

    std::cout << "Remaining others : "
              << nOther
              << std::endl;

    
    double TOFCut_min = 4.0;
    double TOFCut_max = 5.6;
    int binCut_min =
      hPi->FindBin(TOFCut_min);
    int binCut_max = hPi->FindBin(TOFCut_max);

    double N_BAC_TOF =
      hBACAll->Integral(binCut_min, binCut_max);

    double N_TOF =
      hK->Integral(binCut_min, binCut_max)
      + hMuNu->Integral(binCut_min, binCut_max)
      + hPiPi0->Integral(binCut_min, binCut_max)
      + hPiPiPi->Integral(binCut_min, binCut_max)
      + hPiPi0Pi0->Integral(binCut_min, binCut_max)
      + hPi0ENu->Integral(binCut_min, binCut_max)
      + hPi0MuNu->Integral(binCut_min, binCut_max);
    double frac =
    100.0 * N_BAC_TOF / N_TOF;
    
    cout<<"fraction : "<<frac<<endl;
 
}

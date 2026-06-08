// delta_ray_prob.C
// root -l -q 'delta_ray_prob.C'

#include <iostream>
#include <cmath>

const double me = 0.51099895;          // MeV
const double mK = 493.677;             // MeV/c2
const double re = 2.8179403262e-13;    // cm
const double NA = 6.02214076e23;

// BAC aerogel
const double nAerogel = 1.115;
const double rhoAerogel = 0.38;        // g/cm3, adjust if known
const double thicknessAerogel = (12.4 + 12.4 + 12.3) / 10.0; // cm

// approximate SiO2: Z/A ~ 0.499
const double ZoverA_Aerogel = 0.499;

// Kaon momentum
const double pK = 735.0; // MeV/c

double Beta(double p, double m)
{
    return p / std::sqrt(p*p + m*m);
}

double Gamma(double beta)
{
    return 1.0 / std::sqrt(1.0 - beta*beta);
}

double ElectronCherenkovThreshold(double n)
{
    double beta_th = 1.0 / n;
    double gamma_th = Gamma(beta_th);
    return (gamma_th - 1.0) * me; // MeV
}

double TmaxHeavyParticle(double p, double M)
{
    double beta = Beta(p, M);
    double gamma = Gamma(beta);

    double massRatio = me / M;

    return (2.0 * me * beta*beta * gamma*gamma)
        / (1.0 + 2.0*gamma*massRatio + massRatio*massRatio);
}

// Integral of delta-ray cross section above Tcut
// dσ/dT ~= 2π re^2 me / beta^2 * 1/T^2
// sigma(Tcut<T<Tmax) = C * (1/Tcut - 1/Tmax)
double DeltaRayCrossSection(double beta, double Tcut, double Tmax)
{
    if (Tmax <= Tcut) return 0.0;

    double C = 2.0 * M_PI * re * re * me / (beta * beta); // cm2 MeV

    return C * (1.0/Tcut - 1.0/Tmax); // cm2 / electron
}

void delta_ray_prob()
{
    double betaK = Beta(pK, mK);
    double gammaK = Gamma(betaK);

    double Tth = ElectronCherenkovThreshold(nAerogel);
    double Tmax = TmaxHeavyParticle(pK, mK);

    double sigma = DeltaRayCrossSection(betaK, Tth, Tmax);

    // electron column density [electrons/cm2]
    double Ne_col = rhoAerogel * thicknessAerogel * NA * ZoverA_Aerogel;

    double meanN = Ne_col * sigma;

    // Poisson probability for at least one delta ray above threshold
    double prob = 1.0 - std::exp(-meanN);

    std::cout << "========== Delta-ray estimate ==========\n";
    std::cout << "Kaon momentum              : " << pK << " MeV/c\n";
    std::cout << "beta_K                     : " << betaK << "\n";
    std::cout << "gamma_K                    : " << gammaK << "\n";
    std::cout << "Aerogel thickness          : " << thicknessAerogel << " cm\n";
    std::cout << "Aerogel density            : " << rhoAerogel << " g/cm3\n";
    std::cout << "Areal density rho*x        : "
              << rhoAerogel * thicknessAerogel << " g/cm2\n";
    std::cout << "Electron Cherenkov T_th    : " << Tth << " MeV\n";
    std::cout << "Delta-ray T_max from K     : " << Tmax << " MeV\n";
    std::cout << "Integrated cross section   : " << sigma << " cm2/electron\n";
    std::cout << "Mean N_delta above T_th    : " << meanN << "\n";
    std::cout << "P(at least one delta ray)  : " << prob*100.0 << " %\n";
}

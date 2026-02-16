#ifndef BACPrimaryGeneratorAction_E70_hh
#define BACPrimaryGeneratorAction_E70_hh

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <TRandom3.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>

class G4ParticleGun;
class G4Event;


class BACPrimaryGeneratorAction_E70 : public G4VUserPrimaryGeneratorAction
{
public:
  BACPrimaryGeneratorAction_E70();
  virtual ~BACPrimaryGeneratorAction_E70();
  virtual void GeneratePrimaries(G4Event *anEvent);
  const G4ParticleGun* GetParticleGun() const {return fParticleGun;}



private:

  G4ParticleGun* fParticleGun;
  G4ParticleTable* particleTable;
  const double mass_kaonm = 0.493677;
  const double mass_pim = 0.139570;
  const double mass_proton = 0.938272;
  G4String particle = "pion";
  G4double energy;


};

#endif

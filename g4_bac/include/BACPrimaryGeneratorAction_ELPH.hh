#ifndef BACPrimaryGeneratorAction_ELPH_hh
#define BACPrimaryGeneratorAction_ELPH_hh

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
class G4Box;

class BACPrimaryGeneratorAction_ELPH : public G4VUserPrimaryGeneratorAction
{
public:
  BACPrimaryGeneratorAction_ELPH(const G4String &parameter1, const G4String &parameter2);
  virtual ~BACPrimaryGeneratorAction_ELPH();
  virtual void GeneratePrimaries(G4Event *anEvent);
  const G4ParticleGun* GetParticleGun() const {return fParticleGun;}



private:

  G4double pa1;
  G4double pa2;
  
  G4ParticleGun* fParticleGun;
  G4ParticleTable* particleTable;
  G4double energy;

  G4String parameter1;
  G4String parameter2;

};

#endif

#ifndef BACPrimaryGeneratorAction_T110_hh
#define BACPrimaryGeneratorAction_T110_hh

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <TRandom3.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>

class G4ParticleGun;
class G4Event;
class G4Box;

class BACPrimaryGeneratorAction_T110 : public G4VUserPrimaryGeneratorAction
{
public:
  BACPrimaryGeneratorAction_T110();
  virtual ~BACPrimaryGeneratorAction_T110();
  virtual void GeneratePrimaries(G4Event *anEvent);
  const G4ParticleGun* GetParticleGun() const {return fParticleGun;}



private:

  //G4double pa1;
  //G4double pa2;
  
  G4ParticleGun* fParticleGun;
  G4ParticleTable* particleTable;
  Int_t bp_file_ndata;
  Int_t bp_nAccess=0;
  const double mass_kaonm = 0.493677;
  const double mass_pim = 0.139570;
  G4String particle = "pion";
  G4double energy;

  TFile* fInputFile;
  TTree* fTree;

  Long64_t fNentries;

  // branches
  double fX;
  double fY;
  double fZ;
  double fU0;
  double fV0;
  double fP;

  //G4String parameter1;
  //G4String parameter2;

};

#endif

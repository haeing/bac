#include "BACPrimaryGeneratorAction_T110.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"


#include "TFile.h"
#include "TTree.h"
#include "Randomize.hh"

BACPrimaryGeneratorAction_T110::BACPrimaryGeneratorAction_T110()
  : G4VUserPrimaryGeneratorAction(),
    fInputFile(nullptr),
    fTree(nullptr),
    fNentries(0),
    fX(0.), fY(0.), fZ(0.),
    fU0(0.), fV0(0.), fP(0.)
{
  G4int n_particle = 1;
  
  fParticleGun = new G4ParticleGun(n_particle);

  fInputFile = TFile::Open("/Users/ihaein/Work/E72/j-parc-haein/e72/bac/JPARC2025May/t110_beam_735.root", "READ");
  if (!fInputFile || fInputFile->IsZombie()) {
    G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction",
                "FileError", FatalException,
                "Cannot open beam.root");
  }

  fTree = dynamic_cast<TTree*>(fInputFile->Get("tree"));
  if (!fTree) {
    G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction",
                "TreeError", FatalException,
                "Cannot find TTree 'tree'");
  }

  // Set branch addresses
  fTree->SetBranchAddress("x",  &fX);
  fTree->SetBranchAddress("y",  &fY);
  fTree->SetBranchAddress("z",  &fZ);
  fTree->SetBranchAddress("u0", &fU0);
  fTree->SetBranchAddress("v0", &fV0);
  fTree->SetBranchAddress("p",  &fP);

  fNentries = fTree->GetEntries();
  if (fNentries <= 0) {
    G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction",
                "TreeEmpty", FatalException,
                "Input tree has no entries");
  }
}

BACPrimaryGeneratorAction_T110::~BACPrimaryGeneratorAction_T110()
{
  
  delete fParticleGun;
  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
  }
}

void BACPrimaryGeneratorAction_T110::GeneratePrimaries(G4Event* anEvent)
{

  Long64_t evtID = anEvent->GetEventID();

  // 방법 1: event 수보다 tree entry가 적으면 반복 사용
  Long64_t entry = evtID % fNentries;
  
  //G4int event_num = anEvent->GetEventID()+1;
  if(evtID% 1000 == 0)
    {
      G4cout<<"Event# "<<evtID<<G4endl;
    }




  

  
  fTree->GetEntry(entry);

  G4ThreeVector dir(fU0, fV0, 1.0);
  dir = dir.unit();
  
  fParticleGun->SetParticlePosition(G4ThreeVector(fX*mm, fY*mm, fZ*mm));
  
  fParticleGun->SetParticleMomentumDirection(dir);
  double momentum = fP*735/0.9*MeV;
  //fParticleGun->SetParticleMomentum(momentum);
  
  if(particle=="kaon"){
    energy = (sqrt(mass_kaonm*mass_kaonm+momentum*momentum) - mass_kaonm )*GeV;
    fParticleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  }
  if(particle=="pion"){
    energy = (sqrt(mass_pim*mass_pim+momentum*momentum) - mass_pim )*GeV;
    fParticleGun->SetParticleDefinition (particleTable -> FindParticle("pi-"));
  }
  fParticleGun->SetParticleEnergy(energy);
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
}

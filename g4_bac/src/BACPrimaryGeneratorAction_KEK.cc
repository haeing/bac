#include "BACPrimaryGeneratorAction_KEK.hh"
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

BACPrimaryGeneratorAction_KEK::BACPrimaryGeneratorAction_KEK(const G4String &par1, const G4String &par2)
  : G4VUserPrimaryGeneratorAction(),parameter1(par1),parameter2(par2)
{
  G4int n_particle = 1;
  
  fParticleGun = new G4ParticleGun(n_particle);

}

BACPrimaryGeneratorAction_KEK::~BACPrimaryGeneratorAction_KEK()
{
  
  delete fParticleGun;
}

void BACPrimaryGeneratorAction_KEK::GeneratePrimaries(G4Event* anEvent)
{

  G4int event_num = anEvent->GetEventID()+1;
  if(event_num% 1000 == 0)
    {
      G4cout<<"Event# "<<event_num<<G4endl;
    }
  
  //homogeneous test------------------------------------

  pa1 = std::stod(parameter1);
  pa2 = std::stod(parameter2);

  G4double x_moving = pa1*mm;
  G4double y_moving = pa2*mm;
  

  G4double tight_size_x = 115;
  G4double tight_size_y = 115;


  energy = 2010*MeV;
  fParticleGun->SetParticleDefinition (particleTable -> FindParticle("e-"));

    
  
  fParticleGun->SetParticleMomentumDirection ( G4ThreeVector(0,0,1) );
  


  G4double x = tight_size_x*0.5-G4UniformRand()*tight_size_x*mm;
  G4double y = tight_size_y*0.5-G4UniformRand()*tight_size_y*mm;

  fParticleGun->SetParticlePosition(G4ThreeVector(x+x_moving,y+y_moving,-5*cm) );
  //fParticleGun->SetParticlePosition(G4ThreeVector(x_moving,y_moving,-5*cm) );
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);


  


}


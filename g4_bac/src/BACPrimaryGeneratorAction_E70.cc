#include "BACPrimaryGeneratorAction_E70.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"


#include "TFile.h"
#include "TTree.h"
#include "Randomize.hh"
#include "TRandom.h"

double BAC_ra_th = 24.7;
double BAC_z_pos = 3540.3-28.0-110.;
double BAC_width = 115;

double get_x_center(double z_pos){
  double a = (3540.3-3040.3)/(595.4-541.7);
  double b = 3040.3-a*541.7;

  return (z_pos-b)/a;
  
}


BACPrimaryGeneratorAction_E70::BACPrimaryGeneratorAction_E70()
  : G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  
  fParticleGun = new G4ParticleGun(n_particle);

}

BACPrimaryGeneratorAction_E70::~BACPrimaryGeneratorAction_E70()
{
  
  delete fParticleGun;
}

void BACPrimaryGeneratorAction_E70::GeneratePrimaries(G4Event* anEvent)
{

  G4int event_num = anEvent->GetEventID()+1;
  if(event_num% 1000 == 0)
    {
      G4cout<<"Event# "<<event_num<<G4endl;
    }



  
  //homogeneous test------------------------------------

  //G4double momentum = 1.2;
  G4double momentum = 0.735;

  if(particle=="pion"){
    energy = (sqrt(mass_pim*mass_pim+momentum*momentum) - mass_pim )*GeV;
    fParticleGun->SetParticleDefinition (particleTable -> FindParticle("pi+"));
  }
  
  else if(particle=="kaon"){
    energy = (sqrt(mass_kaonm*mass_kaonm+momentum*momentum) - mass_kaonm )*GeV;
    fParticleGun->SetParticleDefinition (particleTable -> FindParticle("kaon+"));
  }

  TFile *beam70_file = new TFile("../run70274_SdcOutTracking.root", "read");
  TTree *sdcout = (TTree*)beam70_file->Get("sdcout");
  Int_t ntrack;
  Double_t x0[10];
  Double_t y0[10];
  Double_t u0[10];
  Double_t v0[10];

  sdcout->SetBranchAddress("ntrack",&ntrack);
  sdcout->SetBranchAddress("x0",x0);
  sdcout->SetBranchAddress("y0",y0);
  sdcout->SetBranchAddress("u0",u0);
  sdcout->SetBranchAddress("v0",v0);

  int total = sdcout->GetEntries();
  TRandom r;

  int entry;

  G4double x;
  G4double y;
  
  entry = std::trunc(G4UniformRand()*(total-1));

  

  sdcout->GetEntry(entry);
    

  x = (x0[0]+u0[0]*(2100.6226+BAC_z_pos-5)-get_x_center(BAC_z_pos-5))*mm;
  y = (y0[0]+v0[0]*(2100.6226+BAC_z_pos-5))*mm;

  if(x<-60||x>60||y<-60||y>60){

    while((x<-60||x>60)||(y<-60||y>60)){
      entry = std::trunc(G4UniformRand()*(total-1));
      sdcout->GetEntry(entry);
      x = (x0[0]+u0[0]*(2100.6226+BAC_z_pos-5)-get_x_center(BAC_z_pos-5))*mm;
      y = (y0[0]+v0[0]*(2100.6226+BAC_z_pos-5))*mm;
      
    }
  }


  G4double tight_size_x = 110;
  G4double tight_size_y = 110;
  fParticleGun->SetParticleMomentumDirection ( G4ThreeVector(0,0,1) );
  

  G4double uniform_x = tight_size_x*0.5-G4UniformRand()*tight_size_x*mm;
  G4double uniform_y = tight_size_y*0.5-G4UniformRand()*tight_size_y*mm;
  //fParticleGun->SetParticleMomentumDirection (G4ThreeVector(u0[0],v0[0],1) );

  //fParticleGun->SetParticlePosition(G4ThreeVector(x,y,-12.3*mm/2-5*mm) );
  fParticleGun->SetParticlePosition(G4ThreeVector(uniform_x,uniform_y,-12.3*mm/2-5*mm) );
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);

  beam70_file->Close();
  


}


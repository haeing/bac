#include "MPPCSD.hh"
#include "MPPCHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include "TGraph.h"

MPPCSD::MPPCSD(G4String name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("MppcCollection");
  SetQETable();
}

MPPCSD::~MPPCSD()
{}

void MPPCSD::Initialize(G4HCofThisEvent *HCTE)
{
  static int HCID = -1;
  MppcCollection = new MPPCHitsCollection(GetName(), collectionName[0]);

  if (HCID<0) HCID = GetCollectionID(0);

  HCTE->AddHitsCollection(HCID, MppcCollection);


  

}

G4bool MPPCSD::ProcessHits(G4Step *astep, G4TouchableHistory *ROhist)
{
  
  const G4StepPoint* preStepPoint = astep-> GetPreStepPoint();
  G4Track* atrack = astep->GetTrack();
  G4int pid = astep->GetTrack()->GetDefinition()-> GetPDGEncoding();
  G4ThreeVector worldPos =preStepPoint->GetPosition();
  G4ThreeVector pos = preStepPoint->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  G4double energy = atrack->GetTotalEnergy();
  G4ThreeVector hitmom = atrack->GetMomentum()/CLHEP::eV;
  G4double E_p = hitmom.mag();

  const G4double h = 6.628e-34;
  const G4double c = 3.0e+8;
  G4double wavelength = ((h*c)/(energy*1.6e-13))*(1e+9); //nm

  G4int copynumber = preStepPoint->GetTouchableHandle()->GetCopyNumber();

  
  
  
  auto hitTime = astep->GetPreStepPoint()->GetGlobalTime();


  atrack->SetTrackStatus(fStopAndKill);

  G4double random = G4UniformRand();

  G4double qe_value = QETable->Eval(E_p);
  //std::cout<<qe_value<<std::endl;
  if (qe_value < 0 || qe_value > 1) qe_value = 0;

  if(random<qe_value){
    MPPCHit* ahit = new MPPCHit(pos,worldPos, hitTime,pid,wavelength, copynumber);
    MppcCollection->insert(ahit);
  }
  return true;
      

}


void MPPCSD::EndOfEvent(G4HCofThisEvent *HCTE)
{
  MppcCollection->PrintAllHits();

}

void MPPCSD::SetQETable(){
  const int num_qe = 44;

  
  G4double mppc_ep1[num_qe] = {1.38,1.43,1.47,1.51,1.56,1.61,1.66,1.7,1.74,1.79,1.84,1.88,1.93,1.97,2,2.06,2.1,2.15,2.19,2.24,2.3,2.33,2.4,2.47,2.57,2.7,2.85,2.96,3.05,3.15,3.22,3.28,3.35,3.41,3.5,3.57,3.65,3.67,3.71,3.73,3.77,3.79,3.83,3.9};
  G4double mppc_effi[num_qe] = {0.035,0.048,0.057,0.07,0.085,0.098,0.113,0.126,0.142,0.158,0.172,0.191,0.206,0.226,0.243,0.258,0.276,0.294,0.308,0.326,0.344,0.356,0.371,0.385,0.395,0.399,0.392,0.376,0.360,0.342,0.321,0.300,0.278,0.251,0.228,0.201,0.175,0.141,0.120,0.098,0.079,0.059,0.039,0.021};

  //for(int i=0;i<num_qe;i++)mppc_effi[i]*=47./40.; //57
  for(int i=0;i<num_qe;i++)mppc_effi[i]*=51./40.; //58

  QETable = new TGraph(num_qe,mppc_ep1,mppc_effi);
  //QETable->SetName("MPPC_eff");
  

  
}
  

#include "BACStackingAction.hh"

#include <TMath.h>
#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ClassificationOfNewTrack.hh"

int gCerenkovCounter;
int check;
int num_check;
double pid_sum;
double decay_check;


BACStackingAction::BACStackingAction()
  : G4UserStackingAction(),
    fScintillationCounter(0), fCerenkovCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BACStackingAction::~BACStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
BACStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{

  G4ClassificationOfNewTrack classification = fWaiting;
  const G4double h = 6.628e-34;
  const G4double c = 3.0e+8;


  //for view
  if (aTrack->GetDefinition() == G4OpticalPhoton::Definition()) {

     
        if (G4UniformRand() > 0.2) {
            return fKill;
        }
    }

  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
    if(aTrack->GetTouchable()->GetVolume()->GetCopyNo() ==123){
      { // particle is optical photon
	
	if(aTrack->GetParentID()>0)
	  { // particle is secondary
	    if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
	      fScintillationCounter++;
	    if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
	      if(((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9)>320&&((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9)<900){
		//if(((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9)>190&&((h*c)/(aTrack->GetTotalEnergy()*1.6e-13))*(1e+9)<900){
		fCerenkovCounter++;
		gCerenkovCounter++;
		check = 1;
		//std::cout<<"Cherenkov"<<std::endl;
	      }
	      else return fKill;
	    }
	  }
      }
    }
  }
  
  if(check ==1){
    if(aTrack->GetTouchable()->GetVolume()->GetCopyNo() ==123){
      if(aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()){
      num_check++;
      pid_sum+=aTrack->GetDefinition()->GetPDGEncoding();
      }
    }
  }
  
  /*
    G4Track* tr = (G4Track*) aTrack;
    G4String volume;
  
  
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  
  if(particle==G4PionMinus::PionMinus()){
    volume = aTrack->GetNextVolume()->GetName();
    //if(volume == "Part2LW")
    //tr->SetTrackStatus(fStopButAlive);
    std::cout<<" "<<particle<<std::endl;
  }
  */

  
  //auto *Change = new BACChange();
  //if(aTrack == nullptr)std::cout<<"error"<<std::endl;
  //Change->Analysis(*aTrack);
  //Change->ShowInfo();
		   
  return fUrgent;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BACStackingAction::NewStage()
{
  //std::cout << "Number of Scintillation photons produced in this event : "
  //	    << fScintillationCounter << std::endl;
  //std::cout << "Number of Cerenkov photons produced in this event : "
  // << fCerenkovCounter << std::endl;
  if(check==0)decay_check = 12283249;
  else if(check==1)decay_check = pid_sum/num_check;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BACStackingAction::PrepareNewEvent()
{
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
  gCerenkovCounter = 0;
  check = 0;
  num_check = 0;
  pid_sum = 0.;

}

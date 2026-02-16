#include "BACRunAction.hh"
#include "BACAnalysisManager.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <stdlib.h>

BACRunAction::BACRunAction(BACAnalysisManager *analysisManager)
  :G4UserRunAction(), anaMan(analysisManager)
{
  
}
	     

BACRunAction::~BACRunAction()
{
}

void BACRunAction::BeginOfRunAction(const G4Run *aRun)
{
  if(G4VVisManager::GetConcreteInstance()){
    G4UImanager *UI = G4UImanager::GetUIpointer();
  }
  if(anaMan) anaMan->BeginOfRun(aRun);
}


  


void BACRunAction::EndOfRunAction(const G4Run *aRun)
{
  if(G4VVisManager::GetConcreteInstance()){
    G4UImanager *UI = G4UImanager::GetUIpointer();
  }
  if (anaMan) anaMan->EndOfRun(aRun);
}



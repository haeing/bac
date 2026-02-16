#include "BACEventAction.hh"
#include "BACAnalysisManager.hh"

BACEventAction::BACEventAction(BACAnalysisManager *analysisManager)
  : G4UserEventAction(), anaMan(analysisManager)
{
}

BACEventAction::~BACEventAction()
{
}

void BACEventAction::BeginOfEventAction(const G4Event *anEvent)
{
  if (anaMan) anaMan->BeginOfEvent(anEvent);
}

void BACEventAction::EndOfEventAction(const G4Event *anEvent)
{
  if (anaMan) anaMan->EndOfEvent(anEvent);
}

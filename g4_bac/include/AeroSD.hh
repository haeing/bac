#ifndef AeroSD_h
#define AeroSD_h 1

#include "G4VSensitiveDetector.hh"
#include "AeroHit.hh"
//#include "TGraph.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class AeroSD : public G4VSensitiveDetector
{
private:
  //int EMFlag;
  AeroHitsCollection *AeroCollection;

public:
  AeroSD(G4String name);
  //AeroSD(G4String name, ParamMan*);
  virtual ~AeroSD();

  virtual void Initialize( G4HCofThisEvent *HCTE);
  virtual G4bool ProcessHits(G4Step* astep, G4TouchableHistory *ROhist);
  virtual void EndOfEvent(G4HCofThisEvent *HCTE);
  //virtual void EndOfEvent(G4HCofThisEvent*);
  
  //void EndOfEvent( G4HCofThisEvent *HCE );
  //G4bool ProcessHits_constStep(const G4Step*, G4TouchableHistory* );
  //G4int QEFlag;
  //TGraph* QETable;
  //void SetQETable(G4int);

  //void PrintAll() const;



};

#endif

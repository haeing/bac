#ifndef MPPCSD_h
#define MPPCSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MPPCHit.hh"
#include "TGraph.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class MPPCSD : public G4VSensitiveDetector
{
public:
  MPPCSD(G4String name);
  //MPPCSD(G4String name, ParamMan*);
  virtual ~MPPCSD();

  virtual void Initialize( G4HCofThisEvent *HCTE );
  virtual G4bool ProcessHits(G4Step* astep, G4TouchableHistory *ROhist);
  virtual void EndOfEvent(G4HCofThisEvent *HCTE);

  void SetQETable();
  TGraph* QETable;
  

private:
  MPPCHitsCollection *MppcCollection;

};

#endif

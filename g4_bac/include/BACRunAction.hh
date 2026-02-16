#ifndef BACRunAction_hh
#define BACRunAction_hh

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"



#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>

class G4Run;
class BACAnalysisManager;


class BACRunAction : public G4UserRunAction
{
public:
  BACRunAction(BACAnalysisManager *analysisManager = 0);
  virtual ~BACRunAction();

  virtual void BeginOfRunAction(const G4Run *aRun);
  virtual void EndOfRunAction(const G4Run *aRun);

private:
  BACAnalysisManager *anaMan;
};

#endif

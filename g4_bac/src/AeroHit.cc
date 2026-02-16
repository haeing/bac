#include "AeroHit.hh"
//#include "BACDetectorConstruction_167_2.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AttDefStore.hh"

#include <iomanip>

G4Allocator<AeroHit> AeroHitAllocator;

AeroHit::AeroHit()
  :xyz(0.,0.,0.),momentum(0.,0.,0.),tof(0.)
{}

AeroHit::AeroHit(G4ThreeVector& axyz, G4ThreeVector& mxyz, G4double t)
  :xyz(axyz),momentum(mxyz),tof(t)
{}

AeroHit::~AeroHit()
{}

AeroHit::AeroHit(G4ThreeVector &axyz, G4ThreeVector &mxyz, G4double t, G4int pid, G4double angle, G4int cn, G4double edep)
  :xyz(axyz),momentum(mxyz), tof(t),particleID(pid),angleget(angle),copynumber(cn),energyloss(edep)
{}
/*
const AeroHit& AeroHit::operator=(const AeroHit &right)
{
  fEdep = right.fEdep;
  fTime = right.fTime;
  fId = right.fId;
  return *this;
}



const std::map<G4String,G4AttDef>* AeroHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("AeroHit",isNew);

  if (isNew){
    //(*store)["Energy"] = G4AttDef("Energy","Energy Deposited","Physics","G4BestUnit","G4double");
    (*store)["Energy"] = G4AttDef("Energy","Energy Deposited","Physics","MeV","G4double");
    (*store)["ID"] = G4AttDef("ID","ID","Physics","","G4int");
    (*store)["Time"] = G4AttDef("Time","Time","Physics","ns","G4double");
  }
  return store;
}

std::vector<G4AttValue>* AeroHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  values ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));
  values ->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
  values ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fId),""));
  return values;
}



*/

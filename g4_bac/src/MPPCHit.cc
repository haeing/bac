#include "MPPCHit.hh"
//#include "BACDetectorConstruction.hh"

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

G4Allocator<MPPCHit> MPPCHitAllocator;

MPPCHit::MPPCHit()
  :xyz(0.,0.,0.),tof(0.)
{}

MPPCHit::MPPCHit(G4ThreeVector& axyz, G4double t)
  :xyz(axyz),tof(t)
{}

MPPCHit::~MPPCHit()
{}

MPPCHit::MPPCHit(G4ThreeVector &axyz, G4ThreeVector &wxyz, G4double t, G4int pid, G4double Wavelength, G4int cn)
  :xyz(axyz), worldxyz(wxyz), tof(t),particleID(pid), wavelengthMP(Wavelength), copynum(cn)
{}

/*
const MPPCHit& MPPCHit::operator=(const MPPCHit &right)
{
  fTime = right.fTime;
  //fId = right.fId;
  fPhotons = right.fPhotons;
  fPos = right.fPos;
  return *this;
}



const std::map<G4String,G4AttDef>* MPPCHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("MPPCHit",isNew);

  if (isNew){
    //(*store)["ID"] = G4AttDef("ID","ID","Physics","","G4int");
    (*store)["Time"] = G4AttDef("Time","Time","Physics","ns","G4double");
    (*store)["Hits"] = G4AttDef("Hit","Photon Hit","Physics","G4BestUnit","G4double");
    (*store)["Pos"] = G4AttDef("Pos","Position","Physics","mm""G4ThreeVector");
  }
  return store;
}

std::vector<G4AttValue>* MPPCHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  values ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));
  values ->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
  //values ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fId),""));
  values->push_back(G4AttValue("Hits",G4BestUnit(fPhotons,"Hits"),""));
  values ->push_back(G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));
  return values;
}

*/




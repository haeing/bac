#ifndef BACDetectorConstruction_edep_hh
#define BACDetectorConstruction_edep_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4VisAttributes;


class BACDetectorConstruction_edep : public G4VUserDetectorConstruction
{
public:

  BACDetectorConstruction_edep(const G4String &parameter4);
  virtual ~BACDetectorConstruction_edep();

  virtual void ConstructSDandField();

  virtual G4VPhysicalVolume* Construct();

private:

  G4double pa4;
  G4String parameter4;

  


  
  G4VPhysicalVolume* physWorld;
  G4VPhysicalVolume* physDetect;


  G4LogicalVolume* Aero1LW;
  G4LogicalVolume* Aero2LW;
  G4LogicalVolume* Aero3LW;
  G4LogicalVolume* Aero_SuppLW;
  G4LogicalVolume* HolderLW;
  G4LogicalVolume* BehindLW;
  G4LogicalVolume* Ae_sideLW;
  G4LogicalVolume* Behind_filmLW;
  G4LogicalVolume* BottomLW;
  G4LogicalVolume* ReflectLW;
  G4LogicalVolume* FilmLW;
  G4LogicalVolume* SideLW;
  G4LogicalVolume* SideLW1;
  G4LogicalVolume* Up_aeroLW;
  G4LogicalVolume* Extra_sideLW;
  G4LogicalVolume* Front_suppLW;
  G4LogicalVolume* Front_mylarLW;
  G4LogicalVolume* MPPCLW;

  G4LogicalVolume* UpstreamLW;
  G4LogicalVolume* DownstreamLW;



  std::vector<G4VisAttributes*> fVisAttributes;



  
};

#endif

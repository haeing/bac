#ifndef BACDetectorConstruction_hh
#define BACDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4VisAttributes;


class BACDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  BACDetectorConstruction(const G4String &parameter4);
  virtual ~BACDetectorConstruction();

  virtual void ConstructSDandField();
  virtual void BuildAerogelAbsLength(
    G4int tileID,
    G4double thickness,
    G4double* photonEnergy,
    G4double* absLength);

  virtual G4VPhysicalVolume* Construct();
  
  
private:
  
  //  extern const G4double kAerogelWavelength_nm[kAerogelNEntries];

  G4double pa4;
  G4String parameter4;

  G4double Aerox_real = 115.0 *mm;
  G4double Aeroy_real = 115.0 *mm;
  G4double Aeroz1 = 12.3*mm;
  G4double Aeroz2 = 12.3*mm;
  G4double Aeroz3 = 12.4*mm;
  G4double Aero_space = 1*mm;
  G4double Aeroz_real = Aeroz1+Aeroz2+Aero_space*4+Aeroz3;

  G4double Aeroy = 125*mm+10*mm;   //For the exact position information of parabola

  G4double aero_suppx = 112*mm;
  G4double aero_suppy = 111*mm;
  G4double aero_suppz = 1*mm;

  
  //G4int numRZ = 150;
  
  G4double win_thick = 1*mm;
  G4double mppc_thick = 1*mm;
  G4double re_thick = 1*mm;
  G4double mppc_place = 30*mm;


  


  
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

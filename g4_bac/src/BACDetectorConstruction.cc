#include "BACDetectorConstruction.hh"
#include "AeroSD.hh"
#include "MPPCSD.hh"


#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnionSolid.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trap.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "TString.h"
#include "TMath.h"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4RotationMatrix.hh"

static const G4int kAerogelNEntries = 61;
static const G4double wavelength_nm[kAerogelNEntries] = {
  800,790,780,770,760,750,740,730,720,710,
  700,690,680,670,660,650,640,630,620,610,
  600,590,580,570,560,550,540,530,520,510,
  500,490,480,470,460,450,440,430,420,410,
  400,390,380,370,360,350,340,330,320,310,
  300,290,280,270,260,250,240,230,220,210,200
};

static const G4double T_tile1[kAerogelNEntries] = {
  0.957,0.956,0.955,0.953,0.953,0.951,0.948,0.948,0.947,0.946,
  0.944,0.942,0.940,0.938,0.936,0.934,0.931,0.928,0.925,0.922,
  0.918,0.915,0.911,0.906,0.901,0.896,0.890,0.883,0.876,0.869,
  0.860,0.851,0.840,0.829,0.816,0.802,0.787,0.771,0.752,0.731,
  0.708,0.681,0.650,0.614,0.579,0.542,0.501,0.455,0.404,0.350,
  0.295,0.234,0.170,0.109,0.055,0.024,0.011,0.002,0.0,0.0,0.0
};

static const G4double T_tile2[kAerogelNEntries] = {
  0.961,0.960,0.959,0.958,0.957,0.954,0.952,0.952,0.951,0.950,
  0.948,0.946,0.944,0.942,0.940,0.937,0.934,0.931,0.928,0.925,
  0.922,0.918,0.914,0.909,0.904,0.899,0.893,0.886,0.879,0.871,
  0.863,0.853,0.843,0.831,0.819,0.804,0.789,0.772,0.753,0.732,
  0.710,0.683,0.651,0.617,0.582,0.543,0.501,0.452,0.400,0.348,
  0.293,0.232,0.170,0.110,0.056,0.025,0.011,0.002,0.0,0.0,0.0
};

static const G4double T_tile3[kAerogelNEntries] = {
  0.954,0.953,0.952,0.950,0.950,0.948,0.945,0.945,0.944,0.943,
  0.941,0.939,0.937,0.935,0.933,0.930,0.928,0.925,0.921,0.919,
  0.915,0.911,0.907,0.903,0.897,0.892,0.886,0.879,0.872,0.864,
  0.856,0.846,0.835,0.824,0.811,0.797,0.782,0.765,0.746,0.725,
  0.702,0.675,0.644,0.610,0.574,0.536,0.494,0.446,0.395,0.343,
  0.288,0.229,0.167,0.108,0.056,0.025,0.011,0.002,0.0,0.0,0.0
};

BACDetectorConstruction::BACDetectorConstruction(const G4String &par4)
  : G4VUserDetectorConstruction(),parameter4(par4)
{
}

BACDetectorConstruction::~BACDetectorConstruction()
{
  for (auto visAttributes: fVisAttributes){
    delete visAttributes;
  }
}

G4VPhysicalVolume* BACDetectorConstruction::Construct()
{
  
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  //Material---------------------------------------------------------------

  G4Material* world_mat = nist-> FindOrBuildMaterial("G4_AIR");
  //G4Material* Al = nist-> FindOrBuildMaterial("G4_Al");

  G4Element* C = new G4Element("Carbon","C",6.,12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
  G4Element* O = new G4Element("Oxygen","O",15.,15.9994*g/mole);
  G4Element* AlE = new G4Element("Aluminum","Al",13.,26.981538*g/mole);
  G4Element* B = new G4Element("Boron","B",5.,10.811*g/mole);
  G4Element* Na = new G4Element("Na","Na",11.,23.0*g/mole);
  G4Element* Si = new G4Element("Silicon","Si",14.,28.0844*g/mole);
  G4Element* Cl = new G4Element("Chlorine","Cl",17.,35.453*g/mole);
  G4Element* K = new G4Element("Potassium","K",19.,39.093*g/mole);
  G4Element* F = new G4Element("Fluorine","F",9.,18.998*g/mole);
  

  G4Material *Aerogel1 = new G4Material("Aerogel1",0.377*g/cm3,2);
  Aerogel1->AddElement(Si,1);
  Aerogel1->AddElement(O,2);

  G4Material *Aerogel2 = new G4Material("Aerogel2",0.379*g/cm3,2);
  Aerogel2->AddElement(Si,1);
  Aerogel2->AddElement(O,2);

  G4Material *Aerogel3 = new G4Material("Aerogel3",0.382*g/cm3,2);
  Aerogel3->AddElement(Si,1);
  Aerogel3->AddElement(O,2);

  G4Material* blacksheet = new G4Material("blacksheet",0.95*g/cm3,2);
  blacksheet->AddElement(C,1);
  blacksheet->AddElement(H,2);

  

  G4Material *Mylar = new G4Material("Mylar", 2.7*g/cm3, 1);
  Mylar->AddElement(AlE,1);
  


  G4Material *Film = new G4Material("Film", 1.39*g/cm3, 3);
  Film->AddElement(C, 5);
  Film->AddElement(H, 4);
  Film->AddElement(O, 2);

  G4Material *Polystyrene = new G4Material("Polystyrene",1*g/cm3,2);
  Polystyrene->AddElement(C,8);
  Polystyrene->AddElement(H,8);
  



  

  G4Material* Epoxi = new G4Material("Epoxi",1.1*g/cm3,4);
  Epoxi->AddElement(C, 21);
  Epoxi->AddElement(H, 25);
  Epoxi->AddElement(O,  5);
  Epoxi->AddElement(Cl, 1);


  G4Material* Acrylic = new G4Material("Acrylic", 1.19*g/cm3, 3);
  Acrylic->AddElement(C, 5);
  Acrylic->AddElement(H, 8);
  Acrylic->AddElement(O, 2);




  //Property--------------------------------------------------------------

  //air property---------------------------------------------------------

   
  G4MaterialPropertiesTable* prop_air = new G4MaterialPropertiesTable();
  G4double air_ep[] = {1.3*eV,7.*eV};
  G4double air_rindex[] = {1.0,1.0};
  prop_air->AddProperty("RINDEX",air_ep,air_rindex,2)->SetSpline(true);
  world_mat->SetMaterialPropertiesTable(prop_air);


  //Aerogel Tiles Property
  

  G4MaterialPropertiesTable* prop_aerogel1 = new G4MaterialPropertiesTable();
  G4double aero1_photonEnergy[kAerogelNEntries];
  G4double aero1_absLength[kAerogelNEntries];
 
  BuildAerogelAbsLength(1,Aeroz1,aero1_photonEnergy,aero1_absLength);
  G4double aero1_ep[] = {1.3*eV,7.*eV};
  G4double aero1_rindex[]={1.115,1.115};
  prop_aerogel1->AddProperty("RINDEX",aero1_ep,aero1_rindex,2)->SetSpline(true);
  prop_aerogel1->AddProperty("ABSLENGTH",aero1_photonEnergy,aero1_absLength,kAerogelNEntries)->SetSpline(true);
  
  Aerogel1->SetMaterialPropertiesTable(prop_aerogel1);

  G4MaterialPropertiesTable* prop_aerogel2 = new G4MaterialPropertiesTable();
  G4double aero2_photonEnergy[kAerogelNEntries];
  G4double aero2_absLength[kAerogelNEntries];
  
  BuildAerogelAbsLength(2,Aeroz2,aero2_photonEnergy,aero2_absLength);
  G4double aero2_ep[] = {1.3*eV,7.*eV};
  G4double aero2_rindex[]={1.115,1.115};
  prop_aerogel2->AddProperty("RINDEX",aero2_ep,aero2_rindex,2)->SetSpline(true);
  prop_aerogel2->AddProperty("ABSLENGTH",aero2_photonEnergy,aero2_absLength,kAerogelNEntries)->SetSpline(true);
  
  Aerogel2->SetMaterialPropertiesTable(prop_aerogel2);

  G4MaterialPropertiesTable* prop_aerogel3 = new G4MaterialPropertiesTable();
  G4double aero3_photonEnergy[kAerogelNEntries];
  G4double aero3_absLength[kAerogelNEntries];
  
  BuildAerogelAbsLength(3,Aeroz3,aero3_photonEnergy,aero3_absLength);
  G4double aero3_ep[] = {1.3*eV,7.*eV};
  G4double aero3_rindex[]={1.117,1.117};
  prop_aerogel3->AddProperty("RINDEX",aero3_ep,aero3_rindex,2)->SetSpline(true);
  prop_aerogel3->AddProperty("ABSLENGTH",aero3_photonEnergy,aero3_absLength,kAerogelNEntries)->SetSpline(true);
  
  Aerogel3->SetMaterialPropertiesTable(prop_aerogel3);

  //Black sheet property-------------------------------------------------
  G4MaterialPropertiesTable* prop_bs = new G4MaterialPropertiesTable();

  G4double bs_ep[] = {1.3*eV,7.*eV};
  G4double bs_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4double bs_rindex[]={1.6,1.6};
  
  prop_bs->AddProperty("RINDEX",bs_ep,bs_rindex,2);
  prop_bs->AddProperty("ABSLENGTH",bs_ep,bs_abs,2);

  //prop_bs->AddProperty("RINDEX",bs_ep,bs_rindex,false,true);
  //prop_bs->AddProperty("ABSLENGTH",bs_ep,bs_abs,false,true);
  blacksheet->SetMaterialPropertiesTable(prop_bs);




  //Checking Scintillation from Mylar
  /*
  G4double film_ep[] = {1.3*eV,7.*eV};
  G4double film_rindex[] = {1.63,1.70};
  G4double film_abs[] = {1*cm,1*cm};
  */

  //G4double film_ep[] = {1.24*eV,7.*eV};
  /*
  G4double film_ep1[] = {2.18*eV,2.21*eV,2.25*eV,2.28*eV,2.32*eV,2.36*eV,2.39*eV,2.41*eV,2.44*eV,2.46*eV,2.49*eV,2.51*eV,
    2.54*eV,2.55*eV,2.58*eV,2.61*eV,2.63*eV,2.64*eV,2.66*eV,2.68*eV,2.7*eV,2.72*eV,2.73*eV,2.75*eV,2.77*eV,2.79*eV,
    2.82*eV,2.85*eV,2.87*eV,2.89*eV,2.91*eV,2.93*eV,2.95*eV,2.96*eV,2.965*eV,2.97*eV,2.99*eV,3*eV,3.02*eV,3.03*eV,
    3.05*eV,3.06*eV,3.08*eV,3.09*eV,3.12*eV,3.15*eV,3.2*eV};
  G4double film_fast[]={0.00014325,0.0015758,0.00329528,0.00587386,0.0080231,0.01060168,0.01862478,0.02134662,
    0.02363911,0.0265042,0.02908321,0.03123203,0.03366778,0.0358166,0.03739239,0.03939838,0.04054442,0.04169045,
    0.04212022,0.0415472,0.04068767,0.03882536,0.03681938,0.03495707,0.03280783,0.03123203,0.02893995,0.02636095,
    0.02406888,0.02163313,0.01991407,0.01833827,0.01776483,0.01590252,0.01446998,0.01318069,0.01146121,0.01017192,
    0.00902588,0.00773659,0.00644688,0.00515759,0.00458457,0.00329528,0.00214882,0.00143255,0.00085953,0.00042976};

  const G4int numentries_scin1 = sizeof(film_ep1)/sizeof(G4double);
  assert (sizeof(film_ep1) == sizeof(film_fast));
  G4MaterialPropertiesTable* prop_film = new G4MaterialPropertiesTable();
  prop_film->AddProperty("RINDEX",film_ep,film_rindex,2);
  prop_film->AddProperty("ABSLENGTH",film_ep,film_abs,2);
  prop_film->AddProperty("FASTCOMPONENT",film_ep1,film_fast,numentries_scin1);
  prop_film->AddProperty("SLOWCOMPONENT",film_ep1,film_fast,numentries_scin1);
  prop_film->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);
  prop_film->AddConstProperty("RESOLUTIONSCALE",1.0);
  prop_film->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  prop_film->AddConstProperty("SLOWTIMECONSTANT",2.8*ns);
  prop_film->AddConstProperty("YIELDRATIO",0.8);
  Film->SetMaterialPropertiesTable(prop_film);
  */

  //Polystyrene property
  /*
  G4double poly_ep[] = {1.3*eV,7.*eV};
  G4double poly_rindex[] = {1.6,1.6};
  G4double poly_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4MaterialPropertiesTable* prop_poly = new G4MaterialPropertiesTable();
  prop_poly->AddProperty("RINDEX",poly_ep,poly_rindex,2);
  prop_poly->AddProperty("ABSLENGTH",poly_ep,poly_abs,2);
  prop_poly->AddProperty("FASTCOMPONENT",film_ep1,film_fast,numentries_scin1);
  prop_poly->AddProperty("SLOWCOMPONENT",film_ep1,film_fast,numentries_scin1);
  prop_poly->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);
  prop_poly->AddConstProperty("RESOLUTIONSCALE",1.0);
  prop_poly->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  prop_poly->AddConstProperty("SLOWTIMECONSTANT",2.8*ns);
  prop_poly->AddConstProperty("YIELDRATIO",0.8);
  Polystyrene->SetMaterialPropertiesTable(prop_poly);
  */
  

  
  //Mylar property
  G4double mylar_ep[]={1.4*eV,1.48*eV,1.52*eV,1.56*eV,1.6*eV,1.7*eV,1.8*eV,1.9*eV,2*eV,2.2*eV,2.4*eV,2.6*eV,2.8*eV,3*eV,3.4*eV,3.8*eV,4*eV,5*eV,6*eV,7*eV};

  G4double mylar_real[]={2.2802,2.6945,2.7668,2.7675,2.6154,2.1606,1.8301,1.5724,1.366,1.0728,0.8734,0.7278,0.6079,0.52135,0.39877,0.31474,0.28003,0.18137,0.12677,0.094236};
  G4double mylar_ima[]={8.1134,8.1878,8.2573,8.3866,8.4914,8.3565,8.0601,7.7354,7.4052,6.7839,6.2418,5.7781,5.3676,5.0008,4.3957,3.9165,3.7081,2.9029,2.3563,1.9519};

  G4double mylar_ep1[] = {1.24*eV,7.*eV};

  assert (sizeof(mylar_ep) == sizeof(mylar_real));
  assert (sizeof(mylar_ep) == sizeof(mylar_ima));
  
  const G4int numentries_mylar = sizeof(mylar_ep)/sizeof(G4double);
  G4double mylar_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4MaterialPropertiesTable* prop_mylar = new G4MaterialPropertiesTable();
  prop_mylar->AddProperty("REALRINDEX",mylar_ep,mylar_real,numentries_mylar)->SetSpline(true);
  prop_mylar->AddProperty("IMAGINARYRINDEX",mylar_ep,mylar_ima,numentries_mylar)->SetSpline(true);
  prop_mylar->AddProperty("ABSLENGTH",mylar_ep1,mylar_abs,2)->SetSpline(true);

  //prop_mylar->AddProperty("REALRINDEX",mylar_ep,mylar_real,false, true);
  //prop_mylar->AddProperty("IMAGINARYRINDEX",mylar_ep,mylar_ima,false,true);
  //prop_mylar->AddProperty("ABSLENGTH",mylar_ep1,mylar_abs,false, true);
  Mylar->SetMaterialPropertiesTable(prop_mylar);


  
  
  //MPPC property
  G4MaterialPropertiesTable* prop_mppc = new G4MaterialPropertiesTable();
  
  G4double mppc_rindex[]={1.,1.};
  G4double mppc_ep[] = {1.6*eV,7.*eV};
  G4double mppc_abs[] = {1.0*cm,1.0*cm};
  //G4double mppc_abs[] = {1.0*cm,1.0*cm};

  prop_mppc->AddProperty("RINDEX",mppc_ep,mppc_rindex,2)->SetSpline(true);
  prop_mppc->AddProperty("ABSLENGTH",mppc_ep,mppc_abs,2)->SetSpline(true);

  //prop_mppc->AddProperty("RINDEX",mppc_ep,mppc_rindex,false,true);
  //prop_mppc->AddProperty("ABSLENGTH",mppc_ep,mppc_abs,false,true);
  //Al->SetMaterialPropertiesTable(prop_mppc);

  //MPPC surface (Epoxi) property
  G4MaterialPropertiesTable* prop_epoxi = new G4MaterialPropertiesTable();
 
  G4double epoxi_rindex[]={1.5,1.5};
  G4double epoxi_ep[] = {1.3*eV,7.*eV};
  G4double epoxi_abs[] = {1.0*cm,1.0*cm};


  prop_epoxi->AddProperty("RINDEX",epoxi_ep,epoxi_rindex,2)->SetSpline(true);
  prop_epoxi->AddProperty("ABSLENGTH",epoxi_ep,epoxi_abs,2)->SetSpline(true);

  //prop_epoxi->AddProperty("RINDEX",epoxi_ep,epoxi_rindex,false,true);
  //prop_epoxi->AddProperty("ABSLENGTH",epoxi_ep,epoxi_abs,false,true);
  Epoxi->SetMaterialPropertiesTable(prop_epoxi);




  
  //Geometry---------------------------------------------------------------

  pa4 = std::stod(parameter4);

  //World------------------------------------------------------------------
  G4double world_size = 1*m;
  G4Box* solidWorld = new G4Box("World",world_size, world_size, world_size); 
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  physWorld = new G4PVPlacement(0,G4ThreeVector(), logicWorld, "World",0,false,0,checkOverlaps);


  //Aerogel
  G4Box* Aero1 = new G4Box("Aero1",Aerox_real/2,Aeroy_real/2,Aeroz1/2);
  G4Box* Aero2 = new G4Box("Aero2",Aerox_real/2,Aeroy_real/2,Aeroz2/2);
  G4Box* Aero3 = new G4Box("Aero3",Aerox_real/2,Aeroy_real/2,Aeroz3/2);

  Aero1LW = new G4LogicalVolume(Aero1,Aerogel1,"Aero1");
  Aero2LW = new G4LogicalVolume(Aero2,Aerogel2,"Aero2");
  Aero3LW = new G4LogicalVolume(Aero3,Aerogel3,"Aero3");


  //Reflector around the aerogel
  G4Box* Behind = new G4Box("Behind",Aerox_real/2,Aeroy_real/2,0.05*mm/2);
  BehindLW = new G4LogicalVolume(Behind,Mylar,"Behind");
  //Behind_filmLW = new G4LogicalVolume(Behind,Film,"Behind_film");  //also for check sincillation from the mylar

  
  if(pa4==2){
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,0*mm),Aero2LW,"Aero2",logicWorld,false,123,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,+Aeroz2/2+Aero_space+Aeroz3/2),Aero3LW,"Aero3",logicWorld,false,123,checkOverlaps);
  
    new G4PVPlacement(0,G4ThreeVector(0,0*mm,-Aeroz2/2-0.05*mm),BehindLW,"Behind",logicWorld,false,0,checkOverlaps);
  }

  else if(pa4==3){
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,-Aeroz2/2-Aero_space-Aeroz1/2),Aero1LW,"Aero1",logicWorld,false,123,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,0*mm),Aero2LW,"Aero2",logicWorld,false,123,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,+Aeroz2/2+Aero_space+Aeroz3/2),Aero3LW,"Aero3",logicWorld,false,123,checkOverlaps);
  
    new G4PVPlacement(0,G4ThreeVector(0,0*mm,-Aeroz2/2-Aero_space-Aeroz1-0.05*mm),BehindLW,"Behind",logicWorld,false,0,checkOverlaps);
  }
    

  


  //Aerogel Supporting
  G4Box *Aero_Suppcut = new G4Box("Aero_Suppcut",aero_suppx/2,aero_suppy/2,aero_suppz);
  G4Box *Aero_Suppcover = new G4Box("Aero_Suppcover",Aerox_real/2,Aeroy_real/2,aero_suppz/2);
  G4SubtractionSolid *Aero_Supp = new G4SubtractionSolid("Aero_Supp",Aero_Suppcover,Aero_Suppcut,0,G4ThreeVector());
  Aero_SuppLW = new G4LogicalVolume(Aero_Supp,blacksheet,"Aero_Supp");

  new G4PVPlacement(0,G4ThreeVector(0,0,Aeroz2/2+aero_suppz+Aeroz3+aero_suppz/2),Aero_SuppLW,"Aero_Supp",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,Aeroz2/2+aero_suppz/2),Aero_SuppLW,"Aero_Supp",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,-Aeroz2/2-aero_suppz/2),Aero_SuppLW,"Aero_Supp",logicWorld,false,0,checkOverlaps);
		    


  
  
  

  


  G4Box *Bottom = new G4Box("Bottom",Aerox_real/2,1*mm/2,Aeroz_real/2);
  BottomLW = new G4LogicalVolume(Bottom,Mylar,"Bottom");

  new G4PVPlacement(0,G4ThreeVector(0,-Aeroy_real/2-0.5*mm-0.05*mm,0),BottomLW,"Bottom",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,Aeroy_real/2+0.5*mm+0.05*mm,0),BottomLW,"Bottom",logicWorld,false,0,checkOverlaps);

  G4Box *Ae_side = new G4Box("Ae_side",1*mm/2,Aeroy_real/2,Aeroz_real/2);
  Ae_sideLW = new G4LogicalVolume(Ae_side,Mylar,"Ae_side");

  new G4PVPlacement(0,G4ThreeVector(-Aerox_real/2-1*mm/2,0,0),Ae_sideLW,"Ae_side",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(Aerox_real/2+1*mm/2,0,0),Ae_sideLW,"Ae_side",logicWorld,false,0,checkOverlaps);

  G4double front_suppx = 7*mm;
  G4double front_suppz = 3.9*mm;
  G4Box *Front_supp = new G4Box("Front_supp",front_suppx/2,Aeroy_real/2,front_suppz/2);
  G4Box *Front_mylar = new G4Box("Front_mylar",front_suppx/2,Aeroy_real/2,1*mm/2);

  Front_suppLW = new G4LogicalVolume(Front_supp,blacksheet,"Front_supp");
  Front_mylarLW = new G4LogicalVolume(Front_mylar,Mylar,"Front_mylar");

  new G4PVPlacement(0,G4ThreeVector(117*mm/2+front_suppx/2,0,Aeroz2/2+Aero_space+Aeroz3+Aero_space+front_suppz/2),Front_suppLW,"Front_supp",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-117*mm/2-front_suppx/2,0,Aeroz2/2+Aero_space+Aeroz3+Aero_space+front_suppz-0.5*mm),Front_suppLW,"Front_mylar",logicWorld,false,0,checkOverlaps);

  /*
  G4Box *Extra_side = new G4Box("Extra_side",7*mm/2,Aeroy_real/2,1*mm);
  Extra_sideLW = new G4LogicalVolume(Extra_side,Mylar,"Extra_side");
  new G4PVPlacement(0,G4ThreeVector(-131*mm/2+3.5*mm,0,Aeroz2/2+Aeroz3+Aero_space+5*mm),Extra_sideLW,"Extra_side",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(131*mm/2-3.5*mm,0,Aeroz2/2+Aeroz3+Aero_space+5*mm),Extra_sideLW,"Extra_side",logicWorld,false,0,checkOverlaps);
  */


  //Parabola
  G4int numRZ = 28;
  G4double x[numRZ];
  G4double y[numRZ];
  G4double x_out[numRZ];


  G4double p = 36;
  for(int i=0;i<numRZ;i++){

    x[i] = i*5;
    x_out[i] = (i*5)+0.05;
    
    y[i] = x[i]*x[i]/(4*p);
      
  }

  std::vector<G4TwoVector> poly(2*numRZ);
  for(int i=0;i<numRZ;i++){
    poly[i].set(x[i]*mm,-y[i]*mm);
    poly[2*numRZ-1-i].set(x_out[i]*mm,-y[i]*mm);
  }


  G4Box* local_x = new G4Box("local_x",100*mm, 0.1*mm, 0.1*mm);
  G4Box* local_y = new G4Box("local_y",0.1*mm, 100*mm, 0.1*mm);

  local_xLW = new G4LogicalVolume(local_x,Mylar,"local_x");
  local_yLW = new G4LogicalVolume(local_y,Mylar,"local_y");
 
    
  G4TwoVector offsetA(0.,0.), offsetB(0.,0.);
  G4double scaleA=1., scaleB=1.;
  G4ExtrudedSolid* Reflect = new G4ExtrudedSolid("Reflect",poly,131*mm/2,offsetA,scaleA, offsetB, scaleB);
  ReflectLW = new G4LogicalVolume(Reflect,Mylar,"Reflect");
  //ReflectLW = new G4LogicalVolume(Reflect,ESR,"Reflect");
  FilmLW = new G4LogicalVolume(Reflect,Film,"Film");
  G4RotationMatrix *rotY = new G4RotationMatrix();
  rotY->rotateY(+90*degree);
  rotY->rotateZ(+65*degree);

  G4double pary = Aeroy/2+5*cm;
  G4double parz = Aeroz_real/2+mppc_place+3.5*cm;

  //new G4PVPlacement(rotY,G4ThreeVector(0,pary,parz),FilmLW,"Film",logicWorld,false,0,checkOverlaps); //to check the scintillation of the mylar film
  new G4PVPlacement(rotY,G4ThreeVector(0,pary,parz+0.2*mm),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);
  //new G4PVPlacement(rotY,G4ThreeVector(0,pary,parz+0.2*mm),local_xLW,"local_x",logicWorld,false,0,checkOverlaps);
  //new G4PVPlacement(rotY,G4ThreeVector(0,pary,parz+0.2*mm),local_yLW,"local_y",logicWorld,false,0,checkOverlaps);

  //To check focus
  G4double focusR = 1.0*mm;
  G4double focusDz = 2.0*mm;
  
  auto solidFocus = new G4Tubs("FocusMarker", 0., focusR, focusDz/2, 0.*deg, 360.*deg);
  auto logicFocus = new G4LogicalVolume(solidFocus, world_mat, "FocusMarkerLV"); // material 아무거나 OK
  
  G4ThreeVector focusLocal(0.*mm, -p*mm, 0.*mm);
  G4ThreeVector focusWorld = (*rotY) * focusLocal + G4ThreeVector(0, pary, parz+0.2*mm);

  //new G4PVPlacement(nullptr, focusWorld, logicFocus, "FocusMarker", logicWorld, false, 0, checkOverlaps);
  

    



  //side reflector
  G4Tubs* LED = new G4Tubs("LED",0,3*mm,2*mm,0,2*TMath::Pi());
  G4Box* Side_full = new G4Box("Side",1*mm/2,14*cm,20*cm/2);
  G4RotationMatrix *rotLED = new G4RotationMatrix();
  rotLED->rotateY(+90*degree);
  G4SubtractionSolid *Side = new G4SubtractionSolid("Side",Side_full,LED,rotLED,G4ThreeVector(0,83.5*mm,75.4*mm));
  //SideLW = new G4LogicalVolume(Side,Mylar,"Side");
  SideLW = new G4LogicalVolume(Side_full,Mylar,"Side");
  SideLW1 = new G4LogicalVolume(Side_full,Mylar,"Side1");

  new G4PVPlacement(0,G4ThreeVector(-131*mm/2-1*mm/2,0,0),SideLW1,"Side1",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(131*mm/2+1*mm/2,0,0),SideLW,"Side",logicWorld,false,0,checkOverlaps);

      

  //reflector upper of aerogel
  G4Box* Up_aero = new G4Box("Up_aero",131*mm/2,12*mm/2,1*mm/2);
  Up_aeroLW = new G4LogicalVolume(Up_aero,Mylar,"Up_aero");
  //new G4PVPlacement(0,G4ThreeVector(0,67*mm,27.2*mm),Up_aeroLW,"Up_aero",logicWorld,false,0,checkOverlaps);
   
  //MPPC---------------------------------------------------------------------------

  G4double mppc_theta = 40*degree;
  G4RotationMatrix *rotM = new G4RotationMatrix();
  rotM->rotateX(90*degree+mppc_theta);

  
  //G4Box* solidMPPC = new G4Box("MPPCWorld",Aerox_real/2+5*mm/2,40*mm,mppc_thick/2);
  //G4Box* solidMPPC = new G4Box("MPPCWorld",131*mm/2,40*mm,mppc_thick/2);
  //G4Box* solidMPPC = new G4Box("MPPCWorld",131*mm/2,67*mm/2,mppc_thick/2);
  G4Box *solidFrameOuter = new G4Box("FrameOuterWorld",131*mm/2,120*mm/2,mppc_thick);
  G4Box* solidMPPC = new G4Box("MPPCWorld",131*mm/2,67*mm/2,mppc_thick);
  auto solidFrame = new G4SubtractionSolid("FrameWorld",solidFrameOuter,solidMPPC);
  mppcworld = new G4LogicalVolume(solidMPPC,Mylar,"MPPCWorld");
  frameworld = new G4LogicalVolume(solidFrame,blacksheet,"FrameWorld");
  //G4LogicalVolume* mppcworld = new G4LogicalVolume(solidMPPC,ESR,"MPPCWorld");
  new G4PVPlacement(rotM,G4ThreeVector(0*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta),Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)),mppcworld,"MPPCWorld",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(rotM,G4ThreeVector(0*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta),Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)),frameworld,"FrameWorld",logicWorld,false,0,checkOverlaps);
  
  G4Box* MPPC = new G4Box("MPPC",12*mm,12*mm,mppc_thick/2);
  MPPCLW = new G4LogicalVolume(MPPC,Epoxi,"MPPC");
  G4double ml = 30;
  //Four MPPCs

  for(int i=0;i<4;i++){
    new G4PVPlacement(0,G4ThreeVector(-(ml*1.5)*mm+(ml*i)*mm,0,0),MPPCLW,"MPPC",mppcworld,false,i+1,checkOverlaps);
  }
  

  G4Box *bottomframe = new G4Box("bottomframe",300*mm/2,1*mm,300*mm/2);
  bottomframeLW = new G4LogicalVolume(bottomframe,blacksheet,"bottomframe");
  new G4PVPlacement(0,G4ThreeVector(0*mm,-1*Aeroy_real/2-5*mm,0*mm),bottomframeLW,"bottomframe",logicWorld,0,checkOverlaps);

  //Check edep

  G4Box* dummy = new G4Box("dummy",15*cm,15*cm,1*mm);
  UpstreamLW = new G4LogicalVolume(dummy,world_mat,"upstream");
  DownstreamLW = new G4LogicalVolume(dummy,world_mat,"downstream");
  //new G4PVPlacement(0,G4ThreeVector(0,0,-Aeroz_real/2-0.30*mm-1*mm),UpstreamLW,"upstream",logicWorld,false,1000,checkOverlaps);
  //new G4PVPlacement(0,G4ThreeVector(0,0,10*cm),DownstreamLW,"downstream",logicWorld,false,2000,checkOverlaps);

  
  


  //visattributes------------------------------------------------

  auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  visAttributes -> SetVisibility(false);
  logicWorld->SetVisAttributes(visAttributes);
  //SideLW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  
  visAttributes = new G4VisAttributes(G4Color::Black()); 
  Aero1LW->SetVisAttributes(visAttributes);
  Aero2LW->SetVisAttributes(visAttributes);
  Aero3LW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(false);
  BehindLW->SetVisAttributes(visAttributes);
  Aero_SuppLW->SetVisAttributes(visAttributes);
  Ae_sideLW->SetVisAttributes(visAttributes);
  BottomLW->SetVisAttributes(visAttributes);
  Front_suppLW->SetVisAttributes(visAttributes);
  Front_mylarLW->SetVisAttributes(visAttributes);
  SideLW->SetVisAttributes(visAttributes);
  SideLW1->SetVisAttributes(visAttributes);
  Up_aeroLW->SetVisAttributes(visAttributes);
  frameworld->SetVisAttributes(visAttributes);
  mppcworld->SetVisAttributes(visAttributes);
  bottomframeLW->SetVisAttributes(visAttributes);
  
  fVisAttributes.push_back(visAttributes);

  

  visAttributes = new G4VisAttributes(G4Color::Black()); 
  MPPCLW->SetVisAttributes(visAttributes);
  local_xLW->SetVisAttributes(visAttributes);
  local_yLW->SetVisAttributes(visAttributes);
  
  fVisAttributes.push_back(visAttributes);

  logicFocus->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  

  visAttributes = new G4VisAttributes(G4Color::Black()); 
  ReflectLW->SetVisAttributes(visAttributes);
  
  BehindLW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  








  //Surface--------------------------------------------------------------
  //black sheet surface------------------
  
  G4OpticalSurface* surface_bs = new G4OpticalSurface("surface_bs");
  surface_bs->SetType(dielectric_dielectric);
  surface_bs->SetModel(unified);
  surface_bs->SetFinish(groundtyvekair);

  G4MaterialPropertiesTable* sp_bs = new G4MaterialPropertiesTable();

  G4double bs_effi[] = {0.0,0.0};
  G4double bs_reflec[]= {0.05,0.05};
  G4double bs_specularLobe[] = {0.3,0.3};
  G4double bs_specularSpike[] = {0.2,0.2};
  G4double bs_backScatter[] = {0,0};

  sp_bs->AddProperty("EFFICIENCY",bs_ep,bs_effi,2)->SetSpline(true);
  sp_bs->AddProperty("REFLECTIVITY",bs_ep,bs_reflec,2)->SetSpline(true);
  sp_bs->AddProperty("SPECULARLOBECONSTANT",bs_ep,bs_specularLobe,2)->SetSpline(true);
  sp_bs->AddProperty("SPECULARSPIKECONSTANT",bs_ep,bs_specularSpike,2)->SetSpline(true);
  sp_bs->AddProperty("BACKSCATTERCONSTANT",bs_ep,bs_backScatter,2)->SetSpline(true);

  //sp_bs->AddProperty("EFFICIENCY",bs_ep,bs_effi,false,true);
  //sp_bs->AddProperty("REFLECTIVITY",bs_ep,bs_reflec,false,true);
  //sp_bs->AddProperty("SPECULARLOBECONSTANT",bs_ep,bs_specularLobe,false,true);
  //sp_bs->AddProperty("SPECULARSPIKECONSTANT",bs_ep,bs_specularSpike,false,true);
  //sp_bs->AddProperty("BACKSCATTERCONSTANT",bs_ep,bs_backScatter,false,true);
  surface_bs->SetMaterialPropertiesTable(sp_bs);

  new G4LogicalSkinSurface("bs_surface",HolderLW,surface_bs);

  //mylar_al surface-------------------
  G4OpticalSurface* surface_mylar = new G4OpticalSurface("surface_mylar");
  surface_mylar->SetType(dielectric_metal);
  surface_mylar->SetFinish(polished);    
  surface_mylar->SetModel(unified);
  
  G4MaterialPropertiesTable* sp_mylar = new G4MaterialPropertiesTable();
  G4double mylar_reflec[] = {0.98,0.98};  //for metal, reflectivity is calculated using rindex, they use polarization, angle, energy
  //G4double mylar_effi[] = {0.9,0.9}; //not used
  //G4double mylar_effi[] = {1,1}; //not used
  //G4double mylar_effi[] = {0.01,0.01};
  G4double mylar_specularLobe[] = {0.0,0.0};
  G4double mylar_specularSpike[]={1.0,1.0};

  //G4double mylar_specularLobe[] = {0.3,0.3};
  //G4double mylar_specularSpike[]={0.7,0.7};
  G4double mylar_backScatter[] = {0,0};
  //sp_mylar->AddProperty("EFFICIENCY",mylar_ep1,mylar_effi,2);
  //sp_mylar->AddProperty("REFLECTIVITY",air_ep,mylar_reflec,2)->SetSpline(true);
  
  sp_mylar->AddProperty("SPECULARLOBECONSTANT",mylar_ep1,mylar_specularLobe,2)->SetSpline(true);
  sp_mylar->AddProperty("SPECULARSPIKECONSTANT",mylar_ep1,mylar_specularSpike,2)->SetSpline(true);
  sp_mylar->AddProperty("BACKSCATTERCONSTANT",mylar_ep1,mylar_backScatter,2)->SetSpline(true);

  //sp_mylar->AddProperty("SPECULARLOBECONSTANT",mylar_ep1,mylar_specularLobe,false,true);
  //sp_mylar->AddProperty("SPECULARSPIKECONSTANT",mylar_ep1,mylar_specularSpike,false,true);
  //sp_mylar->AddProperty("BACKSCATTERCONSTANT",mylar_ep1,mylar_backScatter,false,true);
  surface_mylar->SetMaterialPropertiesTable(sp_mylar);
  

  new G4LogicalSkinSurface("mylar_surface",ReflectLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",SideLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",SideLW1,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",BottomLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",BehindLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",mppcworld,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",Ae_sideLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",Front_mylarLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",Up_aeroLW,surface_mylar);



  return physWorld;
}		    


void BACDetectorConstruction::BuildAerogelAbsLength(
			   G4int tileID,
			   G4double thickness,
			   G4double* photonEnergy,
			   G4double* absLength
			   ){
  const G4double* T = nullptr;

  switch(tileID){
  case 1: T = T_tile1; break;
  case 2: T = T_tile2; break;
  case 3: T = T_tile3; break;
    
  default:
    G4Exception("BuildAerogelAbsLength",
		"InvalidTile",
		FatalException,
		"Unknown AerogelTileID");
  }

  for(int i=0;i<kAerogelNEntries;i++){
    photonEnergy[i] = (1240.0 / wavelength_nm[i]) * eV;

    if(T[i] > 0.0)
      absLength[i] = - thickness / std::log(T[i]);
    else
      absLength[i] = 0.001*mm; 
  }
}

void BACDetectorConstruction::ConstructSDandField()
{
  
  auto aeroSD = new AeroSD("aeroSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(aeroSD);
  //UpstreamLW->SetSensitiveDetector(aeroSD);
  //DownstreamLW->SetSensitiveDetector(aeroSD);
  Aero1LW->SetSensitiveDetector(aeroSD);
  Aero2LW->SetSensitiveDetector(aeroSD);
  Aero3LW->SetSensitiveDetector(aeroSD);
  

  


  auto mppcSD = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
  MPPCLW->SetSensitiveDetector(mppcSD);




}



  

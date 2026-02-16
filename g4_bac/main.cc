#include "BACDetectorConstruction.hh"
#include "BACDetectorConstruction_edep.hh"
#include "BACDetectorConstruction_ELPH.hh"
#include "BACPrimaryGeneratorAction_ELPH.hh"
#include "BACPrimaryGeneratorAction_E70.hh"
#include "BACPrimaryGeneratorAction_E72.hh"
#include "BACPrimaryGeneratorAction_KEK.hh"
#include "BACRunAction.hh"
#include "BACEventAction.hh"
#include "BACStackingAction.hh"
#include "BACAnalysisManager.hh"
#include "AeroSD.hh"
#include "AeroHit.hh"
#include "MPPCSD.hh"
#include "MPPCHit.hh"

//#ifdef G4Multithreded
//#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4OpticalPhysics.hh"
//#include "G4OpticalParameters.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"
//#include "G4RunManagerFactory.hh"

#include "G4VisExecutive.hh"
#include "G4ScoringManager.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  CLHEP::HepRandom::setTheSeed((unsigned)time(NULL));

  G4String histname;

  G4String par1_put;
  G4String par2_put;
  G4String par3_put;
  G4String par4_put;
  G4String par5_put;


  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {


    ui = new G4UIExecutive(argc, argv);
    histname = "geant4_test.root";
    par1_put ="0";  //x position
    par2_put = "0";  //y position
    par3_put = "735"; //momentum
    par4_put = "3"; //Number of aerogel
    par5_put = "4"; //0 : ELPH, 1 : J-PARC E70, 2 : J-PARC Edep check 3: J-PARC E72, 4: KEK test
  }



  
  

  else if(argc>=3){
    histname = argv[2];
    par1_put = argv[3];
    par2_put = argv[4];
    par3_put = argv[5];
    par4_put = argv[6];
    par5_put = argv[7];
  }

  G4int geo_version = stod(par4_put);
  

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;
  
  
  //    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Set mandatory initialization classes
  //
  // Detector construction
  //J-PARC
  
  if(par5_put=="1"||par5_put=="3" || par5_put == "4")
    runManager->SetUserInitialization(new BACDetectorConstruction(par4_put));
  else if(par5_put=="2")
    runManager->SetUserInitialization(new BACDetectorConstruction_edep(par4_put));
  //ELPH
  else if(par5_put=="0")
    runManager->SetUserInitialization(new BACDetectorConstruction_ELPH(par2_put));
  







  

  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  //G4VModularPhysicsList* physicsList = new QGSP_BERT;

  //delete for the edep check

  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  //auto opticalParams = G4OpticalParameters::Instance();
  //opticalParams->SetCerenkovTrackSecondariesFirst(true);
  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  physicsList->RegisterPhysics(opticalPhysics);
  


  runManager->SetUserInitialization(physicsList);

  BACAnalysisManager *anaMan = new BACAnalysisManager(histname);
  //BACPrimaryGeneratorAction *priGen = new BACPrimaryGeneratorAction();
  if(par5_put=="0"){
    BACPrimaryGeneratorAction_ELPH *priGen = new BACPrimaryGeneratorAction_ELPH(par1_put,par2_put);
    runManager->SetUserAction(priGen);
  }
  else if(par5_put=="1"){
    BACPrimaryGeneratorAction_E70 *priGen = new BACPrimaryGeneratorAction_E70();
    runManager->SetUserAction(priGen);
  }

  else if(par5_put=="3"){
    BACPrimaryGeneratorAction_E72 *priGen = new BACPrimaryGeneratorAction_E72();
    runManager->SetUserAction(priGen);
  }
  else if(par5_put=="4"){
    BACPrimaryGeneratorAction_KEK *priGen = new BACPrimaryGeneratorAction_KEK(par1_put,par2_put);
    runManager->SetUserAction(priGen);
  }
  BACRunAction *runAction = new BACRunAction(anaMan);
  BACEventAction *eventAction = new BACEventAction(anaMan);
  BACStackingAction *stackAction = new BACStackingAction();

  //runManager->SetUserAction(priGen);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(stackAction);
    
  // User action initialization
  //runManager->SetUserInitialization(new BACActionInitialization());
  runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);

  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;

  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !


  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "BACAnalysisManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "AeroHit.hh"
#include "MPPCHit.hh"



#include "Randomize.hh"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

#include <string>
#include <sstream>

extern int gCerenkovCounter;
extern double decay_check;

BACAnalysisManager::BACAnalysisManager(const G4String & histname)
  :outfile(histname), fActive_(true)
{}

BACAnalysisManager::~BACAnalysisManager()
{
  SaveFile();
}

void BACAnalysisManager::SaveFile(void) const
{
  if (fActive_)
    hfile->Write();
}

void BACAnalysisManager::Terminate(void) const
{
  if(fActive_)
    {
      hfile->Write();
      hfile->Close();
    }
}

void BACAnalysisManager::BeginOfRun(const G4Run*)
{

  

  G4SDManager* SDManager = G4SDManager::GetSDMpointer();

  hfile = new TFile(outfile, "RECREATE");
  tree = new TTree("tree","EvtGen tree");

  tree->Branch("event",&event, "event/I");
  tree->Branch("nEvt",&nEvt, "nEvt/I");
  //tree->Branch("evtid",evtid,"evtid[nEvt]/I");
  tree->Branch("evtpid",&evtpid,"evtpid/I");
  tree->Branch("evtposx",&evtposx,"evtposx/D");
  tree->Branch("evtposy",&evtposy,"evtposy/D");
  tree->Branch("evtposz",&evtposz,"evtposz/D");
  tree->Branch("evtmomx",&evtmomx,"evtmomx/D");
  tree->Branch("evtmomy",&evtmomy,"evtmomy/D");
  tree->Branch("evtmomz",&evtmomz,"evtmomz/D");
  
  //tree->Branch("evtpid",evtpid,"evtpid[nEvt]/I");
  tree->Branch("evtnumce",&evtnumce,"evtnumce/I");


  //Aerogel
  /*
  tree->Branch("nhAero",&nhAero,"nhAero/I");
  //tree->Branch("aeroangle",aeroangle,"aeroangle[nhAero]/D");
  
  tree->Branch("aeropid",aeropid,"aeropid[nhAero]/I");
  //tree->Branch("aeroposx",aeroposx,"aeroposx[nhAero]/D");
  //tree->Branch("aeroposy",aeroposy,"aeroposy[nhAero]/D");
  //tree->Branch("aeroposz",aeroposz,"aeroposz[nhAero]/D");
  tree->Branch("aeromomx",aeromomx,"aeromomx[nhAero]/D");
  tree->Branch("aeromomy",aeromomy,"aeromomy[nhAero]/D");
  tree->Branch("aeromomz",aeromomz,"aeromomz[nhAero]/D");
  tree->Branch("aeroedep",aeroedep,"aeroedep[nhAero]/D");
  tree->Branch("aeronum",aeronum,"aeronum[nhAero]/I");
  //tree->Branch("aerotime",aerotime,"aerotime[nhAero]/D");

  */

  
  
  //MPPC
  tree->Branch("nhMppc",&nhMppc,"nhMppc/I");
  //tree->Branch("mppcmulti"&mppcmulti,"mppcmulti/I");
  tree->Branch("mppcpid",mppcpid,"mppcpid[nhMppc]/I");
  //tree->Branch("mppcposx",mppcposx,"mppcposx[nhMppc]/D");
  //tree->Branch("mppcposy",mppcposy,"mppcposy[nhMppc]/D");
  //tree->Branch("mppcposz",mppcposz,"mppcposz[nhMppc]/D");
  //tree->Branch("mppctime",mppctime,"mppctime[nhMppc]/D");
  //tree->Branch("mppcsignal",mppcsignal,"mppcsignal/D");
  tree->Branch("mppcwavelength",mppcwavelength,"mppcwavelength[nhMppc]/D");
  tree->Branch("mppcnum",mppcnum,"mppcnum[nhMppc]/I");




  event = 0;
  nEvt = 0;


}

void BACAnalysisManager::EndOfRun(const G4Run*)
{
  tree->Write();
  hfile->Write();
  hfile->Close();
}

void BACAnalysisManager::BeginOfEvent(const G4Event* anEvent)
{

}

void BACAnalysisManager::EndOfEvent(const G4Event* anEvent)
{

  G4HCofThisEvent* HCTE = anEvent->GetHCofThisEvent();
  if(!HCTE) return;
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhmppc = 0;
  G4int nhaero = 0;
  G4int multiplicity = 0;

  G4int pdg = anEvent->GetPrimaryVertex(0)->GetPrimary(0)->GetPDGcode();
  evtpid = pdg;
  
  G4ThreeVector pripos = anEvent->GetPrimaryVertex(0)->GetPosition();
  G4ThreeVector primom = anEvent->GetPrimaryVertex(0)->GetPrimary(0)->GetMomentum();
  
  G4double prix = pripos.x();
  G4double priy = pripos.y();
  G4double priz = pripos.z();

  G4double pripx = primom.x();
  G4double pripy = primom.y();
  G4double pripz = primom.z();
  
  evtposx = prix;
  evtposy = priy;
  evtposz = priz;

  evtmomx = pripx;
  evtmomy = pripy;
  evtmomz = pripz;

  evtnumce = gCerenkovCounter;
  

  
  MPPCHitsCollection *MPPCHC = 0;
  G4int ColIdMPPC = SDMan->GetCollectionID("MppcCollection");
  if(ColIdMPPC>=0)
    {
      MPPCHC=dynamic_cast<MPPCHitsCollection *>(HCTE->GetHC( ColIdMPPC ));
      if(MPPCHC)
	{
	  nhmppc = MPPCHC->entries();
	}
    }

  for(int i=0;i<nhmppc;i++)
    {
      MPPCHit* aHit = (*MPPCHC)[i];
      mppctime[i] = aHit->GetTOF();
      mppcposx[i] = aHit->GetPosition().x();
      mppcposy[i] = aHit->GetPosition().y();
      //if(aHit->GetWorldPosition().x()<=0)mppcposz[i]=-1.0;
      //else if(aHit->GetWorldPosition().x()>0)mppcposz[i]=1.0;
      mppcposz[i] = aHit->GetPosition().z();
      //mppcposz[i] = aHit->GetWorldPosition().y();
      mppcpid[i] = aHit->GetParticleID();
      mppcwavelength[i] = aHit->GetWavelength();
      mppcnum[i] = aHit-> GetCopyNum();


      
    }
  nhMppc = nhmppc;
  std::cout<<nhMppc<<std::endl;
  //mppcsignal = nhmppc+

  

  

  AeroHitsCollection *AEROHC = 0;
  G4int ColIdAERO = SDMan->GetCollectionID("AeroCollection");
  if(ColIdAERO>=0)
    {
      AEROHC=dynamic_cast<AeroHitsCollection *>(HCTE->GetHC(ColIdAERO));
      if(AEROHC)
	{
	  nhaero = AEROHC->entries();
	}
    }

  for(int i=0;i<nhaero;i++)
    {
      AeroHit* aHit = (*AEROHC)[i];
      aeroangle[i] = aHit->GetAngle();
      aerotime[i] = aHit->GetTOF();
      aeroposx[i] = aHit->GetPosition().x();
      aeroposy[i] = aHit->GetPosition().y();
      aeroposz[i] = aHit->GetPosition().z();
      aeromomx[i] = aHit->GetMomentum().x();
      aeromomy[i] = aHit->GetMomentum().y();
      aeromomz[i] = aHit->GetMomentum().z();
      aeropid[i] = aHit->GetParticleID();
      aeronum[i] = aHit-> GetCopyNum();
      aeroedep[i] = aHit->GetEdep();

    }
  nhAero = nhaero;
  //std::cout<<decay_check<<std::endl;
  //if(pdg==decay_check){
    tree->Fill();
    event++;
    // }

  nEvt=0;
  nhaero = 0;
  nhmppc = 0;


  
}

void BACAnalysisManager::SetEvtGen(int j,int partnum)
{
  if(nEvt<j) nEvt=j;
  //evtid[j-1]=j;
  //evtpid[j-1]=partnum;
}


  
      




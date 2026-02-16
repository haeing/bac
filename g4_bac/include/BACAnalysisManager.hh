#ifndef BACAnalysisManager_h
#define BACAnalysisManager_h

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <cmath>
#include <vector>

class G4Run;
class G4Event;
class TFile;
class TTree;
class EvtVector4R;

const int num_evtgen = 1000000;
const int num_mppchit = 1000000;
const int num_aerohit = 1000000;

class BACAnalysisManager
{
public:
  void BeginOfRun(const G4Run*);
  void EndOfRun(const G4Run*);
  void BeginOfEvent(const G4Event*);
  void EndOfEvent(const G4Event*);

private:
  G4String outfile;
  TFile *hfile;
  TTree *tree;

public:
  BACAnalysisManager(const G4String &histname);
  virtual ~BACAnalysisManager();
  BACAnalysisManager(const BACAnalysisManager&);
  BACAnalysisManager& operator=(const BACAnalysisManager&);

  void SetEvtGen(int j, int partnum);
  void Terminate(void) const;
  void SaveFile(void) const;

private:
  G4bool fActive_;
  int event;

  int nEvt;
  int evtid[num_evtgen];
  //int evtpid[num_evtgen];
  int evtpid;
  double evtposx;
  double evtposy;
  double evtposz;

  double evtmomx;
  double evtmomy;
  double evtmomz;
  
  int evtnumce;
  

  int nhMppc;
  //int mppcmulti;
  int mppcpid[num_mppchit];
  double mppctime[num_mppchit];
  double mppcposx[num_mppchit];
  double mppcposy[num_mppchit];
  double mppcposz[num_mppchit];
  double mppcwavelength[num_mppchit];
  int mppcnum[num_mppchit];

  int nhAero;
  double aeroangle[num_aerohit];
  int aeropid[num_aerohit];
  double aerotime[num_aerohit];
  double aeroposx[num_aerohit];
  double aeroposy[num_aerohit];
  double aeroposz[num_aerohit];
  double aeromomx[num_aerohit];
  double aeromomy[num_aerohit];
  double aeromomz[num_aerohit];
  double aeroedep[num_aerohit];
  int aeronum[num_aerohit];
  



  
};

#endif

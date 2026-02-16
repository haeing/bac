#ifndef AeroHit_hh
#define AeroHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"

class G4AttDef;
class G4AttValue;

class AeroHit : public G4VHit
{
private:
  G4ThreeVector xyz;
  G4ThreeVector momentum;
  G4int particleID;
  G4double tof;
  G4double angleget;
  G4int copynumber;
  G4double energyloss;

  
public:
  AeroHit();
  AeroHit(G4ThreeVector& axyz, G4ThreeVector& mxyz, G4double t);
  AeroHit(G4ThreeVector& axyz, G4ThreeVector& mxyz, G4double t, G4int pid, G4double angle, G4int cn,G4double edep);
  virtual ~AeroHit();

  AeroHit(const AeroHit& right);
  const AeroHit& operator=(const AeroHit &right);
  //G4bool operator==(const AeroHit &right) const;

  void *operator new(size_t);
  void operator delete(void *aHit);

  //virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  //virtual std::vector<G4AttValue>* CreateAttValues() const;

  const G4ThreeVector& GetPosition() const {return xyz;}
  const G4ThreeVector& GetMomentum() const {return momentum;}
  G4double GetTOF() const {return tof;}
  
  G4int GetParticleID() const {return particleID;}
  G4double GetAngle() const {return angleget;}

  G4int GetCopyNum() const {return copynumber;}
  G4double GetEdep() const{return energyloss;}
  //G4int GetNum() const {return count_ce;}

  //void SetTime(G4double dt) {fTime = dt;}


  //void IncPhotonCount() {fPhotons++;}
  //G4double GetPhotonCount() {return fPhotons;}
  //void ClearPhotonCount() {fPhotons=0.0;}

  //void SetPos(G4ThreeVector xyz) {fPos = xyz;}
  //G4ThreeVector GetPos() const {return fPos;}
  
  
};

inline AeroHit::AeroHit(const AeroHit& right)
  : G4VHit()
{
  xyz = right.xyz;
  momentum = right.momentum;
  tof = right.tof;
  particleID = right.particleID;
  angleget = right.angleget;
  copynumber = right.copynumber;
  energyloss = right.energyloss;
}

inline const AeroHit& AeroHit::operator=
(const AeroHit& right)
{
  xyz = right.xyz;
  momentum = right.momentum;
  tof = right.tof;
  particleID = right.particleID;
  angleget = right.angleget;
  copynumber = right.copynumber;
  energyloss = right.energyloss;

  return *this;
}

//using AeroHitsCollection = G4THitsCollection<AeroHit>;
typedef G4THitsCollection<AeroHit> AeroHitsCollection;
//extern G4ThreadLocal G4Allocator<AeroHit>* AeroHitAllocator;
extern G4Allocator<AeroHit> AeroHitAllocator;

inline void* AeroHit::operator new(size_t)
{
  //return (void*)AeroHitAllocator-> MallocSingle();
  return static_cast<void*>(AeroHitAllocator.MallocSingle());

}

inline void AeroHit::operator delete(void* aHit)
{
  //AeroHitAllocator->FreeSingle((AeroHit*) aHit);
  AeroHitAllocator.FreeSingle(static_cast<AeroHit*> (aHit));
}

#endif
  

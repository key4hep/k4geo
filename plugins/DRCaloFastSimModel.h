#include "G4OpBoundaryProcess.hh"
#include "G4GenericMessenger.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpAbsorption.hh"
#include "G4OpWLS.hh"
#include "G4Material.hh"

struct FastFiberData
{
public:
  FastFiberData(G4int, G4double, G4double, G4double, G4ThreeVector, G4ThreeVector, G4ThreeVector, G4int status = G4OpBoundaryProcessStatus::Undefined);
  ~FastFiberData() {}

  void reset();

  G4double GetAbsorptionNILL() { return mOpAbsorptionNumIntLenLeft; }
  void SetAbsorptionNILL(G4double in) { mOpAbsorptionNumIntLenLeft = in; }

  G4double GetWLSNILL() { return mOpWLSNumIntLenLeft; }
  void SetWLSNILL(G4double in) { mOpWLSNumIntLenLeft = in; }

  G4int GetOpBoundaryStatus() { return mOpBoundaryStatus; }
  void SetOpBoundaryStatus(G4int in) { mOpBoundaryStatus = in; }

  G4double GetStepLengthInterval() { return mStepLengthInterval; }
  void AddStepLengthInterval(G4double in) { mStepLengthInterval += in; }

  G4bool checkRepetitive(const FastFiberData, G4bool checkInterval = true);

  G4int trackID;
  G4double kineticEnergy;
  G4double globalTime;
  G4double pathLength;
  G4ThreeVector globalPosition;
  G4ThreeVector momentumDirection;
  G4ThreeVector polarization;

private:
  G4int mOpBoundaryStatus;
  G4double mOpAbsorptionNumIntLenLeft;
  G4double mOpWLSNumIntLenLeft;
  G4double mStepLengthInterval;
};

FastFiberData::FastFiberData(G4int id, G4double en, G4double globTime, G4double path, G4ThreeVector pos, G4ThreeVector mom, G4ThreeVector pol, G4int status)
{
  trackID = id;
  kineticEnergy = en;
  globalTime = globTime;
  pathLength = path;
  globalPosition = pos;
  momentumDirection = mom;
  polarization = pol;
  mOpBoundaryStatus = status;
  mOpAbsorptionNumIntLenLeft = DBL_MAX;
  mOpWLSNumIntLenLeft = DBL_MAX;
  mStepLengthInterval = 0.;
}

G4bool FastFiberData::checkRepetitive(const FastFiberData theData, G4bool checkInterval)
{
  if (this->trackID != theData.trackID)
    return false;
  if (this->mOpBoundaryStatus != theData.mOpBoundaryStatus)
    return false;
  if (checkInterval && std::abs(this->mStepLengthInterval - theData.mStepLengthInterval) > G4GeometryTolerance::GetInstance()->GetSurfaceTolerance())
    return false;

  return true;
}

void FastFiberData::reset()
{
  this->mOpBoundaryStatus = G4OpBoundaryProcessStatus::Undefined;
  this->mOpAbsorptionNumIntLenLeft = DBL_MAX;
  this->mOpWLSNumIntLenLeft = DBL_MAX;
  this->mStepLengthInterval = 0.;
}
#ifndef DRCaloFastSimModel_h
#define DRCaloFastSimModel_h

#include "G4GenericMessenger.hh"
#include "G4Material.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"
#include <G4GeometryTolerance.hh>
#include <G4OpProcessSubType.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4ProcessManager.hh>
#include <G4Tubs.hh>

struct FastFiberData {
public:
  FastFiberData(G4int id, G4double en, G4double globTime, G4double path, G4ThreeVector pos, G4ThreeVector mom,
                G4ThreeVector pol, G4int status = G4OpBoundaryProcessStatus::Undefined) {
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

  void reset() {
    this->mOpBoundaryStatus = G4OpBoundaryProcessStatus::Undefined;
    this->mOpAbsorptionNumIntLenLeft = DBL_MAX;
    this->mOpWLSNumIntLenLeft = DBL_MAX;
    this->mStepLengthInterval = 0.;
  }

  G4double GetAbsorptionNILL() { return mOpAbsorptionNumIntLenLeft; }
  void SetAbsorptionNILL(G4double in) { mOpAbsorptionNumIntLenLeft = in; }

  G4double GetWLSNILL() { return mOpWLSNumIntLenLeft; }
  void SetWLSNILL(G4double in) { mOpWLSNumIntLenLeft = in; }

  G4int GetOpBoundaryStatus() { return mOpBoundaryStatus; }
  void SetOpBoundaryStatus(G4int in) { mOpBoundaryStatus = in; }

  G4double GetStepLengthInterval() { return mStepLengthInterval; }
  void AddStepLengthInterval(G4double in) { mStepLengthInterval += in; }

  G4bool checkRepetitive(const FastFiberData theData, G4bool checkInterval = true) {
    if (this->trackID != theData.trackID)
      return false;
    if (this->mOpBoundaryStatus != theData.mOpBoundaryStatus)
      return false;
    if (checkInterval && std::abs(this->mStepLengthInterval - theData.mStepLengthInterval) >
                             G4GeometryTolerance::GetInstance()->GetSurfaceTolerance())
      return false;

    return true;
  }

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

class DRCFiberModel {
public:
  DRCFiberModel() = default;
  ~DRCFiberModel() = default;

  G4OpBoundaryProcess* pOpBoundaryProc = nullptr;
  G4OpAbsorption* pOpAbsorption = nullptr;
  G4OpWLS* pOpWLS = nullptr;
  G4bool fProcAssigned = false;

  FastFiberData mDataPrevious = FastFiberData(0, 0., 0., 0., G4ThreeVector(0), G4ThreeVector(0), G4ThreeVector(0));
  FastFiberData mDataCurrent = FastFiberData(0, 0., 0., 0., G4ThreeVector(0), G4ThreeVector(0), G4ThreeVector(0));

  G4int fSafety = 1;
  G4double mNtransport = 0.;
  G4double mTransportUnit = 0.;
  G4ThreeVector mFiberPos = G4ThreeVector(0);
  G4ThreeVector mFiberAxis = G4ThreeVector(0);
  G4bool fKill = false;
  G4bool fTransported = false;
  G4bool fSwitch = true;
  G4int fVerbose = 0;

  G4bool checkTotalInternalReflection(const G4Track* track) {
    if (!fProcAssigned)
      setPostStepProc(track); // locate OpBoundaryProcess only once

    G4int theStatus = pOpBoundaryProc->GetStatus();

    if (fVerbose > 1) {
      std::cout << "DRCFiberModel::checkTotalInternalReflection | TrackID = " << std::setw(4) << track->GetTrackID();
      std::cout << " | G4OpBoundaryProcessStatus = " << std::setw(2) << theStatus;
      std::cout << " | Track status = " << track->GetTrackStatus();
      std::cout << " | StepLength = " << std::setw(9) << track->GetStepLength() << std::endl;
    }

    if (track->GetTrackStatus() == fStopButAlive || track->GetTrackStatus() == fStopAndKill)
      return false;

    // accumulate step length
    mDataCurrent.AddStepLengthInterval(track->GetStepLength());

    // skip exceptional iteration with FresnelReflection
    if (theStatus == G4OpBoundaryProcessStatus::FresnelReflection)
      mDataCurrent.SetOpBoundaryStatus(theStatus);

    // some cases have a status StepTooSmall when the reflection happens
    // between the boundary of cladding & outer volume (outside->cladding)
    // since the outer volume is not a G4Region
    if (theStatus == G4OpBoundaryProcessStatus::TotalInternalReflection ||
        theStatus == G4OpBoundaryProcessStatus::StepTooSmall) {
      if (theStatus != G4OpBoundaryProcessStatus::TotalInternalReflection) {
        // skip StepTooSmall if the track already has TotalInternalReflection
        if (mDataCurrent.GetOpBoundaryStatus() == G4OpBoundaryProcessStatus::TotalInternalReflection)
          return false;
        if (mDataPrevious.GetOpBoundaryStatus() == G4OpBoundaryProcessStatus::TotalInternalReflection)
          return false;
      }

      G4int trackID = track->GetTrackID();
      G4double kineticEnergy = track->GetKineticEnergy();
      G4double globalTime = track->GetGlobalTime();
      G4double pathLength = track->GetStepLength();
      G4ThreeVector globalPosition = track->GetPosition();
      G4ThreeVector momentumDirection = track->GetMomentumDirection();
      G4ThreeVector polarization = track->GetPolarization();

      auto candidate = FastFiberData(trackID, kineticEnergy, globalTime, pathLength, globalPosition, momentumDirection,
                                     polarization, theStatus);

      if (pOpAbsorption != nullptr)
        candidate.SetAbsorptionNILL(pOpAbsorption->GetNumberOfInteractionLengthLeft());
      if (pOpWLS != nullptr)
        candidate.SetWLSNILL(pOpWLS->GetNumberOfInteractionLengthLeft());

      G4bool repetitive = false;

      if (candidate.checkRepetitive(mDataCurrent, false) && mDataCurrent.checkRepetitive(mDataPrevious))
        repetitive = true;

      mDataPrevious = mDataCurrent;
      mDataCurrent = candidate;

      return repetitive;
    }

    return false;
  } // checkTotalInternalReflection

  G4bool checkAbsorption(const G4double prevNILL, const G4double currentNILL) {
    if (prevNILL < 0. || currentNILL < 0.)
      return false; // the number of interaction length left has to be reset
    if (prevNILL == currentNILL)
      return false; // no absorption
    if (prevNILL == DBL_MAX || currentNILL == DBL_MAX)
      return false; // NILL is re-initialized

    G4double deltaNILL = prevNILL - currentNILL;

    if (currentNILL - deltaNILL * (mNtransport + fSafety) < 0.)
      return true; // absorbed before reaching fiber end

    return false;
  } // checkAbsorption

  G4bool checkNILL() {
    if (!fTransported)
      return true; // do nothing if the track is not already transported

    G4double wlsNILL = DBL_MAX;
    G4double absorptionNILL = DBL_MAX;

    if (pOpWLS != nullptr) {
      wlsNILL = pOpWLS->GetNumberOfInteractionLengthLeft();
      if (mDataPrevious.GetWLSNILL() == DBL_MAX || mDataCurrent.GetWLSNILL() == DBL_MAX)
        return true; // NILL is re-initialized
    }

    if (pOpAbsorption != nullptr) {
      absorptionNILL = pOpAbsorption->GetNumberOfInteractionLengthLeft();
      if (mDataPrevious.GetAbsorptionNILL() == DBL_MAX || mDataCurrent.GetAbsorptionNILL() == DBL_MAX)
        return true; // NILL is re-initialized
    }

    if (wlsNILL < 0. || absorptionNILL < 0.)
      return true; // let GEANT4 to reset them

    G4double deltaWlsNILL = mDataPrevious.GetWLSNILL() - mDataCurrent.GetWLSNILL();
    G4double deltaAbsorptionNILL = mDataPrevious.GetAbsorptionNILL() - mDataCurrent.GetAbsorptionNILL();

    G4double finalWlsNILL = wlsNILL - deltaWlsNILL * fSafety;
    G4double finalAbsorptionNILL = absorptionNILL - deltaAbsorptionNILL * fSafety;

    // prevent double counting of the probability of getting absorbed
    // (which already estimated before transportation)
    // reset NILL again
    if (finalWlsNILL < 0. || finalAbsorptionNILL < 0.)
      return false;

    return true;
  } // checkNILL

  bool check_trigger(const G4Track* track) {
    if (fTransported) {
      // different propagation direction (e.g. mirror)
      if (mFiberAxis.dot(track->GetMomentumDirection()) * mFiberAxis.dot(mDataCurrent.momentumDirection) < 0)
        reset();

      // track is already transported and did pass NILL check, nothing to do
      return false;
    }

    if (!checkTotalInternalReflection(track))
      return false; // nothing to do if the track has no repetitive total internal reflection

    auto theTouchable = track->GetTouchableHandle();
    auto solid = theTouchable->GetSolid();

    if (fVerbose > 0)
      print(); // at this point, the track should have passed all prerequisites before entering computationally heavy
               // operations

    if (solid->GetEntityType() != "G4Tubs")
      return false; // only works for G4Tubs at the moment

    G4Tubs* tubs = static_cast<G4Tubs*>(solid);
    G4double fiberLen = 2. * tubs->GetZHalfLength();

    mFiberPos = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0., 0., 0.));
    mFiberAxis = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformAxis(G4ThreeVector(0., 0., 1.));

    auto delta = mDataCurrent.globalPosition - mDataPrevious.globalPosition;
    mTransportUnit = delta.dot(mFiberAxis);

    // estimate the number of expected total internal reflections before reaching fiber end
    auto fiberEnd =
        (mTransportUnit > 0.) ? mFiberPos + mFiberAxis * fiberLen / 2. : mFiberPos - mFiberAxis * fiberLen / 2.;
    auto toEnd = fiberEnd - track->GetPosition();
    G4double toEndAxis = toEnd.dot(mFiberAxis);
    G4double maxTransport = std::floor(toEndAxis / mTransportUnit);
    mNtransport = maxTransport - fSafety;

    if (mNtransport < 1.)
      return false; // require at least n = fSafety of total internal reflections at the end

    if (checkAbsorption(mDataPrevious.GetWLSNILL(), mDataCurrent.GetWLSNILL()))
      return false; // do nothing if WLS happens before reaching fiber end
    if (checkAbsorption(mDataPrevious.GetAbsorptionNILL(), mDataCurrent.GetAbsorptionNILL()))
      fKill = true; // absorbed before reaching fiber end

    return true;
  }

  void setPostStepProc(const G4Track* track) {
    G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
    auto postStepProcessVector = pm->GetPostStepProcessVector();

    for (unsigned int np = 0; np < postStepProcessVector->entries(); np++) {
      auto theProcess = (*postStepProcessVector)[np];

      auto theType = theProcess->GetProcessType();

      if (theType != fOptical)
        continue;

      if (theProcess->GetProcessSubType() == G4OpProcessSubType::fOpBoundary)
        pOpBoundaryProc = dynamic_cast<G4OpBoundaryProcess*>(theProcess);
      else if (theProcess->GetProcessSubType() == G4OpProcessSubType::fOpAbsorption)
        pOpAbsorption = dynamic_cast<G4OpAbsorption*>(theProcess);
      else if (theProcess->GetProcessSubType() == G4OpProcessSubType::fOpWLS)
        pOpWLS = dynamic_cast<G4OpWLS*>(theProcess);
    }

    fProcAssigned = true;

    return;
  } // setPostStepProc

  void reset() {
    mNtransport = 0.;
    mTransportUnit = 0.;
    mFiberPos = G4ThreeVector(0);
    mFiberAxis = G4ThreeVector(0);
    fKill = false;
    fTransported = false;
    mDataCurrent.reset();
    mDataPrevious.reset();
  }

  void print() {
    if (fVerbose > 1) {
      G4cout << G4endl;

      G4cout << "mDataPrevious.trackID = " << mDataPrevious.trackID;
      G4cout << " | .mOpBoundaryStatus = " << std::setw(4) << mDataPrevious.GetOpBoundaryStatus();
      G4cout << " | .mStepLengthInterval = " << mDataPrevious.GetStepLengthInterval() << G4endl;

      if (fVerbose > 2) {
        G4cout << "  | globalPosition    = (" << std::setw(9) << mDataPrevious.globalPosition.x();
        G4cout << "," << std::setw(9) << mDataPrevious.globalPosition.y();
        G4cout << "," << std::setw(9) << mDataPrevious.globalPosition.z() << ")" << G4endl;

        G4cout << "  | momentumDirection = (" << std::setw(9) << mDataPrevious.momentumDirection.x();
        G4cout << "," << std::setw(9) << mDataPrevious.momentumDirection.y();
        G4cout << "," << std::setw(9) << mDataPrevious.momentumDirection.z() << ")" << G4endl;

        G4cout << "  | polarization      = (" << std::setw(9) << mDataPrevious.polarization.x();
        G4cout << "," << std::setw(9) << mDataPrevious.polarization.y();
        G4cout << "," << std::setw(9) << mDataPrevious.polarization.z() << ")" << G4endl;

        G4cout << "  | globalTime        =  " << std::setw(9) << mDataPrevious.globalTime << G4endl;
        G4cout << "  | WLSNILL           =  " << std::setw(9) << mDataPrevious.GetWLSNILL() << G4endl;
        G4cout << "  | AbsorptionNILL    =  " << std::setw(9) << mDataPrevious.GetAbsorptionNILL() << G4endl;
      }

      G4cout << "mDataCurrent.trackID  = " << mDataCurrent.trackID;
      G4cout << " | .mOpBoundaryStatus  = " << std::setw(4) << mDataCurrent.GetOpBoundaryStatus() << G4endl;

      if (fVerbose > 2) {
        G4cout << "  | globalPosition    = (" << std::setw(9) << mDataCurrent.globalPosition.x();
        G4cout << "," << std::setw(9) << mDataCurrent.globalPosition.y();
        G4cout << "," << std::setw(9) << mDataCurrent.globalPosition.z() << ")" << G4endl;

        G4cout << "  | momentumDirection = (" << std::setw(9) << mDataCurrent.momentumDirection.x();
        G4cout << "," << std::setw(9) << mDataCurrent.momentumDirection.y();
        G4cout << "," << std::setw(9) << mDataCurrent.momentumDirection.z() << ")" << G4endl;

        G4cout << "  | polarization      = (" << std::setw(9) << mDataCurrent.polarization.x();
        G4cout << "," << std::setw(9) << mDataCurrent.polarization.y();
        G4cout << "," << std::setw(9) << mDataCurrent.polarization.z() << ")" << G4endl;

        G4cout << "  | globalTime        =  " << std::setw(9) << mDataCurrent.globalTime << G4endl;
        G4cout << "  | WLSNILL           =  " << std::setw(9) << mDataCurrent.GetWLSNILL() << G4endl;
        G4cout << "  | AbsorptionNILL    =  " << std::setw(9) << mDataCurrent.GetAbsorptionNILL() << G4endl;
      }

      G4cout << G4endl;
    }
  } // print
};  // end of the class

#endif

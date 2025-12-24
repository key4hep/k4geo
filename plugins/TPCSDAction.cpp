#include "DD4hep/Version.h"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "G4VProcess.hh"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

/// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
namespace sim {

  /**
   *  Geant4SensitiveAction<TPCSdData> sensitive detector for the special case of
   *  of a TPC, where every pad row is devided into two halfs in order to get
   *  the position from the crossing of the middle of the pad row from
   *  geant4 via volume boundary. Ported of Mokka/TPCSD04.cc
   *
   *  \author  F.Gaede ( ported from Mokka/TPCSD04.cc )
   *  \version 1.0
   *  \ingroup DD4HEP_SIMULATION
   */
  struct TPCSDData {

    // helper struct with Mokka controls ...
    struct {
      double TPCLowPtCut{};
      bool TPCLowPtStepLimit{};
      double TPCLowPtMaxHitSeparation{};
    } Control{};

    typedef Geant4HitCollection HitCollection;
    Geant4Sensitive* sensitive{};
    const BitFieldElement* layerField{};

    G4double fThresholdEnergyDeposit{};
    Geant4HitCollection* fHitCollection{};
    Geant4HitCollection* fSpaceHitCollection{};
    Geant4HitCollection* fLowPtHitCollection{};
    G4int fHCID{};
    G4int fSpaceHitCollectionID{};
    G4int fLowPtHitCollectionID{};

    G4ThreeVector CrossingOfPadRingCentre{};
    G4ThreeVector MomentumAtPadRingCentre{};
    G4double dEInPadRow{};
    G4double globalTimeAtPadRingCentre{};
    G4double pathLengthInPadRow{};

    G4double CumulativePathLength{};
    G4double CumulativeEnergyDeposit{};
    G4double CumulativeMeanTime{};
    G4ThreeVector CumulativeMeanPosition{};
    G4ThreeVector CumulativeMeanMomentum{};
    G4int CumulativeNumSteps{};

    std::map<int, G4double> padRowCentralRadii{};

    struct StepInfo {
      bool isValid{false};
      G4ThreeVector PreStepPointPosition{};
      G4ThreeVector PreStepPointMomentum{};
      G4ThreeVector PostStepPointPosition{};
      G4ThreeVector PostStepPointMomentum{};
      G4int TrackID{};
      G4int TrackTime{};
      G4int TrackPDGEncoding{};
      G4int cellID{};
      G4int innercopy{};
    };
    StepInfo previousStepInfo{};
    StepInfo stepAtEntranceToPadRingInfo{};

    TPCSDData() : fThresholdEnergyDeposit(0), fHCID(-1), fSpaceHitCollectionID(-1), fLowPtHitCollectionID(-1) {

      Control.TPCLowPtCut = CLHEP::MeV;
      Control.TPCLowPtStepLimit = false;
      Control.TPCLowPtMaxHitSeparation = 5. * CLHEP::mm;
    }

    /// Clear collected information and restart for new hit
    void clear() {
      // nothing to clear
    }

    /// return the layer number of the volume (either pre or post-position )
    int getCopyNumber(G4Step const* step, bool usePostPos) {

      int cellID = this->volID(step, usePostPos);

      return this->layerField->value(cellID);
    }

    /// Returns the volumeID of sensitive volume corresponding to the step (either pre or post-position )
    long long int volID(G4Step const* step, bool usePostPos = false) {

      Geant4StepHandler h(step);

      Geant4VolumeManager volMgr = Geant4Mapping::instance().volumeManager();

      VolumeID volID = (usePostPos ? volMgr.volumeID(h.postTouchable()) : volMgr.volumeID(h.preTouchable()));

      return volID;
    }

    void dumpStep(Geant4StepHandler& h, G4Step const* s) {
      sensitive->debug(
          "TPCSDAction considering step in detector %s ; "
          "prePos: (%e, %e, %e), "
          "postPos: (%e, %e, %e) ; "
          "preStatus : %s postStatus : %s ; "
          "preVolume : %s postVolume : %s ; "
          "CurrentCopyNumbers : %d, %d ; "
          "momentum : pre-step (%e, %e, %e), post-step (%e, %e, %e) ; "
          "PDG : %d ",
          h.sdName(s->GetPreStepPoint()).c_str(), h.prePos().X() / CLHEP::mm, h.prePos().Y() / CLHEP::mm,
          h.prePos().Z() / CLHEP::mm, h.postPos().X() / CLHEP::mm, h.postPos().Y() / CLHEP::mm,
          h.postPos().Z() / CLHEP::mm, h.preStepStatus(), h.postStepStatus(), h.volName(s->GetPreStepPoint()),
          h.volName(s->GetPostStepPoint()), getCopyNumber(s, false), getCopyNumber(s, true),
          s->GetPreStepPoint()->GetMomentum()[0] / CLHEP::GeV, s->GetPreStepPoint()->GetMomentum()[1] / CLHEP::GeV,
          s->GetPreStepPoint()->GetMomentum()[2] / CLHEP::GeV, s->GetPostStepPoint()->GetMomentum()[0] / CLHEP::GeV,
          s->GetPostStepPoint()->GetMomentum()[1] / CLHEP::GeV, s->GetPostStepPoint()->GetMomentum()[2] / CLHEP::GeV,
          s->GetTrack()->GetDefinition()->GetPDGEncoding());
    }

    /// Method for generating hit(s) using the information of G4Step object.
    G4bool process(G4Step const* step, G4TouchableHistory*) {

      fHitCollection = sensitive->collection(0);
      fSpaceHitCollection = sensitive->collection(1);
      fLowPtHitCollection = sensitive->collection(2);

      Geant4StepHandler h(step);
      dumpStep(h, step);

      // FIXME:
      // a particle that crosses the boundry between two pad-ring halves will have the hit
      // placed on this surface at the last crossing point, and will be assinged the total energy
      // deposited in the whole pad-ring. This is a possible source of bias for the hit

      if (fabs(step->GetTrack()->GetDefinition()->GetPDGCharge()) < 0.01)
        return true;

      // this is the info we may need to make a hit
      StepInfo thisStepInfo;
      thisStepInfo.isValid = true;
      thisStepInfo.PreStepPointPosition = step->GetPreStepPoint()->GetPosition();
      thisStepInfo.PreStepPointMomentum = step->GetPreStepPoint()->GetMomentum();
      thisStepInfo.PostStepPointPosition = step->GetPostStepPoint()->GetPosition();
      thisStepInfo.PostStepPointMomentum = step->GetPostStepPoint()->GetMomentum();
      thisStepInfo.TrackID = step->GetTrack()->GetTrackID();
      thisStepInfo.TrackTime = step->GetTrack()->GetGlobalTime();
      thisStepInfo.TrackPDGEncoding = step->GetTrack()->GetDefinition()->GetPDGEncoding();
      thisStepInfo.cellID = sensitive->cellID(step);
      thisStepInfo.innercopy = getCopyNumber(step, false);

      //=========================================================================================================

      if (thisStepInfo.PostStepPointMomentum.perp2() >= pow(Control.TPCLowPtCut, 2)) { // high pT class

        // ===================== first check we have no left-over low-pt stuff
        // This step does not continue the previous path. Deposit the energy up to the previous step
        if (CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
          DepositLowPtHit(previousStepInfo);
        } else {
          ResetCumulativeVariables(); // set low pt cumulation to zero
        }
        //=== end of low pt cleaning-up

        //=========================================================================================================

        // keep info about the step which enters the pad ring
        if (!stepAtEntranceToPadRingInfo.isValid) { // first step in this padrow
          stepAtEntranceToPadRingInfo = thisStepInfo;
        }

        // Step finishes at a geometric boundry
        if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {

          // accumulate step in current pad row
          dEInPadRow += step->GetTotalEnergyDeposit();
          pathLengthInPadRow += step->GetStepLength();
          if (thisStepInfo.innercopy == getCopyNumber(step, true)) {
            // pre == post; step within the same pair of upper and lower pad ring halves
            // this step must have ended on the boundry between these two pad ring halfves
            // record the tracks coordinates at this position
            CrossingOfPadRingCentre = thisStepInfo.PostStepPointPosition;
            MomentumAtPadRingCentre = thisStepInfo.PostStepPointMomentum;
            globalTimeAtPadRingCentre = thisStepInfo.TrackTime;

            // here we memorise the padrow positions, we need them in some special cases
            if (padRowCentralRadii.find(thisStepInfo.innercopy) == padRowCentralRadii.end()) {
              padRowCentralRadii[thisStepInfo.innercopy] = CrossingOfPadRingCentre.perp();
            }

          } else { // has crossed into new padrow, consider making hit
            DepositHiPtHit(thisStepInfo);
          }

        } else { // case for which the step remains within geometric volume

          // accumulate
          dEInPadRow += step->GetTotalEnergyDeposit();
          pathLengthInPadRow += step->GetStepLength();

          if (step->GetPostStepPoint()->GetKineticEnergy() == 0) { // particle stopped in padring
            DepositHiPtHit(thisStepInfo);
          } else if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() ==
                     "StepLimiter") { // step limited by distance
            // write out a zero energy hit in the spacehitcollection
            Geant4Tracker::Hit* hit =
                new Geant4Tracker::Hit(thisStepInfo.TrackID, thisStepInfo.TrackPDGEncoding,
                                       0.0,                        // dE set to ZERO
                                       globalTimeAtPadRingCentre); // this time may or may not be defined already
            hit->position = 0.5 * (thisStepInfo.PreStepPointPosition + thisStepInfo.PostStepPointPosition);
            hit->momentum = thisStepInfo.PostStepPointMomentum;
            hit->length = step->GetStepLength();
            hit->cellID = thisStepInfo.cellID;
            fSpaceHitCollection->add(hit);
          }
        }

      }
      //=========================================================================================================
      //   ptSQRD < (Control.TPCLowPtCut*Control.TPCLowPtCut)

      else if (Control.TPCLowPtStepLimit) { // low pt tracks will be treated differently as their step length is limited
                                            // by the special low pt steplimiter

        //--------------------------------
        // first check if there is padrow stuff left over to deposit
        if (dEInPadRow > 0) {
          sensitive->debug("left-over padrow energy from high pt tracks: deposit HiPtHit");
          DepositHiPtHit(previousStepInfo); // use stored previous step
        }
        //--------------------------------

        if ((previousStepInfo.PostStepPointPosition - thisStepInfo.PreStepPointPosition).mag() > 1.0e-6 * CLHEP::mm) {
          // This step does not continue the previous path. Deposit the energy and begin a new Pt hit.
          if (CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
            // dumpStep( h , step ) ;
            DepositLowPtHit(previousStepInfo);
          }

          else {
            // reset the cumulative variables if the hit has not been deposited.
            // The previous track has ended and the cumulated energy left at the end
            // was not enough to ionize
            ResetCumulativeVariables();
          }
        }

        CumulateLowPtStep(step);

        // check whether to deposit the hit
        if (step->GetPostStepPoint()->GetKineticEnergy() == 0) {   // particle stopped
          if (CumulativeEnergyDeposit > fThresholdEnergyDeposit) { // enough energy
            DepositLowPtHit(thisStepInfo);                         // make hit ending with this step
          } else {
            ResetCumulativeVariables(); // not enough energy: reset/ignore it
          }
        }

        if ((CumulativePathLength > Control.TPCLowPtMaxHitSeparation)) {

          // hit is deposited because the step limit is reached and there is enough energy
          // to ionize

          if (CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
            // dumpStep( h , step ) ;
            DepositLowPtHit(thisStepInfo);
          }
        }
      }

      // keep track of previous step properties
      previousStepInfo = thisStepInfo;

      return true;
    }

    /// Post-event action callback
    void endEvent(const G4Event* /* event */) {

      ResetPadrowVariables();

      // ===================== finally check we have no left-over low-pt stuff
      if (CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
        DepositLowPtHit(previousStepInfo); // Deposit the energy
      }
      ResetCumulativeVariables();
    }

    void ResetPadrowVariables() // DJ added
    {
      dEInPadRow = 0.0;
      globalTimeAtPadRingCentre = 0.0;
      pathLengthInPadRow = 0.0;
      for (int i = 0; i < 3; i++) {
        CrossingOfPadRingCentre[i] = 0.0;
        MomentumAtPadRingCentre[i] = 0.0;
      }
      stepAtEntranceToPadRingInfo.isValid = false;
    }

    void DepositHiPtHit(StepInfo mainStepInfo) {

      if (!mainStepInfo.isValid || !stepAtEntranceToPadRingInfo.isValid) {
        sensitive->error(
            "ERROR from TPCSDAction DepositHiPtHit, invalid StepInfo (main step %d, step at padrow entrance %d",
            mainStepInfo.isValid, stepAtEntranceToPadRingInfo.isValid);
      }

      if (dEInPadRow > fThresholdEnergyDeposit) {

        bool rareError = false;

        if (CrossingOfPadRingCentre[0] < 0.1 && CrossingOfPadRingCentre[1] < 0.1 && CrossingOfPadRingCentre[2] < 0.1) {
          // series of steps did not cross the centre of pad row; make reasonable estimate
          int innercopy = mainStepInfo.innercopy;
          if (padRowCentralRadii.find(innercopy) != padRowCentralRadii.end()) { // we know radius of this pad row
            // average of first and last points of this series of steps
            const G4ThreeVector AvePos =
                0.5 * (stepAtEntranceToPadRingInfo.PreStepPointPosition + mainStepInfo.PostStepPointPosition);
            G4double radius = AvePos.perp();
            CrossingOfPadRingCentre =
                AvePos * (padRowCentralRadii[innercopy] / radius); // move radially to centre of pad row
            // time and momentum: average of intial and final
            globalTimeAtPadRingCentre = 0.5 * (stepAtEntranceToPadRingInfo.TrackTime + mainStepInfo.TrackTime);
            MomentumAtPadRingCentre =
                0.5 * (stepAtEntranceToPadRingInfo.PreStepPointMomentum + mainStepInfo.PostStepPointMomentum);
          } else { // have not yet seen this pad row in this job
            // in principle we can get this from geometry information, but don't want to add extra dependencies...
            sensitive->warning("WARNING: dont yet know radius of this pad row."
                               "ignoring this energy deposit (should be very rare,"
                               "eg possibly in first track(s) of first event)");
            rareError = true;
          }
        }

        if (!rareError) {
          Geant4Tracker::Hit* hit = new Geant4Tracker::Hit(mainStepInfo.TrackID, mainStepInfo.TrackPDGEncoding,
                                                           dEInPadRow, globalTimeAtPadRingCentre);
          hit->position = CrossingOfPadRingCentre;
          hit->momentum = MomentumAtPadRingCentre;
          hit->length = pathLengthInPadRow;
          hit->cellID = mainStepInfo.cellID;
          fHitCollection->add(hit);
          sensitive->printM2("+++ TrackID:%6d [%s] CREATE TPC hit at pad row crossing :"
                             " %e MeV  Pos:%8.2f %8.2f %8.2f",
                             mainStepInfo.TrackID, sensitive->c_name(), dEInPadRow, hit->position.X() / CLHEP::mm,
                             hit->position.Y() / CLHEP::mm, hit->position.Z() / CLHEP::mm);
        }
      } // en threshold
      ResetPadrowVariables();
    }

    void ResetCumulativeVariables() {
      CumulativeMeanPosition.set(0.0, 0.0, 0.0);
      CumulativeMeanMomentum.set(0.0, 0.0, 0.0);
      CumulativeMeanTime = 0;
      CumulativeNumSteps = 0;
      CumulativeEnergyDeposit = 0;
      CumulativePathLength = 0;
    }

    void DepositLowPtHit(StepInfo stInf) {
      if (!stInf.isValid) {
        sensitive->error("ERROR from TPCSDAction DepositLowPtHit: invalid step info, ignoring");
        return;
      }
      Geant4Tracker::Hit* hit = new Geant4Tracker::Hit(stInf.TrackID, stInf.TrackPDGEncoding, CumulativeEnergyDeposit,
                                                       CumulativeMeanTime / CumulativeNumSteps);

      hit->position = CumulativeMeanPosition / CumulativeNumSteps;
      hit->momentum = CumulativeMeanMomentum / CumulativeNumSteps;
      hit->length = CumulativePathLength;
      hit->cellID = stInf.cellID;
      fLowPtHitCollection->add(hit);
      // reset the cumulative variables after positioning the hit
      ResetCumulativeVariables();
    }

    void CumulateLowPtStep(G4Step const* step) {

      const G4ThreeVector meanPosition =
          (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition()) / 2;
      const G4ThreeVector meanMomentum =
          (step->GetPreStepPoint()->GetMomentum() + step->GetPostStepPoint()->GetMomentum()) / 2;

      ++CumulativeNumSteps;
      CumulativeMeanPosition += (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition()) / 2;
      CumulativeMeanMomentum += (step->GetPreStepPoint()->GetMomentum() + step->GetPostStepPoint()->GetMomentum()) / 2;
      CumulativeEnergyDeposit += step->GetTotalEnergyDeposit();
      CumulativePathLength += step->GetStepLength();
      CumulativeMeanTime += step->GetTrack()->GetGlobalTime();
    }
  };

  /// Initialization overload for specialization
  template <>
  void Geant4SensitiveAction<TPCSDData>::initialize() {
    eventAction().callAtEnd(&m_userData, &TPCSDData::endEvent);

    declareProperty("TPCLowPtCut", m_userData.Control.TPCLowPtCut);
    declareProperty("TPCLowPtStepLimit", m_userData.Control.TPCLowPtStepLimit);
    declareProperty("TPCLowPtMaxHitSeparation", m_userData.Control.TPCLowPtMaxHitSeparation);

    m_userData.fThresholdEnergyDeposit = m_sensitive.energyCutoff();
    m_userData.sensitive = this;

    IDDescriptor dsc = m_sensitive.idSpec();
    m_userData.layerField = dsc.field("layer");
  }

  /// Define collections created by this sensitivie action object
  template <>
  void Geant4SensitiveAction<TPCSDData>::defineCollections() {
    m_collectionID = defineCollection<Geant4Tracker::Hit>(m_sensitive.readout().name());
    m_userData.fSpaceHitCollectionID = defineCollection<Geant4Tracker::Hit>("TPCSpacePointCollection");
    m_userData.fLowPtHitCollectionID = defineCollection<Geant4Tracker::Hit>("TPCLowPtCollection");
  }

  /// Method for generating hit(s) using the information of G4Step object.
  template <>
  void Geant4SensitiveAction<TPCSDData>::clear(G4HCofThisEvent*) {
    m_userData.clear();
  }

  /// Method for generating hit(s) using the information of G4Step object.
  template <>
  G4bool Geant4SensitiveAction<TPCSDData>::process(G4Step const* step, G4TouchableHistory* history) {
    return m_userData.process(step, history);
  }

  typedef Geant4SensitiveAction<TPCSDData> TPCSDAction;

} // namespace sim
} // namespace dd4hep

#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(TPCSDAction)

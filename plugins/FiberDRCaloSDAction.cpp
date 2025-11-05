// DD4hep Framework include files
#include "DD4hep/Segmentations.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4Random.h"
#include "DDG4/Geant4SensDetAction.inl"

#include "G4OpticalPhoton.hh"
// k4geo Framework include files

#include "DRCaloFastSimModel.h"
#include "FiberDRCaloSDAction.h"

// Geant4 include files
#include "G4HCofThisEvent.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpProcessSubType.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHistory.hh"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

/// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
namespace sim {

  /*
   *  Geant4SensitiveAction<DRcaloSDAction> sensitive detector for the Dual-readout calorimeter
   *
   *  \author  Sungwon Kim, Sanghyun Ko
   *  \version 1.0
   *  \ingroup DD4HEP_SIMULATION
   */

  struct DRCData {
  public:
    bool skipScint = true;
    DRCFiberModel fastfiber;

    G4int fWavBin;
    G4int fTimeBin;
    G4float fWavlenStart;
    G4float fWavlenEnd;
    G4float fTimeStart;
    G4float fTimeEnd;
    G4float fWavlenStep;
    G4float fTimeStep;

  public:
    G4double wavToE(G4double wav) const { return CLHEP::h_Planck * CLHEP::c_light / wav; }

    int findWavBin(G4double en) const {
      int i = 0;
      for (; i < fWavBin + 1; i++) {
        if (en < wavToE((fWavlenStart - static_cast<float>(i) * fWavlenStep) * CLHEP::nm))
          break;
      }

      return fWavBin + 1 - i;
    }

    int findTimeBin(G4double stepTime) const {
      int i = 0;
      for (; i < fTimeBin + 1; i++) {
        if (stepTime < ((fTimeStart + static_cast<float>(i) * fTimeStep) * CLHEP::ns))
          break;
      }

      return i;
    }

    // default constructor
    DRCData() : fWavBin(120), fTimeBin(2000), fWavlenStart(900.), fWavlenEnd(300.), fTimeStart(-100.), fTimeEnd(100.) {
      fWavlenStep = (fWavlenStart - fWavlenEnd) / (float)fWavBin;
      fTimeStep = (fTimeEnd - fTimeStart) / (float)fTimeBin;
    }
  }; // struct DRCData

  template <>
  Geant4SensitiveAction<DRCData>::Geant4SensitiveAction(Geant4Context* ctxt, const std::string& nam, DetElement det,
                                                        Detector& desc)
      : Geant4Sensitive(ctxt, nam, det, desc), m_collectionName(), m_collectionID(0) {
    declareProperty("skipScint", m_userData.skipScint = true);
    declareProperty("ReadoutName", m_readoutName);
    declareProperty("CollectionName", m_collectionName);
    declareProperty("timeStart", m_userData.fTimeStart); // in ns
    declareProperty("timeEnd", m_userData.fTimeEnd);     // in ns
    // delegate to user's responsibility to ensure timeStep divides (timeEnd - timeStart)
    declareProperty("timeStep", m_userData.fTimeStep); // in ns
    declareProperty("timeBin", m_userData.fTimeBin);   // in ns
    declareProperty("verbose", m_userData.fastfiber.fVerbose);
    declareProperty("safety", m_userData.fastfiber.fSafety);
    InstanceCount::increment(this);

    m_hitCreationMode = HitCreationFlags::DETAILED_MODE; // always store step position
  }

  /// Define collections created by this sensitive action object
  template <>
  void Geant4SensitiveAction<DRCData>::defineCollections() {
    std::string readout_name = m_sensitive.readout().name();
    m_collectionID = defineCollection<Geant4DRCalorimeter::Hit>(m_sensitive.readout().name());
    std::cout << "defineCollection Geant4DRCalorimeter readout_name   : " << readout_name << std::endl;
    std::cout << "defineCollection Geant4DRCalorimeter m_collectionID : " << m_collectionID << std::endl;

    if (m_userData.skipScint) {
      defineCollection<Geant4Calorimeter::Hit>(std::string(m_sensitive.readout().name()) + "_scint");
      std::cout << "defineCollection Geant4Calorimeter readout_name   : " << readout_name + "_scint" << std::endl;
      std::cout << "defineCollection Geant4Calorimeter m_collectionID : " << m_collectionID + 1 << std::endl;
    }
  }

  /// Method for generating hit(s) using the information of G4Step object.
  template <>
  G4bool Geant4SensitiveAction<DRCData>::process(G4Step const* step, G4TouchableHistory*) {
    // optical photons traveling through the cladding is not interesting
    // but consume lots of computation power
    // let's kill optical photons inside the cladding whose status is not StepTooSmall
    if (step->GetPreStepPoint() &&
        step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetNoDaughters() == 1 &&
        step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
      // 1e-9 is the default tolerance
      // for the warnings, see
      // https://geant4-forum.web.cern.ch/t/error-occurs-when-an-optical-photon-hits-the-edge-of-a-cubic-scintillator/8748
      // we're not interested in the optical photon position
      // but only the number of photons (and timing)
      // so we assume it's fine to ignore the warnings
      if (step->GetStepLength() > 1e-8 * CLHEP::mm)
        step->GetTrack()->SetTrackStatus(fStopAndKill);
    }
    // the fiber itself is the SD
    // we skip the scintillation process and only account for the Cherenkov process
    // remember to turn off the scintillation process!
    if (m_userData.skipScint) {
      // we need to move the touchable to the tower to retrieve volID
      auto* touchable = const_cast<G4VTouchable*>(step->GetPostStepPoint()->GetTouchable());
      const auto* logicalVol = touchable->GetVolume()->GetLogicalVolume();

      // we're not interested in the world or assembly volume
      if (touchable->GetHistoryDepth() < 2)
        return false;

      // we're only interested in the fiber core for optical photons
      if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() &&
          logicalVol->GetNoDaughters() != 0)
        return false;

      // now let's make the touchable points to the tower
      // world -> assembly -> tower
      touchable->MoveUpHistory(touchable->GetHistoryDepth() - 2);

      // now the touchable is the tower
      // get local position and volumeID
      // dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
      // dd4hep::VolumeID volID = volMgr.volumeID(touchable);
      // note the above method started to throw warnings since the PR DD4hep#1390
      // so we switch to the copy number

      // trick: we retrieve the tower's volID by using the fact that
      // in the DRconstructor.cpp, the first 32bits of the volID
      // contain the information relevant to the tower
      // and this 32bit integer becomes the tower's copy number
      auto copyNo = touchable->GetCopyNumber();

      // ideally the below operation should be GridDRcalo_k4geo::convertFirst32to64
      // but this is trivial, so let's avoid doing dynamic cast
      dd4hep::VolumeID volID = static_cast<dd4hep::VolumeID>(copyNo);
      G4ThreeVector global = step->GetPostStepPoint()->GetPosition();
      G4ThreeVector local = touchable->GetHistory()->GetTopTransform().TransformPoint(global);

      // convert G4 position to dd4hep position
      dd4hep::Position loc(local.x() * dd4hep::millimeter / CLHEP::millimeter,
                           local.y() * dd4hep::millimeter / CLHEP::millimeter,
                           local.z() * dd4hep::millimeter / CLHEP::millimeter);
      dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter,
                            global.y() * dd4hep::millimeter / CLHEP::millimeter,
                            global.z() * dd4hep::millimeter / CLHEP::millimeter);

      // retrieve cellID and distinguish Cherenkov & scintillation fibers
      auto cID = m_segmentation->cellID(loc, glob, volID);
      auto ceren = static_cast<dd4hep::DDSegmentation::VolumeID>(m_segmentation->decoder()->get(cID, "c"));
      bool IsCeren = static_cast<bool>(ceren);

      Geant4HitCollection* coll = collection(m_collectionID);

      if (IsCeren) {
        // Cherenkov fiber
        // skip anything else than optical photons
        if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
          return false;

        const auto* track = step->GetTrack();

        // reset when moving to the next track
        if (m_userData.fastfiber.mDataCurrent.trackID != track->GetTrackID())
          m_userData.fastfiber.reset();

        // need repetitive total internal reflection
        if (!m_userData.fastfiber.check_trigger(track))
          return false;

        // absorption
        if (m_userData.fastfiber.fKill) {
          step->GetTrack()->SetTrackStatus(fStopAndKill);

          return false;
        }

        // backward transportation
        if (m_userData.fastfiber.mTransportUnit < 0.) {
          step->GetTrack()->SetTrackStatus(fStopAndKill);

          return false;
        }

        // for timing measurement
        double timeUnit = m_userData.fastfiber.mDataCurrent.globalTime - m_userData.fastfiber.mDataPrevious.globalTime;
        double timeShift = timeUnit * m_userData.fastfiber.mNtransport;
        G4double energy = step->GetTrack()->GetTotalEnergy();

        // default hit (optical photon count)
        Geant4DRCalorimeter::Hit* drHit =
            coll->find<Geant4DRCalorimeter::Hit>(CellIDCompare<Geant4DRCalorimeter::Hit>(cID));

        if (!drHit) {
          drHit = new Geant4DRCalorimeter::Hit(m_userData.fWavlenStep, m_userData.fTimeStep);
          drHit->cellID = cID;
          drHit->position = m_segmentation->position(cID) * CLHEP::mm / dd4hep::mm; // segmentation gives dd4hep unit
          drHit->SetSiPMnum(cID);
          drHit->SetTimeStart(m_userData.fTimeStart);
          drHit->SetTimeEnd(m_userData.fTimeEnd);
          drHit->SetWavlenMax(m_userData.fWavlenStart);
          drHit->SetWavlenMin(m_userData.fWavlenEnd);
          coll->add(cID, drHit);
        }

        // everything should be in the G4 unit
        // (approximate) timing at the end of the fiber
        G4double hitTime = step->GetPostStepPoint()->GetGlobalTime() + timeShift;

        drHit->photonCount();
        int wavBin = m_userData.findWavBin(energy);
        drHit->CountWavlenSpectrum(wavBin);
        int timeBin = m_userData.findTimeBin(hitTime);
        drHit->CountTimeStruct(timeBin);

        // finally kill optical photon
        step->GetTrack()->SetTrackStatus(fStopAndKill);

        return true;
      } else {
        // scintillation calo hit

        // assume nPhoton_scint >> nPhoton_cherenkov
        // kill optical photons (from Cherenkov process)
        if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
          step->GetTrack()->SetTrackStatus(fStopAndKill);

          return false;
        }

        // copy-paste of the dd4hep scintillation calorimeter SD
        Geant4HitCollection* coll_scint = collection(m_collectionID + 1);
        Geant4Calorimeter::Hit* caloHit =
            coll_scint->find<Geant4Calorimeter::Hit>(CellIDCompare<Geant4Calorimeter::Hit>(cID));
        HitContribution contrib = Geant4Calorimeter::Hit::extractContribution(step, true);

        if (!caloHit) {
          caloHit = new Geant4Calorimeter::Hit(glob);
          caloHit->cellID = cID;
          caloHit->position = m_segmentation->position(cID) * CLHEP::mm / dd4hep::mm; // segmentation gives dd4hep unit
          coll_scint->add(cID, caloHit);
        }

        caloHit->truth.emplace_back(contrib);
        caloHit->energyDeposit += contrib.deposit;

        Geant4StepHandler h(step);
        mark(h.track);

        return true;
      } // !IsCeren
    } else {
      // no skipping optical photon propagation
      // SiPM wafers are SD in this case
      if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;

      typedef Geant4DRCalorimeter::Hit Hit;

      Geant4HitCollection* coll = collection(m_collectionID);

      auto* touchable = step->GetPostStepPoint()->GetTouchable();
      dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
      dd4hep::VolumeID volID = volMgr.volumeID(touchable);
      G4ThreeVector global = step->GetPostStepPoint()->GetPosition();
      G4ThreeVector local = touchable->GetHistory()->GetTopTransform().TransformPoint(global);

      dd4hep::Position loc(local.x() * dd4hep::millimeter / CLHEP::millimeter,
                           local.y() * dd4hep::millimeter / CLHEP::millimeter,
                           local.z() * dd4hep::millimeter / CLHEP::millimeter);
      dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter,
                            global.y() * dd4hep::millimeter / CLHEP::millimeter,
                            global.z() * dd4hep::millimeter / CLHEP::millimeter);

      auto cID = m_segmentation->cellID(loc, glob, volID); // This returns cID corresponding to SiPM Wafer
      Hit* hit = coll->find<Hit>(CellIDCompare<Hit>(cID));

      G4double hitTime = step->GetPostStepPoint()->GetGlobalTime();
      G4double energy = step->GetTrack()->GetTotalEnergy();

      if (!hit) {
        hit = new Geant4DRCalorimeter::Hit(m_userData.fWavlenStep, m_userData.fTimeStep);
        hit->cellID = cID;
        hit->position = m_segmentation->position(cID) * CLHEP::mm / dd4hep::mm; // segmentation gives dd4hep unit
        hit->SetSiPMnum(cID);
        hit->SetTimeStart(m_userData.fTimeStart);
        hit->SetTimeEnd(m_userData.fTimeEnd);
        hit->SetWavlenMax(m_userData.fWavlenStart);
        hit->SetWavlenMin(m_userData.fWavlenEnd);
        coll->add(cID, hit);
      }

      hit->photonCount();
      int wavBin = m_userData.findWavBin(energy);
      hit->CountWavlenSpectrum(wavBin);
      int timeBin = m_userData.findTimeBin(hitTime);
      hit->CountTimeStruct(timeBin);

      return true;
    } // !skipScint
  } // Geant4SensitiveAction::process

  typedef Geant4SensitiveAction<DRCData> DRCaloSDAction;
} // namespace sim
} // namespace dd4hep

#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE(DRCaloSDAction)

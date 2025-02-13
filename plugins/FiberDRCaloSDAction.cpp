// DD4hep Framework include files
#include "DD4hep/Segmentations.h"
#include "DDG4/Geant4Random.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Geant4Mapping.h"

#include "G4OpticalPhoton.hh"
// k4geo Framework include files

#include "FiberDRCaloSDAction.h"
#include "DRCaloFastSimModel.h"

// Geant4 include files
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4ThreeVector.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpProcessSubType.hh"

#if DD4HEP_VERSION_GE(1, 21)
#define GEANT4_CONST_STEP const
#else
#define GEANT4_CONST_STEP
#endif

/// Namespace for the AIDA detector description toolkit
namespace dd4hep
{

  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim
  {

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

    private:
      // from 900 nm to 460 nm
      const std::vector<double> fGraph_X = {
        1.37760 * CLHEP::eV,
        1.45864 * CLHEP::eV,
        1.54980 * CLHEP::eV,
        1.65312 * CLHEP::eV,
        1.71013 * CLHEP::eV,
        1.77120 * CLHEP::eV,
        1.83680 * CLHEP::eV,
        1.90745 * CLHEP::eV,
        1.98375 * CLHEP::eV,
        2.06640 * CLHEP::eV,

        2.10143 * CLHEP::eV,
        2.13766 * CLHEP::eV,
        2.17516 * CLHEP::eV,
        2.21400 * CLHEP::eV,
        2.25426 * CLHEP::eV,
        2.29600 * CLHEP::eV,
        2.33932 * CLHEP::eV,
        2.38431 * CLHEP::eV,
        2.43106 * CLHEP::eV,
        2.47968 * CLHEP::eV,

        2.53029 * CLHEP::eV,
        2.58300 * CLHEP::eV,
        2.63796 * CLHEP::eV,
        2.69531 * CLHEP::eV,
        2.75520 * CLHEP::eV,
        2.81782 * CLHEP::eV,
        2.88335 * CLHEP::eV,
        2.95200 * CLHEP::eV,
        3.09960 * CLHEP::eV,
        3.54241 * CLHEP::eV,

        4.13281 * CLHEP::eV
      };
      // filter efficiency of the Kodak Wratten No.9 filter
      const std::vector<double> fKodakEff = {
        0.903,
        0.903,
        0.903,
        0.903,
        0.903,
        0.903,
        0.902,
        0.901,
        0.898,
        0.895,

        0.893,
        0.891,
        0.888,
        0.883,
        0.87,
        0.838,
        0.76,
        0.62,
        0.488,
        0.345,

        0.207,
        0.083,
        0.018,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,

        0.
      };

      // SiPM efficiency Hamamatsu S14160-1310PS
      // TODO migrate this part to the digitization step!
      // Note: Ideally, this should be part of the digitization step.
      // (not the simulation step)
      // But, to do this, we need to store the distribution
      // of the optical photon wavelength.
      // While we can develop another feature to enable this,
      // let's emulate the SiPM efficiency in the simulation step for now.
      // We just need a working code and without this,
      // the number of Cherenkov photon will be order of magnitude higher
      const std::vector<double> fSipmEff = {
        0.02,
        0.025,
        0.045,
        0.06,
        0.0675,
        0.075,
        0.0925,
        0.11,
        0.125,
        0.14,

        0.146,
        0.152,
        0.158,
        0.164,
        0.17,
        0.173,
        0.176,
        0.178,
        0.179,
        0.18,

        0.181,
        0.182,
        0.183,
        0.184,
        0.18,
        0.173,
        0.166,
        0.158,
        0.15,
        0.12,

        0.05
      };

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

      // Linear interpolation function for calculating the efficiency of yellow filter
      // used in rejectedByYellowFilter
      double getFilterEff(const std::vector<double>& yarray, const G4double G4energy) const {
        // If the photon energy <= 1.37760 eV, than return maximum filter efficiency
        if (G4energy <= fGraph_X.at(0))
          return yarray.at(0);

        for(unsigned idx = 1; idx < yarray.size(); ++idx) {
          if (G4energy <= fGraph_X.at(idx)) {
            double x1 = fGraph_X.at(idx-1);
            double x2 = fGraph_X.at(idx);
            double y1 = yarray.at(idx-1);
            double y2 = yarray.at(idx);

            // return linear interpolated filter efficiency
            return (y1 + ((y2 - y1) / (x2 - x1))*(G4energy - x1));
          }
        }

        return 0.;
      }

      // If true, then the photon is rejected by yellow filter
      bool rejectedByYellowFilter(G4double G4energy, double rndVal) const {
        const double FilterEff = getFilterEff(fKodakEff,G4energy); // Get efficiency of filter using photon's energy

        // filter efficiency == probability of photon accepted by filter
        // == Probability of random value (from uniform distribution with range 0 ~ 1) smaller than filter efficiency
        // So if the rndVal is larger than the FilterEff, than the photon is rejected
        return (rndVal > FilterEff);
      }

      // check sipm efficiency
      // TODO migrate this to the digitization step
      bool rejectedBySiPM(double G4energy, double rndVal) const {
        return (rndVal > getFilterEff(fSipmEff,G4energy));
      }

      DRCData()
      : fWavBin(120), fTimeBin(650), fWavlenStart(900.),
        fWavlenEnd(300.), fTimeStart(5.), fTimeEnd(70.)
      {
        fWavlenStep = (fWavlenStart - fWavlenEnd) / (float)fWavBin;
        fTimeStep = (fTimeEnd - fTimeStart) / (float)fTimeBin;
      }
    }; // struct DRCData

    template <>
    Geant4SensitiveAction<DRCData>::Geant4SensitiveAction(Geant4Context* ctxt,
                                                          const std::string& nam,
                                                          DetElement det,
                                                          Detector& desc)
    : Geant4Sensitive(ctxt,nam,det,desc), m_collectionName(), m_collectionID(0)
    {
      declareProperty("skipScint", m_userData.skipScint = true);
      declareProperty("ReadoutName", m_readoutName);
      declareProperty("CollectionName", m_collectionName);
      initialize();
      InstanceCount::increment(this);

      m_userData.fastfiber.fSafety = 1;
      m_userData.fastfiber.fVerbose = 0;
    }

    /// Define collections created by this sensitive action object
    template <>
    void Geant4SensitiveAction<DRCData>::defineCollections() {
      std::string readout_name = m_sensitive.readout().name();
      m_collectionID = defineCollection<Geant4DRCalorimeter::Hit>(m_sensitive.readout().name());
      std::cout << "defineCollection Geant4DRCalorimeter readout_name   : " << readout_name << std::endl;
      std::cout << "defineCollection Geant4DRCalorimeter m_collectionID : " << m_collectionID << std::endl;

      if (m_userData.skipScint) {
        defineCollection<Geant4Calorimeter::Hit>(std::string(m_sensitive.readout().name())+"_scint");
        std::cout << "defineCollection Geant4Calorimeter readout_name   : " << readout_name + "_scint" << std::endl;
        std::cout << "defineCollection Geant4Calorimeter m_collectionID : " << m_collectionID + 1 << std::endl;
      }
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <>
    G4bool Geant4SensitiveAction<DRCData>::process(G4Step GEANT4_CONST_STEP *step, G4TouchableHistory *) {
      // optical photons traveling through the cladding is not interesting
      // but consume lots of computation power
      // let's kill optical photons inside the cladding whose status is not StepTooSmall
      if ( step->GetPreStepPoint() &&
           step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetNoDaughters()==1 &&
           step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
        // 1e-9 is the default tolerance
        // for the warnings, see https://geant4-forum.web.cern.ch/t/error-occurs-when-an-optical-photon-hits-the-edge-of-a-cubic-scintillator/8748
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

        // we're only interested in the fiber core
        if ( logicalVol->GetNoDaughters()!=0 )
          return false;

        // now let's make the touchable points to the tower
        // world -> assembly -> tower
        touchable->MoveUpHistory( touchable->GetHistoryDepth()-2 );

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
        dd4hep::Position loc(local.x() * dd4hep::millimeter / CLHEP::millimeter, local.y() * dd4hep::millimeter / CLHEP::millimeter, local.z() * dd4hep::millimeter / CLHEP::millimeter);
        dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter, global.y() * dd4hep::millimeter / CLHEP::millimeter, global.z() * dd4hep::millimeter / CLHEP::millimeter);

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

          // check wheter the optical photon will be rejected by the SiPM
          // TODO migrate this to the digitization step
          Geant4Event& evt = context()->event();
          Geant4Random& rnd = evt.random(); // Get random generator
          double rndVal = rnd.uniform(0, 1); // Get random number from uniform distribution [0, 1]
          G4double energy = step->GetTrack()->GetTotalEnergy();

          if (m_userData.rejectedBySiPM(energy, rndVal)) {
            step->GetTrack()->SetTrackStatus(fStopAndKill);

            return false;
          }

          // default hit (optical photon count)
          Geant4DRCalorimeter::Hit* drHit = coll->find<Geant4DRCalorimeter::Hit>(CellIDCompare<Geant4DRCalorimeter::Hit>(cID));

          if (!drHit) {
            drHit = new Geant4DRCalorimeter::Hit(m_userData.fWavlenStep, m_userData.fTimeStep);
            drHit->cellID = cID;
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

          drHit->position = glob;

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
          Geant4HitCollection* coll_scint = collection(m_collectionID+1);
          Geant4Calorimeter::Hit* caloHit = coll_scint->find<Geant4Calorimeter::Hit>(CellIDCompare<Geant4Calorimeter::Hit>(cID));
          HitContribution contrib = Geant4Calorimeter::Hit::extractContribution(step,true);

          if ( !caloHit ) {
            caloHit = new Geant4Calorimeter::Hit(glob);
            caloHit->cellID = cID;
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

        Geant4HitCollection *coll = collection(m_collectionID);

        auto* touchable = step->GetPostStepPoint()->GetTouchable();
        dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
        dd4hep::VolumeID volID = volMgr.volumeID(touchable);
        G4ThreeVector global = step->GetPostStepPoint()->GetPosition();
        G4ThreeVector local = touchable->GetHistory()->GetTopTransform().TransformPoint(global);

        dd4hep::Position loc(local.x() * dd4hep::millimeter / CLHEP::millimeter, local.y() * dd4hep::millimeter / CLHEP::millimeter, local.z() * dd4hep::millimeter / CLHEP::millimeter);
        dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter, global.y() * dd4hep::millimeter / CLHEP::millimeter, global.z() * dd4hep::millimeter / CLHEP::millimeter);

        auto cID = m_segmentation->cellID(loc, glob, volID); // This returns cID corresponding to SiPM Wafer
        Hit *hit = coll->find<Hit>(CellIDCompare<Hit>(cID));

        G4double hitTime = step->GetPostStepPoint()->GetGlobalTime();
        G4double energy = step->GetTrack()->GetTotalEnergy();

        // Apply yellow filter here
        // Get random number from (0~1) uniform distribution
        // If the photon is from Scint. process, calculate the filter eff. based on its energy
        // If  (random number) > (Eff), the photon is rejected by yellow filter (= do not make hit (or count photon) for that photon)
        dd4hep::DDSegmentation::VolumeID ceren = static_cast<dd4hep::DDSegmentation::VolumeID>(m_segmentation->decoder()->get(cID, "c"));
        bool IsCeren = static_cast<bool>(ceren);

        Geant4Event& evt = context()->event();
        Geant4Random& rnd = evt.random(); // Get random generator
        double rndVal = rnd.uniform(0, 1); // Get random number from uniform distribution [0, 1]

        if (!IsCeren) {
          if (m_userData.rejectedByYellowFilter(energy, rndVal))
            return false;
        }

        // check wheter the optical photon will be rejected by the SiPM
        // TODO migrate this to the digitization step
        if (m_userData.rejectedBySiPM(energy, rndVal))
          return false;

        if (!hit) {
          hit = new Geant4DRCalorimeter::Hit(m_userData.fWavlenStep, m_userData.fTimeStep);
          hit->cellID = cID;
          hit->SetSiPMnum(cID);
          hit->SetTimeStart(m_userData.fTimeStart);
          hit->SetTimeEnd(m_userData.fTimeEnd);
          hit->SetWavlenMax(m_userData.fWavlenStart);
          hit->SetWavlenMin(m_userData.fWavlenEnd);
          coll->add(cID,hit);
        }

        hit->photonCount();
        int wavBin = m_userData.findWavBin(energy);
        hit->CountWavlenSpectrum(wavBin);
        int timeBin = m_userData.findTimeBin(hitTime);
        hit->CountTimeStruct(timeBin);

        hit->position = glob;

        return true;
      } // !skipScint
    } // Geant4SensitiveAction::process

    typedef Geant4SensitiveAction<DRCData> DRCaloSDAction;
  }
}

#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE(DRCaloSDAction)

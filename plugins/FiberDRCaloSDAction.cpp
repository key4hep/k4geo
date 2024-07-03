// DD4hep Framework include files
#include "DD4hep/Version.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Segmentations.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDG4/Geant4Random.h"
// k4geo Framework include files

#include "FiberDRCaloSDAction.h"

// Geant4 include files
#include "G4THitsCollection.hh"
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
#include "G4VProcess.hh"

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
     *  \author  Sungwon Kim
     *  \version 1.0
     *  \ingroup DD4HEP_SIMULATION
     */

    struct DRCData
    {
      Geant4HitCollection *fHitCollection;
      G4int fHCID;

      G4int fWavBin;
      G4int fTimeBin;
      G4float fWavlenStart;
      G4float fWavlenEnd;
      G4float fTimeStart;
      G4float fTimeEnd;
      G4float fWavlenStep;
      G4float fTimeStep;
      static const int fArraySize = 24;
      double fGraph_X[fArraySize] ={1.37760 * dd4hep::eV,
                                    1.45864 * dd4hep::eV,
                                    1.54980 * dd4hep::eV,
                                    1.65312 * dd4hep::eV,
                                    1.71013 * dd4hep::eV,
                                    1.77120 * dd4hep::eV,
                                    1.83680 * dd4hep::eV,
                                    1.90745 * dd4hep::eV,
                                    1.98375 * dd4hep::eV,
                                    2.06640 * dd4hep::eV,
                                    2.10143 * dd4hep::eV,
                                    2.13766 * dd4hep::eV,
                                    2.17516 * dd4hep::eV,
                                    2.21400 * dd4hep::eV,
                                    2.25426 * dd4hep::eV,
                                    2.29600 * dd4hep::eV,
                                    2.33932 * dd4hep::eV,
                                    2.38431 * dd4hep::eV,
                                    2.43106 * dd4hep::eV,
                                    2.47968 * dd4hep::eV,
                                    2.53029 * dd4hep::eV,
                                    2.58300 * dd4hep::eV,
                                    2.63796 * dd4hep::eV,
                                    2.69531 * dd4hep::eV};
        double fGraph_Y[fArraySize] ={0.903 * 0.903,
                                      0.903 * 0.903,
                                      0.903 * 0.903,
                                      0.903 * 0.903,
                                      0.903 * 0.903,
                                      0.903 * 0.903,
                                      0.902 * 0.902,
                                      0.901 * 0.901,
                                      0.898 * 0.898,
                                      0.895 * 0.895,
                                      0.893 * 0.893,
                                      0.891 * 0.891,
                                      0.888 * 0.888,
                                      0.883 * 0.883,
                                      0.87 * 0.87,
                                      0.838 * 0.838,
                                      0.76 * 0.76,
                                      0.62 * 0.62,
                                      0.488 * 0.488,
                                      0.345 * 0.345,
                                      0.207 * 0.207,
                                      0.083 * 0.083,
                                      0.018 * 0.018,
                                      0.};

      G4double wavToE(G4double wav) { return CLHEP::h_Planck * CLHEP::c_light / wav; }

      int findWavBin(G4double en)
      {
        int i = 0;
        for (; i < fWavBin + 1; i++)
        {
          if (en < wavToE((fWavlenStart - static_cast<float>(i) * fWavlenStep) * CLHEP::nm))
            break;
        }

        return fWavBin + 1 - i;
      }

      int findTimeBin(G4double stepTime)
      {
        int i = 0;
        for (; i < fTimeBin + 1; i++)
        {
          if (stepTime < ((fTimeStart + static_cast<float>(i) * fTimeStep) * CLHEP::ns))
            break;
        }

        return i;
      }

      // Linear interpolation function for calculating the efficiency of yellow filter, used in rejectedByYellowFilter
      double getFilterEff(G4double G4energy) {
        double energy = (G4energy * (dd4hep::MeV / CLHEP::MeV)); // Convert G4 MeV to dd4hep MeV
        if (energy <= fGraph_X[0]) // If the photon energy <= 1.37760 eV, than return maximum filter efficiency
          return fGraph_Y[0];
        
        for(int idx = 0; idx < fArraySize; ++idx) {
          if (energy <= fGraph_X[idx]) {
            double x1 = fGraph_X[idx - 1];
            double x2 = fGraph_X[idx];
            double y1 = fGraph_Y[idx - 1];
            double y2 = fGraph_Y[idx];

            return (y1 + ((y2 - y1) / (x2 - x1))*(energy - x1)); // return linear interpolated filter efficiency
          }
        }

        // This should not happen
        std::cout << "Error in Yellow filter efficiency with photon energy : " << energy << " MeV" << std::endl;
        std::cout << "Cannot find corresponding filter efficiency" << std::endl;
        return 0.;
      }

      // If true, then the photon is rejected by yellow filter
      bool rejectedByYellowFilter(G4double G4energy, double rndVal)
      {
        double energy = (G4energy * (dd4hep::MeV / CLHEP::MeV)); // Convert G4 MeV to dd4hep MeV
        if ( energy >= fGraph_X[fArraySize-1] ) // Photon w/ E > 2.69531 eV rejected
          return true;

        double FilterEff = getFilterEff(G4energy); // Get efficiency of filter using photon's energy

        // filter efficiency == probability of photon accepted by filter
        // == Probability of random value (from uniform distribution with range 0 ~ 1) smaller than filter efficiency
        // So if the rndVal is larger than the FilterEff, than the photon is rejected
        return (rndVal > FilterEff);
      }

      DRCData()
          : fHitCollection(0), fHCID(-1), fWavBin(120), fTimeBin(650),
            fWavlenStart(900.), fWavlenEnd(300.), fTimeStart(5.), fTimeEnd(70.)
      {
        fWavlenStep = (fWavlenStart - fWavlenEnd) / (float)fWavBin;
        fTimeStep = (fTimeEnd - fTimeStart) / (float)fTimeBin;
      }
    };

    /// Define collections created by this sensitive action object
    template <>
    void Geant4SensitiveAction<DRCData>::defineCollections()
    {
      std::string readout_name = m_sensitive.readout().name();
      m_collectionID = defineCollection<Geant4DRCalorimeter::Hit>(m_sensitive.readout().name());
      std::cout << "defineCollection Geant4DRCalorimeter readout_name   : " << readout_name << std::endl;
      std::cout << "defineCollection Geant4DRCalorimeter m_collectionID : " << m_collectionID << std::endl;
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <>
    G4bool Geant4SensitiveAction<DRCData>::process(G4Step GEANT4_CONST_STEP *step, G4TouchableHistory *history)
    {
      if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;

      typedef Geant4DRCalorimeter::Hit Hit;

      Geant4StepHandler h(step);
      Geant4HitCollection *coll = collection(m_collectionID);
      HitContribution contrib = Hit::extractContribution(step);

      auto theTouchable = step->GetPostStepPoint()->GetTouchable();
      dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
      dd4hep::VolumeID volID = volMgr.volumeID(theTouchable);
      G4ThreeVector global = step->GetPostStepPoint()->GetPosition();
      G4ThreeVector local = theTouchable->GetHistory()->GetTopTransform().TransformPoint(global);
      // MoveUpHistory(GetHistoryDepth - 2) -> tower touchable, local position
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
      if (!(IsCeren))
      {
        Geant4Event&  evt = context()->event();
        Geant4Random& rnd = evt.random(); // Get random generator
        double rndVal = rnd.uniform(0, 1); // Get random number from uniform distribution [0, 1]
        if (m_userData.rejectedByYellowFilter(energy, rndVal))
          return true;
      }

      if (!hit)
      {
        hit = new Geant4DRCalorimeter::Hit(m_userData.fWavlenStep, m_userData.fTimeStep);
        hit->cellID = cID;
        hit->SetSiPMnum(cID);
        hit->SetTimeStart(m_userData.fTimeStart);
        hit->SetTimeEnd(m_userData.fTimeEnd);
        hit->SetWavlenMax(m_userData.fWavlenStart);
        hit->SetWavlenMin(m_userData.fWavlenEnd);
        coll->add(hit);
      }

      hit->photonCount();
      int wavBin = m_userData.findWavBin(energy);
      hit->CountWavlenSpectrum(wavBin);
      int timeBin = m_userData.findTimeBin(hitTime);
      hit->CountTimeStruct(timeBin);

      hit->position = glob;
      hit->energyDeposit += contrib.deposit;
      hit->truth.emplace_back(contrib);

      return true;
    }

    typedef Geant4SensitiveAction<DRCData> DRCaloSDAction;
  }
}

#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE(DRCaloSDAction)
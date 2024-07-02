// DD4hep Framework include files
#include "DD4hep/Version.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Segmentations.h"
#include "DD4hep/DD4hepUnits.h"

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

      bool applyYellowFilter(G4double en)
      {
        double energy = (en * 1000000.); // MeV -> eV unit change

        if (energy >= 2.69531)
          return true; // Photon w/ E > 2.69531 eV filtered

        TF1 *uniform = new TF1("f1", "pol0", 0, 1); // TODO :: Pointer -> change to unique pointer
        uniform->SetParameter(0, 1);
        double randVal = uniform->GetRandom(0, 1); // TODO :: Change random engine to DDSim or G4 random!
        const int N = 24;
        double graph_X[N] = {1.37760,
                             1.45864,
                             1.54980,
                             1.65312,
                             1.71013,
                             1.77120,
                             1.83680,
                             1.90745,
                             1.98375,
                             2.06640,
                             2.10143,
                             2.13766,
                             2.17516,
                             2.21400,
                             2.25426,
                             2.29600,
                             2.33932,
                             2.38431,
                             2.43106,
                             2.47968,
                             2.53029,
                             2.58300,
                             2.63796,
                             2.69531};
        double graph_Y[N] = {0.903 * 0.903,
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
        TGraph *FilterEffGraph = new TGraph(N, graph_X, graph_Y); // TODO :: Pointer -> change to unique pointer
        double FilterEff = FilterEffGraph->Eval((double)energy);

        delete uniform;
        delete FilterEffGraph;

        // filter efficiency == probability of photon passing filter
        // == Probability of random value (from uniform distribution with range 0 ~ 1) larger than filter efficiency
        return (randVal > FilterEff);
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
      dd4hep::Position loc(local.x() * dd4hep::millimeter / CLHEP::millimeter, local.y() * dd4hep::millimeter / CLHEP::millimeter, local.z() * dd4hep::millimeter / CLHEP::millimeter);
      dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter, global.y() * dd4hep::millimeter / CLHEP::millimeter, global.z() * dd4hep::millimeter / CLHEP::millimeter);

      auto cID = m_segmentation->cellID(loc, glob, volID); // This returns cID corresponding to SiPM Wafer
      Hit *hit = coll->find<Hit>(CellIDCompare<Hit>(cID));

      G4double hitTime = step->GetPostStepPoint()->GetGlobalTime();
      G4double energy = step->GetTrack()->GetTotalEnergy();

      // Apply yellow filter here
      // Get random number from (0~1) uniform distribution
      // If the photon is from Scint. process, calculate the filter eff. based on its energy
      // If  (random number) > (Eff), apply yellow filter (= do not make hit (or count photon) for that photon)
      dd4hep::DDSegmentation::VolumeID ceren = static_cast<dd4hep::DDSegmentation::VolumeID>(m_segmentation->decoder()->get(cID, "c"));
      bool IsCeren = static_cast<bool>(ceren);
      if (!(IsCeren))
      {
        if (m_userData.applyYellowFilter(energy))
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
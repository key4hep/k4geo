//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "detectorSegmentations/SCEPCal_MainSegmentation_k4geo.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

namespace SCEPCal {
  class SCEPCal_MainSDAction {
    public:
      std::size_t  m_collectionID_scint {0}; 
      std::size_t  m_collectionID_ceren {0}; 
  };
}

namespace dd4hep {
  namespace sim {
    using namespace SCEPCal;
    
    template <> void Geant4SensitiveAction<SCEPCal_MainSDAction>::defineCollections()    {
      m_collectionID                     = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainEdep");
      m_userData.m_collectionID_scint    = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainScounts");
      m_userData.m_collectionID_ceren    = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainCcounts");
    }

    template <> bool 
    Geant4SensitiveAction<SCEPCal_MainSDAction>::process(const G4Step* step,G4TouchableHistory* /*hist*/) {
      G4StepPoint        *thePrePoint =step->GetPreStepPoint();
      G4TouchableHandle  thePreStepTouchable =thePrePoint->GetTouchableHandle();

      auto cellID =thePreStepTouchable->GetCopyNumber(0);
      G4Track *track =step->GetTrack();
      
      dd4hep::Segmentation *_geoSeg=&m_segmentation;
      auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCal_MainSegmentation_k4geo *>(_geoSeg->segmentation());

      DDSegmentation::Vector3D pos =segmentation->myPosition(cellID);
      Position global(pos.x(),pos.y(),pos.z());

      G4double edep =step->GetTotalEnergyDeposit();

      auto newOrExistingHitIn = [&](std::size_t id) {
        Geant4HitCollection* coll =collection(id);
        auto* hit =coll->findByKey<Geant4Calorimeter::Hit>(cellID);
        if(!hit) {
          hit =new Geant4Calorimeter::Hit(global);
          hit->cellID =cellID;
          coll->add(cellID, hit);
        }
        return hit;
      };

      // edep hits
      auto* hitedep =newOrExistingHitIn(m_collectionID);
      hitedep->energyDeposit+=edep;

      // Scintillation and Cerenkov hits
      if(track->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()) {
        auto procName =track->GetCreatorProcess()->GetProcessName();
        bool isCerenkov       =(procName =="CerenkovPhys");
        bool isScintillation  =(procName =="ScintillationPhys");
        if (!isCerenkov && !isScintillation) return true;

        if (track->GetCurrentStepNumber()==1) {
          auto* hitSC =newOrExistingHitIn(isScintillation? m_userData.m_collectionID_scint:
                                                           m_userData.m_collectionID_ceren);
          hitSC->energyDeposit+=(int)1000; // edep is in MeV, readout is in GeV, therefore a count of 1 is 1000
          track->SetTrackStatus(fStopAndKill);
        }
      }
      return true;
    }
  }
}

namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<SCEPCal_MainSDAction> SCEPCal_MainSDAction;
}}

DECLARE_GEANT4SENSITIVE(SCEPCal_MainSDAction)
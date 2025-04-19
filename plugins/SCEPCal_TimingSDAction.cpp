//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "detectorSegmentations/SCEPCal_TimingSegmentation_k4geo.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4OpticalPhoton.hh"

namespace SCEPCal {
  class SCEPCal_TimingSDAction {
    public:
      std::size_t  m_collectionID_scint {0}; 
      std::size_t  m_collectionID_ceren {0}; 
  };
}

namespace dd4hep {
  namespace sim {
    using namespace SCEPCal;
    
    template <> void Geant4SensitiveAction<SCEPCal_TimingSDAction>::defineCollections()    {
      m_collectionID                     = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_TimingEdep");
      m_userData.m_collectionID_scint    = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_TimingScounts");
      m_userData.m_collectionID_ceren    = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_TimingCcounts");
    }

    template <> bool 
    Geant4SensitiveAction<SCEPCal_TimingSDAction>::process(const G4Step* step,G4TouchableHistory* /*hist*/) {
      G4StepPoint        *thePrePoint =step->GetPreStepPoint();
      G4TouchableHandle  thePreStepTouchable =thePrePoint->GetTouchableHandle();

      auto cellID =thePreStepTouchable->GetCopyNumber(0);
      G4Track *track =step->GetTrack();
      
      dd4hep::Segmentation *_geoSeg=&m_segmentation;
      auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCal_TimingSegmentation_k4geo *>(_geoSeg->segmentation());

      DDSegmentation::Vector3D pos =segmentation->position(cellID);
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
        auto* creatorProc =track->GetCreatorProcess();
        if (creatorProc) {
          if (creatorProc->GetProcessType() == fOptical) {
            bool isCerenkov      = (creatorProc->GetProcessSubType() == fCerenkov);
            bool isScintillation = (creatorProc->GetProcessSubType() == fScintillation);
        
        
            if (!isCerenkov && !isScintillation) return true;

            if (track->GetCurrentStepNumber()==1) {
              auto* hitSC =newOrExistingHitIn(isScintillation? m_userData.m_collectionID_scint:
                                                              m_userData.m_collectionID_ceren);
              hitSC->energyDeposit+=1/dd4hep::MeV;
              track->SetTrackStatus(fStopAndKill);
            }
          }
        }
      }
      return true;
    }
  }
}

namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<SCEPCal_TimingSDAction> SCEPCal_TimingSDAction;
}}

DECLARE_GEANT4SENSITIVE(SCEPCal_TimingSDAction)
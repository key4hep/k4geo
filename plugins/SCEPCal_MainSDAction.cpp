//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "DDG4/Factories.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "G4EmProcessSubType.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessType.hh"
#include "G4VProcess.hh"
#include "detectorSegmentations/SCEPCal_MainSegmentation_k4geo.h"
#include "G4Poisson.hh"
#include "Randomize.hh"

// This function takes in input the ionizing energy deposit (in MeV) from a step in the main scepcal section
// and returns the corresponding number of sipm fired cells (or photo-electrons).
// It reproduces a scintillating light yield of 2000 p.e./GeV which is expected from test-beam data.
// The smearing is done by a single poissoninan sampling. De facto one must have two samplings:
// a poissonian sampling for light emission fluctuations and a Binomial sampling for light detection
// fluctuations. Thanks to the Poissonian thinning theorem this is equivalent to a single poissonian
// sampling reproducing on average the expected light yield.
int SmearSsignal(G4double edep){
   return G4Poisson((edep/1000.)*2000); // MeV->GeV
}

// This function is a Bernoully trial to decide wether a Cerenkov photon is detected or not
// in the scepcal. Cerenkov photons are created by Geant4 inside crystals, already taking
// into account their poissonian emission fluctuations, therefore we only need to add a
// Binomial sampling that reproduces the expected light yield of 100 p.e./GeV
bool SmearCsignal(){
  return G4UniformRand()<0.00147; 
}

namespace SCEPCal {
class SCEPCal_MainSDAction {
public:
  std::size_t m_collectionID_scint{0};
  std::size_t m_collectionID_ceren{0};
};
} // namespace SCEPCal

namespace dd4hep {
namespace sim {
  using namespace SCEPCal;

  template <>
  void Geant4SensitiveAction<SCEPCal_MainSDAction>::defineCollections() {
    m_collectionID = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainEdep");
    m_userData.m_collectionID_scint = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainScounts");
    m_userData.m_collectionID_ceren = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_MainCcounts");
  }

  template <>
  bool Geant4SensitiveAction<SCEPCal_MainSDAction>::process(const G4Step* step, G4TouchableHistory* /*hist*/) {
    G4StepPoint* thePrePoint = step->GetPreStepPoint();
    G4TouchableHandle thePreStepTouchable = thePrePoint->GetTouchableHandle();

    auto cellID = thePreStepTouchable->GetCopyNumber(0);
    G4Track* track = step->GetTrack();

    dd4hep::Segmentation* _geoSeg = &m_segmentation;
    auto segmentation = dynamic_cast<dd4hep::DDSegmentation::SCEPCal_MainSegmentation_k4geo*>(_geoSeg->segmentation());

    G4double edep = step->GetTotalEnergyDeposit();

    auto newOrExistingHitIn = [&](std::size_t id) {
      Geant4HitCollection* coll = collection(id);
      auto* hit = coll->findByKey<Geant4Calorimeter::Hit>(cellID);
      if (!hit) {
        DDSegmentation::Vector3D pos = segmentation->position(cellID);
        Position global(pos.x(), pos.y(), pos.z());

        hit = new Geant4Calorimeter::Hit(global / dd4hep::mm);
        hit->cellID = cellID;
        coll->add(cellID, hit);
      }
      return hit;
    };

    // Do not fill edep and S hits if there is no ionizing energy deposit
    //
    if(edep > 0.){
      // edep hits
      auto* hitedep = newOrExistingHitIn(m_collectionID);
      hitedep->energyDeposit += edep;

      // scintillation hits
      auto* hitS = newOrExistingHitIn(m_userData.m_collectionID_scint);
      auto Scount = SmearSsignal(edep);
      if(Scount > 0) hitS->energyDeposit += Scount / dd4hep::MeV;
    }

    // Cerenkov hits
    if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
      auto procSubType = track->GetCreatorProcess()->GetProcessSubType();
      bool isCerenkov = (procSubType == 21); // Cerenkov process sub type is int 21
      if (!isCerenkov){
        track->SetTrackStatus(fStopAndKill);
        return true;
      }
      else if (track->GetCurrentStepNumber() == 1) {
        auto* hitC = newOrExistingHitIn(m_userData.m_collectionID_ceren);
        if(SmearCsignal()) hitC->energyDeposit += 1 / dd4hep::MeV;
        track->SetTrackStatus(fStopAndKill);
      }
      else return true;
    }
    return true;
  }
} // namespace sim
} // namespace dd4hep

namespace dd4hep {
namespace sim {
  typedef Geant4SensitiveAction<SCEPCal_MainSDAction> SCEPCal_MainSDAction;
}
} // namespace dd4hep

DECLARE_GEANT4SENSITIVE(SCEPCal_MainSDAction)

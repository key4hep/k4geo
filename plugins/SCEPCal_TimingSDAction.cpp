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
#include "detectorSegmentations/SCEPCal_TimingSegmentation_k4geo.h"
#include "G4Poisson.hh"

// This function takes in input the ionizing energy deposit (in MeV) from a step in the scepcal timing
// layer and returns the corresponding number of sipm fired cells (or photo-electrons).
// It reproduces the expevted very-high LYSO light yield of 6000 p.e./MeV.
// The smearing is done by a single poissoninan sampling. De facto one must have two samplings:
// a poissonian sampling for light emission fluctuations and a Binomial sampling for light detection
// fluctuations. Thanks to the Poissonian thinning theorem this is equivalent to a single poissonian
// sampling reproducing on average the expected light yield.
int SmearTimingLayersignal(G4double edep){
   return G4Poisson((edep)*6000);
}

namespace SCEPCal {
class SCEPCal_TimingSDAction {
public:
  std::size_t m_collectionID_scint{0};
};
} // namespace SCEPCal

namespace dd4hep {
namespace sim {
  using namespace SCEPCal;

  template <>
  void Geant4SensitiveAction<SCEPCal_TimingSDAction>::defineCollections() {
    m_collectionID = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_TimingEdep");
    m_userData.m_collectionID_scint = defineCollection<Geant4Calorimeter::Hit>("SCEPCal_TimingScounts");
  }

  template <>
  bool Geant4SensitiveAction<SCEPCal_TimingSDAction>::process(const G4Step* step, G4TouchableHistory* /*hist*/) {
    G4StepPoint* thePrePoint = step->GetPreStepPoint();
    G4TouchableHandle thePreStepTouchable = thePrePoint->GetTouchableHandle();

    auto cellID = thePreStepTouchable->GetCopyNumber(0);

    dd4hep::Segmentation* _geoSeg = &m_segmentation;
    auto segmentation =
        dynamic_cast<dd4hep::DDSegmentation::SCEPCal_TimingSegmentation_k4geo*>(_geoSeg->segmentation());

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

    // edep hits
    if(edep > 0.) { // skip non-ionizing steps
      auto* hitedep = newOrExistingHitIn(m_collectionID);
      hitedep->energyDeposit += edep;
 
      auto* hitS = newOrExistingHitIn(m_userData.m_collectionID_scint);
      auto Scount = SmearTimingLayersignal(edep);
      if(Scount > 0) hitS->energyDeposit += Scount / dd4hep::MeV;
    }

    return true;
  }
} // namespace sim
} // namespace dd4hep

namespace dd4hep {
namespace sim {
  typedef Geant4SensitiveAction<SCEPCal_TimingSDAction> SCEPCal_TimingSDAction;
}
} // namespace dd4hep

DECLARE_GEANT4SENSITIVE(SCEPCal_TimingSDAction)

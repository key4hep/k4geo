//**************************************************************************
// \file DRTubesSDAction.cpp
// \brief: Implementation of Geant4SensitiveAction template class for
//         dual-readout-tubes calorimeters
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 11 August 2024
//**************************************************************************

// Includers from DD4HEP
#include "DD4hep/Segmentations.h"
#include "DD4hep/Version.h"
#include "DDG4/Factories.h"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4GeneratorAction.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4RunAction.h"
#include "DDG4/Geant4SensDetAction.inl"

// Includers from Geant4
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4Poisson.hh"
#include "G4ProcessManager.hh"
#include "G4Types.hh"
#include "G4UserSteppingAction.hh"
#include "G4VProcess.hh"
#include "globals.hh"

// Includers from project files
#include "DRTubesSglHpr.hh"

// #define DRTubesSDDebug

namespace dd4hep
{
namespace sim
{
class DRTubesSDData
{
    // Constructor and destructor
    //
  public:
    DRTubesSDData() = default;
    ~DRTubesSDData() = default;

    // Fields
    //
  public:
    Geant4Sensitive* sensitive{};
    int collection_cher_right;
    int collection_cher_left;
    int collection_scin_left;
};
}  // namespace sim
}  // namespace dd4hep

namespace dd4hep
{
namespace sim
{

// Function template specialization of Geant4SensitiveAction class.
// Define actions
template<>
void Geant4SensitiveAction<DRTubesSDData>::initialize()
{
  m_userData.sensitive = this;
}

// Function template specialization of Geant4SensitiveAction class.
// Define collections created by this sensitivie action object
template<>
void Geant4SensitiveAction<DRTubesSDData>::defineCollections()
{
  std::string ROname = m_sensitive.readout().name();
  m_collectionID = defineCollection<Geant4Calorimeter::Hit>(ROname + "ScinRight");
  m_userData.collection_cher_right = defineCollection<Geant4Calorimeter::Hit>(ROname + "CherRight");
  m_userData.collection_scin_left = defineCollection<Geant4Calorimeter::Hit>(ROname + "ScinLeft");
  m_userData.collection_cher_left = defineCollection<Geant4Calorimeter::Hit>(ROname + "CherLeft");
}

// Function template specialization of Geant4SensitiveAction class.
// Method that accesses the G4Step object at each track step.
template<>
bool Geant4SensitiveAction<DRTubesSDData>::process(const G4Step* aStep,
                                                         G4TouchableHistory* /*history*/)
{
  // NOTE: Here we do manipulation of the signal in each fiber (Scintillating and Cherenkov)
  // to compute the calorimeter signal and populate the corresponding hit.
  // Sensitive volumes are associated in the steering file using the DD4hep regexSensitiveDetector
  // and matching the substring DRETS. Usually DD4hep volIDs are retrieved with something like this
  // ---
  // dd4hep::BitFieldCoder
  // decoder("stave:10,tower:6,air:1,col:16,row:16,clad:1,core:1,cherenkov:1");
  // auto VolID = volumeID(aStep); auto CherenkovID = decoder.get(VolID,"cherenkov");
  // ---
  // However using regexSensitiveDetector does not populate the cache that allows using the
  // volumeID() method. One can still access volIDs in this sensitive detector action with something
  // like this
  // ---
  // G4VPhysicalVolume* PhysVol = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  // Geant4Mapping& Mapping = Geant4Mapping::instance();
  // PlacedVolume PlacedVol = Mapping.placement(PhysVol);
  // const PlacedVolumeExtension::VolIDs& TestIDs = PlacedVol.volIDs();
  // auto it = TestIDs.find("name");
  // std::cout<<it->first<<" "<<it->second<<std::endl;
  // ---
  // but this brute force method makes the simulaiton slower by more than two order of magnitudes
  // (see https://github.com/AIDASoft/DD4hep/issues/1319).
  // Therefore we use copynumbers instead of volIDs.

  // NOTE: in this SDAction we apply our custom Birks' Law correction
  // on the G4Step via the method DRTubesSglHpr::ApplyBirks().
  // However it should be known that the dd4hep Geant4StepHandler can apply Birks' Law
  // correction too (using internally the Geant4 class G4EmSaturation).
  // Therefore, if the Geant4StepHandler is used on this subdetector it might
  // lead to a double application of Birks' correction.

#ifdef DRTubesSDDebug
  // Print out some info step-by-step in sensitive volumes
  //
  DRTubesSglHpr::PrintStepInfo(aStep);
#endif

  auto Edep = aStep->GetTotalEnergyDeposit();
  auto cpNo = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
  // The second bit of the CopyNumber corresponds to the "core" entry:
  // 1 if the step is in the fiber core (S or C) and 0 if it is
  // in the fiber cladding (C only) 
  unsigned int CoreID = (cpNo & 0b10) >> 1; // take CpNo 2nd bit
  bool IsCherClad = (CoreID == 0);
  // The first bit of the CopyNumber corresponds to the "cherenkov" entry
  // 1 for C fibers and 0 for S fibers
  unsigned int CherenkovID = cpNo & 0b1;
  bool IsScin = (CherenkovID == 0);
  // Skip this step if edep is 0 and it is a scintillating fiber
  if (IsScin && Edep == 0.) return true;
  // If it is a track inside the cherenkov CLADDING skip this step,
  // if it is an optical photon kill it first
  if (IsCherClad) {
    if (aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition()) {
      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    }
    return true;
  }

  // Now we are inside fibers' core volume (either Scintillating or Cherenkov)

  // We recreate the TubeID from the tube copynumber:
  // fist16 bits for the columnID and second 16 bits for the rowID
  auto TubeID = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(2);
  unsigned int ColumnID = TubeID >> 16;
  unsigned int RawID = TubeID & 0xFFFF;
  auto TowerID = static_cast<unsigned int>(aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(3));
  auto StaveID = static_cast<unsigned int>(aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(4));

  VolumeID VolID = 0; // recreate the 64-bit VolumeID
  BitFieldCoder bc("system:5,stave:10,tower:6,air:1,col:16,row:16,clad:1,core:1,cherenkov:1");
  bc.set(VolID, "system", 25); // this number is set in DectDimensions_IDEA_o2_v01.xml
  bc.set(VolID, "stave" , StaveID);
  bc.set(VolID, "tower" , TowerID);
  bc.set(VolID, "air", 0);
  bc.set(VolID, "col", ColumnID);
  bc.set(VolID, "row", RawID);
  bc.set(VolID, "clad", 1);
  bc.set(VolID, "core", CoreID);
  bc.set(VolID, "cherenkov", CherenkovID);

  /* If you want to compare the 64-bits VolID created here
   * with the original DD4hep volumeID:
   * 1. set in DREndcapTubes_o1_v01.xml clad_C, core_C and core_S
   * volumes as sensitive
   * 2. associate DRTubesSDAction to DREncapTubes subdetector
   * in the steering file (instead of using RegexSD)
   * 3. Uncomment the code below */
  /*std::cout<<"Volume id, created "<<VolID<<" and DD4hep oroginal "<<volumeID(aStep)<<std::endl;
  std::cout<<"system id, created "<<25<<" and DD4hep original "<<bc.get(volumeID(aStep),"system")<<std::endl;
  std::cout<<"stave id, created "<<StaveID<<" and DD4hep original "<<bc.get(volumeID(aStep),"stave")<<std::endl;
  std::cout<<"tower id, created "<<TowerID<<" and DD4hep original "<<bc.get(volumeID(aStep),"tower")<<std::endl;
  std::cout<<"air id, created "<<0<<" and DD4hep original "<<bc.get(volumeID(aStep),"air")<<std::endl;
  std::cout<<"col id, created "<<ColumnID<<" and DD4hep original "<<bc.get(volumeID(aStep),"col")<<std::endl;
  std::cout<<"row id, created "<<RawID<<" and DD4hep original "<<bc.get(volumeID(aStep),"row")<<std::endl;
  std::cout<<"clad id, created "<<1<<" and DD4hep original "<<bc.get(volumeID(aStep),"clad")<<std::endl;
  std::cout<<"core id, created "<<CoreID<<" and DD4hep original "<<bc.get(volumeID(aStep),"core")<<std::endl;
  std::cout<<"cherenkov id, created "<<CherenkovID<<" and DD4hep original "<<bc.get(volumeID(aStep),"cherenkov")<<std::endl;*/

  bool IsRight = (aStep->GetPreStepPoint()->GetPosition().z() > 0.);

  // We now calculate the signal in S and C fiber according to the step contribution
  //
  G4double steplength = aStep->GetStepLength();
  G4int signalhit = 0;
  Geant4HitCollection* coll = (IsRight && IsScin)    ? collection(m_collectionID)
                              : (IsRight && !IsScin) ? collection(m_userData.collection_cher_right)
                              : (!IsRight && IsScin) ? collection(m_userData.collection_scin_left)
                                                     : collection(m_userData.collection_cher_left);

  if (IsScin) {  // it is a scintillating fiber

    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || steplength == 0.) {
      return true;  // not ionizing particle
    }
    G4double distance_to_sipm = DRTubesSglHpr::GetDistanceToSiPM(aStep);
    signalhit =
      DRTubesSglHpr::SmearSSignal(DRTubesSglHpr::ApplyBirks(Edep, steplength));
    signalhit = DRTubesSglHpr::AttenuateSSignal(signalhit, distance_to_sipm);
    if (signalhit == 0) return true;
  }  // end of scintillating fibre sigal calculation

  else {  // it is a Cherenkov fiber
    // calculate the signal in terms of Cherenkov photo-electrons
    if (aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition()) {
      G4OpBoundaryProcessStatus theStatus = Undefined;
      G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

      if (OpManager) {
        G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
        G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

        for (G4int i = 0; i < MAXofPostStepLoops; i++) {
          G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
          G4OpBoundaryProcess* fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
          if (fOpProcess) {
            theStatus = fOpProcess->GetStatus();
            break;
          }
        }
      }

      switch (theStatus) {
        case TotalInternalReflection: {
          // Kill Cherenkov photons inside fibers travelling towards the inner tip
          if (!DRTubesSglHpr::IsReflectedForward(aStep)) return true;
          G4double distance_to_sipm = DRTubesSglHpr::GetDistanceToSiPM(aStep);
          G4int c_signal = DRTubesSglHpr::SmearCSignal();
          signalhit = DRTubesSglHpr::AttenuateCSignal(c_signal, distance_to_sipm);
          if (signalhit == 0) return true;
          aStep->GetTrack()->SetTrackStatus(fStopAndKill);
          break;
        }
        default:
          aStep->GetTrack()->SetTrackStatus(fStopAndKill);
          return true;
      }  // end of swich cases
    }  // end of optical photon
    else {
      return true;
    }
  }  // end of Cherenkov fiber

  // We are going to create an hit per each fiber with a signal above 0
  // Each fiber is identified with a unique volID
  //
  Geant4Calorimeter::Hit* hit = coll->findByKey<Geant4Calorimeter::Hit>(VolID);  // the hit
  if (!hit) {  // if the hit does not exist yet, create it
    hit = new Geant4Calorimeter::Hit();
    hit->cellID = VolID;  // this should be assigned only once
    G4ThreeVector FiberVec = DRTubesSglHpr::CalculateFiberPosition(aStep);
    Position FiberPos(FiberVec.x(), FiberVec.y(), FiberVec.z());
    hit->position = FiberPos;  // this should be assigned only once
    // Note, when the hit is saved in edm4hep format the energyDeposit is
    // divided by 1000, i.e. it translates from MeV (Geant4 unit) to GeV (EDM4hep unit).
    // Here I am using this field to save photo-electrons, so I multiply it by 1000
    hit->energyDeposit = signalhit * 1000;
    coll->add(VolID, hit);  // add the hit to the hit collection
  }
  else {  // if the hit exists already, increment its fields
    hit->energyDeposit += signalhit * 1000;
  }
  return true;
}  // end of Geant4SensitiveAction::process() method specialization

}  // namespace sim
}  // namespace dd4hep

//--- Factory declaration
namespace dd4hep
{
namespace sim
{
typedef Geant4SensitiveAction<DRTubesSDData> DRTubesSDAction;
}
}  // namespace dd4hep
DECLARE_GEANT4SENSITIVE(DRTubesSDAction)

//**************************************************************************

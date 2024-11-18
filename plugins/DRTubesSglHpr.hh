//**************************************************************************
// \file DRTubesSglHpr.hh
// \brief:  definition of DRTubesSglHpr class
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 13 August 2024
//**************************************************************************

// Header-only utility class to store methods needed when computing
// the DREndcapTubes signals in scintillating and Cherenkov fibers

#ifndef DRTubesSglHpr_h
#  define DRTubesSglHpr_h 1

// Includers from Geant4
//
#  include "G4Tubs.hh"
#  include "globals.hh"

#  include "CLHEP/Units/SystemOfUnits.h"

class DRTubesSglHpr
{
  public:
    DRTubesSglHpr() = delete;

  private:
    // Fields
    //
    static constexpr G4double fk_B = 0.126;  // Birks constant
    static constexpr G4double fSAttenuationLength = 1000.0 * CLHEP::m;
    static constexpr G4double fCAttenuationLength = 1000.0 * CLHEP::m;

  public:
    // Methods
    //
    // Apply Briks Law
    static constexpr G4double ApplyBirks(const G4double& de, const G4double& steplength)
    {
      return (de / steplength) / (1 + fk_B * (de / steplength)) * steplength;
    }

    // Smear S signal according to Poissonian fluctuations and light yield
    static G4int SmearSSignal(const G4double& satde) { return G4Poisson(satde * 9.5); }

    // Smear C signal according to Poissonian fluctuations and light yield
    static G4int SmearCSignal() { return G4Poisson(0.177); }

    // Calculate distance from step in fiber to SiPM
    inline static G4double GetDistanceToSiPM(const G4Step* step, bool prestep);
    static G4double GetDistanceToSiPM(const G4Step* step) { return GetDistanceToSiPM(step, true); }

    // Calculate how many photons survived after light attenuation
    inline static G4int AttenuateHelper(const G4int& signal, const G4double& distance,
                                        const G4double& attenuation_length);

    // Attenuate light in fibers
    static G4int AttenuateSSignal(const G4int& signal, const G4double& distance)
    {
      return AttenuateHelper(signal, distance, fSAttenuationLength);
    }

    static G4int AttenuateCSignal(const G4int& signal, const G4double& distance)
    {
      return AttenuateHelper(signal, distance, fCAttenuationLength);
    }

    inline static G4ThreeVector CalculateFiberPosition(const G4Step* step);

    // Check if photon is travelling towards SiPM
    inline static bool IsReflectedForward(const G4Step* step);

    // Print step info for debugging purposes
    inline static void PrintStepInfo(const G4Step* aStep);
};

inline G4double DRTubesSglHpr::GetDistanceToSiPM(const G4Step* step, bool prestep)
{
  // Get the pre-step point
  const G4StepPoint* StepPoint = prestep ? step->GetPreStepPoint() : step->GetPostStepPoint();
  // Get the global position of the step point
  G4ThreeVector globalPos = StepPoint->GetPosition();
  // Get the local position of the step point in the current volume's coordinate system
  G4ThreeVector localPos =
    StepPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(globalPos);

  // Get the logical volume of the current step
  G4LogicalVolume* currentVolume = StepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  // Get the solid associated with the logical volume
  G4Tubs* solid = dynamic_cast<G4Tubs*>(currentVolume->GetSolid());
  // Get the dimensions of the solid (size of the volume)
  G4double size = solid->GetZHalfLength();

  G4double distance_to_sipm = size - localPos.z();
  return distance_to_sipm;
}

inline G4int DRTubesSglHpr::AttenuateHelper(const G4int& signal, const G4double& distance,
                                                  const G4double& attenuation_length)
{
  double probability_of_survival = exp(-distance / attenuation_length);

  G4int survived_photons = 0;
  for (int i = 0; i < signal; i++) {
    // Simulate drawing between 0 and 1 with probability x of getting 1
    if (G4UniformRand() <= probability_of_survival) survived_photons++;
  }
  return survived_photons;
}

inline G4ThreeVector DRTubesSglHpr::CalculateFiberPosition(const G4Step* step)
{
  G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
  G4ThreeVector origin(0., 0., 0.);
  G4ThreeVector zdir(0., 0., 1.);
  G4ThreeVector vectPos =
    theTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);
  G4ThreeVector direction =
    theTouchable->GetHistory()->GetTopTransform().Inverse().TransformAxis(zdir);

  // Get the logical volume of the current step
  G4LogicalVolume* currentVolume =
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  // Get the solid associated with the logical volume
  G4Tubs* solid = dynamic_cast<G4Tubs*>(currentVolume->GetSolid());
  // Get the dimensions of the solid (size of the volume)
  G4double size = solid->GetZHalfLength();
  G4double lengthfiber = size * 2.;
  G4ThreeVector Halffibervect = direction * lengthfiber / 2;
  // Fibre tip position
  G4ThreeVector vectPostip = vectPos - Halffibervect;
  // SiPM position
  // G4ThreeVector SiPMvecPos = vectPos + Halffibervect;

  return vectPostip;
}

// Check if photon is travelling towards SiPM
// If a (Cherenkov) photon is internally reflected in fibers
// it might travel towards the SiPM or the inner finer tip.
// As we do not consider reflections at the inner fiber tip
// photons travelling backwards should be killed.
inline bool DRTubesSglHpr::IsReflectedForward(const G4Step* step)
{
  double PreStepDistance = GetDistanceToSiPM(step, true);
  double PostStepDistance = GetDistanceToSiPM(step, false);
  bool IsReflectedForward = (PostStepDistance < PreStepDistance) ? true : false;
  return IsReflectedForward;
}

inline void DRTubesSglHpr::PrintStepInfo(const G4Step* aStep)
{
  std::cout << "-------------------------------" << std::endl;
  std::cout << "--> DREndcapTubes: track info: " << std::endl;
  std::cout << "----> Track #: " << aStep->GetTrack()->GetTrackID() << " "
            << "Step #: " << aStep->GetTrack()->GetCurrentStepNumber() << " "
            << "Volume: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
            << " " << std::endl;
  std::cout << "--> DREndcapTubes:: position info(mm): " << std::endl;
  std::cout << "----> x: " << aStep->GetPreStepPoint()->GetPosition().x()
            << " y: " << aStep->GetPreStepPoint()->GetPosition().y()
            << " z: " << aStep->GetPreStepPoint()->GetPosition().z() << std::endl;
  std::cout << "--> DREndcapTubes: particle info: " << std::endl;
  std::cout << "----> Particle " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName()
            << " "
            << "Dep(MeV) " << aStep->GetTotalEnergyDeposit() << " "
            << "Mat " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " "
            << "Vol " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
            << " "
            << "CpNo " << aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber() << " "
            << "CpNo1 " << aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(1) << " "
            << "CpNo2 " << aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(2) << " "
            << "CpNo3 " << aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(3) << " "
            << "CpNo4 " << aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(4) << std::endl;
}

#endif  // DRTubesSglHpr_h

//**************************************************************************

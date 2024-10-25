//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include <DD4hep/Printout.h>
#include <DD4hep/InstanceCount.h>
#include <DDG4/Geant4Data.h>
// #include <DDG4/Geant4StepHandler.h>
// #include <DDG4/Geant4FastSimHandler.h>
// #include <G4Step.hh>
// #include <G4Allocator.hh>
// #include <G4OpticalPhoton.hh>
#include "detectorSegmentations/DRCrystalHit.h"
// #include "G4Track.hh"

SCEPCal::DRCrystalHit::DRCrystalHit()
: dd4hep::sim::Geant4HitData(), position(), truth(), energyDeposit(0), nCerenkovProd(0), nScintillationProd(0), tAvgC(0), tAvgS(0) {

  dd4hep::InstanceCount::increment(this);

}

SCEPCal::DRCrystalHit::DRCrystalHit(const Position& pos)
: dd4hep::sim::Geant4HitData(), position(pos), truth(), energyDeposit(0), nCerenkovProd(0), nScintillationProd(0), tAvgC(0), tAvgS(0) {

  dd4hep::InstanceCount::increment(this);

}

SCEPCal::DRCrystalHit::~DRCrystalHit() {
  dd4hep::InstanceCount::decrement(this);
}
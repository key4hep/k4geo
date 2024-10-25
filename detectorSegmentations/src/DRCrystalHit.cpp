//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include <DD4hep/Printout.h>
#include <DD4hep/InstanceCount.h>
#include <DDG4/Geant4Data.h>
#include <DDG4/Geant4StepHandler.h>
#include <DDG4/Geant4FastSimHandler.h>
#include <G4Step.hh>
#include <G4Allocator.hh>
#include <G4OpticalPhoton.hh>
#include "detectorSegmentations/DRCrystalHit.h"
#include "G4Track.hh"

using namespace dd4hep::sim;
using namespace dd4hep;
using namespace std;
using namespace SCEPCal;

DRCrystalHit::DRCrystalHit()
: Geant4HitData(), position(), truth(), energyDeposit(0), eta(0), phi(0), depth(0), system(0), ncerenkov(0), nscintillator(0) {

  InstanceCount::increment(this);

  for( int i=0; i<nsipmbins; i++) {
    nwavelen_cer[i]=0;
    nwavelen_scint[i]=0;
    ntime_cer[i]=0;
    ntime_scint[i]=0;
  }
}

DRCrystalHit::DRCrystalHit(const Position& pos)
: Geant4HitData(), position(pos), truth(), energyDeposit(0), eta(0), phi(0), depth(0), system(0), ncerenkov(0), nscintillator(0) {

  InstanceCount::increment(this);

  for( int i=0; i<nsipmbins; i++) {
    nwavelen_cer[i]=0;
    nwavelen_scint[i]=0;
    ntime_cer[i]=0;
    ntime_scint[i]=0;
  }
}

DRCrystalHit::~DRCrystalHit() {
  InstanceCount::decrement(this);
}
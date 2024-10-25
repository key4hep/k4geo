//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#ifndef DRCrystalHit_h
#define DRCrystalHit_h 1
// #include "G4VHit.hh"
// #include "G4THitsCollection.hh"
// #include "G4Allocator.hh"
// #include "G4ThreeVector.hh"
#include "DDG4/Geant4Data.h"
// #include "G4OpticalPhoton.hh"
// #include "G4VProcess.hh"
// #include "DD4hep/Objects.h"
// #include "DD4hep/Segmentations.h"
// #include "CLHEP/Vector/ThreeVector.h"

namespace SCEPCal {

  typedef ROOT::Math::XYZVector Position;

    class DRCrystalHit : public dd4hep::sim::Geant4HitData {

      public:
        typedef dd4hep::sim::Geant4HitData base_t;

        Position      position;
        Contributions truth;
        double        energyDeposit;
        
        int           nCerenkovProd;
        int           nScintillationProd;
        double        tAvgC;
        double        tAvgS;
        
      public:
        DRCrystalHit();
        DRCrystalHit(DRCrystalHit&& c) = delete;
        DRCrystalHit(const DRCrystalHit& c) = delete;
        DRCrystalHit(const Position& cell_pos);
        virtual ~DRCrystalHit();
        DRCrystalHit& operator=(DRCrystalHit&& c) = delete;
        DRCrystalHit& operator=(const DRCrystalHit& c) = delete;
    };
};

#endif

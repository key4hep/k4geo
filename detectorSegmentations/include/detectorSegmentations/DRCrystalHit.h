//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#ifndef DRCrystalHit_h
#define DRCrystalHit_h 1
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "DDG4/Geant4Data.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "DD4hep/Objects.h"
#include "DD4hep/Segmentations.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace SCEPCal {

  typedef ROOT::Math::XYZVector Position;

    const int nsipmbins=6000;

    class DRCrystalHit : public dd4hep::sim::Geant4HitData {

      public:
        typedef dd4hep::sim::Geant4HitData base_t;

        Position      position;
        Contributions truth;
        double        energyDeposit;

        int eta;
        int phi;
        int depth;
        int system;
        
        int ncerenkov;
        int nscintillator;

        int nbins=nsipmbins;
        
        float wavelen_min=300;
        float wavelen_max=1000;

        float time_min=0;
        float time_max=300;
        
        std::array<int,nsipmbins> nwavelen_cer;
        std::array<int,nsipmbins> nwavelen_scint;
        
        std::array<int,nsipmbins> ntime_cer;
        std::array<int,nsipmbins> ntime_scint;
        
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

#if defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__) || defined(__ROOTCLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ namespace dd4hep;
#pragma link C++ namespace dd4hep::sim;
#pragma link C++ namespace SCEPCal;
#pragma link C++ class     SCEPCal::DRCrystalHit+;
#endif
#endif

#ifndef FiberDRCaloSDAction_H
#define FiberDRCaloSDAction_H

// DD4hep Framework include files
#include "DD4hep/Segmentations.h"
#include "DD4hep/Version.h"

#include "DDG4/Geant4Data.h"

#include <G4Types.hh>

#include <map>

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

/// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
namespace sim {

  class Geant4DRCalorimeter {
  public:
    class Hit : public Geant4HitData {
    public:
      typedef Geant4HitData base_t;
      typedef std::map<int, int> DRsimTimeStruct;
      typedef std::map<int, int> DRsimWavlenSpectrum;

      /// Hit position
      ROOT::Math::XYZVector position;
      /// Hit contributions by individual particles
      std::vector<MonteCarloContrib> truth;
      /// Total energy deposit
      double energyDeposit;

    public:
      /// Default constructor (for ROOT)
      Hit();
      // Constructor
      Hit(float wavSampling, float timeSampling)
          : Geant4HitData(), fSiPMnum(0), fPhotons(0), mWavSampling(wavSampling), mTimeSampling(timeSampling) {}
      /// copy constructor
      Hit(const Hit& right) : Geant4HitData() {
        fSiPMnum = right.fSiPMnum;
        fPhotons = right.fPhotons;
        fWavlenSpectrum = right.fWavlenSpectrum;
        fTimeStruct = right.fTimeStruct;
        mWavSampling = right.mWavSampling;
        mTimeSampling = right.mTimeSampling;
      }
      /// Copy assignment operator
      Hit& operator=(Hit& right) {
        fSiPMnum = right.fSiPMnum;
        fPhotons = right.fPhotons;
        fWavlenSpectrum = right.fWavlenSpectrum;
        fTimeStruct = right.fTimeStruct;
        mWavSampling = right.mWavSampling;
        mTimeSampling = right.mTimeSampling;
        return *this;
      }
      G4bool operator==(const Hit& right) const { return (fSiPMnum == right.fSiPMnum); }
      /// Move constructor
      Hit(Hit&& c) = delete;

      /// Standard constructor
      Hit(const Position& cell_pos);
      /// Default destructor
      virtual ~Hit() {};
      /// Move assignment operator
      Hit& operator=(Hit&& c) = delete;

      void photonCount() { fPhotons++; }
      unsigned long GetPhotonCount() const { return fPhotons; }

      void SetSiPMnum(dd4hep::DDSegmentation::CellID n) { fSiPMnum = n; }
      const dd4hep::DDSegmentation::CellID& GetSiPMnum() const { return fSiPMnum; }

      void CountWavlenSpectrum(int ibin) {
        auto it = fWavlenSpectrum.find(ibin);

        if (it == fWavlenSpectrum.end())
          fWavlenSpectrum.insert(std::make_pair(ibin, 1));
        else
          it->second++;
      };
      const DRsimWavlenSpectrum& GetWavlenSpectrum() const { return fWavlenSpectrum; }

      void CountTimeStruct(int ibin) {
        auto it = fTimeStruct.find(ibin);

        if (it == fTimeStruct.end())
          fTimeStruct.insert(std::make_pair(ibin, 1));
        else
          it->second++;
      };
      const DRsimTimeStruct& GetTimeStruct() const { return fTimeStruct; }

      float GetSamplingTime() const { return mTimeSampling; }
      float GetSamplingWavlen() const { return mWavSampling; }

      float GetTimeStart() const { return fTimeStart; }
      void SetTimeStart(float val) { fTimeStart = val; }

      float GetTimeEnd() const { return fTimeEnd; }
      void SetTimeEnd(float val) { fTimeEnd = val; }

      float GetWavlenMax() const { return fWavlenMax; }
      void SetWavlenMax(float val) { fWavlenMax = val; }

      float GetWavlenMin() const { return fWavlenMin; }
      void SetWavlenMin(float val) { fWavlenMin = val; }

    private:
      dd4hep::DDSegmentation::CellID fSiPMnum;
      unsigned long fPhotons;
      DRsimWavlenSpectrum fWavlenSpectrum;
      DRsimTimeStruct fTimeStruct;
      float mWavSampling;
      float mTimeSampling;
      float fWavlenMax;
      float fWavlenMin;
      float fTimeStart;
      float fTimeEnd;
    };
  };

} // namespace sim
} // namespace dd4hep

#endif

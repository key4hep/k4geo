// Framework include files
#include <DDG4/Geant4FastSimShowerModel.inl.h>

#include "G4FastStep.hh"

// C/C++ include files
#include "DRCaloFastSimModel.h"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep
{
  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim
  {
    template <>
    void Geant4FSShowerModel<DRCFiberModel>::initialize()
    {
      this->m_applicablePartNames.emplace_back("opticalphoton");
    }

    template <>
    void Geant4FSShowerModel<DRCFiberModel>::constructSensitives(Geant4DetectorConstructionContext *ctxt)
    {
      this->Geant4FastSimShowerModel::constructSensitives(ctxt);
    }

    template <>
    void Geant4FSShowerModel<DRCFiberModel>::modelShower(const G4FastTrack &fasttrack, G4FastStep &faststep)
    {
      auto *track = fasttrack.GetPrimaryTrack();

      if (locals.fKill)
      { // absorption
        faststep.ProposeTotalEnergyDeposited(track->GetKineticEnergy());
        faststep.KillPrimaryTrack();

        return;
      }

      if (locals.fTransported)
        return; // reset NILL if the track did not meet NILL check

      double timeUnit = locals.mDataCurrent.globalTime - locals.mDataPrevious.globalTime;
      auto posShift = locals.mTransportUnit * locals.mNtransport * locals.mFiberAxis; // #TODO apply shift for xy direction as well
      double timeShift = timeUnit * locals.mNtransport;

      faststep.ProposePrimaryTrackFinalPosition(track->GetPosition() + posShift, false);
      faststep.ProposePrimaryTrackFinalTime(track->GetGlobalTime() + timeShift);
      faststep.ProposePrimaryTrackFinalKineticEnergy(track->GetKineticEnergy());
      faststep.ProposePrimaryTrackFinalMomentumDirection(track->GetMomentumDirection(), false);
      faststep.ProposePrimaryTrackFinalPolarization(track->GetPolarization(), false);
      locals.fTransported = true;
      return;
    }

    template <>
    bool Geant4FSShowerModel<DRCFiberModel>::check_applicability(const G4ParticleDefinition &particle)
    {
      return &particle == G4OpticalPhoton::OpticalPhotonDefinition();
    }

    template <>
    bool Geant4FSShowerModel<DRCFiberModel>::check_trigger(const G4FastTrack &fasttrack) {
      if (!locals.fSwitch)
        return false; // turn on/off the model

      const G4Track *track = fasttrack.GetPrimaryTrack();

      // reset when moving to the next track
      if (locals.mDataCurrent.trackID != track->GetTrackID())
        locals.reset();

      // make sure that the track does not get absorbed after transportation
      // as number of interaction length left is reset when doing transportation
      if (!locals.checkNILL())
        return true; // track is already transported but did not pass NILL check, attempt to reset NILL

      return locals.check_trigger(track);
    }

    typedef Geant4FSShowerModel<DRCFiberModel> Geant4DRCFiberModel;
  }
}

#include <DDG4/Factories.h>

DECLARE_GEANT4ACTION_NS(dd4hep::sim, Geant4DRCFiberModel)

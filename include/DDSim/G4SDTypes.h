//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//====================================================================
#ifndef G4SDTypes_h
#define G4SDTypes_h 1

#include "G4lcioSimHit.h"

#include "DDG4/Geant4SensDetAction.h"
#include "DDG4/Geant4Data.h"
#include "DDG4/Geant4StepHandler.h"

#include "DD4hep/Printout.h"
#include "DD4hep/InstanceCount.h"

namespace dd4hep {
  using namespace DD4hep::Simulation  ;
  using namespace DD4hep ;
}

namespace DDSim{
  
  
  /** Simple Tracker type - wrapping SimTrkHit ... why is this needed ?
   */
  struct DDSimTracker{
    //    typedef SimTrkHit Hit ;
    typedef lcio::SimTrackerHitImpl Hit ;
  } ;
  
  // copied from Geant4SDActions.cpp (why is this not a public class ??????)
  
  /** Simple SensitiveAction class ...
   */
  template <typename T> class Geant4SensitiveAction : public dd4hep::Geant4Sensitive  {
  protected:
    /// Collection identifiers
    size_t m_collectionID;
    
  public:
    //    typedef SimpleHit::Contribution HitContribution;
    // Standard , initializing constructor
    Geant4SensitiveAction(dd4hep::Geant4Context* context, const std::string& name, DetElement det, LCDD& lcdd)
      : Geant4Sensitive(context,name,det,lcdd), m_collectionID(0) {
      
      defineCollections();
      dd4hep::InstanceCount::increment(this);
    }
    
    /// Default destructor
    virtual ~Geant4SensitiveAction(){
      dd4hep::InstanceCount::decrement(this);
    }
    
    /// Define collections created by this sensitivie action object
    virtual void defineCollections() {}

    /// G4VSensitiveDetector interface: Method invoked at the begining of each event.
    virtual void begin(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::begin(hce);
    }
    
    /// G4VSensitiveDetector interface: Method invoked at the end of each event.
    virtual void end(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::end(hce);
    }
    
    /// G4VSensitiveDetector interface: Method for generating hit(s) using the G4Step object.
    virtual bool process(G4Step* step,G4TouchableHistory* history)  {
      
      return dd4hep::Geant4Sensitive::process(step,history);
    }
    
    /// G4VSensitiveDetector interface: Method invoked if the event was aborted.
    virtual void clear(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::clear(hce);
    }
    
  };
  

} // namespace


#endif

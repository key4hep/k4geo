//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  A.Sailer, CERN
//====================================================================
#ifndef DDBeamCalSD_h
#define DDBeamCalSD_h 1


#include "DDG4/Geant4SensDetAction.h"
// why not: #include "DDG4/Geant4Sensitive.h" ???

#include "DD4hep/InstanceCount.h"

class G4Step;

namespace IMPL {
  class SimCalorimeterHitImpl;
}

namespace dd4hep {
  using namespace DD4hep::Simulation  ;
  using namespace DD4hep::DDSegmentation ;
  using namespace DD4hep ;
  using DD4hep::VolumeID ;
}


namespace DDSim {

  class DDBeamCalSD : public dd4hep::Geant4Sensitive {
  
  protected:
    size_t m_collectionID;
    bool _detailedHitsStoring ;

  public:
    //======================================================================================================
  
    DDBeamCalSD(dd4hep::Geant4Context* mContext, const std::string& mName, DetElement mDet, LCDD& mLcdd):
      Geant4Sensitive(mContext,mName,mDet,mLcdd), m_collectionID(-1), _detailedHitsStoring(false)
    {
      declareProperty("detailedHitsStoring", _detailedHitsStoring =0) ;
      defineCollections();
      dd4hep::InstanceCount::increment(this);
    }
  
    /// Default destructor
    virtual ~DDBeamCalSD(){
      dd4hep::InstanceCount::decrement(this);
    }
  
    /// Define collections created by this sensitive action object
    virtual void defineCollections() {
      
      m_collectionID = dd4hep::Geant4Sensitive::defineCollection<IMPL::SimCalorimeterHitImpl>( name() );
  }

    /// G4VSensitiveDetector interface: Method invoked at the begining of each event.
    virtual void begin(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::begin(hce);
    }
  
    /// G4VSensitiveDetector interface: Method invoked at the end of each event.
    virtual void end(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::end(hce);
    }
  
    /// G4VSensitiveDetector interface: Method for generating hit(s) using the G4Step object.
    virtual bool process(G4Step* step,G4TouchableHistory* history) ;
  
    /// G4VSensitiveDetector interface: Method invoked if the event was aborted.
    virtual void clear(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::clear(hce);
    }
  
    //======================================================================================================

  };

} // namespace
#endif


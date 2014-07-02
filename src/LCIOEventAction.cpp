#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4OutputAction.h"

#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCEventImpl.h"

using namespace DD4hep ;
using namespace DD4hep::Simulation ;


namespace DDSim {
  
  /** Geant4EventAction for dealing with LCIO.
   *  Only functionality for now is adding the MCParticle
   *  collection (created in LCIOGeneratorAction) to the
   *  event.
   * @author F.Gaede, DESY
   * @version $Id:$
   */
  class LCIOEventAction : public Geant4EventAction {
    
  public:
    /// Standard constructor
    LCIOEventAction(Geant4Context* context, const std::string& nam) : Geant4EventAction( context , nam ) {}
    
    /// Default destructor
    virtual ~LCIOEventAction() {  }
    
    /// Begin-of-event callback
    virtual void begin(const G4Event* event){}
    
    /// End-of-event callback
    virtual void end(const G4Event* event);
  };
  
  //=================================================================================================
  
  void LCIOEventAction::end(const G4Event* event){
    
    // get the MCParticle collection from the event:
    EVENT::LCCollection* col  = context()->event().extension<EVENT::LCCollection>(false); 
   
    if( col != 0 ) {
      //      std::cout << " ########## LCIOEventAction::end() : found MCParticle collection " << col <<  std::endl ;

     lcio::LCEventImpl* evt = context()->event().extension<lcio::LCEventImpl>();

     if( evt) {

       evt->addCollection( col , "MCParticle" ) ;

     } else { 

       std::cout << " ########## LCIOEventAction::saveEvent() : invalid event pointer in context()->event().extension<lcio::LCEventImpl>()" << std::endl ;
     }

    } else {
      std::cout << " ########## LCIOEventAction::end() : MCParticle collection not found in event context " << col <<  std::endl ;
    }

  }
  
} // end namespace

//------------------------------------------------------------------------------------------
// need to make this class visible in DD4hep namespace:
namespace DD4hep{
  namespace Simulation{
    typedef DDSim::LCIOEventAction LCIOEventAction ;
  }
}
#include "DDG4/Factories.h"
// to use this macro
DECLARE_GEANT4ACTION(LCIOEventAction)


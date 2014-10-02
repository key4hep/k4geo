//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id:$
//====================================================================

#ifndef DDSim_LCIOGeneratorAction_h
#define DDSim_LCIOGeneratorAction_h

// Framework include files
#include "DDG4/Geant4GeneratorAction.h"

namespace dd4hep {
  using namespace DD4hep::Simulation  ;
  //  using namespace DD4hep::DDSegmentation ;
  using namespace DD4hep ;
}


namespace IO{
  class LCReader ;
}
namespace EVENT{
  class MCParticle ;
}

class G4PrimaryParticle ;

#include <string>
#include <set>
// Forward declarations
// class G4ParticleDefinition;
// class G4ParticleGun;

namespace DDSim {

  /** @class LCIOGeneratorAction LCIOGeneratorAction.h DDG4/LCIOGeneratorAction.h
   *
   * Implementation of a Geant4GeneratorAction that uses the lcio::MCParticle collection
   * to generate G4PrimaryParticles. Ported from Mokka::HepLCIOInterfaceNew.
   *
   * @parameter LCIOFileName  input file name 
   * @parameter MCParticleCollectionName input collection name - default: "MCParticle"
   * @parameter PrimaryVertexSpreadZ  spread of z-vertex position of primary vertez in mm - default  0.0 );
   * @parameter LorentzTransformationAngle  - crossing angle for lorentz boost  - default 0.0 
   * @parameter SkipNEvents number of events to be skipped from input file - default  0.0  );
   *
   * @author F.Gaede, DESY
   * @version $Id:$
   */
  class LCIOGeneratorAction: public dd4hep::Geant4GeneratorAction {

    typedef std::map< EVENT::MCParticle*  , G4PrimaryParticle* > LCIO2Geant4Map ;

  protected:
    
    //---  properties ---
    std::string _MCParticleCollectionName ;
    std::string _LCIOFileName ;
    double _PrimaryVertexSpreadZ ;
    double _LorentzTransformationAngle ;
    int    _SkipNEvents ;
    //-------------------

    IO::LCReader* _rdr ;

    
  private:
    /// helper method to get the MCParticles relevant for the simulation (from Mokka::HepLCIOInterfaceNew)
     std::set<G4PrimaryParticle*> getRelevantParticles(EVENT::MCParticle* p) ;

 
    std::set<EVENT::MCParticle*> _visitedParticles ;

    LCIO2Geant4Map _map ;

  public:
    /// Standard constructor
    LCIOGeneratorAction(dd4hep::Geant4Context* context, const std::string& name);

    /// Default destructor
    virtual ~LCIOGeneratorAction();

    /// Callback to generate primary particles
    virtual void operator()(G4Event* event);
  };

}      // End namespace DD4Sim

#endif //DDSim_LCIOGeneratorAction_h 


#include "DDSim/G4SDTypes.h"

//include "gearimpl/Vector3D.h"

namespace DDSim {
  
  /// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ///               Geant4SensitiveAction<SimpleTracker>
  /// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Define collections created by this sensitivie action object
  template <> 
  void Geant4SensitiveAction<DDSimTracker>::defineCollections() {
    
    m_collectionID = dd4hep::Geant4Sensitive::defineCollection<DDSimTracker::Hit>(name()+"Hits");
    //m_collectionID = defineCollection<SimpleHit>(name()+"Hits");
  }
  
  /// Method for generating hit(s) using the information of G4Step object.
  template <> bool Geant4SensitiveAction<DDSimTracker>::process(G4Step* step,G4TouchableHistory* /*hist*/ ) {
    typedef DDSimTracker::Hit Hit;
    
    StepHandler h(step);
    
    dd4hep::Position prePos    = h.prePos();
    dd4hep::Position postPos   = h.postPos();
    dd4hep::Position direction = postPos - prePos;
    dd4hep::Position position  = mean_direction(prePos,postPos);
    double   hit_len   = direction.R();
    
    if (hit_len > 0) {
      
      double new_len = mean_length(h.preMom(),h.postMom())/hit_len;
      direction *= new_len/hit_len;
    }
    

#if 0
    Hit* g4lhit = new Hit ; 
    lcio::SimTrackerHitImpl* hit = g4lhit->lcio_hit() ;
#else
    Hit* hit = new Hit ;
#endif
    //    (h.track->GetTrackID(),
     // h.track->GetDefinition()->GetPDGEncoding(),
     // step->GetTotalEnergyDeposit(),
     // h.track->GetGlobalTime());
    
    // HitContribution contrib = Hit::extractContribution(step);

    long cellID = volumeID( step ) ;

    hit->setCellID0( cellID & 0xffffffff   ) ;
    hit->setCellID1( ( cellID >> 32 ) & 0xffffffff   ) ;


    printout(dd4hep::INFO,"DDSimTracker","%s> Add hit with deposit:%f  Pos:%f %f %f - cellID0: 0x%x cellID1: 0x%x",
	     c_name(),step->GetTotalEnergyDeposit(),position.X(),position.Y(),position.Z() , hit->getCellID0() ,hit->getCellID1() );

    double pos[3] ;
    pos[0] = position.x() ;
    pos[1] = position.y() ;
    pos[2] = position.z() ;
    hit->setPosition( pos  ) ;

    // hit->energyDeposit = contrib.deposit ;
    // hit->position      = position;
    // hit->momentum      = direction;
    // hit->length        = hit_len;

    collection(m_collectionID)->add(hit);
    return hit != 0;
  }



} // namespace


namespace DD4hep{
  namespace Simulation{
    typedef DDSim::Geant4SensitiveAction<DDSim::DDSimTracker> DDSimTrackerAction;
    
  }}

//using namespace DD4hep

#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(DDSimTrackerAction)

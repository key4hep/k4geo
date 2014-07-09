//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  A.Sailer, CERN
//====================================================================
#ifndef DDBeamCalSD_h
#define DDBeamCalSD_h 1

#include "G4Step.hh"

#include "DDG4/Geant4SensDetAction.h"
// why not: #include "DDG4/Geant4Sensitive.h" ???

// #include "DDG4/Geant4Data.h"
// #include "DDG4/Geant4StepHandler.h"
// #include "DD4hep/Printout.h"
#include "DD4hep/InstanceCount.h"

//#include "lcio.h"
#include "Exceptions.h"
#include "IMPL/SimCalorimeterHitImpl.h"

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
    //fg: original methods - to be replaced w/ DDG4 ones...
    // void Initialize(G4HCofThisEvent*HCE);
    // G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    // void EndOfEvent(G4HCofThisEvent*HCE);
  
    DDBeamCalSD(dd4hep::Geant4Context* context, const std::string& name, DetElement det, LCDD& lcdd)
      : Geant4Sensitive(context,name,det,lcdd), m_collectionID(0), currentCylinder(0), currentPID(-1), 
	currentPDG(-1),
	EntryPoint(0.,0.,0.), ExitPoint(0.,0.,0.),
	EntryMomentum (0.,0.,0.), ExitMomentum (0.,0.,0.),
	DepositedEnergy(0.), HitTime(0.), StepLength(0.) {
      
      declareProperty("detailedHitsStoring", _detailedHitsStoring =0) ;
      
      std::cout << "================================================================================"  << std::endl;
      std::cout << "================================================================================"  << std::endl;
      std::cout << "==============================Constructor======================================="  << std::endl;
      std::cout << "================================================================================"  << std::endl;
      std::cout << "================================================================================"  << std::endl;
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


  private:
    void StartNewHit(G4int aLayerNumber,
		     G4ThreeVector theEntryPoint,
		     G4ThreeVector theEntryMomentum,
		     G4int thePDG,
		     G4int theTrackPID);
  
    void UpdateHit(G4Step *aStep);

    void DumpHit(G4Step *aStep);

    void Clear();

    G4double Threshold;
    G4double PrimaryTPCCut;
    G4int HCID;
    G4int currentCylinder;
    G4int currentPID;
    G4int currentPDG;

    G4ThreeVector EntryPoint;
    G4ThreeVector ExitPoint;
    G4ThreeVector EntryMomentum;
    G4ThreeVector ExitMomentum;
    G4double DepositedEnergy;
    G4double HitTime ;
    G4double StepLength;

    G4bool detailedHitsStoring;
    G4bool defaultHitsStoring;

  };

} // namespace
#endif


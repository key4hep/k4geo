//====================================================================
//  lcgeo - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  SensitiveDetector driver for HcalBarrel
//
//--------------------------------------------------------------------
//  S.Lu, DESY
//====================================================================

#ifndef SDHcalBarrel_h
#define SDHcalBarrel_h 1

#include "DDG4/Geant4SensDetAction.h"

#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4EmSaturation.hh"

#include "DD4hep/InstanceCount.h"

//#include "lcio.h"
#include "Exceptions.h"
#include "IMPL/SimCalorimeterHitImpl.h"

 
/** @class SDHcalBarrel SDHcalBarrel.hh "SDHcalBarrel.hh" 
    \brief Sensitive detector of the HCAL barrel
*/


namespace dd4hep {
  using namespace DD4hep::Simulation  ;
  using namespace DD4hep::DDSegmentation ;
  using namespace DD4hep ;
  using DD4hep::VolumeID ;
}


namespace lcgeo {

  /** Sensitive detector for ILD HcalBarrel - ported from Mokka SDHcalBarrel 
   *  (http://polzope.in2p3.fr:8081/MOKKA)
   *
   *  @author: S.Lu, DESY, Jan 2014
   *  @version: $Id:$
   */

  class SDHcalBarrel : public dd4hep::Geant4Sensitive
  {

    /**Constructor
	@param cellSize vector to specify the cell size in each layer
	@param applyBirksLaw flag to indicate if Birks law should be applied or not
    */

  protected:
    size_t m_collectionID;
    std::vector<double> _cellSize;
    std::string _cellIDEncoding;
    G4bool _applyBirksLawFlag;

  public:

    SDHcalBarrel(dd4hep::Geant4Context* context, const std::string& name, DetElement det, LCDD& lcdd)
      : Geant4Sensitive(context,name,det,lcdd), m_collectionID(0), currentLayer(0), currentPID(-1), 
	currentPDG(-1),
	DepositedEnergy(0.), HitTime(0.) {
      
      declareProperty("cellSize", _cellSize );
      declareProperty("applyBirksLawFlag", _applyBirksLawFlag = true );
     
      defineCollections();
      dd4hep::InstanceCount::increment(this);
    }
  
    /**Destructor*/
    virtual ~SDHcalBarrel(){
      dd4hep::InstanceCount::decrement(this);
    }

    /** Define collections created by this sensitivie action object */
    virtual void defineCollections() {
      
      m_collectionID = dd4hep::Geant4Sensitive::defineCollection<IMPL::SimCalorimeterHitImpl>( name() );
    }
    
    
    /** G4VSensitiveDetector interface: Method invoked at the begining of each event. */
    virtual void begin(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::begin(hce);
    }
    
    /** G4VSensitiveDetector interface: Method invoked at the end of each event. */
    virtual void end(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::end(hce);
    }
  
    /** G4VSensitiveDetector interface: Method for generating hit(s) using the G4Step object. */
    virtual bool process(G4Step* step,G4TouchableHistory*) ;
    
    /** G4VSensitiveDetector interface: Method invoked if the event was aborted. */
    virtual void clear(G4HCofThisEvent* hce) {
      dd4hep::Geant4Sensitive::clear(hce);
    }
  

    //======================================================================================================
    
    void DrawAll();
    void PrintAll();
    

  private:
    
    void StartNewHit(G4int aLayerNumber,
		     G4double theDepositedEnergy,
		     G4int thePDG,
		     G4int theTrackPID);
  
    void UpdateHit(G4Step *aStep);

    void DumpHit(G4Step *aStep);

    void Clear();

    // port from mokka SD
    /**Apply Birks attenuation law (example given by Vladimir Ivantchenko, from CERN)*/
    G4double BirkAttenuation(const G4Step* aStep);

    G4int currentLayer;
    G4int currentPID;
    G4int currentPDG;
    
    G4double DepositedEnergy;
    G4double HitTime ;
    
    G4EmSaturation *emSaturation;


  };

}
#endif


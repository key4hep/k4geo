//
//====================================================================
//  lcgeo - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  SensitiveDetector driver for HcalBarrel
//
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id$
//====================================================================
//
#include "SDHcalBarrel.h"

#include "DDG4/Geant4SensDetAction.h"

#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"


#include "G4SDManager.hh"
#include "G4VSolid.hh"
#include "G4ios.hh"
#include "globals.hh"
#include <cassert>
#include <math.h> 

//****************************
//#define SDHcalBarrel_DEBUG
//****************************

#include "IMPL/SimCalorimeterHitImpl.h"
#include <IMPL/LCFlagImpl.h>
#include <EVENT/LCIO.h>

#include "DDG4/Geant4StepHandler.h"
#include "DDG4/Geant4VolumeManager.h"
#include "DDG4/Geant4Mapping.h"

#include "DDSegmentation/BitField64.h"

namespace lcgeo {


  void SDHcalBarrel::DrawAll()
  {
  } 

  void SDHcalBarrel::PrintAll()
  {
  } 

//==================================================================
//
//  BirkAttenuation: port from mokka SD
//
//==================================================================
  G4double SDHcalBarrel::BirkAttenuation(const G4Step* aStep)
  {
    G4double energyDeposition = aStep->GetTotalEnergyDeposit();
    G4double length = aStep->GetStepLength();
    G4double niel   = aStep->GetNonIonizingEnergyDeposit();
    const G4Track* track = aStep->GetTrack();
    const G4ParticleDefinition* particle = track->GetDefinition();
    const G4MaterialCutsCouple* couple = track->GetMaterialCutsCouple();
    
    G4double engyVis = emSaturation->VisibleEnergyDeposition(particle,
							     couple,
							     length,
							     energyDeposition,
							     niel);
    return engyVis;
  }



//==================================================================
//
//  DumpHit: 
//
//==================================================================
  void SDHcalBarrel::DumpHit(G4Step* aStep)
  {

    dd4hep::Geant4StepHandler stepH( aStep );
    dd4hep::Geant4VolumeManager volMgr = dd4hep::Geant4Mapping::instance().volumeManager();
    dd4hep::BitField64& decoder = *  m_readout.idSpec().decoder() ;
    dd4hep::VolumeID preID  = volMgr.volumeID( stepH.preTouchable() );

#ifdef SDHcalBarrel_DEBUG
    dd4hep::VolumeID postID  = volMgr.volumeID( stepH.postTouchable() );
    G4cout << " ***************** DEBUG information! " << G4endl ;
    G4cout << " ***************** Original :  preID : " << decoder(preID)  << G4endl ; 
    G4cout << " ***************** Original : postID : " << decoder(postID) << G4endl ; 
#endif
       

    // Layer count from 1, and _cellSize vector count from 0
    // Layer 1 cell size is _cellSize.at(currentLayer-1)
    currentLayer =  decoder(  preID )["layer"];

#ifdef SDHcalBarrel_DEBUG
    G4cout << " SDHcalBarrel::DumpHit -  working on --- " << name() <<"  currentLayer = " << currentLayer << G4endl; 
    G4cout << " ---------------------------------------------------------  "    << G4endl;
#endif

    // boundary check for the layers number and valid cellSize setup
    // layer start from 0??
    if ( currentLayer > (int)(_cellSize.size()) ){
      G4cout << " ***************** BOUNDARY CHECK, WARNING! :" <<  G4endl ;
      G4cout << " ***************** currentLayer: "<<currentLayer << G4endl ;
      G4cout << " ***************** WARNING : Please check the cellSize in sequences.xml and the layers number in ILD_o1_v05.xml" <<  G4endl ;
      G4cout << " ***************** WARNING : Missing cell size in sequences.xml!" <<  G4endl ;
      throw std::runtime_error("Invalid cell size: 0.0");
      }else if (_cellSize.at(currentLayer-1) <= 0.) {
      throw std::runtime_error("Invalid cell size: 0.0");
    }


#ifdef SDHcalBarrel_DEBUG
    G4cout << name()
	   << "====>DumpHit, currentLayer = " << currentLayer << G4endl;
#endif
  
 
    //Read the cellEncoding from  compact.xml readout <id></id> 
    _cellIDEncoding = decoder.fieldDescription();
    dd4hep::BitField64 cellID_decoder(_cellIDEncoding);

    // get the hit position information
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle thepreTouchable = preStepPoint->GetTouchableHandle();
    G4ThreeVector preworldPosition = preStepPoint->GetPosition();
    G4ThreeVector prelocalPosition = thepreTouchable->GetHistory()->
      GetTopTransform().TransformPoint(preworldPosition);

    G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
    G4TouchableHandle thepostTouchable = postStepPoint->GetTouchableHandle();
    G4ThreeVector postworldPosition = postStepPoint->GetPosition();
    G4ThreeVector postlocalPosition = thepostTouchable->GetHistory()->
      GetTopTransform().TransformPoint(postworldPosition);

    //G4ThreeVector localMiddlePoint = ( prelocalPosition + postlocalPosition ) * 0.5;

    G4ThreeVector localMiddlePoint = prelocalPosition;

    double offSize  = 0.0;


    int iID = int(floor((localMiddlePoint.x() + 0.5 * _cellSize.at(currentLayer-1) - offSize)/_cellSize.at(currentLayer-1)));
    int jID = int(floor((localMiddlePoint.y() + 0.5 * _cellSize.at(currentLayer-1) - offSize)/_cellSize.at(currentLayer-1)));
    //int jID = int(floor((localMiddlePoint.z() + 0.5 * _cellSize.at(currentLayer-1) - offSize)/_cellSize.at(currentLayer-1)));

    long volume_cellID = volumeID( aStep ) ;
    
    cellID_decoder.setValue(volume_cellID);
    cellID_decoder["I"] = iID;
    cellID_decoder["J"] = jID;

    long cellID = cellID_decoder.getValue();


    localMiddlePoint.setX( double(iID)*_cellSize.at(currentLayer-1)); // local cartensian center of SD
    localMiddlePoint.setY( double(jID)*_cellSize.at(currentLayer-1)); // local cartensian center of SD
    localMiddlePoint.setZ(0.0);                                       // local thickness  center of SD

    G4ThreeVector worldMiddlePoint = thepreTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint( localMiddlePoint );

    HitTime    = aStep->GetTrack()->GetGlobalTime();
    currentPID = aStep->GetTrack()->GetTrackID();
    currentPDG = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();


    // Check the Birks Law
    if (_applyBirksLawFlag == false) {

      DepositedEnergy = aStep->GetTotalEnergyDeposit();

#ifdef SDHcalBarrel_DEBUG
      G4cout << "  DepositedEnergy: " << DepositedEnergy/keV << " keV" <<G4endl;
#endif

    } else if (_applyBirksLawFlag == true) {

      double attenuatedEnergy = BirkAttenuation(aStep);
      DepositedEnergy = attenuatedEnergy;

#ifdef SDHcalBarrel_DEBUG
      G4cout << "   applyBirksLawFlag: "  << _applyBirksLawFlag                            << G4endl;
      G4cout << "     DepositedEnergy: "  << aStep->GetTotalEnergyDeposit()/keV << " keV"  << G4endl;
      G4cout << " response after Birk: "  << attenuatedEnergy/keV << " keV"                << G4endl;
#endif

    }

#ifdef SDHcalBarrel_DEBUG
    if( abs(preworldPosition[0]-postworldPosition[0])> _cellSize.at(currentLayer-1)
	|| abs(preworldPosition[1]-postworldPosition[1])> _cellSize.at(currentLayer-1) )
    {
    G4cout   <<" \n preworldPosition: "  << preworldPosition
	     <<" \n postworldPosition: " << postworldPosition
	     <<" \n prelocalPosition: "  << prelocalPosition
	     <<" \n postlocalPosition: " << postlocalPosition
	     <<" \n worldMiddlePoint: "  << worldMiddlePoint
	     << G4endl;
  }
#endif

#ifdef SDHcalBarrel_DEBUG 
     if( (sqrt(worldMiddlePoint.x()*worldMiddlePoint.x()+worldMiddlePoint.y()*worldMiddlePoint.y()) < 1800)
	  || sqrt(  ( preworldPosition.x()-worldMiddlePoint.x() ) * ( preworldPosition.x()-worldMiddlePoint.x() )
		    + ( preworldPosition.y()-worldMiddlePoint.y() ) * ( preworldPosition.y()-worldMiddlePoint.y() )
		    + ( preworldPosition.z()-worldMiddlePoint.z() ) * ( preworldPosition.z()-worldMiddlePoint.z() )
		    ) > (_cellSize.at(currentLayer-1)*_cellSize.at(currentLayer-1)) ){

	G4cout << " ***************** DEBUG information! " << G4endl ;
	G4cout << " ***************** WARNING : preID : " << decoder(preID) << G4endl ; 
	G4cout << " ***************** WARNING : postID : " << decoder(postID) << G4endl ; 
	
	G4cout   <<" \n preworldPosition: "  << preworldPosition
		 <<" \n postworldPosition: " << postworldPosition
		 <<" \n prelocalPosition: "  << prelocalPosition
		 <<" \n postlocalPosition: " << postlocalPosition
		 <<" \n worldMiddlePoint: "  << worldMiddlePoint
		 << G4endl;

      }
#endif
 

    IMPL::SimCalorimeterHitImpl* hit = new IMPL::SimCalorimeterHitImpl ;

    hit->setCellID0( cellID & 0xffffffff   ) ;

    hit->setCellID1( ( cellID >> 32 ) & 0xffffffff   ) ;
  
    hit->setEnergy( DepositedEnergy ) ;
  
    float pos[3] ;
    pos[0] = worldMiddlePoint.x() ;
    pos[1] = worldMiddlePoint.y() ;
    pos[2] = worldMiddlePoint.z() ;
    hit->setPosition( pos  ) ;


    collection(m_collectionID)->add(hit);
  

#ifdef SDHcalBarrel_DEBUG
    G4cout << name()
      	   << "====>DumpHit, currentLayer = " << currentLayer <<G4endl;
    G4cout << " worldMiddlePoint(0) = " << worldMiddlePoint.x()//(0)
	   << ", worldMiddlePoint(1) = " << worldMiddlePoint.y()//(1)
	   << ", worldMiddlePoint(2) = " << worldMiddlePoint.z()//(2)
	   << "\n pos[0] = " << pos[0]
	   << ", pos[1] = " << pos[1]
	   << ", pos[2] = " << pos[2]
	   << "\n currentPID = " << currentPID
	   << " currentPDG = " << currentPDG
	   << " DepositedEnergy = " << DepositedEnergy
	   << G4endl;
    G4cout << " ref  preID : " << decoder(preID)  << G4endl ; 
    G4cout << " ref postID : " << decoder(postID) << G4endl ; 
    G4cout << " cellID : " << cellID_decoder  << G4endl ;  
#endif
  }



//==================================================================
//
//  Clear: 
//
//==================================================================
  void SDHcalBarrel::Clear()
  {
#ifdef SDHcalBarrel_DEBUG
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4cout << name()
	   << "====> SDHcalBarrel::Clear " << currentLayer
	   << G4endl;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif
  
    currentLayer = 0;
    currentPID = -1; 
    DepositedEnergy = 0;
    HitTime = 0. ;
    currentPDG = 0;
  
  }

 
 
  void SDHcalBarrel::StartNewHit(G4int aLayerNumber,
			      G4double theDepositedEnergy,
			      G4int thePDG,
			      G4int theTrackPID)
  {      
    currentLayer = aLayerNumber;
    currentPID = theTrackPID;
    currentPDG = thePDG;
    DepositedEnergy = theDepositedEnergy;
    HitTime = 0.;

#ifdef SDHcalBarrel_DEBUG
    G4cout << "StartNewHit " << name()
	   << ", currentLayer = " << aLayerNumber
	   << " currentPID = " << currentPID
	   << " currentPDG = " << currentPDG
	   << G4endl;
#endif


  }
 


  void SDHcalBarrel::UpdateHit(G4Step *aStep)
  {
#ifdef SDHcalBarrel_DEBUG
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4cout << name()
	   << "====> SDHcalBarrel::UpdateHit " << currentLayer
	   << G4endl;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif
    DepositedEnergy+=aStep->GetTotalEnergyDeposit();
    HitTime = aStep->GetTrack()->GetGlobalTime() ; 
  }


//==================================================================
//
// process:
//
//==================================================================
  bool SDHcalBarrel::process(G4Step* aStep,G4TouchableHistory* ) 
  {

    
    //process only if energy>0, except neutrinos
    G4double energyDeposition = aStep->GetTotalEnergyDeposit();

#ifdef SDHcalBarrel_DEBUG
    if (energyDeposition >0 )
      {      
	G4cout << "\n\n #####################################################################  "    << G4endl; 
	G4cout << "             SDHcalBarrel::process -  start...........................  "   << G4endl; 
	G4cout << "             SDHcalBarrel::process -  energyDeposition : "<<energyDeposition   << G4endl;
      }
#endif

    if (energyDeposition <= 0 
	&& aStep->GetTrack()->GetDefinition()->GetParticleType() != "geantino")
      return true;

    dd4hep::Geant4StepHandler stepH( aStep );
    dd4hep::Geant4VolumeManager volMgr = dd4hep::Geant4Mapping::instance().volumeManager();

    dd4hep::VolumeID preID  = volMgr.volumeID( stepH.preTouchable() );

   
    if( preID == 0 ) {
      
      G4cout << " #####################################################################  "    << G4endl; 
      G4cout << "  WARNING:   SDHcalBarrel::process -  invalid prestep volumeID : " << preID  << G4endl; 
      G4cout << " #####################################################################  "    << G4endl;
      
      return false; // this should never happen really ??
    }
    

    emSaturation = new G4EmSaturation();

    DumpHit(aStep); 

    Clear();

    return true;

  }


} // namespace lcgeo

//##############################################################################################

namespace DD4hep{
  namespace Simulation{
    typedef lcgeo::SDHcalBarrel SDHcalBarrel;
  }
}
#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(SDHcalBarrel)

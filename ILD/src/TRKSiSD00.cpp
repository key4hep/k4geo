//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id:$
//====================================================================
#include "TRKSiSD00.h"


#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"

#include "G4SDManager.hh"
#include "G4ios.hh"
#include "globals.hh"
#include <assert.h>

//**********************
#define DEBUGTRKSD
//**********************


#include "IMPL/SimTrackerHitImpl.h"

#include "DDG4/Geant4StepHandler.h"
#include "DDG4/Geant4VolumeManager.h"
#include "DDG4/Geant4Mapping.h"

#include "DDSegmentation/BitField64.h"

namespace DDSim {


  //dodao PK
  //#include "UserTrackInformation.hh"
  //fg #include "TrackSummary.hh"

  //fgTRKSiSD00::TRKSiSD00(G4String SDname, 
  //fg		 G4double pThreshold,
  //fg		 G4double thePrimaryTPCCut) 
  //fg  : VSensitiveDetector(SDname),Threshold(pThreshold),
  //fg    PrimaryTPCCut(thePrimaryTPCCut),
  //fg    HCID(-1), currentCylinder(0), currentPID(-1), 
  //fg    currentPDG(-1),
  //fg    EntryPoint(0.,0.,0.), ExitPoint(0.,0.,0.),
  //fg    EntryMomentum (0.,0.,0.), ExitMomentum (0.,0.,0.),
  //fg    DepositedEnergy(0.), HitTime(0.), StepLength(0.), 
  //fg    CalCollection(0)
  //fg{
  //fg  G4String CollName=SDname+"Collection";
  //fg  collectionName.insert(CollName);
  //fg}

  //fgvoid TRKSiSD00::Initialize(G4HCofThisEvent *)
  //fg{
  //fg  CalCollection = new TRKHitsCollection
  //fg    (SensitiveDetectorName,collectionName[0]);
  //fg  currentCylinder = 0;
  //fg  currentPID = -1;
  //fg  
  //fg  detailedHitsStoring = false;
  //fg  defaultHitsStoring = true;
  //fg  for(G4int i = 0; i < (G4int)Control::detailedHitsStoring.size(); i++){ 
  //fg    if(Control::detailedHitsStoring[i] == SensitiveDetectorName){
  //fg      detailedHitsStoring = true;
  //fg      defaultHitsStoring = false;
  //fg      break;
  //fg    }
  //fg  }//i
  //fg  
  //fg}

  bool TRKSiSD00::process(G4Step* aStep,G4TouchableHistory* history) 
  //G4bool TRKSiSD00::ProcessHits(G4Step *aStep,G4TouchableHistory*)
  {
    
    if (fabs(aStep->GetTrack()->GetDefinition()->GetPDGCharge()) < 0.01) return true;
    
    dd4hep::Geant4StepHandler stepH( aStep );
    dd4hep::Geant4VolumeManager volMgr = dd4hep::Geant4Mapping::instance().volumeManager();
    dd4hep::BitField64& decoder = *  m_readout.idSpec().decoder() ;

    dd4hep::VolumeID preID  = volMgr.volumeID( stepH.preTouchable() );

    dd4hep::VolumeID postID = volMgr.volumeID( stepH.postTouchable() );
    

    if( preID == 0 || preID == 0xffffffe1 ) {

      std::cout << " ##################################################################### \n " 
		<< "  WARNING:   TRKSiSD00::process -  invalid prestep volumeID : " << preID  << " \n " 
		<< " ##################################################################### \n "  ;
      
      return ; // this should never happen really ??????
    }
    
    // debug
    std::cout << " --------------------------------------------------- \n "
	      << "    prestep volID "  << ( preID )  << " \n " 
	      << "    poststep volID " << ( postID ) << " \n " 
	      << " --------------------------------------------------- \n " ;
    
    if( !_detailedHitsStoring){
      
      // It's a first approach for hit's processing for all the tracking
      // detectors (VXD, SIT, FTD and TPC). In this first release, if
      // the total energy deposited by a primary and its secondaries
      // is upper then a given threshould it keeps:
      // 
      // o the layer number (the plan number for the FTD)
      // o the mean step position when crossing the layer
      // o the mean momentum when crossing the layer
      // o the primary PID number
      // o the PDG particle code (it can be the secondary one)
      // o the total energy deposited when crossing the layer
      //
      // The current  
      //
      // (Paulo Sept. 2002)
      //
    
      // // // Just for particles with more than the PrimaryTPCCut
      // // // given at constructor call (to control TPC output
      // // // file length).
      // // if(aStep->GetPreStepPoint()->GetKineticEnergy() 
      // // 	 < PrimaryTPCCut) return true;
        

      //      decoder.setValue( preID ) ;
      G4int PrelayerNumber = decoder( preID )["layer"] ;

	// aStep->GetPreStepPoint()->
	// GetPhysicalVolume()->GetCopyNo();


      // // G4StepPoint* postStep = aStep->GetPostStepPoint();
      // // G4VPhysicalVolume * physVol = postStep->GetPhysicalVolume();
      // // if(physVol == 0) 
      // // 	{	
      // // 	  G4cout << "WARNING: TRKSiD00::ProcessHits: post step point physical volume pointer is null!!!\n"
      // // 		 << "It's a Geant4 bug. TRKSD will skip this hit to avoid aborting the job!" << G4endl;
      // // 	  return true;
      // // 	}
    
      G4int PostlayerNumber =  decoder( postID )["layer"] ;

	//physVol->GetCopyNo();

      /*
	G4int PostlayerNumber = 
	aStep->GetPostStepPoint()->
	GetPhysicalVolume()->GetCopyNo();
      */
      //#define DEBUGTRKSD
#ifdef DEBUGTRKSD
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      G4cout << "ProcessHits " << name()
	     << ", PrelayerNumber = " << PrelayerNumber
	     << " PreStepPoint = " << aStep->GetPreStepPoint()->GetPosition()
	     << ", PostlayerNumber = " << PostlayerNumber 
	     << " PostStepPoint = " << aStep->GetPostStepPoint()->GetPosition()
	//	     << " Control::primaryId = " << Control::primaryId
	     << " TrackPID = " << aStep->GetTrack()->GetTrackID()
	     << " TrackPDG = " << aStep->GetTrack()->GetDefinition()->GetPDGEncoding()
	     << " TrackCharge = " << aStep->GetTrack()->GetDefinition()->GetPDGCharge()
	     << " currentPID = " << currentPID
	     << " currentPDG = " << currentPDG
	     << " DepositedEnergy = " << DepositedEnergy
	     << G4endl;
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif

      // If new primary tracking or another secondary
      // dump and reset the counters.

      if ( aStep->GetTrack()->GetTrackID() != currentPID ) 
	{
	  DumpHit(aStep); 
	  Clear(); // forces to start a new hit
	}

    
      // First case, starting traversal
      if(PrelayerNumber != currentCylinder)
	{
	  // dump and start a new hit.
	  DumpHit(aStep); 
	  StartNewHit(PrelayerNumber,
		      aStep->GetPreStepPoint()->GetPosition(),
		      aStep->GetPreStepPoint()->GetMomentum(),
		      aStep->GetTrack()->GetDefinition()->GetPDGEncoding(),
		      aStep->GetTrack()->GetTrackID());
	  // take the energy and the exit momentum
	  UpdateHit(aStep);
	
	  // Perhaps it's already on the next boundary:
	  if(PrelayerNumber != PostlayerNumber)
	    {
	      // So dump and start a new hit if the layers
	      // share surfaces (the case of TPC)
	      DumpHit(aStep);
	      if( PostlayerNumber!= 0)
		StartNewHit(PostlayerNumber,
			    aStep->GetPostStepPoint()->GetPosition(),
			    aStep->GetPostStepPoint()->GetMomentum(),
			    aStep->GetTrack()->GetDefinition()->GetPDGEncoding(),
			    aStep->GetTrack()->GetTrackID());
	      else
		// If the layers don't share surfaces, reset the counters.
		Clear();
	    }
	
	  // PAY ATTENTION TO THE RETURN HERE!
	  return true;
	}
    
      // Second case, traveling and perhaps on the next boundary
      // add energy and update the exit momentum.
      // We test here if currentCylinder !=0 just to be sure, it should 
      // never happens except if the user plugged the TRKSD in a zero
      // numbered layer.
      if( currentCylinder !=0 ) UpdateHit(aStep);
    
      // Is it on the next boundary?
      if(PrelayerNumber != PostlayerNumber) 
	{
	  // Yes, dump the Hit and perhaps start a new one, if the
	  // layers share surfaces (the TPC case)
	  DumpHit(aStep);
	  if( PostlayerNumber!= 0) 
	    StartNewHit(PostlayerNumber,
			aStep->GetPostStepPoint()->GetPosition(),
			aStep->GetPostStepPoint()->GetMomentum(),
			aStep->GetTrack()->GetDefinition()->GetPDGEncoding(),
			aStep->GetTrack()->GetTrackID());
	  else 
	    // Else clear the counters.
	    Clear();
	}
      else if (aStep->GetTrack()->GetTrackStatus() == fStopAndKill) 
	{
	  DumpHit(aStep);
	  Clear();
	}

      return true;



    } else { //  if(detailedHitsStoring){
    
      if(aStep->GetPreStepPoint()->GetKineticEnergy() < PrimaryTPCCut) return true;
    
      G4int PrelayerNumber = 
	aStep->GetPreStepPoint()->
	GetPhysicalVolume()->GetCopyNo();
    
      G4StepPoint* postStep = aStep->GetPostStepPoint();
      G4VPhysicalVolume * physVol = postStep->GetPhysicalVolume();
      if(physVol == 0) 
	{	
	  G4cout << "WARNING: TRKSiD00::ProcessHits: post step point physical volume pointer is null!!!\n"
		 << "It's a Geant4 bug. TRKSD will skip this hit to avoid aborting the job!" << G4endl;
	  return true;
	}
    
      //fg fixme: which number to use here ?
      G4int PostlayerNumber = physVol->GetCopyNo();
      /*
	G4int PostlayerNumber = 
	aStep->GetPostStepPoint()->
	GetPhysicalVolume()->GetCopyNo();
      */
      //#define DEBUGTRKSD
#ifdef DEBUGTRKSD
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      G4cout << name()
	     << ", PrelayerNumber = " << PrelayerNumber
	     << " PreStepPoint = " << aStep->GetPreStepPoint()->GetPosition()
	     << ", PostlayerNumber = " << PostlayerNumber 
	     << " PostStepPoint = " << aStep->GetPostStepPoint()->GetPosition()
	     << " currentPID = " << currentPID
	     << " currentPDG = " << currentPDG
	     << " DepositedEnergy = " << DepositedEnergy
	     << G4endl;
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif

      Clear();
  
      StartNewHit(PrelayerNumber,
		  aStep->GetPreStepPoint()->GetPosition(),
		  aStep->GetPreStepPoint()->GetMomentum(),
		  aStep->GetTrack()->GetDefinition()->GetPDGEncoding(),
		  aStep->GetTrack()->GetTrackID());
      // take the energy and the exit momentum
      UpdateHit(aStep);
      DumpHit(aStep);
    
      return true;
    
    }

    return true;

  }

  //fg void 
  //fg TRKSiSD00::EndOfEvent(G4HCofThisEvent*HCE)
  //fg {
  //fg   if(HCID<0)
  //fg     { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  //fg   HCE->AddHitsCollection( HCID, CalCollection );
  //fg }

  //fg void 
  //fg TRKSiSD00::LoadEvent(FILE* theSubDetectorEventHitsFileInput)
  //fg {
  //fg   TRKHit* newHit = new TRKHit();
  //fg   while (newHit->Load(theSubDetectorEventHitsFileInput))
  //fg     {
  //fg       CalCollection->insert(newHit);
  //fg       newHit = new TRKHit();
  //fg     }
  //fg   delete newHit;
  //fg }


  void TRKSiSD00::DrawAll()
  {
  } 

  void TRKSiSD00::PrintAll()
  {
  } 

  void TRKSiSD00::DumpHit(G4Step* aStep)
  {
#ifdef DEBUGTRKSD
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4cout << name()
	   << "====>DumpHit, currentCylinder = " << currentCylinder << G4endl;
#endif
  
    // If currentCylinder==0 there is nothing to dump
    if(currentCylinder == 0) return;

    // It keeps the hit only if the total deposited energy
    // is upper then the given Threshold.

    //fg  if( DepositedEnergy < Threshold && Control::TrackingPhysicsListELossOn == true ) return;    is this needed ????

    G4ThreeVector MiddlePoint = 
      ( EntryPoint + ExitPoint ) / 2.;
  
    G4ThreeVector MeanMomentum; 

    if(detailedHitsStoring) MeanMomentum = EntryMomentum;
    else MeanMomentum = ( ExitMomentum + EntryMomentum ) / 2.;

    // Bug: The last G4Step length counted twice
    //StepLength += aStep->GetStepLength();

    //PK: fix for delta electrons: all particles causing hits
    // have to be saved in the LCIO file
    //fg UserTrackInformation* theUserTrackInformation =
    //fg (UserTrackInformation*) (aStep->GetTrack()->GetUserInformation());
  
 
    //fg    if(theUserTrackInformation) {
    //fg      
    //fg      theUserTrackInformation->GetTheTrackSummary()->SetToBeSaved();
    //fg    }
  
    //fg CalCollection->
    //fg   insert(new TRKHit (currentCylinder,
    //fg 		       MiddlePoint (0),
    //fg 		       MiddlePoint (1),
    //fg 		       MiddlePoint (2),
    //fg 		       MeanMomentum (0),
    //fg 		       MeanMomentum (1),
    //fg 		       MeanMomentum (2),
    //fg 		       currentPID,
    //fg 		       currentPDG,
    //fg 		       DepositedEnergy, 
    //fg 		       HitTime,
    //fg 		       StepLength));

    //******************  here we create the lcio SimTrackerHit now instead of the old Mokka TRKHit ****************************

    IMPL::SimTrackerHitImpl* hit = new IMPL::SimTrackerHitImpl ;
  
    //    (h.track->GetTrackID(),
    // h.track->GetDefinition()->GetPDGEncoding(),
    // step->GetTotalEnergyDeposit(),
    // h.track->GetGlobalTime());
    // HitContribution contrib = Hit::extractContribution(step);

    long cellID = volumeID( aStep ) ;

    unsigned cellID0 = ( cellID & 0xffffffff   ) ;

    if( cellID0 == 0xffffffe1 ){
      std::cout << " ***************** WARNING :  invalid cellID0 set to 0xffffffe1 ! " << std::endl ;

      dd4hep::Geant4StepHandler stepH( aStep );
      dd4hep::Geant4VolumeManager volMgr = dd4hep::Geant4Mapping::instance().volumeManager();
      dd4hep::BitField64& decoder = *  m_readout.idSpec().decoder() ;
      dd4hep::VolumeID preID  = volMgr.volumeID( stepH.preTouchable() );
      dd4hep::VolumeID postID  = volMgr.volumeID( stepH.postTouchable() );
      std::cout << " ***************** WARNING : preID : " << decoder(preID) << std::endl ; 
      std::cout << " ***************** WARNING : postID : " << decoder(postID) << std::endl ; 
    }

    hit->setCellID0( cellID & 0xffffffff   ) ;

    hit->setCellID1( ( cellID >> 32 ) & 0xffffffff   ) ;
  
    hit->setEDep( DepositedEnergy ) ;
    // printout(dd4hep::INFO,"DDSimTracker","%s> Add hit with deposit:%f  Pos:%f %f %f - cellID0: 0x%x cellID1: 0x%x",
    // 	     c_name(),step->GetTotalEnergyDeposit(),position.X(),position.Y(),position.Z() , hit->getCellID0() ,hit->getCellID1() );
  
    double pos[3] ;
    pos[0] = MiddlePoint.x() ;
    pos[1] = MiddlePoint.y() ;
    pos[2] = MiddlePoint.z() ;
    hit->setPosition( pos  ) ;
  
    // hit->energyDeposit = contrib.deposit ;
    // hit->position      = position;
    // hit->momentum      = direction;
    // hit->length        = hit_len;
  
    collection(m_collectionID)->add(hit);
  
    //******************************************************************************************************************************

#ifdef DEBUGTRKSD
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4cout << name()
	   << "====>DumpHit, currentCylinder = " << currentCylinder
	   << ", MiddlePoint(0) = " << MiddlePoint(0)
	   << ", MiddlePoint(1) = " << MiddlePoint(1)
	   << ", EntryPoint = " << EntryPoint
	   << ", ExitPoint = " << ExitPoint
	   << " currentPID = " << currentPID
	   << " currentPDG = " << currentPDG
	   << " DepositedEnergy = " << DepositedEnergy
	   << G4endl;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif
  }

  void TRKSiSD00::StartNewHit(G4int aLayerNumber,
			      G4ThreeVector theEntryPoint,
			      G4ThreeVector theEntryMomentum,
			      G4int thePDG,
			      G4int theTrackPID)
  {      
    currentCylinder = aLayerNumber;
    currentPID = theTrackPID;
    currentPDG = thePDG;
    EntryPoint = ExitPoint = theEntryPoint;
    EntryMomentum = ExitMomentum = theEntryMomentum;
    DepositedEnergy = 0.;
    HitTime = 0.;
    StepLength = 0.;

#ifdef DEBUGTRKSD
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4cout << "StartNewHit " << name()
	   << ", currentCylinder = " << aLayerNumber
	   << " theEntryPoint = " << theEntryPoint
	   << " currentPID = " << currentPID
	   << " currentPDG = " << currentPDG
	   << G4endl;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif


  }

  void TRKSiSD00::Clear()
  {
#ifdef DEBUGTRKSD
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4cout << name()
	   << "====> TRKSiSD00::Clear " << currentCylinder
	   << G4endl;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif
  
    currentCylinder = 0;
    currentPID = -1; 
    DepositedEnergy = 0;
    HitTime = 0. ;
    StepLength  = 0. ;
    currentPDG = 0;
  
  }

  void TRKSiSD00::UpdateHit(G4Step *aStep)
  {
#ifdef DEBUGTRKSD
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4cout << name()
	   << "====> TRKSiSD00::UpdateHit " << currentCylinder
	   << G4endl;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif
    DepositedEnergy+=aStep->GetTotalEnergyDeposit();
    HitTime = aStep->GetTrack()->GetGlobalTime() ; 
    ExitPoint=aStep->GetPostStepPoint()->GetPosition();
    ExitMomentum=aStep->GetPostStepPoint()->GetMomentum();
    StepLength += aStep->GetStepLength();
  }


} // namespace DDSim

//##############################################################################################

namespace DD4hep{
  namespace Simulation{
    typedef DDSim::TRKSiSD00 TRKSiSD00;
  }
}
#include "DDG4/Factories.h"

DECLARE_GEANT4SENSITIVE(TRKSiSD00)

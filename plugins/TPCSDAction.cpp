#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4Mapping.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
  
  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim   {
    
    /**
     *  Geant4SensitiveAction<TPCSdData> sensitive detector for the special case of 
     *  of a TPC, where every pad row is devided into two halfs in order to get
     *  the position from the crossing of the middle of the pad row from
     *  geant4 via volume boundary. Ported of Mokka/TPCSD04.cc
     * 
     *  \author  F.Gaede ( ported from Mokka/TPCSD04.cc )
     *  \version 1.0
     *  \ingroup DD4HEP_SIMULATION
     */
    struct TPCSDData {

      // helper struct with Mokka controls ...
      struct {
	double TPCLowPtCut {};
	bool   TPCLowPtStepLimit {};
        double TPCLowPtMaxHitSeparation {};
      } Control {};

      typedef Geant4HitCollection HitCollection;
      Geant4Sensitive*  sensitive{};
      const BitFieldElement* layerField {};

      G4double fThresholdEnergyDeposit{};
      Geant4HitCollection *fHitCollection{};
      Geant4HitCollection *fSpaceHitCollection{};
      Geant4HitCollection *fLowPtHitCollection{};
      G4int fHCID{};
      G4int fSpaceHitCollectionID{};
      G4int fLowPtHitCollectionID{};
      
      G4ThreeVector CrossingOfPadRingCentre{};
      G4ThreeVector MomentumAtPadRingCentre{};
      G4double dEInPadRow{};
      G4double globalTimeAtPadRingCentre{};
      G4double pathLengthInPadRow{};
      G4double CumulativePathLength{};
      G4double CumulativeEnergyDeposit{};
      G4ThreeVector CumulativeMeanPosition{};
      G4ThreeVector CumulativeMeanMomentum{};
      G4int CumulativeNumSteps{};
      
      G4ThreeVector PreviousPostStepPosition{}; //< the end point of the previous step
      G4int CurrentPDGEncoding{}; //< the PDG encoding of the particle causing the cumulative energy deposit
      G4int CurrentTrackID{}; //< the TrackID of the particle causing the cumulative energy deposit
      G4double CurrentGlobalTime{}; ///< the global time of the track causing the cumulative energy deposit
      G4int CurrentCopyNumber{}; ///< copy number of the preStepPoint's TouchableHandle for the cumulative energy deposit
      

      TPCSDData() : 
	fThresholdEnergyDeposit(0),
	// fHitCollection(0),
	// fSpaceHitCollection(0),
	// fLowPtHitCollection(0),
	fHCID(-1),
	fSpaceHitCollectionID(-1),
	fLowPtHitCollectionID(-1),
	CurrentTrackID(-1) {


	Control.TPCLowPtCut = CLHEP::MeV ;
	Control.TPCLowPtStepLimit = false ;
	Control.TPCLowPtMaxHitSeparation = 5. * CLHEP::mm ;

      }


     /// Clear collected information and restart for new hit
      void clear()  {
	// nothing to clear
      }


      /// return the layer number of the volume (either pre or post-position )
      int getCopyNumber(G4Step* s, bool usePostPos ){

	int cellID = this->volID( s , usePostPos) ;

	return this->layerField->value( cellID ) ;
      }


      /// Returns the volumeID of sensitive volume corresponding to the step (either pre or post-position )
      long long int volID( G4Step* s, bool usePostPos=false ) {

	Geant4StepHandler h(s);

	Geant4VolumeManager volMgr = Geant4Mapping::instance().volumeManager();

	VolumeID volID = ( usePostPos ?  volMgr.volumeID(h.postTouchable()) : volMgr.volumeID(h.preTouchable()) );

	return volID;
      }


      void dumpStep( Geant4StepHandler h, G4Step* s){

	std::cout << " ----- step in detector " << h.sdName( s->GetPreStepPoint() ) 
		  << " prePos  " << h.prePos() 
		  << " postPos " << h.postPos()
		  << " preStatus  " << h.preStepStatus()   
		  << " postStatus  " << h.postStepStatus()   
		  << " preVolume " << h.volName( s->GetPreStepPoint() ) 
		  << " postVolume " << h.volName( s->GetPostStepPoint() ) 
		  << " CurrentCopyNumbers " <<  getCopyNumber( s, false )
		  << " -  " <<  getCopyNumber( s, true )
		  << std::endl
		  << "     momentum : "  << std::scientific
		  <<  s->GetPreStepPoint()->GetMomentum()[0] << ", " <<  s->GetPreStepPoint()->GetMomentum()[1]<< ", " <<  s->GetPreStepPoint()->GetMomentum()[2]
		  << " / "
		  << s->GetPostStepPoint()->GetMomentum()[0] << ", " <<  s->GetPostStepPoint()->GetMomentum()[1]<< ", " <<  s->GetPostStepPoint()->GetMomentum()[2]
		  << ", PDG: " << s->GetTrack()->GetDefinition()->GetPDGEncoding()
		  << std::endl ;
      }




      /// Method for generating hit(s) using the information of G4Step object.
      G4bool process(G4Step* step, G4TouchableHistory* ) {


	fHitCollection = sensitive->collection(0) ;
	fSpaceHitCollection = sensitive->collection(1) ;
	fLowPtHitCollection = sensitive->collection(2) ;

	Geant4StepHandler h(step);
	//	dumpStep( h , step ) ;

	// FIXME: 
	// (i) in the following algorithm if a particle "appears" within a pad-ring half and 
	// leaves without passing through the middle of the pad-ring it will not create a hit in 
	// this ring
	// (ii) a particle that crosses the boundry between two pad-ring halves will have the hit 
	// placed on this surface at the last crossing point, and will be assinged the total energy 
	// deposited in the whole pad-ring. This is a possible source of bias for the hit
  
  
	if (fabs(step->GetTrack()->GetDefinition()->GetPDGCharge()) < 0.01) return true;
  
	const G4ThreeVector PrePosition = step->GetPreStepPoint()->GetPosition();
	const G4ThreeVector PostPosition = step->GetPostStepPoint()->GetPosition();
	const G4ThreeVector thisMomentum = step->GetPostStepPoint()->GetMomentum();
  
	float ptSQRD = thisMomentum[0]*thisMomentum[0]+thisMomentum[1]*thisMomentum[1];


	//=========================================================================================================

	if( ptSQRD >= (Control.TPCLowPtCut*Control.TPCLowPtCut) ){

	  //=========================================================================================================
	  // Step finishes at a geometric boundry

	  if(step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {

	    // step within the same pair of upper and lower pad ring halves
	    if( getCopyNumber( step, false ) == getCopyNumber( step, true ) ){

	      //this step must have ended on the boundry between these two pad ring halfs 
	      //record the tracks coordinates at this position 
	      //and return
        
	      CrossingOfPadRingCentre = PostPosition;
	      MomentumAtPadRingCentre = thisMomentum;
	      dEInPadRow += step->GetTotalEnergyDeposit();
	      globalTimeAtPadRingCentre = step->GetTrack()->GetGlobalTime();
	      pathLengthInPadRow += step->GetStepLength();
        
	      //	    G4cout << "step must have ended on the boundry between these two pad ring halfs" << G4endl;
	      //	    G4cout << "CrossingOfPadRingCentre = "   
	      //		   << sqrt( CrossingOfPadRingCentre[0]*CrossingOfPadRingCentre[0]
	      //			    +CrossingOfPadRingCentre[1]*CrossingOfPadRingCentre[1] )
	      //		   << G4endl;
        
	      return true;

	    }
	    else if(!(CrossingOfPadRingCentre[0]==0.0 && CrossingOfPadRingCentre[1]==0.0 && CrossingOfPadRingCentre[2]==0.0)) {
	      // the above IF statment is to catch the case where the particle "appears" in this pad-row half volume and 
	      // leaves with out crossing the pad-ring centre, as mentioned above
        
	      // it is leaving the pad ring couplet
	      // write out a hit

	      //	    G4cout << "step must be leaving the pad ring couplet" << G4endl;
	      //	    G4cout << "write out hit at " 
	      //		   << sqrt( CrossingOfPadRingCentre[0]*CrossingOfPadRingCentre[0]
	      //			    +CrossingOfPadRingCentre[1]*CrossingOfPadRingCentre[1] )
	      //		   << " " << "dEdx = " << step->GetTotalEnergyDeposit()+dEInPadRow 
	      //		   << " " << "step length = " << step->GetStepLength()+pathLengthInPadRow  
	      //		   << G4endl;
        
	      G4double dE = step->GetTotalEnergyDeposit()+dEInPadRow;
        
	      if ( dE > fThresholdEnergyDeposit ) {
          

//xx		fHitCollection->
//xx		  insert(new TRKHit(touchPre->GetCopyNumber(), 
//xx				    CrossingOfPadRingCentre[0],CrossingOfPadRingCentre[1],CrossingOfPadRingCentre[2],
//xx				    MomentumAtPadRingCentre[0],
//xx				    MomentumAtPadRingCentre[1],
//xx				    MomentumAtPadRingCentre[2],
//xx				    info->GetTheTrackSummary()->GetTrackID(),
//xx				    step->GetTrack()->GetDefinition()->GetPDGEncoding(),
//xx				    dE, 
//xx				    globalTimeAtPadRingCentre,
//xx				    step->GetStepLength()+pathLengthInPadRow));
//xx
		Geant4Tracker::Hit* hit = new Geant4Tracker::Hit(step->GetTrack()->GetTrackID(),
								 step->GetTrack()->GetDefinition()->GetPDGEncoding(),
								 dE, globalTimeAtPadRingCentre);
		hit->position = CrossingOfPadRingCentre ;
		hit->momentum = MomentumAtPadRingCentre;
		hit->length   = step->GetStepLength();
		hit->cellID   = sensitive->cellID( step ) ;

		fHitCollection->add(hit);

		sensitive->printM2("+++ TrackID:%6d [%s] CREATE TPC hit at pad row crossing :"
				   " %e MeV  Pos:%8.2f %8.2f %8.2f",
				   step->GetTrack()->GetTrackID(),sensitive->c_name(), dE,
				   hit->position.X()/CLHEP::mm,hit->position.Y()/CLHEP::mm,hit->position.Z()/CLHEP::mm);
          
	      }
        
	      // zero cumulative variables 
	      dEInPadRow = 0.0;
	      globalTimeAtPadRingCentre=0.0;
	      pathLengthInPadRow=0.0;
	      CrossingOfPadRingCentre[0]=0.0;
	      CrossingOfPadRingCentre[1]=0.0;
	      CrossingOfPadRingCentre[2]=0.0;
	      MomentumAtPadRingCentre[0]=0.0;
	      MomentumAtPadRingCentre[1]=0.0;
	      MomentumAtPadRingCentre[2]=0.0;
	      return true;
	    }



	  }
    
	  //=========================================================================================================
	  // case for which the step remains within geometric volume

	  //FIXME: need and another IF case to catch particles which Stop within the padring

	  else {    // else if(step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) {

	    // the step is not limited by the step length
	    if( step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() != "StepLimiter"){

	      // if(particle not stoped){
	      // add the dEdx and return
	      //	    G4cout << "Step ended by Physics Process: Add dEdx and carry on" << G4endl;
	      dEInPadRow += step->GetTotalEnergyDeposit();
	      pathLengthInPadRow += step->GetStepLength();
	      return true;
	      //}
	      //else{
	      //  write out the hit and clear counters
	      //}


	    } else {  // GetProcessName() == "StepLimiter"
        
	      //	    G4cout << "step must have been stopped by the step limiter" << G4endl;
	      //	    G4cout << "write out hit at " 
	      //		   << sqrt( position_x*position_x
	      //			    +position_y*position_y )
	      //		   << " " << "dEdx = " << step->GetTotalEnergyDeposit() 
	      //		   << " " << "step length = " << step->GetStepLength()  
	      //		   << G4endl;
        
	      // write out step limited hit 
	      // these are just space point hits so do not save the dE, which is set to ZERO
	      //	    if ( step->GetTotalEnergyDeposit() > fThresholdEnergyDeposit ) 
        
//xx	      fSpaceHitCollection->
//xx		insert(new TRKHit(touchPre->GetCopyNumber(), 
//xx				  position_x,position_y,position_z,
//xx				  momentum_x,momentum_y,momentum_z,
//xx				  info->GetTheTrackSummary()->GetTrackID(), 
//xx				  step->GetTrack()->GetDefinition()->GetPDGEncoding(),
//xx				  0.0, // dE set to ZERO 
//xx				  step->GetTrack()->GetGlobalTime(),
//xx				  step->GetStepLength()));

	      Geant4Tracker::Hit* hit = new Geant4Tracker::Hit(step->GetTrack()->GetTrackID(),
							       step->GetTrack()->GetDefinition()->GetPDGEncoding(),
							       0.0, globalTimeAtPadRingCentre);  // dE set to ZERO

	      hit->position = 0.5*( PrePosition + PostPosition );
	      hit->momentum = thisMomentum ;
	      hit->length   = step->GetStepLength();
	      hit->cellID   = sensitive->cellID( step ) ;

	      fSpaceHitCollection->add(hit);



              // add dE and pathlegth and return
	      dEInPadRow += step->GetTotalEnergyDeposit();
	      pathLengthInPadRow += step->GetStepLength();
	      return true;
	    }
	  }



	}
        //=========================================================================================================
	//   ptSQRD < (Control.TPCLowPtCut*Control.TPCLowPtCut)

	else if (Control.TPCLowPtStepLimit) { // low pt tracks will be treated differently as their step length is limited by the special low pt steplimiter
    
	  if ( ( PreviousPostStepPosition - step->GetPreStepPoint()->GetPosition() ).mag() > 1.0e-6 * CLHEP::mm ) {
      
	    // This step does not continue the previous path. Deposit the energy and begin a new Pt hit.
      
	    if (CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
	      //dumpStep( h , step ) ;
	      DepositLowPtHit();
	    }
      
	    else {
	      // reset the cumulative variables if the hit has not been deposited.
	      // The previous track has ended and the cumulated energy left at the end 
	      // was not enough to ionize
	      //G4cout << "reset due to new track , discarding " << CumulativeEnergyDeposit / eV << " eV" << std::endl;
	      ResetCumulativeVariables();
	    }

	  }

	  else {
	    //G4cout << "continuing track" << endl;
	  }
    
	  CumulateLowPtStep(step);  

    
	  // check whether to deposit the hit
	  if( ( CumulativePathLength > Control.TPCLowPtMaxHitSeparation )  ) {
      
	    // hit is deposited because the step limit is reached and there is enough energy
	    // to ionize
      
	    if ( CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
	      //dumpStep( h , step ) ;
	      DepositLowPtHit();
	    }
	    //else {
	    //G4cout << "not deposited, energy is " << CumulativeEnergyDeposit/eV << " eV" << std::endl;
	    //}
	  }
	  else { // even if the step lenth has not been reached the hit might
	    // be deposited because the particle track ends

	    if ( step->GetPostStepPoint()->GetKineticEnergy() == 0 ) {


              // only deposit the hit if the energy is high enough
	      if (CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
          
	        //dumpStep( h , step ) ;
		DepositLowPtHit();
	      }
        
	      else { // energy is not enoug to ionize.
		// However, the track has ended and the energy is discarded and not added to the next step
		//G4cout << "reset due to end of track, discarding " << CumulativeEnergyDeposit/eV << " eV" << std::endl;
		ResetCumulativeVariables();
	      }
	    }
	  }
    
	  PreviousPostStepPosition = step->GetPostStepPoint()->GetPosition();    
    
	  return true;
    
	}

	return true;


      }


      /// Post-event action callback
      void endEvent(const G4Event* /* event */)   {
	// // We need to add the possibly last added hit to the collection here.
	// // otherwise the last hit would be assigned to the next event and the
	// // MC truth would be screwed.
	// //
	// // Alternatively the 'update' method would become rather CPU consuming,
	// // beacuse the extract action would have to be recalculated over and over.
	// if ( current > 0 )   {
	//   Geant4HitCollection* coll = sensitive->collection(0);
	//   extractHit(coll);
      }
  
      void ResetCumulativeVariables()
      {
	CumulativeMeanPosition.set(0.0,0.0,0.0);
	CumulativeMeanMomentum.set(0.0,0.0,0.0);
	CumulativeNumSteps = 0;
	CumulativeEnergyDeposit = 0;
	CumulativePathLength = 0;
      }

      void DepositLowPtHit()
      {

//xx	fLowPtHitCollection->
//xx	  insert(new TRKHit(CurrentCopyNumber, 
//xx			    CumulativeMeanPosition[0], CumulativeMeanPosition[1], CumulativeMeanPosition[2],
//xx			    CumulativeMeanMomentum[0], CumulativeMeanMomentum[1], CumulativeMeanMomentum[2],
//xx			    CurrentTrackID, 
//xx			    CurrentPDGEncoding,
//xx			    CumulativeEnergyDeposit, 
//xx			    CurrentGlobalTime,
//xx			    CumulativePathLength));
  
	Geant4Tracker::Hit* hit = new Geant4Tracker::Hit(CurrentTrackID, CurrentPDGEncoding,
							 CumulativeEnergyDeposit, globalTimeAtPadRingCentre);

	hit->position = CumulativeMeanPosition ;
	hit->momentum = CumulativeMeanMomentum ;
	hit->length   = CumulativePathLength ;
	hit->cellID   = CurrentCopyNumber ;

	fLowPtHitCollection->add(hit);

	// reset the cumulative variables after positioning the hit
	ResetCumulativeVariables();
      }
      
      void CumulateLowPtStep(G4Step *step)
      {
	
	const G4ThreeVector meanPosition = (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition()) / 2;
	const G4ThreeVector meanMomentum = (step->GetPreStepPoint()->GetMomentum() + step->GetPostStepPoint()->GetMomentum()) / 2;
	
	++CumulativeNumSteps;    
	CumulativeMeanPosition = ( (CumulativeMeanPosition*(CumulativeNumSteps-1)) + meanPosition ) / CumulativeNumSteps;
	CumulativeMeanMomentum = ( (CumulativeMeanMomentum*(CumulativeNumSteps-1)) + meanMomentum ) / CumulativeNumSteps;
	CumulativeEnergyDeposit += step->GetTotalEnergyDeposit();
	CumulativePathLength += step->GetStepLength();
	CurrentPDGEncoding = step->GetTrack()->GetDefinition()->GetPDGEncoding();
	CurrentTrackID = step->GetTrack()->GetTrackID();
	CurrentGlobalTime = step->GetTrack()->GetGlobalTime();
	CurrentCopyNumber = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
	
      }
      
    };


    /// Initialization overload for specialization
    template <> void Geant4SensitiveAction<TPCSDData>::initialize() {
      eventAction().callAtEnd(&m_userData,&TPCSDData::endEvent);

      declareProperty("TPCLowPtCut",              m_userData.Control.TPCLowPtCut ); 
      declareProperty("TPCLowPtStepLimit",        m_userData.Control.TPCLowPtStepLimit );
      declareProperty("TPCLowPtMaxHitSeparation", m_userData.Control.TPCLowPtMaxHitSeparation );

      m_userData.fThresholdEnergyDeposit = m_sensitive.energyCutoff();
      m_userData.sensitive = this;

      IDDescriptor dsc = m_sensitive.idSpec() ;
      m_userData.layerField = dsc.field( "layer" ) ;

    }

    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<TPCSDData>::defineCollections() {
      m_collectionID = defineCollection<Geant4Tracker::Hit>(m_sensitive.readout().name());
      m_userData.fSpaceHitCollectionID = defineCollection<Geant4Tracker::Hit>("TPCSpacePointCollection");
      m_userData.fLowPtHitCollectionID = defineCollection<Geant4Tracker::Hit>("TPCLowPtCollection");
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> void Geant4SensitiveAction<TPCSDData>::clear(G4HCofThisEvent*) {
      m_userData.clear();
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> G4bool
    Geant4SensitiveAction<TPCSDData>::process(G4Step* step, G4TouchableHistory* history) {
      return m_userData.process(step, history);
    }

    typedef Geant4SensitiveAction<TPCSDData>  TPCSDAction;

  }
}


#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE( TPCSDAction )

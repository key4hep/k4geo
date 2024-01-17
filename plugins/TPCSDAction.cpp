#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4Mapping.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

//DD4hep
#include "DDRec/DetectorData.h"
#include "IMPL/SimTrackerHitImpl.h"


//helper class
#include "InterpolatedHitGenerator.h"

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
struct trackIDExtension : lcrtrel::LCIntExtension<trackIDExtension> {};

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
	G4ThreeVector CumulativeMeanPosition{}, currentMeanPosition{}, firstStepPrePosition{},firstStepPreMomentum{};
	G4ThreeVector CumulativeMeanMomentum{}, currentMeanMomentum{};
	G4int CumulativeNumSteps{};
	G4int firstStepTrackID;

	bool pixelTPC{false}, fastSimulation{false};
	dd4hep::rec::FixedPadSizeTPCData* tpcData{nullptr};

	G4ThreeVector PreviousPostStepPosition{}; //< the end point of the previous step
	G4Step* previousStep{}; //pointer to previous step
	G4int CurrentPDGEncoding{}; //< the PDG encoding of the particle causing the cumulative energy deposit
	G4int CurrentTrackID{}; //< the TrackID of the particle causing the cumulative energy deposit
	G4double CurrentGlobalTime{}; ///< the global time of the track causing the cumulative energy deposit
	G4int CurrentCopyNumber{}; ///< copy number of the preStepPoint's TouchableHandle for the cumulative energy deposit
	G4int eventNumber{0}; //event number


	TPCSDData() :
		fThresholdEnergyDeposit(0),
		// fHitCollection(0),
		// fSpaceHitCollection(0),
		// fLowPtHitCollection(0),
		fHCID(-1),
		pixelTPC(false),
		fSpaceHitCollectionID(-1),
		fLowPtHitCollectionID(-1),
		CurrentTrackID(-1),
		previousStep(new G4Step())
	{
		ResetCumulativeVariables();
		Control.TPCLowPtCut = CLHEP::MeV ;
		Control.TPCLowPtStepLimit = false ;
		Control.TPCLowPtMaxHitSeparation = 5. * CLHEP::mm ;

	}

//	~TPCSDData() {
//		delete previousStep;
//	}

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

		G4cout << " ----- step in detector " << h.sdName( s->GetPreStepPoint() )
		  		<< " prePos  " << h.prePos()
					<< " postPos " << h.postPos()
					<< " preStatus  " << h.preStepStatus()
					<< " postStatus  " << h.postStepStatus()
					<< " preVolume " << h.volName( s->GetPreStepPoint() )
					<< " postVolume " << h.volName( s->GetPostStepPoint() )
					<< " stepLength " << s->GetStepLength()
					<< " EnergyDeposit " << s->GetTotalEnergyDeposit()/CLHEP::eV << " eV"
					<< " CurrentCopyNumbers " <<  getCopyNumber( s, false )
					<< " -  " <<  getCopyNumber( s, true )
					<< "\n"
					<< "     momentum : "  << std::scientific
					<<  s->GetPreStepPoint()->GetMomentum()[0] << ", " <<  s->GetPreStepPoint()->GetMomentum()[1]<< ", " <<  s->GetPreStepPoint()->GetMomentum()[2]
																																																																											 << " / "
																																																																											 << s->GetPostStepPoint()->GetMomentum()[0] << ", " <<  s->GetPostStepPoint()->GetMomentum()[1]<< ", " <<  s->GetPostStepPoint()->GetMomentum()[2]
																																																																																																																																																			<< ", PDG: " << s->GetTrack()->GetDefinition()->GetPDGEncoding()
																																																																																																																																																			<< "\n";
		G4cout<<"current cumulative energy is "<<CumulativeEnergyDeposit/CLHEP::eV<<" eV"<<std::endl;

	}

	/// Method for generating hit(s) using the information of G4Step object.
	G4bool processPad(G4Step* step, G4TouchableHistory* ) {


		fHitCollection = sensitive->collection(0) ;
		fSpaceHitCollection = sensitive->collection(1) ;
		fLowPtHitCollection = sensitive->collection(2) ;

		Geant4StepHandler h(step);
		//	dumpStep( h , step ) ;

		// FIXME:
		// (i) in the following algorithm if a particle "appears" within a pad-ring half and
		// leaves without passing through the middle of the pad-ring it will not create a hit in
		// this ring
		// (ii) a particle that crosses the boundary between two pad-ring halves will have the hit
		// placed on this surface at the last crossing point, and will be assigned the total energy
		// deposited in the whole pad-ring. This is a possible source of bias for the hit

		//check charge
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
						//G4cout << "reset due to end of track, discarding " << CumulativeEnergyDeposit/CLHEP::eV << " eV" << std::endl;
						ResetCumulativeVariables();
					}
				}
			}

			PreviousPostStepPosition = step->GetPostStepPoint()->GetPosition();

			return true;

		}

		return true;


	}


	/// Method for generating hit(s) using the information of G4Step object for a pad readout
	G4bool processPixel(G4Step* step, G4TouchableHistory* ) {

		//do only primary track:
		//    if(step->GetTrack()->GetTrackID()!=1) return true;
		//do all except primary track:
		//    if(step->GetTrack()->GetTrackID()==1) return true;



		//get positions from step
		//    G4cout<<"get position and momentum from step"<<std::endl;
		const G4ThreeVector PrePosition = step->GetPreStepPoint()->GetPosition();
		const G4ThreeVector PostPosition = step->GetPostStepPoint()->GetPosition();

		// This step does not continue the previous path. Deposit the energy at the previous step.
		if ( ( PreviousPostStepPosition - step->GetPreStepPoint()->GetPosition() ).mag() > 1.0e-6 * CLHEP::mm //new hit is far from previous
				|| ( step->GetTrack()->GetTrackID()!=CurrentTrackID ) ) { //or track id is different
			//		G4cout<<"found new track "<<step->GetTrack()->GetTrackID()<<std::endl;
			if (CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
				DepositMultiplePtHits(previousStep);
			}
			// reset the cumulative variables if the hit has not been deposited. The previous track has ended and the accumulated energy left at the end was not enough to ionize
			ResetCumulativeVariables(false); //reset variables, except energy
		}

		//fill histograms, do at this point because last step must be completed, but new step cannot be started yet!
		//	G4cout<<"Fill histograms"<<std::endl;
		//	_deltahist.Fill(step, eventNumber);

		//    dumpStep( Geant4StepHandler(step) , step ) ;

		//increase all counters
		//	G4cout<<"Cumulate step"<<std::endl;
		CumulatePtStep(step);

		//crossed boundary
		if(step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
			if ( CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
				DepositMultiplePtHits(step);
				ResetCumulativeVariables(false); //reset variables, except energy
			}
			//In contrast to pad readout, do not reset after crossing a boundary here.
		}

		// Maximum hit separation reached
		if( ( CumulativePathLength > Control.TPCLowPtMaxHitSeparation )   ) {
			// hit is deposited because the step limit is reached and there is enough energy to ionize
			if ( CumulativeEnergyDeposit > fThresholdEnergyDeposit) {
				DepositMultiplePtHits(step);
			}
			ResetCumulativeVariables(false); //reset variables, except energy
		}


		//save last postStepPosition
		PreviousPostStepPosition = step->GetPostStepPoint()->GetPosition();
		*previousStep = *step; //save previous step

		return true;
	}



	/// Method for generating hit(s) using the information of G4Step object.
	G4bool process(G4Step* step, G4TouchableHistory* th) {
		//    	  G4cout<<"processing step in TPCSDAction!"<<std::endl;
		//dumpStep(Geant4StepHandler(step),step);

		//set collections
		fHitCollection = sensitive->collection(0) ;
		fSpaceHitCollection = sensitive->collection(1) ;
		fLowPtHitCollection = sensitive->collection(2) ;


		//check charge
		if (fabs(step->GetTrack()->GetDefinition()->GetPDGCharge()) < 0.01) return true;

		//call pixel or pad readout
		if(pixelTPC) return processPixel(step,th);
		else return processPad(step,th);
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

	void ResetCumulativeVariables(bool resetEnergy=true)
	{
		CumulativeMeanPosition.set(0.0,0.0,0.0);
		CumulativeMeanMomentum.set(0.0,0.0,0.0);
		currentMeanPosition.set(0.0,0.0,0.0);
		currentMeanMomentum.set(0.0,0.0,0.0);
		firstStepPrePosition.set(0.0,0.0,0.0);
		firstStepPreMomentum.set(0.0,0.0,0.0);
		firstStepTrackID=0;
		CumulativeNumSteps = 0;
		if(resetEnergy) CumulativeEnergyDeposit = 0;
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

	//deposit hit for pixel
	void DepositPtHit(G4Step* step) {
		DepositPtHit(step, CumulativeEnergyDeposit);
	}

	int getCellIDFromFormula() {
		const double toMM=10;
		int system=(currentMeanPosition.z()<0) ? 100 : 36;//FIXME: properly find system id of tpc
		int nrow=(currentMeanPosition.perp()-(tpcData->rMinReadout)*toMM)/(tpcData->padHeight*toMM);
		nrow=nrow<<7;
		return system+nrow;
	}

	void WriteHit(G4ThreeVector position, G4ThreeVector momentum, long long cid, double EDep, double pathLength) {
		IMPL::SimTrackerHitImpl* hit = new IMPL::SimTrackerHitImpl();
		//use SimTrack
		hit->setEDep(EDep);
		double cmposition[3] = {position.getX(), position.getY(), position.getZ()};
		hit->setPosition(cmposition);
		float cmmomentum[3] = {float(momentum.getX()/CLHEP::GeV), float(momentum.getY()/CLHEP::GeV), float(momentum.getZ()/CLHEP::GeV)}; //why is position double and momentum float?!
		hit->setMomentum(cmmomentum);
		hit->setPathLength(pathLength);
		hit->setCellID0(int(cid));
		hit->setCellID1(int(cid>>32));
		hit->ext<trackIDExtension>()=CurrentTrackID;
		//missing:CurrentTrackID,CurrentPDGEncoding,
		fHitCollection->add(hit);

		//histogram
		//    		_deltahist.FillHit(position, momentum);
	}

	int _nprint=0;
	void displayHit(G4Step* step) {
		//display this hit
		long long cid=sensitive->cellID( step );
		//long long cid=getCellIDFromFormula();
		float cmmomentum[3] = {float(currentMeanMomentum.getX()), float(currentMeanMomentum.getY()), float(currentMeanMomentum.getZ())}; //why is position double and momentum float?!
		double cmposition[3] = {currentMeanPosition.getX(), currentMeanPosition.getY(), currentMeanPosition.getZ()};
		const int PRINTEVERYN=100000;
		if(!(_nprint++%PRINTEVERYN))
			printf(
					"+++ TrackID:%6d [%s] CREATE TPC hit at pixel row crossing :"
					" %e MeV  Pos:%8.2f %8.2f %8.2f Mom:%8.1f %8.1f %8.1f ID:%i-%i (%i)\n",
					step->GetTrack()->GetTrackID(),sensitive->c_name(), CumulativeEnergyDeposit,
					cmposition[0]/CLHEP::mm,cmposition[1]/CLHEP::mm,cmposition[2]/CLHEP::mm,
					cmmomentum[0]/CLHEP::MeV,cmmomentum[1]/CLHEP::MeV,cmmomentum[2]/CLHEP::MeV,
					int(cid), int(cid>>32), getCellIDFromFormula());
		//		  	  std::cout<<"tpcData->padheight="<<tpcData->padHeight<<" tpcData->rMinReadout="<<tpcData->rMinReadout<<" r="<<currentMeanPosition.perp()<<std::endl;
	}

	void DepositPtHit(G4Step *step, double EDep) {
		long long cid=sensitive->cellID( step );
		//long long cid=getCellIDFromFormula();
		//    	  G4DynamicParticle mcparticle=step->GetTrack()->GetDynamicParticle();
		WriteHit(CumulativeMeanPosition, CumulativeMeanMomentum, cid, EDep, CumulativePathLength);
		displayHit(step);
	}



	int DepositFastPtHits(G4Step* step) {
		G4StepPoint* postStepPoint = step->GetPostStepPoint();
		std::function<double(int)> hitFunction;//use flat hit function for secundary electrons and triangular for the higher energy primary electrons
		if ( step->GetTrack()->GetParentID() ) hitFunction=FlatMultipleHitFunction(0);
		else hitFunction=TriangularHitFunction(0.1) ;
		InterpolatedHitGenerator hitGenerator(
				hitFunction,
				AngleInterpolation(firstStepPrePosition, postStepPoint->GetPosition(), firstStepPreMomentum, postStepPoint->GetMomentum())
				//			  LinearVectorInterpolation(firstStepPrePosition, postStepPoint->GetPosition())
		);

		int nHits=CumulativeEnergyDeposit/fThresholdEnergyDeposit;
		std::vector<G4ThreeVector> hitpos=hitGenerator.generateHits(nHits),
				hitmom=hitGenerator.generateFromFunction( LinearVectorInterpolation( firstStepPreMomentum, postStepPoint->GetMomentum() ) );

		for(int i=0; i<nHits; i++) {
			long long cid=getCellIDFromFormula();
			double EDep=fThresholdEnergyDeposit; //last hits deposits slightly less than threshold!
			WriteHit(hitpos[i], hitmom[i], cid,EDep, CumulativePathLength/nHits);
		}
		CumulativeEnergyDeposit=fmod(CumulativeEnergyDeposit,fThresholdEnergyDeposit);
		//G4cout<<"Interpolatclupd hits between "<<step->GetPreStepPoint()->GetPosition()<<" and "<<previousStep->GetPostStepPoint()->GetPosition()<<G4endl;
		displayHit(step);
		return nHits;
	}

	int DepositMultiplePtHits(G4Step* step) {
		if(fastSimulation)
			return DepositFastPtHits(step);
		//else regular->
		int nhits=0;
		while(CumulativeEnergyDeposit>fThresholdEnergyDeposit) {
			DepositPtHit(step, fThresholdEnergyDeposit);
			CumulativeEnergyDeposit-=fThresholdEnergyDeposit;
			++nhits;
		}
		//G4cout<<"would have deposited: "<<CumulativeEnergyDeposit/CLHEP::eV<<"/"<<fThresholdEnergyDeposit/CLHEP::eV<<"="<<CumulativeEnergyDeposit/fThresholdEnergyDeposit<<" hits"<<std::endl;
		//DepositPtHit(step,CumulativeEnergyDeposit);
		//	G4cout<<"deposited "<<nhits<<" and cummulativeEnergyDeposit is now "<<CumulativeEnergyDeposit/CLHEP::eV<<" eV"<<std::endl;
		return nhits;
	}

	void CumulateStep(G4Step *step) {
		const G4ThreeVector meanPosition = (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition()) / 2;
		const G4ThreeVector meanMomentum = (step->GetPreStepPoint()->GetMomentum() + step->GetPostStepPoint()->GetMomentum()) / 2;
		currentMeanPosition=meanPosition; currentMeanMomentum=meanMomentum;
		if(!CumulativeNumSteps) { //first step?
			firstStepPrePosition=step->GetPreStepPoint()->GetPosition();
			firstStepPreMomentum=step->GetPreStepPoint()->GetMomentum();
			firstStepTrackID= step->GetTrack()->GetTrackID();
		}
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

	void CumulatePtStep(G4Step *step) {
		CumulateStep(step);
	}

	void CumulateLowPtStep(G4Step *step) {
		CumulateStep(step);
	}

};


/// Initialization overload for specialization
template <> void Geant4SensitiveAction<TPCSDData>::initialize() {
	eventAction().callAtEnd(&m_userData,&TPCSDData::endEvent);

	declareProperty("TPCLowPtCut",              m_userData.Control.TPCLowPtCut );
	declareProperty("TPCLowPtStepLimit",        m_userData.Control.TPCLowPtStepLimit );
	declareProperty("TPCLowPtMaxHitSeparation", m_userData.Control.TPCLowPtMaxHitSeparation );

	m_userData.fThresholdEnergyDeposit = m_sensitive.energyCutoff()/CLHEP::keV; //set manually in geometry builder by kees, I do not know if this is the proper use of energyCutoff()
	m_userData.sensitive = this;

	IDDescriptor dsc = m_sensitive.idSpec() ;
	m_userData.layerField = dsc.field( "layer" ) ;


	G4cout << "limits of TPC are: "<< this->sensitiveDetector().limits().isValid() << std::endl;

	G4cout << "TPCSDAction: Threshold for Energy Deposit = " << m_userData.fThresholdEnergyDeposit / CLHEP::eV << " eV" << G4endl;

	const std::string TPCReadoutType= m_detDesc.constantAsString("TPC_readoutType");
	m_userData.pixelTPC=(TPCReadoutType=="pixel");
	const std::string TPCSimulationMode= m_detDesc.constantAsString("TPC_simulationMode");
	m_userData.fastSimulation=(TPCSimulationMode=="fast");

	if(not m_userData.pixelTPC) m_userData.fThresholdEnergyDeposit=0;

	m_userData.tpcData = m_detector.extension<dd4hep::rec::FixedPadSizeTPCData>();
	if(!m_userData.tpcData) G4cout<<"tpcData not found"<<std::endl;
}

/// Define collections created by this sensitive action object
template <> void Geant4SensitiveAction<TPCSDData>::defineCollections() {
	if(m_userData.pixelTPC)  {
		G4cout<<"defining m_collectionID with SimTrackerHitImpl"<<G4endl;
		m_collectionID = defineCollection<IMPL::SimTrackerHitImpl>(m_sensitive.readout().name());//define with SimTrackerHit directly
	} else {
		m_collectionID = defineCollection<Geant4Tracker::Hit>(m_sensitive.readout().name());//will later be converted to simtrackerhit
	}
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

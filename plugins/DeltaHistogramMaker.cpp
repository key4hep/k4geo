/*
 * DeltaHistogramMaker.cpp
 *
 *  Created on: Nov 8, 2016
 *      Author: cligtenb
 */

#include "DeltaHistogramMaker.h"

#include <functional>
#include <tuple>
#include <algorithm>
#include <exception>

#include "TMath.h"

#include "ParticleMap.h"

std::tuple<double, double, double> tupleFromVec(const G4ThreeVector  v) {
	return std::make_tuple(v.x(), v.y(), v.z());
}
TVector3 TVector3FromG4ThreeVector(const G4ThreeVector v) {
	return TVector3(v.x(), v.y(), v.z());
}

void DeltaTreeEntryType::SetMaxDiff(const G4ThreeVector v) {
	std::tie(_maxDiffx,_maxDiffy, _maxDiffz) = tupleFromVec(v);
}
void DeltaTreeEntryType::SetMaxPos(const G4ThreeVector v) {
	std::tie(_maxx,_maxy, _maxz) = tupleFromVec(v);
}
void DeltaTreeEntryType::SetMinPos(const G4ThreeVector v) {
	std::tie(_minx,_miny, _minz) = tupleFromVec(v);
}


DeltaHistogramMaker::DeltaHistogramMaker() {
	//set branches
	_deltaTree.Branch("deltax", &_deltaTreeEntry.x);
	_deltaTree.Branch("deltay", &_deltaTreeEntry.y);
	_deltaTree.Branch("deltaz", &_deltaTreeEntry.z);
	_deltaTree.Branch("deltapx", &_deltaTreeEntry.px);
	_deltaTree.Branch("deltapy", &_deltaTreeEntry.py);
	_deltaTree.Branch("deltapz", &_deltaTreeEntry.pz);
	_deltaTree.Branch("deltap", &_deltaTreeEntry.p);
	_deltaTree.Branch("deltaE", &_deltaTreeEntry.E);
	_deltaTree.Branch("EventNumber", &_deltaTreeEntry.eventNumber);
	_deltaTree.Branch("maxDiffx", &_deltaTreeEntry._maxDiffx);
	_deltaTree.Branch("maxDiffy", &_deltaTreeEntry._maxDiffy);
	_deltaTree.Branch("maxDiffz", &_deltaTreeEntry._maxDiffz);
	_deltaTree.Branch("maxPosx", &_deltaTreeEntry._maxx);
	_deltaTree.Branch("maxPosy", &_deltaTreeEntry._maxy);
	_deltaTree.Branch("maxPosz", &_deltaTreeEntry._maxz);
	_deltaTree.Branch("minPosx", &_deltaTreeEntry._minx);
	_deltaTree.Branch("minPosy", &_deltaTreeEntry._miny);
	_deltaTree.Branch("minPosz", &_deltaTreeEntry._minz);
	_deltaTree.Branch("numberOfHits", &_deltaTreeEntry.numberOfHits);
	_deltaTree.Branch("numberOfSteps", &_deltaTreeEntry.numSteps);
	_deltaTree.Branch("totalStepLength", &_deltaTreeEntry.cummulativePathLength);
	_deltaTree.Branch("hasChilds", &_deltaTreeEntry.hasChilds);
	_deltaTree.Branch("energyWithChilds", &_deltaTreeEntry.energyWithChilds);

	_deltaTree.Branch("distanceStartToEnd", &_deltaTreeEntry.distanceStartToEnd);
	_deltaTree.Branch("distanceTrackToEnd", &_deltaTreeEntry.distanceTrackToEnd);
	_deltaTree.Branch("distanceStartToTrackPoint", &_deltaTreeEntry.distanceStartToTrackPoint);

	_deltaTree.Branch("trackID", &_deltaTreeEntry.trackID);
	_deltaTree.Branch("parentID", &_deltaTreeEntry.parentID);


}

DeltaHistogramMaker::~DeltaHistogramMaker() {
//	  _deltaTree.Fill();
	ProcessCache();
	  _rootFile.Write();
}

//make a new vector from two vectors using the given function
G4ThreeVector CombineVectorsByComponent(const G4ThreeVector& v1, const G4ThreeVector& v2, std::function<double(double,double)> combiner) {
	G4ThreeVector result;
	result.setX(combiner(v1.x(), v2.x()));
	result.setY(combiner(v1.y(), v2.y()));
	result.setZ(combiner(v1.z(), v2.z()));
	return result;
}



void DeltaHistogramMaker::UpdateDifferenceFromStart(G4ThreeVector currentPoint) {
	G4ThreeVector distanceToStart(currentPoint-_startOfTrack);
	_maxDiff=CombineVectorsByComponent(_maxDiff,distanceToStart,[](double a, double b){a=fabs(a); b=fabs(b); return a>b ? a : b; });
	_deltaTreeEntry.SetMaxDiff(_maxDiff);
	_maxPos=CombineVectorsByComponent(_maxPos, distanceToStart, [](double a, double b){return a>b ? a : b; });
	_deltaTreeEntry.SetMaxPos(_maxPos);
	_minPos=CombineVectorsByComponent(_minPos, distanceToStart, [](double a, double b){return a<b ? a : b; });
	_deltaTreeEntry.SetMinPos(_minPos);

	if( distanceToStart.mag() > _deltaTreeEntry.distanceStartToEnd.Mag() ){
			_deltaTreeEntry.distanceStartToEnd=TVector3FromG4ThreeVector(distanceToStart);
			_deltaTreeEntry.distanceTrackToEnd=_deltaTreeEntry.distanceStartToEnd+_deltaTreeEntry.distanceStartToTrackPoint;
	}

}

void DeltaHistogramMaker::StartNewTrack(G4Step* step) {

	if(!first) _CachedEntries.push_back(_deltaTreeEntry);
//		if(!first) _deltaTree.Fill(); //fill at the start of new track
//		G4cout<<"deposited "<<_deltaTreeEntry.numberOfHits<<" hits"<<std::endl;

	_firstStepTrackID=step->GetTrack()->GetTrackID();
	_startOfTrack=step->GetPreStepPoint()->GetPosition();
	_maxDiff.set(0,0,0);
	_maxPos.set(0,0,0);
	_minPos.set(0,0,0);
	_deltaTreeEntry.SetPosition(step->GetPreStepPoint()->GetPosition());
	_deltaTreeEntry.setMomentum(step->GetPreStepPoint()->GetMomentum());
//		_deltaTreeEntry.E=step->GetTrack()->GetTotalEnergy(); wrong old energy!
	_deltaTreeEntry.E=step->GetPreStepPoint()->GetTotalEnergy();
	_deltaTreeEntry.eventNumber=_currentEventNumber;
	_deltaTreeEntry.numberOfHits=0;
	_deltaTreeEntry.numSteps=0;
	_deltaTreeEntry.cummulativePathLength=0;
	_deltaTreeEntry.trackID=step->GetTrack()->GetTrackID();
	_deltaTreeEntry.parentID=step->GetTrack()->GetParentID();
	_deltaTreeEntry.hasChilds=false;

	if(step->GetTrack()->GetTrackID()==1)  {
		_trackHits.clear();
		_closestPointOnTrack=trackHit();
	} else {
		struct unexpectedZeroTrackHits : std::exception {};
		if(!_trackHits.size()) throw unexpectedZeroTrackHits();
		_closestPointOnTrack=*std::min_element(_trackHits.begin(), _trackHits.end(),
				[this](trackHit& a, trackHit& b){ return (a.x-_startOfTrack).mag() <  (b.x-_startOfTrack).mag(); } );
	}
	_deltaTreeEntry.distanceStartToTrackPoint=TVector3FromG4ThreeVector(_startOfTrack-_closestPointOnTrack.x);//start - trackpoint
	_deltaTreeEntry.distanceStartToEnd.SetXYZ(0,0,0); //point - start
	_deltaTreeEntry.distanceTrackToEnd.SetXYZ(0,0,0);

	    G4cout<<"new track ID: "<<step->GetTrack()->GetTrackID()
	    		<<" for "<<getNameFromPID(step->GetTrack()->GetParticleDefinition()->GetPDGEncoding())
				<<" at "<<step->GetPreStepPoint()->GetPosition()
				<<" with momentum "<<int(step->GetPreStepPoint()->GetMomentum().mag()*1E3)<<" KeV"
				<<" created in "<<step->GetTrack()->GetCreatorModelName()
				<<" created by "<<step->GetTrack()->GetParentID()
				<<std::endl;

	//keep track of childs
	_childsOfTrack.insert(step->GetTrack()->GetParentID());
	summedEnergy[step->GetTrack()->GetTrackID()]=std::pair<int, double>(step->GetTrack()->GetParentID(), step->GetPreStepPoint()->GetKineticEnergy());

}

void DeltaHistogramMaker::Fill(G4Step* step, int eventNumber) {

	if(step->GetTrack()->GetTrackID()!=_firstStepTrackID ) {
		if(eventNumber!=_currentEventNumber) {
						ProcessCache();//process cache at new event
						_currentEventNumber=eventNumber;
		}
		StartNewTrack(step);
	}

	UpdateDifferenceFromStart(step->GetPostStepPoint()->GetPosition());
	_deltaTreeEntry.numSteps++;
	_deltaTreeEntry.cummulativePathLength+=step->GetStepLength();

	first=false;
}

void DeltaHistogramMaker::FillHit(G4ThreeVector position, G4ThreeVector momentum){
	++_deltaTreeEntry.numberOfHits;
	if(_deltaTreeEntry.trackID==1) {
		_trackHits.push_back( {position, momentum} );
	} else {
		G4ThreeVector distanceToStart(position-_startOfTrack);
		_distanceFromDeltaStart.Fill(distanceToStart.mag());
		_distanceFromDeltaStartZ.Fill(distanceToStart.z());
		_distanceFromDeltaStartXY.Fill(distanceToStart.perp());

		G4ThreeVector distanceToTrackPoint(position-_closestPointOnTrack.x);
		_distanceFromTrackPoint.Fill(distanceToTrackPoint.mag());
		_distanceFromStartToTrackPoint.Fill((_closestPointOnTrack.x-_startOfTrack).mag());

		G4ThreeVector distanceToTrack=distanceToTrackPoint.perpPart(_closestPointOnTrack.p);
		_distanceFromTrack.Fill(distanceToTrack.mag());
		_distanceFromTrackXY.Fill(distanceToTrack.perp());
		_distanceFromTrackZ.Fill(distanceToTrack.z());
	}
}

void DeltaHistogramMaker::ProcessCache() {
	//sum energy from childs
	for(auto rit=summedEnergy.rbegin(); rit!=summedEnergy.rend(); ++rit) {
		if(rit->second.first>1) {
			if(summedEnergy.find(rit->second.first)!=summedEnergy.end()) {
				summedEnergy.at(rit->second.first).second+=rit->second.second; //sum energy on parent
			} else {
				G4cout<<"Track "<<rit->second.first<<" was not found!"<<std::endl;
			}
		}
	}
	//loop over cached entries
	while( _CachedEntries.size() ) {
		_deltaTreeEntry= _CachedEntries.front();
		if(_childsOfTrack.count(_deltaTreeEntry.trackID)) //fill only entries with childs
			_deltaTreeEntry.hasChilds=true;
		else
			_deltaTreeEntry.hasChilds=false;

		if(_deltaTreeEntry.trackID>1)
			if(summedEnergy.find(_deltaTreeEntry.trackID)!=summedEnergy.end()) {
				_deltaTreeEntry.energyWithChilds=summedEnergy.at(_deltaTreeEntry.trackID).second;
			} else {
				G4cout<<"Track "<<_deltaTreeEntry.trackID<<" was not found!"<<std::endl;
			}
		_deltaTree.Fill();
		_CachedEntries.pop_front();
	}
	summedEnergy.clear();


}

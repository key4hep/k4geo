/*
 * DeltaHistogramMaker.h
 *
 *  Created on: Nov 8, 2016
 *      Author: cligtenb
 */

#ifndef SENSITIVE_DELTAHISTOGRAMMAKER_H_
#define SENSITIVE_DELTAHISTOGRAMMAKER_H_
#include <list>
#include <unordered_set>
#include <map>
#include <vector>
//geant4
#include "G4ThreeVector.hh"
#include "G4Step.hh"
//root
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TVector3.h"

struct DeltaTreeEntryType  { /*: TObject*/
public:
	double x,y,z;
	double px, py, pz, p;
	double E;
	int eventNumber;
    int numberOfHits;
	double _maxDiffx, _maxDiffy, _maxDiffz;//distance to start
	double _maxx, _maxy, _maxz;
	double _minx, _miny, _minz;
	double cummulativePathLength;
	int numSteps;
	bool hasChilds;
	double energyWithChilds;
	TVector3 distanceStartToTrackPoint, distanceStartToEnd, distanceTrackToEnd;

	int trackID, parentID;

	void inline SetPosition(const G4ThreeVector v) {x=v.x(); y=v.y();z=v.z();}
	void inline setMomentum(const G4ThreeVector v) {px=v.x(); py=v.y(); pz=v.z(); p=v.mag();}
	void SetMaxDiff(const G4ThreeVector);
	void SetMaxPos(const G4ThreeVector);
	void SetMinPos(const G4ThreeVector);


	//ClassDef(DeltaTreeEntryType, 1) not working???
};

// Work in progress, define a nice structure to do hists in directions
//struct OrthogonalHistSet {
//	OrthogonalHistSet(std::string name, std::string title, std::vector<double> range={0,1}, int nbins=100);
//	void Fill(const G4ThreeVector);
//	TH1D histmag;
//	TH1D histx, histy, histz;
//	TH1D histxy;
//};
//
//OrthogonalHistSet::OrthogonalHistSet(std::string name, std::string title, std::vector<double> range, int nbins) :
//	histmag( (name+"mag").c_str(), (title+"mag").c_str(), nbins, range.at(0), range.at(1) ) ,
//	histx( (name+"x").c_str(), (title+"x").c_str(), nbins, range.at(0), range.at(1) ),
//	histy( (name+"y").c_str(), (title+"y").c_str(), nbins, range.at(0), range.at(1) ),
//	histz( (name+"z").c_str(), (title+"z").c_str(), nbins, range.at(0), range.at(1) ),
//	histxy( (name+"xy").c_str(), (title+"xy").c_str(), nbins, range.at(0), range.at(1) )
//	{
//}
//
//void OrthogonalHistSet::Fill() {
//
//}

class DeltaHistogramMaker {
public:
	DeltaHistogramMaker();
	virtual ~DeltaHistogramMaker();

    //histograms
    TFile _rootFile {"TPCSDActionHistograms.root","RECREATE"};
    TTree _deltaTree {"deltaTree", "Tree with delta particles"};
    DeltaTreeEntryType _deltaTreeEntry;

    TH1D _distanceFromStartToTrackPoint {"distanceFromStartToTrackPoint", "Distance to the track point closest to start of the delta; distance [mm]; Hits", 100,0,1};
    TH1D _distanceFromTrackPoint {"distanceFromTrackPoint", "Distance to the track point closest to start of the delta; distance [mm]; Hits", 100,0,10};
    TH1D _distanceFromTrack {"distanceFromTrack", "Distance to the track; Distance to the track [mm]; Hits", 100,0,10};
    TH1D _distanceFromTrackXY {"distanceFromTrackXY", "Distance to the track in XY plane; Distance to the track in XY plane [mm]; Hits", 100,0,10};
    TH1D _distanceFromTrackZ {"distanceFromTrackZ", "Distance to the track in Z direction; Distance to the track in Z direction [mm]; Hits", 100,0,10};
    TH1D _distanceFromDeltaStart {"distanceFromDeltaStart", "Distance to the first step of delta; distance [mm]; Hits", 100,0,10};
    TH1D _distanceFromDeltaStartZ {"distanceFromDeltaStartZ", "Distance to the first step of delta in Z; distance [mm]; Hits", 100,0,10};
    TH1D _distanceFromDeltaStartXY {"distanceFromDeltaStartXY", "Distance to the first step of delta in XY; distance [mm]; Hits", 100,0,3};

    void Fill(G4Step* step, int eventNumber);
    void FillHit(G4ThreeVector position, G4ThreeVector momentum);
    void UpdateDifferenceFromStart(G4ThreeVector position);


private:
	int _firstStepTrackID=0;
	G4ThreeVector _startOfTrack;
	G4ThreeVector _maxDiff, _maxPos, _minPos;
	bool first=true;
	int _currentEventNumber=0;
	void StartNewTrack(G4Step* step);


	std::list<DeltaTreeEntryType> _CachedEntries; //cache entries in order to check them for childs
	std::unordered_set<int> _childsOfTrack; //keep an unordered set with all parents to count later
	void ProcessCache();

	std::map<int, std::pair<int, double> > summedEnergy; //map containing energy and parent

	struct trackHit { G4ThreeVector x, p; };
	trackHit _closestPointOnTrack;
	std::vector<trackHit> _trackHits;
};


#endif /* SENSITIVE_DELTAHISTOGRAMMAKER_H_ */

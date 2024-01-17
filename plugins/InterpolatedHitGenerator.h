/*
 * InterpolatedHitGenerator.h
 *
 *  Created on: Oct 12, 2016
 *      Author: cligtenb
 */

#ifndef INTERPOLATEDHITGENERATOR_H_
#define INTERPOLATEDHITGENERATOR_H_

#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>

//geant4
#include "G4ThreeVector.hh" //is actually a typedef of CLHEP::ThreeVector
#include "Randomize.hh"//also CLHEP

//root
#include "TRandom.h"
#include "TMath.h"


struct LinearVectorInterpolation {
	LinearVectorInterpolation( G4ThreeVector x1, G4ThreeVector x2) :
		_x1(x1),
		_x2(x2)
	{}
	G4ThreeVector _x1, _x2;
	G4ThreeVector operator()(double t) { //t is parameter between 0 and 1
		return _x1+t*(_x2-_x1);//interpolate linearly
	}
};

struct QuadraticVectorInterpolation {
	//interpolate using momentum at second point
	QuadraticVectorInterpolation( G4ThreeVector x1, G4ThreeVector x2, G4ThreeVector momentum) :
		_x1(x1),
		_x2(x2),
		_mom(momentum)
	{}
	G4ThreeVector _x1, _x2, _mom;
	G4ThreeVector operator()(double t) { //t is parameter between 0 and 1
		G4ThreeVector vT=(_x2-_x1).mag()*_mom.unit();
		return _x1+t*(2*(_x2-_x1)-vT)+t*t*(_x1-_x2+vT);//interpolate quadratically
	}

};

struct AngleInterpolation {
	//interpolate using momentum at angle between momenta
	AngleInterpolation( G4ThreeVector x1, G4ThreeVector x2, G4ThreeVector p1, G4ThreeVector p2) :
		_x1(x1),
		_x2(x2),
		_u1(p1.unit()),
		_u2(p2.unit())
	{}
	G4ThreeVector _x1, _x2, _u1, _u2;
	int sign(int x) {return x < 0 ? -1 : 1; }
	G4ThreeVector operator()(double t) { //t is parameter between 0 and 1
		G4ThreeVector pos {_x1};
		G4ThreeVector difference {_x2-_x1};
		//start with linear base
		G4ThreeVector relativePos=t*difference;
		//add dS factor in XY plane:
		double dphi=_u1.deltaPhi(_u2); //phi difference between vectors
		double L=difference.perp(); //point to point distance in XY plane
		double dS=L/2*sin(dphi/2)/2;
		//G4cout<<"ds="<<dS<<" L="<<L<<G4endl;
		G4ThreeVector dSPos;
		dSPos.setRhoPhiZ(
				4*abs(dS)*t*(1-t), //quadratic function in t
				difference.phi()-sign(dS)*M_PI/2, //orthogonal to difference
				0);
		return pos+relativePos+dSPos;
	}

};

//Produce a new hit at each location with a given chance;
struct FlatMultipleHitFunction {
	std::vector<double> recurrenceChance={0.1};
	double previous=CLHEP::HepRandom::getTheEngine()->flat();
	int multiplicity=0;
	FlatMultipleHitFunction(double p) : recurrenceChance(1,p) {}
	FlatMultipleHitFunction(std::vector<double> p) : recurrenceChance(p) {}
	double operator() (int) {
		if(CLHEP::HepRandom::getTheEngine()->flat()<recurrenceChance.at(multiplicity)) {
			if(unsigned(multiplicity)<recurrenceChance.size()-1) ++multiplicity;
			return previous;
		} else {
			multiplicity=0;
			return previous=CLHEP::HepRandom::getTheEngine()->flat();
		}
	}
};

struct FlatHitFunction : FlatMultipleHitFunction {
	FlatHitFunction() : FlatMultipleHitFunction(0) {};
};

struct LandauBasedHitFunction {
	double mu, sigma;
	int nhit=0;
	double t;
	LandauBasedHitFunction(double mu, double sigma) : mu(mu), sigma(sigma) {};
	double operator() (int) {
		if(nhit<=0) {
			while ( !(nhit=int(gRandom->Landau(mu, sigma))) ) {};//draw number until different from zero
			t=CLHEP::HepRandom::getTheEngine()->flat();
		}
		--nhit;
		return t;
	}
};

//generate a y=x like distribution for the number of hits
struct TriangularHitFunction {
	double triangleFraction;
	int nHitAtSame=0;
	double t=0;
	TriangularHitFunction(double triangleFraction=0.1) : triangleFraction(triangleFraction) {}
	double operator() (int nHitRemaining) {
			if(nHitAtSame<=0) {
				double y=0;//under treshold is poisson=always one hit
//				if(nHitRemaining<treshold) {
//				if(nHitRemaining>10) {
//					nHitAtSame=1;
//				}
//				const int startTrans=10; //treshold = slope length
//				if(gRandom->Rndm() < double(nHitRemaining-startTrans)/treshold ) {
				if(gRandom->Rndm() < triangleFraction) { // chance on generating a triangular dist
//				if(nHitRemaining>treshold){
//					y=1;
					y=TMath::Sqrt( gRandom->Rndm() );//y is a triangle between 0 and 1
//					y= gRandom->Rndm() ;//flat instead!
				}
				const int minHits=0; //fix to 1?
				nHitAtSame=int(minHits+(nHitRemaining-minHits)*y);
				t=gRandom->Rndm();
			}
			--nHitAtSame;
			return t;
	}
};

class InterpolatedHitGenerator {
public:
	InterpolatedHitGenerator(
			std::function<double(int)> hitDistributionFunction, //e.g. [](int) { return CLHEP::HepRandom::getTheEngine()->flat(); }
			std::function<G4ThreeVector(double)> hitLocationFunction); //e.g. LinearVectorInterpolation(x1,x2);
	std::vector<G4ThreeVector> generateHits(int nHits);
	std::vector<G4ThreeVector> interpolateHits(std::function<G4ThreeVector>);
	std::vector<G4ThreeVector> generateFromFunction(
			const std::function<G4ThreeVector(double)>& vectorFunction);

protected:
	std::function<double(int)> _getHitDistribution; //return parameter between 0.-1, for location on track, takes number of remaining hits as parameter
	std::function<G4ThreeVector(double)> _getHitLocation;
	std::vector<double> _param; //param saves parameter used to generate hits such that the corresponding momentum can later be found.

};


InterpolatedHitGenerator::InterpolatedHitGenerator(
		std::function<double(int)> hitDistributionFunction,
		std::function<G4ThreeVector(double)> hitLocationFunction) :
		_getHitDistribution(hitDistributionFunction),
		_getHitLocation(hitLocationFunction)
	{}

std::vector<G4ThreeVector> InterpolatedHitGenerator::generateFromFunction(
		const std::function<G4ThreeVector(double)>& vectorFunction) {
	//generate hits
	std::vector<G4ThreeVector> hits;
	hits.reserve(_param.size());
	for (double t : _param) {
		hits.push_back(vectorFunction(t));
	}
	return hits;
}

std::vector<G4ThreeVector> InterpolatedHitGenerator::generateHits(int nHits) {
	//first generate the parameter, such that hits can be easily placed in order
	_param.reserve(nHits);
	for(int i=0; i<nHits; i++) {
		_param.push_back( _getHitDistribution(nHits-i) );
	}
	//sort
	std::sort(_param.begin(), _param.end());
	//generate hits
	std::vector<G4ThreeVector> hits = generateFromFunction(_getHitLocation);
	return hits;
}

#endif /* INTERPOLATEDHITGENERATOR_H_ */

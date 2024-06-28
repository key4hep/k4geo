#include "detectorSegmentations/DRparamBarrel_k4geo.h"

#include "Math/GenVector/RotationZYX.h"

#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

DRparamBarrel_k4geo::DRparamBarrel_k4geo() {
  fIsRHS = 0;
  fPhiZRot = 0.;
  fInnerX = 0.;
  fTowerH = 0.;
  fNumZRot = 0;
  fDeltaTheta = 0.;
  fThetaOfCenter = 0.;
  fCurrentInnerR = 0.;
  fPhiZRot = 0;
  fCurrentCenter = TVector3();
  fV1 = TVector3();
  fV2 = TVector3();
  fV3 = TVector3();
  fV4 = TVector3();
  fSipmHeight = 0.;
  fCurrentInnerHalf = 0.;
  fCurrentOuterHalf = 0.;
  fFilled = false;
  fFinalized = false;
}

DRparamBarrel_k4geo::~DRparamBarrel_k4geo() {}

void DRparamBarrel_k4geo::init() {
  fCurrentInnerR = fInnerX/std::cos(fThetaOfCenter);
  double trnsLength = fTowerH/2.+fCurrentInnerR;
  fCurrentCenter =  TVector3(std::cos(fThetaOfCenter)*trnsLength,0.,std::sin(fThetaOfCenter)*trnsLength);

  fCurrentInnerHalf = fCurrentInnerR*std::tan(fDeltaTheta/2.);
  fCurrentOuterHalf = (fCurrentInnerR+fTowerH)*std::tan(fDeltaTheta/2.);
  fCurrentOuterHalfSipm = (fCurrentInnerR+fTowerH+fSipmHeight)*std::tan(fDeltaTheta/2.);

  fV1 = TVector3(
    std::cos(fThetaOfCenter)*fCurrentInnerR+std::sin(fThetaOfCenter)*fCurrentInnerR*std::tan(fDeltaTheta/2.),
    0.,
    std::sin(fThetaOfCenter)*fCurrentInnerR-std::cos(fThetaOfCenter)*fCurrentInnerR*std::tan(fDeltaTheta/2.)
  );

  fV2 = TVector3(
    std::cos(fThetaOfCenter)*(fCurrentInnerR+fTowerH)+std::sin(fThetaOfCenter)*(fCurrentInnerR+fTowerH)*std::tan(fDeltaTheta/2.),
    0.,
    std::sin(fThetaOfCenter)*(fCurrentInnerR+fTowerH)-std::cos(fThetaOfCenter)*(fCurrentInnerR+fTowerH)*std::tan(fDeltaTheta/2.)
  );

  fV3 = TVector3(
    std::cos(fThetaOfCenter)*fCurrentInnerR-std::sin(fThetaOfCenter)*fCurrentInnerR*std::tan(fDeltaTheta/2.),
    0.,
    std::sin(fThetaOfCenter)*fCurrentInnerR+std::cos(fThetaOfCenter)*fCurrentInnerR*std::tan(fDeltaTheta/2.)
  );

  fV4 = TVector3(
    std::cos(fThetaOfCenter)*(fCurrentInnerR+fTowerH)-std::sin(fThetaOfCenter)*(fCurrentInnerR+fTowerH)*std::tan(fDeltaTheta/2.),
    0.,
    std::sin(fThetaOfCenter)*(fCurrentInnerR+fTowerH)+std::cos(fThetaOfCenter)*(fCurrentInnerR+fTowerH)*std::tan(fDeltaTheta/2.)
  );

  fV2sipm = TVector3(
    std::cos(fThetaOfCenter)*(fCurrentInnerR+fTowerH+fSipmHeight)+std::sin(fThetaOfCenter)*(fCurrentInnerR+fTowerH+fSipmHeight)*std::tan(fDeltaTheta/2.),
    0.,
    std::sin(fThetaOfCenter)*(fCurrentInnerR+fTowerH+fSipmHeight)-std::cos(fThetaOfCenter)*(fCurrentInnerR+fTowerH+fSipmHeight)*std::tan(fDeltaTheta/2.)
  );

  fV4sipm = TVector3(
    std::cos(fThetaOfCenter)*(fCurrentInnerR+fTowerH+fSipmHeight)-std::sin(fThetaOfCenter)*(fCurrentInnerR+fTowerH+fSipmHeight)*std::tan(fDeltaTheta/2.),
    0.,
    std::sin(fThetaOfCenter)*(fCurrentInnerR+fTowerH+fSipmHeight)+std::cos(fThetaOfCenter)*(fCurrentInnerR+fTowerH+fSipmHeight)*std::tan(fDeltaTheta/2.)
  );

  if (!fFilled) {
    fDeltaThetaVec.push_back(fDeltaTheta);
    fThetaOfCenterVec.push_back(fThetaOfCenter);
  }
}

void DRparamBarrel_k4geo::SetDeltaThetaByTowerNo(int signedTowerNo, int) {
  if (!fFilled) throw std::runtime_error("Attempt to set by tower num while barrel parameter is not filled!");

  fDeltaTheta = fDeltaThetaVec.at( unsignedTowerNo(signedTowerNo) );
}

void DRparamBarrel_k4geo::SetThetaOfCenterByTowerNo(int signedTowerNo, int) {
  if (!fFilled) throw std::runtime_error("Attempt to set by tower num while barrel parameter is not filled!");

  fThetaOfCenter = fThetaOfCenterVec.at( unsignedTowerNo(signedTowerNo) );
}

}
}

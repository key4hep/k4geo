#include "DRparamBase.h"

#include "Math/GenVector/RotationZYX.h"

#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

DRparamBase::DRparamBase() {
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
  fCurrentTowerNum = 0;
  fFilled = false;
  fFinalized = false;
}

DRparamBase::~DRparamBase() {}

dd4hep::RotationZYX DRparamBase::GetRotationZYX(int numPhi) {
  double numPhi_ = (double)numPhi;
  double xRot = fIsRHS ? -fThetaOfCenter : fThetaOfCenter;
  double zRot = fIsRHS ? -M_PI/2. : M_PI/2.;
  dd4hep::RotationZYX rot = dd4hep::RotationZYX(zRot, M_PI/2.+xRot, 0.);
  ROOT::Math::RotationZ rotZ = ROOT::Math::RotationZ(numPhi_*fPhiZRot);
  rot = rotZ*rot;

  return rot;
}

dd4hep::Position DRparamBase::GetTowerPos(int numPhi) {
  double numPhi_ = (double)numPhi;
  double x = std::cos(numPhi_*fPhiZRot)*fCurrentCenter.X();
  double y = std::sin(numPhi_*fPhiZRot)*fCurrentCenter.X();
  double z = fIsRHS ? fCurrentCenter.Z() : -fCurrentCenter.Z();
  dd4hep::Position pos = dd4hep::Position(x,y,z);

  return pos;
}

dd4hep::Position DRparamBase::GetAssemblePos(int numPhi) {
  double numPhi_ = (double)numPhi;
  double x = std::cos(numPhi_*fPhiZRot)*fCurrentCenter.X()*(fCurrentCenter.Mag()+fSipmHeight/2.)/fCurrentCenter.Mag();
  double y = std::sin(numPhi_*fPhiZRot)*fCurrentCenter.X()*(fCurrentCenter.Mag()+fSipmHeight/2.)/fCurrentCenter.Mag();
  double z_abs = fCurrentCenter.Z()*(fCurrentCenter.Mag()+fSipmHeight/2.)/fCurrentCenter.Mag();
  double z = fIsRHS ? z_abs : -z_abs;
  dd4hep::Position pos = dd4hep::Position(x,y,z);

  return pos;
}

dd4hep::Position DRparamBase::GetSipmLayerPos(int numPhi) {
  double numPhi_ = (double)numPhi;
  double x = std::cos(numPhi_*fPhiZRot)*fCurrentCenter.X()*(fCurrentCenter.Mag()+fTowerH/2.+fSipmHeight/2.)/fCurrentCenter.Mag();
  double y = std::sin(numPhi_*fPhiZRot)*fCurrentCenter.X()*(fCurrentCenter.Mag()+fTowerH/2.+fSipmHeight/2.)/fCurrentCenter.Mag();
  double z_abs = fCurrentCenter.Z()*(fCurrentCenter.Mag()+fTowerH/2.+fSipmHeight/2.)/fCurrentCenter.Mag();
  double z = fIsRHS ? z_abs : -z_abs;
  dd4hep::Position pos = dd4hep::Position(x,y,z);

  return pos;
}

dd4hep::Transform3D DRparamBase::GetTransform3D(int numPhi) {
  auto rot = GetRotationZYX(numPhi);
  auto pos = GetTowerPos(numPhi);

  return dd4hep::Transform3D(rot,pos);
}

dd4hep::Transform3D DRparamBase::GetAssembleTransform3D(int numPhi) {
  auto rot = GetRotationZYX(numPhi);
  auto pos = GetAssemblePos(numPhi);

  return dd4hep::Transform3D(rot,pos);
}

dd4hep::Transform3D DRparamBase::GetSipmTransform3D(int numPhi) {
  auto rot = GetRotationZYX(numPhi);
  auto pos = GetSipmLayerPos(numPhi);

  return dd4hep::Transform3D(rot,pos);
}

}
}

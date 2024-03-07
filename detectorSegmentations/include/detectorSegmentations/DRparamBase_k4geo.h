#ifndef DETSEGMENTATION_DRPARAMBASE_H
#define DETSEGMENTATION_DRPARAMBASE_H

#include "TVector3.h"
#include "DD4hep/DetFactoryHelper.h"

#include <vector>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {
  class DRparamBase_k4geo {
  public:
    DRparamBase_k4geo();
    virtual ~DRparamBase_k4geo();

    void SetIsRHS(bool isRHS) { fIsRHS = isRHS; }
    void SetInnerX(double innerX) { fInnerX = innerX; }
    void SetTowerH(double towerH) { fTowerH = towerH; }
    void SetNumZRot(int num) { fNumZRot = num; fPhiZRot = 2*M_PI/(double)num; }
    void SetDeltaTheta(double theta) { fDeltaTheta = theta; }
    void SetThetaOfCenter(double theta) { fThetaOfCenter = theta; }
    void SetSipmHeight(double SipmHeight) { fSipmHeight = SipmHeight; }

    bool GetIsRHS() { return fIsRHS; }
    double GetCurrentInnerR() { return fCurrentInnerR; }
    double GetTowerH() { return fTowerH; }
    double GetSipmHeight() { return fSipmHeight; }
    double GetH1() { return fCurrentInnerHalf; }
    double GetBl1() { return fV3.X()*std::tan(fPhiZRot/2.); }
    double GetTl1() { return fV1.X()*std::tan(fPhiZRot/2.); }
    double GetH2() { return fCurrentOuterHalf; }
    double GetBl2() { return fV4.X()*std::tan(fPhiZRot/2.); }
    double GetTl2() { return fV2.X()*std::tan(fPhiZRot/2.); }

    double GetH2sipm() { return fCurrentOuterHalfSipm; }
    double GetBl2sipm() { return fV4sipm.X()*std::tan(fPhiZRot/2.); }
    double GetTl2sipm() { return fV2sipm.X()*std::tan(fPhiZRot/2.); }

    dd4hep::RotationZYX GetRotationZYX(int numPhi);
    dd4hep::Position GetTowerPos(int numPhi);
    dd4hep::Position GetAssemblePos(int numPhi);
    dd4hep::Position GetSipmLayerPos(int numPhi);

    dd4hep::Transform3D GetTransform3D(int numPhi);
    dd4hep::Transform3D GetAssembleTransform3D(int numPhi);
    dd4hep::Transform3D GetSipmTransform3D(int numPhi);

    int signedTowerNo(int unsignedTowerNo) { return fIsRHS ? unsignedTowerNo : -unsignedTowerNo-1; }
    int unsignedTowerNo(int signedTowerNo) { return signedTowerNo >= 0 ? signedTowerNo : -signedTowerNo-1; }

    virtual void SetDeltaThetaByTowerNo(int , int ) {}
    virtual void SetThetaOfCenterByTowerNo(int , int ) {}
    void SetIsRHSByTowerNo(int signedTowerNo) { fIsRHS = ( signedTowerNo >=0 ? true : false ); }

    int GetTotTowerNum() { return fTotNum; }
    void SetTotTowerNum(int totNum) { fTotNum = totNum; }

    int GetCurrentTowerNum() { return fCurrentTowerNum; }
    void SetCurrentTowerNum(int numEta) { fCurrentTowerNum = numEta; }

    virtual void init() {};
    void filled() { fFilled = true; }
    void finalized() { fFinalized = true; }
    bool IsFinalized() { return fFinalized; }

  protected:
    bool fIsRHS;
    double fPhiZRot;
    double fInnerX;
    double fTowerH;
    int fNumZRot;
    double fDeltaTheta;
    double fThetaOfCenter;
    double fCurrentInnerR;
    TVector3 fCurrentCenter;
    TVector3 fV1;
    TVector3 fV2;
    TVector3 fV3;
    TVector3 fV4;
    TVector3 fV2sipm;
    TVector3 fV4sipm;
    double fSipmHeight;

    double fCurrentInnerHalf;
    double fCurrentOuterHalf;
    double fCurrentOuterHalfSipm;

    int fTotNum;
    int fCurrentTowerNum;
    std::vector<double> fDeltaThetaVec;
    std::vector<double> fThetaOfCenterVec;
    bool fFilled;
    bool fFinalized;
  };
}
}

#endif

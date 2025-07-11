#ifndef DETSEGMENTATION_DRPARAMBASE_H
#define DETSEGMENTATION_DRPARAMBASE_H

#include "DD4hep/DetFactoryHelper.h"
#include "TVector3.h"

#include <cmath>
#include <vector>

namespace dd4hep {
namespace DDSegmentation {
  class DRparamBase_k4geo {
  public:
    DRparamBase_k4geo();
    virtual ~DRparamBase_k4geo();

    void SetIsRHS(bool isRHS) { fIsRHS = isRHS; }
    void SetInnerX(double innerX) { fInnerX = innerX; }
    void SetTowerH(double towerH) { fTowerH = towerH; }
    void SetNumZRot(int num) {
      fNumZRot = num;
      fPhiZRot = 2 * M_PI / (double)num;
    }
    void SetDeltaTheta(double theta) { fDeltaTheta = theta; }
    void SetThetaOfCenter(double theta) { fThetaOfCenter = theta; }
    void SetSipmHeight(double SipmHeight) { fSipmHeight = SipmHeight; }

    bool GetIsRHS() { return fIsRHS; }
    int GetNumZRot() { return fNumZRot; }
    double GetCurrentInnerR() { return fCurrentInnerR; }
    double GetTowerH() { return fTowerH; }
    double GetSipmHeight() { return fSipmHeight; }
    double GetH1() { return fCurrentInnerHalf; }
    double GetBl1() { return fV3.X() * std::tan(fPhiZRot / 2.); }
    double GetTl1() { return fV1.X() * std::tan(fPhiZRot / 2.); }
    double GetH2() { return fCurrentOuterHalf; }
    double GetBl2() { return fV4.X() * std::tan(fPhiZRot / 2.); }
    double GetTl2() { return fV2.X() * std::tan(fPhiZRot / 2.); }

    double GetH2sipm() { return fCurrentOuterHalfSipm; }
    double GetBl2sipm() { return fV4sipm.X() * std::tan(fPhiZRot / 2.); }
    double GetTl2sipm() { return fV2sipm.X() * std::tan(fPhiZRot / 2.); }

    dd4hep::RotationZYX GetRotationZYX(int numPhi);
    dd4hep::Position GetTowerPos(int numPhi);
    dd4hep::Position GetAssemblePos(int numPhi);
    dd4hep::Position GetSipmLayerPos(int numPhi);

    dd4hep::Transform3D GetTransform3D(int numPhi);
    dd4hep::Transform3D GetAssembleTransform3D(int numPhi);
    dd4hep::Transform3D GetSipmTransform3D(int numPhi);

    int signedTowerNo(int unsignedTowerNo) const { return fIsRHS ? unsignedTowerNo : -unsignedTowerNo - 1; }
    int unsignedTowerNo(int signedTowerNo) const { return signedTowerNo >= 0 ? signedTowerNo : -signedTowerNo - 1; }

    virtual void SetDeltaThetaByTowerNo(int, int) {}
    virtual void SetThetaOfCenterByTowerNo(int, int) {}
    void SetIsRHSByTowerNo(int signedTowerNo) { fIsRHS = (signedTowerNo >= 0 ? true : false); }

    int GetTotTowerNum() { return fTotNum; }
    void SetTotTowerNum(int totNum) { fTotNum = totNum; }

    int GetCurrentTowerNum() { return fCurrentTowerNum; }
    void SetCurrentTowerNum(int numEta) { fCurrentTowerNum = numEta; }

    // the size of the neighborhood (radius) in units of fiber-fiber distance (i.e. grid size)
    void SetNeighborSize(double nsize) { fNeighborSize = nsize; }
    double GetNeighborSize() { return fNeighborSize; }
    // the margin in units of fiber-fiber distance when incorporating satelite cells at the edge of the tower
    void SetMargin(int margin) { fMargin = margin; }
    int GetMargin() { return fMargin; }

    virtual void init() {};
    void filled() { fFilled = true; }
    void finalized() { fFinalized = true; }
    bool IsFinalized() { return fFinalized; }

    // store information of which fibers have full length or not
    // full length fibers within rmin <= n_row <= rmax and cmin <= n_column <= cmax
    struct fullLengthFibers {
    public:
      fullLengthFibers(int rmin_, int rmax_, int cmin_, int cmax_)
          : rmin(rmin_), rmax(rmax_), cmin(cmin_), cmax(cmax_) {}

      fullLengthFibers() // default constructor
          : rmin(0), rmax(0), cmin(0), cmax(0) {}

      int rmin; // min n_row with full length fibers
      int rmax; // max n_row
      int cmin; // min n_column with full length fibers
      int cmax; // max n_column
    };

    fullLengthFibers GetFullLengthFibers(int numEta) const { return fFullLengthFibers.at(unsignedTowerNo(numEta)); }
    void SetFullLengthFibers(int rmin, int rmax, int cmin, int cmax);

    // length of fibers that don't have full length
    // it is a map containing (row, col) as a key and the length of the fiber as a value
    // and each tower in eta has one of it
    struct shortFibers {
    public:
      shortFibers(const double towerH) : fTowerH(towerH) {}
      ~shortFibers() = default;

      void addShortFibers(const int row, const int col, const double len) {
        m_fiberLengths_.insert(std::make_pair(std::make_pair(row, col), len));
      }
      double retrieveFiberLength(const int row, const int col) const;

    private:
      std::map<std::pair<int, int>, double> m_fiberLengths_;
      const double fTowerH;
    };

    shortFibers GetShortFibers(int numEta) const { return fShortFibers.at(unsignedTowerNo(numEta)); }
    void SetShortFibers(const shortFibers& input);

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

    // parameters for neighbour finding algorithm
    double fNeighborSize;
    int fMargin;

    int fTotNum;
    int fCurrentTowerNum;
    std::vector<double> fDeltaThetaVec;
    std::vector<double> fThetaOfCenterVec;
    std::map<int, fullLengthFibers> fFullLengthFibers;
    std::map<int, shortFibers> fShortFibers;
    bool fFilled;
    bool fFinalized;
  };
} // namespace DDSegmentation
} // namespace dd4hep

#endif

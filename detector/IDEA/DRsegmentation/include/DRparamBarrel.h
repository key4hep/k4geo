#ifndef DRparamBarrel_h
#define DRparamBarrel_h 1

#include "DRparamBase.h"

#include "TVector3.h"
#include "DD4hep/DetFactoryHelper.h"

#include <vector>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {
  class DRparamBarrel : public DRparamBase {
  public:
    DRparamBarrel();
    virtual ~DRparamBarrel();

    virtual void SetDeltaThetaByTowerNo(int signedTowerNo, int) override;
    virtual void SetThetaOfCenterByTowerNo(int signedTowerNo, int) override;

    virtual void init() override;
  };
}
}

#endif

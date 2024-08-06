#ifndef DETSEGMENTATION_DRPARAMBARREL_H
#define DETSEGMENTATION_DRPARAMBARREL_H

#include "detectorSegmentations/DRparamBase_k4geo.h"

#include "TVector3.h"
#include "DD4hep/DetFactoryHelper.h"

#include <vector>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {
  class DRparamBarrel_k4geo : public DRparamBase_k4geo {
  public:
    DRparamBarrel_k4geo();
    virtual ~DRparamBarrel_k4geo();

    virtual void SetDeltaThetaByTowerNo(int signedTowerNo, int) override;
    virtual void SetThetaOfCenterByTowerNo(int signedTowerNo, int) override;

    virtual void init() override;
  };
}
}

#endif

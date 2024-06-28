#ifndef DETSEGMENTATION_DRPARAMENDCAP_H
#define DETSEGMENTATION_DRPARAMENDCAP_H

#include "detectorSegmentations/DRparamBase_k4geo.h"

#include "TVector3.h"
#include "DD4hep/DetFactoryHelper.h"

#include <vector>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {
  class DRparamEndcap_k4geo : public DRparamBase_k4geo {
  public:
    DRparamEndcap_k4geo();
    virtual ~DRparamEndcap_k4geo();

    virtual void SetDeltaThetaByTowerNo(int signedTowerNo, int BEtrans) override;
    virtual void SetThetaOfCenterByTowerNo(int signedTowerNo, int BEtrans) override;

    virtual void init() override;
  };
}
}

#endif

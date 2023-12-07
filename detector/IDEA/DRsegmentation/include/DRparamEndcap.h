#ifndef DRparamEndcap_h
#define DRparamEndcap_h 1

#include "DRparamBase.h"

#include "TVector3.h"
#include "DD4hep/DetFactoryHelper.h"

#include <vector>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {
  class DRparamEndcap : public DRparamBase {
  public:
    DRparamEndcap();
    virtual ~DRparamEndcap();

    virtual void SetDeltaThetaByTowerNo(int signedTowerNo, int BEtrans) override;
    virtual void SetThetaOfCenterByTowerNo(int signedTowerNo, int BEtrans) override;

    virtual void init() override;
  };
}
}

#endif

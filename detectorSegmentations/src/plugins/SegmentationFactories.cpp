#include "DD4hep/Factories.h"
#include "DD4hep/detail/SegmentationsInterna.h"

namespace {
template <typename T>
dd4hep::SegmentationObject* create_segmentation(const dd4hep::BitFieldCoder* decoder) {
  return new dd4hep::SegmentationWrapper<T>(decoder);
}
}

#include "detectorSegmentations/GridEta_k4geo.h"
DECLARE_SEGMENTATION(GridEta_k4geo, create_segmentation<dd4hep::DDSegmentation::GridEta_k4geo>)

#include "detectorSegmentations/GridTheta_k4geo.h"
DECLARE_SEGMENTATION(GridTheta_k4geo, create_segmentation<dd4hep::DDSegmentation::GridTheta_k4geo>)

#include "detectorSegmentations/FCCSWGridPhiTheta.h"
DECLARE_SEGMENTATION(FCCSWGridPhiTheta, create_segmentation<dd4hep::DDSegmentation::FCCSWGridPhiTheta>)

#include "detectorSegmentations/FCCSWGridPhiEta.h"
DECLARE_SEGMENTATION(FCCSWGridPhiEta, create_segmentation<dd4hep::DDSegmentation::FCCSWGridPhiEta>)

#include "detectorSegmentations/GridRPhiEta.h"
DECLARE_SEGMENTATION(GridRPhiEta, create_segmentation<dd4hep::DDSegmentation::GridRPhiEta>)

#include "detectorSegmentations/GridSimplifiedDriftChamber.h"
DECLARE_SEGMENTATION(GridSimplifiedDriftChamber, create_segmentation<dd4hep::DDSegmentation::GridSimplifiedDriftChamber>)


#include "DD4hep/Factories.h"
#include "DD4hep/detail/SegmentationsInterna.h"

namespace {
template <typename T>
dd4hep::SegmentationObject* create_segmentation(const dd4hep::BitFieldCoder* decoder) {
  return new dd4hep::SegmentationWrapper<T>(decoder);
}
} // namespace

#include "detectorSegmentations/GridEta_k4geo.h"
DECLARE_SEGMENTATION(GridEta_k4geo, create_segmentation<dd4hep::DDSegmentation::GridEta_k4geo>)

#include "detectorSegmentations/GridTheta_k4geo.h"
DECLARE_SEGMENTATION(GridTheta_k4geo, create_segmentation<dd4hep::DDSegmentation::GridTheta_k4geo>)

#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
DECLARE_SEGMENTATION(FCCSWGridPhiTheta_k4geo, create_segmentation<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo>)

#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"
DECLARE_SEGMENTATION(FCCSWGridModuleThetaMerged_k4geo,
                     create_segmentation<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo>)

#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"
DECLARE_SEGMENTATION(FCCSWGridPhiEta_k4geo, create_segmentation<dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo>)

#include "detectorSegmentations/GridRPhiEta_k4geo.h"
DECLARE_SEGMENTATION(GridRPhiEta_k4geo, create_segmentation<dd4hep::DDSegmentation::GridRPhiEta_k4geo>)

#include "detectorSegmentations/GridSimplifiedDriftChamber_k4geo.h"
DECLARE_SEGMENTATION(GridSimplifiedDriftChamber_k4geo,
                     create_segmentation<dd4hep::DDSegmentation::GridSimplifiedDriftChamber_k4geo>)

#include "detectorSegmentations/GridDRcalo_k4geo.h"
DECLARE_SEGMENTATION(GridDRcalo_k4geo, create_segmentation<dd4hep::DDSegmentation::GridDRcalo_k4geo>)

#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
DECLARE_SEGMENTATION(FCCSWEndcapTurbine_k4geo, create_segmentation<dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo>)

#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
DECLARE_SEGMENTATION(FCCSWHCalPhiTheta_k4geo, create_segmentation<dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo>)

#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"
DECLARE_SEGMENTATION(FCCSWHCalPhiRow_k4geo, create_segmentation<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo>)

#include "detectorSegmentations/SCEPCal_MainSegmentation_k4geo.h"
DECLARE_SEGMENTATION(SCEPCal_MainSegmentation_k4geo,
                     create_segmentation<dd4hep::DDSegmentation::SCEPCal_MainSegmentation_k4geo>)

#include "detectorSegmentations/SCEPCal_TimingSegmentation_k4geo.h"
DECLARE_SEGMENTATION(SCEPCal_TimingSegmentation_k4geo,
                     create_segmentation<dd4hep::DDSegmentation::SCEPCal_TimingSegmentation_k4geo>)
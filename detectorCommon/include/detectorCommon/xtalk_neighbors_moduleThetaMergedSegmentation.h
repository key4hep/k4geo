#ifndef DETECTORCOMMON_XTALKMODULETHETAMERGED_H
#define DETECTORCOMMON_XTALKMODULETHETAMERGED_H

// k4geo
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

// DD4hep
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Segmentations.h"
#include "DDSegmentation/BitFieldCoder.h"

namespace det {
namespace crosstalk {

  /** The crosstalk neighbours finding function for the readout with module and theta merged cells
   *  Compared to the `neighbours` function from DetUtils, it needs a reference to the segmentation class to
   *  access the number of merged cells per layer. The input parameters and return value are the following:
   *   @param[in] aSeg Reference to the segmentation object.
   *   @param[in] aDecoder Handle to the bitfield decoder.
   *   @param[in] aFieldNames Names of the fields for which crosstalk neighbours are found.
   *   @param[in] aFieldExtremes Minimal and maximal values for the fields.
   *   @param[in] aCellId ID of cell.
   *   @param[in] aXtalkCoefRadial radial crosstalk coefficient.
   *   @param[in] aXtalkCoefTheta theta crosstalk coefficient.
   *   @param[in] aXtalkCoefDiagonal diagonal crosstalk coefficient.
   *   @param[in] aXtalkCoefTower crosstalk coefficient in the same theta tower.
   *   return Vector of neighbours and their crosstalk coefficients.
   */
  std::vector<std::pair<uint64_t, double>>
  getNeighboursModuleThetaMerged(const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo& aSeg,
                                 const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                 const std::vector<std::string>& aFieldNames,
                                 const std::vector<std::vector<std::pair<int, int>>>& aFieldExtremes, uint64_t aCellId,
                                 const double aXtalkCoefRadial = 0.7e-2, const double aXtalkCoefTheta = 0.2e-2,
                                 const double aXtalkCoefDiagonal = 0.04e-2, const double aXtalkCoefTower = 0.1e-2);

  // debug: return cell layer/module/theta indices
  std::vector<int> getCellIndices(const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo& aSeg,
                                  const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                  const std::vector<std::string>& aFieldNames, uint64_t aCellId);
} // namespace crosstalk
} // namespace det
#endif /* DETCOMMON_XTALKMODULETHETAMERGED_H */

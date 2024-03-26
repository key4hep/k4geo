#ifndef DETECTORCOMMON_XTALK_H
#define DETECTORCOMMON_XTALK_H

// k4geo
#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

// DD4hep
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Segmentations.h"
#include "DDSegmentation/BitFieldCoder.h"
// Include GridPhiEta from dd4hep once eta calculation is fixed
//#include "DDSegmentation/GridPhiEta.h"
#include "DDSegmentation/CartesianGridXY.h"
#include "DDSegmentation/CartesianGridXYZ.h"
#include "DDSegmentation/PolarGridRPhi.h"

// Geant
#include "G4Step.hh"

// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

#include "TGeoManager.h"

/** Given a XML element with several daughters with the same name, e.g.
 <detector> <layer name="1" /> <layer name="2"> </detector>
 this method returns the first daughter of type nodeName whose attribute has a given value
 e.g. returns <layer name="2"/> when called with (detector, "layer", "name", "1") */
namespace det {
namespace xtalk {

 /** Special version of the crosstalk neighbours function for the readout with module and theta merged cells
 *  Compared to the standard version, it needs a reference to the segmentation class to
 *  access the number of merged cells per layer. The other parameters and return value are the same
 *   @param[in] aSeg Reference to the segmentation object.
 *   @param[in] aDecoder Handle to the bitfield decoder.
 *   @param[in] aFieldNames Names of the fields for which neighbours are found.
 *   @param[in] aFieldExtremes Minimal and maximal values for the fields.
 *   @param[in] aCellId ID of cell.
 *   return Vector of neighbours and cross talk coefficients.
 */
std::vector<std::pair<uint64_t, double>> xtalk_neighbours_ModuleThetaMerged(const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo& aSeg,
                                 const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                 const std::vector<std::string>& aFieldNames,
                                 const std::vector<std::vector<std::pair<int, int>>>& aFieldExtremes,
                                 uint64_t aCellId);
// debug: return cell position
std::vector<int> xtalk_get_cell_position(const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo& aSeg,
                                 const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                 const std::vector<std::string>& aFieldNames,
                                 uint64_t aCellId);
}
}
#endif /* DETCOMMON_XTALK_H */

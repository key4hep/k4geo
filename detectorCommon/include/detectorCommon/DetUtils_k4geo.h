#ifndef DETECTORCOMMON_DETUTILS_H
#define DETECTORCOMMON_DETUTILS_H

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
namespace utils {
dd4hep::xml::Component getNodeByStrAttr(const dd4hep::xml::Handle_t& mother, const std::string& nodeName,
                                        const std::string& attrName, const std::string& attrValue);

/// try to get attribute with double value, return defaultValue if attribute not found
double getAttrValueWithFallback(const dd4hep::xml::Component& node, const std::string& attrName,
                                const double& defaultValue);

/** Retrieves the cellID based on the position of the step and the detector segmentation.
 *  @param aSeg Handle to the segmentation of the volume.
 *  @param aStep Step in which particle deposited the energy.
 *  @param aPreStepPoint Flag indicating if the position of the deposit is the beginning of the step (default)
 *         or the middle of the step.
 */

uint64_t cellID(const dd4hep::Segmentation& aSeg, const G4Step& aStep, bool aPreStepPoint = true);

/** Get number of possible combinations of bit fields for determination of neighbours.
 *   @param[in] aN number of field names.
 *   @param[in] aK length of bit fields included for index search.
 *   return vector of possible combinations of field values.
 */

std::vector<std::vector<uint>> combinations(int N, int K);

/** Get number of possible permutations for certain combination of bit field indeces.
 *   @param[in] aN number of field names.
 *   return vector of permuations for certain field values.
 */

std::vector<std::vector<int>> permutations(int K);

/** Get true field value of neighbour in cyclic bit field
 *   @param[in] aCyclicId field value of neighbour
 *   @param[in] aFieldExtremes Minimal and Maximal values of the fields.
 *   return field value for neighbour in cyclic bit field
 */

int cyclicNeighbour(int aCyclicId, std::pair<int, int> aFieldExtremes);

/**  Get neighbours in many dimensions.
 *   @param[in] aDecoder Handle to the bitfield decoder.
 *   @param[in] aFieldNames Names of the fields for which neighbours are found.
 *   @param[in] aFieldExtremes Minimal and maximal values for the fields.
 *   @param[in] aCellId ID of cell.
 *   @param[in] aDiagonal If diagonal neighbours should be included (all combinations of fields).
 *   return Vector of neighbours.
 */
std::vector<uint64_t> neighbours(const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                 const std::vector<std::string>& aFieldNames,
                                 const std::vector<std::pair<int, int>>& aFieldExtremes,
                                 uint64_t aCellId,
                                 const std::vector<bool>& aFieldCyclic = {false, false, false, false},
                                 bool aDiagonal = true);

/** Special version of the neighbours function for the readout with module and theta merged cells
 *  Compared to the standard version, it needs a reference to the segmentation class to
 *  access the number of merged cells per layer. The other parameters and return value are the same
 *   @param[in] aSeg Reference to the segmentation object.
 *   @param[in] aDecoder Handle to the bitfield decoder.
 *   @param[in] aFieldNames Names of the fields for which neighbours are found.
 *   @param[in] aFieldExtremes Minimal and maximal values for the fields.
 *   @param[in] aCellId ID of cell.
 *   @param[in] aDiagonal If diagonal neighbours should be included (all combinations of fields).
 *   return Vector of neighbours.
 */
std::vector<uint64_t> neighbours_ModuleThetaMerged(const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo& aSeg,
                                 const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                 const std::vector<std::string>& aFieldNames,
                                 const std::vector<std::pair<int, int>>& aFieldExtremes,
                                 uint64_t aCellId,
                                 bool aDiagonal = false);

/** Get minimal and maximal values that can be decoded in the fields of the bitfield.
 *   @param[in] aDecoder Handle to the bitfield decoder.
 *   @param[in] aFieldNames Names of the fields for which extremes are found.
 *   return Vector of pairs (min,max)
 */
std::vector<std::pair<int, int>> bitfieldExtremes(const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                                  const std::vector<std::string>& aFieldNames);

/** Get the half widths of the box envelope (TGeoBBox).
 *   @param[in] aVolumeId The volume ID.
 *   return Half-widths of the volume (x,y,z).
 */
CLHEP::Hep3Vector envelopeDimensions(uint64_t aVolumeId);

/** Get the dimensions of a tube (TGeoConeSeg).
 *   @param[in] aVolumeId The volume ID.
 *   return Dimensions of the tube (rmin, rmax, z(half-length)).
 */
CLHEP::Hep3Vector tubeDimensions(uint64_t aVolumeId);

/** Get the dimensions of a cone (TGeoCone).
 *   @param[in] aVolumeId The volume ID.
 *   return Dimensions of the cone (rmin, rmax, z(half-length)).
 */
CLHEP::Hep3Vector coneDimensions(uint64_t aVolumeId);

/** Get the extrema in pseudorapidity of a tube or cone volume.
 *   @param[in] aVolumeId The volume ID.
 *   return Pseudorapidity extrema (eta_min, eta_max).
 */
std::array<double, 2> tubeEtaExtremes(uint64_t aVolumeId);

/** Get the extrema in pseudorapidity of an envelope.
 *   @param[in] aVolumeId The volume ID.
 *   return Pseudorapidity extrema (eta_min, eta_max).
 */
std::array<double, 2> envelopeEtaExtremes(uint64_t aVolumeId);

/** Get the extrema in pseudorapidity of a volume. First try to match tube or cone, if it fails use an envelope shape.
 *   @param[in] aVolumeId The volume ID.
 *   return Pseudorapidity extrema (eta_min, eta_max).
 */
std::array<double, 2> volumeEtaExtremes(uint64_t aVolumeId);

/** Get the number of cells for the volume and a given Cartesian XY segmentation.
 *   For an example see: Test/TestReconstruction/tests/options/testcellcountingXYZ.py.
 *   @warning No offset in segmentation is currently taken into account.
 *   @param[in] aVolumeId The volume for which the cells are counted.
 *   @param[in] aSeg Handle to the segmentation of the volume.
 *   return Array of the number of cells in (X, Y).
 */
std::array<uint, 2> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::CartesianGridXY& aSeg);

/** Get the number of cells for the volume and a given Cartesian XYZ segmentation.
 *   For an example see: Test/TestReconstruction/tests/options/testcellcountingXYZ.py.
 *   @warning No offset in segmentation is currently taken into account.
 *   @param[in] aVolumeId The volume for which the cells are counted.
 *   @param[in] aSeg Handle to the segmentation of the volume.
 *   return Array of the number of cells in (X, Y, Z).
 */
std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::CartesianGridXYZ& aSeg);

/** Get the number of cells for the volume and a given Phi-Eta / Phi-Theta / Module-Theta segmentation.
 *   It is assumed that the volume has a cylindrical shape (and full azimuthal coverage)
 *   and that it is centred at (0,0,0).
 *   For an example see: Test/TestReconstruction/tests/options/testcellcountingPhiEta.py.
 *   @warning No offset in segmentation is currently taken into account.
 *   @param[in] aVolumeId The volume for which the cells are counted.
 *   @param[in] aSeg Handle to the segmentation of the volume.
 *   return Array of the number of cells in (phi, eta) / (phi, theta) / (module, theta) and the minimum eta / theta ID.
 */
std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo& aSeg);
std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo& aSeg);
std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo& aSeg);

/** Get the number of cells for the volume and a given R-phi segmentation.
 *   It is assumed that the volume has a cylindrical shape - TGeoTube (and full azimuthal coverage)
 *   and that it is centred at (0,0,0).
 *   For an example see: Test/TestReconstruction/tests/options/testcellcountingRPhi.py.
 *   @warning No offset in segmentation is currently taken into account.
 *   @param[in] aVolumeId The volume for which the cells are counted.
 *   @param[in] aSeg Handle to the segmentation of the volume.
 *   return Array of the number of cells in (r, phi).
 */
std::array<uint, 2> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::PolarGridRPhi& aSeg);

/** Get the number of the volumes containing a given name.
 *   For an example see: Test/TestReconstruction/tests/options/testcellcountingXYZ.py.
 *   @param[in] aHighestVolume The top volume in the geometry.
 *   @param[in] aMatchName Name (or its part) of the volume.
 *   return Number of the volumes.
 */
unsigned int countPlacedVolumes(TGeoVolume* aHighestVolume, const std::string& aMatchName);
}
}
#endif /* DETCOMMON_DETUTILS_H */

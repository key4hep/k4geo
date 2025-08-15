#include "detectorSegmentations/GridRPhiEta_k4geo.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  GridRPhiEta_k4geo::GridRPhiEta_k4geo(const std::string& cellEncoding) : FCCSWGridPhiEta_k4geo(cellEncoding) {
    commonSetup();
  }

  GridRPhiEta_k4geo::GridRPhiEta_k4geo(const BitFieldCoder* decoder) : FCCSWGridPhiEta_k4geo(decoder) { commonSetup(); }

  /// Initialization common to all ctors.
  void GridRPhiEta_k4geo::commonSetup() {
    // define type and description
    _type = "GridRPhiEta_k4geo";
    _description = "R-phi-eta segmentation in the global coordinates";

    // register all necessary parameters (additional to those registered in FCCSWGridPhiEta_k4geo)
    registerParameter("grid_size_r", "Cell size in radial distance", m_gridSizeR, 1.,
                      SegmentationParameter::LengthUnit);
    registerParameter("offset_r", "Angular offset in radial distance", m_offsetR, 0., SegmentationParameter::LengthUnit,
                      true);
    registerIdentifier("identifier_r", "Cell ID identifier for R", m_rID, "r");

    m_etaIndex = decoder()->index(fieldNameEta());
    m_phiIndex = decoder()->index(fieldNamePhi());
    m_rIndex = decoder()->index(m_rID);
  }

  /// determine the local based on the cell ID
  Vector3D GridRPhiEta_k4geo::position(const CellID& cID) const {
    return positionFromREtaPhi(r(cID), eta(cID), phi(cID));
  }

  /// determine the cell ID based on the position
  CellID GridRPhiEta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                   const VolumeID& vID) const {
    CellID cID = vID;
    double lRadius = radiusFromXYZ(globalPosition);
    double lEta = etaFromXYZ(globalPosition);
    double lPhi = phiFromXYZ(globalPosition);
    decoder()->set(cID, m_etaIndex, positionToBin(lEta, gridSizeEta(), offsetEta()));
    decoder()->set(cID, m_phiIndex, positionToBin(lPhi, 2 * M_PI / (double)phiBins(), offsetPhi()));
    decoder()->set(cID, m_rIndex, positionToBin(lRadius, m_gridSizeR, m_offsetR));
    return cID;
  }

  /// determine the radial distance R based on the cell ID
  double GridRPhiEta_k4geo::r(const CellID cID) const {
    CellID rValue = decoder()->get(cID, m_rIndex);
    return binToPosition(rValue, m_gridSizeR, m_offsetR);
  }
} // namespace DDSegmentation
} // namespace dd4hep

#include "detectorSegmentations/GridRPhiEta_k4geo.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
GridRPhiEta_k4geo::GridRPhiEta_k4geo(const std::string& cellEncoding) : FCCSWGridPhiEta_k4geo(cellEncoding) {
  // define type and description
  _type = "GridRPhiEta_k4geo";
  _description = "R-phi-eta segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in FCCSWGridPhiEta_k4geo)
  registerParameter("grid_size_r", "Cell size in radial distance", m_gridSizeR, 1., SegmentationParameter::LengthUnit);
  registerParameter("offset_r", "Angular offset in radial distance", m_offsetR, 0., SegmentationParameter::LengthUnit,
                    true);
  registerIdentifier("identifier_r", "Cell ID identifier for R", m_rID, "r");
}

GridRPhiEta_k4geo::GridRPhiEta_k4geo(const BitFieldCoder* decoder) : FCCSWGridPhiEta_k4geo(decoder) {
  // define type and description
  _type = "GridRPhiEta_k4geo";
  _description = "R-phi-eta segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in FCCSWGridPhiEta_k4geo)
  registerParameter("grid_size_r", "Cell size in radial distance", m_gridSizeR, 1., SegmentationParameter::LengthUnit);
  registerParameter("offset_r", "Angular offset in radial distance", m_offsetR, 0., SegmentationParameter::LengthUnit,
                    true);
  registerIdentifier("identifier_r", "Cell ID identifier for R", m_rID, "r");
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
  _decoder->set(cID, m_etaID, positionToBin(lEta, m_gridSizeEta, m_offsetEta));
  _decoder->set(cID, m_phiID, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
  _decoder->set(cID, m_rID, positionToBin(lRadius, m_gridSizeR, m_offsetR));
  return cID;
}

/// determine the radial distance R based on the current cell ID
//double GridRPhiEta_k4geo::r() const {
//  CellID rValue = (*_decoder)[m_rID].value();
//  return binToPosition(rValue, m_gridSizeR, m_offsetR);
//}

/// determine the radial distance R based on the cell ID
double GridRPhiEta_k4geo::r(const CellID& cID) const {
  CellID rValue = _decoder->get(cID, m_rID);
  return binToPosition(rValue, m_gridSizeR, m_offsetR);
}
} // namespace DDSegmentation
} // namespace dd4hep


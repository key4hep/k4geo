#include "detectorSegmentations/GridTheta_k4geo.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
GridTheta_k4geo::GridTheta_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  // define type and description
  _type = "GridTheta_k4geo";
  _description = "Theta segmentation in the global coordinates";

  // register all necessary parameters
  registerParameter("grid_size_theta", "Cell size in Theta", m_gridSizeTheta, 1., SegmentationParameter::LengthUnit);
  registerParameter("offset_theta", "Angular offset in theta", m_offsetTheta, 0., SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_theta", "Cell ID identifier for theta", m_thetaID, "theta");
}

GridTheta_k4geo::GridTheta_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
  // define type and description
  _type = "GridTheta_k4geo";
  _description = "Etheta segmentation in the global coordinates";

  // register all necessary parameters
  registerParameter("grid_size_theta", "Cell size in Theta", m_gridSizeTheta, 1., SegmentationParameter::LengthUnit);
  registerParameter("offset_theta", "Angular offset in theta", m_offsetTheta, 0., SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_theta", "Cell ID identifier for theta", m_thetaID, "theta");
}

/// determine the local based on the cell ID
Vector3D GridTheta_k4geo::position(const CellID& cID) const {
  return positionFromRThetaPhi(1.0, theta(cID), 0.);
}

/// determine the cell ID based on the position
CellID GridTheta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition, const VolumeID& vID) const {
  CellID cID = vID;
  double lTheta = thetaFromXYZ(globalPosition);
  _decoder->set(cID, m_thetaID, positionToBin(lTheta, m_gridSizeTheta, m_offsetTheta));
  return cID;
}

/// determine the pseudorapidity based on the current cell ID
//double GridTheta_k4geo::theta() const {
//  CellID thetaValue = (*_decoder)[m_thetaID].value();
//  return binToPosition(thetaValue, m_gridSizeTheta, m_offsetTheta);
//}

/// determine the polar angle theta based on the cell ID
double GridTheta_k4geo::theta(const CellID& cID) const {
  CellID thetaValue = _decoder->get(cID, m_thetaID);
  return binToPosition(thetaValue, m_gridSizeTheta, m_offsetTheta);
}

} // namespace DDSegmentation
} // namespace dd4hep


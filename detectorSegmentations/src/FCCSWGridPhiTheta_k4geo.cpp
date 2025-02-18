#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
FCCSWGridPhiTheta_k4geo::FCCSWGridPhiTheta_k4geo(const std::string& cellEncoding) : GridTheta_k4geo(cellEncoding) {
  // define type and description
  _type = "FCCSWGridPhiTheta_k4geo";
  _description = "Phi-theta segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in GridTheta_k4geo)
  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
}

FCCSWGridPhiTheta_k4geo::FCCSWGridPhiTheta_k4geo(const BitFieldCoder* decoder) : GridTheta_k4geo(decoder) {
  // define type and description
  _type = "FCCSWGridPhiTheta_k4geo";
  _description = "Phi-theta segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in GridTheta_k4geo)
  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
}

/// determine the local based on the cell ID
Vector3D FCCSWGridPhiTheta_k4geo::position(const CellID& cID) const {
  return positionFromRThetaPhi(1.0, theta(cID), phi(cID));
}

/// determine the cell ID based on the position
CellID FCCSWGridPhiTheta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {
  CellID cID = vID;
  double lTheta = thetaFromXYZ(globalPosition);
  double lPhi = phiFromXYZ(globalPosition);
  _decoder->set(cID, m_thetaID, positionToBin(lTheta, m_gridSizeTheta, m_offsetTheta));
  _decoder->set(cID, m_phiID, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
  return cID;
}

/// determine the azimuthal angle phi based on the current cell ID
//double FCCSWGridPhiTheta_k4geo::phi() const {
//  CellID phiValue = (*_decoder)[m_phiID].value();
//  return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
//}

/// determine the azimuthal angle phi based on the cell ID
double FCCSWGridPhiTheta_k4geo::phi(const CellID& cID) const {
  CellID phiValue = _decoder->get(cID, m_phiID);
  return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
}
}
}

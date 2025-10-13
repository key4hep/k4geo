#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWGridPhiTheta_k4geo::FCCSWGridPhiTheta_k4geo(const std::string& cellEncoding) : GridTheta_k4geo(cellEncoding) {
    commonSetup();
  }

  FCCSWGridPhiTheta_k4geo::FCCSWGridPhiTheta_k4geo(const BitFieldCoder* decoder) : GridTheta_k4geo(decoder) {
    commonSetup();
  }

  /// Initialization common to all ctors.
  void FCCSWGridPhiTheta_k4geo::commonSetup() {
    // define type and description
    _type = "FCCSWGridPhiTheta_k4geo";
    _description = "Phi-theta segmentation in the global coordinates";

    // register all necessary parameters (additional to those registered in GridTheta_k4geo)
    registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
    registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
    registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");

    m_thetaIndex = decoder()->index(fieldNameTheta());
    m_phiIndex = decoder()->index(m_phiID);
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
    decoder()->set(cID, m_thetaIndex, positionToBin(lTheta, gridSizeTheta(), offsetTheta()));
    decoder()->set(cID, m_phiIndex, positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi));
    return cID;
  }

  /// determine the azimuthal angle phi based on the cell ID
  double FCCSWGridPhiTheta_k4geo::phi(const CellID cID) const {
    CellID phiValue = decoder()->get(cID, m_phiIndex);
    return binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);
  }
} // namespace DDSegmentation
} // namespace dd4hep

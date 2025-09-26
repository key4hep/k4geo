#include "detectorSegmentations/GridEta_k4geo.h"

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  GridEta_k4geo::GridEta_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) { commonSetup(); }

  GridEta_k4geo::GridEta_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) { commonSetup(); }

  /// Initialization common to all ctors.
  void GridEta_k4geo::commonSetup() {
    // define type and description
    _type = "GridEta_k4geo";
    _description = "Eta segmentation in the global coordinates";

    // register all necessary parameters
    registerParameter("grid_size_eta", "Cell size in eta", m_gridSizeEta, 1., SegmentationParameter::LengthUnit);
    registerParameter("offset_eta", "Angular offset in eta", m_offsetEta, 0., SegmentationParameter::AngleUnit, true);
    registerIdentifier("identifier_eta", "Cell ID identifier for eta", m_etaID, "eta");

    m_etaIndex = decoder()->index(m_etaID);
  }

  /// determine the local based on the cell ID
  Vector3D GridEta_k4geo::position(const CellID& cID) const { return positionFromREtaPhi(1.0, eta(cID), 0.); }

  /// determine the cell ID based on the position
  CellID GridEta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                               const VolumeID& vID) const {
    CellID cID = vID;
    double lEta = etaFromXYZ(globalPosition);
    decoder()->set(cID, m_etaIndex, positionToBin(lEta, m_gridSizeEta, m_offsetEta));
    return cID;
  }

  /// determine the pseudorapidity based on the cell ID
  double GridEta_k4geo::eta(const CellID cID) const {
    CellID etaValue = decoder()->get(cID, m_etaIndex);
    return binToPosition(etaValue, m_gridSizeEta, m_offsetEta);
  }

} // namespace DDSegmentation
} // namespace dd4hep

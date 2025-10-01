//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "detectorSegmentations/SCEPCal_TimingSegmentation_k4geo.h"
#include "DD4hep/Printout.h"
#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

  SCEPCal_TimingSegmentation_k4geo::SCEPCal_TimingSegmentation_k4geo(const std::string& cellEncoding)
      : Segmentation(cellEncoding) {
    commonSetup();
  }

  SCEPCal_TimingSegmentation_k4geo::SCEPCal_TimingSegmentation_k4geo(const BitFieldCoder* decoder)
      : Segmentation(decoder) {
    commonSetup();
  }

  /// Initialization common to all ctors.
  void SCEPCal_TimingSegmentation_k4geo::commonSetup() {
    _type = "SCEPCal_TimingSegmentation_k4geo";
    _description = "SCEPCal timing layer segmentation";
    registerIdentifier("identifier_system", "Cell ID identifier for System", m_systemId, "system");
    registerIdentifier("identifier_phi", "Cell ID identifier for Phi", m_phiId, "phi");
    registerIdentifier("identifier_theta", "Cell ID identifier for Theta", m_thetaId, "theta");
    registerIdentifier("identifier_gamma", "Cell ID identifier for Gamma", m_gammaId, "gamma");

    m_systemIndex = decoder()->index(m_systemId);
    m_phiIndex = decoder()->index(m_phiId);
    m_thetaIndex = decoder()->index(m_thetaId);
    m_gammaIndex = decoder()->index(m_gammaId);
  }

  Vector3D SCEPCal_TimingSegmentation_k4geo::position(const CellID& cellId) const {
    auto it = m_positionOf.find(cellId);
    if (it != m_positionOf.end()) {
      return it->second;
    }

    dd4hep::printout(dd4hep::WARNING, "SCEPCal_TimingSegmentation_k4geo::position",
                     "cellID %lld not found, returning (0,0,0)!", cellId);
    return Vector3D(0, 0, 0);
  }
} // namespace DDSegmentation
} // namespace dd4hep

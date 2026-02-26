//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "detectorSegmentations/SCEPCal_MainSegmentation_k4geo.h"
#include "DD4hep/Printout.h"
#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

  SCEPCal_MainSegmentation_k4geo::SCEPCal_MainSegmentation_k4geo(const std::string& cellEncoding)
      : Segmentation(cellEncoding) {
    commonSetup();
  }

  SCEPCal_MainSegmentation_k4geo::SCEPCal_MainSegmentation_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) {
    commonSetup();
  }

  /// Initialization common to all ctors.
  void SCEPCal_MainSegmentation_k4geo::commonSetup() {
    _type = "SCEPCal_MainSegmentation_k4geo";
    _description = "SCEPCal main layer segmentation";
    registerIdentifier("identifier_system", "Cell ID identifier for System", m_systemId, "system");
    registerIdentifier("identifier_phi", "Cell ID identifier for Phi", m_phiId, "phi");
    registerIdentifier("identifier_theta", "Cell ID identifier for Theta", m_thetaId, "theta");
    registerIdentifier("identifier_gamma", "Cell ID identifier for Gamma", m_gammaId, "gamma");
    registerIdentifier("identifier_epsilon", "Cell ID identifier for Epsilon", m_epsilonId, "epsilon");
    registerIdentifier("identifier_depth", "Cell ID identifier for Depth", m_depthId, "depth");
    registerIdentifier("identifier_isCheren", "Cell ID identifier for isCherenkov flag", m_isCherenId, "isCheren");

    m_systemIndex = decoder()->index(m_systemId);
    m_phiIndex = decoder()->index(m_phiId);
    m_thetaIndex = decoder()->index(m_thetaId);
    m_gammaIndex = decoder()->index(m_gammaId);
    m_epsilonIndex = decoder()->index(m_epsilonId);
    m_depthIndex = decoder()->index(m_depthId);
    m_isCherenIndex = decoder()->index(m_isCherenId);
  }

  Vector3D SCEPCal_MainSegmentation_k4geo::position(const CellID& cellId) const {
    auto it = m_positionOf.find(cellId);
    if (it != m_positionOf.end()) {
      return it->second;
    }
    dd4hep::printout(dd4hep::WARNING, "SCEPCal_MainSegmentation_k4geo::position",
                     "cellID %lld not found, returning (0,0,0)!", cellId);
    return Vector3D(0, 0, 0);
  }

  // function that retrieves the std::set of neighbouring cell IDs
  // overrides the DDSegmentation::Segmentation::neighbours method
  void SCEPCal_MainSegmentation_k4geo::neighbours(const CellID& cID, std::set<CellID>& neighbours) const {
    // cherenkov and scintillation channel shares the same physical crystal
    neighbours.insert(
        setCellID(System(cID), Phi(cID), Theta(cID), Gamma(cID), Epsilon(cID), Depth(cID), !isCherenkov(cID)));
    // depth is trivial in this segmentation (only two layers)
    int neighborDepth = Depth(cID) == 0 ? 1 : 0;
    neighbours.insert(
        setCellID(System(cID), Phi(cID), Theta(cID), Gamma(cID), Epsilon(cID), neighborDepth, isCherenkov(cID)));

    // then let's do phi direction first (easier one)
    // lambda function to get modulo with wrap-around for negative numbers
    auto modulo = [](int in, int n) -> int { return (in % n + n) % n; };
    // get gamma
    int igamma = Gamma(cID);
    int iphi = Phi(cID);

    if (igamma > 0 && igamma < m_nGamma_ - 1) {
      // middle crystals, just add +/-1 in gamma
      CellID cID_low =
          setCellID(System(cID), Phi(cID), Theta(cID), igamma - 1, Epsilon(cID), Depth(cID), isCherenkov(cID));
      CellID cID_high =
          setCellID(System(cID), Phi(cID), Theta(cID), igamma + 1, Epsilon(cID), Depth(cID), isCherenkov(cID));
      neighbours.insert(cID_low);
      neighbours.insert(cID_high);
    } else if (igamma == 0) {
      // lowest gamma, add +1 in gamma and -1 in phi (with wrap-around)
      CellID cID_high =
          setCellID(System(cID), Phi(cID), Theta(cID), igamma + 1, Epsilon(cID), Depth(cID), isCherenkov(cID));
      int iphi_low = modulo(iphi - 1, m_nPhi_);
      CellID cID_low =
          setCellID(System(cID), iphi_low, Theta(cID), m_nGamma_ - 1, Epsilon(cID), Depth(cID), isCherenkov(cID));
      neighbours.insert(cID_low);
      neighbours.insert(cID_high);
    } else if (igamma == m_nGamma_ - 1) {
      // highest gamma, add -1 in gamma and +1 in phi (with wrap-around)
      CellID cID_low =
          setCellID(System(cID), Phi(cID), Theta(cID), igamma - 1, Epsilon(cID), Depth(cID), isCherenkov(cID));
      int iphi_high = modulo(iphi + 1, m_nPhi_);
      CellID cID_high = setCellID(System(cID), iphi_high, Theta(cID), 0, Epsilon(cID), Depth(cID), isCherenkov(cID));
      neighbours.insert(cID_low);
      neighbours.insert(cID_high);
    }

    // now theta direction
    int itheta = Theta(cID);
    int iDetId = System(cID);
    int differentDetId = (iDetId == m_detId_barrel_) ? m_detId_endcap_ : m_detId_barrel_;
    bool isMiddleOfBarrel = (itheta > m_iTheta_barrel_start_ && itheta < m_iTheta_barrel_end_);
    bool isMiddleOfEndcap = (itheta < m_iTheta_barrel_start_ - 1 || itheta > m_iTheta_barrel_end_ + 1);

    if (isMiddleOfBarrel || isMiddleOfEndcap) {
      // middle crystals, just add +/-1 in theta
      CellID cID_low =
          setCellID(System(cID), Phi(cID), itheta - 1, Gamma(cID), Epsilon(cID), Depth(cID), isCherenkov(cID));
      CellID cID_high =
          setCellID(System(cID), Phi(cID), itheta + 1, Gamma(cID), Epsilon(cID), Depth(cID), isCherenkov(cID));
      neighbours.insert(cID_low);
      neighbours.insert(cID_high);
    } else if (itheta == m_iTheta_barrel_start_ || itheta == m_iTheta_barrel_end_ + 1) {
      // beginning of barrel/endcap, add +1 in theta
      // and -1 in theta to endcap/barrel
      CellID cID_high =
          setCellID(System(cID), Phi(cID), itheta + 1, Gamma(cID), Epsilon(cID), Depth(cID), isCherenkov(cID));
      CellID cID_low =
          setCellID(differentDetId, Phi(cID), itheta - 1, Gamma(cID), Epsilon(cID), Depth(cID), isCherenkov(cID));
      neighbours.insert(cID_low);
      neighbours.insert(cID_high);
    } else if (itheta == m_iTheta_barrel_end_ || itheta == m_iTheta_barrel_start_ - 1) {
      // end of barrel/endcap, add -1 in theta
      // and +1 in theta to endcap/barrel
      CellID cID_low =
          setCellID(System(cID), Phi(cID), itheta - 1, Gamma(cID), Epsilon(cID), Depth(cID), isCherenkov(cID));
      CellID cID_high =
          setCellID(differentDetId, Phi(cID), itheta + 1, Gamma(cID), Epsilon(cID), Depth(cID), isCherenkov(cID));
      neighbours.insert(cID_low);
      neighbours.insert(cID_high);
    }

    return;
  }
} // namespace DDSegmentation
} // namespace dd4hep

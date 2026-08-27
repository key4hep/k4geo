#include "detectorSegmentations/FCCSWGridRhoPhiTheta_k4geo.h"
#include "DD4hep/Detector.h"

#include <algorithm>
#include <cmath>

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWGridRhoPhiTheta_k4geo::FCCSWGridRhoPhiTheta_k4geo(const std::string& cellEncoding) : FCCSWGridPhiTheta_k4geo(cellEncoding) {
    commonSetup();
  }

  FCCSWGridRhoPhiTheta_k4geo::FCCSWGridRhoPhiTheta_k4geo(const BitFieldCoder* decoder) : FCCSWGridPhiTheta_k4geo(decoder) {
    commonSetup();
  }

  /// Initialization common to all ctors.
  void FCCSWGridRhoPhiTheta_k4geo::commonSetup() {
    // define type and description
    _type = "FCCSWGridRhoPhiTheta_k4geo";
    _description = "Rho-Phi-theta segmentation in the global coordinates";
    m_rhoBins.clear();

    // register all necessary parameters (additional to those registered in GridTheta_k4geo)
    registerParameter("grid_rho", "Grid size in rho", m_grid_rho, 0., SegmentationParameter::LengthUnit, true);
    registerParameter("rho_bins", "Non-uniform rho bins", m_rhoBins, std::vector<double>{}, SegmentationParameter::NoUnit, true);
    registerParameter("rho_bin_count", "Number of bins for a uniform rho grid", m_rhoBinCount, 0,
                      SegmentationParameter::NoUnit, true);
    registerParameter("offset_R", "Offset in R", m_offsetR, 0., SegmentationParameter::LengthUnit, true);
    registerParameter("min_theta", "Minimum physical theta", m_minTheta, 0., SegmentationParameter::AngleUnit, true);
    registerParameter("max_theta", "Maximum physical theta", m_maxTheta, M_PI, SegmentationParameter::AngleUnit, true);
    registerParameter("clear_extra_fields", "Clear CellID fields other than system, phi, theta and rho",
                      m_clearExtraFields, false);
    registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");

    m_thetaIndex = decoder()->index(fieldNameTheta());
    m_phiIndex = decoder()->index(fieldNamePhi());
    m_rhoIndex = decoder()->index(m_rhoID);

  }

  void FCCSWGridRhoPhiTheta_k4geo::clearExtraFields(CellID& cID) const {
    if (!m_clearExtraFields)
      return;

    for (const auto& field : decoder()->fields()) {
      const std::string& name = field.name();
      if (name == "system" || name == fieldNamePhi() || name == fieldNameTheta() || name == m_rhoID)
        continue;
      field.set(cID, 0);
    }
  }


  /// determine the local based on the cell ID
  Vector3D FCCSWGridRhoPhiTheta_k4geo::position(const CellID& cID) const {
    return positionFromRThetaPhi(rho(cID)*std::sin(theta(cID)) , theta(cID), phi(cID));
  }


  /// determine the cell ID based on the position
  CellID FCCSWGridRhoPhiTheta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                         const VolumeID& vID) const {
    CellID cID = vID;
    double lTheta = thetaFromXYZ(globalPosition);
    double lPhi = phiFromXYZ(globalPosition);
    double lrho = sqrt( radiusFromXYZ(globalPosition)*radiusFromXYZ(globalPosition) + globalPosition.Z*globalPosition.Z );

    //dd4hep::Segmentation::positionToBin: int(floor((position + 0.5 * cellSize - offset) / cellSize));
    decoder()->set(cID, m_thetaIndex, positionToBin(lTheta, gridSizeTheta(), offsetTheta()));
    decoder()->set(cID, m_phiIndex, positionToBin(lPhi, gridSizePhi(), offsetPhi() ));

    //std::cout<<" Position to cellID "<<std::endl;
    //std::cout<<" Print rho bins: "<<std::endl;
    //for(size_t i=0; i<m_rhoBins.size(); i++) std::cout<<m_rhoBins[i]<<'\t';
    //std::cout<<std::endl;
    //std::cout<<" Global position: ("<<globalPosition.X<<", "<<globalPosition.Y<<", "<<globalPosition.Z<<") "<<std::endl;
    //std::cout<<" Input: rho "<<lrho<<", grid rho: "<<m_grid_rho<<", offsetR "<<m_offsetR<<", theta "<<lTheta<<", offset rho "<<m_offsetR/std::sin(lTheta)<<std::endl;
    //std::cout<<" rho ID:  "<<positionToBin(lrho, m_grid_rho, m_offsetR/std::sin(lTheta)+m_grid_rho/2. )<<std::endl;

    //Two ways for rho segmentation:
    if(m_rhoBins.size()<2)
      // positionToBin(position, cellsize, offset): int(floor((position + 0.5 * cellSize - offset) / cellSize));
      decoder()->set(cID, m_rhoIndex, positionToBin(lrho, m_grid_rho, m_offsetR/std::sin(lTheta)+m_grid_rho/2. ));
    else{
      // Note: a bug was found in dd4hep::Segmentation, and fixed in PR: https://github.com/AIDASoft/DD4hep/pull/1625
      // Following function can be used after updating dd4hep to latest ver.
      // int rhoID = Segmentation::positionToBin(lrho, m_rhoBins, m_offsetR/std::sin(lTheta));  //Fixed in DD4hep main branch, no tag yet.

      std::vector<double>::const_iterator bin = std::upper_bound(m_rhoBins.begin(),
                                                                 m_rhoBins.end(),
                                                                 lrho-m_offsetR/std::sin(lTheta));
      int rhoID = bin - m_rhoBins.begin() - 1 ;

      //Add a protection for cell ID, if the position is out of bin range
      //WARNING: may introduce bug. Better to well-define bin edges.
      if(rhoID<0) rhoID = 0;
      if(rhoID>m_rhoBins.size()-2) rhoID = m_rhoBins.size()-2;
      decoder()->set(cID, m_rhoIndex, rhoID);
    }

    clearExtraFields(cID);
    return cID;
  }


  /// determine the distance to IP rho based on the cell ID
  double FCCSWGridRhoPhiTheta_k4geo::rho(const CellID cID) const {
    CellID rhoValue = decoder()->get(cID, m_rhoIndex);
    double m_theta = theta(cID);

    // dd4hep::Segmentation::binToPosition: bin * cellSize + offset;

    //std::cout<<" CellID to position"<<std::endl;
    //std::cout<<" Print rho bins: "<<std::endl;
    //for(size_t i=0; i<m_rhoBins.size(); i++) std::cout<<m_rhoBins[i]<<'\t';
    //std::cout<<std::endl;
    //std::cout<<" rho ID val: "<<rhoValue<<", theta val: "<<m_theta<<", offsetR "<<m_offsetR<<", rho = "<<0.5*(m_rhoBins[rhoValue] + m_rhoBins[rhoValue+1])<<std::endl;

    if(m_rhoBins.size()<2)
      return binToPosition(rhoValue, m_grid_rho, m_offsetR/std::sin(m_theta)+m_grid_rho/2. );
    else
      return m_offsetR/std::sin(m_theta) + 0.5*(m_rhoBins[rhoValue] + m_rhoBins[rhoValue+1]);

  }

  // overrides the DDSegmentation::Segmentation::neighbours method
  void FCCSWGridRhoPhiTheta_k4geo::neighbours(const CellID& cID, std::set<CellID>& neighbours) const {
    const int phiBin = static_cast<int>(decoder()->get(cID, m_phiIndex));
    const int thetaBin = static_cast<int>(decoder()->get(cID, m_thetaIndex));
    const int rhoBin = static_cast<int>(decoder()->get(cID, m_rhoIndex));

    auto withBin = [&](int fieldIndex, int bin) {
      CellID neighbour = cID;
      decoder()->set(neighbour, fieldIndex, bin);
      clearExtraFields(neighbour);
      return neighbour;
    };

    // phiFromXYZ() uses atan2(), hence phi bins are generally signed.  In
    // particular, they are not [0, phiBins()-1] when offset_phi is zero.
    if (phiBins() > 1) {
      const int firstPhiBin = positionToBin(-M_PI, gridSizePhi(), offsetPhi());
      const auto wrapPhi = [&](int bin) {
        const int relativeBin = bin - firstPhiBin;
        return firstPhiBin + (relativeBin % phiBins() + phiBins()) % phiBins();
      };
      neighbours.insert(withBin(m_phiIndex, wrapPhi(phiBin - 1)));
      neighbours.insert(withBin(m_phiIndex, wrapPhi(phiBin + 1)));
    }

    // Theta is not periodic.  Use the bins reached by points at the detector
    // boundaries rather than testing bin centres: the first/last physical
    // cells can have centres just outside the boundary.
    const auto& thetaField = (*decoder())[m_thetaIndex];
    const int firstThetaBin = std::max(positionToBin(m_minTheta, gridSizeTheta(), offsetTheta()),
                                       static_cast<int>(thetaField.minValue()));
    const int lastThetaBin = std::min(positionToBin(m_maxTheta, gridSizeTheta(), offsetTheta()),
                                      static_cast<int>(thetaField.maxValue()));
    for (const int candidate : {thetaBin - 1, thetaBin + 1}) {
      if (candidate >= firstThetaBin && candidate <= lastThetaBin)
        neighbours.insert(withBin(m_thetaIndex, candidate));
    }

    // Non-uniform rho bins define their count directly.  A uniform grid needs
    // rho_bin_count from the compact file; the encoded field width is only a
    // storage limit and is not normally the number of physical layers.
    const auto& rhoField = (*decoder())[m_rhoIndex];
    const int minRhoBin = std::max(0, static_cast<int>(rhoField.minValue()));
    int maxRhoBin = static_cast<int>(rhoField.maxValue());
    if (m_rhoBins.size() >= 2)
      maxRhoBin = std::min(static_cast<int>(m_rhoBins.size()) - 2, maxRhoBin);
    else if (m_rhoBinCount > 0)
      maxRhoBin = std::min(m_rhoBinCount - 1, maxRhoBin);
    for (const int candidate : {rhoBin - 1, rhoBin + 1}) {
      if (candidate >= minRhoBin && candidate <= maxRhoBin)
        neighbours.insert(withBin(m_rhoIndex, candidate));
    }
  }


} // namespace DDSegmentation
} // namespace dd4hep

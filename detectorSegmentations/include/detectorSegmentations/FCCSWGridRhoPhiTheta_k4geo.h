#ifndef DETECTORSEGMENTATIONS_GRIDRHOPHITHETA_K4GEO_H
#define DETECTORSEGMENTATIONS_GRIDRHOPHITHETA_K4GEO_H

// FCCSW
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"

/** FCCSWGridRhoPhiTheta_k4geo Detector/detectorSegmentations/detectorSegmentations/FCCSWGridRhoPhiTheta_k4geo.h
 * FCCSWGridRhoPhiTheta_k4geo.h
 *
 *  Segmentation in rho, theta and phi.
 *  Based on FCCSWGridPhiTheta_k4geo, addition of Rho segmentation
 *
 */

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWGridRhoPhiTheta_k4geo : public FCCSWGridPhiTheta_k4geo {
  public:
    /// default constructor using an arbitrary type
    FCCSWGridRhoPhiTheta_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    FCCSWGridRhoPhiTheta_k4geo(const BitFieldCoder* decoder);

    /// destructor
    virtual ~FCCSWGridRhoPhiTheta_k4geo() = default;

    /**  Determine the global position based on the cell ID.
     *   @warning This segmentation has no knowledge of radius, so radius = 1 is taken into calculations.
     *   @param[in] aCellId ID of a cell.
     *   return Position (radius = 1).
     */
    virtual Vector3D position(const CellID& aCellID) const override;
    /**  Determine the cell ID based on the position.
     *   @param[in] aLocalPosition (not used).
     *   @param[in] aGlobalPosition position in the global coordinates.
     *   @param[in] aVolumeId ID of a volume.
     *   return Cell ID.
     */
    virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                          const VolumeID& aVolumeID) const override;

    /** Add the face-sharing neighbours in rho, phi and theta.
     *  Phi neighbours wrap around at 2 pi; rho and theta do not.
     */
    virtual void neighbours(const CellID& cellID, std::set<CellID>& neighbours) const override;


    /**  Find neighbours of the cell.
     *   Definition of neighbours is explained on slide 9:
     * https://indico.cern.ch/event/1475808/contributions/6219554/attachments/2966253/5218774/FCC_FullSim_HCal_slides.pdf
     *   @param[in] cID ID of a cell.
     *   @param[in] aDiagonal if true, will include neighbours from diagonal positions in the next and previous layers.
     *   return vector of neighbour cellIDs.
     */
    //std::vector<uint64_t> neighbours(const CellID cID, bool aDiagonal) const;

    /**  Find neighbours of the cell.
     *   Implement the signature from the Segmentation base class.
     */
    //virtual void neighbours(const CellID& cellID, std::set<CellID>& neighbours) const override;

  protected:
    /// Get rho from cellID
    double rho(const CellID cID) const;

  private:
    /// Clear all CellID fields except system, phi, theta and rho when enabled.
    void clearExtraFields(CellID& cellID) const;

    /// the grid size in rho
    double m_grid_rho;
    std::vector<double> m_rhoBins;
    /// number of bins for a uniform rho grid (0 keeps the legacy field-range fallback)
    int m_rhoBinCount;

    /// physical theta coverage of the segmented detector
    double m_minTheta;
    double m_maxTheta;

    /// whether placement and other non-grid CellID fields are cleared
    bool m_clearExtraFields;

    /// the coordinate offset in rho / R
    double m_offsetR;   // In cylinder case: offset depends on R not rho.
    double m_offsetrho;

    /// the field name used for rho
    std::string m_rhoID;

    /// Initialization common to all ctors.
    void commonSetup();
    /// the field index used for theta
    int m_thetaIndex = -1;
    /// the field index used for phi
    int m_phiIndex = -1;
    /// the field index used for rho
    int m_rhoIndex = -1;



  };
} // namespace DDSegmentation
} // namespace dd4hep
#endif /* DETSEGMENTATION_GRIDRHOPHITHETA_H */

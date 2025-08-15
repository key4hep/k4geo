#ifndef DETECTORSEGMENTATIONS_GRIDPHITHETA_K4GEO_H
#define DETECTORSEGMENTATIONS_GRIDPHITHETA_K4GEO_H

// FCCSW
#include "detectorSegmentations/GridTheta_k4geo.h"

/** FCCSWGridPhiTheta_k4geo Detector/detectorSegmentations/detectorSegmentations/FCCSWGridPhiTheta_k4geo.h
 * FCCSWGridPhiTheta_k4geo.h
 *
 *  Segmentation in theta and phi.
 *  Based on GridTheta_k4geo, addition of azimuthal angle coordinate.
 *
 */

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWGridPhiTheta_k4geo : public GridTheta_k4geo {
  public:
    /// default constructor using an arbitrary type
    FCCSWGridPhiTheta_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    FCCSWGridPhiTheta_k4geo(const BitFieldCoder* decoder);

    /// destructor
    virtual ~FCCSWGridPhiTheta_k4geo() = default;

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
    /**  Determine the azimuthal angle based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   return Phi.
     */
    double phi(const CellID aCellID) const;
    /**  Get the grid size in phi.
     *   return Grid size in phi.
     */
    inline double gridSizePhi() const { return 2 * M_PI / static_cast<double>(m_phiBins); }
    /**  Get the number of bins in azimuthal angle.
     *   return Number of bins in phi.
     */
    inline int phiBins() const { return m_phiBins; }
    /**  Get the coordinate offset in azimuthal angle.
     *   return The offset in phi.
     */
    inline double offsetPhi() const { return m_offsetPhi; }
    /**  Get the field name for azimuthal angle.
     *   return The field name for phi.
     */
    inline const std::string& fieldNamePhi() const { return m_phiID; }
    /**  Set the number of bins in azimuthal angle.
     *   @param[in] aNumberBins Number of bins in phi.
     */
    inline void setPhiBins(int bins) { m_phiBins = bins; }
    /**  Set the coordinate offset in azimuthal angle.
     *   @param[in] aOffset Offset in phi.
     */
    inline void setOffsetPhi(double offset) { m_offsetPhi = offset; }
    /**  Set the field name used for azimuthal angle.
     *   @param[in] aFieldName Field name for phi.
     */
    inline void setFieldNamePhi(const std::string& fieldName) { m_phiID = fieldName; }
    /** Returns a std::vector<double> of the cellDimensions of the given cell ID
     *  in natural order of dimensions (dPhi, dTheta)
     *  @param[in] cellID
     *  return a std::vector of size 2 with the cellDimensions of the given cell ID (phi, theta)
     */
    virtual std::vector<double> cellDimensions(const CellID& /* id */) const override {
      return {gridSizePhi(), gridSizeTheta()};
    }

  private:
    /// determine the azimuthal angle phi based on the current cell ID
    double phi() const;
    /// the number of bins in phi
    int m_phiBins;
    /// the coordinate offset in phi
    double m_offsetPhi;
    /// the field name used for phi
    std::string m_phiID;

    /// Initialization common to all ctors.
    void commonSetup();
    /// the field index used for theta
    int m_thetaIndex = -1;
    /// the field index used for phi
    int m_phiIndex = -1;
  };
} // namespace DDSegmentation
} // namespace dd4hep
#endif /* DETSEGMENTATION_GRIDPHITHETA_H */

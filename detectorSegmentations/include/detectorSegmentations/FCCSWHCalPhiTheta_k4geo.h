#ifndef DETECTORSEGMENTATIONS_HCALPHITHETA_K4GEO_H
#define DETECTORSEGMENTATIONS_HCALPHITHETA_K4GEO_H

// FCCSW
#include "detectorSegmentations/GridTheta_k4geo.h"

/** FCCSWHCalPhiTheta_k4geo Detector/detectorSegmentations/detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h
 * FCCSWHCalPhiTheta_k4geo.h
 *
 *  Segmentation in theta and phi.
 *  Based on GridTheta_k4geo, addition of azimuthal angle coordinate.
 *
 *  Rectangular shape cells are defined in r-z plan and all the hits within a defined cell boundaries are assigned the
 * same cellID. Cell borders are defined such that closely follow the theta projective towers. More details:
 * https://indico.cern.ch/event/1475808/contributions/6219554/attachments/2966253/5218774/FCC_FullSim_HCal_slides.pdf
 *
 */

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWHCalPhiTheta_k4geo : public GridTheta_k4geo {
  public:
    /// default constructor using an arbitrary type
    FCCSWHCalPhiTheta_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    FCCSWHCalPhiTheta_k4geo(const BitFieldCoder* decoder);

    /// destructor
    virtual ~FCCSWHCalPhiTheta_k4geo() = default;

    /**  Get the postion of the geometric center of the cell based on the cellID
     *   @param[in] aCellID
     *   return the global coordinates of cell center
     */
    virtual Vector3D position(const CellID& aCellID) const;

    /**  Assign a cellID based on the global position.
     *  @param[in] aLocalPosition (not used).
     *  @param[in] aGlobalPosition position in the global coordinates.
     *  @param[in] aVolumeId ID of a volume.
     *  return Cell ID.
     */
    virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                          const VolumeID& aVolumeID) const;

    /**  Find neighbours of the cell.
     *   Definition of neighbours is explained on slide 9:
     * https://indico.cern.ch/event/1475808/contributions/6219554/attachments/2966253/5218774/FCC_FullSim_HCal_slides.pdf
     *   @param[in] aCellId ID of a cell.
     *   @param[in] aDiagonal if true, will include neighbours from diagonal positions in the next and previous layers.
     *   return vector of neighbour cellIDs.
     */
    std::vector<uint64_t> neighbours(const CellID& cID, bool aDiagonal) const;

    /**  Calculate layer radii and edges in z-axis, then define cell edges in each layer using defineCellEdges().
     *    Following member variables are calculated:
     *      m_radii
     *      m_layerEdges
     *      m_layerDepth
     *      m_thetaBins (updated through defineCellEdges())
     *      m_cellEdges (updated through defineCellEdges())
     */
    void defineCellsInRZplan() const;

    /**  Define cell edges in z-axis for the given layer.
     *   Logic:
     *      1) Find theta bin centers that fit within the given layer;
     *      2) Define a cell edge in z-axis as the middle of each pair of theta bin centers
     *   @param[in] layer index
     */
    void defineCellEdges(const unsigned int layer) const;

    /**  Determine the azimuthal angle of HCal cell based on the cellID.
     *   @param[in] aCellId ID of a cell.
     *   return Phi.
     */
    double phi(const CellID& aCellID) const;

    /**  Get the grid size in phi.
     *   return Grid size in phi.
     */
    inline double gridSizePhi() const { return 2 * M_PI / static_cast<double>(m_phiBins); }

    /**  Get the number of bins in azimuthal angle.
     *   return Number of bins in phi.
     */
    inline int phiBins() const { return m_phiBins; }

    /**  Get the coordinate offset in azimuthal angle, which is the center of cell with m_phiID=0
     *   return The offset in phi.
     */
    inline double offsetPhi() const { return m_offsetPhi; }

    /**  Get the vector of theta bins (cells) in a given layer.
     */
    inline std::vector<int> thetaBins(const uint layer) const {
      if (m_radii.empty())
        defineCellsInRZplan();
      if (!m_thetaBins.empty())
        return m_thetaBins[layer];
      else
        return std::vector<int>();
    }

    /**  Get the coordinate offset in z-axis.
     *   Offset is the middle position of the Barrel or each section of the Endcap.
     *   For the Barrel, the vector size is 1, while for the Endcap - number of section.
     *   return The offset in z.
     */
    inline std::vector<double> offsetZ() const { return m_offsetZ; }

    /**  Get the z length of the layer.
     *   return the z length.
     */
    inline std::vector<double> widthZ() const { return m_widthZ; }

    /**  Get the coordinate offset in radius.
     *   Offset is the inner radius of the first layer in the Barrel or in each section of the Endcap.
     *   For the Barrel, the vector size is 1, while for the Endcap - number of sections.
     *   return the offset in radius.
     */
    inline std::vector<double> offsetR() const { return m_offsetR; }

    /**  Get the number of layers for each different thickness retrieved with dRlayer().
     *   For the Barrel, the vector size equals to the number of different thicknesses used to form the layers.
     *   For the Endcap, the vector size equals to the number of sections in the Endcap times the number of different
     * thicknesses used to form the layers. return the number of layers.
     */
    inline std::vector<int> numLayers() const { return m_numLayers; }

    /**  Get the dR (thickness) of layers.
     *   The size of the vector equals to the number of different thicknesses used to form the layers.
     *   return the dR.
     */
    inline std::vector<double> dRlayer() const { return m_dRlayer; }

    /**  Get the field name for azimuthal angle.
     *   return The field name for phi.
     */
    inline const std::string& fieldNamePhi() const { return m_phiID; }

    /**  Get the field name for layer.
     *   return The field name for layer.
     */
    inline const std::string& fieldNameLayer() const { return m_layerID; }

    /**  Determine the minimum and maximum polar angle of HCal cell based on the cellID.
     *   @param[in] aCellId ID of a cell.
     *   return Theta.
     */
    std::array<double, 2> cellTheta(const CellID& cID) const;

    /**  Get the min and max layer indexes of each HCal part.
     * For Endcap, returns the three elements vector, while for Barrel - single element vector.
     */
    std::vector<std::pair<uint, uint>> getMinMaxLayerId() const;

    /**  Set the number of bins in azimuthal angle.
     *   @param[in] aNumberBins Number of bins in phi.
     */
    inline void setPhiBins(int bins) { m_phiBins = bins; }

    /**  Set the coordinate offset in azimuthal angle.
     *   @param[in] aOffset Offset in phi.
     */
    inline void setOffsetPhi(double offset) { m_offsetPhi = offset; }

    /**  Set the detector layout (0 = Barrel; 1 = Endcap).
     *   @param[in] detLayout
     */
    inline void setDetLayout(int detLayout) { m_detLayout = detLayout; }

    /**  Set the coordinate offset in z-axis.
     *   @param[in] offset in z (centre of the layer).
     */
    inline void setOffsetZ(std::vector<double> const& offset) { m_offsetZ = offset; }

    /**  Set the z width.
     *   @param[in] width in z.
     */
    inline void setWidthZ(std::vector<double> const& width) { m_widthZ = width; }

    /**  Set the coordinate offset in radius.
     *   @param[in] offset in radius.
     */
    inline void setOffsetR(std::vector<double> const& offset) { m_offsetR = offset; }

    /**  Set the number of layers.
     *   @param[in] number of layers
     */
    inline void setNumLayers(std::vector<int> const& num) { m_numLayers = num; }

    /**  Set the dR of layers.
     *   @param[in] dR of layers.
     */
    inline void setdRlayer(std::vector<double> const& dRlayer) { m_dRlayer = dRlayer; }

    /**  Set the field name used for azimuthal angle.
     *   @param[in] aFieldName Field name for phi.
     */
    inline void setFieldNamePhi(const std::string& fieldName) { m_phiID = fieldName; }

  protected:
    /// determine the azimuthal angle phi based on the current cell ID
    double phi() const;
    /// the number of bins in phi
    int m_phiBins;
    /// the coordinate offset in phi
    double m_offsetPhi;
    /// the field name used for phi
    std::string m_phiID;
    /// the field name used for layer
    std::string m_layerID;
    /// the detector layout (0 = Barrel; 1 = Endcap)
    int m_detLayout;
    /// the z offset of middle of the layer
    std::vector<double> m_offsetZ;
    /// the z width of the layer
    std::vector<double> m_widthZ;
    /// the radius of the layer
    std::vector<double> m_offsetR;
    /// number of layers
    std::vector<int> m_numLayers;
    /// dR of the layer
    std::vector<double> m_dRlayer;
    /// radius of each layer
    mutable std::vector<double> m_radii;
    /// z-min and z-max of each layer
    mutable std::vector<std::pair<double, double>> m_layerEdges;
    /// dR of each layer
    mutable std::vector<double> m_layerDepth;
    /// theta bins (cells) in each layer
    mutable std::vector<std::vector<int>> m_thetaBins;
    /// z-min and z-max of each cell (theta bin) in each layer
    mutable std::vector<std::unordered_map<int, std::pair<double, double>>> m_cellEdges;
  };
} // namespace DDSegmentation
} // namespace dd4hep
#endif /* DETSEGMENTATION_HCALPHITHETA_H */

#ifndef DETECTORSEGMENTATIONS_HCALPHIROW_K4GEO_H
#define DETECTORSEGMENTATIONS_HCALPHIROW_K4GEO_H

#include "DDSegmentation/Segmentation.h"
#include "DDSegmentation/SegmentationUtil.h"

#include <array>
#include <atomic>
#include <string>
#include <unordered_map>
#include <vector>

/** FCCSWHCalPhiRow_k4geo Detector/detectorSegmentations/detectorSegmentations/FCCSWHCalPhiRow_k4geo.h
 * FCCSWHCalPhiRow_k4geo.h
 *
 *  Segmentation in row and phi.
 *
 *  Cells are defined in r-z plan by merging the row/sequence of scintillator/absorber.
 *  Each row consists of 2 * Master_Plate + 1 * Spacer_Plate + 1 * Scintillator.
 *  Considering a single row as a single cell gives the highest possible granularity.
 *  More details:
 * https://indico.cern.ch/event/1475808/contributions/6219554/attachments/2966253/5218774/FCC_FullSim_HCal_slides.pdf
 *
 */

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWHCalPhiRow_k4geo : public Segmentation {
  public:
    /// default constructor using an arbitrary type
    FCCSWHCalPhiRow_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    FCCSWHCalPhiRow_k4geo(const BitFieldCoder* decoder);

    /**  Determine the position of HCal cell based on the cellID.
     *   @param[in] aCellID ID of a cell.
     *   return Position (radius = 1).
     */
    virtual Vector3D position(const CellID& aCellID) const override;

    /**  Assign a cellID based on the global position.
     *  @param[in] aLocalPosition (not used).
     *  @param[in] aGlobalPosition position in the global coordinates.
     *  @param[in] aVolumeID ID of a volume.
     *  return Cell ID.
     */
    virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                          const VolumeID& aVolumeID) const override;

    /**  Find neighbours of the cell
     *   Definition of neighbours is explained on slide 7:
     * https://indico.cern.ch/event/1475808/contributions/6219554/attachments/2966253/5218774/FCC_FullSim_HCal_slides.pdf
     *   @param[in] cID ID of a cell.
     *   return vector of neighbour cellIDs.
     */
    std::vector<uint64_t> neighbours(const CellID cID) const;

    /**  Find neighbours of the cell.
     *   Implement the signature from the Segmentation base class.
     */
    virtual void neighbours(const CellID& cellID, std::set<CellID>& neighbours) const override;

    /**  Calculate layer radii and edges in z-axis.
     *    Following member variables are calculated:
     *      m_radii
     *      m_layerEdges
     *      m_layerDepth
     */
    void calculateLayerRadii() const;

    /**  Define cell indexes for the given layer.
     *  This function fills the m_cellIndexes vector per layer with the cell indexes.
     *  The cell index is encoded in CellID with "row" field.
     *  In case of a cell with single row/sequence, the index is directly the number of row in the layer.
     *  In case of a cell with several rows/sequences merged, the index is the number of cell in the layer.
     *  For the layers of negative-z Endcap, indexes of cells are negative.
     *  Following member variables are calculated:
     *   m_cellIndexes
     *   m_cellEdges
     *   @param[in] layer index
     */
    void defineCellIndexes(const unsigned int layer) const;

    /**  Determine the azimuthal angle of HCal cell based on the cellID.
     *   @param[in] aCellID ID of a cell.
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

    /**  Get the coordinate offset in azimuthal angle, which is the center of cell with m_phiID=0.
     *   return The offset in phi.
     */
    inline double offsetPhi() const { return m_offsetPhi; }

    /**  Get the grid size in row for each layer.
     *   return Grid size in row.
     */
    inline const std::vector<int>& gridSizeRow() const { return m_gridSizeRow; }

    /**  Determine the polar angle of HCal cell center based on the cellID.
     *   @param[in] aCellID ID of a cell.
     *   return Theta.
     */
    inline double theta(const CellID aCellID) const {
      return dd4hep::DDSegmentation::Util::thetaFromXYZ(position(aCellID));
    }

    /**  Determine the minimum and maximum polar angle of HCal cell based on the cellID.
     *   @param[in] cID ID of a cell.
     *   return Theta.
     */
    std::array<double, 2> cellTheta(const CellID cID) const;

    /**  Get the vector of cell indexes in a given layer.
     */
    inline std::vector<int> cellIndexes(const uint layer) const {
      const LayerInfo& li = getLayerInfo(layer);
      return li.cellIndexes;
    }

    /**  Get the thetaMax needed for SW clustering
     *   return max theta value of the detector
     */
    double thetaMax() const;

    /**  Get the min and max layer indexes of each HCal part.
     * For Endcap, returns the three elements vector, while for Barrel - single element vector.
     */
    std::vector<std::pair<uint, uint>> getMinMaxLayerId() const;

    /**  Get the coordinate offset in z-axis.
     *   Offset is the middle position of the Barrel or each section of the Endcap.
     *   For the Barrel, the vector size is 1, while for the Endcap - number of section.
     *   return The offset in z.
     */
    inline const std::vector<double>& offsetZ() const { return m_offsetZ; }

    /**  Get the z length of the layer.
     *   return the z length.
     */
    inline const std::vector<double>& widthZ() const { return m_widthZ; }

    /**  Get the coordinate offset in radius.
     *   Offset is the inner radius of the first layer in the Barrel or in each section of the Endcap.
     *   For the Barrel, the vector size is 1, while for the Endcap - number of sections.
     *   return the offset in radius.
     */
    inline const std::vector<double>& offsetR() const { return m_offsetR; }

    /**  Get the number of layers for each different thickness retrieved with dRlayer().
     *   For the Barrel, the vector size equals to the number of different thicknesses used to form the layers.
     *   For the Endcap, the vector size equals to the number of sections in the Endcap times the number of different
     * thicknesses used to form the layers. return the number of layers.
     */
    inline const std::vector<int>& numLayers() const { return m_numLayers; }

    /**  Get the dR (thickness) of layers.
     *   The size of the vector equals to the number of different thicknesses used to form the layers.
     *   return the dR.
     */
    inline const std::vector<double>& dRlayer() const { return m_dRlayer; }

    /**  Get the field name for azimuthal angle.
     *   return The field name for phi.
     */
    inline const std::string& fieldNamePhi() const { return m_phiID; }

    /**  Get the field name for layer.
     *   return The field name for layer.
     */
    inline const std::string& fieldNameLayer() const { return m_layerID; }

    /**  Get the field name for row number.
     *   return The field name for row.
     */
    inline const std::string& fieldNameRow() const { return m_rowID; }

    /**  Get the layer for a cell given cell ID
     *   @param[in] aCellID the cell ID
     *   return The layer number
     */
    inline int layer(const CellID aCellID) const { return _decoder->get(aCellID, fieldNameLayer()); }

    /**  Set the number of bins in azimuthal angle.
     *   @param[in] bins Number of bins in phi.
     */
    inline void setPhiBins(int bins) { m_phiBins = bins; }

    /**  Set the coordinate offset in azimuthal angle.
     *   @param[in] offset Offset in phi.
     */
    inline void setOffsetPhi(double offset) { m_offsetPhi = offset; }

    /**  Set the grid size in theta angle.
     *   @param[in] size Cell size in theta.
     */
    inline void setGridSizeRow(std::vector<int> const& size) { m_gridSizeRow = size; }

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
     *   @param[in] num number of layers
     */
    inline void setNumLayers(std::vector<int> const& num) { m_numLayers = num; }

    /**  Set the dR of layers.
     *   @param[in] dRlayer dR of layers.
     */
    inline void setdRlayer(std::vector<double> const& dRlayer) { m_dRlayer = dRlayer; }

    /**  Set the field name used for azimuthal angle.
     *   @param[in] fieldName Field name for phi.
     */
    inline void setFieldNamePhi(const std::string& fieldName) { m_phiID = fieldName; }

    /**  Set the field name used for row number.
     *   @param[in] fieldName Field name for row.
     */
    inline void setFieldNameRow(const std::string& fieldName) { m_rowID = fieldName; }

    /** Returns a std::vector<double> of the cellDimensions of the given cell ID
     *  in natural order of dimensions (phi, z)
     *  @param[in] id
     *  return a std::vector of size 2 with the cellDimensions of the given cell ID (phi, z)
     */
    virtual std::vector<double> cellDimensions(const CellID& id) const override {
      const int aLayer = layer(id);
      return {gridSizePhi(), m_gridSizeRow[aLayer] * m_dz_row};
    }

  private:
    /// the number of bins in phi
    int m_phiBins;
    /// the coordinate offset in phi
    double m_offsetPhi;
    /// the field name used for phi
    std::string m_phiID;
    /// the field name used for layer
    std::string m_layerID;
    /// the grid size in row for each layer
    std::vector<int> m_gridSizeRow;
    /// dz of row
    double m_dz_row;
    /// the field name used for row
    std::string m_rowID;
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

    /// Initialization common to all ctors.
    void commonSetup();
    /// the field index used for layer
    int m_layerIndex = -1;
    /// the field index used for row
    int m_rowIndex = -1;
    /// the field index used for type
    int m_typeIndex = -1;
    /// the field index used for phi
    int m_phiIndex = -1;

    // Derived geometrical information about each layer.
    struct LayerInfo {
      /// Radius of the layer.
      double radius = 1;

      /// Half the layer depth (dR).
      double halfDepth = 0;

      /// z-min and z-max of the layer
      double zmin = 0;
      double zmax = 0;

      /// cell indexes in each layer
      std::vector<int> cellIndexes{};

      /// z-min and z-max of each cell in each layer
      std::unordered_map<int, std::pair<double, double>> cellEdges{};
    };

    // The vector of tabulated values, indexed by layer number.
    // We can't build this in the constructor --- the volumes won't have
    // been created yet.  Instead, build it lazily the first time it's needed.
    // Since that's in a const method, make it thread-safe.
    mutable std::atomic<const std::vector<LayerInfo>*> m_layerInfo = nullptr;

    // Retrieve the derived geometrical information for a given layer.
    const LayerInfo& getLayerInfo(const unsigned layer) const;

    /**  Construct the derived geometrical information.
     * Calculate layer radii and edges in z-axis, then define cell edges in each layer using defineCellEdges().
     *    Following member variables are calculated:
     *      radius
     *      layerEdges
     *      layerDepth
     *      thetaBins (updated through defineCellEdges())
     *      cellEdges* (updated through defineCellEdges())
     */
    std::vector<LayerInfo> initLayerInfo() const;

    /**  Define cell indexes for the given layer.
     *  This function fills the cellIndexes vector per layer with the cell indexes.
     *  The cell index is encoded in CellID with "row" field.
     *  In case of a cell with single row/sequence, the index is directly the number of row in the layer.
     *  In case of a cell with several rows/sequences merged, the index is the number of cell in the layer.
     *  For the layers of negative-z Endcap, indexes of cells are negative.
     *  Following member variables are calculated:
     *   cellIndexes
     *   cellEdges
     *   @param[in] li Layer info entry corresponding to layer.
     *   @param[in] layer index
     */
    void defineCellIndexes(LayerInfo& li, const unsigned int layer) const;

    // Check consistency of input geometric variables.
    bool checkParameters() const;
  };
} // namespace DDSegmentation
} // namespace dd4hep
#endif /* DETSEGMENTATION_HCALPHITHETA_H */

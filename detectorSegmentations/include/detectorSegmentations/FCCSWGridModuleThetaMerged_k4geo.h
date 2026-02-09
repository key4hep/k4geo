#ifndef DETECTORSEGMENTATIONS_FCCSWGridModuleThetaMerged_k4geo_H
#define DETECTORSEGMENTATIONS_FCCSWGridModuleThetaMerged_k4geo_H

// FCCSW
#include "detectorSegmentations/GridTheta_k4geo.h"
#include <atomic>

/** FCCSWGridModuleThetaMerged_k4geo Detector/DetSegmentation/DetSegmentation/FCCSWGridModuleThetaMerged_k4geo.h
 * FCCSWGridModuleThetaMerged_k4geo.h
 *
 *  Segmentation in theta and module.
 *  Based on GridTheta, merges modules and theta cells based on layer number
 *
 */

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWGridModuleThetaMerged_k4geo : public GridTheta_k4geo {
  public:
    /// default constructor using an arbitrary type
    FCCSWGridModuleThetaMerged_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    FCCSWGridModuleThetaMerged_k4geo(const BitFieldCoder* decoder);

    /// destructor
    virtual ~FCCSWGridModuleThetaMerged_k4geo();

    /// read n(modules) from detector metadata
    void GetNModulesFromGeom();

    /// read n(layers) from detector metadata
    void GetNLayersFromGeom();

    /**  Determine the local position based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   The position is in the local coordinate system of the associated
     *   dd4hep volume.
     */
    virtual Vector3D position(const CellID& aCellID) const override;
    /**  Determine the cell ID based on the position.
     *   @param[in] aLocalPosition (not used).
     *   @param[in] aGlobalPosition
     *   @param[in] aVolumeId ID of the Geant4 volume
     *   return Cell ID.
     */
    virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                          const VolumeID& aVolumeID) const override;
    /**  Determine the azimuthal angle (relative to the G4 volume) based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   return Phi.
     */
    double phi(const CellID aCellID) const;
    /**  Determine the polar angle (relative to the G4 volume) based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   return Theta.
     */
    double theta(const CellID aCellID) const;
    /**  Determine the radius based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   return Radius.
     */
    // double radius(const CellID aCellID) const;
    /**  Get the number of merged cells in theta for given layer
     *   @param[in] layer
     *   return The number of merged cells in theta
     */
    inline int mergedThetaCells(const int layer) const {
      if (layer < (int)m_mergedCellsTheta.size())
        return m_mergedCellsTheta[layer];
      else
        return 1;
    }
    /**  Get the number of merged modules (inclined in phi)
     *   @param[in] layer
     *   return The number of merged modules
     */
    inline int mergedModules(const int layer) const {
      if (layer < (int)m_mergedModules.size())
        return m_mergedModules[layer];
      else
        return 1;
    }
    /**  Get the total number of modules of detector
     *   return The number of modules (as it was set by the user in the xml file..)
     */
    inline int nModules() const { return m_nModules; }
    /**  Get the number of layers
     *   return The number of layers
     */
    inline int nLayers() const { return m_nLayers; }
    /**  Get the field name used for the layer
     *   return The field name for the layer.
     */
    inline const std::string& fieldNameLayer() const { return m_layerID; }
    /**  Get the field name used for the module number
     *   return The field name for the module.
     */
    inline const std::string& fieldNameModule() const { return m_moduleID; }

    /// Extract the layer index fom a cell ID.
    int layer(const CellID aCellID) const;

    /// Determine the volume ID from the full cell ID by removing all local fields
    virtual VolumeID volumeID(const CellID& cellID) const override;

    /// Return true if this segmentation can have cells that span multiple
    /// volumes.  That is, points from multiple distinct volumes may
    /// be assigned to the same cell.
    virtual bool cellsSpanVolumes() const override { return true; }

    /** Returns a std::vector<double> of the cellDimensions of the given cell ID
     *  in natural order of dimensions (nModules, dTheta)
     *  @param[in] cellID
     *  return a std::vector of size 2 with the cellDimensions of the given cell ID(modules, theta)
     */
    virtual std::vector<double> cellDimensions(const CellID& id) const override {
      const int aLayer = layer(id);
      return {(double)mergedModules(aLayer), gridSizeTheta() * mergedThetaCells(aLayer)};
    }

  protected:
    /// the field name used for the cryostat
    std::string m_cryoID;
    /// the field name used for layer
    std::string m_layerID;
    /// the field name used for the read-out module (can differ from module due to merging)
    std::string m_moduleID;
    /// vector of number of cells to be merged along theta for each layer (typically between 1 and 4)
    std::vector<int> m_mergedCellsTheta;
    /// vector of number of modules to be merged for each layer (typically 1 or 2)
    std::vector<int> m_mergedModules;

    /// number of modules (or, equivalently, the deltaPhi between adjacent modules)
    int m_nModules;

    /// number of layers (from the geometry)
    int m_nLayers;

  private:
    /// Tabulate the cylindrical radii of all layers, as well as the
    /// local x and z components needed for the proper phi offset.
    struct LayerInfo {
      LayerInfo(double the_rho, double the_xloc, double the_zloc) : rho(the_rho), xloc(the_xloc), zloc(the_zloc) {}
      double rho;
      double xloc;
      double zloc;
    };
    std::vector<LayerInfo> initLayerInfo(const CellID& cID) const;

    // The vector of tabulated values, indexed by layer number.
    // We can't build this in the constructor --- the volumes won't have
    // been created yet.  Instead, build it lazily the first time it's needed.
    // Since that's in a const method, make it thread-safe.
    mutable std::atomic<const std::vector<LayerInfo>*> m_layerInfo = nullptr;

    /// Initialization common to all ctors.
    void commonSetup();
    /// the field index used for the cryostat
    int m_cryoIndex = -1;
    /// the field index used for layer
    int m_layerIndex = -1;
    /// the field index used for theta
    int m_thetaIndex = -1;
    /// the field index used for module
    int m_moduleIndex = -1;
  };
} // namespace DDSegmentation
} // namespace dd4hep
#endif /* DETSEGMENTATION_FCCSWGridModuleThetaMerged_k4geo_H */

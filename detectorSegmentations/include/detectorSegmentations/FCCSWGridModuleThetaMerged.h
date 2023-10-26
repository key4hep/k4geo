#ifndef DETSEGMENTATION_FCCSWGRIDMODULETHETAMERGED_H
#define DETSEGMENTATION_FCCSWGRIDMODULETHETAMERGED_H

// FCCSW
#include "detectorSegmentations/GridTheta_k4geo.h"
#include "DD4hep/VolumeManager.h"

/** FCCSWGridModuleThetaMerged Detector/DetSegmentation/DetSegmentation/FCCSWGridModuleThetaMerged.h FCCSWGridModuleThetaMerged.h
 *
 *  Segmentation in theta and module.
 *  Based on GridTheta, merges modules and theta cells based on layer number
 *
 */

namespace dd4hep {
namespace DDSegmentation {
class FCCSWGridModuleThetaMerged : public GridTheta_k4geo {
public:
  /// default constructor using an arbitrary type
  FCCSWGridModuleThetaMerged(const std::string& aCellEncoding);
  /// Default constructor used by derived classes passing an existing decoder
  FCCSWGridModuleThetaMerged(const BitFieldCoder* decoder);

  /// destructor
  virtual ~FCCSWGridModuleThetaMerged() = default;

  /// read n(modules) from detector metadata
  void GetNModulesFromGeom();

  /// read n(layers) from detector metadata
  void GetNLayersFromGeom();
  
  /**  Determine the local position based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Position (relative to R, phi of Geant4 volume it belongs to, scaled for R=1).
   */
  virtual Vector3D position(const CellID& aCellID) const;
  /**  Determine the cell ID based on the position.
   *   @param[in] aLocalPosition (not used).
   *   @param[in] aGlobalPosition
   *   @param[in] aVolumeId ID of the Geant4 volume
   *   return Cell ID.
   */
  virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;
  /**  Determine the azimuthal angle (relative to the G4 volume) based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Phi.
   */
  double phi(const CellID& aCellID) const;
  /**  Determine the polar angle (relative to the G4 volume) based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Theta.
   */
  double theta(const CellID& aCellID) const;
  /**  Determine the radius based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Radius.
   */
  // double radius(const CellID& aCellID) const;
  /**  Get the number of merged cells in theta for given layer
   *   @param[in] layer
   *   return The number of merged cells in theta
   */
  inline int mergedThetaCells(const int layer) const {
    if (layer < (int) m_mergedCellsTheta.size())
      return m_mergedCellsTheta[layer];
    else
      return 1;
  }
  /**  Get the number of merged modules (inclined in phi)
   *   @param[in] layer
   *   return The number of merged modules
   */
  inline int mergedModules(const int layer) const {
    if (layer < (int) m_mergedModules.size())
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

protected:
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
  
};
}
}
#endif /* DETSEGMENTATION_FCCSWGRIDMODULETHETAMERGED_H */

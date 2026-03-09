#ifndef DETSEGMENTATION_ENDCAPTURBINE_H
#define DETSEGMENTATION_ENDCAPTURBINE_H

// FCCSW
// #include "detectorSegmentations/GridTheta_k4geo.h"
#include "DDSegmentation/Segmentation.h"
#include "TVector3.h"

/** FCCSWEndcapTurbine_k4geo
 *
 *  Segmentation for turbine-style endcap calorimeter, allowing for
 *  modules to be merged to reduce the readout channel count
 *
 */

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWEndcapTurbine_k4geo : public Segmentation {
  public:
    /// default constructor using an arbitrary type
    FCCSWEndcapTurbine_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    FCCSWEndcapTurbine_k4geo(const BitFieldCoder* decoder);
    /// common setup that needs to be done for either of the above constructors
    void commonSetup();
    /// destructor
    virtual ~FCCSWEndcapTurbine_k4geo() = default;

    /**  Determine the global position based on the cell ID. **/
    virtual Vector3D position(const CellID& aCellID) const override;
    /**  Determine the cell ID based on the position.
     *   @param[in] aLocalPosition (not used).
     *   @param[in] aGlobalPosition position in the global coordinates.
     *   @param[in] aVolumeId ID of a volume.
     *   return Cell ID.
     */
    virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                          const VolumeID& aVolumeID) const override;

    /**  Determine the transverse distance from the beamline (rho) based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   return rho.
     */
    double rho(const CellID aCellID) const;
    /** Get the grid size in rho for a given wheel
     * return grid size in rho
     */
    inline double gridSizeRho(int iWheel) const { return m_gridSizeRho[iWheel]; }
    /** Get the number of cells in rho for a given wheel
     * @param[in] iWheel wheel index
     * return number of cells in rho for the specified wheel
     */
    inline int numCellsRho(int iWheel) const { return m_numReadoutRhoLayers[iWheel]; }

    /** Get the number of calibration cells in rho for a given wheel
     * @param[in] iWheel wheel index
     * return number of calibration cells in rho for the specified wheel
     */
    inline int numCellsRhoCalib(int iWheel) const { return m_numCalibRhoLayers[iWheel]; }

    /**  Get the coordinate offset in rho for a given wheel.
     *   return The offset in rho.
     */
    inline double offsetRho(int iWheel) const { return m_offsetRho[iWheel]; }
    /**  Get the field name for rho.
     *   return The field name for rho.
     */

    inline const std::string& fieldNameRho() const { return m_rhoID; }
    /**  Determine the azimuthal angle based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   return Phi.
     */
    double phi(const CellID aCellID) const;

    /**  Get the coordinate offset in azimuthal angle.
     *   return The offset in phi.
     */
    inline double offsetPhi() const { return m_offsetPhi; }
    /**  Get the coordinate offset in theta angle.
     *   return The offset in theta.
     */
    inline double offsetTheta() const { return m_offsetTheta; }
    /**  Get the field name for azimuthal angle.
     *   return The field name for phi.
     */
    inline const std::string& fieldNamePhi() const { return m_phiID; }
    /** Get the angle of the turbine blades in a given wheel
     *   @param[in] iWheel index of wheel.
     *  return the blade angle for the requested wheel
     */
    double bladeAngle(unsigned iWheel) const { return m_bladeAngle[iWheel]; }
    /**  Set the number of bins in azimuthal angle.
     *   @param[in] aNumberBins Number of bins in phi.
     */
    inline void setPhiBins(int bins) { m_phiBins = bins; }
    /**  Set the coordinate offset in azimuthal angle.
     *   @param[in] aOffset Offset in phi.
     */
    inline void setOffsetPhi(double offset) { m_offsetPhi = offset; }
    /**  Set the coordinate offset in theta angle.
     *   @param[in] aOffset Offset in theta.
     */
    inline void setOffsetTheta(double offset) { m_offsetTheta = offset; }
    /**  Set the field name used for phi.
     *   @param[in] aFieldName Field name for phi.
     */
    inline void setFieldNamePhi(const std::string& fieldName) { m_phiID = fieldName; }
    /**  Set the field name used for the wheel ID.
     *   @param[in] aFieldName Field name for wheel.
     */
    inline void setFieldNameWheel(const std::string& fieldName) { m_wheelID = fieldName; }
    /**  Determine the x coordinate based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   return x.
     */
    double x(const CellID aCellID) const;
    /**  Determine the z coordinate based on the cell ID.
     *   @param[in] aCellId ID of a cell.
     *   return z.
     */
    double z(const CellID aCellID) const;
    /** Get the grid size in z for a given wheel
     * return grid size in z
     */
    inline double gridSizeZ(int iWheel) const { return m_gridSizeZ[iWheel]; }
    /**  Get the coordinate offset in z.
     *   return The offset in z.
     */
    /** Get the number of cells in z for a given wheel
     * @param[in] iWheel wheel index
     * return number of cells in z for the specified wheel
     */
    inline int numCellsZ(int iWheel) const { return m_numReadoutZLayers[iWheel]; }
    /** Get the number of calibration cells in z for a given wheel
     * @param[in] iWheel wheel index
     * return number of calibration cells in z for the specified wheel
     */
    inline int numCellsZCalib(int iWheel) const {
      return m_numCalibZLayers[iWheel];
    } /** Get the offset in z for a given wheel
       * @param[in] iWheel wheel index
       * return offset in z for the specified wheel
       */
    inline double offsetZ(int iWheel) const { return m_offsetZ[iWheel]; }
    /**  Get the field name for z.
     *   return The field name for z.
     */
    inline const std::string& fieldNameZ() const { return m_zID; }
    /**  Set the number of bins in z.
     *   @param[in] aNumberBins Number of bins in z.
     */
    inline void setZBins(int bins) { m_zBins = bins; }
    /**  Set the coordinate offset in z for the specified wheel.
     *   @param[in] iWheel wheel index
     *   @param[in] aOffset Offset in z.
     */
    inline void setOffsetZ(int iWheel, double offset) { m_offsetZ[iWheel] = offset; }
    /**  Set the field name used for z.
     *   @param[in] aFieldName Field name for z.
     */
    inline void setFieldNameZ(const std::string& fieldName) { m_zID = fieldName; }
    inline double rhoFromXYZ(const Vector3D& aposition) const {
      TVector3 vec(aposition.X, aposition.Y, aposition.Z);
      return vec.Perp();
    }
    inline double phiFromXYZ(const Vector3D& aposition) const {
      TVector3 vec(aposition.X, aposition.Y, aposition.Z);
      return vec.Phi();
    }

    /** return the number of unit cells in each wheel
     * @param[in] iWheel wheel identifier
     */
    inline int nModules(int iWheel) const { return m_nUnitCells[iWheel] / m_mergedModules[iWheel]; }
    /** return the expected value of the layer index
     * @param[in] iWheel wheel identifier
     * @param[in] iRho rho readout cell identifier
     * @param[in] iZ z readout cell identifier
     */
    unsigned expLayer(unsigned iWheel, unsigned iRho, unsigned iZ) const;

    /**  Get the field name for side.
     *   return The field name for side.
     */
    inline const std::string& fieldNameSide() const { return m_sideID; }
    /**  Get the field name for the wheel.
     *   return The field name for wheel.
     */
    inline const std::string& fieldNameWheel() const { return m_wheelID; }
    /**  Get the field name for the module.
     *   return The field name for module.
     */
    inline const std::string& fieldNameModule() const { return m_moduleID; }
    /**  Get the field name for the layer.
     *   return The field name for layer.
     */
    inline const std::string& fieldNameLayer() const { return m_layerID; }

  private:
    /// turbine blade angle in each wheel
    std::vector<double> m_bladeAngle;
    /// number of unit cells in each wheel
    std::vector<int> m_nUnitCells;
    /// number of merged modules in each wheel
    std::vector<int> m_mergedModules;
    /// the number of cells in rho for each wheel
    std::vector<int> m_numReadoutRhoLayers;
    /// the number of cells in z for each wheel
    std::vector<int> m_numReadoutZLayers;
    /// the number of calibration cells in rho for each wheel
    std::vector<int> m_numCalibRhoLayers;
    /// the number of calibration cells in z for each wheel
    std::vector<int> m_numCalibZLayers;
    /// the number of bins in phi
    int m_phiBins;
    /// the coordinate offset in phi
    double m_offsetPhi;
    /// the coordinate offset in theta
    double m_offsetTheta; /// the field name used for phi
    std::string m_phiID;
    /// the number of bins in rho
    int m_rhoBins;
    ////grid size in rho
    std::vector<double> m_gridSizeRho;
    /// the coordinate offset in rho
    std::vector<double> m_offsetRho;
    /// the field name used for rho
    std::string m_rhoID;
    /// the field name used for wheel
    std::string m_wheelID;
    /// the field name used for module
    std::string m_moduleID;
    /// the number of bins in z
    int m_zBins;
    /// grid size in z
    std::vector<double> m_gridSizeZ;
    /// the coordinate offset in z
    std::vector<double> m_offsetZ;
    /// the field name used for z
    std::string m_zID;
    std::string m_sideID;
    std::string m_layerID;

    /// the field index used for rho
    int m_rhoIndex = -1;
    /// the field index used for wheel
    int m_wheelIndex = -1;
    /// the field index used for module
    int m_moduleIndex = -1;
    /// the field index used for z
    int m_zIndex = -1;
    /// the field index used for side
    int m_sideIndex = -1;
    /// the field index used for layer
    int m_layerIndex = -1;
  };
} // namespace DDSegmentation
} // namespace dd4hep
#endif /* DETSEGMENTATION_ENDCAPTURBINE_H */

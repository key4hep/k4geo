#ifndef DETSEGMENTATION_GRIDDRCALO_H
#define DETSEGMENTATION_GRIDDRCALO_H

#include "detectorSegmentations/DRparamBarrel_k4geo.h"
#include "detectorSegmentations/DRparamEndcap_k4geo.h"

#include "DDSegmentation/Segmentation.h"

namespace dd4hep {
namespace DDSegmentation {
  class GridDRcalo_k4geo : public Segmentation {
  public:
    /// default constructor using an arbitrary type
    GridDRcalo_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    GridDRcalo_k4geo(const BitFieldCoder* decoder);
    /// destructor
    virtual ~GridDRcalo_k4geo() override;

    // Determine the global (local) position based on the cell ID.
    // position of the front end of the fiber
    virtual Vector3D position(const CellID& aCellID) const override;
    Vector3D localPosition(const CellID aCellID) const;
    Vector3D localPosition(int numx, int numy, int x_, int y_) const;

    // position of the rear end (readout location) of the fiber
    Vector3D sipmPosition(const CellID& aCellID) const;

    virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                          const VolumeID& aVolumeID) const override;

    VolumeID setVolumeID(int systemId, int numEta, int numPhi) const;
    CellID setCellID(bool isRHS, int systemId, int numEta, int numPhi, int x, int y) const;

    virtual void neighbours(const CellID& cellID, std::set<CellID>& neighbours) const override;

    void setGridSize(double grid) { m_gridSize = grid; }
    void setSipmSize(double sipm) { m_sipmSize = sipm; }

    // remove different type of channels in the neighborhood set if set to true
    void setRemoveDifferentCh(bool rmDiffCh) { m_removeDifferentCh = rmDiffCh; }

    // Get the identifier number of a mother tower in eta or phi direction
    int numEta(const CellID aCellID) const;
    int numPhi(const CellID aCellID) const;

    // Get the total number of SiPMs of the mother tower in x or y direction (local coordinate)
    int numX(const CellID aCellID) const;
    int numY(const CellID aCellID) const;

    // Get the identifier number of a SiPM in x or y direction (local coordinate)
    int x(const CellID aCellID) const; // approx phi direction
    int y(const CellID aCellID) const; // approx eta direction

    bool IsCerenkov(const CellID aCellID) const;
    bool IsCerenkov(int col, int row) const;

    bool IsTower(const CellID aCellID) const;
    bool IsSiPM(const CellID aCellID) const;
    bool IsRHS(const CellID aCellID) const;

    int getFirst32bits(const CellID aCellID) const { return (int)aCellID; }
    int getLast32bits(const CellID aCellID) const;
    CellID convertFirst32to64(const int aId32) const { return (CellID)aId32; }
    CellID convertLast32to64(const int aId32) const;

    // Methods for 32bit to 64bit en/decoder
    int numEta(const int aId32) const { return numEta(convertFirst32to64(aId32)); }
    int numPhi(const int aId32) const { return numPhi(convertFirst32to64(aId32)); }

    int numX(const int aId32) const { return numX(convertFirst32to64(aId32)); }
    int numY(const int aId32) const { return numY(convertFirst32to64(aId32)); }

    int x(const int aId32) const { return x(convertLast32to64(aId32)); }
    int y(const int aId32) const { return y(convertLast32to64(aId32)); }

    bool IsCerenkov(const int aId32) const { return IsCerenkov(convertLast32to64(aId32)); }

    bool IsTower(const int aId32) const { return IsTower(convertLast32to64(aId32)); }
    bool IsSiPM(const int aId32) const { return IsSiPM(convertLast32to64(aId32)); }

    inline const std::string& fieldNameAssembly() const { return m_assemblyID; }
    inline const std::string& fieldNameNumEta() const { return m_numEtaID; }
    inline const std::string& fieldNameNumPhi() const { return m_numPhiID; }
    inline const std::string& fieldNameX() const { return m_xID; }
    inline const std::string& fieldNameY() const { return m_yID; }
    inline const std::string& fieldNameIsCerenkov() const { return m_isCerenkovID; }
    inline const std::string& fieldNameModule() const { return m_moduleID; }

    inline void setFieldNameAssembly(const std::string& fieldName) { m_assemblyID = fieldName; }
    inline void setFieldNameNumEta(const std::string& fieldName) { m_numEtaID = fieldName; }
    inline void setFieldNameNumPhi(const std::string& fieldName) { m_numPhiID = fieldName; }
    inline void setFieldNameX(const std::string& fieldName) { m_xID = fieldName; }
    inline void setFieldNameY(const std::string& fieldName) { m_yID = fieldName; }
    inline void setFieldNameIsCerenkov(const std::string& fieldName) { m_isCerenkovID = fieldName; }
    inline void setFieldNameModule(const std::string& fieldName) { m_moduleID = fieldName; }

    DRparamBarrel_k4geo* paramBarrel() { return m_paramBarrel; }
    DRparamEndcap_k4geo* paramEndcap() { return m_paramEndcap; }

    DRparamBase_k4geo* setParamBase(int noEta) const;

  private:
    std::string m_assemblyID;
    std::string m_numEtaID;
    std::string m_numPhiID;
    std::string m_xID;
    std::string m_yID;
    std::string m_isCerenkovID;
    std::string m_moduleID;

    double m_gridSize;
    double m_sipmSize;

    // switch to remove different type cells from the neighborhood
    bool m_removeDifferentCh;

    DRparamBarrel_k4geo* m_paramBarrel;
    DRparamEndcap_k4geo* m_paramEndcap;

    /// Initialization common to all ctors.
    void commonSetup();
    /// the field index used for system
    int m_systemIndex = -1;
    /// the field index used for numEta
    int m_numEtaIndex = -1;
    /// the field index used for numPhi
    int m_numPhiIndex = -1;
    /// the field index used for module
    int m_moduleIndex = -1;
    /// the field index used for x
    int m_xIndex = -1;
    /// the field index used for y
    int m_yIndex = -1;
    /// the field index used for isCerenkov
    int m_isCerenkovIndex = -1;
    /// the field index used for assembly
    int m_assemblyIndex = -1;
  };
} // namespace DDSegmentation
} // namespace dd4hep

#endif

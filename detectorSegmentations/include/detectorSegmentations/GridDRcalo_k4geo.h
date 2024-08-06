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

  //  Determine the global(local) position based on the cell ID.
  virtual Vector3D position(const CellID& aCellID) const;
  Vector3D localPosition(const CellID& aCellID) const;
  Vector3D localPosition(int numx, int numy, int x_, int y_) const;

  virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;

  VolumeID setVolumeID(int numEta, int numPhi) const;
  CellID setCellID(int numEta, int numPhi, int x, int y) const;

  void setGridSize(double grid) { fGridSize = grid; }
  void setSipmSize(double sipm) { fSipmSize = sipm; }

  // Get the identifier number of a mother tower in eta or phi direction
  int numEta(const CellID& aCellID) const;
  int numPhi(const CellID& aCellID) const;

  // Get the total number of SiPMs of the mother tower in x or y direction (local coordinate)
  int numX(const CellID& aCellID) const;
  int numY(const CellID& aCellID) const;

  // Get the identifier number of a SiPM in x or y direction (local coordinate)
  int x(const CellID& aCellID) const; // approx eta direction
  int y(const CellID& aCellID) const; // approx phi direction

  bool IsCerenkov(const CellID& aCellID) const;
  bool IsCerenkov(int col, int row) const;

  bool IsTower(const CellID& aCellID) const;
  bool IsSiPM(const CellID& aCellID) const;

  int getFirst32bits(const CellID& aCellID) const { return (int)aCellID; }
  int getLast32bits(const CellID& aCellID) const;
  CellID convertFirst32to64(const int aId32) const { return (CellID)aId32; }
  CellID convertLast32to64(const int aId32) const;

  // Methods for 32bit to 64bit en/decoder
  int numEta(const int& aId32) const { return numEta( convertFirst32to64(aId32) ); }
  int numPhi(const int& aId32) const { return numPhi( convertFirst32to64(aId32) ); }

  int numX(const int& aId32) const { return numX( convertFirst32to64(aId32) ); }
  int numY(const int& aId32) const { return numY( convertFirst32to64(aId32) ); }

  int x(const int& aId32) const { return x( convertLast32to64(aId32) ); }
  int y(const int& aId32) const { return y( convertLast32to64(aId32) ); }

  bool IsCerenkov(const int& aId32) const { return IsCerenkov( convertLast32to64(aId32) ); }

  bool IsTower(const int& aId32) const { return IsTower( convertLast32to64(aId32) ); }
  bool IsSiPM(const int& aId32) const { return IsSiPM( convertLast32to64(aId32) ); }

  inline const std::string& fieldNameAssembly() const { return fAssemblyId; }
  inline const std::string& fieldNameNumEta() const { return fNumEtaId; }
  inline const std::string& fieldNameNumPhi() const { return fNumPhiId; }
  inline const std::string& fieldNameX() const { return fXId; }
  inline const std::string& fieldNameY() const { return fYId; }
  inline const std::string& fieldNameIsCerenkov() const { return fIsCerenkovId; }
  inline const std::string& fieldNameModule() const { return fModule; }

  inline void setFieldNameAssembly(const std::string& fieldName) { fAssemblyId = fieldName; }
  inline void setFieldNameNumEta(const std::string& fieldName) { fNumEtaId = fieldName; }
  inline void setFieldNameNumPhi(const std::string& fieldName) { fNumPhiId = fieldName; }
  inline void setFieldNameX(const std::string& fieldName) { fXId = fieldName; }
  inline void setFieldNameY(const std::string& fieldName) { fYId = fieldName; }
  inline void setFieldNameIsCerenkov(const std::string& fieldName) { fIsCerenkovId = fieldName; }
  inline void setFieldNameModule(const std::string& fieldName) { fModule = fieldName; }

  DRparamBarrel_k4geo* paramBarrel() { return fParamBarrel; }
  DRparamEndcap_k4geo* paramEndcap() { return fParamEndcap; }

  DRparamBase_k4geo* setParamBase(int noEta) const;

protected:
  std::string fAssemblyId;
  std::string fNumEtaId;
  std::string fNumPhiId;
  std::string fXId;
  std::string fYId;
  std::string fIsCerenkovId;
  std::string fModule;

  double fGridSize;
  double fSipmSize;

private:
  DRparamBarrel_k4geo* fParamBarrel;
  DRparamEndcap_k4geo* fParamEndcap;
};
}
}

#endif

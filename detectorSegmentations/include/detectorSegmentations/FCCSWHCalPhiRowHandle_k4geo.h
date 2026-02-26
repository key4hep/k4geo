#ifndef DETECTORSEGMENTATIONS_HCALPHIROWHANDLE_K4GEO_H
#define DETECTORSEGMENTATIONS_HCALPHIROWHANDLE_K4GEO_H

// FCCSW
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"

// DD4hep
#include "DD4hep/Segmentations.h"
#include "DD4hep/detail/SegmentationsInterna.h"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
/// Namespace for base segmentations

// Forward declarations
class Segmentation;
template <typename T>
class SegmentationWrapper;

/// We need some abbreviation to make the code more readable.
typedef Handle<SegmentationWrapper<DDSegmentation::FCCSWHCalPhiRow_k4geo>> FCCSWHCalPhiRowHandle_k4geo;

/// Implementation class for the HCal phi-row segmentation.
/**
 *  Concrete user handle to serve specific needs of client code
 *  which requires access to the base functionality not served
 *  by the super-class Segmentation.
 *
 *  Note:
 *  We only check the validity of the underlying handle.
 *  If for whatever reason the implementation object is not valid
 *  This is not checked.
 *  In principle this CANNOT happen unless some brain-dead has
 *  fiddled with the handled object directly.....
 *
 *  Note:
 *  The handle base corrsponding to this object in for
 *  conveniance reasons instantiated in DD4hep/src/Segmentations.cpp.
 *
 */
class FCCSWHCalPhiRow_k4geo : public FCCSWHCalPhiRowHandle_k4geo {
public:
  /// Defintion of the basic handled object
  typedef FCCSWHCalPhiRowHandle_k4geo::Object Object;

public:
  /// Default constructor
  FCCSWHCalPhiRow_k4geo() = default;
  /// Copy constructor
  FCCSWHCalPhiRow_k4geo(const FCCSWHCalPhiRow_k4geo& e) = default;
  /// Copy Constructor from segmentation base object
  FCCSWHCalPhiRow_k4geo(const Segmentation& e) : Handle<Object>(e) {}
  /// Copy constructor from handle
  FCCSWHCalPhiRow_k4geo(const Handle<Object>& e) : Handle<Object>(e) {}
  /// Copy constructor from other polymorph/equivalent handle
  template <typename Q>
  FCCSWHCalPhiRow_k4geo(const Handle<Q>& e) : Handle<Object>(e) {}
  /// Assignment operator
  FCCSWHCalPhiRow_k4geo& operator=(const FCCSWHCalPhiRow_k4geo& seg) = default;
  /// Equality operator
  bool operator==(const FCCSWHCalPhiRow_k4geo& seg) const { return m_element == seg.m_element; }

  /// determine the position based on the cell ID
  inline Position position(const CellID id) const { return Position(access()->implementation->position(id)); }

  /// determine the cell ID based on the position
  inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID volID) const {
    return access()->implementation->cellID(local, global, volID);
  }

  /// access the grid size in row for each layer
  inline const std::vector<int>& gridSizeRow() const { return access()->implementation->gridSizeRow(); }

  /// access the number of rows grouped in the pseudo-layers of the Endcap
  inline const std::vector<int>& groupedRows() const { return access()->implementation->groupedRows(); }

  /// access the grid size in Phi
  inline int phiBins() const { return access()->implementation->phiBins(); }

  /// access the coordinate offset in Phi
  inline double offsetPhi() const { return access()->implementation->offsetPhi(); }

  /// set the coordinate offset in Phi
  inline void setOffsetPhi(double offset) const { access()->implementation->setOffsetPhi(offset); }

  /// set the detector layout
  inline void setDetLayout(int detLayout) const { access()->implementation->setDetLayout(detLayout); }

  /// set the coordinate offset in z-axis
  inline void setOffsetZ(std::vector<double> const& offset) const { access()->implementation->setOffsetZ(offset); }

  /// set the z width
  inline void setWidthZ(std::vector<double> const& width) const { access()->implementation->setWidthZ(width); }

  /// set the offset in radius
  inline void setOffsetR(std::vector<double> const& offset) const { access()->implementation->setOffsetR(offset); }

  /// set the number of layers with different dR
  inline void setNumLayers(std::vector<int> const& num) const { access()->implementation->setNumLayers(num); }

  /// set the dR of each layer
  inline void setdRlayer(std::vector<double> const& dRlayer) const { access()->implementation->setdRlayer(dRlayer); }

  /// set the grid size in theta
  inline void setGridSizeRow(std::vector<int> const& cellSize) const {
    access()->implementation->setGridSizeRow(cellSize);
  }

  /// set the number of rows grouped in the pseudo-layer of the Endcap
  inline void setGroupedRows(std::vector<int> const& cellSize) const {
    access()->implementation->setGroupedRows(cellSize);
  }

  /// set the grid size in Phi
  inline void setPhiBins(int cellSize) const { access()->implementation->setPhiBins(cellSize); }

  /// access the field name used for layer
  inline const std::string& fieldNameLayer() const { return access()->implementation->fieldNameLayer(); }

  /// access the field name used for theta
  inline const std::string& fieldNameRow() const { return access()->implementation->fieldNameRow(); }

  /// access the field name used for Phi
  inline const std::string& fieldNamePhi() const { return access()->implementation->fieldNamePhi(); }
};

} /* End namespace dd4hep */
#endif // DETECTORSEGMENTATIONS_HCALPHITHETAHANDLE_K4GEO_H

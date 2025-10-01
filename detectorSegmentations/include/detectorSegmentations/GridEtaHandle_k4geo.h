#ifndef DETECTORSEGMENTATIONS_GRIDETAHANDLE_K4GEO_H
#define DETECTORSEGMENTATIONS_GRIDETAHANDLE_K4GEO_H

// FCCSW
#include "detectorSegmentations/GridEta_k4geo.h"

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
typedef Handle<SegmentationWrapper<DDSegmentation::GridEta_k4geo>> GridEtaHandle_k4geo;

/// Implementation class for the grid phi-eta segmentation.
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
 *  The handle base corresponding to this object is for
 *  convenience reasons instantiated in DD4hep/src/Segmentations.cpp.
 *
 *  \author  A. Zaborowska
 *  \version 1.0
 */
class GridEta_k4geo : public GridEtaHandle_k4geo {
public:
  /// Defintiion of the basic handled object
  typedef GridEtaHandle_k4geo::Object Object;

public:
  /// Default constructor
  GridEta_k4geo() = default;
  /// Copy constructor
  GridEta_k4geo(const GridEta_k4geo& e) = default;
  /// Copy Constructor from segmentation base object
  GridEta_k4geo(const Segmentation& e) : Handle<Object>(e) {}
  /// Copy constructor from handle
  GridEta_k4geo(const Handle<Object>& e) : Handle<Object>(e) {}
  /// Copy constructor from other polymorph/equivalent handle
  template <typename Q>
  GridEta_k4geo(const Handle<Q>& e) : Handle<Object>(e) {}
  /// Assignment operator
  GridEta_k4geo& operator=(const GridEta_k4geo& seg) = default;
  /// Equality operator
  bool operator==(const GridEta_k4geo& seg) const { return m_element == seg.m_element; }
  /// determine the position based on the cell ID
  inline Position position(const CellID id) const { return Position(access()->implementation->position(id)); }

  /// determine the cell ID based on the position
  inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID volID) const {
    return access()->implementation->cellID(local, global, volID);
  }

  /// access the grid size in eta
  inline double gridSizeEta() const { return access()->implementation->gridSizeEta(); }

  /// access the coordinate offset in eta
  inline double offsetEta() const { return access()->implementation->offsetEta(); }

  /// set the coordinate offset in eta
  inline void setOffsetEta(double offset) const { access()->implementation->setOffsetEta(offset); }

  /// set the grid size in eta
  inline void setGridSizeEta(double cellSize) const { access()->implementation->setGridSizeEta(cellSize); }

  /// access the field name used for eta
  inline const std::string& fieldNameEta() const { return access()->implementation->fieldNameEta(); }

  /** \brief Returns a std::vector<double> of the cellDimensions of the given cell ID
      in natural order of dimensions (dEta)

      Returns a std::vector of the cellDimensions of the given cell ID
      \param cellID is ignored as all cells have the same dimension
      \return std::vector<double> size 1:
      -# size in eta
  */
  inline std::vector<double> cellDimensions(const CellID /*id*/) const {
    return {access()->implementation->gridSizeEta()};
  }
};

} /* End namespace dd4hep                */
#endif // DD4HEP_DDCORE_GRIDETA_H

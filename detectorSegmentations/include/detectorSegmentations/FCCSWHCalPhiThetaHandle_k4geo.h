#ifndef DETECTORSEGMENTATIONS_HCALPHITHETAHANDLE_K4GEO_H
#define DETECTORSEGMENTATIONS_HCALPHITHETAHANDLE_K4GEO_H

// FCCSW
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"

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
typedef Handle<SegmentationWrapper<DDSegmentation::FCCSWHCalPhiTheta_k4geo>> FCCSWHCalPhiThetaHandle_k4geo;

/// Implementation class for the HCal phi-theta segmentation.
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
class FCCSWHCalPhiTheta_k4geo : public FCCSWHCalPhiThetaHandle_k4geo {
public:
  /// Defintion of the basic handled object
  typedef FCCSWHCalPhiThetaHandle_k4geo::Object Object;

public:
  /// Default constructor
  FCCSWHCalPhiTheta_k4geo() = default;
  /// Copy constructor
  FCCSWHCalPhiTheta_k4geo(const FCCSWHCalPhiTheta_k4geo& e) = default;
  /// Copy Constructor from segmentation base object
  FCCSWHCalPhiTheta_k4geo(const Segmentation& e) : Handle<Object>(e) {}
  /// Copy constructor from handle
  FCCSWHCalPhiTheta_k4geo(const Handle<Object>& e) : Handle<Object>(e) {}
  /// Copy constructor from other polymorph/equivalent handle
  template <typename Q>
  FCCSWHCalPhiTheta_k4geo(const Handle<Q>& e) : Handle<Object>(e) {}
  /// Assignment operator
  FCCSWHCalPhiTheta_k4geo& operator=(const FCCSWHCalPhiTheta_k4geo& seg) = default;
  /// Equality operator
  bool operator==(const FCCSWHCalPhiTheta_k4geo& seg) const { return m_element == seg.m_element; }

  /// determine the position based on the cell ID
  inline Position position(const CellID id) const { return Position(access()->implementation->position(id)); }

  /// determine the cell ID based on the position
  inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID volID) const {
    return access()->implementation->cellID(local, global, volID);
  }

  /// access the grid size in theta
  inline double gridSizeTheta() const { return access()->implementation->gridSizeTheta(); }

  /// access the grid size in Phi
  inline int phiBins() const { return access()->implementation->phiBins(); }

  /// access the coordinate offset in theta
  inline double offsetTheta() const { return access()->implementation->offsetTheta(); }

  /// access the coordinate offset in Phi
  inline double offsetPhi() const { return access()->implementation->offsetPhi(); }

  /// set the coordinate offset in theta
  inline void setOffsetTheta(double offset) const { access()->implementation->setOffsetTheta(offset); }

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
  inline void setGridSizeTheta(double cellSize) const { access()->implementation->setGridSizeTheta(cellSize); }

  /// set the grid size in Phi
  inline void setPhiBins(int cellSize) const { access()->implementation->setPhiBins(cellSize); }

  /// access the field name used for theta
  inline const std::string& fieldNameTheta() const { return access()->implementation->fieldNameTheta(); }

  /// access the field name used for Phi
  inline const std::string& fieldNamePhi() const { return access()->implementation->fieldNamePhi(); }

  /// access the field name used for layer
  inline const std::string& fieldNameLayer() const { return access()->implementation->fieldNameLayer(); }

  /** \brief Returns a std::vector<double> of the cellDimensions of the given cell ID
      in natural order of dimensions (dPhi, dTheta)

      Returns a std::vector of the cellDimensions of the given cell ID
      \param cellID is ignored as all cells have the same dimension
      \return std::vector<double> size 2:
      -# size in phi
      -# size in theta
  */
  inline std::vector<double> cellDimensions(const CellID /*id*/) const {
    return {access()->implementation->gridSizePhi(), access()->implementation->gridSizeTheta()};
  }
};

} /* End namespace dd4hep */
#endif // DETECTORSEGMENTATIONS_HCALPHITHETAHANDLE_K4GEO_H

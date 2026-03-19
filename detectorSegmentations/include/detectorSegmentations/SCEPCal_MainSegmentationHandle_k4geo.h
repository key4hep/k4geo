//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#ifndef SCEPCal_MainSegmentationHandle_k4geo_h
#define SCEPCal_MainSegmentationHandle_k4geo_h 1
#include "DD4hep/Segmentations.h"
#include "DD4hep/detail/SegmentationsInterna.h"
#include "Math/Vector3D.h"
#include "SCEPCal_MainSegmentation_k4geo.h"

namespace dd4hep {

class Segmentation;
template <typename T>
class SegmentationWrapper;

typedef Handle<SegmentationWrapper<DDSegmentation::SCEPCal_MainSegmentation_k4geo>>
    SCEPCal_MainSegmentationHandle_k4geo;

class SCEPCal_MainSegmentation_k4geo : public SCEPCal_MainSegmentationHandle_k4geo {
public:
  typedef SCEPCal_MainSegmentationHandle_k4geo::Object Object;

public:
  SCEPCal_MainSegmentation_k4geo() = default;
  SCEPCal_MainSegmentation_k4geo(const SCEPCal_MainSegmentation_k4geo& e) = default;
  SCEPCal_MainSegmentation_k4geo(const Segmentation& e) : Handle<Object>(e) {}
  SCEPCal_MainSegmentation_k4geo(const Handle<Object>& e) : Handle<Object>(e) {}
  template <typename Q>
  SCEPCal_MainSegmentation_k4geo(const Handle<Q>& e) : Handle<Object>(e) {}
  SCEPCal_MainSegmentation_k4geo& operator=(const SCEPCal_MainSegmentation_k4geo& seg) = default;
  bool operator==(const SCEPCal_MainSegmentation_k4geo& seg) const { return m_element == seg.m_element; }

  inline Position position(const CellID& id) const { return Position(access()->implementation->position(id)); }

  inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID& volID) const {
    return access()->implementation->cellID(local, global, volID);
  }

  inline VolumeID setVolumeID(int System, int Phi, int Theta, int Gamma, int Epsilon, int Depth) const {
    return access()->implementation->setVolumeID(System, Phi, Theta, Gamma, Epsilon, Depth);
  }

  inline CellID setCellID(int System, int Phi, int Theta, int Gamma, int Epsilon, int Depth, bool isCheren) const {
    return access()->implementation->setCellID(System, Phi, Theta, Gamma, Epsilon, Depth, isCheren);
  }

  inline void neighbours(const CellID& cellID, std::set<CellID>& neighbours) const {
    return access()->implementation->neighbours(cellID, neighbours);
  };

  inline int System(const CellID& aCellID) const { return access()->implementation->System(aCellID); }
  inline int Phi(const CellID& aCellID) const { return access()->implementation->Phi(aCellID); }
  inline int Theta(const CellID& aCellID) const { return access()->implementation->Theta(aCellID); }
  inline int Gamma(const CellID& aCellID) const { return access()->implementation->Gamma(aCellID); }
  inline int Epsilon(const CellID& aCellID) const { return access()->implementation->Epsilon(aCellID); }
  inline int Depth(const CellID& aCellID) const { return access()->implementation->Depth(aCellID); }
  inline bool isCherenkov(const CellID& aCellID) const { return access()->implementation->isCherenkov(aCellID); }
  inline int getFirst32bits(const CellID& aCellID) const { return access()->implementation->getFirst32bits(aCellID); }
  inline int getLast32bits(const CellID& aCellID) const { return access()->implementation->getLast32bits(aCellID); }

  inline CellID convertFirst32to64(const int aId32) const {
    return access()->implementation->convertFirst32to64(aId32);
  }
  inline CellID convertLast32to64(const int aId32) const { return access()->implementation->convertLast32to64(aId32); }

  inline int System(const int& aId32) const { return access()->implementation->System(aId32); }
  inline int Phi(const int& aId32) const { return access()->implementation->Phi(aId32); }
  inline int Theta(const int& aId32) const { return access()->implementation->Theta(aId32); }
  inline int Gamma(const int& aId32) const { return access()->implementation->Gamma(aId32); }
  inline int Epsilon(const int& aId32) const { return access()->implementation->Epsilon(aId32); }
  inline int Depth(const int& aId32) const { return access()->implementation->Depth(aId32); }
};

} // namespace dd4hep
#endif

//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#ifndef SCEPCal_TimingSegmentationHandle_k4geo_h
#define SCEPCal_TimingSegmentationHandle_k4geo_h 1
#include "SCEPCal_TimingSegmentation_k4geo.h"
#include "DD4hep/Segmentations.h"
#include "DD4hep/detail/SegmentationsInterna.h"
#include "Math/Vector3D.h"

namespace dd4hep {
    
class Segmentation;
template <typename T>
class SegmentationWrapper;

typedef Handle<SegmentationWrapper<DDSegmentation::SCEPCal_TimingSegmentation_k4geo>> SCEPCal_TimingSegmentationHandle_k4geo;

class SCEPCal_TimingSegmentation_k4geo : public SCEPCal_TimingSegmentationHandle_k4geo {
public:
    typedef SCEPCal_TimingSegmentationHandle_k4geo::Object Object;

public:
    SCEPCal_TimingSegmentation_k4geo() = default;
    SCEPCal_TimingSegmentation_k4geo(const SCEPCal_TimingSegmentation_k4geo& e) = default;
    SCEPCal_TimingSegmentation_k4geo(const Segmentation& e) : Handle<Object>(e) {}
    SCEPCal_TimingSegmentation_k4geo(const Handle<Object>& e) : Handle<Object>(e) {}
    template <typename Q>
    SCEPCal_TimingSegmentation_k4geo(const Handle<Q>& e) : Handle<Object>(e) {}
    SCEPCal_TimingSegmentation_k4geo& operator=(const SCEPCal_TimingSegmentation_k4geo& seg) = default;
    bool operator==(const SCEPCal_TimingSegmentation_k4geo& seg) const { return m_element == seg.m_element; }

    inline Position position(const CellID& id) const {
        return Position(access()->implementation->position(id));
    }

    inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID& volID) const {
        return access()->implementation->cellID(local, global, volID);
    }

    inline VolumeID setVolumeID(int System, int Phi, int Eta, int Gamma) const {
        return access()->implementation->setVolumeID(System, Phi, Eta, Gamma);
    }

    inline CellID setCellID(int System, int Phi, int Eta, int Gamma) const {
        return access()->implementation->setCellID(System, Phi, Eta, Gamma);
    }

    inline int System(const CellID& aCellID) const { return access()->implementation->System(aCellID); }
    inline int Phi(const CellID& aCellID) const { return access()->implementation->Phi(aCellID); }
    inline int Eta(const CellID& aCellID) const { return access()->implementation->Eta(aCellID); }
    inline int Gamma(const CellID& aCellID) const { return access()->implementation->Gamma(aCellID); }

    inline int getFirst32bits(const CellID& aCellID) const { return access()->implementation->getFirst32bits(aCellID); }
    inline int getLast32bits(const CellID& aCellID) const { return access()->implementation->getLast32bits(aCellID); }

    inline CellID convertFirst32to64(const int aId32) const { return access()->implementation->convertFirst32to64(aId32); }
    inline CellID convertLast32to64(const int aId32) const { return access()->implementation->convertLast32to64(aId32); }

    inline int System(const int& aId32) const { return access()->implementation->System(aId32); }
    inline int Phi(const int& aId32) const { return access()->implementation->Phi(aId32); }
    inline int Eta(const int& aId32) const { return access()->implementation->Eta(aId32); }
    inline int Gamma(const int& aId32) const { return access()->implementation->Gamma(aId32); }
};

}
#endif

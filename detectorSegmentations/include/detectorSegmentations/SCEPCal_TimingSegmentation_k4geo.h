//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#ifndef SCEPCal_TimingSegmentation_k4geo_h
#define SCEPCal_TimingSegmentation_k4geo_h 1
#include "DD4hep/DetFactoryHelper.h"
#include "DDSegmentation/Segmentation.h"
#include "Math/Vector3D.h"
#include <cmath>
#include <vector>

namespace dd4hep {
namespace DDSegmentation {

  class SCEPCal_TimingSegmentation_k4geo : public Segmentation {
  public:
    SCEPCal_TimingSegmentation_k4geo(const std::string& aCellEncoding);
    SCEPCal_TimingSegmentation_k4geo(const BitFieldCoder* decoder);
    virtual ~SCEPCal_TimingSegmentation_k4geo() = default;

    virtual Vector3D position(const CellID& aCellID) const override;

    CellID cellID(const Vector3D& /*localPosition*/, const Vector3D& /*globalPosition*/,
                  const VolumeID& vID) const override {
      return setCellID(System(vID), Phi(vID), Theta(vID), Gamma(vID));
    }

    VolumeID setVolumeID(int System, int Phi, int Theta, int Gamma) const {
      VolumeID SystemId = static_cast<VolumeID>(System);
      VolumeID PhiId = static_cast<VolumeID>(Phi);
      VolumeID ThetaId = static_cast<VolumeID>(Theta);
      VolumeID GammaId = static_cast<VolumeID>(Gamma);
      VolumeID vID = 0;
      decoder()->set(vID, m_systemIndex, SystemId);
      decoder()->set(vID, m_phiIndex, PhiId);
      decoder()->set(vID, m_thetaIndex, ThetaId);
      decoder()->set(vID, m_gammaIndex, GammaId);
      return vID;
    }

    CellID setCellID(int System, int Phi, int Theta, int Gamma) const {
      VolumeID SystemId = static_cast<VolumeID>(System);
      VolumeID PhiId = static_cast<VolumeID>(Phi);
      VolumeID ThetaId = static_cast<VolumeID>(Theta);
      VolumeID GammaId = static_cast<VolumeID>(Gamma);
      VolumeID vID = 0;
      decoder()->set(vID, m_systemIndex, SystemId);
      decoder()->set(vID, m_phiIndex, PhiId);
      decoder()->set(vID, m_thetaIndex, ThetaId);
      decoder()->set(vID, m_gammaIndex, GammaId);
      return vID;
    }

    int System(const CellID aCellID) const {
      VolumeID System = static_cast<VolumeID>(decoder()->get(aCellID, m_systemIndex));
      return static_cast<int>(System);
    }

    int Phi(const CellID aCellID) const {
      VolumeID Phi = static_cast<VolumeID>(decoder()->get(aCellID, m_phiIndex));
      return static_cast<int>(Phi);
    }

    int Theta(const CellID aCellID) const {
      VolumeID Theta = static_cast<VolumeID>(decoder()->get(aCellID, m_thetaIndex));
      return static_cast<int>(Theta);
    }

    int Gamma(const CellID aCellID) const {
      VolumeID Gamma = static_cast<VolumeID>(decoder()->get(aCellID, m_gammaIndex));
      return static_cast<int>(Gamma);
    }

    int getFirst32bits(const CellID aCellID) const { return static_cast<int>(aCellID); }
    int getLast32bits(const CellID aCellID) const {
      CellID aId64 = aCellID >> sizeof(int) * CHAR_BIT;
      int aId32 = static_cast<int>(aId64);
      return aId32;
    }

    CellID convertFirst32to64(const int aId32) const { return static_cast<CellID>(aId32); }
    CellID convertLast32to64(const int aId32) const {
      CellID aId64 = static_cast<CellID>(aId32);
      aId64 <<= sizeof(int) * CHAR_BIT;
      return aId64;
    }

    int System(const int aId32) const { return System(convertFirst32to64(aId32)); }
    int Phi(const int aId32) const { return Phi(convertFirst32to64(aId32)); }
    int Theta(const int aId32) const { return Theta(convertFirst32to64(aId32)); }
    int Gamma(const int aId32) const { return Gamma(convertFirst32to64(aId32)); }

    inline void savePosition(int volID_32, Vector3D pos) { m_positionOf.emplace(volID_32, pos); }

  private:
    /// Initialization common to all ctors.
    void commonSetup();

    std::string m_systemId;
    std::string m_phiId;
    std::string m_thetaId;
    std::string m_gammaId;

    int m_systemIndex = -1;
    int m_phiIndex = -1;
    int m_thetaIndex = -1;
    int m_gammaIndex = -1;

    std::unordered_map<int, Vector3D> m_positionOf;
  };
} // namespace DDSegmentation
} // namespace dd4hep

#endif

//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#ifndef SCEPCal_MainSegmentation_k4geo_h
#define SCEPCal_MainSegmentation_k4geo_h 1
#include "DD4hep/DetFactoryHelper.h"
#include "DDSegmentation/Segmentation.h"
#include "Math/Vector3D.h"
#include <climits> // For CHAR_BIT
#include <cmath>
#include <vector>

namespace dd4hep {
namespace DDSegmentation {

  class SCEPCal_MainSegmentation_k4geo : public Segmentation {
  public:
    SCEPCal_MainSegmentation_k4geo(const std::string& aCellEncoding);
    SCEPCal_MainSegmentation_k4geo(const BitFieldCoder* decoder);
    virtual ~SCEPCal_MainSegmentation_k4geo() = default;

    virtual Vector3D position(const CellID& aCellID) const override;

    CellID cellID(const Vector3D& /*localPosition*/, const Vector3D& /*globalPosition*/,
                  const VolumeID& vID) const override {
      // crystal cell ID returns scintillation channel by default
      return setCellID(System(vID), Phi(vID), Theta(vID), Gamma(vID), Epsilon(vID), Depth(vID), false);
    }

    VolumeID setVolumeID(int System, int Phi, int Theta, int Gamma, int Epsilon, int Depth) const {
      VolumeID SystemId = static_cast<VolumeID>(System);
      VolumeID PhiId = static_cast<VolumeID>(Phi);
      VolumeID ThetaId = static_cast<VolumeID>(Theta);
      VolumeID GammaId = static_cast<VolumeID>(Gamma);
      VolumeID EpsilonId = static_cast<VolumeID>(Epsilon);
      VolumeID DepthId = static_cast<VolumeID>(Depth);
      VolumeID vID = 0;
      decoder()->set(vID, m_systemIndex, SystemId);
      decoder()->set(vID, m_phiIndex, PhiId);
      decoder()->set(vID, m_thetaIndex, ThetaId);
      decoder()->set(vID, m_gammaIndex, GammaId);
      decoder()->set(vID, m_epsilonIndex, EpsilonId);
      decoder()->set(vID, m_depthIndex, DepthId);
      return vID;
    }

    CellID setCellID(int System, int Phi, int Theta, int Gamma, int Epsilon, int Depth, bool isCheren) const {
      VolumeID SystemId = static_cast<VolumeID>(System);
      VolumeID PhiId = static_cast<VolumeID>(Phi);
      VolumeID ThetaId = static_cast<VolumeID>(Theta);
      VolumeID GammaId = static_cast<VolumeID>(Gamma);
      VolumeID EpsilonId = static_cast<VolumeID>(Epsilon);
      VolumeID DepthId = static_cast<VolumeID>(Depth);
      VolumeID vID = 0;
      decoder()->set(vID, m_systemIndex, SystemId);
      decoder()->set(vID, m_phiIndex, PhiId);
      decoder()->set(vID, m_thetaIndex, ThetaId);
      decoder()->set(vID, m_gammaIndex, GammaId);
      decoder()->set(vID, m_epsilonIndex, EpsilonId);
      decoder()->set(vID, m_depthIndex, DepthId);
      decoder()->set(vID, m_isCherenIndex, static_cast<VolumeID>(isCheren));

      return vID;
    }

    CellID setCellID(CellID cID, bool isCheren) const {
      return setCellID(System(cID), Phi(cID), Theta(cID), Gamma(cID), Epsilon(cID), Depth(cID), isCheren);
    }

    // function that retrieves the std::set of neighbouring cell IDs
    // overrides the DDSegmentation::Segmentation::neighbours method
    virtual void neighbours(const CellID& cellID, std::set<CellID>& neighbours) const override;

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

    int Epsilon(const CellID aCellID) const {
      VolumeID Epsilon = static_cast<VolumeID>(decoder()->get(aCellID, m_epsilonIndex));
      return static_cast<int>(Epsilon);
    }

    int Depth(const CellID aCellID) const {
      VolumeID Depth = static_cast<VolumeID>(decoder()->get(aCellID, m_depthIndex));
      return static_cast<int>(Depth);
    }

    bool isCherenkov(const CellID aCellID) const {
      VolumeID isCheren = static_cast<VolumeID>(decoder()->get(aCellID, m_isCherenIndex));
      return static_cast<bool>(isCheren);
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
    int Epsilon(const int aId32) const { return Epsilon(convertFirst32to64(aId32)); }
    int Depth(const int aId32) const { return Depth(convertFirst32to64(aId32)); }

    inline void savePosition(int volID_32, Vector3D pos) { m_positionOf.emplace(volID_32, pos); }

    // setters for geometry parameters
    void setDetIdBarrel(int detId) { m_detId_barrel_ = detId; }
    void setDetIdEndcap(int detId) { m_detId_endcap_ = detId; }
    void setIThetaBarrelStart(int iThetaStart) { m_iTheta_barrel_start_ = iThetaStart; }
    void setIThetaBarrelEnd(int iThetaEnd) { m_iTheta_barrel_end_ = iThetaEnd; }
    void setNPhi(int nPhi) { m_nPhi_ = nPhi; }
    void setNGamma(int nGamma) { m_nGamma_ = nGamma; }

  private:
    /// Initialization common to all ctors.
    void commonSetup();

    std::string m_systemId;
    std::string m_phiId;
    std::string m_thetaId;
    std::string m_gammaId;
    std::string m_epsilonId;
    std::string m_depthId;
    std::string m_isCherenId;

    int m_systemIndex = -1;
    int m_phiIndex = -1;
    int m_thetaIndex = -1;
    int m_gammaIndex = -1;
    int m_epsilonIndex = -1;
    int m_depthIndex = -1;
    int m_isCherenIndex = -1;

    // geometry parameters
    int m_detId_barrel_ = 0;         // system id for barrel
    int m_detId_endcap_ = 0;         // system id for endcap
    int m_iTheta_barrel_start_ = -1; // theta index start for barrel
    int m_iTheta_barrel_end_ = -1;   // theta index end for barrel
    int m_nPhi_ = -1;                // number of phi segments
    int m_nGamma_ = -1;              // number of gamma segments
    int m_nEpsilon_ = 0;             // note: always zero in this segmentation

    std::unordered_map<int, Vector3D> m_positionOf;
  };
} // namespace DDSegmentation
} // namespace dd4hep

#endif

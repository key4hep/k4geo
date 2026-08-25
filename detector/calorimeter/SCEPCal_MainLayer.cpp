//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorTools.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Volumes.h"
#include "DDRec/DetectorData.h"
#include "TGeoTrd2.h"
#include "detectorSegmentations/SCEPCal_MainSegmentationHandle_k4geo.h"
#include <bitset>
#include <unordered_map>

#include "XML/Utilities.h"

using dd4hep::Position;
using dd4hep::RotationZYX;
using dd4hep::Transform3D;
using ROOT::Math::RotationY;
using ROOT::Math::RotationZ;
using ROOT::Math::XYZVector;
using namespace dd4hep;

// -----------------------------------------------------------------------
// Free helper functions
// -----------------------------------------------------------------------

static void bilinearInterpolateTower(const double* vf, double u, double v, std::array<double, 4>& result) {
  result[0] = (1 - u) * (1 - v) * vf[0] + u * (1 - v) * vf[2] + u * v * vf[4] + (1 - u) * v * vf[6];
  result[1] = (1 - u) * (1 - v) * vf[1] + u * (1 - v) * vf[3] + u * v * vf[5] + (1 - u) * v * vf[7];
  result[2] = (1 - u) * (1 - v) * vf[8] + u * (1 - v) * vf[10] + u * v * vf[12] + (1 - u) * v * vf[14];
  result[3] = (1 - u) * (1 - v) * vf[9] + u * (1 - v) * vf[11] + u * v * vf[13] + (1 - u) * v * vf[15];
}

static void getSingleCrystalVertices(int i, int j, int xtalDiv, const double* vertices, bool reflected,
                                     std::array<double, 16>& result) {
  double u0 = double(i) / double(xtalDiv);
  double u1 = double(i + 1) / double(xtalDiv);
  double v0 = double(j) / double(xtalDiv);
  double v1 = double(j + 1) / double(xtalDiv);

  std::array<double, 4> P00, P10, P11, P01;
  bilinearInterpolateTower(vertices, u0, v0, P00);
  bilinearInterpolateTower(vertices, u1, v0, P10);
  bilinearInterpolateTower(vertices, u1, v1, P11);
  bilinearInterpolateTower(vertices, u0, v1, P01);

  if (reflected)
    result = {P00[2], P00[3], P10[2], P10[3], P11[2], P11[3], P01[2], P01[3],
              P00[0], P00[1], P10[0], P10[1], P11[0], P11[1], P01[0], P01[1]};
  else
    result = {P00[0], P00[1], P10[0], P10[1], P11[0], P11[1], P01[0], P01[1],
              P00[2], P00[3], P10[2], P10[3], P11[2], P11[3], P01[2], P01[3]};
}

static void getSingleCrystalCenter(const std::array<double, 16>& vSub, std::array<double, 2>& result) {
  result[0] = 0.125 * (vSub[0] + vSub[2] + vSub[4] + vSub[6] + vSub[8] + vSub[10] + vSub[12] + vSub[14]);
  result[1] = 0.125 * (vSub[1] + vSub[3] + vSub[5] + vSub[7] + vSub[9] + vSub[11] + vSub[13] + vSub[15]);
}

static void placeCrystal(dd4hep::Detector& theDetector, dd4hep::SensitiveDetector& sens,
                         dd4hep::DDSegmentation::SCEPCal_MainSegmentation_k4geo* segmentation,
                         dd4hep::DetElement& ScepcalDetElement, bool USE_OPTICAL_SURFACES,
                         dd4hep::OpticalSurface& PbWO4_to_ESR, const std::string& volName, double dz,
                         const std::array<double, 16>& vertices, const std::string& materialStr,
                         const std::string& visStr, const dd4hep::Transform3D& transform,
                         const dd4hep::Volume& assemblyVol, int nSystem, int nTheta, int nGamma, int nEpsilon,
                         int nDepth, const XYZVector& posGlobal, int phi_start, int phi_end,
                         const std::unordered_map<int, RotationZ>& phi_rotations) {
  dd4hep::EightPointSolid theShape(dz, vertices.data());
  dd4hep::Volume theVolume(volName, theShape, theDetector.material(materialStr));
  theVolume.setVisAttributes(theDetector, visStr);

  if (USE_OPTICAL_SURFACES) {
    dd4hep::SkinSurface(theDetector, ScepcalDetElement,
                        volName + "Surface_" + std::to_string(nTheta) + "_" + std::to_string(nGamma) + "_" +
                            std::to_string(nDepth),
                        PbWO4_to_ESR, theVolume);
  }

  theVolume.setSensitiveDetector(sens);
  auto volID_0 = segmentation->setVolumeID(nSystem, 0, nTheta, nGamma, nEpsilon, nDepth);
  int volID_0_32 = segmentation->getFirst32bits(volID_0);
  dd4hep::PlacedVolume thePlacedVol = assemblyVol.placeVolume(theVolume, volID_0_32, transform);

  thePlacedVol.addPhysVolID("system", nSystem);
  thePlacedVol.addPhysVolID("theta", nTheta);
  thePlacedVol.addPhysVolID("gamma", nGamma);
  thePlacedVol.addPhysVolID("epsilon", nEpsilon);
  thePlacedVol.addPhysVolID("depth", nDepth);

  for (int iPhi = phi_start; iPhi < phi_end; iPhi++) {
    XYZVector posGlobalPhi = phi_rotations.at(iPhi) * posGlobal;
    auto volID = segmentation->setVolumeID(nSystem, iPhi, nTheta, nGamma, nEpsilon, nDepth);
    int volID_32 = segmentation->getFirst32bits(volID);
    segmentation->savePosition(volID_32, posGlobalPhi);
  }
}

// -----------------------------------------------------------------------
// Shared geometry parameters
// -----------------------------------------------------------------------

struct SCEPCalParams {
  // Basic dimensions
  double BARREL_HALF_Z, BARREL_INNER_R;
  double XTAL_LEN_F, XTAL_LEN_R;
  int XTAL_DIV_F, XTAL_DIV_R;
  double XTAL_TH_WIDTH;
  double PROJ_OFFSET_R, PROJ_OFFSET_X;
  int PHI_SEGMENTS;
  double D_PHI_GLOBAL;

  // Barrel
  int N_THETA_BARREL, N_GAMMA_BARREL;
  double D_THETA_BARREL, D_GAMMA_BARREL;
  double THETA_SIZE_BARREL, THETA_SIZE_ENDCAP;
  int BARREL_SYSTEM_NO, BARREL_PHI_START, BARREL_PHI_END;
  bool CONSTRUCT_BARREL;

  // Endcap
  int N_THETA_ENDCAP, ENDCAP_THETA_START;
  double D_THETA_ENDCAP;
  int ENDCAP_SYSTEM_NO, ENDCAP_PHI_START, ENDCAP_PHI_END;
  bool CONSTRUCT_ENDCAP;

  // Partial visualisation
  int PHI_LOAD_START, PHI_LOAD_END;
  int THETA_LOAD_START, THETA_LOAD_END;
  int GAMMA_LOAD_START, GAMMA_LOAD_END;

  // Crystal material and vis strings (extracted from XML at setup time)
  std::string crystalFMaterial, crystalFVis;
  std::string crystalRMaterial, crystalRVis;

  // Assembly vis strings
  std::string barrelAssemblyPhiVis, barrelAssemblyThetaVis;
  std::string endcapAssemblyPhiVis, endcapAssemblyThetaVis;

  // Crystal placement transforms
  Transform3D trans_dispFsub, trans_dispRsub, trans_dispFsub_1, trans_dispRsub_1;
  XYZVector DISP_PROJ_R;

  // Barrel phi-slice envelope polyhedra
  std::vector<double> zBarrelPolyhedra, rminBarrelPolyhedra, rmaxBarrelPolyhedra;

  // Endcap (+z) phi-slice envelope polyhedra
  std::vector<double> zEndcapPolyhedra, rminEndcapPolyhedra, rmaxEndcapPolyhedra;
  // Endcap (-z) phi-slice envelope polyhedra
  std::vector<double> zEndcapPolyhedra_1, rminEndcapPolyhedra_1, rmaxEndcapPolyhedra_1;
};

// -----------------------------------------------------------------------
// Barrel builder
// -----------------------------------------------------------------------

static int buildBarrel(dd4hep::Detector& theDetector, dd4hep::SensitiveDetector& sens,
                       dd4hep::DDSegmentation::SCEPCal_MainSegmentation_k4geo* segmentation,
                       dd4hep::DetElement& ScepcalDetElement, bool USE_OPTICAL_SURFACES,
                       dd4hep::OpticalSurface& PbWO4_to_ESR, const SCEPCalParams& p,
                       const std::unordered_map<int, RotationZ>& phi_barrel_rotations,
                       dd4hep::Volume& barrelGlobalAssemblyVol) {
  int numCrystals = 0;

  const double tan_d_theta_barrel2 = tan(p.D_THETA_BARREL / 2.);
  const double tan_d_phi_global2 = tan(p.D_PHI_GLOBAL / 2.);
  const double tan_m_d_phi_global2 = -tan_d_phi_global2;
  const double tan_d_gamma_barrel2 = tan(p.D_GAMMA_BARREL / 2);

  std::vector<double> tan_gamma_m_dgb2s, tan_gamma_p_dgb2s;
  tan_gamma_m_dgb2s.reserve(p.N_GAMMA_BARREL);
  tan_gamma_p_dgb2s.reserve(p.N_GAMMA_BARREL);
  for (int nGamma = 0; nGamma < p.N_GAMMA_BARREL; nGamma++) {
    double gamma = -p.D_PHI_GLOBAL / 2 + p.D_GAMMA_BARREL / 2 + p.D_GAMMA_BARREL * nGamma;
    tan_gamma_m_dgb2s.push_back(tan(gamma - p.D_GAMMA_BARREL / 2));
    tan_gamma_p_dgb2s.push_back(tan(gamma + p.D_GAMMA_BARREL / 2));
  }

  dd4hep::Polyhedra barrelPhiAssemblyShape(1, -p.D_PHI_GLOBAL / 2, p.D_PHI_GLOBAL, p.zBarrelPolyhedra,
                                           p.rminBarrelPolyhedra, p.rmaxBarrelPolyhedra);
  dd4hep::Volume barrelPhiAssemblyVolume("barrelPhiAssembly", barrelPhiAssemblyShape, theDetector.material("Vacuum"));
  barrelPhiAssemblyVolume.setVisAttributes(theDetector, p.barrelAssemblyPhiVis);

  for (int iTheta = 0; iTheta < p.N_THETA_BARREL; iTheta++) {
    if ((iTheta < p.THETA_LOAD_START) || (iTheta > p.THETA_LOAD_END))
      continue;

    double thC = p.THETA_SIZE_ENDCAP + p.D_THETA_BARREL / 2 + (iTheta * p.D_THETA_BARREL);
    double sin_thC = sin(thC);
    double cos_thC = cos(thC);
    RotationY rotYthGlobal(thC);

    double r0e = (p.BARREL_INNER_R + p.PROJ_OFFSET_R) / sin_thC;
    double r1e = r0e + p.XTAL_LEN_F;
    double r2e = r1e + p.XTAL_LEN_R;
    double y0e = r0e * tan_d_theta_barrel2;
    double y1e = r1e * tan_d_theta_barrel2;
    double y2e = r2e * tan_d_theta_barrel2;

    double x0y0 = r0e * sin_thC - y0e * cos_thC - p.PROJ_OFFSET_R;
    double x1y0 = r0e * sin_thC + y0e * cos_thC - p.PROJ_OFFSET_R;
    double x0y1 = r1e * sin_thC - y1e * cos_thC - p.PROJ_OFFSET_R;
    double x1y1 = r1e * sin_thC + y1e * cos_thC - p.PROJ_OFFSET_R;
    double x0y2 = r2e * sin_thC - y2e * cos_thC - p.PROJ_OFFSET_R;
    double x1y2 = r2e * sin_thC + y2e * cos_thC - p.PROJ_OFFSET_R;

    double x0y0l_E = x0y0 * tan_m_d_phi_global2;
    double x0y0r_E = x0y0 * tan_d_phi_global2;
    double x1y0l_E = x1y0 * tan_m_d_phi_global2;
    double x1y0r_E = x1y0 * tan_d_phi_global2;
    double x0y2l_E = x0y2 * tan_m_d_phi_global2;
    double x0y2r_E = x0y2 * tan_d_phi_global2;
    double x1y2l_E = x1y2 * tan_m_d_phi_global2;
    double x1y2r_E = x1y2 * tan_d_phi_global2;

    double verticesE[] = {x0y0r_E, y0e, x1y0r_E, -y0e, x1y0l_E, -y0e, x0y0l_E, y0e,
                          x0y2r_E, y2e, x1y2r_E, -y2e, x1y2l_E, -y2e, x0y2l_E, y2e};

    double rE = r0e + (p.XTAL_LEN_F + p.XTAL_LEN_R) / 2.;
    RotationZYX rotE(M_PI / 2, thC, 0);
    Position dispE(rE * sin_thC - p.PROJ_OFFSET_R, 0, rE * cos_thC);

    dd4hep::EightPointSolid barrelThetaAssemblyShape((p.XTAL_LEN_F + p.XTAL_LEN_R) / 2, verticesE);
    dd4hep::Volume barrelThetaAssemblyVolume("barrelThetaAssembly", barrelThetaAssemblyShape,
                                             theDetector.material("Vacuum"));
    barrelThetaAssemblyVolume.setVisAttributes(theDetector, p.barrelAssemblyThetaVis);
    barrelPhiAssemblyVolume.placeVolume(barrelThetaAssemblyVolume, Transform3D(rotE, dispE));

    for (int nGamma = 0; nGamma < p.N_GAMMA_BARREL; nGamma++) {
      if ((nGamma < p.GAMMA_LOAD_START) || (nGamma > p.GAMMA_LOAD_END))
        continue;

      double projOffsetXmax = std::min(r0e * tan_d_gamma_barrel2, p.PROJ_OFFSET_X);
      double r1_x_gamma_shift = projOffsetXmax * (r1e - r0e) / r0e;
      double r2_x_gamma_shift = projOffsetXmax * (r2e - r0e) / r0e;

      int left = 0, right = 0;
      if (p.N_GAMMA_BARREL % 2 == 0) {
        if (nGamma < (p.N_GAMMA_BARREL / 2 - 1)) {
          left = 1;
          right = 1;
        }
        if (nGamma == (p.N_GAMMA_BARREL / 2 - 1)) {
          left = 1;
          right = 0;
        }
        if (nGamma == p.N_GAMMA_BARREL / 2) {
          left = 0;
          right = -1;
        }
        if (nGamma > p.N_GAMMA_BARREL / 2) {
          left = -1;
          right = -1;
        }
      } else {
        if (nGamma < p.N_GAMMA_BARREL / 2) {
          left = 1;
          right = 1;
        }
        if (nGamma == p.N_GAMMA_BARREL / 2) {
          left = 1;
          right = -1;
        }
        if (nGamma > p.N_GAMMA_BARREL / 2) {
          left = -1;
          right = -1;
        }
      }
      if (nGamma == 0)
        left = 0;
      if (nGamma == p.N_GAMMA_BARREL - 1)
        right = 0;

      double x0y0l = x0y0 * tan_gamma_m_dgb2s[nGamma];
      double x0y0r = x0y0 * tan_gamma_p_dgb2s[nGamma];
      double x1y0l = x1y0 * tan_gamma_m_dgb2s[nGamma];
      double x1y0r = x1y0 * tan_gamma_p_dgb2s[nGamma];

      double x0y1l = x0y1 * tan_gamma_m_dgb2s[nGamma] + left * r1_x_gamma_shift;
      double x0y1r = x0y1 * tan_gamma_p_dgb2s[nGamma] + right * r1_x_gamma_shift;
      double x1y1l = x1y1 * tan_gamma_m_dgb2s[nGamma] + left * r1_x_gamma_shift;
      double x1y1r = x1y1 * tan_gamma_p_dgb2s[nGamma] + right * r1_x_gamma_shift;

      double x0y2l = x0y2 * tan_gamma_m_dgb2s[nGamma] + left * r2_x_gamma_shift;
      double x0y2r = x0y2 * tan_gamma_p_dgb2s[nGamma] + right * r2_x_gamma_shift;
      double x1y2l = x1y2 * tan_gamma_m_dgb2s[nGamma] + left * r2_x_gamma_shift;
      double x1y2r = x1y2 * tan_gamma_p_dgb2s[nGamma] + right * r2_x_gamma_shift;

      double verticesF[] = {x0y0r, y0e, x1y0r, -y0e, x1y0l, -y0e, x0y0l, y0e,
                            x0y1r, y1e, x1y1r, -y1e, x1y1l, -y1e, x0y1l, y1e};
      double verticesR[] = {x0y1r, y1e, x1y1r, -y1e, x1y1l, -y1e, x0y1l, y1e,
                            x0y2r, y2e, x1y2r, -y2e, x1y2l, -y2e, x0y2l, y2e};

      std::array<double, 16> vsub;
      std::array<double, 2> center;

      double rGlobal = r0e + p.XTAL_LEN_F / 2;
      for (int i = 0; i < p.XTAL_DIV_F; i++) {
        for (int j = 0; j < p.XTAL_DIV_F; j++) {
          int nEpsilon = i * p.XTAL_DIV_F + j;
          getSingleCrystalVertices(i, j, p.XTAL_DIV_F, verticesF, false, vsub);
          getSingleCrystalCenter(vsub, center);
          XYZVector dispGlobal(-center[1], center[0], rGlobal);
          XYZVector posGlobal = (rotYthGlobal * dispGlobal + p.DISP_PROJ_R);
          placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, PbWO4_to_ESR,
                       "BarrelCrystalF", p.XTAL_LEN_F / 2, vsub, p.crystalFMaterial, p.crystalFVis, p.trans_dispFsub,
                       barrelThetaAssemblyVolume, p.BARREL_SYSTEM_NO, p.N_THETA_ENDCAP + iTheta, nGamma, nEpsilon, 0,
                       posGlobal, p.BARREL_PHI_START, p.BARREL_PHI_END, phi_barrel_rotations);
          numCrystals++;
        }
      }

      rGlobal = r0e + p.XTAL_LEN_F + p.XTAL_LEN_R / 2;
      for (int i = 0; i < p.XTAL_DIV_R; i++) {
        for (int j = 0; j < p.XTAL_DIV_R; j++) {
          int nEpsilon = i * p.XTAL_DIV_R + j;
          getSingleCrystalVertices(i, j, p.XTAL_DIV_R, verticesR, false, vsub);
          getSingleCrystalCenter(vsub, center);
          XYZVector dispGlobal(-center[1], center[0], rGlobal);
          XYZVector posGlobal = (rotYthGlobal * dispGlobal + p.DISP_PROJ_R);
          placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, PbWO4_to_ESR,
                       "BarrelCrystalR", p.XTAL_LEN_R / 2, vsub, p.crystalRMaterial, p.crystalRVis, p.trans_dispRsub,
                       barrelThetaAssemblyVolume, p.BARREL_SYSTEM_NO, p.N_THETA_ENDCAP + iTheta, nGamma, nEpsilon, 1,
                       posGlobal, p.BARREL_PHI_START, p.BARREL_PHI_END, phi_barrel_rotations);
          numCrystals++;
        }
      }
    }
  }

  unsigned int nPhiSlicePlaced = 0;
  for (int iPhi = p.BARREL_PHI_START; iPhi < p.BARREL_PHI_END; iPhi++) {
    if (p.PHI_LOAD_START <= p.PHI_LOAD_END) {
      if ((iPhi < p.PHI_LOAD_START) || (iPhi > p.PHI_LOAD_END))
        continue;
    } else {
      if ((iPhi < p.PHI_LOAD_START) && (iPhi > p.PHI_LOAD_END))
        continue;
    }
    nPhiSlicePlaced++;
    double phiGlobal = iPhi * p.D_PHI_GLOBAL;
    RotationZ rotZphiGlobal(phiGlobal);
    auto pv = barrelGlobalAssemblyVol.placeVolume(barrelPhiAssemblyVolume, Transform3D(rotZphiGlobal));
    pv.addPhysVolID("phi", iPhi);
  }
  return numCrystals * static_cast<int>(nPhiSlicePlaced);
}

// -----------------------------------------------------------------------
// Endcap builder
// -----------------------------------------------------------------------

static int buildEndcap(dd4hep::Detector& theDetector, dd4hep::SensitiveDetector& sens,
                       dd4hep::DDSegmentation::SCEPCal_MainSegmentation_k4geo* segmentation,
                       dd4hep::DetElement& ScepcalDetElement, bool USE_OPTICAL_SURFACES,
                       dd4hep::OpticalSurface& PbWO4_to_ESR, const SCEPCalParams& p,
                       const std::unordered_map<int, RotationZ>& phi_endcap_rotations,
                       dd4hep::Volume& endcapGlobalAssemblyVol, dd4hep::Volume& endcapGlobalAssemblyVol_1) {
  int numCrystals = 0;

  const double tan_d_theta_endcap2 = tan(p.D_THETA_ENDCAP / 2.);
  const double tan_d_phi_global2 = tan(p.D_PHI_GLOBAL / 2.);
  const double tan_m_d_phi_global2 = -tan_d_phi_global2;

  dd4hep::Polyhedra endcapPhiAssemblyShape(1, -p.D_PHI_GLOBAL / 2, p.D_PHI_GLOBAL, p.zEndcapPolyhedra,
                                           p.rminEndcapPolyhedra, p.rmaxEndcapPolyhedra);
  dd4hep::Polyhedra endcapPhiAssemblyShape_1(1, -p.D_PHI_GLOBAL / 2, p.D_PHI_GLOBAL, p.zEndcapPolyhedra_1,
                                             p.rminEndcapPolyhedra_1, p.rmaxEndcapPolyhedra_1);

  dd4hep::Volume endcapPhiAssemblyVolume("endcapPhiVol", endcapPhiAssemblyShape, theDetector.material("Vacuum"));
  endcapPhiAssemblyVolume.setVisAttributes(theDetector, p.endcapAssemblyPhiVis);
  dd4hep::Volume endcapPhiAssemblyVolume_1("endcapPhiVol_1", endcapPhiAssemblyShape_1, theDetector.material("Vacuum"));
  endcapPhiAssemblyVolume_1.setVisAttributes(theDetector, p.endcapAssemblyPhiVis);

  for (int iTheta = p.ENDCAP_THETA_START; iTheta < p.N_THETA_ENDCAP; iTheta++) {
    double thC = p.D_THETA_ENDCAP / 2 + iTheta * p.D_THETA_ENDCAP;
    double sin_thC = sin(thC);
    double cos_thC = cos(thC);
    RotationY rotYthGlobal(thC);
    RotationY rotYthGlobal_1(-thC);

    double RinEndcap = p.BARREL_HALF_Z * tan(thC);
    int nGammaEndcap = std::max(int(2 * M_PI * RinEndcap / (p.PHI_SEGMENTS * p.XTAL_TH_WIDTH)), 1);
    double dGammaEndcap = p.D_PHI_GLOBAL / nGammaEndcap;
    double tan_dGamma_Endcap2 = tan(dGammaEndcap / 2);

    double r0e = RinEndcap / sin_thC;
    double r1e = r0e + p.XTAL_LEN_F;
    double r2e = r1e + p.XTAL_LEN_R;
    double y0e = r0e * tan_d_theta_endcap2;
    double y1e = r1e * tan_d_theta_endcap2;
    double y2e = r2e * tan_d_theta_endcap2;

    double x0y0 = r0e * sin_thC - y0e * cos_thC - p.PROJ_OFFSET_R;
    double x1y0 = r0e * sin_thC + y0e * cos_thC - p.PROJ_OFFSET_R;
    double x0y1 = r1e * sin_thC - y1e * cos_thC - p.PROJ_OFFSET_R;
    double x1y1 = r1e * sin_thC + y1e * cos_thC - p.PROJ_OFFSET_R;
    double x0y2 = r2e * sin_thC - y2e * cos_thC - p.PROJ_OFFSET_R;
    double x1y2 = r2e * sin_thC + y2e * cos_thC - p.PROJ_OFFSET_R;

    double x0y0l_E = x0y0 * tan_m_d_phi_global2;
    double x0y0r_E = x0y0 * tan_d_phi_global2;
    double x1y0l_E = x1y0 * tan_m_d_phi_global2;
    double x1y0r_E = x1y0 * tan_d_phi_global2;
    double x0y2l_E = x0y2 * tan_m_d_phi_global2;
    double x0y2r_E = x0y2 * tan_d_phi_global2;
    double x1y2l_E = x1y2 * tan_m_d_phi_global2;
    double x1y2r_E = x1y2 * tan_d_phi_global2;

    double verticesE[] = {x0y0r_E, y0e, x1y0r_E, -y0e, x1y0l_E, -y0e, x0y0l_E, y0e,
                          x0y2r_E, y2e, x1y2r_E, -y2e, x1y2l_E, -y2e, x0y2l_E, y2e};
    double verticesE_1[] = {x0y2r_E, y2e, x1y2r_E, -y2e, x1y2l_E, -y2e, x0y2l_E, y2e,
                            x0y0r_E, y0e, x1y0r_E, -y0e, x1y0l_E, -y0e, x0y0l_E, y0e};

    double rE = r0e + (p.XTAL_LEN_F + p.XTAL_LEN_R) / 2.;
    RotationZYX rotE(M_PI / 2, thC, 0);
    RotationZYX rotE_1(M_PI / 2, -thC, 0);
    Position dispE(rE * sin_thC - p.PROJ_OFFSET_R, 0, rE * cos_thC);
    Position dispE_1(rE * sin_thC - p.PROJ_OFFSET_R, 0, -rE * cos_thC);

    dd4hep::EightPointSolid endcapThetaAssemblyShape((p.XTAL_LEN_F + p.XTAL_LEN_R) / 2, verticesE);
    dd4hep::EightPointSolid endcapThetaAssemblyShape_1((p.XTAL_LEN_F + p.XTAL_LEN_R) / 2, verticesE_1);

    dd4hep::Volume endcapThetaAssemblyVolume("endcapThetaAssembly", endcapThetaAssemblyShape,
                                             theDetector.material("Vacuum"));
    endcapThetaAssemblyVolume.setVisAttributes(theDetector, p.endcapAssemblyThetaVis);
    endcapPhiAssemblyVolume.placeVolume(endcapThetaAssemblyVolume, Transform3D(rotE, dispE));

    dd4hep::Volume endcapThetaAssemblyVolume_1("endcapThetaAssembly_1", endcapThetaAssemblyShape_1,
                                               theDetector.material("Vacuum"));
    endcapThetaAssemblyVolume_1.setVisAttributes(theDetector, p.endcapAssemblyThetaVis);
    endcapPhiAssemblyVolume_1.placeVolume(endcapThetaAssemblyVolume_1, Transform3D(rotE_1, dispE_1));

    for (int nGamma = 0; nGamma < nGammaEndcap; nGamma++) {
      double gamma = -p.D_PHI_GLOBAL / 2 + dGammaEndcap / 2 + dGammaEndcap * nGamma;
      double tan_g_m_dg2 = tan(gamma - dGammaEndcap / 2);
      double tan_g_p_dg2 = tan(gamma + dGammaEndcap / 2);
      double projOffsetXmax = std::min(r0e * tan_dGamma_Endcap2, p.PROJ_OFFSET_X);

      double r1_x_gamma_shift = projOffsetXmax * (r1e - r0e) / r0e;
      double r2_x_gamma_shift = projOffsetXmax * (r2e - r0e) / r0e;

      int left = 0, right = 0;
      if (nGammaEndcap % 2 == 0) {
        if (nGamma < (nGammaEndcap / 2 - 1)) {
          left = 1;
          right = 1;
        }
        if (nGamma == (nGammaEndcap / 2 - 1)) {
          left = 1;
          right = 0;
        }
        if (nGamma == nGammaEndcap / 2) {
          left = 0;
          right = -1;
        }
        if (nGamma > nGammaEndcap / 2) {
          left = -1;
          right = -1;
        }
      } else {
        if (nGamma < nGammaEndcap / 2) {
          left = 1;
          right = 1;
        }
        if (nGamma == nGammaEndcap / 2) {
          left = 1;
          right = -1;
        }
        if (nGamma > nGammaEndcap / 2) {
          left = -1;
          right = -1;
        }
      }
      if (nGamma == 0)
        left = 0;
      if (nGamma == nGammaEndcap - 1)
        right = 0;

      double x0y0l = x0y0 * tan_g_m_dg2;
      double x0y0r = x0y0 * tan_g_p_dg2;
      double x1y0l = x1y0 * tan_g_m_dg2;
      double x1y0r = x1y0 * tan_g_p_dg2;

      double x0y1l = x0y1 * tan_g_m_dg2 + left * r1_x_gamma_shift;
      double x0y1r = x0y1 * tan_g_p_dg2 + right * r1_x_gamma_shift;
      double x1y1l = x1y1 * tan_g_m_dg2 + left * r1_x_gamma_shift;
      double x1y1r = x1y1 * tan_g_p_dg2 + right * r1_x_gamma_shift;

      double x0y2l = x0y2 * tan_g_m_dg2 + left * r2_x_gamma_shift;
      double x0y2r = x0y2 * tan_g_p_dg2 + right * r2_x_gamma_shift;
      double x1y2l = x1y2 * tan_g_m_dg2 + left * r2_x_gamma_shift;
      double x1y2r = x1y2 * tan_g_p_dg2 + right * r2_x_gamma_shift;

      double verticesF[] = {x0y0r, y0e, x1y0r, -y0e, x1y0l, -y0e, x0y0l, y0e,
                            x0y1r, y1e, x1y1r, -y1e, x1y1l, -y1e, x0y1l, y1e};
      double verticesR[] = {x0y1r, y1e, x1y1r, -y1e, x1y1l, -y1e, x0y1l, y1e,
                            x0y2r, y2e, x1y2r, -y2e, x1y2l, -y2e, x0y2l, y2e};

      std::array<double, 16> vsub, vsub1;
      std::array<double, 2> center, center1;
      const int iTheta_neg = 2 * p.N_THETA_ENDCAP + p.N_THETA_BARREL - 1 - iTheta;

      for (int i = 0; i < p.XTAL_DIV_F; i++) {
        for (int j = 0; j < p.XTAL_DIV_F; j++) {
          int nEpsilon = i * p.XTAL_DIV_F + j;

          getSingleCrystalVertices(i, j, p.XTAL_DIV_F, verticesF, false, vsub);
          getSingleCrystalCenter(vsub, center);
          getSingleCrystalVertices(i, j, p.XTAL_DIV_F, verticesF, true, vsub1);
          getSingleCrystalCenter(vsub1, center1);

          double rGlobal = r0e + p.XTAL_LEN_F / 2;
          XYZVector dispGlobal(-center[1], center[0], rGlobal);
          XYZVector posGlobal = (rotYthGlobal * dispGlobal + p.DISP_PROJ_R);
          XYZVector dispGlobal_1(-center1[1], center1[0], -rGlobal);
          XYZVector posGlobal_1 = (rotYthGlobal_1 * dispGlobal_1 + p.DISP_PROJ_R);

          placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, PbWO4_to_ESR,
                       "EndcapCrystalF", p.XTAL_LEN_F / 2, vsub, p.crystalFMaterial, p.crystalFVis, p.trans_dispFsub,
                       endcapThetaAssemblyVolume, p.ENDCAP_SYSTEM_NO, iTheta, nGamma, nEpsilon, 0, posGlobal,
                       p.ENDCAP_PHI_START, p.ENDCAP_PHI_END, phi_endcap_rotations);
          numCrystals++;

          placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, PbWO4_to_ESR,
                       "EndcapCrystalF_1", p.XTAL_LEN_F / 2, vsub1, p.crystalFMaterial, p.crystalFVis,
                       p.trans_dispFsub_1, endcapThetaAssemblyVolume_1, p.ENDCAP_SYSTEM_NO, iTheta_neg, nGamma,
                       nEpsilon, 0, posGlobal_1, p.ENDCAP_PHI_START, p.ENDCAP_PHI_END, phi_endcap_rotations);
          numCrystals++;
        }
      }

      for (int i = 0; i < p.XTAL_DIV_R; i++) {
        for (int j = 0; j < p.XTAL_DIV_R; j++) {
          int nEpsilon = i * p.XTAL_DIV_R + j;

          getSingleCrystalVertices(i, j, p.XTAL_DIV_R, verticesR, false, vsub);
          getSingleCrystalCenter(vsub, center);
          getSingleCrystalVertices(i, j, p.XTAL_DIV_R, verticesR, true, vsub1);
          getSingleCrystalCenter(vsub1, center1);

          double rGlobal = r0e + p.XTAL_LEN_F + p.XTAL_LEN_R / 2;
          XYZVector dispGlobal(-center[1], center[0], rGlobal);
          XYZVector posGlobal = (rotYthGlobal * dispGlobal + p.DISP_PROJ_R);
          XYZVector dispGlobal_1(-center1[1], center1[0], -rGlobal);
          XYZVector posGlobal_1 = (rotYthGlobal_1 * dispGlobal_1 + p.DISP_PROJ_R);

          placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, PbWO4_to_ESR,
                       "EndcapCrystalR", p.XTAL_LEN_R / 2, vsub, p.crystalRMaterial, p.crystalRVis, p.trans_dispRsub,
                       endcapThetaAssemblyVolume, p.ENDCAP_SYSTEM_NO, iTheta, nGamma, nEpsilon, 1, posGlobal,
                       p.ENDCAP_PHI_START, p.ENDCAP_PHI_END, phi_endcap_rotations);
          numCrystals++;

          placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, PbWO4_to_ESR,
                       "EndcapCrystalR_1", p.XTAL_LEN_R / 2, vsub1, p.crystalRMaterial, p.crystalRVis,
                       p.trans_dispRsub_1, endcapThetaAssemblyVolume_1, p.ENDCAP_SYSTEM_NO, iTheta_neg, nGamma,
                       nEpsilon, 1, posGlobal_1, p.ENDCAP_PHI_START, p.ENDCAP_PHI_END, phi_endcap_rotations);
          numCrystals++;
        }
      }
    }
  }

  unsigned int nPhiSlice = 0;
  for (int iPhi = p.ENDCAP_PHI_START; iPhi < p.ENDCAP_PHI_END; iPhi++) {
    if (p.PHI_LOAD_START <= p.PHI_LOAD_END) {
      if ((iPhi >= p.PHI_LOAD_START) && (iPhi <= p.PHI_LOAD_END)) {
        auto pv =
            endcapGlobalAssemblyVol.placeVolume(endcapPhiAssemblyVolume, Transform3D(RotationZ(iPhi * p.D_PHI_GLOBAL)));
        auto pv1 = endcapGlobalAssemblyVol_1.placeVolume(endcapPhiAssemblyVolume_1,
                                                         Transform3D(RotationZ(iPhi * p.D_PHI_GLOBAL)));
        nPhiSlice++;
        pv.addPhysVolID("phi", iPhi);
        pv1.addPhysVolID("phi", iPhi);
      }
    } else {
      if ((iPhi >= p.PHI_LOAD_START) || (iPhi <= p.PHI_LOAD_END)) {
        auto pv =
            endcapGlobalAssemblyVol.placeVolume(endcapPhiAssemblyVolume, Transform3D(RotationZ(iPhi * p.D_PHI_GLOBAL)));
        auto pv1 = endcapGlobalAssemblyVol_1.placeVolume(endcapPhiAssemblyVolume_1,
                                                         Transform3D(RotationZ(iPhi * p.D_PHI_GLOBAL)));
        nPhiSlice++;
        pv.addPhysVolID("phi", iPhi);
        pv1.addPhysVolID("phi", iPhi);
      }
    }
  }
  return numCrystals * static_cast<int>(nPhiSlice);
}

// -----------------------------------------------------------------------
// Detector factory
// -----------------------------------------------------------------------

static dd4hep::Ref_t create_detector_SCEPCal_MainLayer(dd4hep::Detector& theDetector, xml_h xmlElement,
                                                       dd4hep::SensitiveDetector sens) {
  xml_det_t detectorXML = xmlElement;
  xml_comp_t dimXML = detectorXML.child(_Unicode(dim));
  xml_comp_t barrelXML = detectorXML.child(_Unicode(barrel));
  xml_comp_t endcapXML = detectorXML.child(_Unicode(endcap));
  xml_comp_t crystalFXML = detectorXML.child(_Unicode(crystalF));
  xml_comp_t crystalRXML = detectorXML.child(_Unicode(crystalR));

  xml_comp_t barrelAssemblyGlobalVisXML = detectorXML.child(_Unicode(barrelAssemblyGlobalVis));
  xml_comp_t barrelAssemblyPhiVisXML = detectorXML.child(_Unicode(barrelAssemblyPhiVis));
  xml_comp_t barrelAssemblyThetaVisXML = detectorXML.child(_Unicode(barrelAssemblyThetaVis));

  xml_comp_t endcapAssemblyGlobalVisXML = detectorXML.child(_Unicode(endcapAssemblyGlobalVis));
  xml_comp_t endcapAssemblyPhiVisXML = detectorXML.child(_Unicode(endcapAssemblyPhiVis));
  xml_comp_t endcapAssemblyThetaVisXML = detectorXML.child(_Unicode(endcapAssemblyThetaVis));

  std::string detName = detectorXML.nameStr();

  const double BARREL_HALF_Z = dimXML.attr<double>(_Unicode(barrelHalfZ));
  const double BARREL_INNER_R = dimXML.attr<double>(_Unicode(barrelInnerR));

  const double XTAL_TH_WIDTH = dimXML.attr<double>(_Unicode(crystalThetaWidth));
  const int XTAL_DIV_F = dimXML.attr<int>(_Unicode(crystalDivisionsF));
  const int XTAL_DIV_R = dimXML.attr<int>(_Unicode(crystalDivisionsR));
  const double XTAL_LEN_F = dimXML.attr<double>(_Unicode(crystalFlength));
  const double XTAL_LEN_R = dimXML.attr<double>(_Unicode(crystalRlength));

  const double BEAMPIPE_OPENING = dimXML.attr<double>(_Unicode(beampipe_opening));
  const double REAR_GAP = dimXML.attr<double>(_Unicode(rear_gap));

  const int PHI_SEGMENTS = dimXML.attr<int>(_Unicode(phiSegments));
  const double PROJ_OFFSET_R = dimXML.attr<double>(_Unicode(projectiveOffsetR));
  const double PROJ_OFFSET_X = dimXML.attr<double>(_Unicode(projectiveOffsetX));

  const bool USE_OPTICAL_SURFACES = dimXML.attr<bool>(_Unicode(useOpticalSurfaces));

  const int PHI_LOAD_START = dimXML.attr<int>(_Unicode(phi_load_start));
  const int PHI_LOAD_END = dimXML.attr<int>(_Unicode(phi_load_end));
  const int THETA_LOAD_START = dimXML.attr<int>(_Unicode(theta_load_start));
  const int THETA_LOAD_END = dimXML.attr<int>(_Unicode(theta_load_end));
  const int GAMMA_LOAD_START = dimXML.attr<int>(_Unicode(gamma_load_start));
  const int GAMMA_LOAD_END = dimXML.attr<int>(_Unicode(gamma_load_end));

  const int BARREL_SYSTEM_NO = barrelXML.attr<int>(_Unicode(system));
  const bool CONSTRUCT_BARREL = barrelXML.attr<bool>(_Unicode(construct));
  const int BARREL_PHI_START = 0;
  const int BARREL_PHI_END = PHI_SEGMENTS;

  const int ENDCAP_SYSTEM_NO = endcapXML.attr<int>(_Unicode(system));
  const bool CONSTRUCT_ENDCAP = endcapXML.attr<bool>(_Unicode(construct));
  const int ENDCAP_PHI_START = 0;
  const int ENDCAP_PHI_END = PHI_SEGMENTS;

  const double D_PHI_GLOBAL = 2 * M_PI / PHI_SEGMENTS;

  double THETA_SIZE_BARREL = atan(BARREL_HALF_Z / (BARREL_INNER_R + PROJ_OFFSET_R));
  double THETA_SIZE_ENDCAP = atan((BARREL_INNER_R + PROJ_OFFSET_R) / BARREL_HALF_Z);

  int N_THETA_BARREL = 2 * floor(BARREL_HALF_Z / XTAL_TH_WIDTH);
  int N_THETA_ENDCAP = floor(BARREL_INNER_R / XTAL_TH_WIDTH);
  double D_THETA_BARREL = 2 * THETA_SIZE_BARREL / N_THETA_BARREL;
  double D_THETA_ENDCAP = THETA_SIZE_ENDCAP / N_THETA_ENDCAP;

  int N_GAMMA_BARREL = std::max(int(2 * M_PI * BARREL_INNER_R / (PHI_SEGMENTS * XTAL_TH_WIDTH)), 1);
  double D_GAMMA_BARREL = D_PHI_GLOBAL / N_GAMMA_BARREL;

  XYZVector DISP_PROJ_R(-PROJ_OFFSET_R, 0, 0);

  int ENDCAP_THETA_START = 0;
  for (int iTheta = 0; iTheta < N_THETA_ENDCAP; iTheta++) {
    double thC = D_THETA_ENDCAP / 2 + iTheta * D_THETA_ENDCAP;
    double RinEndcap = BARREL_HALF_Z * tan(thC);
    if (RinEndcap < (BEAMPIPE_OPENING + PROJ_OFFSET_R))
      ENDCAP_THETA_START++;
  }

  // Barrel phi-slice envelope geometry
  double thC_br_end = THETA_SIZE_BARREL - D_THETA_BARREL / 2;
  double thC_br_beg = D_THETA_BARREL / 2;

  double r0_br_end = (BARREL_INNER_R + PROJ_OFFSET_R) / cos(thC_br_end);
  double r0_br_beg = (BARREL_INNER_R + PROJ_OFFSET_R) / cos(thC_br_beg);

  double r0_proj_arm_br_end = r0_br_end / cos(D_THETA_BARREL / 2);
  double r2_proj_arm_br_beg = (r0_br_beg + XTAL_LEN_F + XTAL_LEN_R) / cos(D_THETA_BARREL / 2);

  double br_phislice_8pa_z1 = r0_proj_arm_br_end * cos(thC_br_end + D_THETA_BARREL / 2) - PROJ_OFFSET_R;
  double br_phislice_8pa_z2 = r2_proj_arm_br_beg * cos(thC_br_beg - D_THETA_BARREL / 2) - PROJ_OFFSET_R + REAR_GAP;
  double br_phislice_8pa_y0 = (br_phislice_8pa_z1 + PROJ_OFFSET_R) * tan(thC_br_end + D_THETA_BARREL / 2);
  double br_phislice_8pa_y2 = (br_phislice_8pa_z2 + PROJ_OFFSET_R) * tan(thC_br_end + D_THETA_BARREL / 2);

  // Endcap phi-slice envelope geometry
  double thC_ec_end = THETA_SIZE_ENDCAP - D_THETA_ENDCAP / 2;
  double thC_ec_beg = D_THETA_ENDCAP / 2 + ENDCAP_THETA_START * D_THETA_ENDCAP;

  double r0_ec_end = BARREL_HALF_Z / cos(thC_ec_end);
  double r0_ec_beg = BARREL_HALF_Z / cos(thC_ec_beg);

  double r0_proj_arm_ec_end = r0_ec_end / cos(D_THETA_ENDCAP / 2);
  double r2_proj_arm_ec_beg = (r0_ec_beg + XTAL_LEN_F + XTAL_LEN_R) / cos(D_THETA_ENDCAP / 2);

  double ec_phislice_8pa_z1 = r0_proj_arm_ec_end * cos(thC_ec_end + D_THETA_ENDCAP / 2);
  double ec_phislice_8pa_z2 = r2_proj_arm_ec_beg * cos(thC_ec_beg - D_THETA_ENDCAP / 2) + REAR_GAP;

  double ec_end_phislice_8pa_y0 = ec_phislice_8pa_z1 * tan(thC_ec_end + D_THETA_ENDCAP / 2) - PROJ_OFFSET_R;
  double ec_end_phislice_8pa_y2 = ec_phislice_8pa_z2 * tan(thC_ec_end + D_THETA_ENDCAP / 2) - PROJ_OFFSET_R;
  double ec_beg_phislice_8pa_y0 = ec_phislice_8pa_z1 * tan(thC_ec_beg - D_THETA_ENDCAP / 2) - PROJ_OFFSET_R;
  double ec_beg_phislice_8pa_y2 = ec_phislice_8pa_z2 * tan(thC_ec_beg - D_THETA_ENDCAP / 2) - PROJ_OFFSET_R;

  Transform3D trans_dispFsub(Position(0, 0, -XTAL_LEN_R / 2));
  Transform3D trans_dispRsub(Position(0, 0, XTAL_LEN_F / 2));
  Transform3D trans_dispFsub_1(Position(0, 0, XTAL_LEN_R / 2));
  Transform3D trans_dispRsub_1(Position(0, 0, -XTAL_LEN_F / 2));

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "=SCEPCAL MAIN LAYER INPUTS=======================" << std::endl;
  std::cout << std::endl;
  std::cout << "BARREL_HALF_Z:        " << BARREL_HALF_Z << std::endl;
  std::cout << "BARREL_INNER_R:       " << BARREL_INNER_R << std::endl;
  std::cout << "PHI_SEGMENTS:         " << PHI_SEGMENTS << std::endl;
  std::cout << std::endl;
  std::cout << "XTAL_TH_WIDTH:        " << XTAL_TH_WIDTH << std::endl;
  std::cout << "XTAL_DIV_F:           " << XTAL_DIV_F << std::endl;
  std::cout << "XTAL_DIV_R:           " << XTAL_DIV_R << std::endl;
  std::cout << std::endl;
  std::cout << "XTAL_LEN_F:           " << XTAL_LEN_F << std::endl;
  std::cout << "XTAL_LEN_R:           " << XTAL_LEN_R << std::endl;
  std::cout << std::endl;
  std::cout << "PROJECTIVE_OFFSET_R:  " << PROJ_OFFSET_R << std::endl;
  std::cout << "PROJECTIVE_OFFSET_X:  " << PROJ_OFFSET_X << std::endl;
  std::cout << "BEAMPIPE_OPENING:     " << BEAMPIPE_OPENING << std::endl;
  std::cout << "REAR_GAP:             " << REAR_GAP << std::endl;
  std::cout << std::endl;
  std::cout << "=PARTIAL VISUALISATION PARAMETERS================" << std::endl;
  std::cout << "PHI_LOAD_START:       " << PHI_LOAD_START << std::endl;
  std::cout << "PHI_LOAD_END:         " << PHI_LOAD_END << std::endl;
  std::cout << "THETA_LOAD_START:     " << THETA_LOAD_START << std::endl;
  std::cout << "THETA_LOAD_END:       " << THETA_LOAD_END << std::endl;
  std::cout << "GAMMA_LOAD_START:     " << GAMMA_LOAD_START << std::endl;
  std::cout << "GAMMA_LOAD_END:       " << GAMMA_LOAD_END << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "=CONTROL=========================================" << std::endl;
  std::cout << std::endl;
  std::cout << "CONSTRUCT_BARREL:     " << CONSTRUCT_BARREL << std::endl;
  std::cout << "CONSTRUCT_ENDCAP:     " << CONSTRUCT_ENDCAP << std::endl;
  std::cout << "USE_OPTICAL_SURFACES: " << USE_OPTICAL_SURFACES << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "=CALCULATED PARAMETERS===========================" << std::endl;
  std::cout << std::endl;
  std::cout << "N_THETA_BARREL:       " << N_THETA_BARREL << std::endl;
  std::cout << "N_GAMMA_BARREL:       " << N_GAMMA_BARREL << std::endl;
  std::cout << std::endl;
  std::cout << "N_THETA_ENDCAP:       " << N_THETA_ENDCAP << std::endl;
  std::cout << "ENDCAP_THETA_START:   " << ENDCAP_THETA_START << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  // Initialization
  dd4hep::DetElement ScepcalDetElement(detName, detectorXML.id());
  dd4hep::Volume experimentalHall = theDetector.pickMotherVolume(ScepcalDetElement);
  dd4hep::xml::Dimension sdType = detectorXML.child(_Unicode(sensitive));
  sens.setType(sdType.typeStr());
  dd4hep::Readout readout = sens.readout();
  dd4hep::Segmentation geomseg = readout.segmentation();
  dd4hep::Segmentation* _geoSeg = &geomseg;
  auto segmentation = dynamic_cast<dd4hep::DDSegmentation::SCEPCal_MainSegmentation_k4geo*>(_geoSeg->segmentation());

  dd4hep::OpticalSurfaceManager surfMgr = theDetector.surfaceManager();
  dd4hep::OpticalSurface PbWO4_to_ESR = surfMgr.opticalSurface("/world/" + detName + "#PbWO4_to_ESR");

  segmentation->setDetIdBarrel(BARREL_SYSTEM_NO);
  segmentation->setDetIdEndcap(ENDCAP_SYSTEM_NO);
  segmentation->setIThetaBarrelStart(N_THETA_ENDCAP);
  segmentation->setIThetaBarrelEnd(N_THETA_BARREL + N_THETA_ENDCAP - 1);
  segmentation->setNPhi(PHI_SEGMENTS);
  segmentation->setNGamma(N_GAMMA_BARREL);

  std::map<int, int> nGammaPerTheta;
  for (int iTheta = ENDCAP_THETA_START; iTheta < N_THETA_ENDCAP; ++iTheta) {
    double thC = D_THETA_ENDCAP / 2 + iTheta * D_THETA_ENDCAP;
    double RinEndcap = BARREL_HALF_Z * std::tan(thC);
    int nGammaEndcap = std::max(int(2 * M_PI * RinEndcap / (PHI_SEGMENTS * XTAL_TH_WIDTH)), 1);
    nGammaPerTheta[iTheta] = nGammaEndcap;
    nGammaPerTheta[2 * N_THETA_ENDCAP + N_THETA_BARREL - 1 - iTheta] = nGammaEndcap;
  }
  segmentation->setNGammaPerTheta(nGammaPerTheta);

  dd4hep::xml::setDetectorTypeFlag(detectorXML, ScepcalDetElement);

  auto extensionData = new dd4hep::rec::LayeredCalorimeterData;
  const double XTAL_LEN = XTAL_LEN_F + XTAL_LEN_R;
  const double thEcInner = thC_ec_beg - D_THETA_ENDCAP / 2;

  extensionData->extent[0] = BARREL_INNER_R;
  extensionData->extent[1] = BARREL_INNER_R + XTAL_LEN;
  extensionData->extent[2] = BARREL_HALF_Z;
  extensionData->extent[3] = BARREL_HALF_Z + XTAL_LEN * std::cos(thEcInner);
  extensionData->extent[4] = BARREL_HALF_Z * std::tan(thEcInner) - PROJ_OFFSET_R;
  extensionData->extent[5] = BARREL_INNER_R + XTAL_LEN * std::sin(THETA_SIZE_ENDCAP);
  extensionData->layers.resize(2);
  extensionData->layers.at(0).distance = BARREL_INNER_R;
  extensionData->layers.at(1).distance = BARREL_INNER_R + XTAL_LEN_F;
  extensionData->layers.at(0).sensitive_thickness = XTAL_LEN_F;
  extensionData->layers.at(1).sensitive_thickness = XTAL_LEN_R;
  extensionData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;
  ScepcalDetElement.addExtension<dd4hep::rec::LayeredCalorimeterData>(extensionData);

  // Global assembly volumes
  std::vector<double> zBarrelPolyhedra = {-br_phislice_8pa_y2, -br_phislice_8pa_y0, br_phislice_8pa_y0,
                                          br_phislice_8pa_y2};
  std::vector<double> rminBarrelPolyhedra = {br_phislice_8pa_z2, br_phislice_8pa_z1, br_phislice_8pa_z1,
                                             br_phislice_8pa_z2};
  std::vector<double> rmaxBarrelPolyhedra = {br_phislice_8pa_z2, br_phislice_8pa_z2, br_phislice_8pa_z2,
                                             br_phislice_8pa_z2};

  std::vector<double> zEndcapPolyhedra = {ec_phislice_8pa_z1, ec_phislice_8pa_z2};
  std::vector<double> rminEndcapPolyhedra = {ec_beg_phislice_8pa_y0, ec_beg_phislice_8pa_y2};
  std::vector<double> rmaxEndcapPolyhedra = {ec_end_phislice_8pa_y0, ec_end_phislice_8pa_y2};

  std::vector<double> zEndcapPolyhedra_1 = {-ec_phislice_8pa_z2, -ec_phislice_8pa_z1};
  std::vector<double> rminEndcapPolyhedra_1 = {ec_beg_phislice_8pa_y2, ec_beg_phislice_8pa_y0};
  std::vector<double> rmaxEndcapPolyhedra_1 = {ec_end_phislice_8pa_y2, ec_end_phislice_8pa_y0};

  dd4hep::Polyhedra barrelGlobalAssemblyShape(PHI_SEGMENTS, D_PHI_GLOBAL / 2, 2 * M_PI, zBarrelPolyhedra,
                                              rminBarrelPolyhedra, rmaxBarrelPolyhedra);
  dd4hep::Volume barrelGlobalAssemblyVol("barrelGlobalAssemblyVol", barrelGlobalAssemblyShape,
                                         theDetector.material("Vacuum"));
  barrelGlobalAssemblyVol.setVisAttributes(theDetector, barrelAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume barrelAssemblyPlacedVol = experimentalHall.placeVolume(barrelGlobalAssemblyVol);
  barrelAssemblyPlacedVol.addPhysVolID("system", BARREL_SYSTEM_NO);

  dd4hep::Polyhedra endcapGlobalAssemblyShape(PHI_SEGMENTS, D_PHI_GLOBAL / 2, 2 * M_PI, zEndcapPolyhedra,
                                              rminEndcapPolyhedra, rmaxEndcapPolyhedra);
  dd4hep::Volume endcapGlobalAssemblyVol("endcapGlobalAssemblyVol", endcapGlobalAssemblyShape,
                                         theDetector.material("Vacuum"));
  endcapGlobalAssemblyVol.setVisAttributes(theDetector, endcapAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume endcapAssemblyPlacedVol = experimentalHall.placeVolume(endcapGlobalAssemblyVol);
  endcapAssemblyPlacedVol.addPhysVolID("system", ENDCAP_SYSTEM_NO);

  dd4hep::Polyhedra endcapGlobalAssemblyShape_1(PHI_SEGMENTS, D_PHI_GLOBAL / 2, 2 * M_PI, zEndcapPolyhedra_1,
                                                rminEndcapPolyhedra_1, rmaxEndcapPolyhedra_1);
  dd4hep::Volume endcapGlobalAssemblyVol_1("endcapGlobalAssemblyVol_1", endcapGlobalAssemblyShape_1,
                                           theDetector.material("Vacuum"));
  endcapGlobalAssemblyVol_1.setVisAttributes(theDetector, endcapAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume endcapAssemblyPlacedVol_1 = experimentalHall.placeVolume(endcapGlobalAssemblyVol_1);
  endcapAssemblyPlacedVol_1.addPhysVolID("system", ENDCAP_SYSTEM_NO);
  endcapAssemblyPlacedVol_1.addPhysVolID("theta", 1);

  ScepcalDetElement.setPlacement(barrelAssemblyPlacedVol);

  // Pre-build phi rotation maps
  std::unordered_map<int, RotationZ> phi_barrel_rotations, phi_endcap_rotations;
  for (int iPhi = BARREL_PHI_START; iPhi < BARREL_PHI_END; iPhi++)
    phi_barrel_rotations[iPhi] = RotationZ(iPhi * D_PHI_GLOBAL);
  for (int iPhi = ENDCAP_PHI_START; iPhi < ENDCAP_PHI_END; iPhi++)
    phi_endcap_rotations[iPhi] = RotationZ(iPhi * D_PHI_GLOBAL);

  // Pack shared parameters
  SCEPCalParams p;
  p.BARREL_HALF_Z = BARREL_HALF_Z;
  p.BARREL_INNER_R = BARREL_INNER_R;
  p.XTAL_LEN_F = XTAL_LEN_F;
  p.XTAL_LEN_R = XTAL_LEN_R;
  p.XTAL_DIV_F = XTAL_DIV_F;
  p.XTAL_DIV_R = XTAL_DIV_R;
  p.XTAL_TH_WIDTH = XTAL_TH_WIDTH;
  p.PROJ_OFFSET_R = PROJ_OFFSET_R;
  p.PROJ_OFFSET_X = PROJ_OFFSET_X;
  p.PHI_SEGMENTS = PHI_SEGMENTS;
  p.D_PHI_GLOBAL = D_PHI_GLOBAL;

  p.N_THETA_BARREL = N_THETA_BARREL;
  p.N_GAMMA_BARREL = N_GAMMA_BARREL;
  p.D_THETA_BARREL = D_THETA_BARREL;
  p.D_GAMMA_BARREL = D_GAMMA_BARREL;
  p.THETA_SIZE_BARREL = THETA_SIZE_BARREL;
  p.THETA_SIZE_ENDCAP = THETA_SIZE_ENDCAP;
  p.BARREL_SYSTEM_NO = BARREL_SYSTEM_NO;
  p.BARREL_PHI_START = BARREL_PHI_START;
  p.BARREL_PHI_END = BARREL_PHI_END;
  p.CONSTRUCT_BARREL = CONSTRUCT_BARREL;

  p.N_THETA_ENDCAP = N_THETA_ENDCAP;
  p.ENDCAP_THETA_START = ENDCAP_THETA_START;
  p.D_THETA_ENDCAP = D_THETA_ENDCAP;
  p.ENDCAP_SYSTEM_NO = ENDCAP_SYSTEM_NO;
  p.ENDCAP_PHI_START = ENDCAP_PHI_START;
  p.ENDCAP_PHI_END = ENDCAP_PHI_END;
  p.CONSTRUCT_ENDCAP = CONSTRUCT_ENDCAP;

  p.PHI_LOAD_START = PHI_LOAD_START;
  p.PHI_LOAD_END = PHI_LOAD_END;
  p.THETA_LOAD_START = THETA_LOAD_START;
  p.THETA_LOAD_END = THETA_LOAD_END;
  p.GAMMA_LOAD_START = GAMMA_LOAD_START;
  p.GAMMA_LOAD_END = GAMMA_LOAD_END;

  p.crystalFMaterial = crystalFXML.materialStr();
  p.crystalFVis = crystalFXML.visStr();
  p.crystalRMaterial = crystalRXML.materialStr();
  p.crystalRVis = crystalRXML.visStr();
  p.barrelAssemblyPhiVis = barrelAssemblyPhiVisXML.visStr();
  p.barrelAssemblyThetaVis = barrelAssemblyThetaVisXML.visStr();
  p.endcapAssemblyPhiVis = endcapAssemblyPhiVisXML.visStr();
  p.endcapAssemblyThetaVis = endcapAssemblyThetaVisXML.visStr();

  p.trans_dispFsub = trans_dispFsub;
  p.trans_dispRsub = trans_dispRsub;
  p.trans_dispFsub_1 = trans_dispFsub_1;
  p.trans_dispRsub_1 = trans_dispRsub_1;
  p.DISP_PROJ_R = DISP_PROJ_R;

  p.zBarrelPolyhedra = zBarrelPolyhedra;
  p.rminBarrelPolyhedra = rminBarrelPolyhedra;
  p.rmaxBarrelPolyhedra = rmaxBarrelPolyhedra;

  p.zEndcapPolyhedra = zEndcapPolyhedra;
  p.rminEndcapPolyhedra = rminEndcapPolyhedra;
  p.rmaxEndcapPolyhedra = rmaxEndcapPolyhedra;
  p.zEndcapPolyhedra_1 = zEndcapPolyhedra_1;
  p.rminEndcapPolyhedra_1 = rminEndcapPolyhedra_1;
  p.rmaxEndcapPolyhedra_1 = rmaxEndcapPolyhedra_1;

  int numCrystalsBarrel = 0;
  if (CONSTRUCT_BARREL)
    numCrystalsBarrel = buildBarrel(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES,
                                    PbWO4_to_ESR, p, phi_barrel_rotations, barrelGlobalAssemblyVol);

  int numCrystalsEndcap = 0;
  if (CONSTRUCT_ENDCAP)
    numCrystalsEndcap = buildEndcap(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES,
                                    PbWO4_to_ESR, p, phi_endcap_rotations, endcapGlobalAssemblyVol,
                                    endcapGlobalAssemblyVol_1);

  std::cout << std::endl;
  std::cout << "NUM_CRYSTALS_BARREL:  " << numCrystalsBarrel << std::endl;
  std::cout << "NUM_CRYSTALS_ENDCAP:  " << numCrystalsEndcap << std::endl;
  std::cout << std::endl;

  return ScepcalDetElement;
}

DECLARE_DETELEMENT(SCEPCal_MainLayer, create_detector_SCEPCal_MainLayer)

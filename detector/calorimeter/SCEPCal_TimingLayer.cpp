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
#include "detectorSegmentations/SCEPCal_TimingSegmentationHandle_k4geo.h"
#include <bitset>
#include <unordered_map>

using dd4hep::Position;
using dd4hep::RotationZYX;
using dd4hep::Transform3D;
using ROOT::Math::RotationY;
using ROOT::Math::RotationZ;
using ROOT::Math::XYZVector;
using namespace dd4hep;

// -----------------------------------------------------------------------
// Crystal placement helper
// -----------------------------------------------------------------------

static void placeCrystal(dd4hep::Detector& theDetector, dd4hep::SensitiveDetector& sens,
                         dd4hep::DDSegmentation::SCEPCal_TimingSegmentation_k4geo* segmentation,
                         dd4hep::DetElement& ScepcalDetElement, bool USE_OPTICAL_SURFACES,
                         dd4hep::OpticalSurface& LYSO_to_ESR, const std::string& volName, double dz,
                         const double* vertices, const std::string& materialStr, const std::string& visStr,
                         const dd4hep::Transform3D& transform, const dd4hep::Volume& assemblyVol, int nSystem,
                         int nTheta, int nGamma, const XYZVector& posGlobal, int phi_start, int phi_end,
                         const std::vector<RotationZ>& phi_rotations, int phi0) {
  dd4hep::EightPointSolid theShape(dz, vertices);
  dd4hep::Volume theVolume(volName, theShape, theDetector.material(materialStr));
  theVolume.setVisAttributes(theDetector, visStr);

  if (USE_OPTICAL_SURFACES) {
    dd4hep::SkinSurface(theDetector, ScepcalDetElement,
                        volName + "Surface_" + std::to_string(nTheta) + "_" + std::to_string(nGamma), LYSO_to_ESR,
                        theVolume);
  }

  theVolume.setSensitiveDetector(sens);

  auto volID_0 = segmentation->setVolumeID(nSystem, 0, nTheta, nGamma);
  int volID_0_32 = segmentation->getFirst32bits(volID_0);
  dd4hep::PlacedVolume thePlacedVol = assemblyVol.placeVolume(theVolume, volID_0_32, transform);

  thePlacedVol.addPhysVolID("system", nSystem);
  thePlacedVol.addPhysVolID("theta", nTheta);
  thePlacedVol.addPhysVolID("gamma", nGamma);

  for (int iPhi = phi_start; iPhi < phi_end; iPhi++) {
    XYZVector posGlobalPhi = phi_rotations[iPhi + phi0] * posGlobal;
    auto volID = segmentation->setVolumeID(nSystem, iPhi, nTheta, nGamma);
    int volID_32 = segmentation->getFirst32bits(volID);
    segmentation->savePosition(volID_32, posGlobalPhi);
  }
}

// -----------------------------------------------------------------------
// Shared geometry parameters
// -----------------------------------------------------------------------

struct SCEPCalTimingParams {
  // Basic dimensions
  double BARREL_HALF_Z, BARREL_INNER_R;
  double XTAL_DEPTH_T, XTAL_LENGTH_T, XTAL_TH_WIDTH;
  int PHI_SEGMENTS;
  double D_PHI_GLOBAL;

  // Timing barrel
  double rT; // inner radius of timing barrel crystals
  int N_THETA_TLBARREL, N_GAMMA_TLBARREL;
  double D_THETA_TLBARREL, D_GAMMA_TLBARREL;
  double THETA_SIZE_TLBARREL, THETA_SIZE_TLENDCAP;
  int TLBARREL_SYSTEM_NO, TIMING_PHI_START, TIMING_PHI_END;
  bool CONSTRUCT_TLBARREL;

  // Timing endcap
  double zT; // z-face of timing endcap crystals
  int N_THETA_TLENDCAP, TLENDCAP_THETA_START;
  double D_THETA_TLENDCAP;
  int TLENDCAP_SYSTEM_NO, TLENDCAP_PHI_START, TLENDCAP_PHI_END;
  bool CONSTRUCT_TLENDCAP;

  // Crystal material/vis strings
  std::string crystalTMaterial, crystalTVis;

  // Assembly vis strings
  std::string tlbarrelAssemblyPhiVis, tlbarrelAssemblyThetaVis;
  std::string tlendcapAssemblyPhiVis, tlendcapAssemblyThetaVis;

  // Timing barrel polyhedra
  std::vector<double> zTimingPolyhedra, rminTimingPolyhedra, rmaxTimingPolyhedra;

  // Timing endcap (+z) polyhedra
  std::vector<double> zTlEndcapPolyhedra, rminTlEndcapPolyhedra, rmaxTlEndcapPolyhedra;
  // Timing endcap (-z) polyhedra
  std::vector<double> zTlEndcapPolyhedra_1, rminTlEndcapPolyhedra_1, rmaxTlEndcapPolyhedra_1;
};

// -----------------------------------------------------------------------
// Timing barrel builder
// -----------------------------------------------------------------------

static int buildTimingBarrel(dd4hep::Detector& theDetector, dd4hep::SensitiveDetector& sens,
                             dd4hep::DDSegmentation::SCEPCal_TimingSegmentation_k4geo* segmentation,
                             dd4hep::DetElement& ScepcalDetElement, bool USE_OPTICAL_SURFACES,
                             dd4hep::OpticalSurface& LYSO_to_ESR, const SCEPCalTimingParams& p,
                             const std::vector<RotationZ>& phi_barrel_rotations,
                             dd4hep::Volume& timingGlobalAssemblyVol) {
  int numCrystals = 0;

  const double tan_d_phi_global2 = tan(p.D_PHI_GLOBAL / 2);
  const double tan_m_d_phi_global2 = -tan_d_phi_global2;
  const Transform3D trans_dispT(Position(0, 0, 0));

  dd4hep::Polyhedra timingPhiAssemblyShape(1, -p.D_PHI_GLOBAL / 2, p.D_PHI_GLOBAL, p.zTimingPolyhedra,
                                           p.rminTimingPolyhedra, p.rmaxTimingPolyhedra);
  dd4hep::Volume timingPhiAssemblyVolume("timingPhiAssembly", timingPhiAssemblyShape, theDetector.material("Vacuum"));
  timingPhiAssemblyVolume.setVisAttributes(theDetector, p.tlbarrelAssemblyPhiVis);

  for (int iTheta = 0; iTheta < p.N_THETA_TLBARREL; iTheta++) {
    double thC = p.THETA_SIZE_TLENDCAP + p.D_THETA_TLBARREL / 2 + (iTheta * p.D_THETA_TLBARREL);
    double sin_thC = sin(thC);
    double cos_thC = cos(thC);

    double r0e = p.rT / sin_thC;
    double r2e = r0e + p.XTAL_DEPTH_T;
    double y0e = r0e * tan(p.D_THETA_TLBARREL / 2.);
    double y2e = r2e * tan(p.D_THETA_TLBARREL / 2.);

    double x0y0 = r0e * sin_thC - y0e * cos_thC;
    double x1y0 = r0e * sin_thC + y0e * cos_thC;
    double x0y2 = r2e * sin_thC - y2e * cos_thC;
    double x1y2 = r2e * sin_thC + y2e * cos_thC;

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

    double rE = r0e + p.XTAL_DEPTH_T / 2.;
    RotationZYX rotE(M_PI / 2, thC, 0);
    Position dispE(rE * sin_thC, 0, rE * cos_thC);

    dd4hep::EightPointSolid timingThetaAssemblyShape(p.XTAL_DEPTH_T / 2, verticesE);
    dd4hep::Volume timingThetaAssemblyVolume("timingThetaAssembly", timingThetaAssemblyShape,
                                             theDetector.material("Vacuum"));
    timingThetaAssemblyVolume.setVisAttributes(theDetector, p.tlbarrelAssemblyThetaVis);
    timingPhiAssemblyVolume.placeVolume(timingThetaAssemblyVolume, Transform3D(rotE, dispE));

    for (int nGamma = 0; nGamma < p.N_GAMMA_TLBARREL; nGamma++) {
      double gamma = -p.D_PHI_GLOBAL / 2 + p.D_GAMMA_TLBARREL / 2 + p.D_GAMMA_TLBARREL * nGamma;
      double tan_g_m_dg = tan(gamma - p.D_GAMMA_TLBARREL / 2);
      double tan_g_p_dg = tan(gamma + p.D_GAMMA_TLBARREL / 2);

      double x0y0l = x0y0 * tan_g_m_dg;
      double x0y0r = x0y0 * tan_g_p_dg;
      double x1y0l = x1y0 * tan_g_m_dg;
      double x1y0r = x1y0 * tan_g_p_dg;
      double x0y2l = x0y2 * tan_g_m_dg;
      double x0y2r = x0y2 * tan_g_p_dg;
      double x1y2l = x1y2 * tan_g_m_dg;
      double x1y2r = x1y2 * tan_g_p_dg;

      double verticesT[] = {x0y0r, y0e, x1y0r, -y0e, x1y0l, -y0e, x0y0l, y0e,
                            x0y2r, y2e, x1y2r, -y2e, x1y2l, -y2e, x0y2l, y2e};

      XYZVector posGlobal(rE * sin_thC, rE * sin_thC * tan(gamma), rE * cos_thC);

      placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, LYSO_to_ESR,
                   "TimingCrystal", p.XTAL_DEPTH_T / 2, verticesT, p.crystalTMaterial, p.crystalTVis, trans_dispT,
                   timingThetaAssemblyVolume, p.TLBARREL_SYSTEM_NO, p.N_THETA_TLENDCAP + iTheta, nGamma, posGlobal,
                   p.TIMING_PHI_START, p.TIMING_PHI_END, phi_barrel_rotations, p.TIMING_PHI_START);
      numCrystals++;
    }
  }

  unsigned int nPhiPlaced = 0;
  for (int iPhi = p.TIMING_PHI_START; iPhi < p.TIMING_PHI_END; iPhi++) {
    double phiGlobal = iPhi * p.D_PHI_GLOBAL;
    RotationZ rotZphiGlobal(phiGlobal);
    auto pv = timingGlobalAssemblyVol.placeVolume(timingPhiAssemblyVolume, iPhi, Transform3D(rotZphiGlobal));
    pv.addPhysVolID("phi", iPhi);
    nPhiPlaced++;
  }
  return numCrystals * static_cast<int>(nPhiPlaced);
}

// -----------------------------------------------------------------------
// Timing endcap builder
// -----------------------------------------------------------------------

static int buildTimingEndcap(dd4hep::Detector& theDetector, dd4hep::SensitiveDetector& sens,
                             dd4hep::DDSegmentation::SCEPCal_TimingSegmentation_k4geo* segmentation,
                             dd4hep::DetElement& ScepcalDetElement, bool USE_OPTICAL_SURFACES,
                             dd4hep::OpticalSurface& LYSO_to_ESR, const SCEPCalTimingParams& p,
                             const std::vector<RotationZ>& phi_endcap_rotations,
                             dd4hep::Volume& tlendcapGlobalAssemblyVol, dd4hep::Volume& tlendcapGlobalAssemblyVol_1) {
  int numCrystals = 0;

  const double tan_d_theta_tl2 = tan(p.D_THETA_TLENDCAP / 2.);
  const double tan_d_phig2 = tan(p.D_PHI_GLOBAL / 2);
  const double tan_m_d_phig2 = -tan_d_phig2;

  dd4hep::Polyhedra tlendcapPhiAssemblyShape(1, -p.D_PHI_GLOBAL / 2, p.D_PHI_GLOBAL, p.zTlEndcapPolyhedra,
                                             p.rminTlEndcapPolyhedra, p.rmaxTlEndcapPolyhedra);
  dd4hep::Polyhedra tlendcapPhiAssemblyShape_1(1, -p.D_PHI_GLOBAL / 2, p.D_PHI_GLOBAL, p.zTlEndcapPolyhedra_1,
                                               p.rminTlEndcapPolyhedra_1, p.rmaxTlEndcapPolyhedra_1);

  dd4hep::Volume tlendcapPhiAssemblyVolume("tlendcapPhiVol", tlendcapPhiAssemblyShape, theDetector.material("Vacuum"));
  tlendcapPhiAssemblyVolume.setVisAttributes(theDetector, p.tlendcapAssemblyPhiVis);

  dd4hep::Volume tlendcapPhiAssemblyVolume_1("tlendcapPhiVol_1", tlendcapPhiAssemblyShape_1,
                                             theDetector.material("Vacuum"));
  tlendcapPhiAssemblyVolume_1.setVisAttributes(theDetector, p.tlendcapAssemblyPhiVis);

  for (int iTheta = p.TLENDCAP_THETA_START; iTheta < p.N_THETA_TLENDCAP; iTheta++) {
    double thC = p.D_THETA_TLENDCAP / 2 + iTheta * p.D_THETA_TLENDCAP;
    double sin_thC = sin(thC);
    double cos_thC = cos(thC);

    double RinEndcap = p.zT * tan(thC);
    int nGammaEndcap = std::max(int(2 * M_PI * RinEndcap / (p.PHI_SEGMENTS * p.XTAL_LENGTH_T)), 1);
    double dGammaEndcap = p.D_PHI_GLOBAL / nGammaEndcap;

    double r0e = RinEndcap / sin_thC;
    double r2e = r0e + p.XTAL_DEPTH_T;
    double y0e = r0e * tan_d_theta_tl2;
    double y2e = r2e * tan_d_theta_tl2;

    double x0y0 = r0e * sin_thC - y0e * cos_thC;
    double x1y0 = r0e * sin_thC + y0e * cos_thC;
    double x0y2 = r2e * sin_thC - y2e * cos_thC;
    double x1y2 = r2e * sin_thC + y2e * cos_thC;

    double x0y0l_E = x0y0 * tan_m_d_phig2;
    double x0y0r_E = x0y0 * tan_d_phig2;
    double x1y0l_E = x1y0 * tan_m_d_phig2;
    double x1y0r_E = x1y0 * tan_d_phig2;
    double x0y2l_E = x0y2 * tan_m_d_phig2;
    double x0y2r_E = x0y2 * tan_d_phig2;
    double x1y2l_E = x1y2 * tan_m_d_phig2;
    double x1y2r_E = x1y2 * tan_d_phig2;

    double verticesE[] = {x0y0r_E, y0e, x1y0r_E, -y0e, x1y0l_E, -y0e, x0y0l_E, y0e,
                          x0y2r_E, y2e, x1y2r_E, -y2e, x1y2l_E, -y2e, x0y2l_E, y2e};
    double verticesE_1[] = {x0y2r_E, y2e, x1y2r_E, -y2e, x1y2l_E, -y2e, x0y2l_E, y2e,
                            x0y0r_E, y0e, x1y0r_E, -y0e, x1y0l_E, -y0e, x0y0l_E, y0e};

    double rE = r0e + p.XTAL_DEPTH_T / 2.;
    RotationZYX rotE(M_PI / 2, thC, 0);
    RotationZYX rotE_1(M_PI / 2, -thC, 0);
    Position dispE(rE * sin_thC, 0, rE * cos_thC);
    Position dispE_1(rE * sin_thC, 0, -rE * cos_thC);

    dd4hep::EightPointSolid tlendcapThetaAssemblyShape(p.XTAL_DEPTH_T / 2, verticesE);
    dd4hep::EightPointSolid tlendcapThetaAssemblyShape_1(p.XTAL_DEPTH_T / 2, verticesE_1);

    dd4hep::Volume tlendcapThetaAssemblyVolume("tlendcapThetaAssembly", tlendcapThetaAssemblyShape,
                                               theDetector.material("Vacuum"));
    tlendcapThetaAssemblyVolume.setVisAttributes(theDetector, p.tlendcapAssemblyThetaVis);
    tlendcapPhiAssemblyVolume.placeVolume(tlendcapThetaAssemblyVolume, Transform3D(rotE, dispE));

    dd4hep::Volume tlendcapThetaAssemblyVolume_1("tlendcapThetaAssembly_1", tlendcapThetaAssemblyShape_1,
                                                 theDetector.material("Vacuum"));
    tlendcapThetaAssemblyVolume_1.setVisAttributes(theDetector, p.tlendcapAssemblyThetaVis);
    tlendcapPhiAssemblyVolume_1.placeVolume(tlendcapThetaAssemblyVolume_1, Transform3D(rotE_1, dispE_1));

    const int iTheta_neg = 2 * p.N_THETA_TLENDCAP + p.N_THETA_TLBARREL - 1 - iTheta;

    for (int nGamma = 0; nGamma < nGammaEndcap; nGamma++) {
      double gamma = -p.D_PHI_GLOBAL / 2 + dGammaEndcap / 2 + dGammaEndcap * nGamma;
      double tan_g_m_dg2 = tan(gamma - dGammaEndcap / 2);
      double tan_g_p_dg2 = tan(gamma + dGammaEndcap / 2);

      double x0y0l = x0y0 * tan_g_m_dg2;
      double x0y0r = x0y0 * tan_g_p_dg2;
      double x1y0l = x1y0 * tan_g_m_dg2;
      double x1y0r = x1y0 * tan_g_p_dg2;
      double x0y2l = x0y2 * tan_g_m_dg2;
      double x0y2r = x0y2 * tan_g_p_dg2;
      double x1y2l = x1y2 * tan_g_m_dg2;
      double x1y2r = x1y2 * tan_g_p_dg2;

      double verticesT[] = {x0y0r, y0e, x1y0r, -y0e, x1y0l, -y0e, x0y0l, y0e,
                            x0y2r, y2e, x1y2r, -y2e, x1y2l, -y2e, x0y2l, y2e};
      double verticesT_1[] = {x0y2r, y2e, x1y2r, -y2e, x1y2l, -y2e, x0y2l, y2e,
                              x0y0r, y0e, x1y0r, -y0e, x1y0l, -y0e, x0y0l, y0e};

      Position dispT(0, 0, 0);

      XYZVector posGlobal(rE * sin_thC, rE * sin_thC * tan(gamma), rE * cos_thC);
      XYZVector posGlobal_1(rE * sin_thC, rE * sin_thC * tan(gamma), -rE * cos_thC);

      placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, LYSO_to_ESR,
                   "TlEndcapCrystal", p.XTAL_DEPTH_T / 2, verticesT, p.crystalTMaterial, p.crystalTVis,
                   Transform3D(dispT), tlendcapThetaAssemblyVolume, p.TLENDCAP_SYSTEM_NO, iTheta, nGamma, posGlobal,
                   p.TLENDCAP_PHI_START, p.TLENDCAP_PHI_END, phi_endcap_rotations, p.TLENDCAP_PHI_START);
      numCrystals++;

      placeCrystal(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, LYSO_to_ESR,
                   "TlEndcapCrystal_1", p.XTAL_DEPTH_T / 2, verticesT_1, p.crystalTMaterial, p.crystalTVis,
                   Transform3D(dispT), tlendcapThetaAssemblyVolume_1, p.TLENDCAP_SYSTEM_NO, iTheta_neg, nGamma,
                   posGlobal_1, p.TLENDCAP_PHI_START, p.TLENDCAP_PHI_END, phi_endcap_rotations, p.TLENDCAP_PHI_START);
      numCrystals++;
    }
  }

  unsigned int nPhiPlaced = 0;
  for (int iPhi = p.TLENDCAP_PHI_START; iPhi < p.TLENDCAP_PHI_END; iPhi++) {
    double phiGlobal = iPhi * p.D_PHI_GLOBAL;
    RotationZ rotZphiGlobal(phiGlobal);
    auto pv = tlendcapGlobalAssemblyVol.placeVolume(tlendcapPhiAssemblyVolume, iPhi, Transform3D(rotZphiGlobal));
    pv.addPhysVolID("phi", iPhi);
    auto pv_1 = tlendcapGlobalAssemblyVol_1.placeVolume(tlendcapPhiAssemblyVolume_1, iPhi, Transform3D(rotZphiGlobal));
    pv_1.addPhysVolID("phi", iPhi);

    nPhiPlaced++;
  }
  return numCrystals * static_cast<int>(nPhiPlaced);
}

// -----------------------------------------------------------------------
// Detector factory
// -----------------------------------------------------------------------

static dd4hep::Ref_t create_detector_SCEPCal_TimingLayer(dd4hep::Detector& theDetector, xml_h xmlElement,
                                                         dd4hep::SensitiveDetector sens) {
  xml_det_t detectorXML = xmlElement;
  xml_comp_t dimXML = detectorXML.child(_Unicode(dim));
  xml_comp_t tlbarrelXML = detectorXML.child(_Unicode(tlbarrel));
  xml_comp_t tlendcapXML = detectorXML.child(_Unicode(tlendcap));
  xml_comp_t crystalTXML = detectorXML.child(_Unicode(crystalT));

  xml_comp_t tlbarrelAssemblyGlobalVisXML = detectorXML.child(_Unicode(tlbarrelAssemblyGlobalVis));
  xml_comp_t tlbarrelAssemblyThetaVisXML = detectorXML.child(_Unicode(tlbarrelAssemblyThetaVis));
  xml_comp_t tlbarrelAssemblyPhiVisXML = detectorXML.child(_Unicode(tlbarrelAssemblyPhiVis));

  xml_comp_t tlendcapAssemblyGlobalVisXML = detectorXML.child(_Unicode(tlendcapAssemblyGlobalVis));
  xml_comp_t tlendcapAssemblyPhiVisXML = detectorXML.child(_Unicode(tlendcapAssemblyPhiVis));
  xml_comp_t tlendcapAssemblyThetaVisXML = detectorXML.child(_Unicode(tlendcapAssemblyThetaVis));

  std::string detName = detectorXML.nameStr();

  const double BARREL_HALF_Z = dimXML.attr<double>(_Unicode(barrelHalfZ));
  const double BARREL_INNER_R = dimXML.attr<double>(_Unicode(barrelInnerR));

  const double XTAL_TH_WIDTH = dimXML.attr<double>(_Unicode(crystalThetaWidth));
  const double XTAL_DEPTH_T = dimXML.attr<double>(_Unicode(crystalTdepth));
  const double XTAL_LENGTH_T = dimXML.attr<double>(_Unicode(crystalTlength));

  const double BEAMPIPE_OPENING = dimXML.attr<double>(_Unicode(beampipe_opening));
  const double TIMING_GAP = dimXML.attr<double>(_Unicode(timing_gap));

  const int PHI_SEGMENTS = dimXML.attr<int>(_Unicode(phiSegments));
  const double PROJ_OFFSET_R = dimXML.attr<double>(_Unicode(projectiveOffsetR));

  const bool USE_OPTICAL_SURFACES = dimXML.attr<bool>(_Unicode(useOpticalSurfaces));

  const int TLBARREL_SYSTEM_NO = tlbarrelXML.attr<int>(_Unicode(system));
  const bool CONSTRUCT_TLBARREL = tlbarrelXML.attr<bool>(_Unicode(construct));
  const int TIMING_PHI_START = 0;
  const int TIMING_PHI_END = PHI_SEGMENTS;

  const int TLENDCAP_SYSTEM_NO = tlendcapXML.attr<int>(_Unicode(system));
  const bool CONSTRUCT_TLENDCAP = tlendcapXML.attr<bool>(_Unicode(construct));
  const int TLENDCAP_PHI_START = 0;
  const int TLENDCAP_PHI_END = PHI_SEGMENTS;

  const double D_PHI_GLOBAL = 2 * M_PI / PHI_SEGMENTS;

  double THETA_SIZE_BARREL = atan(BARREL_HALF_Z / (BARREL_INNER_R + PROJ_OFFSET_R));
  double THETA_SIZE_ENDCAP = atan((BARREL_INNER_R + PROJ_OFFSET_R) / BARREL_HALF_Z);

  double THETA_SIZE_TLBARREL = atan(BARREL_HALF_Z / BARREL_INNER_R);
  double THETA_SIZE_TLENDCAP = atan(BARREL_INNER_R / BARREL_HALF_Z);

  int N_THETA_BARREL = 2 * floor(BARREL_HALF_Z / XTAL_TH_WIDTH);
  int N_THETA_ENDCAP = floor(BARREL_INNER_R / XTAL_TH_WIDTH);
  int N_THETA_TLBARREL = 2 * floor(BARREL_HALF_Z / XTAL_DEPTH_T);
  int N_THETA_TLENDCAP = floor(BARREL_INNER_R / XTAL_DEPTH_T);

  double D_THETA_BARREL = 2 * THETA_SIZE_BARREL / N_THETA_BARREL;
  double D_THETA_ENDCAP = THETA_SIZE_ENDCAP / N_THETA_ENDCAP;
  double D_THETA_TLBARREL = 2 * THETA_SIZE_TLBARREL / N_THETA_TLBARREL;
  double D_THETA_TLENDCAP = THETA_SIZE_TLENDCAP / N_THETA_TLENDCAP;

  int N_GAMMA_TLBARREL = std::max(int(2 * M_PI * BARREL_INNER_R / (PHI_SEGMENTS * XTAL_LENGTH_T)), 1);
  double D_GAMMA_TLBARREL = D_PHI_GLOBAL / N_GAMMA_TLBARREL;

  // Inner radius/z of timing crystals (just inside the main crystal layer)
  double thC_br_end = THETA_SIZE_BARREL - D_THETA_BARREL / 2;
  double r0_br_end = (BARREL_INNER_R + PROJ_OFFSET_R) / cos(thC_br_end);
  double r0_proj_arm_br_end = r0_br_end / cos(D_THETA_BARREL / 2);
  double br_phislice_8pa_z1 = r0_proj_arm_br_end * cos(thC_br_end + D_THETA_BARREL / 2) - PROJ_OFFSET_R;
  double rT = br_phislice_8pa_z1 - TIMING_GAP - XTAL_DEPTH_T;

  double thC_ec_end = THETA_SIZE_ENDCAP - D_THETA_ENDCAP / 2;
  double r0_ec_end = BARREL_HALF_Z / cos(thC_ec_end);
  double r0_proj_arm_ec_end = r0_ec_end / cos(D_THETA_ENDCAP / 2);
  double ec_phislice_8pa_z1 = r0_proj_arm_ec_end * cos(thC_ec_end + D_THETA_ENDCAP / 2);
  double zT = ec_phislice_8pa_z1 - TIMING_GAP - XTAL_DEPTH_T;

  // Endcap theta-start (skip rings obscured by beampipe)
  int TLENDCAP_THETA_START = 0;
  for (int iTheta = 0; iTheta < N_THETA_TLENDCAP; iTheta++) {
    double thC = D_THETA_TLENDCAP / 2 + iTheta * D_THETA_TLENDCAP;
    double RinEndcap = zT * tan(thC);
    if (RinEndcap < BEAMPIPE_OPENING)
      TLENDCAP_THETA_START++;
  }

  // Timing barrel phi-slice envelope geometry
  double thC_tl_end = THETA_SIZE_TLBARREL - D_THETA_TLBARREL / 2;
  double thC_tl_beg = D_THETA_TLBARREL / 2;

  double r0_tl_end = rT / cos(thC_tl_end);
  double r0_tl_beg = rT / cos(thC_tl_beg);

  double r0_proj_arm_tl_end = r0_tl_end / cos(D_THETA_TLBARREL / 2);
  double r2_proj_arm_tl_beg = (r0_tl_beg + XTAL_DEPTH_T) / cos(D_THETA_TLBARREL / 2);

  double tl_phislice_8pa_z1 = r0_proj_arm_tl_end * cos(thC_tl_end + D_THETA_TLBARREL / 2);
  double tl_phislice_8pa_z2 = r2_proj_arm_tl_beg * cos(thC_tl_beg - D_THETA_TLBARREL / 2) + XTAL_DEPTH_T;
  double tl_phislice_8pa_y0 = tl_phislice_8pa_z1 * tan(thC_tl_end + D_THETA_TLBARREL / 2);
  double tl_phislice_8pa_y2 = tl_phislice_8pa_z2 * tan(thC_tl_end + D_THETA_TLBARREL / 2);

  // Timing endcap phi-slice envelope geometry
  double thC_tlec_end = THETA_SIZE_TLENDCAP - D_THETA_TLENDCAP / 2;
  double thC_tlec_beg = D_THETA_TLENDCAP / 2 + TLENDCAP_THETA_START * D_THETA_TLENDCAP;

  double r0_tlec_end = zT / cos(thC_tlec_end);
  double r0_tlec_beg = zT / cos(thC_tlec_beg);

  double r0_proj_arm_tlec_end = r0_tlec_end / cos(D_THETA_TLENDCAP / 2);
  double r2_proj_arm_tlec_beg = (r0_tlec_beg + XTAL_DEPTH_T) / cos(D_THETA_TLENDCAP / 2);

  double tlec_phislice_8pa_z1 = r0_proj_arm_tlec_end * cos(thC_tlec_end + D_THETA_TLENDCAP / 2);
  double tlec_phislice_8pa_z2 = r2_proj_arm_tlec_beg * cos(thC_tlec_beg - D_THETA_TLENDCAP / 2) + XTAL_DEPTH_T;

  double tlec_end_phislice_8pa_y0 = tlec_phislice_8pa_z1 * tan(thC_tlec_end + D_THETA_TLENDCAP / 2);
  double tlec_end_phislice_8pa_y2 = tlec_phislice_8pa_z2 * tan(thC_tlec_end + D_THETA_TLENDCAP / 2);
  double tlec_beg_phislice_8pa_y0 = tlec_phislice_8pa_z1 * tan(thC_tlec_beg - D_THETA_TLENDCAP / 2);
  double tlec_beg_phislice_8pa_y2 = tlec_phislice_8pa_z2 * tan(thC_tlec_beg - D_THETA_TLENDCAP / 2);

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "=SCEPCAL TIMING LAYER INPUTS=====================" << std::endl;
  std::cout << std::endl;
  std::cout << "BARREL_HALF_Z:        " << BARREL_HALF_Z << std::endl;
  std::cout << "BARREL_INNER_R:       " << BARREL_INNER_R << std::endl;
  std::cout << "PHI_SEGMENTS:         " << PHI_SEGMENTS << std::endl;
  std::cout << std::endl;
  std::cout << "XTAL_DEPTH_T:         " << XTAL_DEPTH_T << std::endl;
  std::cout << "XTAL_LENGTH_T:        " << XTAL_LENGTH_T << std::endl;
  std::cout << std::endl;
  std::cout << "PROJECTIVE_OFFSET_R:  " << PROJ_OFFSET_R << std::endl;
  std::cout << "BEAMPIPE_OPENING:     " << BEAMPIPE_OPENING << std::endl;
  std::cout << "TIMING_LAYER_GAP:     " << TIMING_GAP << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "=CONTROL=========================================" << std::endl;
  std::cout << std::endl;
  std::cout << "CONSTRUCT_TLBARREL:   " << CONSTRUCT_TLBARREL << std::endl;
  std::cout << "CONSTRUCT_TLENDCAP:   " << CONSTRUCT_TLENDCAP << std::endl;
  std::cout << "USE_OPTICAL_SURFACES: " << USE_OPTICAL_SURFACES << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "=CALCULATED PARAMETERS===========================" << std::endl;
  std::cout << std::endl;
  std::cout << "N_THETA_TLBARREL:     " << N_THETA_TLBARREL << std::endl;
  std::cout << "N_GAMMA_TIMING:       " << N_GAMMA_TLBARREL << std::endl;
  std::cout << std::endl;
  std::cout << "N_THETA_TLENDCAP:     " << N_THETA_TLENDCAP << std::endl;
  std::cout << "TLENDCAP_THETA_START: " << TLENDCAP_THETA_START << std::endl;
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
  auto segmentation = dynamic_cast<dd4hep::DDSegmentation::SCEPCal_TimingSegmentation_k4geo*>(_geoSeg->segmentation());

  dd4hep::OpticalSurfaceManager surfMgr = theDetector.surfaceManager();
  dd4hep::OpticalSurface LYSO_to_ESR = surfMgr.opticalSurface("/world/" + detName + "#LYSO_to_ESR");

  // Global assembly volumes
  std::vector<double> zTimingPolyhedra = {-tl_phislice_8pa_y2, -tl_phislice_8pa_y0, tl_phislice_8pa_y0,
                                          tl_phislice_8pa_y2};
  std::vector<double> rminTimingPolyhedra = {tl_phislice_8pa_z2, tl_phislice_8pa_z1, tl_phislice_8pa_z1,
                                             tl_phislice_8pa_z2};
  std::vector<double> rmaxTimingPolyhedra = {tl_phislice_8pa_z2, tl_phislice_8pa_z2, tl_phislice_8pa_z2,
                                             tl_phislice_8pa_z2};

  std::vector<double> zTlEndcapPolyhedra = {tlec_phislice_8pa_z1, tlec_phislice_8pa_z2};
  std::vector<double> rminTlEndcapPolyhedra = {tlec_beg_phislice_8pa_y0, tlec_beg_phislice_8pa_y2};
  std::vector<double> rmaxTlEndcapPolyhedra = {tlec_end_phislice_8pa_y0, tlec_end_phislice_8pa_y2};

  std::vector<double> zTlEndcapPolyhedra_1 = {-tlec_phislice_8pa_z2, -tlec_phislice_8pa_z1};
  std::vector<double> rminTlEndcapPolyhedra_1 = {tlec_beg_phislice_8pa_y2, tlec_beg_phislice_8pa_y0};
  std::vector<double> rmaxTlEndcapPolyhedra_1 = {tlec_end_phislice_8pa_y2, tlec_end_phislice_8pa_y0};

  dd4hep::Polyhedra timingGlobalAssemblyShape(PHI_SEGMENTS, D_PHI_GLOBAL / 2, 2 * M_PI, zTimingPolyhedra,
                                              rminTimingPolyhedra, rmaxTimingPolyhedra);
  dd4hep::Volume timingGlobalAssemblyVol("timingGlobalAssemblyVol", timingGlobalAssemblyShape,
                                         theDetector.material("Vacuum"));
  timingGlobalAssemblyVol.setVisAttributes(theDetector, tlbarrelAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume tlbarrelAssemblyPlacedVol = experimentalHall.placeVolume(timingGlobalAssemblyVol);
  tlbarrelAssemblyPlacedVol.addPhysVolID("system", TLBARREL_SYSTEM_NO);

  dd4hep::Polyhedra tlendcapGlobalAssemblyShape(PHI_SEGMENTS, D_PHI_GLOBAL / 2, 2 * M_PI, zTlEndcapPolyhedra,
                                                rminTlEndcapPolyhedra, rmaxTlEndcapPolyhedra);
  dd4hep::Volume tlendcapGlobalAssemblyVol("tlendcapGlobalAssemblyVol", tlendcapGlobalAssemblyShape,
                                           theDetector.material("Vacuum"));
  tlendcapGlobalAssemblyVol.setVisAttributes(theDetector, tlendcapAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume tlendcapAssemblyPlacedVol = experimentalHall.placeVolume(tlendcapGlobalAssemblyVol);
  tlendcapAssemblyPlacedVol.addPhysVolID("system", TLENDCAP_SYSTEM_NO);

  dd4hep::Polyhedra tlendcapGlobalAssemblyShape_1(PHI_SEGMENTS, D_PHI_GLOBAL / 2, 2 * M_PI, zTlEndcapPolyhedra_1,
                                                  rminTlEndcapPolyhedra_1, rmaxTlEndcapPolyhedra_1);
  dd4hep::Volume tlendcapGlobalAssemblyVol_1("tlendcapGlobalAssemblyVol_1", tlendcapGlobalAssemblyShape_1,
                                             theDetector.material("Vacuum"));
  tlendcapGlobalAssemblyVol_1.setVisAttributes(theDetector, tlendcapAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume tlendcapAssemblyPlacedVol_1 = experimentalHall.placeVolume(tlendcapGlobalAssemblyVol_1);
  tlendcapAssemblyPlacedVol_1.addPhysVolID("system", TLENDCAP_SYSTEM_NO);
  tlendcapAssemblyPlacedVol_1.addPhysVolID("theta", 1);

  ScepcalDetElement.setPlacement(tlbarrelAssemblyPlacedVol);

  // Pre-build phi rotation vectors
  std::vector<RotationZ> phi_barrel_rotations, phi_endcap_rotations;
  phi_barrel_rotations.reserve(TIMING_PHI_END - TIMING_PHI_START);
  phi_endcap_rotations.reserve(TLENDCAP_PHI_END - TLENDCAP_PHI_START);
  for (int iPhi = TIMING_PHI_START; iPhi < TIMING_PHI_END; iPhi++)
    phi_barrel_rotations.push_back(RotationZ(iPhi * D_PHI_GLOBAL));
  for (int iPhi = TLENDCAP_PHI_START; iPhi < TLENDCAP_PHI_END; iPhi++)
    phi_endcap_rotations.push_back(RotationZ(iPhi * D_PHI_GLOBAL));

  // Pack shared parameters
  SCEPCalTimingParams p;
  p.BARREL_HALF_Z = BARREL_HALF_Z;
  p.BARREL_INNER_R = BARREL_INNER_R;
  p.XTAL_DEPTH_T = XTAL_DEPTH_T;
  p.XTAL_LENGTH_T = XTAL_LENGTH_T;
  p.XTAL_TH_WIDTH = XTAL_TH_WIDTH;
  p.PHI_SEGMENTS = PHI_SEGMENTS;
  p.D_PHI_GLOBAL = D_PHI_GLOBAL;

  p.rT = rT;
  p.N_THETA_TLBARREL = N_THETA_TLBARREL;
  p.N_GAMMA_TLBARREL = N_GAMMA_TLBARREL;
  p.D_THETA_TLBARREL = D_THETA_TLBARREL;
  p.D_GAMMA_TLBARREL = D_GAMMA_TLBARREL;
  p.THETA_SIZE_TLBARREL = THETA_SIZE_TLBARREL;
  p.THETA_SIZE_TLENDCAP = THETA_SIZE_TLENDCAP;
  p.TLBARREL_SYSTEM_NO = TLBARREL_SYSTEM_NO;
  p.TIMING_PHI_START = TIMING_PHI_START;
  p.TIMING_PHI_END = TIMING_PHI_END;
  p.CONSTRUCT_TLBARREL = CONSTRUCT_TLBARREL;

  p.zT = zT;
  p.N_THETA_TLENDCAP = N_THETA_TLENDCAP;
  p.TLENDCAP_THETA_START = TLENDCAP_THETA_START;
  p.D_THETA_TLENDCAP = D_THETA_TLENDCAP;
  p.TLENDCAP_SYSTEM_NO = TLENDCAP_SYSTEM_NO;
  p.TLENDCAP_PHI_START = TLENDCAP_PHI_START;
  p.TLENDCAP_PHI_END = TLENDCAP_PHI_END;
  p.CONSTRUCT_TLENDCAP = CONSTRUCT_TLENDCAP;

  p.crystalTMaterial = crystalTXML.materialStr();
  p.crystalTVis = crystalTXML.visStr();
  p.tlbarrelAssemblyPhiVis = tlbarrelAssemblyPhiVisXML.visStr();
  p.tlbarrelAssemblyThetaVis = tlbarrelAssemblyThetaVisXML.visStr();
  p.tlendcapAssemblyPhiVis = tlendcapAssemblyPhiVisXML.visStr();
  p.tlendcapAssemblyThetaVis = tlendcapAssemblyThetaVisXML.visStr();

  p.zTimingPolyhedra = zTimingPolyhedra;
  p.rminTimingPolyhedra = rminTimingPolyhedra;
  p.rmaxTimingPolyhedra = rmaxTimingPolyhedra;

  p.zTlEndcapPolyhedra = zTlEndcapPolyhedra;
  p.rminTlEndcapPolyhedra = rminTlEndcapPolyhedra;
  p.rmaxTlEndcapPolyhedra = rmaxTlEndcapPolyhedra;
  p.zTlEndcapPolyhedra_1 = zTlEndcapPolyhedra_1;
  p.rminTlEndcapPolyhedra_1 = rminTlEndcapPolyhedra_1;
  p.rmaxTlEndcapPolyhedra_1 = rmaxTlEndcapPolyhedra_1;

  int numCrystalsTiming = 0;
  if (CONSTRUCT_TLBARREL)
    numCrystalsTiming = buildTimingBarrel(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES,
                                          LYSO_to_ESR, p, phi_barrel_rotations, timingGlobalAssemblyVol);

  int numCrystalsTlEndcap = 0;
  if (CONSTRUCT_TLENDCAP)
    numCrystalsTlEndcap =
        buildTimingEndcap(theDetector, sens, segmentation, ScepcalDetElement, USE_OPTICAL_SURFACES, LYSO_to_ESR, p,
                          phi_endcap_rotations, tlendcapGlobalAssemblyVol, tlendcapGlobalAssemblyVol_1);

  std::cout << std::endl;
  std::cout << "NUM_CRYSTALS_TLBARREL:" << numCrystalsTiming << std::endl;
  std::cout << "NUM_CRYSTALS_TLENDCAP:" << numCrystalsTlEndcap << std::endl;
  std::cout << std::endl;

  return ScepcalDetElement;
}

DECLARE_DETELEMENT(SCEPCal_TimingLayer, create_detector_SCEPCal_TimingLayer)

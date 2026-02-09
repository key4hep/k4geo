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

using dd4hep::Position;
using dd4hep::RotationZYX;
using dd4hep::Transform3D;
using ROOT::Math::RotationY;
using ROOT::Math::RotationZ;
using ROOT::Math::XYZVector;
using namespace dd4hep;

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

  // Start of partial visualisation parameters
  const int PHI_LOAD_START = dimXML.attr<int>(_Unicode(phi_load_start));
  const int PHI_LOAD_END = dimXML.attr<int>(_Unicode(phi_load_end));
  const int THETA_LOAD_START = dimXML.attr<int>(_Unicode(theta_load_start));
  const int THETA_LOAD_END = dimXML.attr<int>(_Unicode(theta_load_end));
  const int GAMMA_LOAD_START = dimXML.attr<int>(_Unicode(gamma_load_start));
  const int GAMMA_LOAD_END = dimXML.attr<int>(_Unicode(gamma_load_end));
  // End of partial visualisation parameters

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
    double RinEndcap = (BARREL_HALF_Z)*tan(thC);
    if (RinEndcap < (BEAMPIPE_OPENING + PROJ_OFFSET_R))
      ENDCAP_THETA_START++;
  }

  // Barrel Phi envelope with protrusions for tilted crystals
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

  // Endcap Phi envelope with protrusions for tilted crystals
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

  int numCrystalsBarrel = 0;
  int numCrystalsEndcap = 0;

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

  // Lambda for crystals
  auto CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId =
      [&theDetector, &sens, &segmentation, &ScepcalDetElement, &USE_OPTICAL_SURFACES,
       &PbWO4_to_ESR](std::string volName, double dz, const std::array<double, 16> vertices, xml_comp_t compXml,
                      dd4hep::Transform3D transform, dd4hep::Volume assemblyVol, int nSystem, int nPhi, int nTheta,
                      int nGamma, int nEpsilon, int nDepth, XYZVector posGlobal) {
        dd4hep::EightPointSolid theShape(dz, vertices.data());
        dd4hep::Volume theVolume(volName, theShape, theDetector.material(compXml.materialStr()));
        theVolume.setVisAttributes(theDetector, compXml.visStr());

        if (USE_OPTICAL_SURFACES) {
          dd4hep::SkinSurface(theDetector, ScepcalDetElement,
                              volName + "Surface_" + std::to_string(nPhi) + "_" + std::to_string(nTheta) + "_" +
                                  std::to_string(nGamma) + "_" + std::to_string(nDepth),
                              PbWO4_to_ESR, theVolume);
        }

        theVolume.setSensitiveDetector(sens);

        auto volID = segmentation->setVolumeID(nSystem, nPhi, nTheta, nGamma, nEpsilon, nDepth);
        int volID_32 = segmentation->getFirst32bits(volID);
        dd4hep::PlacedVolume thePlacedVol = assemblyVol.placeVolume(theVolume, volID_32, transform);

        thePlacedVol.addPhysVolID("system", nSystem);
        thePlacedVol.addPhysVolID("phi", nPhi);
        thePlacedVol.addPhysVolID("theta", nTheta);
        thePlacedVol.addPhysVolID("gamma", nGamma);
        thePlacedVol.addPhysVolID("epsilon", nEpsilon);
        thePlacedVol.addPhysVolID("depth", nDepth);

        segmentation->savePosition(volID_32, posGlobal);
      };

  // Lambdas for crystal tower divisions
  auto bilinearInterpolateTower = [&](const double* vf, double u, double v) {
    double x0 = (1 - u) * (1 - v) * vf[0] + u * (1 - v) * vf[2] + u * v * vf[4] + (1 - u) * v * vf[6];

    double y0 = (1 - u) * (1 - v) * vf[1] + u * (1 - v) * vf[3] + u * v * vf[5] + (1 - u) * v * vf[7];

    double x1 = (1 - u) * (1 - v) * vf[8] + u * (1 - v) * vf[10] + u * v * vf[12] + (1 - u) * v * vf[14];

    double y1 = (1 - u) * (1 - v) * vf[9] + u * (1 - v) * vf[11] + u * v * vf[13] + (1 - u) * v * vf[15];

    return std::array<double, 4>{x0, y0, x1, y1};
  };

  auto getSingleCrystalVertices = [&](int i, int j, int xtalDiv, const double* vertices,
                                      bool reflected) -> std::array<double, 16> {
    double u0 = double(i) / double(xtalDiv);
    double u1 = double(i + 1) / double(xtalDiv);
    double v0 = double(j) / double(xtalDiv);
    double v1 = double(j + 1) / double(xtalDiv);

    auto P00 = bilinearInterpolateTower(vertices, u0, v0);
    auto P10 = bilinearInterpolateTower(vertices, u1, v0);
    auto P11 = bilinearInterpolateTower(vertices, u1, v1);
    auto P01 = bilinearInterpolateTower(vertices, u0, v1);

    std::array<double, 16> verticesSub = {P00[0], P00[1], P10[0], P10[1], P11[0], P11[1], P01[0], P01[1],
                                          P00[2], P00[3], P10[2], P10[3], P11[2], P11[3], P01[2], P01[3]};

    std::array<double, 16> verticesSub_1 = {P00[2], P00[3], P10[2], P10[3], P11[2], P11[3], P01[2], P01[3],
                                            P00[0], P00[1], P10[0], P10[1], P11[0], P11[1], P01[0], P01[1]};

    if (reflected)
      return verticesSub_1;
    return verticesSub;
  };

  auto getSingleCrystalCenter = [&](const std::array<double, 16> vSub) -> std::array<double, 2> {
    std::array<double, 2> center = {
        0.125 * (vSub[0] + vSub[2] + vSub[4] + vSub[6] + vSub[8] + vSub[10] + vSub[12] + vSub[14]),
        0.125 * (vSub[1] + vSub[3] + vSub[5] + vSub[7] + vSub[9] + vSub[11] + vSub[13] + vSub[15])};
    return center;
  };

  //////////////////////////////
  // Barrel
  // Phi - divide barrel into phi slices
  // Theta - divide phi slice into slabs in theta
  // Gamma - divide a theta slab into square NxN crystal towers in phi
  // Epsilon - index of single crystal in a NxN tower
  // Projective offsets - tilt the rear face of the trapezoids to point at specified offset away from IP
  //////////////////////////////

  for (int iPhi = CONSTRUCT_BARREL ? BARREL_PHI_START : BARREL_PHI_END; iPhi < BARREL_PHI_END; iPhi++) {
    double phiGlobal = iPhi * D_PHI_GLOBAL;
    RotationZ rotZphiGlobal(phiGlobal);
    dd4hep::Polyhedra barrelPhiAssemblyShape(1, -D_PHI_GLOBAL / 2, D_PHI_GLOBAL, zBarrelPolyhedra, rminBarrelPolyhedra,
                                             rmaxBarrelPolyhedra);
    dd4hep::Volume barrelPhiAssemblyVolume("barrelPhiAssembly", barrelPhiAssemblyShape, theDetector.material("Vacuum"));
    barrelPhiAssemblyVolume.setVisAttributes(theDetector, barrelAssemblyPhiVisXML.visStr());
    if (PHI_LOAD_START <= PHI_LOAD_END) {
      if ((iPhi >= PHI_LOAD_START) && (iPhi <= PHI_LOAD_END)) {
        barrelGlobalAssemblyVol.placeVolume(barrelPhiAssemblyVolume, Transform3D(rotZphiGlobal));
      }
    } else if (PHI_LOAD_START > PHI_LOAD_END) {
      if ((iPhi >= PHI_LOAD_START) || (iPhi <= PHI_LOAD_END)) {
        barrelGlobalAssemblyVol.placeVolume(barrelPhiAssemblyVolume, Transform3D(rotZphiGlobal));
      }
    }

    for (int iTheta = 0; iTheta < N_THETA_BARREL; iTheta++) {
      double thC = THETA_SIZE_ENDCAP + D_THETA_BARREL / 2 + (iTheta * D_THETA_BARREL);
      RotationY rotYthGlobal(thC);

      double r0e = (BARREL_INNER_R + PROJ_OFFSET_R) / sin(thC);
      double r1e = r0e + XTAL_LEN_F;
      double r2e = r1e + XTAL_LEN_R;
      double y0e = r0e * tan(D_THETA_BARREL / 2.);
      double y1e = r1e * tan(D_THETA_BARREL / 2.);
      double y2e = r2e * tan(D_THETA_BARREL / 2.);

      double x0y0 = r0e * sin(thC) - y0e * cos(thC) - PROJ_OFFSET_R;
      double x1y0 = r0e * sin(thC) + y0e * cos(thC) - PROJ_OFFSET_R;

      double x0y1 = r1e * sin(thC) - y1e * cos(thC) - PROJ_OFFSET_R;
      double x1y1 = r1e * sin(thC) + y1e * cos(thC) - PROJ_OFFSET_R;

      double x0y2 = r2e * sin(thC) - y2e * cos(thC) - PROJ_OFFSET_R;
      double x1y2 = r2e * sin(thC) + y2e * cos(thC) - PROJ_OFFSET_R;

      double x0y0l_E = x0y0 * tan(-D_PHI_GLOBAL / 2);
      double x0y0r_E = x0y0 * tan(D_PHI_GLOBAL / 2);
      double x1y0l_E = x1y0 * tan(-D_PHI_GLOBAL / 2);
      double x1y0r_E = x1y0 * tan(D_PHI_GLOBAL / 2);

      double x0y2l_E = x0y2 * tan(-D_PHI_GLOBAL / 2);
      double x0y2r_E = x0y2 * tan(D_PHI_GLOBAL / 2);
      double x1y2l_E = x1y2 * tan(-D_PHI_GLOBAL / 2);
      double x1y2r_E = x1y2 * tan(D_PHI_GLOBAL / 2);

      double verticesE[] = {x0y0r_E, y0e, x1y0r_E, -y0e, x1y0l_E, -y0e, x0y0l_E, y0e,
                            x0y2r_E, y2e, x1y2r_E, -y2e, x1y2l_E, -y2e, x0y2l_E, y2e};

      double rE = r0e + (XTAL_LEN_F + XTAL_LEN_R) / 2.;
      RotationZYX rotE(M_PI / 2, thC, 0);
      Position dispE(rE * sin(thC) - PROJ_OFFSET_R, 0, rE * cos(thC));

      dd4hep::EightPointSolid barrelThetaAssemblyShape((XTAL_LEN_F + XTAL_LEN_R) / 2, verticesE);
      dd4hep::Volume barrelThetaAssemblyVolume("barrelThetaAssembly", barrelThetaAssemblyShape,
                                               theDetector.material("Vacuum"));
      barrelThetaAssemblyVolume.setVisAttributes(theDetector, barrelAssemblyThetaVisXML.visStr());

      if ((iTheta >= THETA_LOAD_START) && (iTheta <= THETA_LOAD_END)) {
        barrelPhiAssemblyVolume.placeVolume(barrelThetaAssemblyVolume, Transform3D(rotE, dispE));
      }

      for (int nGamma = 0; nGamma < N_GAMMA_BARREL; nGamma++) {
        double gamma = -D_PHI_GLOBAL / 2 + D_GAMMA_BARREL / 2 + D_GAMMA_BARREL * nGamma;

        double projOffsetXmax =
            std::min(r0e * tan(N_GAMMA_BARREL % 2 == 0 ? D_GAMMA_BARREL / 2 : D_GAMMA_BARREL / 2), PROJ_OFFSET_X);

        double r1_x_gamma_shift = projOffsetXmax * (r1e - r0e) / r0e;
        double r2_x_gamma_shift = projOffsetXmax * (r2e - r0e) / r0e;

        int left = 0;
        int right = 0;

        if (N_GAMMA_BARREL % 2 == 0) {
          if (nGamma < (N_GAMMA_BARREL / 2 - 1)) {
            left = 1;
            right = 1;
          }
          if (nGamma == (N_GAMMA_BARREL / 2 - 1)) {
            left = 1;
            right = 0;
          }
          if (nGamma == N_GAMMA_BARREL / 2) {
            left = 0;
            right = -1;
          }
          if (nGamma > N_GAMMA_BARREL / 2) {
            left = -1;
            right = -1;
          }
        } else {
          if (nGamma < N_GAMMA_BARREL / 2) {
            left = 1;
            right = 1;
          }
          if (nGamma == N_GAMMA_BARREL / 2) {
            left = 1;
            right = -1;
          }
          if (nGamma > N_GAMMA_BARREL / 2) {
            left = -1;
            right = -1;
          }
        }

        if (nGamma == 0) {
          left = 0;
        }
        if (nGamma == N_GAMMA_BARREL - 1) {
          right = 0;
        }

        double x0y0l = x0y0 * tan(gamma - D_GAMMA_BARREL / 2);
        double x0y0r = x0y0 * tan(gamma + D_GAMMA_BARREL / 2);
        double x1y0l = x1y0 * tan(gamma - D_GAMMA_BARREL / 2);
        double x1y0r = x1y0 * tan(gamma + D_GAMMA_BARREL / 2);

        double x0y1l = x0y1 * tan(gamma - D_GAMMA_BARREL / 2) + left * r1_x_gamma_shift;
        double x0y1r = x0y1 * tan(gamma + D_GAMMA_BARREL / 2) + right * r1_x_gamma_shift;
        double x1y1l = x1y1 * tan(gamma - D_GAMMA_BARREL / 2) + left * r1_x_gamma_shift;
        double x1y1r = x1y1 * tan(gamma + D_GAMMA_BARREL / 2) + right * r1_x_gamma_shift;

        double x0y2l = x0y2 * tan(gamma - D_GAMMA_BARREL / 2) + left * r2_x_gamma_shift;
        double x0y2r = x0y2 * tan(gamma + D_GAMMA_BARREL / 2) + right * r2_x_gamma_shift;
        double x1y2l = x1y2 * tan(gamma - D_GAMMA_BARREL / 2) + left * r2_x_gamma_shift;
        double x1y2r = x1y2 * tan(gamma + D_GAMMA_BARREL / 2) + right * r2_x_gamma_shift;

        double verticesF[] = {x0y0r, y0e, x1y0r, -y0e, x1y0l, -y0e, x0y0l, y0e,
                              x0y1r, y1e, x1y1r, -y1e, x1y1l, -y1e, x0y1l, y1e};
        double verticesR[] = {x0y1r, y1e, x1y1r, -y1e, x1y1l, -y1e, x0y1l, y1e,
                              x0y2r, y2e, x1y2r, -y2e, x1y2l, -y2e, x0y2l, y2e};

        for (int i = 0; i < XTAL_DIV_F; i++) {
          for (int j = 0; j < XTAL_DIV_F; j++) {
            int nEpsilon = i * XTAL_DIV_F + j;

            auto vFsub = getSingleCrystalVertices(i, j, XTAL_DIV_F, verticesF, false);
            auto center = getSingleCrystalCenter(vFsub);

            Position dispFsub(0, 0, -XTAL_LEN_R / 2);

            double rGlobal = r0e + XTAL_LEN_F / 2;
            XYZVector dispGlobal(-center[1], center[0], rGlobal);
            XYZVector posGlobal = (rotZphiGlobal * (rotYthGlobal * dispGlobal + DISP_PROJ_R));

            if ((nGamma >= GAMMA_LOAD_START) && (nGamma <= GAMMA_LOAD_END)) {

              CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
                  "BarrelCrystalF", XTAL_LEN_F / 2, vFsub, crystalFXML, Transform3D(dispFsub),
                  barrelThetaAssemblyVolume, BARREL_SYSTEM_NO, iPhi, N_THETA_ENDCAP + iTheta, nGamma, nEpsilon, 0,
                  posGlobal);
              numCrystalsBarrel += 1;
            }
          }
        }

        for (int i = 0; i < XTAL_DIV_R; i++) {
          for (int j = 0; j < XTAL_DIV_R; j++) {
            int nEpsilon = i * XTAL_DIV_R + j;

            auto vRsub = getSingleCrystalVertices(i, j, XTAL_DIV_R, verticesR, false);
            auto center = getSingleCrystalCenter(vRsub);

            Position dispRsub(0, 0, XTAL_LEN_F / 2);

            double rGlobal = r0e + XTAL_LEN_F + XTAL_LEN_R / 2;
            XYZVector dispGlobal(-center[1], center[0], rGlobal);
            XYZVector posGlobal = (rotZphiGlobal * (rotYthGlobal * dispGlobal + DISP_PROJ_R));

            if ((nGamma >= GAMMA_LOAD_START) && (nGamma <= GAMMA_LOAD_END)) {
              CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
                  "BarrelCrystalR", XTAL_LEN_R / 2, vRsub, crystalRXML, Transform3D(dispRsub),
                  barrelThetaAssemblyVolume, BARREL_SYSTEM_NO, iPhi, N_THETA_ENDCAP + iTheta, nGamma, nEpsilon, 1,
                  posGlobal);
              numCrystalsBarrel += 1;
            }
          }
        }
      }
    }
  }

  //////////////////////////////
  // Endcap
  // Phi - divide endcap into phi slices
  // Theta - divide phi slice into slabs in theta
  // Gamma - divide a theta slab into square NxN crystal towers in phi
  // Epsilon - index of single crystal in a NxN tower
  // Projective offsets - tilt the rear face of the trapezoids to point at specified offset away from IP
  //////////////////////////////

  for (int iPhi = CONSTRUCT_ENDCAP ? ENDCAP_PHI_START : ENDCAP_PHI_END; iPhi < ENDCAP_PHI_END; iPhi++) {
    double phiGlobal = iPhi * D_PHI_GLOBAL;
    RotationZ rotZphiGlobal(phiGlobal);

    dd4hep::Polyhedra endcapPhiAssemblyShape(1, -D_PHI_GLOBAL / 2, D_PHI_GLOBAL, zEndcapPolyhedra, rminEndcapPolyhedra,
                                             rmaxEndcapPolyhedra);
    dd4hep::Polyhedra endcapPhiAssemblyShape_1(1, -D_PHI_GLOBAL / 2, D_PHI_GLOBAL, zEndcapPolyhedra_1,
                                               rminEndcapPolyhedra_1, rmaxEndcapPolyhedra_1);

    dd4hep::Volume endcapPhiAssemblyVolume("endcapPhiVol", endcapPhiAssemblyShape, theDetector.material("Vacuum"));
    endcapPhiAssemblyVolume.setVisAttributes(theDetector, endcapAssemblyPhiVisXML.visStr());
    dd4hep::Volume endcapPhiAssemblyVolume_1("endcapPhiVol_1", endcapPhiAssemblyShape_1,
                                             theDetector.material("Vacuum"));
    endcapPhiAssemblyVolume_1.setVisAttributes(theDetector, endcapAssemblyPhiVisXML.visStr());

    if (PHI_LOAD_START <= PHI_LOAD_END) {
      if ((iPhi >= PHI_LOAD_START) && (iPhi <= PHI_LOAD_END)) {
        endcapGlobalAssemblyVol.placeVolume(endcapPhiAssemblyVolume, Transform3D(rotZphiGlobal));
        endcapGlobalAssemblyVol_1.placeVolume(endcapPhiAssemblyVolume_1, Transform3D(rotZphiGlobal));
      }
    } else if (PHI_LOAD_START > PHI_LOAD_END) {
      if ((iPhi >= PHI_LOAD_START) || (iPhi <= PHI_LOAD_END)) {
        endcapGlobalAssemblyVol.placeVolume(endcapPhiAssemblyVolume, Transform3D(rotZphiGlobal));
        endcapGlobalAssemblyVol_1.placeVolume(endcapPhiAssemblyVolume_1, Transform3D(rotZphiGlobal));
      }
    }

    for (int iTheta = ENDCAP_THETA_START; iTheta < N_THETA_ENDCAP; iTheta++) {

      double thC = D_THETA_ENDCAP / 2 + iTheta * D_THETA_ENDCAP;
      RotationY rotYthGlobal(thC);
      RotationY rotYthGlobal_1(-thC);

      double RinEndcap = (BARREL_HALF_Z)*tan(thC);

      int nGammaEndcap = std::max(int(2 * M_PI * RinEndcap / (PHI_SEGMENTS * XTAL_TH_WIDTH)), 1);
      double dGammaEndcap = D_PHI_GLOBAL / nGammaEndcap;

      double r0e = RinEndcap / sin(thC);
      double r1e = r0e + XTAL_LEN_F;
      double r2e = r1e + XTAL_LEN_R;
      double y0e = r0e * tan(D_THETA_ENDCAP / 2.);
      double y1e = r1e * tan(D_THETA_ENDCAP / 2.);
      double y2e = r2e * tan(D_THETA_ENDCAP / 2.);

      double x0y0 = r0e * sin(thC) - y0e * cos(thC) - PROJ_OFFSET_R;
      double x1y0 = r0e * sin(thC) + y0e * cos(thC) - PROJ_OFFSET_R;

      double x0y1 = r1e * sin(thC) - y1e * cos(thC) - PROJ_OFFSET_R;
      double x1y1 = r1e * sin(thC) + y1e * cos(thC) - PROJ_OFFSET_R;

      double x0y2 = r2e * sin(thC) - y2e * cos(thC) - PROJ_OFFSET_R;
      double x1y2 = r2e * sin(thC) + y2e * cos(thC) - PROJ_OFFSET_R;

      double x0y0l_E = x0y0 * tan(-D_PHI_GLOBAL / 2);
      double x0y0r_E = x0y0 * tan(D_PHI_GLOBAL / 2);
      double x1y0l_E = x1y0 * tan(-D_PHI_GLOBAL / 2);
      double x1y0r_E = x1y0 * tan(D_PHI_GLOBAL / 2);

      double x0y2l_E = x0y2 * tan(-D_PHI_GLOBAL / 2);
      double x0y2r_E = x0y2 * tan(D_PHI_GLOBAL / 2);
      double x1y2l_E = x1y2 * tan(-D_PHI_GLOBAL / 2);
      double x1y2r_E = x1y2 * tan(D_PHI_GLOBAL / 2);

      double verticesE[] = {x0y0r_E, y0e, x1y0r_E, -y0e, x1y0l_E, -y0e, x0y0l_E, y0e,
                            x0y2r_E, y2e, x1y2r_E, -y2e, x1y2l_E, -y2e, x0y2l_E, y2e};
      double verticesE_1[] = {x0y2r_E, y2e, x1y2r_E, -y2e, x1y2l_E, -y2e, x0y2l_E, y2e,
                              x0y0r_E, y0e, x1y0r_E, -y0e, x1y0l_E, -y0e, x0y0l_E, y0e};

      double rE = r0e + (XTAL_LEN_F + XTAL_LEN_R) / 2.;

      RotationZYX rotE(M_PI / 2, thC, 0);
      RotationZYX rotE_1(M_PI / 2, -thC, 0);

      Position dispE(rE * sin(thC) - PROJ_OFFSET_R, 0, rE * cos(thC));
      Position dispE_1(rE * sin(thC) - PROJ_OFFSET_R, 0, -rE * cos(thC));

      dd4hep::EightPointSolid endcapThetaAssemblyShape((XTAL_LEN_F + XTAL_LEN_R) / 2, verticesE);
      dd4hep::EightPointSolid endcapThetaAssemblyShape_1((XTAL_LEN_F + XTAL_LEN_R) / 2, verticesE_1);

      dd4hep::Volume endcapThetaAssemblyVolume("endcapThetaAssembly", endcapThetaAssemblyShape,
                                               theDetector.material("Vacuum"));
      endcapThetaAssemblyVolume.setVisAttributes(theDetector, endcapAssemblyThetaVisXML.visStr());
      endcapPhiAssemblyVolume.placeVolume(endcapThetaAssemblyVolume, Transform3D(rotE, dispE));
      dd4hep::Volume endcapThetaAssemblyVolume_1("endcapThetaAssembly_1", endcapThetaAssemblyShape_1,
                                                 theDetector.material("Vacuum"));
      endcapThetaAssemblyVolume_1.setVisAttributes(theDetector, endcapAssemblyThetaVisXML.visStr());
      endcapPhiAssemblyVolume_1.placeVolume(endcapThetaAssemblyVolume_1, Transform3D(rotE_1, dispE_1));

      for (int nGamma = 0; nGamma < nGammaEndcap; nGamma++) {
        double gamma = -D_PHI_GLOBAL / 2 + dGammaEndcap / 2 + dGammaEndcap * nGamma;

        double projOffsetXmax =
            std::min(r0e * tan(nGammaEndcap % 2 == 0 ? dGammaEndcap / 2 : dGammaEndcap / 2), PROJ_OFFSET_X);

        double r1_x_gamma_shift = projOffsetXmax * (r1e - r0e) / r0e;
        double r2_x_gamma_shift = projOffsetXmax * (r2e - r0e) / r0e;

        int left = 0;
        int right = 0;

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

        if (nGamma == 0) {
          left = 0;
        }
        if (nGamma == nGammaEndcap - 1) {
          right = 0;
        }

        double x0y0l = x0y0 * tan(gamma - dGammaEndcap / 2);
        double x0y0r = x0y0 * tan(gamma + dGammaEndcap / 2);
        double x1y0l = x1y0 * tan(gamma - dGammaEndcap / 2);
        double x1y0r = x1y0 * tan(gamma + dGammaEndcap / 2);

        double x0y1l = x0y1 * tan(gamma - dGammaEndcap / 2) + left * r1_x_gamma_shift;
        double x0y1r = x0y1 * tan(gamma + dGammaEndcap / 2) + right * r1_x_gamma_shift;
        double x1y1l = x1y1 * tan(gamma - dGammaEndcap / 2) + left * r1_x_gamma_shift;
        double x1y1r = x1y1 * tan(gamma + dGammaEndcap / 2) + right * r1_x_gamma_shift;

        double x0y2l = x0y2 * tan(gamma - dGammaEndcap / 2) + left * r2_x_gamma_shift;
        double x0y2r = x0y2 * tan(gamma + dGammaEndcap / 2) + right * r2_x_gamma_shift;
        double x1y2l = x1y2 * tan(gamma - dGammaEndcap / 2) + left * r2_x_gamma_shift;
        double x1y2r = x1y2 * tan(gamma + dGammaEndcap / 2) + right * r2_x_gamma_shift;

        double verticesF[] = {x0y0r, y0e, x1y0r, -y0e, x1y0l, -y0e, x0y0l, y0e,
                              x0y1r, y1e, x1y1r, -y1e, x1y1l, -y1e, x0y1l, y1e};
        double verticesR[] = {x0y1r, y1e, x1y1r, -y1e, x1y1l, -y1e, x0y1l, y1e,
                              x0y2r, y2e, x1y2r, -y2e, x1y2l, -y2e, x0y2l, y2e};

        for (int i = 0; i < XTAL_DIV_F; i++) {
          for (int j = 0; j < XTAL_DIV_F; j++) {
            int nEpsilon = i * XTAL_DIV_F + j;

            auto vFsub = getSingleCrystalVertices(i, j, XTAL_DIV_F, verticesF, false);
            auto center = getSingleCrystalCenter(vFsub);
            Position dispFsub(0, 0, -XTAL_LEN_R / 2);

            auto vFsub_1 = getSingleCrystalVertices(i, j, XTAL_DIV_F, verticesF, true);
            auto center_1 = getSingleCrystalCenter(vFsub_1);
            Position dispFsub_1(0, 0, XTAL_LEN_R / 2);

            double rGlobal = r0e + XTAL_LEN_F / 2;
            XYZVector dispGlobal(-center[1], center[0], rGlobal);
            XYZVector posGlobal = (rotZphiGlobal * (rotYthGlobal * dispGlobal + DISP_PROJ_R));

            XYZVector dispGlobal_1(-center_1[1], center_1[0], -rGlobal);
            XYZVector posGlobal_1 = (rotZphiGlobal * (rotYthGlobal_1 * dispGlobal_1 + DISP_PROJ_R));

            CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
                "EndcapCrystalF", XTAL_LEN_F / 2, vFsub, crystalFXML, Transform3D(dispFsub), endcapThetaAssemblyVolume,
                ENDCAP_SYSTEM_NO, iPhi, iTheta, nGamma, nEpsilon, 0, posGlobal);
            numCrystalsEndcap += 1;

            CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
                "EndcapCrystalF_1", XTAL_LEN_F / 2, vFsub_1, crystalFXML, Transform3D(dispFsub_1),
                endcapThetaAssemblyVolume_1, ENDCAP_SYSTEM_NO, iPhi, 2 * N_THETA_ENDCAP + N_THETA_BARREL - iTheta,
                nGamma, nEpsilon, 0, posGlobal_1);
            numCrystalsEndcap += 1;
          }
        }

        for (int i = 0; i < XTAL_DIV_R; i++) {
          for (int j = 0; j < XTAL_DIV_R; j++) {
            int nEpsilon = i * XTAL_DIV_R + j;

            auto vRsub = getSingleCrystalVertices(i, j, XTAL_DIV_R, verticesR, false);
            auto center = getSingleCrystalCenter(vRsub);
            Position dispRsub(0, 0, XTAL_LEN_F / 2);

            auto vRsub_1 = getSingleCrystalVertices(i, j, XTAL_DIV_R, verticesR, true);
            auto center_1 = getSingleCrystalCenter(vRsub_1);
            Position dispRsub_1(0, 0, -XTAL_LEN_F / 2);

            double rGlobal = r0e + XTAL_LEN_F + XTAL_LEN_R / 2;
            XYZVector dispGlobal(-center[1], center[0], rGlobal);
            XYZVector posGlobal = (rotZphiGlobal * (rotYthGlobal * dispGlobal + DISP_PROJ_R));

            XYZVector dispGlobal_1(-center_1[1], center_1[0], -rGlobal);
            XYZVector posGlobal_1 = (rotZphiGlobal * (rotYthGlobal_1 * dispGlobal_1 + DISP_PROJ_R));

            CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
                "EndcapCrystalR", XTAL_LEN_R / 2, vRsub, crystalRXML, Transform3D(dispRsub), endcapThetaAssemblyVolume,
                ENDCAP_SYSTEM_NO, iPhi, iTheta, nGamma, nEpsilon, 1, posGlobal);
            numCrystalsEndcap += 1;

            CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
                "EndcapCrystalR_1", XTAL_LEN_R / 2, vRsub_1, crystalRXML, Transform3D(dispRsub_1),
                endcapThetaAssemblyVolume_1, ENDCAP_SYSTEM_NO, iPhi, 2 * N_THETA_ENDCAP + N_THETA_BARREL - iTheta,
                nGamma, nEpsilon, 1, posGlobal_1);
            numCrystalsEndcap += 1;
          }
        }
      }
    }
  }

  std::cout << std::endl;
  std::cout << "NUM_CRYSTALS_BARREL:  " << numCrystalsBarrel << std::endl;
  std::cout << "NUM_CRYSTALS_ENDCAP:  " << numCrystalsEndcap << std::endl;
  std::cout << std::endl;

  return ScepcalDetElement;
}

DECLARE_DETELEMENT(SCEPCal_MainLayer, create_detector_SCEPCal_MainLayer)

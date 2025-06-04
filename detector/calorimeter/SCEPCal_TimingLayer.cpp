//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================

#include "detectorSegmentations/SCEPCal_TimingSegmentationHandle_k4geo.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetectorTools.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/Volumes.h"
#include "TGeoTrd2.h"
#include <bitset>

using dd4hep::Transform3D;
using dd4hep::RotationZYX;
using ROOT::Math::RotationY;
using ROOT::Math::RotationZ;
using dd4hep::Position;
using ROOT::Math::XYZVector;
using namespace dd4hep;

static dd4hep::Ref_t
create_detector_SCEPCal_TimingLayer(dd4hep::Detector &theDetector,xml_h xmlElement,dd4hep::SensitiveDetector sens) {
  xml_det_t  detectorXML                 =xmlElement;
  xml_comp_t dimXML                      =detectorXML.child(_Unicode(dim));
  xml_comp_t tlbarrelXML                 =detectorXML.child(_Unicode(tlbarrel));
  xml_comp_t tlendcapXML                 =detectorXML.child(_Unicode(tlendcap));
  xml_comp_t crystalTXML                 =detectorXML.child(_Unicode(crystalT));
  
  xml_comp_t tlbarrelAssemblyGlobalVisXML=detectorXML.child(_Unicode(tlbarrelAssemblyGlobalVis));
  xml_comp_t tlbarrelAssemblyThetaVisXML =detectorXML.child(_Unicode(tlbarrelAssemblyThetaVis));
  xml_comp_t tlbarrelAssemblyPhiVisXML   =detectorXML.child(_Unicode(tlbarrelAssemblyPhiVis));

  xml_comp_t tlendcapAssemblyGlobalVisXML=detectorXML.child(_Unicode(tlendcapAssemblyGlobalVis));
  xml_comp_t tlendcapAssemblyPhiVisXML   =detectorXML.child(_Unicode(tlendcapAssemblyPhiVis));
  xml_comp_t tlendcapAssemblyThetaVisXML =detectorXML.child(_Unicode(tlendcapAssemblyThetaVis));

  std::string detName                    =detectorXML.nameStr();

  const double  BARREL_HALF_Z            =dimXML.attr<double>(_Unicode(barrelHalfZ));
  const double  BARREL_INNER_R           =dimXML.attr<double>(_Unicode(barrelInnerR));

  const double  XTAL_TH_WIDTH            =dimXML.attr<double>(_Unicode(crystalThetaWidth));
  const double  XTAL_DEPTH_T             =dimXML.attr<double>(_Unicode(crystalTdepth));
  const double  XTAL_LENGTH_T            =dimXML.attr<double>(_Unicode(crystalTlength));

  const double  BEAMPIPE_OPENING         =dimXML.attr<double>(_Unicode(beampipe_opening));
  const double  TIMING_GAP               =dimXML.attr<double>(_Unicode(timing_gap));
  
  const int     PHI_SEGMENTS             =dimXML.attr<int>(_Unicode(phiSegments));
  const double  PROJ_OFFSET_R            =dimXML.attr<double>(_Unicode(projectiveOffsetR));

  const bool    USE_OPTICAL_SURFACES     =dimXML.attr<bool>(_Unicode(useOpticalSurfaces));
  
  const int     TLBARREL_SYSTEM_NO       =tlbarrelXML.attr<int>(_Unicode(system));
  const bool    CONSTRUCT_TLBARREL       =tlbarrelXML.attr<bool>(_Unicode(construct));
  const int     TIMING_PHI_START         =0;
  const int     TIMING_PHI_END           =PHI_SEGMENTS;

  const int     TLENDCAP_SYSTEM_NO       =tlendcapXML.attr<int>(_Unicode(system));
  const bool    CONSTRUCT_TLENDCAP       =tlendcapXML.attr<bool>(_Unicode(construct));
  const int     TLENDCAP_PHI_START       =0;
  const int     TLENDCAP_PHI_END         =PHI_SEGMENTS;

  const double  D_PHI_GLOBAL             =2*M_PI/PHI_SEGMENTS;

  double        THETA_SIZE_BARREL        =atan(BARREL_HALF_Z/(BARREL_INNER_R+PROJ_OFFSET_R));
  double        THETA_SIZE_ENDCAP        =atan((BARREL_INNER_R+PROJ_OFFSET_R)/BARREL_HALF_Z);
  
  double        THETA_SIZE_TLBARREL      =atan(BARREL_HALF_Z/(BARREL_INNER_R));
  double        THETA_SIZE_TLENDCAP      =atan(BARREL_INNER_R/BARREL_HALF_Z);
    
  int           N_THETA_BARREL           =2*floor(BARREL_HALF_Z/XTAL_TH_WIDTH);
  int           N_THETA_ENDCAP           =floor(BARREL_INNER_R/XTAL_TH_WIDTH);
  
  int           N_THETA_TLBARREL         =2*floor(BARREL_HALF_Z/XTAL_DEPTH_T);
  int           N_THETA_TLENDCAP         =floor(BARREL_INNER_R/XTAL_DEPTH_T);
  
  double        D_THETA_BARREL           =2*THETA_SIZE_BARREL/N_THETA_BARREL;
  double        D_THETA_ENDCAP           =THETA_SIZE_ENDCAP/N_THETA_ENDCAP;

  double        D_THETA_TLBARREL         =2*THETA_SIZE_TLBARREL/N_THETA_TLBARREL;
  double        D_THETA_TLENDCAP         =THETA_SIZE_TLENDCAP/N_THETA_TLENDCAP;
  
  int           N_GAMMA_TLBARREL         =std::max(int(2*M_PI*BARREL_INNER_R/(PHI_SEGMENTS*XTAL_LENGTH_T)),1);
  double        D_GAMMA_TLBARREL         =D_PHI_GLOBAL/N_GAMMA_TLBARREL;

  // Main layer Barrel - minimal set needed for timing layer calculations
  double        thC_br_end               =THETA_SIZE_BARREL-D_THETA_BARREL/2;
  double        r0_br_end                =(BARREL_INNER_R+PROJ_OFFSET_R)/cos(thC_br_end);
  double        r0_proj_arm_br_end       =r0_br_end/cos(D_THETA_BARREL/2);
  double        br_phislice_8pa_z1       =r0_proj_arm_br_end*cos(thC_br_end+D_THETA_BARREL/2)-PROJ_OFFSET_R;

  // Main layer Endcap - minimal set needed for timing layer calculations
  double        thC_ec_end               =THETA_SIZE_ENDCAP-D_THETA_ENDCAP/2;
  double        r0_ec_end                =BARREL_HALF_Z/cos(thC_ec_end);
  double        r0_proj_arm_ec_end       =r0_ec_end/cos(D_THETA_ENDCAP/2);
  double        ec_phislice_8pa_z1       =r0_proj_arm_ec_end*cos(thC_ec_end+D_THETA_ENDCAP/2);

  // Timing Barrel Phi envelope with protrusions for tilted crystals
  double        rT                       =br_phislice_8pa_z1-TIMING_GAP-XTAL_DEPTH_T;

  double        thC_tl_end               =THETA_SIZE_TLBARREL-D_THETA_TLBARREL/2;
  double        thC_tl_beg               =D_THETA_TLBARREL/2;

  double        r0_tl_end                =(rT)/cos(thC_tl_end);
  double        r0_tl_beg                =(rT)/cos(thC_tl_beg);

  double        r0_proj_arm_tl_end       =r0_tl_end/cos(D_THETA_TLBARREL/2);
  double        r2_proj_arm_tl_beg       =(r0_tl_beg+XTAL_DEPTH_T)/cos(D_THETA_TLBARREL/2);
  
  double        tl_phislice_8pa_z1       =r0_proj_arm_tl_end*cos(thC_tl_end+D_THETA_TLBARREL/2);
  double        tl_phislice_8pa_z2       =r2_proj_arm_tl_beg*cos(thC_tl_beg-D_THETA_TLBARREL/2)+XTAL_DEPTH_T;

  double        tl_phislice_8pa_y0       =tl_phislice_8pa_z1*tan(thC_tl_end+D_THETA_TLBARREL/2);
  double        tl_phislice_8pa_y2       =tl_phislice_8pa_z2*tan(thC_tl_end+D_THETA_TLBARREL/2);

  // Timing Endcap Phi envelope with protrusions for tilted crystals
  double        zT                       =ec_phislice_8pa_z1-TIMING_GAP-XTAL_DEPTH_T;

  int           TLENDCAP_THETA_START     =0;

  for (int iTheta=0; iTheta<N_THETA_TLENDCAP; iTheta++) {
    double thC=D_THETA_TLENDCAP/2+iTheta*D_THETA_TLENDCAP;
    double RinEndcap=(zT)*tan(thC);
    if (RinEndcap<(BEAMPIPE_OPENING)) TLENDCAP_THETA_START++;
  }

  double        thC_tlec_end             =THETA_SIZE_TLENDCAP-D_THETA_TLENDCAP/2;
  double        thC_tlec_beg             =D_THETA_TLENDCAP/2+TLENDCAP_THETA_START*D_THETA_TLENDCAP;

  double        r0_tlec_end              =(zT)/cos(thC_tlec_end);
  double        r0_tlec_beg              =(zT)/cos(thC_tlec_beg);

  double        r0_proj_arm_tlec_end     =r0_tlec_end/cos(D_THETA_TLENDCAP/2);
  double        r2_proj_arm_tlec_beg     =(r0_tlec_beg+XTAL_DEPTH_T)/cos(D_THETA_TLENDCAP/2);

  double        tlec_phislice_8pa_z1     =r0_proj_arm_tlec_end*cos(thC_tlec_end+D_THETA_TLENDCAP/2);
  double        tlec_phislice_8pa_z2     =r2_proj_arm_tlec_beg*cos(thC_tlec_beg-D_THETA_TLENDCAP/2)+XTAL_DEPTH_T;

  double        tlec_end_phislice_8pa_y0 =tlec_phislice_8pa_z1*tan(thC_tlec_end+D_THETA_TLENDCAP/2);
  double        tlec_end_phislice_8pa_y2 =tlec_phislice_8pa_z2*tan(thC_tlec_end+D_THETA_TLENDCAP/2);

  double        tlec_beg_phislice_8pa_y0 =tlec_phislice_8pa_z1*tan(thC_tlec_beg-D_THETA_TLENDCAP/2);
  double        tlec_beg_phislice_8pa_y2 =tlec_phislice_8pa_z2*tan(thC_tlec_beg-D_THETA_TLENDCAP/2);

  int numCrystalsTiming = 0;
  int numCrystalsTlEndcap = 0;

  std::cout                                                         << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "=SCEPCAL TIMING LAYER INPUTS====================="  << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "BARREL_HALF_Z:        " << BARREL_HALF_Z            << std::endl;
  std::cout << "BARREL_INNER_R:       " << BARREL_INNER_R           << std::endl;
  std::cout << "PHI_SEGMENTS:         " << PHI_SEGMENTS             << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "XTAL_DEPTH_T:         " << XTAL_DEPTH_T             << std::endl;
  std::cout << "XTAL_LENGTH_T:        " << XTAL_LENGTH_T            << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "PROJECTIVE_OFFSET_R:  " << PROJ_OFFSET_R            << std::endl;
  std::cout << "BEAMPIPE_OPENING:     " << BEAMPIPE_OPENING         << std::endl;
  std::cout << "TIMING_LAYER_GAP:     " << TIMING_GAP               << std::endl;
  std::cout                                                         << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "=CONTROL========================================="  << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "CONSTRUCT_TLBARREL:   " << CONSTRUCT_TLBARREL       << std::endl;
  std::cout << "CONSTRUCT_TLENDCAP:   " << CONSTRUCT_TLENDCAP       << std::endl;
  std::cout << "USE_OPTICAL_SURFACES: " << USE_OPTICAL_SURFACES     << std::endl;
  std::cout                                                         << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "=CALCULATED PARAMETERS==========================="  << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "N_THETA_TLBARREL:     " << N_THETA_TLBARREL         << std::endl;
  std::cout << "N_GAMMA_TIMING:       " << N_GAMMA_TLBARREL         << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "N_THETA_TLENDCAP:     " << N_THETA_TLENDCAP         << std::endl;
  std::cout << "TLENDCAP_THETA_START: " << TLENDCAP_THETA_START     << std::endl;
  std::cout                                                         << std::endl;
  std::cout                                                         << std::endl;

  // Initialization
  dd4hep::DetElement ScepcalDetElement(detName,detectorXML.id());
  dd4hep::Volume experimentalHall=theDetector.pickMotherVolume(ScepcalDetElement);
  dd4hep::xml::Dimension sdType=detectorXML.child(_Unicode(sensitive));
  sens.setType(sdType.typeStr());
  dd4hep::Readout readout=sens.readout();
  dd4hep::Segmentation geomseg=readout.segmentation();
  dd4hep::Segmentation* _geoSeg=&geomseg;
  auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCal_TimingSegmentation_k4geo *>(_geoSeg->segmentation());

  dd4hep::OpticalSurfaceManager surfMgr = theDetector.surfaceManager();
  dd4hep::OpticalSurface LYSO_to_ESR = surfMgr.opticalSurface("/world/"+detName+"#LYSO_to_ESR");

  // Global assembly volumes
  std::vector<double> zTimingPolyhedra   ={-tl_phislice_8pa_y2,-tl_phislice_8pa_y0,tl_phislice_8pa_y0,tl_phislice_8pa_y2};
  std::vector<double> rminTimingPolyhedra={tl_phislice_8pa_z2,tl_phislice_8pa_z1,tl_phislice_8pa_z1,tl_phislice_8pa_z2};
  std::vector<double> rmaxTimingPolyhedra={tl_phislice_8pa_z2,tl_phislice_8pa_z2,tl_phislice_8pa_z2,tl_phislice_8pa_z2};
  
  std::vector<double> zTlEndcapPolyhedra   ={tlec_phislice_8pa_z1,tlec_phislice_8pa_z2};
  std::vector<double> rminTlEndcapPolyhedra={tlec_beg_phislice_8pa_y0,tlec_beg_phislice_8pa_y2};
  std::vector<double> rmaxTlEndcapPolyhedra={tlec_end_phislice_8pa_y0,tlec_end_phislice_8pa_y2};

  std::vector<double> zTlEndcapPolyhedra_1   ={-tlec_phislice_8pa_z2,-tlec_phislice_8pa_z1};
  std::vector<double> rminTlEndcapPolyhedra_1={tlec_beg_phislice_8pa_y2,tlec_beg_phislice_8pa_y0};
  std::vector<double> rmaxTlEndcapPolyhedra_1={tlec_end_phislice_8pa_y2,tlec_end_phislice_8pa_y0};

  dd4hep::Polyhedra timingGlobalAssemblyShape(PHI_SEGMENTS,D_PHI_GLOBAL/2,2*M_PI,zTimingPolyhedra,rminTimingPolyhedra,rmaxTimingPolyhedra);
  dd4hep::Volume timingGlobalAssemblyVol("timingGlobalAssemblyVol", timingGlobalAssemblyShape, theDetector.material("Vacuum"));
  timingGlobalAssemblyVol.setVisAttributes(theDetector,tlbarrelAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume tlbarrelAssemblyPlacedVol=experimentalHall.placeVolume(timingGlobalAssemblyVol);
  tlbarrelAssemblyPlacedVol.addPhysVolID("system",TLBARREL_SYSTEM_NO);

  dd4hep::Polyhedra tlendcapGlobalAssemblyShape(PHI_SEGMENTS,D_PHI_GLOBAL/2,2*M_PI,zTlEndcapPolyhedra,rminTlEndcapPolyhedra,rmaxTlEndcapPolyhedra);
  dd4hep::Volume tlendcapGlobalAssemblyVol("tlendcapGlobalAssemblyVol", tlendcapGlobalAssemblyShape, theDetector.material("Vacuum"));
  tlendcapGlobalAssemblyVol.setVisAttributes(theDetector,tlendcapAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume tlendcapAssemblyPlacedVol=experimentalHall.placeVolume(tlendcapGlobalAssemblyVol);
  tlendcapAssemblyPlacedVol.addPhysVolID("system",TLENDCAP_SYSTEM_NO);

  dd4hep::Polyhedra tlendcapGlobalAssemblyShape_1(PHI_SEGMENTS,D_PHI_GLOBAL/2,2*M_PI,zTlEndcapPolyhedra_1,rminTlEndcapPolyhedra_1,rmaxTlEndcapPolyhedra_1);
  dd4hep::Volume tlendcapGlobalAssemblyVol_1("tlendcapGlobalAssemblyVol_1", tlendcapGlobalAssemblyShape_1, theDetector.material("Vacuum"));
  tlendcapGlobalAssemblyVol_1.setVisAttributes(theDetector,tlendcapAssemblyGlobalVisXML.visStr());
  dd4hep::PlacedVolume tlendcapAssemblyPlacedVol_1=experimentalHall.placeVolume(tlendcapGlobalAssemblyVol_1);
  tlendcapAssemblyPlacedVol_1.addPhysVolID("system",TLENDCAP_SYSTEM_NO);
  tlendcapAssemblyPlacedVol_1.addPhysVolID("theta",1);

  ScepcalDetElement.setPlacement(tlbarrelAssemblyPlacedVol);

  // Lambda for crystals
  auto CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId = [
    &theDetector, 
    &sens, 
    &segmentation, 
    &ScepcalDetElement,
    &USE_OPTICAL_SURFACES,
    &LYSO_to_ESR
  ](
    std::string volName,
    double dz,
    const double *vertices,
    xml_comp_t compXml,
    dd4hep::Transform3D transform,
    dd4hep::Volume assemblyVol,
    int nSystem, int nPhi, int nTheta, int nGamma,
    XYZVector posGlobal
  ) {
      dd4hep::EightPointSolid theShape(dz,vertices);
      dd4hep::Volume theVolume(volName, theShape, theDetector.material(compXml.materialStr()));
      theVolume.setVisAttributes(theDetector,compXml.visStr());

      if (USE_OPTICAL_SURFACES) {
        dd4hep::SkinSurface(theDetector, ScepcalDetElement, 
          volName+"Surface_"+std::to_string(nPhi)+"_"+std::to_string(nTheta)+"_"+std::to_string(nGamma), 
          LYSO_to_ESR, theVolume);
      }

      theVolume.setSensitiveDetector(sens);
      auto volID   =segmentation->setVolumeID(nSystem,nPhi,nTheta,nGamma);
      int  volID_32=segmentation->getFirst32bits(volID);
      dd4hep::PlacedVolume thePlacedVol=assemblyVol.placeVolume(theVolume,volID_32,transform);

      thePlacedVol.addPhysVolID("system",nSystem);
      thePlacedVol.addPhysVolID("phi",nPhi);
      thePlacedVol.addPhysVolID("theta",nTheta);
      thePlacedVol.addPhysVolID("gamma",nGamma);

      segmentation->savePosition(volID_32,posGlobal);
  };

  //////////////////////////////
  // Timing Barrel
  // Phi - divide barrel into phi slices
  // Theta - divide phi slice into slabs in theta
  // Gamma - divide a theta slab into shorter segments in phi
  //////////////////////////////
  
  for (int iPhi=CONSTRUCT_TLBARREL? TIMING_PHI_START:TIMING_PHI_END; iPhi<TIMING_PHI_END; iPhi++) {
    double phiGlobal=iPhi*D_PHI_GLOBAL;
    RotationZ rotZphiGlobal(phiGlobal);

    dd4hep::Polyhedra timingPhiAssemblyShape(1,-D_PHI_GLOBAL/2,D_PHI_GLOBAL,
      zTimingPolyhedra,rminTimingPolyhedra,rmaxTimingPolyhedra);
    dd4hep::Volume timingPhiAssemblyVolume("timingPhiAssembly", timingPhiAssemblyShape, theDetector.material("Vacuum"));
    timingPhiAssemblyVolume.setVisAttributes(theDetector,tlbarrelAssemblyPhiVisXML.visStr());
    timingGlobalAssemblyVol.placeVolume(timingPhiAssemblyVolume,Transform3D(rotZphiGlobal));

    for (int iTheta=0; iTheta<N_THETA_TLBARREL; iTheta++) {

      double thC =THETA_SIZE_TLENDCAP+D_THETA_TLBARREL/2+(iTheta*D_THETA_TLBARREL);

      double r0e  =(rT)/sin(thC);
      double r2e  =r0e+XTAL_DEPTH_T;
      double y0e  =r0e*tan(D_THETA_TLBARREL/2.);
      double y2e  =r2e*tan(D_THETA_TLBARREL/2.);

      double x0y0 =r0e*sin(thC)-y0e*cos(thC);
      double x1y0 =r0e*sin(thC)+y0e*cos(thC);

      double x0y2 =r2e*sin(thC)-y2e*cos(thC);
      double x1y2 =r2e*sin(thC)+y2e*cos(thC);
      
      double x0y0l_E=x0y0*tan(-D_PHI_GLOBAL/2);
      double x0y0r_E=x0y0*tan( D_PHI_GLOBAL/2);
      double x1y0l_E=x1y0*tan(-D_PHI_GLOBAL/2);
      double x1y0r_E=x1y0*tan( D_PHI_GLOBAL/2);

      double x0y2l_E=x0y2*tan(-D_PHI_GLOBAL/2);
      double x0y2r_E=x0y2*tan( D_PHI_GLOBAL/2);
      double x1y2l_E=x1y2*tan(-D_PHI_GLOBAL/2);
      double x1y2r_E=x1y2*tan( D_PHI_GLOBAL/2);

      double verticesE[]={x0y0r_E,y0e,x1y0r_E,-y0e,x1y0l_E,-y0e,x0y0l_E,y0e,
                          x0y2r_E,y2e,x1y2r_E,-y2e,x1y2l_E,-y2e,x0y2l_E,y2e};

      double rE=r0e+XTAL_DEPTH_T/2.;
      RotationZYX rotE(M_PI/2,thC,0);
      Position    dispE(rE*sin(thC),0,rE*cos(thC));

      dd4hep::EightPointSolid timingThetaAssemblyShape(XTAL_DEPTH_T/2, verticesE);
      dd4hep::Volume timingThetaAssemblyVolume("timingThetaAssembly", timingThetaAssemblyShape, theDetector.material("Vacuum"));
      timingThetaAssemblyVolume.setVisAttributes(theDetector,tlbarrelAssemblyThetaVisXML.visStr());
      timingPhiAssemblyVolume.placeVolume(timingThetaAssemblyVolume,Transform3D(rotE,dispE));

      for (int nGamma=0; nGamma<N_GAMMA_TLBARREL; nGamma++) {
        double gamma =-D_PHI_GLOBAL/2+D_GAMMA_TLBARREL/2+D_GAMMA_TLBARREL*nGamma;

        double x0y0l=x0y0*tan(gamma-D_GAMMA_TLBARREL/2);
        double x0y0r=x0y0*tan(gamma+D_GAMMA_TLBARREL/2);
        double x1y0l=x1y0*tan(gamma-D_GAMMA_TLBARREL/2);
        double x1y0r=x1y0*tan(gamma+D_GAMMA_TLBARREL/2);
      
        double x0y2l=x0y2*tan(gamma-D_GAMMA_TLBARREL/2);
        double x0y2r=x0y2*tan(gamma+D_GAMMA_TLBARREL/2);
        double x1y2l=x1y2*tan(gamma-D_GAMMA_TLBARREL/2);
        double x1y2r=x1y2*tan(gamma+D_GAMMA_TLBARREL/2);

        double verticesT[]={x0y0r,y0e,x1y0r,-y0e,x1y0l,-y0e,x0y0l,y0e,
                            x0y2r,y2e,x1y2r,-y2e,x1y2l,-y2e,x0y2l,y2e};

        Position dispT(0,0,0);

        XYZVector dispGlobal(rE*sin(thC),rE*sin(thC)*tan(gamma),rE*cos(thC));
        XYZVector posGlobal=(rotZphiGlobal*dispGlobal)/dd4hep::mm;

        CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
          "TimingCrystal", XTAL_DEPTH_T/2, verticesT, crystalTXML,
          Transform3D(dispT),
          timingThetaAssemblyVolume,
          TLBARREL_SYSTEM_NO, iPhi, N_THETA_TLENDCAP+iTheta, nGamma,
          posGlobal
        );
        numCrystalsTiming+=1;
      }
    }
  }

  //////////////////////////////
  // Timing Endcap
  // Phi - divide barrel into phi slices
  // Theta - divide phi slice into slabs in theta
  // Gamma - divide a theta slab into shorter segments in phi
  //////////////////////////////

  for (int iPhi=CONSTRUCT_TLENDCAP? TLENDCAP_PHI_START:TLENDCAP_PHI_END;iPhi<TLENDCAP_PHI_END;iPhi++) {
    double phiGlobal=iPhi*D_PHI_GLOBAL;
    RotationZ rotZphiGlobal(phiGlobal);
    RotationY rotMirror(M_PI);

    dd4hep::Polyhedra tlendcapPhiAssemblyShape(1,-D_PHI_GLOBAL/2,D_PHI_GLOBAL,zTlEndcapPolyhedra,rminTlEndcapPolyhedra,rmaxTlEndcapPolyhedra);
    
    dd4hep::Volume tlendcapPhiAssemblyVolume("tlendcapPhiVol", tlendcapPhiAssemblyShape, theDetector.material("Vacuum"));
    tlendcapPhiAssemblyVolume.setVisAttributes(theDetector,tlendcapAssemblyPhiVisXML.visStr());
    tlendcapGlobalAssemblyVol.placeVolume(tlendcapPhiAssemblyVolume,Transform3D(rotZphiGlobal));

    dd4hep::Volume tlendcapPhiAssemblyVolume_1("tlendcapPhiVol_1", tlendcapPhiAssemblyShape, theDetector.material("Vacuum"));
    tlendcapPhiAssemblyVolume_1.setVisAttributes(theDetector,tlendcapAssemblyPhiVisXML.visStr());
    tlendcapGlobalAssemblyVol_1.placeVolume(tlendcapPhiAssemblyVolume_1,Transform3D(rotMirror*rotZphiGlobal));

    for (int iTheta=TLENDCAP_THETA_START; iTheta<N_THETA_TLENDCAP; iTheta++) {
      
      double thC=D_THETA_TLENDCAP/2+iTheta*D_THETA_TLENDCAP;
      double RinEndcap=(zT)*tan(thC);
      
      int    nGammaEndcap=std::max(int(2*M_PI*RinEndcap/(PHI_SEGMENTS*XTAL_LENGTH_T)),1);
      double dGammaEndcap=D_PHI_GLOBAL/nGammaEndcap;

      double r0e=RinEndcap/sin(thC);
      double r2e=r0e+XTAL_DEPTH_T;
      double y0e=r0e*tan(D_THETA_TLENDCAP/2.);
      double y2e=r2e*tan(D_THETA_TLENDCAP/2.);

      double x0y0 =r0e*sin(thC)-y0e*cos(thC);
      double x1y0 =r0e*sin(thC)+y0e*cos(thC);

      double x0y2 =r2e*sin(thC)-y2e*cos(thC);
      double x1y2 =r2e*sin(thC)+y2e*cos(thC);
      
      double x0y0l_E=x0y0*tan(-D_PHI_GLOBAL/2);
      double x0y0r_E=x0y0*tan( D_PHI_GLOBAL/2);
      double x1y0l_E=x1y0*tan(-D_PHI_GLOBAL/2);
      double x1y0r_E=x1y0*tan( D_PHI_GLOBAL/2);

      double x0y2l_E=x0y2*tan(-D_PHI_GLOBAL/2);
      double x0y2r_E=x0y2*tan( D_PHI_GLOBAL/2);
      double x1y2l_E=x1y2*tan(-D_PHI_GLOBAL/2);
      double x1y2r_E=x1y2*tan( D_PHI_GLOBAL/2);

      double verticesE[]={x0y0r_E,y0e,x1y0r_E,-y0e,x1y0l_E,-y0e,x0y0l_E,y0e,
                          x0y2r_E,y2e,x1y2r_E,-y2e,x1y2l_E,-y2e,x0y2l_E,y2e};

      double rE=r0e+XTAL_DEPTH_T/2.;
      RotationZYX rotE(M_PI/2,thC,0);
      Position    dispE(rE*sin(thC),0,rE*cos(thC));

      dd4hep::EightPointSolid tlendcapThetaAssemblyShape(XTAL_DEPTH_T/2, verticesE);

      dd4hep::Volume tlendcapThetaAssemblyVolume("tlendcapThetaAssembly", tlendcapThetaAssemblyShape, theDetector.material("Vacuum"));
      tlendcapThetaAssemblyVolume.setVisAttributes(theDetector,tlendcapAssemblyThetaVisXML.visStr());
      tlendcapPhiAssemblyVolume.placeVolume(tlendcapThetaAssemblyVolume,Transform3D(rotE,dispE));      

      dd4hep::Volume tlendcapThetaAssemblyVolume_1("tlendcapThetaAssembly_1", tlendcapThetaAssemblyShape, theDetector.material("Vacuum"));
      tlendcapThetaAssemblyVolume_1.setVisAttributes(theDetector,tlendcapAssemblyThetaVisXML.visStr());
      tlendcapPhiAssemblyVolume_1.placeVolume(tlendcapThetaAssemblyVolume_1,Transform3D(rotE,dispE));   

      for (int nGamma=0;nGamma<nGammaEndcap;nGamma++) {
        double gamma=-D_PHI_GLOBAL/2+dGammaEndcap/2+dGammaEndcap*nGamma;
        
        double x0y0l=x0y0*tan(gamma-dGammaEndcap/2);
        double x0y0r=x0y0*tan(gamma+dGammaEndcap/2);
        double x1y0l=x1y0*tan(gamma-dGammaEndcap/2);
        double x1y0r=x1y0*tan(gamma+dGammaEndcap/2);

        double x0y2l=x0y2*tan(gamma-dGammaEndcap/2);
        double x0y2r=x0y2*tan(gamma+dGammaEndcap/2);
        double x1y2l=x1y2*tan(gamma-dGammaEndcap/2);
        double x1y2r=x1y2*tan(gamma+dGammaEndcap/2);

        double verticesT[]={x0y0r,y0e,x1y0r,-y0e,x1y0l,-y0e,x0y0l,y0e,
                            x0y2r,y2e,x1y2r,-y2e,x1y2l,-y2e,x0y2l,y2e};

        Position dispT(0,0,0);

        XYZVector dispGlobal(rE*sin(thC),rE*sin(thC)*tan(gamma),rE*cos(thC));
        XYZVector posGlobal=(rotZphiGlobal*dispGlobal)/dd4hep::mm;

        CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
          "TlEndcapCrystal", XTAL_DEPTH_T/2, verticesT, crystalTXML,
          Transform3D(dispT),
          tlendcapThetaAssemblyVolume,
          TLENDCAP_SYSTEM_NO, iPhi, iTheta, nGamma,
          posGlobal
        );
        numCrystalsTlEndcap+=1;

        CreateEightPointShapeVolume_SetVolAttributes_Place_SetCellId(
          "TlEndcapCrystal_1", XTAL_DEPTH_T/2, verticesT, crystalTXML,
          Transform3D(dispT),
          tlendcapThetaAssemblyVolume_1,
          TLENDCAP_SYSTEM_NO, iPhi, 2*N_THETA_TLENDCAP+N_THETA_TLBARREL-iTheta, nGamma,
          rotMirror*posGlobal
        );
        numCrystalsTlEndcap+=1;
      }
    }
  }

  std::cout                                                         << std::endl;
  std::cout << "NUM_CRYSTALS_TLBARREL:" << numCrystalsTiming        << std::endl; 
  std::cout << "NUM_CRYSTALS_TLENDCAP:" << numCrystalsTlEndcap      << std::endl;
  std::cout                                                         << std::endl;

  return ScepcalDetElement;
}

DECLARE_DETELEMENT(SCEPCal_TimingLayer,create_detector_SCEPCal_TimingLayer)
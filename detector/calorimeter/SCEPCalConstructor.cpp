//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================

#include "detectorSegmentations/SCEPCalSegmentation_k4geo.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetectorTools.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "TGeoTrd2.h"
#include <bitset>

using dd4hep::Transform3D;
using dd4hep::RotationZYX;
using dd4hep::RotationY;
using ROOT::Math::RotationZ;
using dd4hep::Rotation3D;
using dd4hep::Position;

static dd4hep::Ref_t
create_detector_SCEPCal(dd4hep::Detector &theDetector,xml_h xmlElement,dd4hep::SensitiveDetector sens) {
  xml_det_t  detectorXML                 =xmlElement;
  xml_comp_t dimXML                      =detectorXML.child(_Unicode(dim));
  xml_comp_t timingXML                   =detectorXML.child(_Unicode(timing));
  xml_comp_t barrelXML                   =detectorXML.child(_Unicode(barrel));
  xml_comp_t endcapXML                   =detectorXML.child(_Unicode(endcap));
  xml_comp_t projFXML                    =detectorXML.child(_Unicode(projF));
  xml_comp_t projRXML                    =detectorXML.child(_Unicode(projR));
  xml_comp_t crystalFXML                 =detectorXML.child(_Unicode(crystalF));
  xml_comp_t crystalRXML                 =detectorXML.child(_Unicode(crystalR));
  xml_comp_t timingTrXML                 =detectorXML.child(_Unicode(timingLayerTr));
  xml_comp_t timingLgXML                 =detectorXML.child(_Unicode(timingLayerLg));
  xml_comp_t sipmLgXML                   =detectorXML.child(_Unicode(sipmLg));
  xml_comp_t sipmTrXML                   =detectorXML.child(_Unicode(sipmTr));
  // xml_comp_t instXML                     =detectorXML.child(_Unicode(inst));
  xml_comp_t timingAssemblyGlobalVisXML  =detectorXML.child(_Unicode(timingAssemblyGlobalVis));
  xml_comp_t barrelAssemblyGlobalVisXML  =detectorXML.child(_Unicode(barrelAssemblyGlobalVis));
  xml_comp_t endcapAssemblyGlobalVisXML  =detectorXML.child(_Unicode(endcapAssemblyGlobalVis));
  xml_comp_t scepcalAssemblyXML          =detectorXML.child(_Unicode(scepcalAssembly));
  dd4hep::Material crystalFMat           =theDetector.material(crystalFXML.materialStr());
  dd4hep::Material crystalRMat           =theDetector.material(crystalRXML.materialStr());
  dd4hep::Material timingTrMat           =theDetector.material(timingTrXML.materialStr());
  dd4hep::Material timingLgMat           =theDetector.material(timingLgXML.materialStr());
  dd4hep::Material sipmLgMat             =theDetector.material(sipmLgXML.materialStr());
  dd4hep::Material sipmTrMat             =theDetector.material(sipmTrXML.materialStr());
  // dd4hep::Material instMat               =theDetector.material(instXML.materialStr());
  const double  EBz                      =dimXML.attr<double>(_Unicode(barrelHalfZ));
  const double  Rin                      =dimXML.attr<double>(_Unicode(barrelInnerR));
  const double  nomfw                    =dimXML.attr<double>(_Unicode(crystalFaceWidthNominal));
  const double  Fdz                      =dimXML.attr<double>(_Unicode(crystalFlength));
  const double  Rdz                      =dimXML.attr<double>(_Unicode(crystalRlength));
  const double  nomth                    =dimXML.attr<double>(_Unicode(crystalTimingThicknessNominal));
  const double  sipmth                   =dimXML.attr<double>(_Unicode(sipmThickness));
  const int     PHI_SEGMENTS             =dimXML.attr<int>(_Unicode(phiSegments));
  const int     N_PROJECTIVE_FILL        =dimXML.attr<int>(_Unicode(projectiveFill));
  const bool    CONSTRUCT_TIMING         =timingXML.attr<bool>(_Unicode(construct));
  // const int     TIMING_PHI_START         =timingXML.attr<int>(_Unicode(phistart));
  // const int     TIMING_PHI_END           =timingXML.attr<int>(_Unicode(phiend));
  const bool    CONSTRUCT_BARREL         =barrelXML.attr<bool>(_Unicode(construct));
  const int     BARREL_PHI_START         =barrelXML.attr<int>(_Unicode(phistart));
  const int     BARREL_PHI_END           =barrelXML.attr<int>(_Unicode(phiend));
  const bool    CONSTRUCT_ENDCAP         =endcapXML.attr<bool>(_Unicode(construct));
  const int     ENDCAP_PHI_START         =endcapXML.attr<int>(_Unicode(phistart));
  const int     ENDCAP_PHI_END           =endcapXML.attr<int>(_Unicode(phiend));
  const int     ENDCAP_THETA_START       =endcapXML.attr<int>(_Unicode(thetastart));
  const double  D_PHI_GLOBAL             =2*M_PI/PHI_SEGMENTS;
  const double  PROJECTIVE_GAP           =(N_PROJECTIVE_FILL*nomfw)/2;
  double        THETA_SIZE_BARREL        =atan(EBz/Rin);
  double        THETA_SIZE_ENDCAP        =atan(Rin/EBz);
  int           N_THETA_BARREL           =2*floor(EBz/nomfw);
  int           N_THETA_ENDCAP           =floor(Rin/nomfw);
  double        D_THETA_BARREL           =(M_PI-2*THETA_SIZE_ENDCAP)/(N_THETA_BARREL);
  double        D_THETA_ENDCAP           =THETA_SIZE_ENDCAP/N_THETA_ENDCAP;
  int           N_PHI_BARREL_CRYSTAL     =floor(2*M_PI*Rin/(PHI_SEGMENTS*nomfw));
  double        D_PHI_BARREL_CRYSTAL     =D_PHI_GLOBAL/N_PHI_BARREL_CRYSTAL;
  double        thC_end                  =THETA_SIZE_ENDCAP+D_THETA_BARREL/2;
  double        r0slice_end              =Rin/sin(thC_end);
  double        z0slice_end              =r0slice_end*cos(thC_end)+PROJECTIVE_GAP;
  double        y0slice_end              =r0slice_end*tan(D_THETA_BARREL/2.);
  double        slice_front_jut          =y0slice_end*sin(M_PI/2-thC_end);
  double        slice_side_jut           =y0slice_end*cos(M_PI/2-thC_end);
  double        z1slice                  =Rin-slice_front_jut;
  double        z2slice                  =Rin+Fdz+Rdz+slice_front_jut;
  double        zheight_slice            =(z2slice-z1slice)/2;
  double        y1slice                  =z1slice*tan(M_PI/2-THETA_SIZE_ENDCAP)+PROJECTIVE_GAP;
  double        y2slice                  =z2slice*tan(M_PI/2-THETA_SIZE_ENDCAP)+PROJECTIVE_GAP;
  double        rT                       =z1slice-2*nomth;
  double        wT                       =rT *tan(D_PHI_GLOBAL/2);
  int           nTiles                   =ceil(y1slice/wT);
  double        lT                       =2*y1slice/nTiles;
  int           nCy                      =floor(lT/nomth);
  double        actY                     =lT/nCy;
  double        actX                     =2*wT/nCy;
  // double        r2slice_end              =r0slice_end+Fdz+Rdz;
  // double        z2slice_end              =r2slice_end*cos(thC_end)+PROJECTIVE_GAP;
  // double        Rin2slice_end            =r2slice_end*sin(thC_end);
  // double        y2slice_end              =r2slice_end*tan(D_THETA_BARREL/2.);
  // double        slice_front_jut2         =y2slice_end*sin(M_PI/2-thC_end);
  // double        slice_side_jut2          =y2slice_end*cos(M_PI/2-thC_end);
  double        barrelSlice_z1           =z0slice_end+slice_side_jut;
  double        barrelSlice_z2           =y2slice;
  // double        barrelSlice_rmin2        =Rin2slice_end-slice_front_jut2;
  double        thCEnd                   =THETA_SIZE_ENDCAP-D_THETA_ENDCAP/2;
  double        thCBeg                   =D_THETA_ENDCAP/2+ENDCAP_THETA_START*D_THETA_ENDCAP;
  double        r0eEnd                   =EBz/cos(thCEnd);
  // double        r2eEnd                   =r0eEnd+Fdz+Rdz;
  // double        y2eEnd                   =r2eEnd*tan(D_THETA_ENDCAP/2.);
  double        r0eBeg                   =EBz/cos(thCBeg);
  double        r2eBeg                   =r0eBeg+Fdz+Rdz;
  double        y2eBeg                   =r2eBeg*tan(D_THETA_ENDCAP/2.);
  double        aEnd                     =r0eEnd/cos(D_THETA_ENDCAP/2);
  // double        bEnd                     =sqrt(r2eEnd*r2eEnd+y2eEnd*y2eEnd);
  double        z1End                    =aEnd*cos(thCEnd+D_THETA_ENDCAP/2);
  // double        z2End                    =bEnd*cos(thCEnd-D_THETA_ENDCAP/2);
  // double        aBeg                     =r0eBeg/cos(D_THETA_ENDCAP/2);
  double        bBeg                     =sqrt(r2eBeg*r2eBeg+y2eBeg*y2eBeg);
  // double        z1Beg                    =aBeg*cos(thCBeg+D_THETA_ENDCAP/2);
  double        z2Beg                    =bBeg*cos(thCBeg-D_THETA_ENDCAP/2);
  double        z1rmaxE                  =z1End*tan(thCEnd+D_THETA_ENDCAP/2);
  double        z2rmaxE                  =z2Beg*tan(thCEnd+D_THETA_ENDCAP/2);
  double        z1rminB                  =z1End*tan(thCBeg-D_THETA_ENDCAP/2);
  double        z2rminB                  =z2Beg*tan(thCBeg-D_THETA_ENDCAP/2);

  std::cout                                                         << std::endl;
  std::cout << "=GEOMETRY INPUTS="                                  << std::endl;
  std::cout << "BARREL_HALF_Z:        " << EBz                      << std::endl;
  std::cout << "BARREL_INNER_R:       " << Rin                      << std::endl;
  std::cout << "PHI_SEGMENTS:         " << PHI_SEGMENTS             << std::endl;
  std::cout << "N_PROJECTIVE_FILL:    " << N_PROJECTIVE_FILL        << std::endl;
  std::cout << "F_CRYSTAL_LENGTH:     " << Fdz                      << std::endl;
  std::cout << "R_CRYSTAL_LENGTH:     " << Rdz                      << std::endl;
  std::cout << "XTAL_FACE_NOM:        " << nomfw                    << std::endl;
  std::cout << "TIMING_THICK_NOM:     " << nomth                    << std::endl;
  std::cout << "SIPM_THICKNESS:       " << sipmth                   << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "=CONTROL="                                          << std::endl;
  std::cout << "CONSTRUCT_TIMING:     " << CONSTRUCT_TIMING         << std::endl;
  std::cout << "CONSTRUCT_BARREL:     " << CONSTRUCT_BARREL         << std::endl;
  std::cout << "CONSTRUCT_ENDCAP:     " << CONSTRUCT_ENDCAP         << std::endl;
  std::cout << "BARREL_PHI_START:     " << BARREL_PHI_START         << std::endl;
  std::cout << "BARREL_PHI_END  :     " << BARREL_PHI_END           << std::endl;
  std::cout << "ENDCAP_PHI_START:     " << ENDCAP_PHI_START         << std::endl;
  std::cout << "ENDCAP_PHI_END  :     " << ENDCAP_PHI_END           << std::endl;
  std::cout << "ENDCAP_THETA_START:   " << ENDCAP_THETA_START       << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "=CALCULATED PARAMETERS="                            << std::endl;
  std::cout << "D_PHI_GLOBAL:         " << D_PHI_GLOBAL             << std::endl;
  std::cout << "PROJECTIVE_GAP:       " << PROJECTIVE_GAP           << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "=PROJECTIVE LAYER="                                 << std::endl;
  std::cout << "THETA_SIZE_BARREL:    " << THETA_SIZE_BARREL        << std::endl;
  std::cout << "N_THETA_BARREL:       " << N_THETA_BARREL           << std::endl;
  std::cout << "D_THETA_BARREL:       " << D_THETA_BARREL           << std::endl;
  std::cout << "THETA_SIZE_ENDCAP:    " << THETA_SIZE_ENDCAP        << std::endl;
  std::cout << "N_THETA_ENDCAP:       " << N_THETA_ENDCAP           << std::endl;
  std::cout << "D_THETA_ENDCAP:       " << D_THETA_ENDCAP           << std::endl;
  std::cout << "N_PHI_BARREL_CRYSTAL: " << N_PHI_BARREL_CRYSTAL     << std::endl;
  std::cout << "D_PHI_BARREL_CRYSTAL: " << D_PHI_BARREL_CRYSTAL     << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "=TIMING LAYER="                                     << std::endl;
  std::cout << "N_TILES:              " << nTiles                   << std::endl;
  std::cout << "N_XTALS_TILE:         " << nCy                      << std::endl;
  std::cout << "LENGTH_TILE:          " << lT                       << std::endl;
  std::cout << "WIDTH_TILE:           " << 2*wT                     << std::endl;
  std::cout << "ACTUAL_X:             " << actX                     << std::endl;
  std::cout << "ACTUAL_Y:             " << actY                     << std::endl;
  std::cout                                                         << std::endl;

  dd4hep::DetElement ScepcalDetElement(detectorXML.nameStr(),detectorXML.id());
  dd4hep::Volume experimentalHall=theDetector.pickMotherVolume(ScepcalDetElement);
  dd4hep::xml::Dimension sdType=detectorXML.child(_Unicode(sensitive));
  sens.setType(sdType.typeStr());
  dd4hep::Readout readout=sens.readout();
  dd4hep::Segmentation geomseg=readout.segmentation();
  dd4hep::Segmentation* _geoSeg=&geomseg;
  auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCalSegmentation_k4geo *>(_geoSeg->segmentation());
  segmentation->setGeomParams(Fdz,Rdz,nomfw,nomth,EBz,Rin,sipmth,PHI_SEGMENTS,N_PROJECTIVE_FILL);

  std::vector<double> zTimingPolyhedra   ={-barrelSlice_z1,barrelSlice_z1};
  std::vector<double> rminTimingPolyhedra={rT,rT,};
  std::vector<double> rmaxTimingPolyhedra={z1slice,z1slice};
  std::vector<double> zBarrelPolyhedra   ={-barrelSlice_z2,-barrelSlice_z1,barrelSlice_z1,barrelSlice_z2};
  std::vector<double> rminBarrelPolyhedra={z2slice,z1slice,z1slice,z2slice};
  std::vector<double> rmaxBarrelPolyhedra={z2slice,z2slice,z2slice,z2slice};
  std::vector<double> zEndcapPolyhedra   ={z1End+PROJECTIVE_GAP,z2Beg+PROJECTIVE_GAP};
  std::vector<double> rminEndcapPolyhedra={z1rminB,z2rminB};
  std::vector<double> rmaxEndcapPolyhedra={z1rmaxE,z2rmaxE};
  std::vector<double> zEndcap1Polyhedra   ={-(z2Beg+PROJECTIVE_GAP),-(z1End+PROJECTIVE_GAP)};
  std::vector<double> rminEndcap1Polyhedra={z2rminB,z1rminB};
  std::vector<double> rmaxEndcap1Polyhedra={z2rmaxE,z1rmaxE};

  dd4hep::Polyhedra timingAssemblyShape(PHI_SEGMENTS,D_PHI_GLOBAL/2,2*M_PI,zTimingPolyhedra,rminTimingPolyhedra,rmaxTimingPolyhedra);
  dd4hep::Volume    timingAssemblyVol("timingAssemblyVol",timingAssemblyShape,theDetector.material("Vacuum"));
  timingAssemblyVol.setVisAttributes(theDetector,timingAssemblyGlobalVisXML.visStr());
  dd4hep::Polyhedra barrelAssemblyShape(PHI_SEGMENTS,D_PHI_GLOBAL/2,2*M_PI,zBarrelPolyhedra,rminBarrelPolyhedra,rmaxBarrelPolyhedra);
  dd4hep::Volume    barrelAssemblyVol("barrelAssemblyVol",barrelAssemblyShape,theDetector.material("Vacuum"));
  barrelAssemblyVol.setVisAttributes(theDetector,barrelAssemblyGlobalVisXML.visStr());
  dd4hep::Polyhedra endcapAssemblyShape(PHI_SEGMENTS,D_PHI_GLOBAL/2,2*M_PI,zEndcapPolyhedra,rminEndcapPolyhedra,rmaxEndcapPolyhedra);
  dd4hep::Volume    endcapAssemblyVol("endcapAssemblyVol",endcapAssemblyShape,theDetector.material("Vacuum"));
  endcapAssemblyVol.setVisAttributes(theDetector,endcapAssemblyGlobalVisXML.visStr());
  dd4hep::Polyhedra endcap1AssemblyShape(PHI_SEGMENTS,D_PHI_GLOBAL/2,2*M_PI,zEndcap1Polyhedra,rminEndcap1Polyhedra,rmaxEndcap1Polyhedra);
  dd4hep::Volume    endcap1AssemblyVol("endcap1AssemblyVol",endcap1AssemblyShape,theDetector.material("Vacuum"));
  endcap1AssemblyVol.setVisAttributes(theDetector,endcapAssemblyGlobalVisXML.visStr());
  auto timingAssemblyVolId  =segmentation->setVolumeID(6,0,0,0);
  int  timingAssemblyVolId32=segmentation->getFirst32bits(timingAssemblyVolId);
  auto barrelAssemblyVolId  =segmentation->setVolumeID(4,0,0,0);
  int  barrelAssemblyVolId32=segmentation->getFirst32bits(barrelAssemblyVolId);
  auto endcapAssemblyVolId  =segmentation->setVolumeID(5,0,0,0);
  int  endcapAssemblyVolId32=segmentation->getFirst32bits(endcapAssemblyVolId);
  auto endcap1AssemblyVolId  =segmentation->setVolumeID(5,1,0,0);
  int  endcap1AssemblyVolId32=segmentation->getFirst32bits(endcap1AssemblyVolId);
  experimentalHall.placeVolume(timingAssemblyVol,timingAssemblyVolId32);
  dd4hep::PlacedVolume barrelPlacedVol =experimentalHall.placeVolume(barrelAssemblyVol,barrelAssemblyVolId32);
  experimentalHall.placeVolume(endcapAssemblyVol,endcapAssemblyVolId32);
  experimentalHall.placeVolume(endcap1AssemblyVol,endcap1AssemblyVolId32);

  barrelPlacedVol.addPhysVolID("system", detectorXML.id());

  ScepcalDetElement.setPlacement(barrelPlacedVol);

  int numCrystalsBarrel = 0;
  int numCrystalsEndcap = 0;
  int numCrystalsTiming = 0;

  for (int iPhi=CONSTRUCT_BARREL? BARREL_PHI_START:BARREL_PHI_END;iPhi<BARREL_PHI_END;iPhi++) {
    double phiEnvBarrel=iPhi*D_PHI_GLOBAL;
    
    dd4hep::Box timingPhiAssemblyShape(nomth,wT,y1slice);
    dd4hep::Volume timingPhiAssemblyVolume("timingPhiAssembly",timingPhiAssemblyShape,theDetector.material("Vacuum"));
    timingPhiAssemblyVolume.setVisAttributes(theDetector,scepcalAssemblyXML.visStr());
    RotationZ rotZPhi(phiEnvBarrel);
    double rTimingAssembly=rT+nomth;
    Position dispTimingAssembly(rTimingAssembly*cos(phiEnvBarrel),rTimingAssembly*sin(phiEnvBarrel),0);
    timingAssemblyVol.placeVolume(timingPhiAssemblyVolume,Transform3D(rotZPhi,dispTimingAssembly));

    dd4hep::Box timingCrystalLg(nomth/2,actX/2,lT/2-sipmth);
    dd4hep::Box timingCrystalTr(nomth/2,wT-sipmth,actY/2);
    dd4hep::Volume timingCrystalLgVol("TimingCrystalLg",timingCrystalLg,timingLgMat);
    dd4hep::Volume timingCrystalTrVol("TimingCrystalTr",timingCrystalTr,timingTrMat);
    timingCrystalLgVol.setVisAttributes(theDetector,timingLgXML.visStr());
    timingCrystalTrVol.setVisAttributes(theDetector,timingTrXML.visStr());
    timingCrystalLgVol.setSensitiveDetector(sens);
    timingCrystalTrVol.setSensitiveDetector(sens);
    dd4hep::Box sipmBoxLg(nomth/2,actX/2,sipmth/2);
    dd4hep::Box sipmBoxTr(nomth/2,sipmth/2,actY/2);
    dd4hep::Volume sipmBoxLgVol("sipmBoxLg",sipmBoxLg,sipmLgMat);
    dd4hep::Volume sipmBoxTrVol("sipmBoxTr",sipmBoxTr,sipmTrMat);
    sipmBoxLgVol.setVisAttributes(theDetector,sipmLgXML.visStr());
    sipmBoxTrVol.setVisAttributes(theDetector,sipmTrXML.visStr());
    sipmBoxLgVol.setSensitiveDetector(sens);
    sipmBoxTrVol.setSensitiveDetector(sens);

    for (int nTile=CONSTRUCT_TIMING? 0:nTiles;nTile<nTiles;nTile++) {
      dd4hep::Box tileAssemblyShape(nomth,wT,lT/2);
      dd4hep::Volume tileAssemblyVolume("tileAssembly",tileAssemblyShape,theDetector.material("Vacuum"));
      tileAssemblyVolume.setVisAttributes(theDetector,scepcalAssemblyXML.visStr());
      Position dispTileAssembly(0,0,-y1slice+nTile*lT+lT/2);
      timingPhiAssemblyVolume.placeVolume(tileAssemblyVolume,dispTileAssembly);

      for (int nC=0;nC<nCy;nC++) {
        int phiEnvBarrelSign=iPhi%2==0? 1:-1;
        int sign            =nTile%2==0? 1:-1;
        Position dispLg(sign*phiEnvBarrelSign*(nomth/2),-wT+actX/2+ nC*actX,0);
        Position dispTr(sign*phiEnvBarrelSign*(-nomth/2),0,-lT/2+actY/2+nC*actY);
        Position dispSipmLg(0,0,lT/2-sipmth/2);
        Position dispSipmTr(0,wT-sipmth/2,0);

        auto timingLgId64=segmentation->setVolumeID(6,nTile*nCy+nC ,iPhi,3);
        auto timingTrId64=segmentation->setVolumeID(6,nTile*nCy+nC ,iPhi,6);
        int  timingLgId32=segmentation->getFirst32bits(timingLgId64);
        int  timingTrId32=segmentation->getFirst32bits(timingTrId64);
        dd4hep::PlacedVolume timingLgp=tileAssemblyVolume.placeVolume(timingCrystalLgVol,timingLgId32,dispLg);
        dd4hep::PlacedVolume timingTrp=tileAssemblyVolume.placeVolume(timingCrystalTrVol,timingTrId32,dispTr);
        timingLgp.addPhysVolID("system",6);
        timingLgp.addPhysVolID("eta",nTile*nCy+nC);
        timingLgp.addPhysVolID("phi",iPhi);
        timingLgp.addPhysVolID("depth",3);
        timingTrp.addPhysVolID("system",6);
        timingTrp.addPhysVolID("eta",nTile*nCy+nC);
        timingTrp.addPhysVolID("phi",iPhi);
        timingTrp.addPhysVolID("depth",6);
        auto sipmLgId64_1=segmentation->setVolumeID(6,nTile*nCy+nC ,iPhi,4);
        auto sipmLgId64_2=segmentation->setVolumeID(6,nTile*nCy+nC ,iPhi,5);
        auto sipmTrId64_1=segmentation->setVolumeID(6,nTile*nCy+nC ,iPhi,7);
        auto sipmTrId64_2=segmentation->setVolumeID(6,nTile*nCy+nC ,iPhi,8);
        int  sipmLgId32_1=segmentation->getFirst32bits(sipmLgId64_1);
        int  sipmLgId32_2=segmentation->getFirst32bits(sipmLgId64_2);
        int  sipmTrId32_1=segmentation->getFirst32bits(sipmTrId64_1);
        int  sipmTrId32_2=segmentation->getFirst32bits(sipmTrId64_2);
        dd4hep::PlacedVolume sipmLgp1=tileAssemblyVolume.placeVolume(sipmBoxLgVol,sipmLgId32_1,dispLg+dispSipmLg);
        dd4hep::PlacedVolume sipmLgp2=tileAssemblyVolume.placeVolume(sipmBoxLgVol,sipmLgId32_2,dispLg-dispSipmLg);
        dd4hep::PlacedVolume sipmTrp1=tileAssemblyVolume.placeVolume(sipmBoxTrVol,sipmTrId32_1,dispTr+dispSipmTr);
        dd4hep::PlacedVolume sipmTrp2=tileAssemblyVolume.placeVolume(sipmBoxTrVol,sipmTrId32_2,dispTr-dispSipmTr);
        sipmLgp1.addPhysVolID("system",6);
        sipmLgp1.addPhysVolID("eta",nTile*nCy+nC);
        sipmLgp1.addPhysVolID("phi",iPhi);
        sipmLgp1.addPhysVolID("depth",4);
        sipmLgp2.addPhysVolID("system",6);
        sipmLgp2.addPhysVolID("eta",nTile*nCy+nC);
        sipmLgp2.addPhysVolID("phi",iPhi);
        sipmLgp2.addPhysVolID("depth",5);
        sipmTrp1.addPhysVolID("system",6);
        sipmTrp1.addPhysVolID("eta",nTile*nCy+nC);
        sipmTrp1.addPhysVolID("phi",iPhi);
        sipmTrp1.addPhysVolID("depth",7);
        sipmTrp2.addPhysVolID("system",6);
        sipmTrp2.addPhysVolID("eta",nTile*nCy+nC);
        sipmTrp2.addPhysVolID("phi",iPhi);
        sipmTrp2.addPhysVolID("depth",8);

        numCrystalsTiming+=2;
      }
    }

    for (int nGamma=0; nGamma<N_PHI_BARREL_CRYSTAL; nGamma++) {
      double gamma =-D_PHI_GLOBAL/2+D_PHI_BARREL_CRYSTAL/2+D_PHI_BARREL_CRYSTAL*nGamma;
      double x0y0le=z1slice*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
      double x0y0re=z1slice*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
      double x1y1le=z2slice*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
      double x1y1re=z2slice*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
      double verticesS[]={-y1slice,x0y0le,-y1slice,x0y0re,y1slice,x0y0re,y1slice,x0y0le,
                         -y2slice,x1y1le,-y2slice,x1y1re,y2slice,x1y1re,y2slice,x1y1le};
      double rSlice =(z1slice+z2slice)/2;
      RotationZYX rotSlice(0,M_PI/2,0);RotationZ rotZSlice(phiEnvBarrel);rotSlice=rotZSlice*rotSlice;
      Position dispSlice(rSlice*cos(phiEnvBarrel),rSlice*sin(phiEnvBarrel),0);

      dd4hep::EightPointSolid barrelPhiAssemblyShape(zheight_slice,verticesS);
      dd4hep::Volume barrelPhiAssemblyVolume("barrelPhiAssembly",barrelPhiAssemblyShape,theDetector.material("Vacuum"));
      barrelPhiAssemblyVolume.setVisAttributes(theDetector,scepcalAssemblyXML.visStr());
      barrelAssemblyVol.placeVolume(barrelPhiAssemblyVolume,Transform3D(rotSlice,dispSlice));

      for (int iTheta=0; iTheta<N_THETA_BARREL; iTheta++) {
        double thC  =THETA_SIZE_ENDCAP+D_THETA_BARREL/2+(iTheta*D_THETA_BARREL);
        double r0e  =Rin/sin(thC);
        double r1e  =r0e+Fdz;
        double r2e  =r1e+Rdz;
        double y0e  =r0e*tan(D_THETA_BARREL/2.);
        double y1e  =r1e*tan(D_THETA_BARREL/2.);
        double y2e  =r2e*tan(D_THETA_BARREL/2.);
        double x0y0 =r0e*sin(thC)-y0e*cos(thC);
        double x1y0 =r0e*sin(thC)+y0e*cos(thC);
        double x0y0l=x0y0*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y0r=x0y0*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x1y0l=x1y0*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y0r=x1y0*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x0y1 =r1e*sin(thC)-y1e*cos(thC);
        double x1y1 =r1e*sin(thC)+y1e*cos(thC);
        double x0y1l=x0y1*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y1r=x0y1*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x1y1l=x1y1*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y1r=x1y1*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x0y2 =r2e*sin(thC)-y2e*cos(thC);
        double x1y2 =r2e*sin(thC)+y2e*cos(thC);
        double x0y2l=x0y2*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y2r=x0y2*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x1y2l=x1y2*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y2r=x1y2*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double verticesF[]={x0y0r,y0e,x1y0r,-y0e,x1y0l,-y0e,x0y0l,y0e,
                            x0y1r,y1e,x1y1r,-y1e,x1y1l,-y1e,x0y1l,y1e};
        double verticesR[]={x0y1r,y1e,x1y1r,-y1e,x1y1l,-y1e,x0y1l,y1e,
                            x0y2r,y2e,x1y2r,-y2e,x1y2l,-y2e,x0y2l,y2e};
        double rF=r0e+Fdz/2.;
        double rR=r1e+Rdz/2.;
        int projective_sign=cos(thC)>0?-1:1;
        RotationZYX rot(M_PI/2,-M_PI/2+thC,0);
        Position dispF(-rF*cos(thC)+projective_sign*PROJECTIVE_GAP,0,-(rSlice-rF*sin(thC)));
        Position dispR(-rR*cos(thC)+projective_sign*PROJECTIVE_GAP,0,-(rSlice-rR*sin(thC)));

        dd4hep::EightPointSolid crystalFShape(Fdz/2,verticesF);
        dd4hep::EightPointSolid crystalRShape(Rdz/2,verticesR);
        dd4hep::Volume crystalFVol("BarrelCrystalF",crystalFShape,crystalFMat);
        dd4hep::Volume crystalRVol("BarrelCrystalR",crystalRShape,crystalRMat);
        crystalFVol.setVisAttributes(theDetector,crystalFXML.visStr());
        crystalRVol.setVisAttributes(theDetector,crystalRXML.visStr());
        crystalFVol.setSensitiveDetector(sens);
        crystalRVol.setSensitiveDetector(sens);
        auto crystalFId64=segmentation->setVolumeID(4,N_THETA_ENDCAP+iTheta ,iPhi*N_PHI_BARREL_CRYSTAL+nGamma,1);
        auto crystalRId64=segmentation->setVolumeID(4,N_THETA_ENDCAP+iTheta ,iPhi*N_PHI_BARREL_CRYSTAL+nGamma,2);
        int  crystalFId32=segmentation->getFirst32bits(crystalFId64);
        int  crystalRId32=segmentation->getFirst32bits(crystalRId64);
        dd4hep::PlacedVolume crystalFp=barrelPhiAssemblyVolume.placeVolume(crystalFVol,crystalFId32,Transform3D(rot,dispF));
        dd4hep::PlacedVolume crystalRp=barrelPhiAssemblyVolume.placeVolume(crystalRVol,crystalRId32,Transform3D(rot,dispR));
        crystalFp.addPhysVolID("system",4);
        crystalFp.addPhysVolID("eta",N_THETA_ENDCAP+iTheta);
        crystalFp.addPhysVolID("phi",iPhi*N_PHI_BARREL_CRYSTAL+nGamma);
        crystalFp.addPhysVolID("depth",1);
        crystalRp.addPhysVolID("system",4);
        crystalRp.addPhysVolID("eta",N_THETA_ENDCAP+iTheta);
        crystalRp.addPhysVolID("phi",iPhi*N_PHI_BARREL_CRYSTAL+nGamma);
        crystalRp.addPhysVolID("depth",2);

        numCrystalsBarrel+=2;

      }

      for (int iTheta=0; iTheta<N_PROJECTIVE_FILL; iTheta++) {
        double thC  =M_PI/2;
        double r0e  =Rin/sin(thC);
        double r1e  =r0e+Fdz;
        double r2e  =r1e+Rdz;
        double y0e  =nomfw/2;
        double y1e  =nomfw/2;
        double y2e  =nomfw/2;
        double x0y0 =r0e*sin(thC)-y0e*cos(thC);
        double x1y0 =r0e*sin(thC)+y0e*cos(thC);
        double x0y0l=x0y0*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y0r=x0y0*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x1y0l=x1y0*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y0r=x1y0*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x0y1 =r1e*sin(thC)-y1e*cos(thC);
        double x1y1 =r1e*sin(thC)+y1e*cos(thC);
        double x0y1l=x0y1*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y1r=x0y1*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x1y1l=x1y1*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y1r=x1y1*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x0y2 =r2e*sin(thC)-y2e*cos(thC);
        double x1y2 =r2e*sin(thC)+y2e*cos(thC);
        double x0y2l=x0y2*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x0y2r=x0y2*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double x1y2l=x1y2*tan(gamma-D_PHI_BARREL_CRYSTAL/2);
        double x1y2r=x1y2*tan(gamma+D_PHI_BARREL_CRYSTAL/2);
        double verticesF[]={x0y0r,y0e,x1y0r,-y0e,x1y0l,-y0e,x0y0l,y0e,
                            x0y1r,y1e,x1y1r,-y1e,x1y1l,-y1e,x0y1l,y1e};
        double verticesR[]={x0y1r,y1e,x1y1r,-y1e,x1y1l,-y1e,x0y1l,y1e,
                            x0y2r,y2e,x1y2r,-y2e,x1y2l,-y2e,x0y2l,y2e};
        double rF=r0e+Fdz/2.;
        double rR=r1e+Rdz/2.;
        RotationZYX rot(M_PI/2,-M_PI/2+thC,0);
        Position dispF(-rF*cos(thC)-nomfw*(N_PROJECTIVE_FILL-1)/2+iTheta*nomfw,0,-(rSlice-rF*sin(thC)));
        Position dispR(-rR*cos(thC)-nomfw*(N_PROJECTIVE_FILL-1)/2+iTheta*nomfw,0,-(rSlice-rR*sin(thC)));

        dd4hep::EightPointSolid crystalFShape(Fdz/2,verticesF);
        dd4hep::EightPointSolid crystalRShape(Rdz/2,verticesR);
        dd4hep::Volume crystalFVol("BarrelCrystalF",crystalFShape,crystalFMat);
        dd4hep::Volume crystalRVol("BarrelCrystalR",crystalRShape,crystalRMat);
        crystalFVol.setVisAttributes(theDetector,projFXML.visStr());
        crystalRVol.setVisAttributes(theDetector,projRXML.visStr());
        crystalFVol.setSensitiveDetector(sens);
        crystalRVol.setSensitiveDetector(sens);
        auto crystalFId64=segmentation->setVolumeID(7,iTheta,iPhi*N_PHI_BARREL_CRYSTAL+nGamma,1);
        auto crystalRId64=segmentation->setVolumeID(7,iTheta,iPhi*N_PHI_BARREL_CRYSTAL+nGamma,2);
        int  crystalFId32=segmentation->getFirst32bits(crystalFId64);
        int  crystalRId32=segmentation->getFirst32bits(crystalRId64);
        dd4hep::PlacedVolume crystalFp=barrelPhiAssemblyVolume.placeVolume(crystalFVol,crystalFId32,Transform3D(rot,dispF));
        dd4hep::PlacedVolume crystalRp=barrelPhiAssemblyVolume.placeVolume(crystalRVol,crystalRId32,Transform3D(rot,dispR));
        crystalFp.addPhysVolID("system",7);
        crystalFp.addPhysVolID("eta",iTheta);
        crystalFp.addPhysVolID("phi",iPhi*N_PHI_BARREL_CRYSTAL+nGamma);
        crystalFp.addPhysVolID("depth",1);
        crystalRp.addPhysVolID("system",7);
        crystalRp.addPhysVolID("eta",iTheta);
        crystalRp.addPhysVolID("phi",iPhi*N_PHI_BARREL_CRYSTAL+nGamma);
        crystalRp.addPhysVolID("depth",2);

        numCrystalsBarrel+=2;
      }
    }
  }

  for (int iTheta=CONSTRUCT_ENDCAP? ENDCAP_THETA_START:N_THETA_ENDCAP;iTheta<N_THETA_ENDCAP;iTheta++) {
    double thC=D_THETA_ENDCAP/2+iTheta*D_THETA_ENDCAP;
    double RinEndcap=EBz*tan(thC);
    int    nPhiEndcapCrystal=floor(2*M_PI*RinEndcap/(PHI_SEGMENTS*nomfw));
    double dPhiEndcapCrystal=D_PHI_GLOBAL/nPhiEndcapCrystal;
    double r0e=RinEndcap/sin(thC);
    double r1e=r0e+Fdz;
    double r2e=r1e+Rdz;
    double y0e=r0e*tan(D_THETA_ENDCAP/2.);
    double y1e=r1e*tan(D_THETA_ENDCAP/2.);
    double y2e=r2e*tan(D_THETA_ENDCAP/2.);
    double a=r0e/cos(D_THETA_ENDCAP/2);
    double z1=a*cos(thC+D_THETA_ENDCAP/2);
    double r1min=z1*tan(thC-D_THETA_ENDCAP/2);
    double r1max=z1*tan(thC+D_THETA_ENDCAP/2);
    double b=sqrt(r2e*r2e+y2e*y2e);
    double z2=b*cos(thC-D_THETA_ENDCAP/2);
    double r2min=z2*tan(thC-D_THETA_ENDCAP/2);
    double r2max=z2*tan(thC+D_THETA_ENDCAP/2);

    std::vector<double> zPolyhedra={z1+PROJECTIVE_GAP,z2+PROJECTIVE_GAP};
    std::vector<double> rminPolyhedra={r1min,r2min};
    std::vector<double> rmaxPolyhedra={r1max,r2max};
    dd4hep::Polyhedra endcapRingAssemblyShape(PHI_SEGMENTS,D_PHI_GLOBAL/2,2*M_PI,zPolyhedra,rminPolyhedra,rmaxPolyhedra);
    dd4hep::Volume    endcapRingAssemblyVolume("endcapRingAssembly",endcapRingAssemblyShape,theDetector.material("Vacuum"));
    endcapRingAssemblyVolume.setVisAttributes(theDetector,scepcalAssemblyXML.visStr());
    dd4hep::Volume    endcap1RingAssemblyVolume("endcapRingAssembly",endcapRingAssemblyShape,theDetector.material("Vacuum"));
    endcap1RingAssemblyVolume.setVisAttributes(theDetector,scepcalAssemblyXML.visStr());
    RotationY rotMirror(M_PI);
    endcapAssemblyVol.placeVolume(endcapRingAssemblyVolume);
    endcap1AssemblyVol.placeVolume(endcap1RingAssemblyVolume,Transform3D(rotMirror));

    for (int iPhi=ENDCAP_PHI_START;iPhi<ENDCAP_PHI_END;iPhi++) {
      double phiEnvEndcap=iPhi*D_PHI_GLOBAL;
      for (int nGamma=0;nGamma<nPhiEndcapCrystal;nGamma++) {
        double gamma=-D_PHI_GLOBAL/2+dPhiEndcapCrystal/2+dPhiEndcapCrystal*nGamma;
        double x0y0 =r0e*sin(thC)-y0e*cos(thC);
        double x1y0 =r0e*sin(thC)+y0e*cos(thC);
        double x0y0l=x0y0*tan(gamma-dPhiEndcapCrystal/2);
        double x0y0r=x0y0*tan(gamma+dPhiEndcapCrystal/2);
        double x1y0l=x1y0*tan(gamma-dPhiEndcapCrystal/2);
        double x1y0r=x1y0*tan(gamma+dPhiEndcapCrystal/2);
        double x0y1 =r1e*sin(thC)-y1e*cos(thC);
        double x1y1 =r1e*sin(thC)+y1e*cos(thC);
        double x0y1l=x0y1*tan(gamma-dPhiEndcapCrystal/2);
        double x0y1r=x0y1*tan(gamma+dPhiEndcapCrystal/2);
        double x1y1l=x1y1*tan(gamma-dPhiEndcapCrystal/2);
        double x1y1r=x1y1*tan(gamma+dPhiEndcapCrystal/2);
        double x0y2 =r2e*sin(thC)-y2e*cos(thC);
        double x1y2 =r2e*sin(thC)+y2e*cos(thC);
        double x0y2l=x0y2*tan(gamma-dPhiEndcapCrystal/2);
        double x0y2r=x0y2*tan(gamma+dPhiEndcapCrystal/2);
        double x1y2l=x1y2*tan(gamma-dPhiEndcapCrystal/2);
        double x1y2r=x1y2*tan(gamma+dPhiEndcapCrystal/2);
        double verticesF[]={x0y0r,y0e,x1y0r,-y0e,x1y0l,-y0e,x0y0l,y0e,
                            x0y1r,y1e,x1y1r,-y1e,x1y1l,-y1e,x0y1l,y1e};
        double verticesR[]={x0y1r,y1e,x1y1r,-y1e,x1y1l,-y1e,x0y1l,y1e,
                            x0y2r,y2e,x1y2r,-y2e,x1y2l,-y2e,x0y2l,y2e};
        double rF=r0e+Fdz/2.;
        double rR=r1e+Rdz/2.;
        RotationZYX rot(M_PI/2,thC,0);RotationZ rotZ(phiEnvEndcap);rot=rotZ*rot;
        Position dispF(rF*sin(thC)*cos(phiEnvEndcap),rF*sin(thC)*sin(phiEnvEndcap),rF*cos(thC)+PROJECTIVE_GAP);
        Position dispR(rR*sin(thC)*cos(phiEnvEndcap),rR*sin(thC)*sin(phiEnvEndcap),rR*cos(thC)+PROJECTIVE_GAP);

        dd4hep::EightPointSolid crystalFShape(Fdz/2,verticesF);
        dd4hep::EightPointSolid crystalRShape(Rdz/2,verticesR);
        dd4hep::Volume crystalFVol("EndcapCrystalF",crystalFShape,crystalFMat);
        dd4hep::Volume crystalRVol("EndcapCrystalR",crystalRShape,crystalRMat);
        crystalFVol.setVisAttributes(theDetector,crystalFXML.visStr());
        crystalRVol.setVisAttributes(theDetector,crystalRXML.visStr());
        crystalFVol.setSensitiveDetector(sens);
        crystalRVol.setSensitiveDetector(sens);
        auto crystalFId64=segmentation->setVolumeID(5,iTheta,iPhi*nPhiEndcapCrystal+nGamma,1);
        auto crystalRId64=segmentation->setVolumeID(5,iTheta,iPhi*nPhiEndcapCrystal+nGamma,2);
        int  crystalFId32=segmentation->getFirst32bits(crystalFId64);
        int  crystalRId32=segmentation->getFirst32bits(crystalRId64);
        dd4hep::PlacedVolume crystalFp=endcapRingAssemblyVolume.placeVolume(crystalFVol,crystalFId32,Transform3D(rot,dispF));
        dd4hep::PlacedVolume crystalRp=endcapRingAssemblyVolume.placeVolume(crystalRVol,crystalRId32,Transform3D(rot,dispR));
        crystalFp.addPhysVolID("system",5);
        crystalFp.addPhysVolID("eta",iTheta);
        crystalFp.addPhysVolID("phi",iPhi*nPhiEndcapCrystal+nGamma);
        crystalFp.addPhysVolID("depth",1);
        crystalRp.addPhysVolID("system",5);
        crystalRp.addPhysVolID("eta",iTheta);
        crystalRp.addPhysVolID("phi",iPhi*nPhiEndcapCrystal+nGamma);
        crystalRp.addPhysVolID("depth",2);
        auto crystalFId641=segmentation->setVolumeID(5,N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP-iTheta,iPhi*nPhiEndcapCrystal+nGamma,1);
        auto crystalRId641=segmentation->setVolumeID(5,N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP-iTheta,iPhi*nPhiEndcapCrystal+nGamma,2);
        int  crystalFId321=segmentation->getFirst32bits(crystalFId641);
        int  crystalRId321=segmentation->getFirst32bits(crystalRId641);
        dd4hep::PlacedVolume crystalFp1=endcap1RingAssemblyVolume.placeVolume(crystalFVol,crystalFId321,Transform3D(rot,dispF));
        dd4hep::PlacedVolume crystalRp1=endcap1RingAssemblyVolume.placeVolume(crystalRVol,crystalRId321,Transform3D(rot,dispR));
        crystalFp1.addPhysVolID("system",5);
        crystalFp1.addPhysVolID("eta",N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP-iTheta);
        crystalFp1.addPhysVolID("phi",iPhi*nPhiEndcapCrystal+nGamma);
        crystalFp1.addPhysVolID("depth",1);
        crystalRp1.addPhysVolID("system",5);
        crystalRp1.addPhysVolID("eta",N_THETA_ENDCAP+N_THETA_BARREL+N_THETA_ENDCAP-iTheta);
        crystalRp1.addPhysVolID("phi",iPhi*nPhiEndcapCrystal+nGamma);
        crystalRp1.addPhysVolID("depth",2);

        numCrystalsEndcap+=2;
      }
    }
  }

  std::cout                                                         << std::endl;
  std::cout                                                         << std::endl;
  std::cout << "NUM_CRYSTALS_BARREL:  " << numCrystalsBarrel        << std::endl;
  std::cout << "NUM_CRYSTALS_ENDCAP:  " << numCrystalsEndcap        << std::endl;
  std::cout << "NUM_CRYSTALS_TIMING:  " << numCrystalsTiming        << std::endl;
  std::cout                                                         << std::endl;
  std::cout                                                         << std::endl;

  return ScepcalDetElement;
}

DECLARE_DETELEMENT(SegmentedCrystalECAL,create_detector_SCEPCal)

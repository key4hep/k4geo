#include "DD4hep/DetFactoryHelper.h"
#include "TMatrixT.h"

// todo: remove gaudi logging and properly capture output
#define endmsg std::endl
#define lLog std::cout
namespace MSG {
const std::string ERROR = " Error: ";
const std::string DEBUG = " Debug: ";
const std::string INFO  = " Info: ";
}

namespace det {

  double tForArcLength(double s, double plateangle, double delZ, double r) {

    // some intermediate constants
    double zpos = delZ/2.;
    double zp = zpos/TMath::Tan(plateangle);
    double b = zp/(TMath::Sqrt(r*r-zp*zp));
    double c = (TMath::Tan(s/r) +b)/(1.-b*TMath::Tan(s/r));
    // double A = (1-c*c)/(TMath::Sin(plateangle)*TMath::Sin(plateangle));
    //double B = 2*zp*(1-c*c)/TMath::Sin(plateangle);
    //double C = zp*zp-c*c*r*r;
    double D = c*c*r*r/(1+c*c);
    return (TMath::Sqrt(D)-zp)*TMath::Sin(plateangle);

  }

  dd4hep::Solid buildOneBlade(double thickness_inner,
			      double thickness_outer,
			      double width,
			      double ro, double ri,
			      double plateangle,
			      double delZ) 
  {

    dd4hep::Trd2 tmp1(thickness_inner/2., thickness_outer/2., width/2., width/2., ro/2. );
    dd4hep::Tube outerTube(ro, ro*100, delZ);
    dd4hep::Tube innerTube(0., ri, delZ); 
    dd4hep::SubtractionSolid tmp2(tmp1, outerTube, dd4hep::Transform3D(dd4hep::RotationZYX( 0, TMath::Pi()/2.-plateangle, TMath::Pi()/2.),dd4hep::Position(0, 0, ro/2.)));
    //    dd4hep::SubtractionSolid electrodePlate(electrodePlateOuterTrimmed, innerTube, dd4hep::RotationZYX(PlateAngle, 0., 0.));
    return dd4hep::SubtractionSolid(tmp2, innerTube, dd4hep::Transform3D(dd4hep::RotationZYX( 0 , TMath::Pi()/2.-plateangle, TMath::Pi()/2.),dd4hep::Position(0, 0, ro/2.)));

  }
			      
  void buildSubCylinder(dd4hep::Detector& aLcdd,
			dd4hep::SensitiveDetector& aSensDet,
			dd4hep::Volume& aEnvelope,
			dd4hep::xml::Handle_t& aXmlElement,
			float ri, float ro,
			unsigned &iModule,
			int sign) {

    dd4hep::xml::DetElement genericBladeElem = aXmlElement.child(_Unicode(turbineBlade));
    dd4hep::xml::DetElement absBladeElem = genericBladeElem.child(_Unicode(absorberBlade));
    dd4hep::xml::DetElement electrodeBladeElem = genericBladeElem.child(_Unicode(electrodeBlade));
    dd4hep::xml::DetElement nobleLiquidElem = genericBladeElem.child(_Unicode(nobleLiquidGap));
    

    float PlateAngle = genericBladeElem.attr<float>(_Unicode(angle));
    bool decreaseAnglePerSubCylinder = genericBladeElem.attr<bool>(_Unicode(decreaseAnglePerSubCylinder));
    std::cout << "Plate angle is " << PlateAngle << "; decrease angle per subcylinder? " << decreaseAnglePerSubCylinder << std::endl;
    dd4hep::xml::Dimension dim(aXmlElement.child(_Unicode(dimensions)));
    double grmin = dim.rmin1();
    double delZ = dim.dz()*2;
    std::cout << "delZ is " << delZ << std::endl;
    if (decreaseAnglePerSubCylinder) {
      float tubeFracCovered = delZ/(2*grmin*TMath::Tan(PlateAngle));
      PlateAngle = TMath::ATan(delZ/(2*ri*tubeFracCovered));
    }
    
    if (TMath::Abs(TMath::Tan(PlateAngle)) < delZ/(2.*ri)) {
      std::cout << "The requested plate angle is too small for the given delZ and ri values.  Please adjust to at least " << TMath::ATan(delZ/(2.*ri))*180./TMath::Pi() << " degrees!" << std::endl;
      return;
    }

    Float_t xRange = delZ/(TMath::Sin(PlateAngle));

    int nPlates;
    double delrPhi, delrPhiNoGap;
    
    float AbsThicki = absBladeElem.attr<float>(_Unicode(thickness));
    float ElectrodeThick = electrodeBladeElem.attr<float>(_Unicode(thickness));
    float LArgapi = nobleLiquidElem.attr<float>(_Unicode(gap));
    
    bool sameNplates = genericBladeElem.attr<bool>(_Unicode(sameNplates));;

    bool scalePlateThickness = absBladeElem.attr<bool>(_Unicode(scaleThickness));
    float AbsThicko;
    if (scalePlateThickness) {
      AbsThicko = AbsThicki*(ro/ri);
    } else {
      AbsThicko = AbsThicki;
    }

    if (sameNplates) {
      nPlates = 2*TMath::Pi()*grmin/(2*LArgapi+AbsThicki+ElectrodeThick)*TMath::Sin(PlateAngle);
    } else {
      double x1 = delZ/2./TMath::Tan(PlateAngle);
      double x2 = delZ/2./TMath::Tan(PlateAngle)+(2*LArgapi+AbsThicki+ElectrodeThick)/TMath::Sin(PlateAngle);
      double y1 = TMath::Sqrt(ri*ri-x1*x1);
      double y2 = TMath::Sqrt(ri*ri-x2*x2);
      double rPhi1 = ri*TMath::ATan(x1/y1);
      double rPhi2 = ri*TMath::ATan(x2/y2);
      delrPhi = TMath::Abs(rPhi1-rPhi2);
      nPlates = 2*TMath::Pi()*ri/delrPhi;

    }

    // adjust gap thickness at inner layer
    double circ = 2*TMath::Pi()*ri;
    double x1 = delZ/2./TMath::Tan(PlateAngle);
    double x2 = delZ/2./TMath::Tan(PlateAngle)+(AbsThicki+ElectrodeThick)/TMath::Sin(PlateAngle);
    double y1 = TMath::Sqrt(ri*ri-x1*x1);
    double y2 = TMath::Sqrt(ri*ri-x2*x2);
    double rPhi1 = ri*TMath::ATan(x1/y1);
    double rPhi2 = ri*TMath::ATan(x2/y2);    
    delrPhiNoGap = TMath::Abs(rPhi1-rPhi2);
    double leftoverS = (circ - nPlates*delrPhiNoGap);
    double delrPhiGapOnly = leftoverS/(2*nPlates);
    std::cout << "LAr gap was " << LArgapi ;
    LArgapi = tForArcLength(delrPhiGapOnly, PlateAngle, delZ, ri);
    std::cout << " but is now " << LArgapi << " since s = " << delrPhiGapOnly << " s for full unit cell is " << delrPhi << std::endl;

    // now find gap at outer layer
    circ = 2*TMath::Pi()*ro;
    x1 = delZ/2./TMath::Tan(PlateAngle);
    x2 = delZ/2./TMath::Tan(PlateAngle)+(AbsThicko+ElectrodeThick)/TMath::Sin(PlateAngle);
    y1 = TMath::Sqrt(ro*ro-x1*x1);
    y2 = TMath::Sqrt(ro*ro-x2*x2);
    rPhi1 = ro*TMath::ATan(x1/y1);
    rPhi2 = ro*TMath::ATan(x2/y2);    
    delrPhiNoGap = TMath::Abs(rPhi1-rPhi2);
    leftoverS = (circ - nPlates*delrPhiNoGap);
    delrPhiGapOnly = leftoverS/(2*nPlates);   
    double LArgapo = tForArcLength(delrPhiGapOnly, PlateAngle, delZ, ro);
    std::cout << "Outer LAr gap is " << LArgapo << std::endl ;
    
    dd4hep::Solid  absPlate = buildOneBlade(AbsThicki, AbsThicko, xRange, ro, ri, PlateAngle, delZ );
    dd4hep::Volume absPlateVol("absPlate", absPlate, aLcdd.material("Lead"));
    dd4hep::Solid electrodePlateAndGap = buildOneBlade(ElectrodeThick+LArgapo*2, ElectrodeThick+LArgapi*2, xRange, ro, ri, PlateAngle, delZ);
    dd4hep::Volume electrodePlateAndGapVol("electrodePlateAndGap", electrodePlateAndGap, aLcdd.material("Air"));
     dd4hep::Solid electrodePlate = buildOneBlade(ElectrodeThick, ElectrodeThick, xRange, ro, ri, PlateAngle, delZ);
    dd4hep::Volume electrodePlateVol("electrodePlate", electrodePlate, aLcdd.material("PCB"));
    dd4hep::SubtractionSolid LArShape(electrodePlateAndGap, electrodePlate);
    dd4hep::Volume LArVol("nobleLiquid", LArShape, aLcdd.material("LAr"));
    LArVol.setSensitiveDetector(aSensDet);

    //    electrodePlateVol.setSensitiveDetector(aSensDet);
    
    int    nPlatesToDraw = nPlates;
    // nPlatesToDraw = 1;

    for (int iPlate = 0; iPlate < nPlatesToDraw; iPlate++) {
      //    pl->SetInvisible();
      float phi = iPlate*2*TMath::Pi()/nPlates;
      float delPhi = 2*TMath::Pi()/nPlates;
      
      float x = (-ro/2.)*TMath::Cos(phi);
      float y = (-ro/2.)*TMath::Sin(phi); //ri*TMath::Sin(phi)/6.;
      float z =  0.;
      
      //      x = (ri+(ro-ri)/2.)*TMath::Cos(phi+delPhi);
      // y = (ri+(ro-ri)/2.)*TMath::Sin(phi+delPhi);


      TGeoRotation tgr;
      tgr.RotateZ(PlateAngle*180/TMath::Pi());
      tgr.RotateX(-phi*180/TMath::Pi());   
      tgr.RotateY(90);

      const Double_t *rotMatPtr;

      rotMatPtr = tgr.GetRotationMatrix();
      TMatrixT<Double_t> rotMat(3,3, rotMatPtr);
      dd4hep::Rotation3D r3d;
      r3d.SetComponents(rotMat(0,0), rotMat(0,1), rotMat(0,2),
			rotMat(1,0), rotMat(1,1), rotMat(1,2),
			rotMat(2,0), rotMat(2,1), rotMat(2,2));

      //      dd4hep::Transform3D com(dd4hep::RotationZYX(PlateAngle, TMath::Pi()/2.,  -phi), dd4hep::Translation3D(x,y,z));
      dd4hep::Transform3D com(r3d, dd4hep::Translation3D(x,y,z));
      //      dd4hep::PlacedVolume absPlateVol_pv = aEnvelope.placeVolume(absPlateVol, com);
      //      absPlateVol_pv.addPhysVolID("module", iModule);
      //absPlateVol_pv.addPhysVolID("type", 1);

      tgr.Clear();
      tgr.RotateZ(PlateAngle*180/TMath::Pi());
      tgr.RotateX(-(phi+delPhi/2.)*180/TMath::Pi());   
      tgr.RotateY(90);

      rotMatPtr = tgr.GetRotationMatrix();
      TMatrixT<Double_t> rotMat2(3,3, rotMatPtr);

      r3d.SetComponents(rotMat2(0,0), rotMat2(0,1), rotMat2(0,2),
			rotMat2(1,0), rotMat2(1,1), rotMat2(1,2),
			rotMat2(2,0), rotMat2(2,1), rotMat2(2,2));

      //  x = (ri+(ro-ri)/2.)*TMath::Cos(phi+delPhi/2.);
      // y = (ri+(ro-ri)/2.)*TMath::Sin(phi+delPhi/2.);
      x = (-ro/2.)*TMath::Cos(phi+delPhi/2.);
      y = (-ro/2.)*TMath::Sin(phi+delPhi/2.);

      dd4hep::Transform3D com2(r3d, dd4hep::Translation3D(x,y,z));
      //      dd4hep::PlacedVolume electrodePlateVol_pv = aEnvelope.placeVolume(electrodePlateVol, com2);
      dd4hep::PlacedVolume LArVol_pv = aEnvelope.placeVolume(LArVol, com2);
      LArVol_pv.addPhysVolID("layer", 0);
      LArVol_pv.addPhysVolID("module", iModule);
      //      electrodePlateVol_pv.addPhysVolID("layer", 0);
      if (sign > 0) {
	LArVol_pv.addPhysVolID("cryo", 1);
      } else {
	LArVol_pv.addPhysVolID("cryo", 0);
      }
      LArVol_pv.addPhysVolID("type", 2);

      iModule++;

    }

    return;

  }

  void buildOneSide_Turbine(dd4hep::Detector& aLcdd, dd4hep::SensitiveDetector& aSensDet,
                  dd4hep::Volume& aEnvelope, dd4hep::DetElement& aEnvelopeDetElem, dd4hep::xml::Handle_t& aXmlElement,
                  int sign) {

  dd4hep::xml::Dimension dim(aXmlElement.child(_Unicode(dimensions)));

  dd4hep::xml::DetElement active = aXmlElement.child(_Unicode(active));
  std::string activeMaterial = active.materialStr();
  // double activeThickness = active.thickness();

  dd4hep::xml::DetElement readout = aXmlElement.child(_Unicode(readout));
  std::string readoutMaterial = readout.materialStr();
  //  double readoutThickness = readout.thickness();

  dd4hep::xml::DetElement passive = aXmlElement.child(_Unicode(passive));
  dd4hep::xml::DetElement passiveInner = passive.child(_Unicode(inner));
  dd4hep::xml::DetElement passiveOuter = passive.child(_Unicode(outer));
  dd4hep::xml::DetElement passiveGlue = passive.child(_Unicode(glue));
  double passiveInnerThickness = passiveInner.thickness();
  double passiveOuterThickness = passiveOuter.thickness();
  double passiveGlueThickness = passiveGlue.thickness();
  double passiveThickness = passiveInnerThickness + passiveOuterThickness + passiveGlueThickness;
  std::string passiveInnerMaterial = passiveInner.materialStr();
  std::string passiveOuterMaterial = passiveOuter.materialStr();
  std::string passiveGlueMaterial = passiveGlue.materialStr();
  std::string passiveMaterial;
  if (passiveInnerThickness < passiveThickness) {
    passiveMaterial = "Air";
  } else {
    passiveMaterial = passiveInnerMaterial;
  }
  
  double length = dim.dz() * 2.;
  double zOffsetEnvelope = length / 2. * -sign;
 
  //  dd4hep::Tube supportTube(100., 110., 24);
  //  dd4hep::DetElement supportTubeElem("supportTube", idDet);
  dd4hep::xml::DetElement supportTubeElem = aXmlElement.child(_Unicode(supportTube));
  float radiusRatio = supportTubeElem.attr<float>(_Unicode(radiusRatio));  
  double rmin = dim.rmin1();
  double rmax = dim.rmax();
  std::cout << "Min, max radii are " << rmin << " " << rmax << std::endl;
  double ro = rmin*radiusRatio;
  double ri = rmin;

  float supportTubeThickness=supportTubeElem.thickness();
  unsigned iSupportTube = 0;
  // std::vector<float> plateAngle_list;
  //std::vector<int> nPlates_list;

  unsigned iModule = 0;
  while (ri < rmax) {
    dd4hep::Tube supportTube(ro, ro+supportTubeThickness, dim.dz() );
  
    dd4hep::Volume supportTubeVol("supportTube", supportTube, aLcdd.material("Steel235"));
    dd4hep::PlacedVolume supportTube_pv = aEnvelope.placeVolume(supportTubeVol, dd4hep::Position(0,0,zOffsetEnvelope + sign * (dim.dz() )));
    supportTube_pv.addPhysVolID("supportTube", iSupportTube);
    buildSubCylinder(aLcdd, aSensDet, aEnvelope, aXmlElement, ri+supportTubeThickness, ro, 
		     iModule, sign);
    ri = ro;
    ro *= radiusRatio;
    if (ro > rmax) ro = rmax;
    iSupportTube++;
  }

  std::cout << "Total number of modules: " << iModule << std::endl;

  return;
}
  


  static dd4hep::Ref_t
createECalEndcapTurbine(dd4hep::Detector& aLcdd, dd4hep::xml::Handle_t aXmlElement, dd4hep::SensitiveDetector aSensDet) {

  dd4hep::xml::DetElement xmlDetElem = aXmlElement;
  std::string nameDet = xmlDetElem.nameStr();
  int idDet = xmlDetElem.id();
  dd4hep::xml::Dimension dim(xmlDetElem.dimensions());
  dd4hep::DetElement caloDetElem(nameDet, idDet);
  dd4hep::xml::Dimension sdType = xmlDetElem.child(_U(sensitive));
  aSensDet.setType(sdType.typeStr());
  
 
  // Create air envelope for the whole endcap
  dd4hep::Cone envelopePositive(dim.dz(), dim.rmin1(), dim.rmax(), dim.rmin2(), dim.rmax());
  dd4hep::Cone envelopeNegative(dim.dz(), dim.rmin2(), dim.rmax(), dim.rmin1(), dim.rmax());
  dd4hep::UnionSolid envelopeShape(envelopePositive, envelopeNegative, dd4hep::Position(0, 0, -2 * dim.z_offset()));
  dd4hep::Volume envelopeVol(nameDet + "_vol", envelopeShape, aLcdd.material("Air"));
  dd4hep::Volume envelopePositiveVol(nameDet + "_positive_vol", envelopePositive, aLcdd.material("Air"));
  dd4hep::Volume envelopeNegativeVol(nameDet + "_negative_vol", envelopeNegative, aLcdd.material("Air"));

  dd4hep::DetElement caloPositiveDetElem(caloDetElem, "positive", 0);
  dd4hep::DetElement caloNegativeDetElem(caloDetElem, "negative", 0);

  //  xml_det_t x_det = aXmlElement;
  // if (x_det.isSensitive()) {
  //  std::cout << "Yes I'm sensitive " << std::endl;
  //  envelopePositiveVol.setSensitiveDetector(aSensDet);
  //}

  lLog << MSG::DEBUG << "Placing dector on the positive side: (cm) " << dim.z_offset() << " with min, max radii " << dim.rmin1() << " " << dim.rmax() << endmsg;
  buildOneSide_Turbine(aLcdd, aSensDet, envelopePositiveVol, caloPositiveDetElem, aXmlElement, 1);
  lLog << MSG::DEBUG << "Placing dector on the negative side: (cm) " << -dim.z_offset()  << " with min, max radii " << dim.rmin1() << " " << dim.rmax() << endmsg;
  buildOneSide_Turbine(aLcdd, aSensDet, envelopeNegativeVol, caloNegativeDetElem, aXmlElement, -1);

  // Place the envelope
  dd4hep::PlacedVolume envelopePositivePhysVol = envelopeVol.placeVolume(envelopePositiveVol);
  //  envelopePositivePhysVol.addPhysVolID("subsystem", 0);
  caloPositiveDetElem.setPlacement(envelopePositivePhysVol);
  dd4hep::PlacedVolume envelopeNegativePhysVol =
      envelopeVol.placeVolume(envelopeNegativeVol, dd4hep::Position(0, 0, -2 * dim.z_offset()));
  //  envelopeNegativePhysVol.addPhysVolID("subsystem", 1);

  caloNegativeDetElem.setPlacement(envelopeNegativePhysVol);

  dd4hep::Volume motherVol = aLcdd.pickMotherVolume(caloDetElem);
  dd4hep::PlacedVolume envelopePhysVol = motherVol.placeVolume(envelopeVol, dd4hep::Position(0., 0., dim.z_offset()));
  caloDetElem.setPlacement(envelopePhysVol);
  envelopePhysVol.addPhysVolID("system", idDet);
  return caloDetElem;
}
}  // namespace det

DECLARE_DETELEMENT(ECalEndcap_Turbine_o1_v01, det::createECalEndcapTurbine)

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

  void buildSubCylinder(dd4hep::Detector& aLcdd,
			dd4hep::SensitiveDetector& aSensDet,
			dd4hep::Volume& aEnvelope,
			dd4hep::xml::Handle_t& aXmlElement,
			float ri, float ro,
			std::vector<float> &plateAngle_list,
			std::vector<int> &nPlates_list,
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

    int nSubSteps = genericBladeElem.attr<int>(_Unicode(nSubSteps));;
    Float_t xInit = -delZ/(2*TMath::Sin(PlateAngle));
    Float_t xRange = delZ/(TMath::Sin(PlateAngle));
    Float_t xStep = delZ/(TMath::Sin(PlateAngle))/nSubSteps;
    std::vector<Float_t> xVals, yVals, zVals;
    xVals.reserve(2*nSubSteps+8);
    yVals.reserve(2*nSubSteps+8);
    zVals.reserve(2*nSubSteps+8);

    int nPlates;
    double delrPhi;
    
    float AbsThicki = absBladeElem.attr<float>(_Unicode(thickness));
    float ElectrodeThick = electrodeBladeElem.attr<float>(_Unicode(thickness));
    float LArgapi = nobleLiquidElem.attr<float>(_Unicode(gap));
    
    bool sameNplates = genericBladeElem.attr<bool>(_Unicode(sameNplates));;
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

    bool scalePlateThickness = absBladeElem.attr<bool>(_Unicode(scaleThickness));
    float AbsThicko;
    if (scalePlateThickness) {
      AbsThicko = AbsThicki*(ro/ri);
    } else {
      AbsThicko = AbsThicki;
    }

    dd4hep::Trd2 absPlate(AbsThicki/2., AbsThicko/2., xRange/2., xRange/2., (ro-ri)/2. );
    dd4hep::Volume absPlateVol("absPlate", absPlate, aLcdd.material("Lead"));
    dd4hep::Trd2 electrodePlate(ElectrodeThick/2.+LArgapi, ElectrodeThick/2.+LArgapi, xRange/2., xRange/2., (ro-ri)/2. );
    dd4hep::Volume electrodePlateVol("electrodePlate", electrodePlate, aLcdd.material("Air"));
    
    electrodePlateVol.setSensitiveDetector(aSensDet);

    for (int iStep = 0; iStep < nSubSteps; iStep++) {
      Float_t z = 0;
      
      // convert to cartesian
      xVals[2*iStep] = xInit+iStep*xStep;
      yVals[2*iStep] = TMath::Sqrt(ri*ri-xVals[2*iStep]*xVals[2*iStep]*TMath::Cos(PlateAngle)*TMath::Cos(PlateAngle));
      zVals[2*iStep] = 0.;
      //    xVals[2*iStep+1] = ro*TMath::Cos(phi_o);
      //    yVals[2*iStep+1] = ro*TMath::Sin(phi_o);
      xVals[2*iStep+1] =  xVals[2*iStep];
      yVals[2*iStep+1] = TMath::Sqrt(ro*ro-xVals[2*iStep]*xVals[2*iStep]*TMath::Cos(PlateAngle)*TMath::Cos(PlateAngle));
      zVals[2*iStep+1] = 0.;

      char name[10];
      sprintf(name, "b%d", iStep);
      
      dd4hep::Trd2 b(AbsThicki/2., AbsThicko/2., xStep/2., xStep/2., (yVals[2*iStep+1]- yVals[2*iStep])/2.);
      dd4hep::Volume bVol("bar", b, aLcdd.material("Lead"));
      dd4hep::Transform3D tr1(dd4hep::RotationZYX(0.,0.,0), dd4hep::Translation3D(0, xVals[2*iStep+1], (yVals[2*iStep+1]+ yVals[2*iStep])/2.-(ro+ri)/2.));
      dd4hep::PlacedVolume b_pv = absPlateVol.placeVolume(bVol, tr1);
      //      b_pv.addPhysVolID("b", iStep);

      dd4hep::Trd2 b2(ElectrodeThick/2., ElectrodeThick/2., xStep/2., xStep/2., (yVals[2*iStep+1]- yVals[2*iStep])/2.);
      dd4hep::Volume b2Vol("bar2", b2, aLcdd.material("PCB"));     
      dd4hep::PlacedVolume b2_pv = electrodePlateVol.placeVolume(b2Vol, tr1);
      // b2_pv.addPhysVolID("b2", iStep); 

      dd4hep::Trd2 b3(ElectrodeThick/2.+LArgapi, ElectrodeThick/2.+LArgapi, xStep/2., xStep/2., (yVals[2*iStep+1]- yVals[2*iStep])/2.);
      dd4hep::SubtractionSolid b4(b3, b2);
      dd4hep::Volume b4Vol("bar4", b4, aLcdd.material("LAr"));     
      dd4hep::PlacedVolume b4_pv = electrodePlateVol.placeVolume(b4Vol, tr1);

    }

    int    nPlatesToDraw = nPlates;
    //    nPlatesToDraw = 1;

    for (int iPlate = 0; iPlate < nPlatesToDraw; iPlate++) {
      //    pl->SetInvisible();
      float phi = iPlate*2*TMath::Pi()/nPlates;
      float delPhi = 2*TMath::Pi()/nPlates;
      
      float x = (ri+(ro-ri)/2.)*TMath::Cos(phi);
      float y = (ri+(ro-ri)/2.)*TMath::Sin(phi); //ri*TMath::Sin(phi)/6.;
      float z = xStep/2.; //ri*TMath::Cos(phi)/6.;
      z = 0.;
      
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
      dd4hep::PlacedVolume absPlateVol_pv = aEnvelope.placeVolume(absPlateVol, com);
      absPlateVol_pv.addPhysVolID("module", iModule);
      absPlateVol_pv.addPhysVolID("type", 1);

      tgr.Clear();
      tgr.RotateZ(PlateAngle*180/TMath::Pi());
      tgr.RotateX(-(phi+delPhi/2.)*180/TMath::Pi());   
      tgr.RotateY(90);

      rotMatPtr = tgr.GetRotationMatrix();
      TMatrixT<Double_t> rotMat2(3,3, rotMatPtr);

      r3d.SetComponents(rotMat2(0,0), rotMat2(0,1), rotMat2(0,2),
			rotMat2(1,0), rotMat2(1,1), rotMat2(1,2),
			rotMat2(2,0), rotMat2(2,1), rotMat2(2,2));

      x = (ri+(ro-ri)/2.)*TMath::Cos(phi+delPhi/2.);
      y = (ri+(ro-ri)/2.)*TMath::Sin(phi+delPhi/2.);
      

      dd4hep::Transform3D com2(r3d, dd4hep::Translation3D(x,y,z));
      dd4hep::PlacedVolume electrodePlateVol_pv = aEnvelope.placeVolume(electrodePlateVol, com2);
      //      electrodePlateVol_pv.addPhysVolID("electrodePlate", iPlate);
      electrodePlateVol_pv.addPhysVolID("module", iModule);
      electrodePlateVol_pv.addPhysVolID("layer", 0);
      if (sign > 0) {
	electrodePlateVol_pv.addPhysVolID("cryo", 1);
      } else {
	electrodePlateVol_pv.addPhysVolID("cryo", 0);
      }
      electrodePlateVol_pv.addPhysVolID("type", 2);

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
  double activeThickness = active.thickness();

  dd4hep::xml::DetElement readout = aXmlElement.child(_Unicode(readout));
  std::string readoutMaterial = readout.materialStr();
  double readoutThickness = readout.thickness();

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
  double ro = rmin*radiusRatio;
  double ri = rmin;

  float supportTubeThickness=supportTubeElem.thickness();
  unsigned iSupportTube = 0;
  std::vector<float> plateAngle_list;
  std::vector<int> nPlates_list;

  unsigned iModule = 0;
  while (ri < rmax) {
    dd4hep::Tube supportTube(ro, ro+supportTubeThickness, dim.dz() );
  
    dd4hep::Volume supportTubeVol("supportTube", supportTube, aLcdd.material("Steel235"));
    dd4hep::PlacedVolume supportTube_pv = aEnvelope.placeVolume(supportTubeVol, dd4hep::Position(0,0,zOffsetEnvelope + sign * (dim.dz() )));
    supportTube_pv.addPhysVolID("supportTube", iSupportTube);
    buildSubCylinder(aLcdd, aSensDet, aEnvelope, aXmlElement, ri, ro, 
		     plateAngle_list, nPlates_list, iModule, sign);
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

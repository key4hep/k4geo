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

  double tForArcLength(double s, double bladeangle, double delZ, double r) {

    // some intermediate constants
    double zpos = delZ/2.;
    double zp = zpos/TMath::Tan(bladeangle);
    double b = zp/(TMath::Sqrt(r*r-zp*zp));
    double c = (TMath::Tan(s/r) +b)/(1.-b*TMath::Tan(s/r));
    double d = c*c*r*r/(1+c*c);
    return (TMath::Sqrt(d)-zp)*TMath::Sin(bladeangle);

  }

  dd4hep::Solid buildOneBlade(double thickness_inner,
			      double thickness_outer,
			      double width,
			      double ro, double ri,
			      double bladeangle,
			      double delZ) 
  {

    // extrapolate inner thickness down to r = 0
    double t0 = (thickness_outer*ri-thickness_inner*ro)/(ri-ro);
    double r0 = 0.;
    dd4hep::Solid shapeBeforeSubtraction;
    if (t0 < 0.) { // now we need to adjust what we call r0...
      r0 = (thickness_outer*ri-thickness_inner*ro)/(thickness_outer-thickness_inner);
      t0 = 0.;
      dd4hep::Trd2 tmp1(t0/2.+0.01, thickness_outer/2., width/2., width/2., (ro-r0)/2. );
      dd4hep::Trd2 tmp2(thickness_inner/2., thickness_inner/2., width/2., width/2., ro/2. );
      dd4hep::UnionSolid tmp3(tmp2, tmp1, dd4hep::Transform3D(dd4hep::RotationZYX( 0, 0, 0),dd4hep::Position(0, 0, r0)));
      shapeBeforeSubtraction = tmp3;
      lLog << MSG::INFO << "extrapolated thickness is " << t0 <<  " and r0 is " << r0 << endmsg;
    } else {
      dd4hep::Trd2 tmp1(t0/2., thickness_outer/2., width/2., width/2., ro/2. );
      shapeBeforeSubtraction = tmp1;
    }

    dd4hep::Tube outerTube(ro, ro*100, delZ);
    dd4hep::Tube innerTube(0., ri, delZ); 
    dd4hep::SubtractionSolid tmp4(shapeBeforeSubtraction, outerTube, dd4hep::Transform3D(dd4hep::RotationZYX( 0, TMath::Pi()/2.-bladeangle, TMath::Pi()/2.),dd4hep::Position(0, 0, -ro/2.)));
    return dd4hep::SubtractionSolid(tmp4, innerTube, dd4hep::Transform3D(dd4hep::RotationZYX( 0 , TMath::Pi()/2.-bladeangle, TMath::Pi()/2.),dd4hep::Position(0, 0, -ro/2.)));
    //return shapeBeforeSubtraction;

  }
			      
  void buildSubCylinder(dd4hep::Detector& aLcdd,
			dd4hep::SensitiveDetector& aSensDet,
			dd4hep::Volume& aEnvelope,
			dd4hep::xml::Handle_t& aXmlElement,
			dd4hep::DetElement& bathDetElem,
			float ri, float ro,
			unsigned &iModule) {


    dd4hep::xml::DetElement genericBladeElem = aXmlElement.child(_Unicode(turbineBlade));
    dd4hep::xml::DetElement absBladeElem = genericBladeElem.child(_Unicode(absorberBlade));
    dd4hep::xml::DetElement electrodeBladeElem = genericBladeElem.child(_Unicode(electrodeBlade));
    dd4hep::xml::DetElement nobleLiquidElem = genericBladeElem.child(_Unicode(nobleLiquidGap));
    

    float BladeAngle = genericBladeElem.attr<float>(_Unicode(angle));
    bool decreaseAnglePerSubCylinder = genericBladeElem.attr<bool>(_Unicode(decreaseAnglePerSubCylinder));
    std::cout << "Blade angle is " << BladeAngle << "; decrease angle per subcylinder? " << decreaseAnglePerSubCylinder << std::endl;
    dd4hep::xml::Dimension dim(aXmlElement.child(_Unicode(dimensions)));
    double grmin = dim.rmin1();
    double delZ = dim.dz()*2;
    std::cout << "delZ is " << delZ << std::endl;
    if (decreaseAnglePerSubCylinder) {
      float tubeFracCovered = delZ/(2*grmin*TMath::Tan(BladeAngle));
      BladeAngle = TMath::ATan(delZ/(2*ri*tubeFracCovered));
    }
    
    if (TMath::Abs(TMath::Tan(BladeAngle)) < delZ/(2.*ri)) {
      std::cout << "The requested blade angle is too small for the given delZ and ri values.  Please adjust to at least " << TMath::ATan(delZ/(2.*ri))*180./TMath::Pi() << " degrees!" << std::endl;
      return;
    }

    Float_t xRange = delZ/(TMath::Sin(BladeAngle));

    int nBlades;
    double delrPhi, delrPhiNoGap;
    
    float AbsThicki = absBladeElem.attr<float>(_Unicode(thickness));
    float ElectrodeThick = electrodeBladeElem.attr<float>(_Unicode(thickness));
    float LArgapi = nobleLiquidElem.attr<float>(_Unicode(gap));
    
    bool sameNblades = genericBladeElem.attr<bool>(_Unicode(sameNblades));;

    bool scaleBladeThickness = absBladeElem.attr<bool>(_Unicode(scaleThickness));
    float AbsThicko;
    if (scaleBladeThickness) {
      AbsThicko = AbsThicki*(ro/ri);
    } else {
      AbsThicko = AbsThicki;
    }

    if (sameNblades) {
      nBlades = 2*TMath::Pi()*grmin/(2*LArgapi+AbsThicki+ElectrodeThick)*TMath::Sin(BladeAngle);
    } else {
      double x1 = delZ/2./TMath::Tan(BladeAngle);
      double x2 = delZ/2./TMath::Tan(BladeAngle)+(2*LArgapi+AbsThicki+ElectrodeThick)/TMath::Sin(BladeAngle);
      double y1 = TMath::Sqrt(ri*ri-x1*x1);
      double y2 = TMath::Sqrt(ri*ri-x2*x2);
      double rPhi1 = ri*TMath::ATan(x1/y1);
      double rPhi2 = ri*TMath::ATan(x2/y2);
      delrPhi = TMath::Abs(rPhi1-rPhi2);
      nBlades = 2*TMath::Pi()*ri/delrPhi;

    }

    std::cout << "nBlades: " << nBlades << std::endl;
    // adjust gap thickness at inner layer
    double circ = 2*TMath::Pi()*ri;
    double x1 = delZ/2./TMath::Tan(BladeAngle);
    double x2 = delZ/2./TMath::Tan(BladeAngle)+(AbsThicki+ElectrodeThick)/TMath::Sin(BladeAngle);
    double y1 = TMath::Sqrt(ri*ri-x1*x1);
    double y2 = TMath::Sqrt(ri*ri-x2*x2);
    double rPhi1 = ri*TMath::ATan(x1/y1);
    double rPhi2 = ri*TMath::ATan(x2/y2);    
    delrPhiNoGap = TMath::Abs(rPhi1-rPhi2);
    double leftoverS = (circ - nBlades*delrPhiNoGap);
    double delrPhiGapOnly = leftoverS/(2*nBlades);
    std::cout << "LAr gap was " << LArgapi ;
    std::cout << "fixed s is " << delrPhiNoGap << std::endl;
    LArgapi = tForArcLength(delrPhiGapOnly, BladeAngle, delZ, ri);
    std::cout << " but is now " << LArgapi << " since s = " << delrPhiGapOnly << " s for full unit cell is " << delrPhi << std::endl;

    // now find gap at outer layer
    circ = 2*TMath::Pi()*ro;
    x1 = delZ/2./TMath::Tan(BladeAngle);
    x2 = delZ/2./TMath::Tan(BladeAngle)+(AbsThicko+ElectrodeThick)/TMath::Sin(BladeAngle);
    y1 = TMath::Sqrt(ro*ro-x1*x1);
    y2 = TMath::Sqrt(ro*ro-x2*x2);
    rPhi1 = ro*TMath::ATan(x1/y1);
    rPhi2 = ro*TMath::ATan(x2/y2);    
    delrPhiNoGap = TMath::Abs(rPhi1-rPhi2);
    leftoverS = (circ - nBlades*delrPhiNoGap);
    delrPhiGapOnly = leftoverS/(2*nBlades);   
    double LArgapo = tForArcLength(delrPhiGapOnly, BladeAngle, delZ, ro);
    std::cout << "Outer LAr gap is " << LArgapo << std::endl ;
    
    dd4hep::Solid  absBlade = buildOneBlade(AbsThicki, AbsThicko, xRange, ro, ri, BladeAngle, delZ );
    dd4hep::Volume absBladeVol("absBlade", absBlade, aLcdd.material(absBladeElem.materialStr()));
    lLog << MSG::INFO << "Inner and outer thicknesses of noble liquid volume " << ElectrodeThick+LArgapi*2 <<  " " <<  ElectrodeThick+LArgapo*2 << endmsg;
    dd4hep::Solid electrodeBladeAndGap = buildOneBlade(ElectrodeThick+LArgapi*2, ElectrodeThick+LArgapo*2, xRange, ro, ri, BladeAngle, delZ);
    dd4hep::Volume electrodeBladeAndGapVol("electrodeBladeAndGap", electrodeBladeAndGap, aLcdd.material("Air"));
     dd4hep::Solid electrodeBlade = buildOneBlade(ElectrodeThick, ElectrodeThick, xRange, ro, ri, BladeAngle, delZ);
     dd4hep::Volume electrodeBladeVol("electrodeBlade", electrodeBlade, aLcdd.material(electrodeBladeElem.materialStr()));
    dd4hep::SubtractionSolid LArShapeTotal(electrodeBladeAndGap, electrodeBlade);
    unsigned nLayers = aLcdd.constant<unsigned>("ECalEndcapNumLayers");

    std::vector<dd4hep::Volume> LArVolLayers;
    double xRangePerLayer = xRange/nLayers;
    dd4hep::Trd2 bigDumbBlock(ElectrodeThick/2.+LArgapo, ElectrodeThick/2.+LArgapo, xRange/2, xRange/2., ro/2. );
    for (unsigned iLayer = 0; iLayer < nLayers; iLayer++) {
      // make a layer by subtracting off the unwanted bits
      dd4hep::SubtractionSolid tmp1(LArShapeTotal, bigDumbBlock, dd4hep::Position( 0., -xRange+xRangePerLayer*iLayer, 0.));
      dd4hep::SubtractionSolid tmp2(tmp1, bigDumbBlock,  dd4hep::Position(0., xRangePerLayer*(1+iLayer), 0.));
      dd4hep::Volume LArVolThisLayer("nobleLiquidLayer", tmp2, aLcdd.material(nobleLiquidElem.materialStr()));
      LArVolThisLayer.setSensitiveDetector(aSensDet);
      LArVolLayers.push_back(LArVolThisLayer);
    }

    lLog << MSG::INFO << "ECal endcap materials:  nobleLiquid: " << nobleLiquidElem.materialStr() << " absorber: " << absBladeElem.materialStr() << " electrode: " << electrodeBladeElem.materialStr() << endmsg; 
  //build cryostat
// Retrieve cryostat data
    int    nBladesToDraw = nBlades;
    //    nBladesToDraw = 1;

    for (int iBlade = 0; iBlade < nBladesToDraw; iBlade++) {
      //    pl->SetInvisible();
      float phi = iBlade*2*TMath::Pi()/nBlades;
      float delPhi = 2*TMath::Pi()/nBlades;
      
      float x = (ro/2.)*TMath::Cos(phi);
      float y = (ro/2.)*TMath::Sin(phi); //ri*TMath::Sin(phi)/6.;
      float z =  0.;
      
      TGeoRotation tgr;
      tgr.RotateZ(BladeAngle*180/TMath::Pi());
      tgr.RotateX(-phi*180/TMath::Pi());   
      tgr.RotateY(90);

      const Double_t *rotMatPtr;

      rotMatPtr = tgr.GetRotationMatrix();
      TMatrixT<Double_t> rotMat(3,3, rotMatPtr);
      dd4hep::Rotation3D r3d;
      r3d.SetComponents(rotMat(0,0), rotMat(0,1), rotMat(0,2),
			rotMat(1,0), rotMat(1,1), rotMat(1,2),
			rotMat(2,0), rotMat(2,1), rotMat(2,2));

      dd4hep::Transform3D com(r3d, dd4hep::Translation3D(x,y,z));
      dd4hep::PlacedVolume absBladeVol_pv = aEnvelope.placeVolume(absBladeVol, com);
      
      absBladeVol_pv.addPhysVolID("module", iModule);
      absBladeVol_pv.addPhysVolID("type", 1);  // 0 = active, 1 = passive, 2 = readout
      dd4hep::DetElement absBladeDetElem(bathDetElem, "absorber"+std::to_string(iModule), iModule);
      absBladeDetElem.setPlacement(absBladeVol_pv);

      tgr.Clear();
      tgr.RotateZ(BladeAngle*180/TMath::Pi());
      tgr.RotateX(-(phi+delPhi/2.)*180/TMath::Pi());   
      tgr.RotateY(90);

      rotMatPtr = tgr.GetRotationMatrix();
      TMatrixT<Double_t> rotMat2(3,3, rotMatPtr);

      r3d.SetComponents(rotMat2(0,0), rotMat2(0,1), rotMat2(0,2),
			rotMat2(1,0), rotMat2(1,1), rotMat2(1,2),
			rotMat2(2,0), rotMat2(2,1), rotMat2(2,2));

      x = (ro/2.)*TMath::Cos(phi+delPhi/2.);
      y = (ro/2.)*TMath::Sin(phi+delPhi/2.);

      dd4hep::Transform3D com2(r3d, dd4hep::Translation3D(x,y,z));
      dd4hep::PlacedVolume electrodeBladeVol_pv = aEnvelope.placeVolume(electrodeBladeVol, com2);
      electrodeBladeVol_pv.addPhysVolID("module", iModule);
      electrodeBladeVol_pv.addPhysVolID("type", 2);  // 0 = active, 1 = passive, 2 = readout
      dd4hep::DetElement electrodeBladeDetElem(bathDetElem, "electrode"+std::to_string(iModule), iModule);
      electrodeBladeDetElem.setPlacement(electrodeBladeVol_pv);

      for (unsigned iLayer = 0; iLayer<nLayers; iLayer++) {
	dd4hep::PlacedVolume LArVol_pv = aEnvelope.placeVolume(LArVolLayers[iLayer], com2);
	LArVol_pv.addPhysVolID("layer", iLayer);
	LArVol_pv.addPhysVolID("module", iModule);
	LArVol_pv.addPhysVolID("type", 0);  // 0 = active, 1 = passive, 2 = readout
	dd4hep::DetElement LArDetElem(bathDetElem, "active"+std::to_string(iModule)+"layer"+std::to_string(iLayer), (iModule<<8)+iLayer);
	LArDetElem.setPlacement(LArVol_pv);

      }
      iModule++;

    }

    return;

  }

  void buildOneSide_Turbine(dd4hep::Detector& aLcdd, dd4hep::SensitiveDetector& aSensDet,
                  dd4hep::Volume& aEnvelope, dd4hep::xml::Handle_t& aXmlElement,
			    int sign,
			    unsigned& iModule) {

    dd4hep::xml::DetElement calo = aXmlElement.child(_Unicode(calorimeter));
    dd4hep::xml::Dimension caloDim(calo.dimensions());

    dd4hep::xml::DetElement blade = aXmlElement.child(_Unicode(turbineBlade));
    dd4hep::xml::DetElement nobleLiquid = blade.child(_Unicode(nobleLiquidGap));

    dd4hep::xml::DetElement xmlDetElem = aXmlElement;
    std::string nameDet = xmlDetElem.nameStr();
    dd4hep::DetElement caloDetElem(nameDet, xmlDetElem.id());

    dd4hep::xml::Dimension dim(aXmlElement.child(_Unicode(dimensions)));

  //build cryostat
// Retrieve cryostat data
  dd4hep::xml::DetElement cryostat = aXmlElement.child(_Unicode(cryostat));
  dd4hep::xml::Dimension cryoDim(cryostat.dimensions());
  double cryoThicknessFront = cryoDim.rmin2() - cryoDim.rmin1();

  dd4hep::xml::DetElement cryoFront = cryostat.child(_Unicode(front));
  dd4hep::xml::DetElement cryoBack = cryostat.child(_Unicode(back));
  dd4hep::xml::DetElement cryoSide = cryostat.child(_Unicode(side));
  bool cryoFrontSensitive = cryoFront.isSensitive();
  bool cryoBackSensitive = cryoBack.isSensitive();
  bool cryoSideSensitive = cryoSide.isSensitive();

  double bathRmin = caloDim.rmin(); // - margin for inclination
  double bathRmax = caloDim.rmax(); // + margin for inclination
  dd4hep::Tube bathOuterShape(bathRmin, bathRmax, caloDim.dz()); // make it 4 volumes + 5th for detector envelope
  dd4hep::Tube bathAndServicesOuterShape(cryoDim.rmin2(), cryoDim.rmax1(), caloDim.dz()); // make it 4 volumes + 5th for detector envelope

  lLog << MSG::INFO << "Cryostat front thickness is " <<  cryoDim.rmin2() << endmsg;
  if (cryoThicknessFront > 0) {
    // 1. Create cryostat
    dd4hep::Tube cryoFrontShape(cryoDim.rmin1(), cryoDim.rmin2(), cryoDim.dz());
    dd4hep::Tube cryoBackShape(cryoDim.rmax1(), cryoDim.rmax2(), cryoDim.dz());
    dd4hep::Tube cryoSideOuterShape(cryoDim.rmin2(), cryoDim.rmax1(), cryoDim.dz());
    dd4hep::SubtractionSolid cryoSideShape(cryoSideOuterShape, bathAndServicesOuterShape);
    lLog << MSG::INFO << "ECAL endcap cryostat: front: rmin (cm) = " << cryoDim.rmin1() << " rmax (cm) = " << cryoDim.rmin2() << " dz (cm) = " << cryoDim.dz()  << endmsg;
    lLog << MSG::INFO << "ECAL encdap cryostat: back: rmin (cm) = " << cryoDim.rmax1() << " rmax (cm) = " << cryoDim.rmax2() << " dz (cm) = " << cryoDim.dz() << endmsg;
    lLog << MSG::INFO << "ECAL endcap cryostat: side: rmin (cm) = " << cryoDim.rmin2() << " rmax (cm) = " << cryoDim.rmax1() << " dz (cm) = " << cryoDim.dz() - caloDim.dz()  << endmsg;
    dd4hep::Volume cryoFrontVol(cryostat.nameStr()+"_front", cryoFrontShape, aLcdd.material(cryostat.materialStr()));
    dd4hep::Volume cryoBackVol(cryostat.nameStr()+"_back", cryoBackShape, aLcdd.material(cryostat.materialStr()));
    dd4hep::Volume cryoSideVol(cryostat.nameStr()+"_side", cryoSideShape, aLcdd.material(cryostat.materialStr()));
    dd4hep::PlacedVolume cryoFrontPhysVol = aEnvelope.placeVolume(cryoFrontVol);
    dd4hep::PlacedVolume cryoBackPhysVol = aEnvelope.placeVolume(cryoBackVol);
    dd4hep::PlacedVolume cryoSidePhysVol = aEnvelope.placeVolume(cryoSideVol);
    unsigned sidetype = sign > 0 ? 0x4 : 0x0; 
    if (cryoFrontSensitive) {
      cryoFrontVol.setSensitiveDetector(aSensDet);
      cryoFrontPhysVol.addPhysVolID("cryo", 1);
      cryoFrontPhysVol.addPhysVolID("type", sidetype+1);
      lLog << MSG::INFO << "Cryostat front volume set as sensitive" << endmsg;
    }
    if (cryoBackSensitive) {
      cryoBackVol.setSensitiveDetector(aSensDet);
      cryoBackPhysVol.addPhysVolID("cryo", 1);
      cryoBackPhysVol.addPhysVolID("type", sidetype+2);
      lLog << MSG::INFO << "Cryostat back volume set as sensitive" << endmsg;
    }
    if (cryoSideSensitive) {
      cryoSideVol.setSensitiveDetector(aSensDet);
      cryoSidePhysVol.addPhysVolID("cryo", 1);
      cryoSidePhysVol.addPhysVolID("type", sidetype+3);
      lLog << MSG::INFO << "Cryostat front volume set as sensitive" << endmsg;
    }
    dd4hep::DetElement cryoFrontDetElem(caloDetElem, "cryo_front", 0);
    cryoFrontDetElem.setPlacement(cryoFrontPhysVol);
    dd4hep::DetElement cryoBackDetElem(caloDetElem, "cryo_back", 0);
    cryoBackDetElem.setPlacement(cryoBackPhysVol);
    dd4hep::DetElement cryoSideDetElem(caloDetElem, "cryo_side", 0);
    cryoSideDetElem.setPlacement(cryoSidePhysVol);
  }
						     
  // 2. Create noble liquid bath
  std::string nobleLiquidMaterial = nobleLiquid.materialStr();
  dd4hep::Volume bathVol(nobleLiquidMaterial + "_bath", bathOuterShape, aLcdd.material(nobleLiquidMaterial));
  lLog << MSG::INFO << "ECAL endcap bath: material = " << nobleLiquidMaterial << " rmin (cm) =  " << bathRmin
       << " rmax (cm) = " << bathRmax << " thickness in front of ECal (cm) = " << caloDim.rmin() - cryoDim.rmin2()
       << " thickness behind ECal (cm) = " << cryoDim.rmax1() - caloDim.rmax() << endmsg;
  dd4hep::DetElement bathDetElem(caloDetElem, "bath", 1);

  // 3. Create detector structure
  double length = dim.dz() * 2.;
  double zOffsetEnvelope = length / 2. * -sign;
 
  dd4hep::xml::DetElement supportTubeElem = aXmlElement.child(_Unicode(supportTube));
  float radiusRatio = supportTubeElem.attr<float>(_Unicode(radiusRatio));  
  double rmin = dim.rmin1();
  double rmax = dim.rmax();
  std::cout << "Min, max radii are " << rmin << " " << rmax << std::endl;
  double ro = rmin*radiusRatio;
  double ri = rmin;

  float supportTubeThickness=supportTubeElem.thickness();
  unsigned iSupportTube = 0;

  while (ri < rmax) {
    dd4hep::Tube supportTube(ro, ro+supportTubeThickness, dim.dz() );
  
    dd4hep::Volume supportTubeVol("supportTube", supportTube, aLcdd.material(supportTubeElem.materialStr()));
    dd4hep::PlacedVolume supportTube_pv = bathVol.placeVolume(supportTubeVol, dd4hep::Position(0,0,zOffsetEnvelope + sign * (dim.dz() )));
    buildSubCylinder(aLcdd, aSensDet, bathVol, aXmlElement, bathDetElem, ri+supportTubeThickness, ro, iModule);
    ri = ro;
    ro *= radiusRatio;
    if (ro > rmax) ro = rmax;
    iSupportTube++;
  }


  dd4hep::PlacedVolume bathPhysVol = aEnvelope.placeVolume(bathVol);
  bathDetElem.setPlacement(bathPhysVol);

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

  lLog << MSG::DEBUG << "Placing dector on the positive side: (cm) " << dim.z_offset() << " with min, max radii " << dim.rmin1() << " " << dim.rmax() << endmsg;
  unsigned iModule = 0;
  buildOneSide_Turbine(aLcdd, aSensDet, envelopePositiveVol,  aXmlElement, 1, iModule);
  lLog << MSG::DEBUG << "Placing dector on the negative side: (cm) " << -dim.z_offset()  << " with min, max radii " << dim.rmin1() << " " << dim.rmax() << endmsg;
  buildOneSide_Turbine(aLcdd, aSensDet, envelopeNegativeVol, aXmlElement, -1, iModule);

  // Place the envelope
  dd4hep::PlacedVolume envelopePositivePhysVol = envelopeVol.placeVolume(envelopePositiveVol);

  caloPositiveDetElem.setPlacement(envelopePositivePhysVol);
  dd4hep::PlacedVolume envelopeNegativePhysVol =
      envelopeVol.placeVolume(envelopeNegativeVol, dd4hep::Position(0, 0, -2 * dim.z_offset()));

  caloNegativeDetElem.setPlacement(envelopeNegativePhysVol);

  dd4hep::Volume motherVol = aLcdd.pickMotherVolume(caloDetElem);
  dd4hep::PlacedVolume envelopePhysVol = motherVol.placeVolume(envelopeVol, dd4hep::Position(0., 0., dim.z_offset()));
  caloDetElem.setPlacement(envelopePhysVol);
  envelopePhysVol.addPhysVolID("system", idDet);
  return caloDetElem;
}
}  // namespace det

DECLARE_DETELEMENT(ECalEndcap_Turbine_o1_v01, det::createECalEndcapTurbine)

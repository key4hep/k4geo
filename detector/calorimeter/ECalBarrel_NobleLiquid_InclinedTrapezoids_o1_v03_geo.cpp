#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Handle.h"
#include "DD4hep/Printout.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Vector3D.h"
#include "XML/Utilities.h"
#include <DDRec/DetectorData.h>

#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"

// like v02, but in xml the layer dimensions are along the electrode
// directions, without rescaling by R/L where R = radial extent of ECAL
// and L = electrode length
// in addition, the user can configure material in inner part of absorber
// in first layer

namespace det {
static dd4hep::detail::Ref_t createECalBarrelInclined(dd4hep::Detector& aLcdd, dd4hep::xml::Handle_t aXmlElement,
                                                      dd4hep::SensitiveDetector aSensDet) {

  dd4hep::xml::DetElement xmlDetElem = aXmlElement;
  std::string nameDet = xmlDetElem.nameStr();
  dd4hep::xml::Dimension dim(xmlDetElem.dimensions());
  dd4hep::DetElement caloDetElem(nameDet, xmlDetElem.id());

  // Create air envelope for the whole barrel
  dd4hep::Volume envelopeVol(nameDet + "_vol", dd4hep::Tube(dim.rmin(), dim.rmax(), dim.dz()), aLcdd.material("Air"));
  envelopeVol.setVisAttributes(aLcdd, dim.visStr());

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

  // Retrieve active and passive material data
  dd4hep::xml::DetElement calo = aXmlElement.child(_Unicode(calorimeter));
  dd4hep::xml::Dimension caloDim(calo.dimensions());
  dd4hep::xml::DetElement active = calo.child(_Unicode(active));
  std::string activeMaterial = active.materialStr();
  double activeThickness = active.thickness();

  // Retrieve information about active/passive overlap
  dd4hep::xml::DetElement overlap = active.child(_Unicode(overlap));
  double activePassiveOverlap = overlap.offset();
  if (activePassiveOverlap < 0 || activePassiveOverlap > 0.5) {
    dd4hep::printout(dd4hep::ERROR, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Overlap between active and passive cannot be more than half of passive plane!");
    throw std::runtime_error("Overlap too large!");
  }

  // Retrieve length of layers along electrode
  dd4hep::xml::DetElement layers = calo.child(_Unicode(layers));
  uint numLayers = 0;
  std::vector<double> layerHeight;
  double layersTotalHeight = 0;
  for (dd4hep::xml::Collection_t layer_coll(layers, _Unicode(layer)); layer_coll; ++layer_coll) {
    dd4hep::xml::Component layer = layer_coll;
    numLayers += layer.repeat();
    for (int iLay = 0; iLay < layer.repeat(); iLay++) {
      layerHeight.push_back(layer.thickness());
    }
    layersTotalHeight += layer.repeat() * layer.thickness();
  }
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "Total electrode length from calorimeter xml description (cm): %f ", layersTotalHeight / dd4hep::cm);

  // The following code checks if the xml geometry file contains a constant defining
  // the number of layers the barrel. In that case, it makes the program abort
  // if the number of planes in the xml is different from the one calculated from
  // the geometry. This is because the number of layers is needed
  // in other parts of the code (the readout for the FCC-ee ECAL with
  // inclined modules).
  int nLayers = -1;
  try {
    nLayers = aLcdd.constant<int>("ECalBarrelNumLayers");
  } catch (...) {
    ;
  }
  if (nLayers > 0 && nLayers != int(numLayers)) {
    dd4hep::printout(dd4hep::ERROR, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Incorrect number of layers (ECalBarrelNumLayers = %d) in readout in calorimeter xml description!",
                     nLayers);
    dd4hep::printout(dd4hep::ERROR, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Number of layers should be: %d", numLayers);
    throw std::runtime_error("Incorrect number of layers (ECalBarrelNumLayers) in calorimeter xml description!");
  }

  // Retrieve information about the readout
  dd4hep::xml::DetElement readout = calo.child(_Unicode(readout));
  std::string readoutMaterial = readout.materialStr();
  double readoutThickness = readout.thickness();

  // Retrieve info about passive elements
  // (glue, outer - typically steel, inner - typically Pb, and LAr in presampling layer)
  dd4hep::xml::DetElement passive = calo.child(_Unicode(passive));
  dd4hep::xml::DetElement passiveInner = passive.child(_Unicode(inner));
  dd4hep::xml::DetElement passiveInnerMax = passive.child(_Unicode(innerMax));
  dd4hep::xml::DetElement passiveOuter = passive.child(_Unicode(outer));
  dd4hep::xml::DetElement passiveGlue = passive.child(_Unicode(glue));
  std::string passiveInnerMaterial = passiveInnerMax.materialStr();
  std::string passiveInnerMaterialFirstLayer = passiveInner.materialStr();
  std::string passiveOuterMaterial = passiveOuter.materialStr();
  std::string passiveGlueMaterial = passiveGlue.materialStr();
  double passiveInnerThicknessMin = passiveInner.thickness();
  double passiveInnerThicknessMax = passiveInnerMax.thickness();
  double passiveOuterThickness = passiveOuter.thickness();
  double passiveGlueThickness = passiveGlue.thickness();
  double passiveThickness = passiveInnerThicknessMin + passiveOuterThickness + passiveGlueThickness;
  // inclination angle
  double angle = passive.rotation().angle();

  // Retrieve info about bath
  dd4hep::xml::DetElement bath = aXmlElement.child(_Unicode(bath));
  dd4hep::xml::Dimension bathDim(bath.dimensions());
  double bathRmin = bathDim.rmin();
  double bathRmax = bathDim.rmax();
  dd4hep::Tube bathOuterShape(bathRmin, bathRmax, caloDim.dz()); // make it 4 volumes + 5th for detector envelope
  dd4hep::Tube bathAndServicesOuterShape(cryoDim.rmin2(), cryoDim.rmax1(),
                                         caloDim.dz()); // make it 4 volumes + 5th for detector envelope

  // 1. Create cryostat
  if (cryoThicknessFront > 0) {
    dd4hep::Tube cryoFrontShape(cryoDim.rmin1(), cryoDim.rmin2(), cryoDim.dz());
    dd4hep::Tube cryoBackShape(cryoDim.rmax1(), cryoDim.rmax2(), cryoDim.dz());
    dd4hep::Tube cryoSideOuterShape(cryoDim.rmin2(), cryoDim.rmax1(), cryoDim.dz());
    dd4hep::SubtractionSolid cryoSideShape(cryoSideOuterShape, bathAndServicesOuterShape);

    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "ECAL cryostat: front: rmin (cm) = %f  rmax (cm) = %f  dz (cm) = %f", cryoDim.rmin1() / dd4hep::cm,
                     cryoDim.rmin2() / dd4hep::cm, cryoDim.dz() / dd4hep::cm);
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "ECAL cryostat: back: rmin (cm) = %f  rmax (cm) = %f  dz (cm) = %f", cryoDim.rmax1() / dd4hep::cm,
                     cryoDim.rmax2() / dd4hep::cm, cryoDim.dz() / dd4hep::cm);
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "ECAL cryostat: side: rmin (cm) = %f  rmax (cm) = %f  dz (cm) = %f", cryoDim.rmin2() / dd4hep::cm,
                     cryoDim.rmax1() / dd4hep::cm, (cryoDim.dz() - caloDim.dz()) / dd4hep::cm);

    dd4hep::Volume cryoFrontVol(cryostat.nameStr() + "_front", cryoFrontShape, aLcdd.material(cryostat.materialStr()));
    dd4hep::Volume cryoBackVol(cryostat.nameStr() + "_back", cryoBackShape, aLcdd.material(cryostat.materialStr()));
    dd4hep::Volume cryoSideVol(cryostat.nameStr() + "_side", cryoSideShape, aLcdd.material(cryostat.materialStr()));
    dd4hep::PlacedVolume cryoFrontPhysVol = envelopeVol.placeVolume(cryoFrontVol);
    dd4hep::PlacedVolume cryoBackPhysVol = envelopeVol.placeVolume(cryoBackVol);
    dd4hep::PlacedVolume cryoSidePhysVol = envelopeVol.placeVolume(cryoSideVol);

    if (cryoFrontSensitive) {
      cryoFrontVol.setSensitiveDetector(aSensDet);
      cryoFrontPhysVol.addPhysVolID("cryo", 1);
      cryoFrontPhysVol.addPhysVolID("type", 1);
      dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                       "ECAL Cryostat front volume set as sensitive");
    }

    if (cryoBackSensitive) {
      cryoBackVol.setSensitiveDetector(aSensDet);
      cryoBackPhysVol.addPhysVolID("cryo", 1);
      cryoBackPhysVol.addPhysVolID("type", 2);
      dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                       "ECAL Cryostat back volume set as sensitive");
    }

    if (cryoSideSensitive) {
      cryoSideVol.setSensitiveDetector(aSensDet);
      cryoSidePhysVol.addPhysVolID("cryo", 1);
      cryoSidePhysVol.addPhysVolID("type", 3);
      dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                       "ECAL Cryostat side volume set as sensitive");
    }

    dd4hep::DetElement cryoFrontDetElem(caloDetElem, "cryo_front", 0);
    cryoFrontDetElem.setPlacement(cryoFrontPhysVol);
    dd4hep::DetElement cryoBackDetElem(caloDetElem, "cryo_back", 0);
    cryoBackDetElem.setPlacement(cryoBackPhysVol);
    dd4hep::DetElement cryoSideDetElem(caloDetElem, "cryo_side", 0);
    cryoSideDetElem.setPlacement(cryoSidePhysVol);

    // 1.2. Create place-holder for services
    dd4hep::Tube servicesFrontShape(cryoDim.rmin2(), bathRmin, caloDim.dz());
    dd4hep::Tube servicesBackShape(bathRmax, cryoDim.rmax1(), caloDim.dz());

    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "ECAL services: front: rmin (cm) = %f  rmax (cm) = %f  dz (cm) = %f", cryoDim.rmin2() / dd4hep::cm,
                     bathRmin / dd4hep::cm / dd4hep::cm, caloDim.dz() / dd4hep::cm);
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "ECAL services: back: rmin (cm) = %f  rmax (cm) = %f  dz (cm) = %f", bathRmax / dd4hep::cm,
                     cryoDim.rmax1() / dd4hep::cm, caloDim.dz() / dd4hep::cm);

    dd4hep::Volume servicesFrontVol("services_front", servicesFrontShape, aLcdd.material(activeMaterial));
    servicesFrontVol.setVisAttributes(aLcdd, "service_bath");
    dd4hep::Volume servicesBackVol("services_back", servicesBackShape, aLcdd.material(activeMaterial));
    servicesBackVol.setVisAttributes(aLcdd, "service_bath");
    dd4hep::PlacedVolume servicesFrontPhysVol = envelopeVol.placeVolume(servicesFrontVol);
    dd4hep::PlacedVolume servicesBackPhysVol = envelopeVol.placeVolume(servicesBackVol);

    if (cryoFrontSensitive) {
      servicesFrontVol.setSensitiveDetector(aSensDet);
      servicesFrontPhysVol.addPhysVolID("cryo", 1);
      servicesFrontPhysVol.addPhysVolID("type", 4);
      dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                       "ECAL Services front volume set as sensitive");
    }

    if (cryoBackSensitive) {
      servicesBackVol.setSensitiveDetector(aSensDet);
      servicesBackPhysVol.addPhysVolID("cryo", 1);
      servicesBackPhysVol.addPhysVolID("type", 5);
      dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                       "ECAL Services back volume set as sensitive");
    }

    dd4hep::DetElement servicesFrontDetElem(caloDetElem, "services_front", 0);
    servicesFrontDetElem.setPlacement(servicesFrontPhysVol);
    dd4hep::DetElement servicesBackDetElem(caloDetElem, "services_back", 0);
    servicesBackDetElem.setPlacement(servicesBackPhysVol);
  }

  // 2. Create bath that is inside the cryostat and surrounds the detector
  //    Bath is filled with active material -> but not sensitive
  dd4hep::Volume bathVol(activeMaterial + "_bath", bathOuterShape, aLcdd.material(activeMaterial));
  double Rmin = caloDim.rmin();
  double Rmax = caloDim.rmax();
  double dR = Rmax - Rmin;

  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "ECAL bath: material = %s  rmin (cm) = %f  rmax (cm) = %f", activeMaterial.c_str(),
                   bathRmin / dd4hep::cm, bathRmax / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "ECAL calorimeter volume rmin (cm) = %f  rmax (cm) = %f", Rmin / dd4hep::cm, Rmax / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "ECAL thickness of calorimeter (cm) = %f", dR / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "ECAL thickness in front of calorimeter (between cryostat front and calorimeter) (cm) = %f",
                   (Rmin - cryoDim.rmin2()) / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "ECAL thickness behind calorimeter (between calorimeter and cryostat back) (cm) = %f",
                   (cryoDim.rmax1() - Rmax) / dd4hep::cm);

  // 3. Create the calorimeter by placing the passive material, trapezoid active layers, readout and again trapezoid
  // active layers in the bath.
  // sensitive detector for the layers (and, if desired to study energy deposited in absorbers or up/downstream, in
  // services and cryo)
  dd4hep::SensitiveDetector sd = aSensDet;
  dd4hep::xml::Dimension sdType = xmlDetElem.child(_U(sensitive));
  sd.setType(sdType.typeStr());

  // 3.a. Create the passive planes, readout in between of 2 passive planes and the remaining space filled with active
  // material

  // Start by first calculating geometry parameters and printing out info

  //////////////////////////////
  // PASSIVE PLANES
  //////////////////////////////

  // passive volumes consist of inner part and two outer, joined by glue
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "Passive elements:");
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   material in inner part of absorber (except 1st layer) = %s", passiveInnerMaterial.c_str());
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   material in inner part of absorber (1st layer) = %s", passiveInnerMaterialFirstLayer.c_str());
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   material in outer part of absorber = %s", passiveOuterMaterial.c_str());
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   material in middle part between inner and outer = %s", passiveGlueMaterial.c_str());
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   thickness of inner part at inner radius (cm) = %f", passiveInnerThicknessMin / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   thickness of inner part at outer radius (cm) = ", passiveInnerThicknessMax / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   thickness of outer part (cm) = %f", passiveOuterThickness / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   thickness of middle part (cm) = %f", passiveGlueThickness / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   total thickness of absorber at inner radius (cm) = %f", passiveThickness / dd4hep::cm);

  //////////////////////////////
  // ELECTRODES
  //////////////////////////////

  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "Electrodes:");
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   rotation angle (radians) = %f , (degrees) = %f", angle, angle * 57.295780);

  // calculate number of modules
  uint numPlanes =
      round(M_PI / asin((passiveThickness + activeThickness + readoutThickness) / (2. * Rmin * cos(angle))));

  double dPhi = 2. * M_PI / numPlanes;
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   number of planes (calculated) = %d , azim. angle difference = %f", numPlanes, dPhi);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   distance at inner radius (cm) = %f", 2. * M_PI * Rmin / dd4hep::cm / numPlanes);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   distance at outer radius (cm) = %f", 2. * M_PI * Rmax / dd4hep::cm / numPlanes);

  // The following code checks if the xml geometry file contains a constant defining
  // the number of planes in the barrel. In that case, it makes the program abort
  // if the number of planes in the xml is different from the one calculated from
  // the geometry. This is because the number of plane information (retrieved from the
  // xml) is used in other parts of the code (the readout for the FCC-ee ECAL with
  // inclined modules). In principle the code above should be refactored so that the number
  // of planes is one of the inputs of the calculation and other geometrical parameters
  // are adjusted accordingly. This is left for the future, and we use the workaround
  // below to enforce for the time being that the number of planes is "correct"
  int nModules = -1;
  try {
    nModules = aLcdd.constant<int>("ECalBarrelNumPlanes");
  } catch (...) {
    ;
  }
  if (nModules > 0 && nModules != int(numPlanes)) {
    dd4hep::printout(dd4hep::ERROR, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Incorrect number of planes (ECalBarrelNumPlanes) in calorimeter xml description!");
    throw std::runtime_error("Incorrect number of planes (ECalBarrelNumPlanes) in calorimeter xml description!");
  }

  // Readout is in the middle between two passive planes
  double offsetPassivePhi = caloDim.offset() + dPhi / 2.;
  double offsetReadoutPhi = caloDim.offset() + 0;
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "   readout material = %s",
                   readoutMaterial.c_str());
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   thickness of readout planes (cm) = %f", readoutThickness / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "   number of layers = %d",
                   numLayers);

  // Electrode length, given inclination angle and min/max radius of the active calorimeter volume
  double planeLength = -Rmin * cos(angle) + sqrt(pow(Rmax, 2) - pow(Rmin * sin(angle), 2));

  double runningHeight = 0.;
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   total length from Rmin, Rmax and angle (cm) = %f", planeLength / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "   predicted layer radii:");
  for (uint iLay = 0; iLay < numLayers; iLay++) {
    double Lmin = runningHeight;
    double Lmax = runningHeight + layerHeight[iLay];
    runningHeight += layerHeight[iLay];
    double rMin = sqrt(Rmin * Rmin + Lmin * Lmin + 2 * Rmin * Lmin * cos(angle));
    double rMax = sqrt(Rmin * Rmin + Lmax * Lmax + 2 * Rmin * Lmax * cos(angle));
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "      layer %d (cm) = %f - %f",
                     iLay, rMin / dd4hep::cm, rMax / dd4hep::cm);
  }

  // Check that electrode length is consistent with calorimeter radial extent
  // and inclination angle to within 0.5 mm
  if (fabs(planeLength - layersTotalHeight) > 0.05 * dd4hep::cm) {
    dd4hep::printout(
        dd4hep::ERROR, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
        "   the sum of the electrode lengths per layer in the calorimeter xml file is not consistent with the "
        "length calculated from the calorimeter radial extent and the inclination angle");
    throw std::runtime_error("Incorrect length of electrode layers in calorimeter xml description!");
  }

  // Calculate the thickness of the passive material in each layer
  // it's not constant in case of trapezoidal absorbers  (passiveInnerThicknessMax != passiveInnerThicknessMin)
  // the code calculates the (max) passive thickness per layer i.e. at Rout of layer
  // rescaling by runningHeight / Ltot
  std::vector<double> passiveInnerThicknessLayer(numLayers + 1);
  dd4hep::printout(dd4hep::DEBUG, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "passiveInnerThickness: ");
  runningHeight = 0.;
  for (uint iLay = 0; iLay < numLayers; iLay++) {
    passiveInnerThicknessLayer[iLay] =
        passiveInnerThicknessMin +
        (passiveInnerThicknessMax - passiveInnerThicknessMin) * (runningHeight) / layersTotalHeight;
    dd4hep::printout(dd4hep::DEBUG, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "   layer %d = %f cm", iLay,
                     passiveInnerThicknessLayer[iLay]);
  }
  passiveInnerThicknessLayer[numLayers] =
      passiveInnerThicknessMin +
      (passiveInnerThicknessMax - passiveInnerThicknessMin) * (runningHeight) / layersTotalHeight;
  dd4hep::printout(dd4hep::DEBUG, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "   layer %d = %f cm", numLayers,
                   passiveInnerThicknessLayer[numLayers]);

  // Calculate the angle of the passive elements if trapezoidal absorbers are used (Max != Min)
  // if parallel absorbers are used then passiveAngle = 0 and cosPassiveAngle = 1
  double passiveAngle = atan2((passiveInnerThicknessMax - passiveInnerThicknessMin) / 2., planeLength);
  double cosPassiveAngle = cos(passiveAngle);
  double rotatedOuterThickness = passiveOuterThickness / cosPassiveAngle;
  double rotatedGlueThickness = passiveGlueThickness / cosPassiveAngle;

  // Distance from the center of the electrode to the center of the 1st layer, in the electrode direction
  double layerFirstOffset = -planeLength / 2. + layerHeight[0] / 2.;

  //////////////////////////////
  // ACTIVE ELEMENTS
  //////////////////////////////

  // The active (LAr/LKr) elements are constructed subtracting from their
  // envelope the readout and passive elements.
  // This requires some calculations

  // Thickness of active layers at inner radius and outer ( = distance between passive plane and readout plane)
  // at inner radius: distance projected at plane perpendicular to readout plane
  double activeInThickness = Rmin * sin(dPhi / 2.) * cos(angle);
  activeInThickness -= passiveThickness * (0.5 - activePassiveOverlap);
  // at outer radius: distance projected at plane perpendicular to readout plane
  double activeOutThickness = (Rmin + planeLength) * sin(dPhi / 2.) * cos(angle);
  // make correction for outer readius caused by inclination angle
  // first calculate intersection of readout plane and plane parallel to shifted passive plane
  double xIntersect = (Rmin * (tan(angle) - cos(dPhi / 2.) * tan(angle + dPhi / 2.)) - planeLength * sin(dPhi / 2.)) /
                      (tan(angle) - tan(angle + dPhi / 2.));
  double yIntersect = tan(angle) * xIntersect + Rmin * (sin(dPhi / 2.) - tan(angle)) + planeLength * sin(dPhi / 2.);
  // distance from inner radius to intersection
  double correction =
      planeLength - sqrt(pow(xIntersect - Rmin * cos(dPhi / 2), 2) + pow(yIntersect - Rmin * sin(dPhi / 2), 2));
  // correction to the active thickness
  activeOutThickness += 2. * correction * sin(dPhi / 4.);
  activeOutThickness -= passiveThickness * (0.5 - activePassiveOverlap);
  // print the active layer dimensions
  double activeInThicknessAfterSubtraction =
      2. * activeInThickness - readoutThickness - 2. * activePassiveOverlap * passiveThickness;
  double activeOutThicknessAfterSubtraction =
      2. * activeOutThickness - readoutThickness -
      2. * activePassiveOverlap *
          (passiveThickness + passiveInnerThicknessMax - passiveInnerThicknessMin); // correct thickness for trapezoid
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "Active elements:");
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "   material = %s",
                   activeMaterial.c_str());
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   thickness at inner radius (cm) = %f", activeInThicknessAfterSubtraction / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   thickness at outer radius (cm) = %f", activeOutThicknessAfterSubtraction / dd4hep::cm);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   thickness relative increase from inner to outer radius = %f %%",
                   (activeOutThicknessAfterSubtraction - activeInThicknessAfterSubtraction) * 100 /
                       activeInThicknessAfterSubtraction);
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                   "   active passive initial overlap (before subtraction) (cm) = %f  (%f %%)",
                   passiveThickness / dd4hep::cm * activePassiveOverlap, activePassiveOverlap * 100);

  // Now create the volumes

  //////////////////////////////
  // PASSIVE ELEMENTS
  //////////////////////////////

  // in the first calo layer we might use a different material (LAr instead of Pb) for the inner passive material
  // to sample more uniformly for the upstream correction. So we need to split the volume into
  // passiveInnerShapeFirstLayer (for first layer) and passiveInnerShape (for other layers)
  // passiveInner elements' lengths along electrode are length of layer0 and suf of lengths of other layers
  // for the other passive elements, which are boxes, their length has to be corrected for the rotation
  // needed in case of trapezoidal absorbers, leading to the 1 / cosPassiveAngle factor
  // For the same reason, their thickness in the direction perpendicular to the electrode
  // has to be scaled by 1 / cosPassiveAngle, and thus we use rotatedOuterThickness and rotatedGlueThickness
  // here
  dd4hep::Trd1 passiveShape(passiveInnerThicknessMin / 2. + rotatedOuterThickness / 2. + rotatedGlueThickness / 2.,
                            passiveInnerThicknessMax / 2. + rotatedOuterThickness / 2. + rotatedGlueThickness / 2.,
                            caloDim.dz(), planeLength / 2.);
  dd4hep::Trd1 passiveInnerShapeFirstLayer(passiveInnerThicknessMin / 2., passiveInnerThicknessLayer[1] / 2.,
                                           caloDim.dz(), layerHeight[0] / 2.);
  dd4hep::Trd1 passiveInnerShape(passiveInnerThicknessLayer[1] / 2., passiveInnerThicknessMax / 2., caloDim.dz(),
                                 planeLength / 2. - layerHeight[0] / 2.);
  dd4hep::Box passiveOuterShape(passiveOuterThickness / 4., caloDim.dz(), planeLength / 2. / cosPassiveAngle);
  dd4hep::Box passiveGlueShape(passiveGlueThickness / 4., caloDim.dz(), planeLength / 2. / cosPassiveAngle);
  dd4hep::Volume passiveVol("passive", passiveShape, aLcdd.material("Air"));
  passiveVol.setVisAttributes(aLcdd, "absorbers");
  dd4hep::Volume passiveInnerVol(passiveInnerMaterial + "_passive", passiveInnerShape,
                                 aLcdd.material(passiveInnerMaterial));
  passiveInnerVol.setVisAttributes(aLcdd, "absorbers");
  dd4hep::Volume passiveInnerVolFirstLayer(passiveInnerMaterialFirstLayer + "_passive", passiveInnerShapeFirstLayer,
                                           aLcdd.material(passiveInnerMaterialFirstLayer));
  passiveInnerVolFirstLayer.setVisAttributes(aLcdd, "absorbers");
  dd4hep::Volume passiveOuterVol(passiveOuterMaterial + "_passive", passiveOuterShape,
                                 aLcdd.material(passiveOuterMaterial));
  passiveOuterVol.setVisAttributes(aLcdd, "absorbers");
  dd4hep::Volume passiveGlueVol(passiveGlueMaterial + "_passive", passiveGlueShape,
                                aLcdd.material(passiveGlueMaterial));
  passiveGlueVol.setVisAttributes(aLcdd, "absorbers");

  // translate and rotate the elements of the module appropriately
  // the Pb absorber does not need to be rotated, but the glue and
  // the outer absorbers have to, in case of trapezoidal absorbers.
  // The glue and steel absorbers also have to be translated in x-y
  // to the proper positions on the two sides of the inner absorber

  // - inner part of absorber, 2-N layers
  //   Their center is shifted wrt center of fulle absorber by length of 1st layer/2.0
  dd4hep::PlacedVolume passiveInnerPhysVol =
      passiveVol.placeVolume(passiveInnerVol, dd4hep::Position(0, 0, layerHeight[0] / 2.));
  // - inner part of absorber, first layer
  //   its center is shifted wrt center of full absorber by -(length of plane - length of 1st layer)/2.0
  dd4hep::PlacedVolume passiveInnerPhysVolFirstLayer =
      passiveVol.placeVolume(passiveInnerVolFirstLayer, dd4hep::Position(0, 0, layerFirstOffset));
  // - outer part of absorber, all layers, left side
  dd4hep::PlacedVolume passiveOuterPhysVolBelow = passiveVol.placeVolume(
      passiveOuterVol,
      dd4hep::Transform3D(dd4hep::RotationY(-passiveAngle),
                          dd4hep::Position(-(passiveInnerThicknessMin + passiveInnerThicknessMax) / 4. -
                                               rotatedGlueThickness / 2. - rotatedOuterThickness / 4.,
                                           0, 0)));
  // - outer part of absorber, all layers, right side
  dd4hep::PlacedVolume passiveOuterPhysVolAbove = passiveVol.placeVolume(
      passiveOuterVol, dd4hep::Transform3D(dd4hep::RotationY(passiveAngle),
                                           dd4hep::Position((passiveInnerThicknessMin + passiveInnerThicknessMax) / 4. +
                                                                rotatedGlueThickness / 2. + rotatedOuterThickness / 4.,
                                                            0, 0)));
  // - glue in absorber, all layers, left side
  dd4hep::PlacedVolume passiveGluePhysVolBelow = passiveVol.placeVolume(
      passiveGlueVol, dd4hep::Transform3D(dd4hep::RotationY(-passiveAngle),
                                          dd4hep::Position(-(passiveInnerThicknessMin + passiveInnerThicknessMax) / 4. -
                                                               rotatedGlueThickness / 4.,
                                                           0, 0)));
  // - glue in absorber, all layers, right side
  dd4hep::PlacedVolume passiveGluePhysVolAbove = passiveVol.placeVolume(
      passiveGlueVol, dd4hep::Transform3D(dd4hep::RotationY(passiveAngle),
                                          dd4hep::Position((passiveInnerThicknessMin + passiveInnerThicknessMax) / 4. +
                                                               rotatedGlueThickness / 4.,
                                                           0, 0)));
  passiveInnerPhysVol.addPhysVolID("subtype", 0);
  passiveInnerPhysVolFirstLayer.addPhysVolID("subtype", 0);
  passiveOuterPhysVolBelow.addPhysVolID("subtype", 1);
  passiveOuterPhysVolAbove.addPhysVolID("subtype", 2);
  passiveGluePhysVolBelow.addPhysVolID("subtype", 3);
  passiveGluePhysVolAbove.addPhysVolID("subtype", 4);

  // if the inner part of the absorber is sensitive (to study energy deposited in it, for calculation of per-layer
  // sampling fraction), then it is divided in layer volumes, and each layer volume is set as sensitive
  // first layer
  if (passiveInner.isSensitive()) {
    dd4hep::printout(dd4hep::DEBUG, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Passive inner volume (1st layer) set as sensitive");
    passiveInnerVolFirstLayer.setSensitiveDetector(aSensDet);
    passiveInnerPhysVolFirstLayer.addPhysVolID("layer", 0);
    dd4hep::DetElement passiveInnerDetElemFirstLayer("layer", 0);
    passiveInnerDetElemFirstLayer.setPlacement(passiveInnerPhysVolFirstLayer);
  }
  // other layers
  if (passiveInnerMax.isSensitive()) {
    dd4hep::printout(dd4hep::DEBUG, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Passive inner volume (2-N layers) set as sensitive");
    double layerOffset = layerFirstOffset + layerHeight[1] / 2.;
    for (uint iLayer = 1; iLayer < numLayers; iLayer++) {
      dd4hep::Trd1 layerPassiveInnerShape(passiveInnerThicknessLayer[iLayer] / 2.,
                                          passiveInnerThicknessLayer[iLayer + 1] / 2., caloDim.dz(),
                                          layerHeight[iLayer] / 2.);
      dd4hep::Volume layerPassiveInnerVol(passiveInnerMaterial, layerPassiveInnerShape,
                                          aLcdd.material(passiveInnerMaterial));
      layerPassiveInnerVol.setSensitiveDetector(aSensDet);
      dd4hep::PlacedVolume layerPassiveInnerPhysVol =
          passiveInnerVol.placeVolume(layerPassiveInnerVol, dd4hep::Position(0, 0, layerOffset));
      layerPassiveInnerPhysVol.addPhysVolID("layer", iLayer);
      dd4hep::DetElement layerPassiveInnerDetElem("layer", iLayer);
      layerPassiveInnerDetElem.setPlacement(layerPassiveInnerPhysVol);
      if (iLayer != numLayers - 1) {
        layerOffset += layerHeight[iLayer] / 2. + layerHeight[iLayer + 1] / 2.;
      }
    }
  }

  if (passiveOuter.isSensitive()) {
    // if the outer part of the absorber is sensitive (to study energy deposited in it, for calculation of per-layer
    // sampling fraction), then it is divided in layer volumes, and each layer volume is set as sensitive
    dd4hep::printout(dd4hep::DEBUG, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Passive outer volume set as sensitive");
    double layerOffset = layerFirstOffset / cosPassiveAngle;
    for (uint iLayer = 0; iLayer < numLayers; iLayer++) {
      dd4hep::Box layerPassiveOuterShape(passiveOuterThickness / 4., caloDim.dz(),
                                         layerHeight[iLayer] / 2. / cosPassiveAngle);
      dd4hep::Volume layerPassiveOuterVol(passiveOuterMaterial, layerPassiveOuterShape,
                                          aLcdd.material(passiveOuterMaterial));
      layerPassiveOuterVol.setSensitiveDetector(aSensDet);
      dd4hep::PlacedVolume layerPassiveOuterPhysVol =
          passiveOuterVol.placeVolume(layerPassiveOuterVol, dd4hep::Position(0, 0, layerOffset));
      layerPassiveOuterPhysVol.addPhysVolID("layer", iLayer);
      dd4hep::DetElement layerPassiveOuterDetElem("layer", iLayer);
      layerPassiveOuterDetElem.setPlacement(layerPassiveOuterPhysVol);
      if (iLayer != numLayers - 1) {
        layerOffset += (layerHeight[iLayer] / 2. + layerHeight[iLayer + 1] / 2.) / cosPassiveAngle;
      }
    }
  }

  if (passiveGlue.isSensitive()) {
    // if the glue is sensitive (to study energy deposited in it, for calculation of per-layer
    // sampling fraction), then it is divided in layer volumes, and each layer volume is set as sensitive
    dd4hep::printout(dd4hep::DEBUG, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Passive glue volume set as sensitive");
    double layerOffset = layerFirstOffset / cosPassiveAngle;
    for (uint iLayer = 0; iLayer < numLayers; iLayer++) {
      dd4hep::Box layerPassiveGlueShape(passiveGlueThickness / 4., caloDim.dz(),
                                        layerHeight[iLayer] / 2. / cosPassiveAngle);
      dd4hep::Volume layerPassiveGlueVol(passiveGlueMaterial, layerPassiveGlueShape,
                                         aLcdd.material(passiveGlueMaterial));
      layerPassiveGlueVol.setSensitiveDetector(aSensDet);
      dd4hep::PlacedVolume layerPassiveGluePhysVol =
          passiveGlueVol.placeVolume(layerPassiveGlueVol, dd4hep::Position(0, 0, layerOffset));
      layerPassiveGluePhysVol.addPhysVolID("layer", iLayer);
      dd4hep::DetElement layerPassiveGlueDetElem("layer", iLayer);
      layerPassiveGlueDetElem.setPlacement(layerPassiveGluePhysVol);
      if (iLayer != numLayers - 1) {
        layerOffset += (layerHeight[iLayer] / 2. + layerHeight[iLayer + 1] / 2.) / cosPassiveAngle;
      }
    }
  }

  //////////////////////////////
  // READOUT PLANES
  //////////////////////////////
  dd4hep::Box readoutShape(readoutThickness / 2., caloDim.dz(), planeLength / 2.);
  dd4hep::Volume readoutVol(readoutMaterial, readoutShape, aLcdd.material(readoutMaterial));
  readoutVol.setVisAttributes(aLcdd, "pcb");
  // if the readout is sensitive (to study energy deposited in it, for calculation of per-layer
  // sampling fraction), then it is divided in layer volumes, and each layer volume is set as sensitive
  if (readout.isSensitive()) {
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "ECAL readout volume set as sensitive");
    double layerOffset = layerFirstOffset;
    for (uint iLayer = 0; iLayer < numLayers; iLayer++) {
      dd4hep::Box layerReadoutShape(readoutThickness / 2., caloDim.dz(), layerHeight[iLayer] / 2.);
      dd4hep::Volume layerReadoutVol(readoutMaterial, layerReadoutShape, aLcdd.material(readoutMaterial));
      layerReadoutVol.setSensitiveDetector(aSensDet);
      dd4hep::PlacedVolume layerReadoutPhysVol =
          readoutVol.placeVolume(layerReadoutVol, dd4hep::Position(0, 0, layerOffset));
      layerReadoutPhysVol.addPhysVolID("layer", iLayer);
      dd4hep::DetElement layerReadoutDetElem("layer", iLayer);
      layerReadoutDetElem.setPlacement(layerReadoutPhysVol);
      if (iLayer != numLayers - 1) {
        layerOffset += layerHeight[iLayer] / 2. + layerHeight[iLayer + 1] / 2.;
      }
    }
  }

  //////////////////////////////
  // ACTIVE ELEMENTS
  //////////////////////////////

  // creating shape for rows of layers (active material between two passive planes, with readout in the middle)

  // - first define area between two passive planes, area can reach up to the symmetry axis of passive plane
  dd4hep::Trd1 activeOuterShape(activeInThickness, activeOutThickness, caloDim.dz(), planeLength / 2.);

  // - subtract readout shape from the middle
  dd4hep::SubtractionSolid activeShapeNoReadout(activeOuterShape, readoutShape);

  // - make calculation for active plane that is inclined with 0 deg (= offset + angle)
  double Cx = Rmin * cos(-angle) + planeLength / 2.;
  double Cy = Rmin * sin(-angle);
  double Ax = Rmin * cos(-angle + dPhi / 2.) + planeLength / 2. * cos(dPhi / 2.);
  double Ay = Rmin * sin(-angle + dPhi / 2.) + planeLength / 2. * sin(dPhi / 2.);
  double CAx = fabs(Ax - Cx);
  double CAy = fabs(Ay - Cy);
  double zprim, xprim;
  zprim = CAx;
  xprim = CAy;

  double Bx = Rmin * cos(-angle - dPhi / 2.) + planeLength / 2. * cos(-dPhi / 2.);
  double By = Rmin * sin(-angle - dPhi / 2.) + planeLength / 2. * sin(-dPhi / 2.);
  double CBx = fabs(Bx - Cx);
  double CBy = fabs(By - Cy);
  double zprimB, xprimB;
  zprimB = CBx;
  xprimB = CBy;

  // - subtract readout shape from the middle
  dd4hep::SubtractionSolid activeShapeNoPassiveAbove(
      activeShapeNoReadout, passiveShape,
      dd4hep::Transform3D(dd4hep::RotationY(-dPhi / 2.), dd4hep::Position(-fabs(xprim), 0, fabs(zprim))));

  // - subtract passive volume below
  dd4hep::SubtractionSolid activeShape(
      activeShapeNoPassiveAbove, passiveShape,
      dd4hep::Transform3D(dd4hep::RotationY(dPhi / 2.), dd4hep::Position(fabs(xprimB), 0, -fabs(zprimB))));

  // - create the active volume, which will contain the layers filled with LAr
  dd4hep::Volume activeVol("active", activeShape, aLcdd.material("Air"));
  activeVol.setVisAttributes(aLcdd, "sensitive_bath");

  // place layers within active volume
  std::vector<dd4hep::PlacedVolume> layerPhysVols;

  // - first, calculate the layer widths at inner and outer radii
  std::vector<double> layerInThickness;
  std::vector<double> layerOutThickness;
  double layerIncreasePerUnitThickness = (activeOutThickness - activeInThickness) / layersTotalHeight;
  for (uint iLay = 0; iLay < numLayers; iLay++) {
    if (iLay == 0) {
      layerInThickness.push_back(activeInThickness);
    } else {
      layerInThickness.push_back(layerOutThickness[iLay - 1]);
    }
    layerOutThickness.push_back(layerInThickness[iLay] + layerIncreasePerUnitThickness * layerHeight[iLay]);
  }

  // - then, loop on the layers to create and place the volumes
  double layerOffset = layerFirstOffset;
  for (uint iLayer = 0; iLayer < numLayers; iLayer++) {
    // define the layer envelope
    dd4hep::Trd1 layerOuterShape(layerInThickness[iLayer], layerOutThickness[iLayer], caloDim.dz(),
                                 layerHeight[iLayer] / 2.);

    // subtract readout shape from the middle
    dd4hep::SubtractionSolid layerShapeNoReadout(layerOuterShape, readoutShape);

    // subtract readout shape from the middle
    dd4hep::SubtractionSolid layerShapeNoPassiveAbove(
        layerShapeNoReadout, passiveShape,
        dd4hep::Transform3D(dd4hep::RotationY(-dPhi / 2.),
                            dd4hep::Position(-fabs(xprim), 0, fabs(zprim) - layerOffset)));

    // subtract passive volume below
    dd4hep::SubtractionSolid layerShape(
        layerShapeNoPassiveAbove, passiveShape,
        dd4hep::Transform3D(dd4hep::RotationY(dPhi / 2.),
                            dd4hep::Position(fabs(xprimB), 0, -fabs(zprimB) - layerOffset)));

    // create the volume, filled with active material, set as sensitive, and position it properly within active volume
    dd4hep::Volume layerVol("layer", layerShape, aLcdd.material(activeMaterial));
    layerVol.setVisAttributes(aLcdd, "sensitive_bath");
    layerVol.setSensitiveDetector(aSensDet);
    layerPhysVols.push_back(activeVol.placeVolume(layerVol, dd4hep::Position(0, 0, layerOffset)));
    layerPhysVols.back().addPhysVolID("layer", iLayer);
    if (iLayer != numLayers - 1) {
      layerOffset += layerHeight[iLayer] / 2. + layerHeight[iLayer + 1] / 2.;
    }
  }

  // Place elements in bath: passive planes, readout planes and rows of layers
  dd4hep::DetElement bathDetElem(caloDetElem, "bath", 1);
  for (uint iPlane = 0; iPlane < numPlanes; iPlane++) {

    //
    // PASSIVE
    //

    // calculate centre position of the plane without plane rotation
    double phi = offsetPassivePhi + iPlane * dPhi;
    double xRadial = (Rmin + planeLength / 2.) * cos(phi);
    double yRadial = (Rmin + planeLength / 2.) * sin(phi);
    // calculate position of the beginning of plane
    double xRmin = Rmin * cos(phi);
    double yRmin = Rmin * sin(phi);
    // rotate centre by angle wrt beginning of plane
    double xRotated = xRmin + (xRadial - xRmin) * cos(angle) - (yRadial - yRmin) * sin(angle);
    double yRotated = yRmin + (xRadial - xRmin) * sin(angle) + (yRadial - yRmin) * cos(angle);
    dd4hep::Transform3D transform(dd4hep::RotationX(-M_PI / 2.)     // to get in XY plane
                                      * dd4hep::RotationY(M_PI / 2. // to get pointed towards centre
                                                          - phi - angle),
                                  dd4hep::Position(xRotated, yRotated, 0));
    // create and place volume in bath
    dd4hep::PlacedVolume passivePhysVol = bathVol.placeVolume(passiveVol, transform);
    passivePhysVol.addPhysVolID("module", iPlane);
    passivePhysVol.addPhysVolID("type", 1); // 0 = active, 1 = passive, 2 = readout
    dd4hep::DetElement passiveDetElem(bathDetElem, "passive" + std::to_string(iPlane), iPlane);
    passiveDetElem.setPlacement(passivePhysVol);

    //
    // READOUT
    //

    // calculate centre position of the plane without plane rotation
    double phiRead = offsetReadoutPhi + iPlane * dPhi;
    double xRadialRead = (Rmin + planeLength / 2.) * cos(phiRead);
    double yRadialRead = (Rmin + planeLength / 2.) * sin(phiRead);
    // calculate position of the beginning of plane
    double xRminRead = Rmin * cos(phiRead);
    double yRminRead = Rmin * sin(phiRead);
    // rotate centre by angle wrt beginning of plane
    double xRotatedRead = xRminRead + (xRadialRead - xRminRead) * cos(angle) - (yRadialRead - yRminRead) * sin(angle);
    double yRotatedRead = yRminRead + (xRadialRead - xRminRead) * sin(angle) + (yRadialRead - yRminRead) * cos(angle);
    dd4hep::Transform3D transformRead(dd4hep::RotationX(-M_PI / 2.)     // to get in XY plane
                                          * dd4hep::RotationY(M_PI / 2. // to get pointed towards centre
                                                              - phiRead - angle),
                                      dd4hep::Position(xRotatedRead, yRotatedRead, 0));
    // create and place volume in bath
    dd4hep::PlacedVolume readoutPhysVol = bathVol.placeVolume(readoutVol, transformRead);
    readoutPhysVol.addPhysVolID("module", iPlane);
    readoutPhysVol.addPhysVolID("type", 2); // 0 = active, 1 = passive, 2 = readout
    dd4hep::DetElement readoutDetElem(bathDetElem, "readout" + std::to_string(iPlane), iPlane);
    readoutDetElem.setPlacement(readoutPhysVol);

    //
    // ACTIVE
    //

    // same positioning as readout
    dd4hep::Transform3D transformActive(dd4hep::RotationX(-M_PI / 2) * dd4hep::RotationY(M_PI / 2 - phiRead - angle),
                                        dd4hep::Position(xRotatedRead, yRotatedRead, 0));
    // create and place volume in bath
    dd4hep::PlacedVolume activePhysVol = bathVol.placeVolume(activeVol, transformActive);
    activePhysVol.addPhysVolID("module", iPlane);
    activePhysVol.addPhysVolID("type", 0); // 0 = active, 1 = passive, 2 = readout
    dd4hep::DetElement activeDetElem(bathDetElem, "active" + std::to_string(iPlane), iPlane);
    activeDetElem.setPlacement(activePhysVol);
    // place the layers inside the active element
    for (uint iLayer = 0; iLayer < numLayers; iLayer++) {
      dd4hep::DetElement layerDetElem(activeDetElem, "layer" + std::to_string(iLayer), iLayer);
      layerDetElem.setPlacement(layerPhysVols[iLayer]);
    }
  }

  // Place bath in envelope
  dd4hep::PlacedVolume bathPhysVol = envelopeVol.placeVolume(bathVol);
  bathDetElem.setPlacement(bathPhysVol);

  // Place the envelope
  dd4hep::Volume motherVol = aLcdd.pickMotherVolume(caloDetElem);
  dd4hep::PlacedVolume envelopePhysVol = motherVol.placeVolume(envelopeVol);
  envelopePhysVol.addPhysVolID("system", xmlDetElem.id());
  caloDetElem.setPlacement(envelopePhysVol);

  // Create caloData object for the reconstruction
  auto caloData = new dd4hep::rec::LayeredCalorimeterData;
  caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;
  caloDetElem.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);
  dd4hep::xml::setDetectorTypeFlag(xmlDetElem, caloDetElem);

  // Fill caloData information

  // Extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in dd4hep units (mm)
  // (for barrel detectors zmin is 0)
  caloData->extent[0] = Rmin;
  caloData->extent[1] = Rmax;
  caloData->extent[2] = 0.;
  caloData->extent[3] = caloDim.dz();

  // Retrieve segmentation, needed to get cell size in theta and phi
  dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo* seg =
      dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(
          aSensDet.readout().segmentation().segmentation());
  if (seg == nullptr) {
    dd4hep::printout(dd4hep::ERROR, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "Incorrect readout, cannot cast to FCCSWGridModuleThetaMerged");
    throw std::runtime_error("Incorrect readout in calorimeter xml description!");
  }
  std::string layerFieldName = seg->fieldNameLayer();
  std::string cellIDEncoding = aSensDet.readout().idSpec().fieldDescription();
  dd4hep::BitFieldCoder encoder(cellIDEncoding);

  // Information about each layer
  // double distance : distance from Origin (or the z-axis) to the inner-most face of the layer
  // double phi0 : phi0 of layer: potential rotation around normal to absorber plane, e.g. if layers are 'staggered' in
  // phi in fwd. calos double absorberThickness : thickness of the absorber part of the layer. Consider using
  // inner/outer_nRadiationLengths and inner/outer_nInteractionLengths double inner_nRadiationLengths : Absorber
  // material in front of sensitive element in the layer, units of radiation lengths double inner_nInteractionLengths :
  // Absorber material in front of sensitive element in the layer, units of radiation lengths double
  // outer_nRadiationLengths : Absorber material in behind of sensitive element in the layer, units of radiation lengths
  // double outer_nInteractionLengths : Absorber material in behind of sensitive element in the layer, units of
  // radiation lengths double inner_thickness : Distance between the innermost face of the layer (closest to IP) and the
  // center of the sensitive element double outer_thickness : Distance between the center of the sensitive element and
  // the outermost face of the layer double sensitive_thickness : Thickness of the sensitive element (e.g. scintillator)
  // double cellSize0 : cell size along the first axis where first is either along the beam (BarrelLayout) or up
  // (EndcapLayout) or the direction closest to that double cellSize1 : second cell size, perpendicular to the first
  // direction cellSize0 and the depth of the layers
  dd4hep::rec::MaterialManager matMgr(envelopeVol);
  dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;
  dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "Layer structure information:");
  runningHeight = 0.0;
  for (uint iLay = 0; iLay < numLayers; iLay++) {
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "  Layer %d", iLay);
    double Lmin = runningHeight;
    double Lmax = runningHeight + layerHeight[iLay];
    double Lmid = runningHeight + layerHeight[iLay] / 2.0;
    runningHeight += layerHeight[iLay];
    double rad_first = sqrt(Rmin * Rmin + Lmin * Lmin + 2 * Rmin * Lmin * cos(angle));
    double rad_last = sqrt(Rmin * Rmin + Lmax * Lmax + 2 * Rmin * Lmax * cos(angle));
    double rad_mid = sqrt(Rmin * Rmin + Lmid * Lmid + 2 * Rmin * Lmid * cos(angle));

    double thickness_sen = 0.;
    double absorberThickness = 0.;

    dd4hep::rec::Vector3D ivr1 =
        dd4hep::rec::Vector3D(0., rad_first, 0); // defining starting vector points of the given layer
    dd4hep::rec::Vector3D ivr2 =
        dd4hep::rec::Vector3D(0., rad_last, 0); // defining end vector points of the given layer

    dd4hep::printout(dd4hep::DEBUG, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "    radius first = %f, radius last = %f, radius middle = %f", rad_first, rad_last, rad_mid);
    const dd4hep::rec::MaterialVec& materials =
        matMgr.materialsBetween(ivr1, ivr2); // calling material manager to get material info between two points
    auto mat = matMgr.createAveragedMaterial(materials); // creating average of all the material between two points to
                                                         // calculate X0 and lambda of averaged material
    const double nRadiationLengths = mat.radiationLength();
    const double nInteractionLengths = mat.interactionLength();
    const double difference_bet_r1r2 = (ivr1 - ivr2).r();
    const double value_of_x0 = difference_bet_r1r2 / nRadiationLengths;
    const double value_of_lambda = difference_bet_r1r2 / nInteractionLengths;
    std::string str1("LAr");

    for (size_t imat = 0; imat < materials.size(); imat++) {

      std::string str2(materials.at(imat).first.name());
      if (str1.compare(str2) == 0) {
        thickness_sen += materials.at(imat).second;
      } else {
        absorberThickness += materials.at(imat).second;
      }
    }
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "    sensitive thickness : %f",
                     thickness_sen);
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03", "    absorber thickness : %f",
                     absorberThickness);
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "    radiation length : %f , interaction length : %f", value_of_x0, value_of_lambda);

    caloLayer.distance = rad_first;
    caloLayer.absorberThickness = absorberThickness;
    caloLayer.sensitive_thickness = thickness_sen;

    caloLayer.inner_nRadiationLengths = value_of_x0 * (rad_mid - rad_first) / (rad_last - rad_first);
    caloLayer.inner_nInteractionLengths = value_of_lambda * (rad_mid - rad_first) / (rad_last - rad_first);
    caloLayer.inner_thickness = rad_mid - rad_first;

    caloLayer.outer_nRadiationLengths = value_of_x0 * (rad_last - rad_mid) / (rad_last - rad_first);
    caloLayer.outer_nInteractionLengths = value_of_lambda * (rad_last - rad_mid) / (rad_last - rad_first);
    caloLayer.outer_thickness = rad_last - rad_mid;

    // retrieve cell dimensions vector from segmentation class
    // set volume ID and layer ID
    // cell size does not depend on theta/module so no need to set moduleID and thetaID
    dd4hep::CellID cID;
    encoder.set(cID, layerFieldName, iLay);
    std::vector<double> cellSizeVector = seg->cellDimensions(cID);
    double cellSizeTheta = cellSizeVector[1];
    double cellSizeModule = cellSizeVector[0];
    double cellSizePhi = dPhi * cellSizeModule;
    dd4hep::printout(dd4hep::INFO, "ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03",
                     "    cell sizes in theta, phi: %f , %f", cellSizeTheta, cellSizePhi);
    caloLayer.cellSize0 = cellSizeTheta;
    caloLayer.cellSize1 = cellSizePhi;

    caloData->layers.push_back(caloLayer);
  }
  return caloDetElem;
}
} // namespace det

DECLARE_DETELEMENT(ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03, det::createECalBarrelInclined)

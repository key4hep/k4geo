#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMatrixT.h"
#include "XML/Utilities.h"
#include <DDRec/DetectorData.h>

namespace det {

namespace ECalEndcap_Turbine_o1_v03 {
  unsigned ECalEndCapElementCounter = 0;
  const unsigned nWheels = 3;

  unsigned ECalEndcapNumCalibRhoLayersArr[nWheels], ECalEndcapNumCalibZLayersArr[nWheels];

  double tForArcLength(double s, double bladeangle, double delZ, double r) {

    // some intermediate constants
    double zpos = delZ / 2.;
    double zp = zpos / TMath::Tan(bladeangle);
    double b = zp / (TMath::Sqrt(r * r - zp * zp));
    double c = (TMath::Tan(s / r) + b) / (1. - b * TMath::Tan(s / r));
    double d = c * c * r * r / (1 + c * c);
    return (TMath::Sqrt(d) - zp) * TMath::Sin(bladeangle);
  }

  // return position of the inner edge of a blade
  double getZmin(double r, double bladeangle, double delZ) {
    // r: distance from the beamline
    // bladeangle: angle of turbine blades wrt xy plane, in radians
    // delZ: z extent of the blades
    return TMath::Sqrt(r * r - ((delZ / 2) / TMath::Tan(bladeangle)) * ((delZ / 2) / TMath::Tan(bladeangle)));
  }

  dd4hep::Solid buildOneBlade(double thickness_inner, double thickness_outer, double width, double ro, double ri,
                              double bladeangle, double delZ, double zStart) {

    dd4hep::Solid shapeBeforeSubtraction;

    // set max and min extent of the blade (along the z axis in the body frame)
    double zmax = ro;
    double zmin = getZmin(ri, bladeangle, delZ);

    dd4hep::Trd2 tmp1(thickness_inner / 2., thickness_outer / 2., width / 2., width / 2., (zmax - zmin) / 2.);
    shapeBeforeSubtraction = tmp1;

    dd4hep::Tube allowedTube(ri, ro, delZ / 2.);

    return dd4hep::IntersectionSolid(
        shapeBeforeSubtraction, allowedTube,
        dd4hep::Transform3D(dd4hep::RotationZYX(0, TMath::Pi() / 2. - bladeangle, TMath::Pi() / 2.),
                            dd4hep::Position(0, -zStart, -(zmin + zmax) / 2.)));
  }

  void buildWheel(dd4hep::Detector& aLcdd, dd4hep::SensitiveDetector& aSensDet, dd4hep::Volume& aEnvelope,
                  dd4hep::xml::Handle_t& aXmlElement, dd4hep::DetElement& bathDetElem, float ri, float ro, float delZ,
                  float offsetZ, unsigned iWheel) {

    dd4hep::xml::DetElement calorimeterElem = aXmlElement.child(_Unicode(calorimeter));
    dd4hep::xml::DetElement genericBladeElem = calorimeterElem.child(_Unicode(turbineBlade));
    dd4hep::xml::DetElement absBladeElem = genericBladeElem.child(_Unicode(absorberBlade));
    dd4hep::xml::DetElement claddingElem = genericBladeElem.child(_Unicode(cladding));
    dd4hep::xml::DetElement glueElem = genericBladeElem.child(_Unicode(glue));
    dd4hep::xml::DetElement electrodeBladeElem = genericBladeElem.child(_Unicode(electrodeBlade));
    dd4hep::xml::DetElement nobleLiquidElem = genericBladeElem.child(_Unicode(nobleLiquidGap));

    float BladeAngle = 0.0, AbsThickMin = 0.0, BladeThicknessScaleFactor = 0.0;
    unsigned nUnitCells = 0;
    // hardcode for three wheels
    unsigned ECalEndcapNumCalibRhoLayers = ECalEndcapNumCalibRhoLayersArr[iWheel],
             ECalEndcapNumCalibZLayers = ECalEndcapNumCalibZLayersArr[iWheel];

    unsigned LayerIndexBaseline = 0;

    if (iWheel == 0) {
      BladeAngle = genericBladeElem.attr<float>(_Unicode(angle1));
      AbsThickMin = absBladeElem.attr<float>(_Unicode(thickness1));
      BladeThicknessScaleFactor = absBladeElem.attr<float>(_Unicode(thicknessScaleFactor1));

      nUnitCells = genericBladeElem.attr<int>(_Unicode(nUnitCells1));
    }
    if (iWheel == 1) {
      BladeAngle = genericBladeElem.attr<float>(_Unicode(angle2));
      AbsThickMin = absBladeElem.attr<float>(_Unicode(thickness2));
      BladeThicknessScaleFactor = absBladeElem.attr<float>(_Unicode(thicknessScaleFactor2));
      nUnitCells = genericBladeElem.attr<int>(_Unicode(nUnitCells2));
      LayerIndexBaseline = ECalEndcapNumCalibRhoLayersArr[0] * ECalEndcapNumCalibZLayersArr[0];
    }
    if (iWheel == 2) {
      BladeAngle = genericBladeElem.attr<float>(_Unicode(angle3));
      AbsThickMin = absBladeElem.attr<float>(_Unicode(thickness3));
      BladeThicknessScaleFactor = absBladeElem.attr<float>(_Unicode(thicknessScaleFactor3));
      nUnitCells = genericBladeElem.attr<int>(_Unicode(nUnitCells3));
      LayerIndexBaseline = ECalEndcapNumCalibRhoLayersArr[0] * ECalEndcapNumCalibZLayersArr[0] +
                           ECalEndcapNumCalibRhoLayersArr[1] * ECalEndcapNumCalibZLayersArr[1];
    }

    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "Making wheel with inner, outer radii %f, %f", ri, ro);
    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "Blade angle is %f ", BladeAngle);
    dd4hep::xml::Dimension dim(aXmlElement.child(_Unicode(dimensions)));
    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "delZ is %f", delZ);

    if (TMath::Abs(TMath::Tan(BladeAngle)) < delZ / (2. * ri)) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "The requested blade angle is too small for the given delZ and ri values.  Please adjust to at "
                       "least %f degrees!",
                       TMath::ATan(delZ / (2. * ri)) * 180. / TMath::Pi());
      return;
    }

    Float_t xRange = delZ / (TMath::Sin(BladeAngle));

    double delrPhiNoGap;

    float GlueThick = glueElem.attr<float>(_Unicode(thickness));
    float CladdingThick = claddingElem.attr<float>(_Unicode(thickness));
    AbsThickMin = AbsThickMin - (GlueThick + CladdingThick);
    if (AbsThickMin < 0.) {
      dd4hep::printout(
          dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
          "Error: requested absorber thickness is negative after accounting for glue and cladding thickness");
    }
    float ElectrodeThick = electrodeBladeElem.attr<float>(_Unicode(thickness));
    float LArgapi = nobleLiquidElem.attr<float>(_Unicode(gap));

    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "nUnitCells: %d", nUnitCells);

    float AbsThicki = AbsThickMin;
    // make volumes for the noble liquid, electrode, and absorber blades
    float AbsThicko;
    AbsThicko = AbsThicki + BladeThicknessScaleFactor * ((ro / ri) - 1.) * AbsThicki;
    // Calculate gap thickness at inner layer
    double circ = 2 * TMath::Pi() * ri;
    double x2 = (AbsThickMin + (GlueThick + CladdingThick) + ElectrodeThick) / TMath::Sin(BladeAngle);
    double y2 = TMath::Sqrt(ri * ri - x2 * x2);
    double rPhi1 = ri * TMath::Pi() / 2.;
    double rPhi2 = ri * TMath::ATan(y2 / x2);
    delrPhiNoGap = TMath::Abs(rPhi1 - rPhi2);
    double leftoverS = (circ - nUnitCells * delrPhiNoGap);
    double delrPhiGapOnly = leftoverS / (2 * nUnitCells);
    LArgapi = delrPhiGapOnly * TMath::Sin(BladeAngle);
    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "LArGap at inner radius is %f", LArgapi);

    // now find gap at outer radius
    circ = 2 * TMath::Pi() * ro;
    x2 = (AbsThicko + GlueThick + CladdingThick + ElectrodeThick) / TMath::Sin(BladeAngle);
    y2 = TMath::Sqrt(ro * ro - x2 * x2);
    rPhi1 = ro * TMath::Pi() / 2.;
    rPhi2 = ro * TMath::ATan(y2 / x2);
    delrPhiNoGap = TMath::Abs(rPhi1 - rPhi2);
    leftoverS = (circ - nUnitCells * delrPhiNoGap);
    delrPhiGapOnly = leftoverS / (2 * nUnitCells);
    float LArgapo = delrPhiGapOnly * TMath::Sin(BladeAngle);

    float riLayer = ri;

    std::vector<dd4hep::Volume> claddingLayerVols;
    std::vector<dd4hep::Volume> glueLayerVols;
    std::vector<dd4hep::Volume> absBladeLayerVols;
    std::vector<dd4hep::Volume> LArTotalLayerVols;
    std::vector<dd4hep::Volume> electrodeBladeLayerVols;

    dd4hep::Solid passiveShape =
        buildOneBlade(AbsThicki + GlueThick + CladdingThick, AbsThicko + GlueThick + CladdingThick, xRange, ro, ri,
                      BladeAngle, delZ, 0);
    dd4hep::Volume passiveVol("passive", passiveShape, aLcdd.material("Air"));

    dd4hep::Solid activeShape =
        buildOneBlade(ElectrodeThick + LArgapi * 2, ElectrodeThick + LArgapo * 2, xRange, ro, ri, BladeAngle, delZ, 0);
    dd4hep::Volume activeVol("active", activeShape, aLcdd.material("Air"));

    unsigned numNonActiveRhoLayers = 1;
    unsigned numNonActiveZLayers = 1;

    // check that either all non-active volumes are set to sensitive (for
    // sampling fraction calculations) or none are (for normal running)
    bool allNonActiveSensitive = (claddingElem.isSensitive() && glueElem.isSensitive() && absBladeElem.isSensitive() &&
                                  electrodeBladeElem.isSensitive());
    bool allNonActiveNotSensitive = (!claddingElem.isSensitive() && !glueElem.isSensitive() &&
                                     !absBladeElem.isSensitive() && !electrodeBladeElem.isSensitive());
    if (allNonActiveSensitive) {
      numNonActiveRhoLayers = ECalEndcapNumCalibRhoLayers;
      numNonActiveZLayers = ECalEndcapNumCalibZLayers;
    } else if (!allNonActiveNotSensitive) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "Some non-active layers are sensitive and others are not -- this is likely a misconfiguration");
    }

    float delrNonActive = (ro - ri) / numNonActiveRhoLayers;
    float delrActive = (ro - ri) / ECalEndcapNumCalibRhoLayers;

    for (unsigned iRhoLayer = 0; iRhoLayer < numNonActiveRhoLayers; iRhoLayer++) {
      float roLayer = riLayer + delrNonActive;
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Making layer with inner, outer radii %f, %f",
                       riLayer, roLayer);

      AbsThicko = AbsThicki + BladeThicknessScaleFactor * ((roLayer / riLayer) - 1.) * AbsThicki;

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "Inner and outer absorber thicknesses %f, %f ",
                       AbsThicki, AbsThicko);

      float zStart = -xRange / 2. + xRange / (2. * numNonActiveZLayers);

      for (unsigned iZLayer = 0; iZLayer < numNonActiveZLayers; iZLayer++) {

        dd4hep::Solid claddingLayer =
            buildOneBlade(AbsThicki + GlueThick + CladdingThick, AbsThicko + GlueThick + CladdingThick,
                          xRange / numNonActiveZLayers, roLayer, riLayer, BladeAngle, delZ, zStart);

        dd4hep::Solid glueLayer =
            buildOneBlade(AbsThicki + GlueThick, AbsThicko + GlueThick, xRange / numNonActiveZLayers, roLayer, riLayer,
                          BladeAngle, delZ, zStart);

        dd4hep::Solid absBladeLayer = buildOneBlade(AbsThicki, AbsThicko, xRange / numNonActiveZLayers, roLayer,
                                                    riLayer, BladeAngle, delZ, zStart);

        dd4hep::Volume claddingLayerVol("claddingLayer", claddingLayer, aLcdd.material(claddingElem.materialStr()));
        if (claddingElem.isSensitive()) {
          claddingLayerVol.setSensitiveDetector(aSensDet);
        }
        claddingLayerVols.push_back(claddingLayerVol);

        dd4hep::Volume glueLayerVol("glueLayer", glueLayer, aLcdd.material(glueElem.materialStr()));
        if (glueElem.isSensitive()) {
          glueLayerVol.setSensitiveDetector(aSensDet);
        }
        glueLayerVols.push_back(glueLayerVol);

        dd4hep::Volume absBladeLayerVol("absBladeLayer", absBladeLayer, aLcdd.material(absBladeElem.materialStr()));
        if (absBladeElem.isSensitive()) {
          absBladeLayerVol.setSensitiveDetector(aSensDet);
        }
        absBladeLayerVols.push_back(absBladeLayerVol);

        zStart += xRange / numNonActiveZLayers;
      }
      riLayer = roLayer;
      AbsThicki = AbsThicko;
    }

    riLayer = ri;

    AbsThicki = AbsThickMin;

    float LArgapiLayer = LArgapi;

    for (unsigned iRhoLayer = 0; iRhoLayer < ECalEndcapNumCalibRhoLayers; iRhoLayer++) {

      float roLayer = riLayer + delrActive;

      AbsThicko = AbsThicki + BladeThicknessScaleFactor * ((roLayer / riLayer) - 1.) * AbsThicki;

      /*
      // now find gap at outer layer
      circ = 2 * TMath::Pi() * roLayer;
      x2 = (AbsThicko + GlueThick + CladdingThick + ElectrodeThick) / TMath::Sin(BladeAngle);
      y2 = TMath::Sqrt(roLayer * roLayer - x2 * x2);
      rPhi1 = roLayer * TMath::Pi() / 2.;
      rPhi2 = roLayer * TMath::ATan(y2 / x2);
      delrPhiNoGap = TMath::Abs(rPhi1 - rPhi2);
      leftoverS = (circ - nUnitCells * delrPhiNoGap);
      delrPhiGapOnly = leftoverS / (2 * nUnitCells);
      LArgapo = delrPhiGapOnly * TMath::Sin(BladeAngle);
      */
      float LArgapoLayer = LArgapi + (LArgapo - LArgapi) * (1.0 * iRhoLayer) / ECalEndcapNumCalibRhoLayers;

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "Outer LAr gap is %f", LArgapoLayer);
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03",
                       "Inner and outer thicknesses of noble liquid volume %f, %f", ElectrodeThick + LArgapiLayer * 2,
                       ElectrodeThick + LArgapoLayer * 2);

      float zStart = -xRange / 2. + xRange / (2. * ECalEndcapNumCalibZLayers);

      for (unsigned iZLayer = 0; iZLayer < ECalEndcapNumCalibZLayers; iZLayer++) {

        dd4hep::Solid electrodeBladeAndGapLayer =
            buildOneBlade(ElectrodeThick + LArgapiLayer * 2, ElectrodeThick + LArgapoLayer * 2,
                          xRange / ECalEndcapNumCalibZLayers, roLayer, riLayer, BladeAngle, delZ, zStart);

        dd4hep::Solid electrodeBladeLayer =
            buildOneBlade(ElectrodeThick, ElectrodeThick, xRange / ECalEndcapNumCalibZLayers, roLayer, riLayer,
                          BladeAngle, delZ, zStart);

        dd4hep::Volume electrodeBladeLayerVol("electrodeBladeLayer", electrodeBladeLayer,
                                              aLcdd.material(electrodeBladeElem.materialStr()));
        if (electrodeBladeElem.isSensitive()) {
          electrodeBladeLayerVol.setSensitiveDetector(aSensDet);
        }
        electrodeBladeLayerVols.push_back(electrodeBladeLayerVol);

        dd4hep::Volume LArTotalLayerVol("LArTotalLayerVol", electrodeBladeAndGapLayer,
                                        aLcdd.material(nobleLiquidElem.materialStr()));

        if (nobleLiquidElem.isSensitive()) {
          LArTotalLayerVol.setSensitiveDetector(aSensDet);
        }
        LArTotalLayerVols.push_back(LArTotalLayerVol);

        zStart += xRange / ECalEndcapNumCalibZLayers;
      }
      riLayer = roLayer;
      LArgapiLayer = LArgapoLayer;
      AbsThicki = AbsThicko;
    }
    dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03",
                     "ECal endcap materials:  nobleLiquid: %s absorber %s electrode %s",
                     nobleLiquidElem.materialStr().c_str(), absBladeElem.materialStr().c_str(),
                     electrodeBladeElem.materialStr().c_str());

    int nUnitCellsToDraw = nUnitCells;

    dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Number of unit cells %d", nUnitCells);

    // place all components of the absorber blade inside passive volume

    unsigned iLayer = 0;

    riLayer = ri;

    double xOffset = -xRange / 2 + xRange / (2 * numNonActiveZLayers);

    for (auto absBladeLayerVol : absBladeLayerVols) {

      float roLayer = riLayer + delrNonActive;

      dd4hep::Position posLayer(0, 0, 0);
      xOffset += xRange / numNonActiveZLayers;
      dd4hep::PlacedVolume absBladeVol_pv = glueLayerVols[iLayer].placeVolume(absBladeLayerVol, posLayer);

      absBladeVol_pv.addPhysVolID("subtype", 1); // 1 = absorber, 2 = glue, 3 = cladding
      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03_geo", "Blade layer, rho is %d, %f, %f", iLayer,
                       absBladeVol_pv.position().Rho(), roLayer / 2.);

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03_geo", "AbsBalde volume %s",
                       absBladeVol_pv.toString().c_str());

      iLayer++;
      if (iLayer % numNonActiveZLayers == 0) {
        riLayer = roLayer;
        xOffset = -xRange / 2 + xRange / (2 * numNonActiveZLayers);
      }
    }

    riLayer = ri;
    iLayer = 0;

    xOffset = -xRange / 2 + xRange / (2 * numNonActiveZLayers);

    for (auto glueLayerVol : glueLayerVols) {

      float roLayer = riLayer + delrNonActive;

      dd4hep::Position posLayer(0, 0, 0);
      dd4hep::PlacedVolume glueVol_pv = claddingLayerVols[iLayer].placeVolume(glueLayerVol, posLayer);
      xOffset += xRange / numNonActiveZLayers;

      glueVol_pv.addPhysVolID("subtype", 2); // 1 = absorber, 2 = glue, 3 = cladding

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03_geo", "Glue volume %s", glueVol_pv.toString().c_str());

      iLayer++;
      if (iLayer % numNonActiveZLayers == 0) {
        riLayer = roLayer;
        xOffset = -xRange / 2 + xRange / (2 * numNonActiveZLayers);
      }
    }

    riLayer = ri;
    iLayer = 0;

    double zminri = getZmin(ri, BladeAngle, delZ);

    xOffset = -xRange / 2 + xRange / (2 * numNonActiveZLayers);

    for (auto claddingLayerVol : claddingLayerVols) {

      float roLayer = riLayer + delrNonActive;

      double zminLayer = getZmin(riLayer, BladeAngle, delZ);

      dd4hep::Position posLayer(0, xOffset, (zminLayer - zminri + roLayer - ro) / 2.);
      xOffset += xRange / numNonActiveZLayers;
      dd4hep::PlacedVolume claddingVol_pv = passiveVol.placeVolume(claddingLayerVol, posLayer);

      claddingVol_pv.addPhysVolID("layer", LayerIndexBaseline + iLayer);

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03_geo", "Cladding volume %s",
                       claddingVol_pv.toString().c_str());

      iLayer++;
      if (iLayer % numNonActiveZLayers == 0) {
        riLayer = roLayer;
        xOffset = -xRange / 2 + xRange / (2 * numNonActiveZLayers);
      }
    }

    riLayer = ri;
    iLayer = 0;

    xOffset = -xRange / 2 + xRange / (2 * ECalEndcapNumCalibZLayers);

    for (auto electrodeBladeLayerVol : electrodeBladeLayerVols) {

      float roLayer = riLayer + delrActive;
      dd4hep::Position posLayer(0, 0, 0);

      dd4hep::PlacedVolume electrodeBladeVol_pv =
          LArTotalLayerVols[iLayer].placeVolume(electrodeBladeLayerVol, posLayer);
      xOffset += xRange / ECalEndcapNumCalibZLayers;
      electrodeBladeVol_pv.addPhysVolID("type", 2); // 0 = active, 1 = passive, 2 = readout

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03_geo", "Electrode volume %s",
                       electrodeBladeVol_pv.toString().c_str());

      iLayer++;
      if (iLayer % ECalEndcapNumCalibZLayers == 0) {
        riLayer = roLayer;
        xOffset = -xRange / 2 + xRange / (2 * ECalEndcapNumCalibZLayers);
      }
    }

    riLayer = ri;
    iLayer = 0;

    std::vector<dd4hep::PlacedVolume> LArVol_pvs;

    xOffset = -xRange / 2 + xRange / (2 * ECalEndcapNumCalibZLayers);

    for (auto LArTotalLayerVol : LArTotalLayerVols) {

      float roLayer = riLayer + delrActive;

      double zminLayer = getZmin(riLayer, BladeAngle, delZ);

      dd4hep::Position posLayer(0, xOffset, (zminLayer - zminri + roLayer - ro) / 2.);

      dd4hep::PlacedVolume LArVol_pv(activeVol.placeVolume(LArTotalLayerVol, posLayer));
      xOffset += xRange / ECalEndcapNumCalibZLayers;

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "LAr layer: %d layer in readout: %d", iLayer,
                       LayerIndexBaseline + iLayer);
      LArVol_pv.addPhysVolID("layer", LayerIndexBaseline + iLayer);
      LArVol_pvs.push_back(LArVol_pv);

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03_geo", "LAr volume %s", LArVol_pv.toString().c_str());

      iLayer++;
      if (iLayer % ECalEndcapNumCalibZLayers == 0) {
        riLayer = roLayer;
        xOffset = -xRange / 2 + xRange / (2 * ECalEndcapNumCalibZLayers);
      }
    }

    for (int iUnitCell = 0; iUnitCell < nUnitCellsToDraw; iUnitCell++) {

      int modIndex = iUnitCell - nUnitCellsToDraw / 2;
      if (modIndex < 0)
        modIndex += nUnitCells;
      float phi = (iUnitCell - nUnitCellsToDraw / 2) * 2 * TMath::Pi() / nUnitCells;
      float delPhi = 2 * TMath::Pi() / nUnitCells;

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "Placing blade, ro, ri = %f %f", ro, ri);

      TGeoRotation tgr;
      tgr.RotateZ(BladeAngle * 180 / TMath::Pi());
      tgr.RotateX(-phi * 180 / TMath::Pi());
      tgr.RotateY(90);

      const Double_t* rotMatPtr;

      rotMatPtr = tgr.GetRotationMatrix();
      TMatrixT<Double_t> rotMat(3, 3, rotMatPtr);
      dd4hep::Rotation3D r3d;
      r3d.SetComponents(rotMat(0, 0), rotMat(0, 1), rotMat(0, 2), rotMat(1, 0), rotMat(1, 1), rotMat(1, 2),
                        rotMat(2, 0), rotMat(2, 1), rotMat(2, 2));

      tgr.Clear();
      tgr.RotateZ(BladeAngle * 180 / TMath::Pi());
      tgr.RotateX(-(phi + delPhi / 2.) * 180 / TMath::Pi());
      tgr.RotateY(90);

      rotMatPtr = tgr.GetRotationMatrix();
      TMatrixT<Double_t> rotMat2(3, 3, rotMatPtr);
      dd4hep::Rotation3D r3d2;
      r3d2.SetComponents(rotMat2(0, 0), rotMat2(0, 1), rotMat2(0, 2), rotMat2(1, 0), rotMat2(1, 1), rotMat2(1, 2),
                         rotMat2(2, 0), rotMat2(2, 1), rotMat2(2, 2));

      riLayer = ri;

      float xCell = ((ro + zminri) / 2.) * TMath::Cos(phi);
      float yCell = ((ro + zminri) / 2.) * TMath::Sin(phi); // ri*TMath::Sin(phi)/6.;
      float zCell = offsetZ;

      dd4hep::Transform3D comCell(r3d, dd4hep::Translation3D(xCell, yCell, zCell));

      // place passive volume in LAr bath
      dd4hep::PlacedVolume passivePhysVol = aEnvelope.placeVolume(passiveVol, comCell);
      passivePhysVol.addPhysVolID("module", modIndex);
      passivePhysVol.addPhysVolID("wheel", iWheel);
      passivePhysVol.addPhysVolID("type", 1); // 0 = active, 1 = passive, 2 = readout
      dd4hep::DetElement passiveDetElem("passive_" + std::to_string(iUnitCell) + "_" + std::to_string(iWheel),
                                        modIndex);
      passiveDetElem.setPlacement(passivePhysVol);
      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03_geo", "Passive volume %s",
                       passivePhysVol.toString().c_str());

      // place active volume in LAr bath

      xCell = ((ro + zminri) / 2.) * TMath::Cos(phi + delPhi / 2.);
      yCell = ((ro + zminri) / 2.) * TMath::Sin(phi + delPhi / 2.); // ri*TMath::Sin(phi)/6.;
      zCell = offsetZ;
      dd4hep::Transform3D comCell2(r3d2, dd4hep::Translation3D(xCell, yCell, zCell));
      dd4hep::PlacedVolume activePhysVol = aEnvelope.placeVolume(activeVol, comCell2);
      activePhysVol.addPhysVolID("module", modIndex);
      activePhysVol.addPhysVolID("wheel", iWheel);
      activePhysVol.addPhysVolID("type", 0); // 0 = active, 1 = passive, 2 = readout

      dd4hep::DetElement activeDetElem(bathDetElem, "active" + std::to_string(iUnitCell) + "_" + std::to_string(iWheel),
                                       modIndex);

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03_geo", "Active volume %s",
                       activePhysVol.toString().c_str());

      activeDetElem.setPlacement(activePhysVol);
      iLayer = 0;
      for (auto LArVol_pv : LArVol_pvs) {
        dd4hep::DetElement LArDetElem(
            activeDetElem,
            "layer" + std::to_string(modIndex) + "_" + std::to_string(iWheel) + "_" + std::to_string(iLayer), iLayer);
        LArDetElem.setPlacement(LArVol_pv);
        iLayer++;
      }

      dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "LArTotalLayerVols.size = %d",
                       LArTotalLayerVols.size());
    }

    return;
  }

  void buildMechanicalSupport(dd4hep::Detector& aLcdd, dd4hep::Volume& aEnvelope, dd4hep::xml::Handle_t& aXmlElement,
                              dd4hep::DetElement& bathDetElem) {

    dd4hep::xml::DetElement calo = aXmlElement.child(_Unicode(calorimeter));
    dd4hep::xml::DetElement mechSupport = calo.child(_Unicode(mechSupport));

    float mechSupportZCenter = aLcdd.constant<float>("ECalEndcapSupportZCenter");

    // make inner ring
    dd4hep::xml::DetElement innerRing = mechSupport.child(_Unicode(innerRing));

    dd4hep::xml::Dimension innerRingDim(innerRing.dimensions());

    double innerRingRmin = innerRingDim.rmin1();
    double innerRingRmax = innerRingDim.rmax1();
    double innerRingdZ = innerRingDim.dz();

    dd4hep::Tube innerRingTube(innerRingRmin, innerRingRmax, innerRingdZ);
    dd4hep::Volume innerRingVol("mechSupportInnerRing", innerRingTube, aLcdd.material(mechSupport.materialStr()));
    dd4hep::PlacedVolume innerRing_pv = aEnvelope.placeVolume(innerRingVol, dd4hep::Position(0, 0, mechSupportZCenter));
    dd4hep::DetElement innerRingDetElem(bathDetElem, "mechSupportInnerRing", 0);
    innerRingDetElem.setPlacement(innerRing_pv);

    // make outer ring
    dd4hep::xml::DetElement outerRing = mechSupport.child(_Unicode(outerRing));

    dd4hep::xml::Dimension outerRingDim(outerRing.dimensions());

    double outerRingRmin = outerRingDim.rmin1();
    double outerRingRmax = outerRingDim.rmax1();
    double outerRingdZ = outerRingDim.dz();

    dd4hep::Tube outerRingTube(outerRingRmin, outerRingRmax, outerRingdZ);
    dd4hep::Volume outerRingVol("mechSupportOuterRing", outerRingTube, aLcdd.material(mechSupport.materialStr()));
    dd4hep::PlacedVolume outerRing_pv = aEnvelope.placeVolume(outerRingVol, dd4hep::Position(0, 0, mechSupportZCenter));
    dd4hep::DetElement outerRingDetElem(bathDetElem, "mechSupportOuterRing", 0);
    outerRingDetElem.setPlacement(outerRing_pv);

    // make intermediate rings
    dd4hep::xml::DetElement supportTubeElem = calo.child(_Unicode(supportTube));
    float supportTubeThickness = supportTubeElem.thickness();
    unsigned nWheelsXML = supportTubeElem.attr<unsigned>(_Unicode(nWheels));
    dd4hep::xml::DetElement cryostat = calo.child(_Unicode(cryostat));
    dd4hep::xml::Dimension cryoDim(cryostat.dimensions());

    float intermedRingdR = aLcdd.constant<float>("ECalEndcapSupportMidRingdR");

    double rmin = cryoDim.rmin2();
    double rmax = cryoDim.rmax1();
    float radiusRatio = pow(rmax / rmin, 1. / nWheelsXML);
    double ro = rmin * radiusRatio;

    for (unsigned iWheel = 0; iWheel < nWheelsXML - 1; iWheel++) {
      double rInner = ro + (supportTubeThickness - intermedRingdR) / 2.;
      double rOuter = rInner + intermedRingdR;
      dd4hep::Tube intermedRingTube(rInner, rOuter, innerRingDim.dz());
      dd4hep::Volume intermedRingVol("mechSupportIntermedRing", intermedRingTube,
                                     aLcdd.material(mechSupport.materialStr()));
      dd4hep::PlacedVolume intermedRing_pv =
          aEnvelope.placeVolume(intermedRingVol, dd4hep::Position(0, 0, mechSupportZCenter));
      dd4hep::DetElement intermedRingDetElem(bathDetElem, "mechSupportInterMedRing" + std::to_string(iWheel), iWheel);
      intermedRingDetElem.setPlacement(intermedRing_pv);
      ro = ro * radiusRatio;
    }
    // make spokes
    unsigned nSpokes = aLcdd.constant<unsigned>("ECalEndcapSupportNSpokes");
    double rminSpoke = innerRingRmax;
    double rmaxSpoke = rmin * radiusRatio + (supportTubeThickness - intermedRingdR) / 2.;
    double spokeWidth = aLcdd.constant<float>("ECalEndcapSupportSpokeWidth");
    for (unsigned iWheel = 0; iWheel < nWheelsXML; iWheel++) {
      if (iWheel == nWheelsXML - 1)
        rmaxSpoke = outerRingRmin;

      dd4hep::Tube allowedSpokeRegion(rminSpoke, rmaxSpoke, innerRingDim.dz());
      dd4hep::Box spokeBox(spokeWidth / 2., rmaxSpoke - rminSpoke / 2., innerRingDim.dz());
      dd4hep::IntersectionSolid spokeShape(spokeBox, allowedSpokeRegion,
                                           dd4hep::Position(0, -(rminSpoke + rmaxSpoke) / 2., 0));
      dd4hep::Volume spokeVol("mechSupportSpoke", spokeShape, aLcdd.material(mechSupport.materialStr()));

      double rmid = (rminSpoke + rmaxSpoke) / 2.;
      for (unsigned iSpoke = 0; iSpoke < nSpokes; iSpoke++) {
        double phiSpoke = iSpoke * 2. * TMath::Pi() / nSpokes;

        dd4hep::PlacedVolume spoke_pv = aEnvelope.placeVolume(
            spokeVol, dd4hep::Transform3D(dd4hep::RotationZYX(-phiSpoke, 0, 0),
                                          dd4hep::Translation3D(rmid * TMath::Sin(phiSpoke),
                                                                rmid * TMath::Cos(phiSpoke), mechSupportZCenter)));
        dd4hep::DetElement spokeDetElem(bathDetElem,
                                        "mechSupportSpoke" + std::to_string(iWheel) + std::to_string(iSpoke), 0);
        spokeDetElem.setPlacement(spoke_pv);
      }
      rminSpoke = rmaxSpoke + intermedRingdR;
      if (iWheel == 0) {
        rmaxSpoke = rmin * radiusRatio * radiusRatio + (supportTubeThickness - intermedRingdR) / 2.;
      } else {
        rmaxSpoke = rminSpoke * radiusRatio + (supportTubeThickness - intermedRingdR) / 2.;
      }
    }
  }

  void buildOneSide_Turbine(dd4hep::Detector& aLcdd, dd4hep::DetElement& caloDetElem,
                            dd4hep::SensitiveDetector& aSensDet, dd4hep::Volume& aEnvelope,
                            dd4hep::xml::Handle_t& aXmlElement, unsigned& iModule) {

    dd4hep::xml::DetElement calo = aXmlElement.child(_Unicode(calorimeter));
    dd4hep::xml::Dimension caloDim(calo.dimensions());

    dd4hep::xml::DetElement blade = calo.child(_Unicode(turbineBlade));
    dd4hep::xml::DetElement nobleLiquid = blade.child(_Unicode(nobleLiquidGap));

    dd4hep::xml::DetElement xmlDetElem = aXmlElement;
    std::string nameDet = xmlDetElem.nameStr();

    dd4hep::xml::Dimension dim(aXmlElement.child(_Unicode(dimensions)));

    // build cryostat
    //  Retrieve cryostat data
    dd4hep::xml::DetElement cryostat = calo.child(_Unicode(cryostat));
    dd4hep::xml::Dimension cryoDim(cryostat.dimensions());
    double cryoThicknessFront = aLcdd.constant<float>("CryoEMECThicknessFront");
    double cryoThicknessBack = aLcdd.constant<float>("CryoEMECThicknessBack");
    double bathThicknessFront = aLcdd.constant<float>("BathThicknessFront");
    double bathThicknessBack = aLcdd.constant<float>("BathThicknessBack");

    dd4hep::xml::DetElement cryoFront = cryostat.child(_Unicode(front));
    dd4hep::xml::DetElement cryoBack = cryostat.child(_Unicode(back));
    dd4hep::xml::DetElement cryoInner = cryostat.child(_Unicode(inner));
    dd4hep::xml::DetElement cryoOuter = cryostat.child(_Unicode(outer));

    bool cryoFrontSensitive = cryoFront.isSensitive();
    bool cryoBackSensitive = cryoBack.isSensitive();
    bool cryoInnerSensitive = cryoInner.isSensitive();
    bool cryoOuterSensitive = cryoOuter.isSensitive();

    double bathRmin = cryoDim.rmin2(); // - margin for inclination
    double bathRmax = cryoDim.rmax1(); // + margin for inclination
    double bathDelZ = cryoDim.dz();
    dd4hep::Tube bathOuterShape(bathRmin, bathRmax, bathDelZ); // make it 4 volumes + 5th for detector envelope

    dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Cryostat front thickness is %f", cryoDim.rmin2());

    float NLcenterZ = (cryoThicknessFront - cryoThicknessBack) / 2.;
    if (cryoThicknessFront > 0) {
      // 1. Create cryostat
      dd4hep::Tube cryoFrontShape(cryoDim.rmin1(), cryoDim.rmax2(), cryoThicknessFront / 2.);
      dd4hep::Tube cryoBackShape(cryoDim.rmin1(), cryoDim.rmax2(), cryoThicknessBack / 2.);
      dd4hep::Tube cryoInnerShape(cryoDim.rmin1(), cryoDim.rmin2(), cryoDim.dz());
      dd4hep::Tube cryoOuterShape(cryoDim.rmax1(), cryoDim.rmax2(), cryoDim.dz());
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03",
                       "ECAL endcap cryostat: front: rmin (cm) = %f rmax (cm) = %f dz (cm) = %f ", cryoDim.rmin1(),
                       cryoDim.rmin2(), cryoDim.dz());
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03",
                       "ECAL encdap cryostat: back: rmin (cm) =  %f rmax (cm) = %f dz (cm) = %f", cryoDim.rmax1(),
                       cryoDim.rmax2(), cryoDim.dz());
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03",
                       "ECAL endcap cryostat: side: rmin (cm) =  %f rmax (cm) = %f dz (cm) = %f", cryoDim.rmin2(),
                       cryoDim.rmax1(), cryoDim.dz() - caloDim.dz());
      dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Cryostat is made out of %s",
                       cryostat.materialStr().c_str());

      dd4hep::Volume cryoFrontVol(cryostat.nameStr() + "_front", cryoFrontShape,
                                  aLcdd.material(cryostat.materialStr()));
      dd4hep::Volume cryoBackVol(cryostat.nameStr() + "_back", cryoBackShape, aLcdd.material(cryostat.materialStr()));
      dd4hep::Volume cryoInnerVol(cryostat.nameStr() + "_inner", cryoInnerShape,
                                  aLcdd.material(cryostat.materialStr()));
      dd4hep::Volume cryoOuterVol(cryostat.nameStr() + "_outer", cryoOuterShape,
                                  aLcdd.material(cryostat.materialStr()));

      dd4hep::Position cryoFrontPos(0, 0, NLcenterZ - cryoDim.dz() - cryoThicknessFront / 2.);
      dd4hep::PlacedVolume cryoFrontPhysVol = aEnvelope.placeVolume(cryoFrontVol, cryoFrontPos);
      dd4hep::Position cryoBackPos(0, 0, NLcenterZ + cryoDim.dz() + cryoThicknessBack / 2.);
      dd4hep::PlacedVolume cryoBackPhysVol = aEnvelope.placeVolume(cryoBackVol, cryoBackPos);
      dd4hep::Position cryoInnerOuterPos(0, 0, NLcenterZ);
      dd4hep::PlacedVolume cryoInnerPhysVol = aEnvelope.placeVolume(cryoInnerVol, cryoInnerOuterPos);
      dd4hep::PlacedVolume cryoOuterPhysVol = aEnvelope.placeVolume(cryoOuterVol, cryoInnerOuterPos);
      unsigned sidetype = 0x4; // probably not needed anymore...
      if (cryoFrontSensitive) {
        cryoFrontVol.setSensitiveDetector(aSensDet);
        cryoFrontPhysVol.addPhysVolID("cryo", 1);
        cryoFrontPhysVol.addPhysVolID("type", sidetype + 1);
        dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Cryostat front volume set as sensitive");
      }
      if (cryoBackSensitive) {
        cryoBackVol.setSensitiveDetector(aSensDet);
        cryoBackPhysVol.addPhysVolID("cryo", 1);
        cryoBackPhysVol.addPhysVolID("type", sidetype + 2);
        dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Cryostat back volume set as sensitive");
      }
      if (cryoInnerSensitive) {
        cryoInnerVol.setSensitiveDetector(aSensDet);
        cryoInnerPhysVol.addPhysVolID("cryo", 1);
        cryoInnerPhysVol.addPhysVolID("type", sidetype + 3);
        dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Cryostat inner volume set as sensitive");
      }
      if (cryoOuterSensitive) {
        cryoOuterVol.setSensitiveDetector(aSensDet);
        cryoOuterPhysVol.addPhysVolID("cryo", 1);
        cryoOuterPhysVol.addPhysVolID("type", sidetype + 4);
        dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Cryostat outer volume set as sensitive");
      }
      dd4hep::DetElement cryoFrontDetElem(caloDetElem, "cryo_front", 0);
      cryoFrontDetElem.setPlacement(cryoFrontPhysVol);
      dd4hep::DetElement cryoBackDetElem(caloDetElem, "cryo_back", 0);
      cryoBackDetElem.setPlacement(cryoBackPhysVol);
      dd4hep::DetElement cryoInnerDetElem(caloDetElem, "cryo_inner", 0);
      cryoInnerDetElem.setPlacement(cryoInnerPhysVol);
      dd4hep::DetElement cryoOuterDetElem(caloDetElem, "cryo_outer", 0);
      cryoOuterDetElem.setPlacement(cryoOuterPhysVol);
    }

    // 2. Create noble liquid bath
    std::string nobleLiquidMaterial = nobleLiquid.materialStr();
    dd4hep::Volume bathVol(nobleLiquidMaterial + "_bath", bathOuterShape, aLcdd.material(nobleLiquidMaterial));
    dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03",
                     "ECAL endcap bath: material = %s rmin (cm) = %f rmax (cm) = %f, dz (cm) = %f, thickness in front "
                     "of ECal (cm) = %f,  thickness behind ECal (cm) = %f",
                     nobleLiquidMaterial.c_str(), bathRmin, bathRmax, caloDim.dz(), caloDim.rmin() - cryoDim.rmin2(),
                     cryoDim.rmax1() - caloDim.rmax());

    //    dd4hep::Position bathPos(0, 0, (cryoThicknessFront - cryoThicknessBack)/2.);
    dd4hep::Position bathPos(0, 0, NLcenterZ);

    dd4hep::PlacedVolume bathPhysVol = aEnvelope.placeVolume(bathVol, bathPos);

    dd4hep::DetElement bathDetElem(caloDetElem, "bath", 1);

    bathDetElem.setPlacement(bathPhysVol);

    // 3. Create detector structure
    double length = dim.dz() * 2.;
    double zOffsetEnvelope = -length / 2.;

    dd4hep::xml::DetElement supportTubeElem = calo.child(_Unicode(supportTube));
    unsigned nWheelsXML = supportTubeElem.attr<unsigned>(_Unicode(nWheels));
    if (nWheelsXML != nWheels) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "Number of wheels in XML (%d) does not match hard-coded number of wheels (%d) ", nWheelsXML,
                       nWheels);
    }

    dd4hep::printout(dd4hep::INFO, "ECalEndcap_Turbine_o1_v03", "Will build %d wheels", nWheels);

    float supportTubeThickness = supportTubeElem.thickness();
    float mechSupportZGap = aLcdd.constant<float>("ECalEndcapSupportZGap");

    double rmin = bathRmin;
    double rmax = bathRmax - supportTubeThickness;
    float radiusRatio = pow(rmax / rmin, 1. / nWheels);
    double ro = rmin * radiusRatio;
    double ri = rmin;

    for (unsigned iWheel = 0; iWheel < nWheels; iWheel++) {

      dd4hep::Tube supportTube(ro, ro + supportTubeThickness, bathDelZ - (bathThicknessBack - mechSupportZGap) / 2.);

      dd4hep::Volume supportTubeVol("supportTube", supportTube, aLcdd.material(supportTubeElem.materialStr()));
      if (supportTubeElem.isSensitive()) {
        supportTubeVol.setSensitiveDetector(aSensDet);
      }
      dd4hep::PlacedVolume supportTube_pv = bathVol.placeVolume(
          supportTubeVol,
          dd4hep::Position(0, 0, zOffsetEnvelope + dim.dz() - (bathThicknessBack - mechSupportZGap) / 2.));
      supportTube_pv.addPhysVolID("cryo", 1);
      supportTube_pv.addPhysVolID("wheel", iWheel);
      dd4hep::DetElement supportTubeDetElem(bathDetElem, "supportTube_" + std::to_string(iWheel), 0);
      supportTubeDetElem.setPlacement(supportTube_pv);

      buildWheel(aLcdd, aSensDet, bathVol, aXmlElement, bathDetElem, ri + supportTubeThickness, ro,
                 bathDelZ * 2 - bathThicknessFront - bathThicknessBack, (bathThicknessFront - bathThicknessBack) / 2.,
                 iWheel);
      ri = ro;
      ro *= radiusRatio;
      if (ro > rmax)
        ro = rmax;
    }

    buildMechanicalSupport(aLcdd, bathVol, aXmlElement, bathDetElem);

    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03", "Total number of modules:  %d", iModule);

    return;
  }

  static dd4hep::Ref_t createECalEndcapTurbine(dd4hep::Detector& aLcdd, dd4hep::xml::Handle_t aXmlElement,
                                               dd4hep::SensitiveDetector aSensDet) {

    dd4hep::xml::DetElement xmlDetElem = aXmlElement;
    std::string nameDet = xmlDetElem.nameStr();
    int idDet = xmlDetElem.id();
    dd4hep::xml::Dimension dim(xmlDetElem.dimensions());
    dd4hep::DetElement caloDetElem(nameDet, idDet);
    dd4hep::xml::Dimension sdType = xmlDetElem.child(_U(sensitive));
    aSensDet.setType(sdType.typeStr());

    unsigned numReadoutRhoLayers, numReadoutZLayers;
    ECalEndcapNumCalibRhoLayersArr[0] = aLcdd.constant<int>("ECalEndcapNumCalibRhoLayersWheel1");
    numReadoutRhoLayers = aLcdd.constant<int>("ECalEndcapNumReadoutRhoLayersWheel1");
    if ((numReadoutRhoLayers % ECalEndcapNumCalibRhoLayersArr[0]) != 0) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "Number of readout layers must be a multiple of number of calibration layers");
    }
    ECalEndcapNumCalibRhoLayersArr[1] = aLcdd.constant<int>("ECalEndcapNumCalibRhoLayersWheel2");
    numReadoutRhoLayers = aLcdd.constant<int>("ECalEndcapNumReadoutRhoLayersWheel2");
    if ((numReadoutRhoLayers % ECalEndcapNumCalibRhoLayersArr[1]) != 0) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "Number of readout layers must be a multiple of number of calibration layers");
    }
    ECalEndcapNumCalibRhoLayersArr[2] = aLcdd.constant<int>("ECalEndcapNumCalibRhoLayersWheel3");
    numReadoutRhoLayers = aLcdd.constant<int>("ECalEndcapNumReadoutRhoLayersWheel3");
    if ((numReadoutRhoLayers % ECalEndcapNumCalibRhoLayersArr[2]) != 0) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "Number of readout layers must be a multiple of number of calibration layers");
    }
    ECalEndcapNumCalibZLayersArr[0] = aLcdd.constant<int>("ECalEndcapNumCalibZLayersWheel1");
    numReadoutZLayers = aLcdd.constant<int>("ECalEndcapNumReadoutZLayersWheel1");
    if ((numReadoutZLayers % ECalEndcapNumCalibZLayersArr[0]) != 0) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "Number of readout layers must be a multiple of number of calibration layers");
    }
    ECalEndcapNumCalibZLayersArr[1] = aLcdd.constant<int>("ECalEndcapNumCalibZLayersWheel2");
    numReadoutZLayers = aLcdd.constant<int>("ECalEndcapNumReadoutZLayersWheel2");
    if ((numReadoutZLayers % ECalEndcapNumCalibZLayersArr[1]) != 0) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "Number of readout layers must be a multiple of number of calibration layers");
    }
    ECalEndcapNumCalibZLayersArr[2] = aLcdd.constant<int>("ECalEndcapNumCalibZLayersWheel3");
    numReadoutZLayers = aLcdd.constant<int>("ECalEndcapNumReadoutZLayersWheel3");
    if ((numReadoutZLayers % ECalEndcapNumCalibZLayersArr[2]) != 0) {
      dd4hep::printout(dd4hep::ERROR, "ECalEndcap_Turbine_o1_v03",
                       "Number of readout layers must be a multiple of number of calibration layers");
    }

    // Create air envelope for one endcap (will be copied to make both endcaps)
    dd4hep::Tube endcapShape(dim.rmin1(), dim.rmax1(), dim.dz());

    dd4hep::Volume envelopeVol(nameDet + "_vol", endcapShape, aLcdd.material("Air"));

    dd4hep::printout(dd4hep::DEBUG, "ECalEndcap_Turbine_o1_v03",
                     "Placing detector on the positive side: (cm) %f  with min, max radii %f %f and depth %f",
                     dim.z_offset(), dim.rmin1(), dim.rmax1(), dim.dz());

    unsigned iModule = 0;
    buildOneSide_Turbine(aLcdd, caloDetElem, aSensDet, envelopeVol, aXmlElement, iModule);

    dd4hep::Assembly endcapsAssembly("ECalEndcaps_turbine");

    // Place the envelope
    dd4hep::Transform3D envelopePositiveVolume_tr(dd4hep::RotationZYX(0, 0, 0),
                                                  dd4hep::Translation3D(0, 0, dim.z_offset()));
    dd4hep::PlacedVolume envelopePositivePhysVol = endcapsAssembly.placeVolume(envelopeVol, envelopePositiveVolume_tr);
    envelopePositivePhysVol.addPhysVolID("side", 1);

    // make another placement for the negative z endcap
    dd4hep::Transform3D envelopeNegativeVolume_tr(dd4hep::RotationZYX(0, 0, 180 * dd4hep::deg),
                                                  dd4hep::Translation3D(0, 0, -dim.z_offset()));
    dd4hep::PlacedVolume envelopeNegativePhysVol = endcapsAssembly.placeVolume(envelopeVol, envelopeNegativeVolume_tr);
    envelopeNegativePhysVol.addPhysVolID("side", -1);

    dd4hep::Volume motherVol = aLcdd.pickMotherVolume(caloDetElem);
    dd4hep::PlacedVolume envelopePhysVol = motherVol.placeVolume(endcapsAssembly);
    caloDetElem.setPlacement(envelopePhysVol);
    envelopePhysVol.addPhysVolID("system", idDet);

    // Create dummy caloData object for PandoraPFA
    // FIXME: fill calo and layer data information
    auto caloData = new dd4hep::rec::LayeredCalorimeterData;
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::EndcapLayout;
    caloDetElem.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

    // save extent information
    caloData->extent[0] = dim.rmin1();
    caloData->extent[1] = dim.rmax1();
    caloData->extent[2] = dim.z_offset() - dim.dz();
    caloData->extent[3] = dim.z_offset() + dim.dz();

    // Set type flags
    dd4hep::xml::setDetectorTypeFlag(xmlDetElem, caloDetElem);

    return caloDetElem;
  }
} // namespace ECalEndcap_Turbine_o1_v03
} // namespace det

DECLARE_DETELEMENT(ECalEndcap_Turbine_o1_v03, det::ECalEndcap_Turbine_o1_v03::createECalEndcapTurbine)

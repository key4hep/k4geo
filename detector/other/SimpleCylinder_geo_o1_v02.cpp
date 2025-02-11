#include "DD4hep/DetFactoryHelper.h"
#include <DDRec/DetectorData.h>
#include "XML/Utilities.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Vector3D.h"
using dd4hep::_toString;

#define endmsg std::endl
#define lLog std::cout
namespace MSG {
const std::string ERROR = "createSimpleCylinder   ERROR  ";
const std::string DEBUG = "createSimpleCylinder   DEBUG  ";
const std::string INFO  = "createSimpleCylinder   INFO   ";
}

namespace det {
  /**
  Simple cylinder using Tube to be used to define cylinder composed of 1 single material
  Based on SimpleCylinder_geo_o1_v01.cpp prepared by Clement Helsens
  When used for an endcap detector, create both endcap physical volumes and places them 
  in a single detector.
  It also allows the possibility to segment the volume into N layers
  @author A. Durglishvili
  @author G. Marchiori
**/
static dd4hep::Ref_t
createSimpleCylinder(dd4hep::Detector& lcdd, xml_h e, dd4hep::SensitiveDetector sensDet) {
  xml_det_t x_det = e;

  // get detector name, ID and dimensions from compact file
  std::string name = x_det.nameStr();
  int detID = x_det.id();
  xml_comp_t cylinderDim(x_det.child(_U(dimensions)));

  // retrieve layer information
  dd4hep::xml::DetElement layers = x_det.child(_Unicode(layers));
  int nLayers = 0;
  std::vector<double> layerDepth;
  double layersTotalDepth = 0;
  for (dd4hep::xml::Collection_t layer_coll(layers, _Unicode(layer)); layer_coll; ++layer_coll) {
    dd4hep::xml::Component layer = layer_coll;
    nLayers += layer.repeat();
    for (int iLay = 0; iLay < layer.repeat(); iLay++) {
      layerDepth.push_back(layer.thickness());
    }
    layersTotalDepth += layer.repeat() * layer.thickness();
  }
  lLog << MSG::DEBUG << "Number of layers: " << nLayers << endmsg;
  lLog << MSG::DEBUG << "Total thickness from sum of layers in xml description (cm): " << layersTotalDepth/dd4hep::cm << endmsg;
  lLog << MSG::INFO << "Ignoring layer thickness from xml description, assuming all layers have same thickness" << endmsg;

  // create the mother Detector element to be returned at the end
  dd4hep::DetElement detMaster(name, detID);

  // get the world volume, where the detector will be placed
  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(detMaster);

  // create caloData object and fill rmin, rmax info
  auto caloData = new dd4hep::rec::LayeredCalorimeterData;
  caloData->extent[0] = cylinderDim.rmin();
  caloData->extent[1] = cylinderDim.rmax();
  
  dd4hep::rec::MaterialManager matMgr(experimentalHall);
  double zoff = cylinderDim.z_offset();
  double zmin = zoff - cylinderDim.dz();
  double zmax = zoff + cylinderDim.dz();
  bool isEndcap = (zmin*zmax > 0.);
  
  if (isEndcap)
  {
    // top volume of endcaps is an assembly
    dd4hep::Assembly endcapAssembly("Endcaps_assembly");

    // top volume element (do I need this?)
    //  dd4hep::DetElement endcapElement(detMaster);

    // loop over the endcaps
    for (int iEndcap = 0; iEndcap<2; iEndcap++) {
      // create DetElement for endcap, as daughter of detMaster
      dd4hep::DetElement endcap(detMaster);

      // define the tranform for positioning the endcap
      double zoffset = iEndcap == 1 ? zoff : -zoff;
      int rot = iEndcap == 1 ? 0 : 180;
      dd4hep::Transform3D endcap_position(dd4hep::RotationZ( rot*dd4hep::deg), dd4hep::Translation3D(0, 0,  zoffset));

      // define the geometrical shape of the endcap
      dd4hep::Tube cylinder(
        cylinderDim.rmin(), cylinderDim.rmax(),cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

      // define the volume (shape + material) of the detector envelope
      dd4hep::Volume cylinderVol(
        x_det.nameStr() + "_SimpleCylinder", cylinder, lcdd.material("Air"));
      detMaster.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);

      // place the endcap
      auto endcap_pv = endcapAssembly.placeVolume( cylinderVol, endcap_position );

      // mark each placed volume (pv) with the proper phys vol ID
      endcap_pv.addPhysVolID("subsystem", iEndcap);

      // link each pv with its corresponding det element
      endcap.setPlacement( endcap_pv );

      // segment the endcap into layers
      double dzLayer = cylinderDim.dz()*2.0/nLayers;
      lLog << MSG::DEBUG << "dZ of each layer : " << dzLayer << endmsg;
      for (int iLayer=0; iLayer<nLayers; iLayer++) {
        // calculate z extent
        double zMiddle = -cylinderDim.dz() + dzLayer/2.0 + iLayer*dzLayer;
        lLog <<  MSG::DEBUG << "Layer : " << iLayer << " , z offset wrt center of detector : " << zMiddle << endmsg;

        // create detector element as daughter of endcap
        dd4hep::DetElement layer(endcap);

        // define the geometrical shape of the detector layer
        dd4hep::Tube cylinderL(cylinderDim.rmin(), cylinderDim.rmax(), dzLayer/2.0, cylinderDim.phi0(), cylinderDim.deltaphi());

        // define the volume (shape + material) of the detector
        dd4hep::Volume cylinderLVol(
        //  x_det.nameStr() + _toString(iLayer, "_layer%d"), cylinderL, lcdd.material(cylinderDim.materialStr()));
          x_det.nameStr() + _toString(iEndcap, "_side%d") + _toString(iLayer, "_layer%d"), cylinderL, lcdd.material(cylinderDim.materialStr()));
        if (x_det.isSensitive()) {
          dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
          cylinderLVol.setSensitiveDetector(sensDet);
          sensDet.setType(sdType.typeStr());
        }
        detMaster.setVisAttributes(lcdd, x_det.visStr(), cylinderLVol);

        // place the layer volume inside the endcap volume
        dd4hep::Transform3D layerPosition(dd4hep::RotationZ( 000*dd4hep::deg), dd4hep::Translation3D(0, 0, zMiddle));
        auto detLayer_pv = cylinderVol.placeVolume( cylinderLVol, layerPosition );

        // link PV with corresponding det element
        layer.setPlacement( detLayer_pv );

        // set the layer ID
        detLayer_pv.addPhysVolID("layer", iLayer);
      }
    }

    // place the assembly volume in the world
    auto endcapAssembly_pv = experimentalHall.placeVolume(endcapAssembly);

    // assign the system ID to the assembly volume
    endcapAssembly_pv.addPhysVolID("system", detID);

    // link volume with top DetElement to be returned
    detMaster.setPlacement(endcapAssembly_pv);

    // fill the caloData info
    caloData->extent[2] = zmin;
    caloData->extent[3] = zmax;
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::EndcapLayout;

    dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;
    double dzLayer = cylinderDim.dz()*2.0/nLayers;
    auto mat = lcdd.material(cylinderDim.materialStr());
    for (int idxLayer = 0; idxLayer < nLayers; ++idxLayer) {
      double zIn = zmin + idxLayer*dzLayer;
      // double zOut = zIn + dzLayer;

      caloLayer.distance                  = zIn; // distance from origin to innermost face of layer
      caloLayer.sensitive_thickness       = dzLayer; // thickness of the sensitive element
      // caloLayer.absorberThickness         = 0.0; // thickness of absorber part of layer. consider using inner/outer_nRadiationLenghts and inner/outer_nInteractionLengths

      caloLayer.inner_thickness           = dzLayer/2.0; // distance between center of sensitive element and innermost face or layer
      caloLayer.inner_nRadiationLengths   = caloLayer.inner_thickness / mat.radLength(); // absorber material in front of sensitive element in layer
      caloLayer.inner_nInteractionLengths = caloLayer.inner_thickness / mat.intLength(); // absorber material in front of sensitive element in layer
      caloLayer.outer_thickness           = dzLayer/2.0; // distance between center of sensitive element and outermost face or layer
      caloLayer.outer_nRadiationLengths   = caloLayer.outer_thickness / mat.radLength(); // absorber material behind sensitive element in layer
      caloLayer.outer_nInteractionLengths = caloLayer.outer_thickness / mat.intLength(); // absorber material behind sensitive element in layer

      // FIXME! we should get the latter from segmentation class but requires also to set
      // the cell geometry to pointing in theta-phi, which does not exist yet in Pandora
      // For the moment, we use an approximation of fixed-size cells with reasonable numbers
      // so that we can make Pandora run and then we can fix all the missing/outdated pieces
      caloLayer.cellSize0 = 40 * dd4hep::mm;  // 0.5 degrees in theta, R~4.5 m
      caloLayer.cellSize1 = 40 * dd4hep::mm;  // 704 bins in phi (2pi),  R~4.5 m

      // attach the layer to the caloData
      caloData->layers.push_back(caloLayer);
    }
  }
  else
  {
    // top volume of barrel is an assembly
    dd4hep::Assembly barrelAssembly("Barrel_assembly");

    // top volume element
    dd4hep::DetElement barrelElement(detMaster);

    // define the geometrical shape of the barrel
    dd4hep::Tube cylinder(
        cylinderDim.rmin(), cylinderDim.rmax(),cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

    // define the volume (shape + material) of the detector envelope
    dd4hep::Volume cylinderVol(
      x_det.nameStr() + "_SimpleCylinder", cylinder, lcdd.material("Air"));
    detMaster.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);

    // place the barrel
    auto barrel_pv = barrelAssembly.placeVolume( cylinderVol );

    // mark each placed volume (pv) with the proper phys vol ID
    barrel_pv.addPhysVolID("subsystem", 0);

    // link each pv with its corresponding det element
    barrelElement.setPlacement( barrel_pv );

    // create the layers
    double drLayer = (cylinderDim.rmax() - cylinderDim.rmin())/nLayers;
    double rIn = cylinderDim.rmin();
    for (int iLayer=0; iLayer<nLayers; iLayer++) {
      // calculate radial extent
      double rOut = rIn + drLayer;

      // define the geometrical shape of the detector layer
      dd4hep::Tube cylinderL(rIn, rOut, cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

      // define the volume (shape + material) of the detector
      dd4hep::Volume cylinderLVol(
        x_det.nameStr() + _toString(iLayer, "_layer%d"), cylinderL, lcdd.material(cylinderDim.materialStr()));
      if (x_det.isSensitive()) {
        dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
        cylinderLVol.setSensitiveDetector(sensDet);
        sensDet.setType(sdType.typeStr());
      }
      detMaster.setVisAttributes(lcdd, x_det.visStr(), cylinderLVol);

      // create DetElement for layer, as daughter of the barrel assembly
      dd4hep::DetElement detLayer(barrelElement);

      // place the layer volume inside the assembly
      auto detLayer_pv = cylinderVol.placeVolume( cylinderLVol );

      // link PV with corresponding det element
      detLayer.setPlacement( detLayer_pv );

      // set the layer ID
      detLayer_pv.addPhysVolID("layer", iLayer);

      rIn += drLayer;
    }

    // place the assembly volume in the world
    auto barrelAssembly_pv = experimentalHall.placeVolume(barrelAssembly);

    // assign the system ID to the assembly volume
    barrelAssembly_pv.addPhysVolID("system", detID);

    // link volume with top DetElement to be returned
    detMaster.setPlacement(barrelAssembly_pv);

    // Fill caloData object
    caloData->extent[2] = 0;
    caloData->extent[3] = cylinderDim.dz();
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;

    auto mat = lcdd.material(cylinderDim.materialStr());
    dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;
    for (int idxLayer = 0; idxLayer < nLayers; ++idxLayer) {
      rIn = cylinderDim.rmin() + idxLayer*drLayer;
      // double rOut = rIn + drLayer;

      caloLayer.distance                  = rIn; // distance from origin to innermost face of layer
      caloLayer.sensitive_thickness       = drLayer; // thickness of the sensitive element
      // caloLayer.absorberThickness         = 0.0; // thickness of absorber part of layer. consider using inner/outer_nRadiationLenghts and inner/outer_nInteractionLengths

      caloLayer.inner_thickness           = drLayer/2.0; // distance between center of sensitive element and innermost face or layer
      caloLayer.inner_nRadiationLengths   = caloLayer.inner_thickness / mat.radLength(); // absorber material in front of sensitive element in layer
      caloLayer.inner_nInteractionLengths = caloLayer.inner_thickness / mat.intLength(); // absorber material in front of sensitive element in layer
      caloLayer.outer_thickness           = drLayer/2.0; // distance between center of sensitive element and outermost face or layer
      caloLayer.outer_nRadiationLengths   = caloLayer.outer_thickness / mat.radLength(); // absorber material behind sensitive element in layer
      caloLayer.outer_nInteractionLengths = caloLayer.outer_thickness / mat.intLength(); // absorber material behind sensitive element in layer

      caloLayer.cellSize0 = 20 * dd4hep::mm; // FIXME! AD: should be corrected from DDGeometryCreatorALLEGRO. GM: get it from segmentation class
      caloLayer.cellSize1 = 20 * dd4hep::mm; // FIXME! AD: should be corrected from DDGeometryCreatorALLEGRO. GM: get it from segmentation class

      // attach the layer to the caloData
      caloData->layers.push_back(caloLayer);
    }

  }

  // attach the calo data to the detector
  detMaster.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

  // Set type flags
  dd4hep::xml::setDetectorTypeFlag(x_det, detMaster);

  return detMaster;
}
}
DECLARE_DETELEMENT(SimpleCylinder_o1_v02, det::createSimpleCylinder)


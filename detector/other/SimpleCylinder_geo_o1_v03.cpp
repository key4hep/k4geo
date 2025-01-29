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
  Based on SimpleCylinder_geo_o1_v02.cpp, but with 2 layers (todo: read N layers from xml)
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

  const int nLayers(2);

  // create the mother Detector element to be returned at the end
  dd4hep::DetElement detMaster(name, detID);

  // get the world volume, where the detector will be placed
  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(detMaster);

/*
  // define the geometrical shape of the detector (barrel or each endcap)
  dd4hep::Tube cylinder(
    cylinderDim.rmin(), cylinderDim.rmax(),cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

  // define the volume (shape + material) of the detector
  dd4hep::Volume cylinderVol(
    x_det.nameStr() + "_SimpleCylinder", cylinder, lcdd.material(cylinderDim.materialStr()));

  if (x_det.isSensitive()) {
    dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
    cylinderVol.setSensitiveDetector(sensDet);
    sensDet.setType(sdType.typeStr());
  }
  detMaster.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);
*/

  // create caloData object and fill rmin, rmax info
  auto caloData = new dd4hep::rec::LayeredCalorimeterData;
  caloData->extent[0] = cylinderDim.rmin();
  caloData->extent[1] = cylinderDim.rmax();
  
  dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;
  dd4hep::rec::MaterialManager matMgr(experimentalHall);
  double zoff = cylinderDim.z_offset();
  double zmin = zoff - cylinderDim.dz();
  double zmax = zoff + cylinderDim.dz();
  bool isEndcap = (zmin*zmax > 0.);
  
  if (isEndcap)
  {
    // top volume of endcaps is an assembly
    dd4hep::Assembly endcapAssembly("Endcaps_assembly");

    // create DetElements for each endcap, as daughters of detMaster
    dd4hep::DetElement endcapPos(detMaster);
    dd4hep::DetElement endcapNeg(detMaster);

    // define the tranforms for positioning the two endcaps
    dd4hep::Transform3D endcapPos_position(dd4hep::RotationZ( 000*dd4hep::deg), dd4hep::Translation3D(0, 0,  zoff));
    dd4hep::Transform3D endcapNeg_position(dd4hep::RotationZ( 180*dd4hep::deg), dd4hep::Translation3D(0, 0, -zoff));

    // define the geometrical shape of the endcap
    dd4hep::Tube cylinder(
      cylinderDim.rmin(), cylinderDim.rmax(),cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

    // define the volume (shape + material) of the detector
    dd4hep::Volume cylinderVol(
      x_det.nameStr() + "_SimpleCylinder", cylinder, lcdd.material(cylinderDim.materialStr()));

    if (x_det.isSensitive()) {
      dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
      cylinderVol.setSensitiveDetector(sensDet);
      sensDet.setType(sdType.typeStr());
    }
    detMaster.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);

    // place the endcap on the right and left
    auto endcapPos_pv = endcapAssembly.placeVolume( cylinderVol, endcapPos_position );
    auto endcapNeg_pv = endcapAssembly.placeVolume( cylinderVol, endcapNeg_position );

    // mark each placed volume (pv) with the proper phys vol ID
    endcapPos_pv.addPhysVolID("subsystem", 1);
    endcapNeg_pv.addPhysVolID("subsystem", 0);

    // link each pv with its corresponding det element
    endcapPos.setPlacement( endcapPos_pv );
    endcapNeg.setPlacement( endcapNeg_pv );

    // segment the endcaps
    float dzLayer = cylinderDim.dz()*2.0/nLayers;
    lLog << MSG::DEBUG << "dZ of each layer : " << dzLayer << endmsg;
    for (int iLayer=0; iLayer<nLayers; iLayer++) {
      // calculate z extent
      float zMiddle = -cylinderDim.dz() + dzLayer/2.0 + iLayer*dzLayer;
      lLog <<  MSG::DEBUG << "Layer : " << iLayer << " , z offset wrt center of detector : " << zMiddle << endmsg;

      // define the geometrical shape of the detector layer
      dd4hep::Tube cylinderL(cylinderDim.rmin(), cylinderDim.rmax(), dzLayer/2.0, cylinderDim.phi0(), cylinderDim.deltaphi());

      // define the volume (shape + material) of the detector
      dd4hep::Volume cylinderLVol(
        x_det.nameStr() + _toString(iLayer, "_layer%d"), cylinderL, lcdd.material(cylinderDim.materialStr()));
      if (x_det.isSensitive()) {
        dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
        cylinderLVol.setSensitiveDetector(sensDet);
        sensDet.setType(sdType.typeStr());
      }
      detMaster.setVisAttributes(lcdd, x_det.visStr(), cylinderLVol);

      // create DetElement for layer, as daughter of the cylindrical volume
      // dd4hep::DetElement detLayer( cylinderVol );

      // place the layer volume inside the endcap volume
      dd4hep::Transform3D layerPosition(dd4hep::RotationZ( 000*dd4hep::deg), dd4hep::Translation3D(0, 0, zMiddle));
      auto detLayer_pv = cylinderVol.placeVolume( cylinderLVol, layerPosition );

      // link PV with corresponding det element (does not work, needs to be called for DetElement)
      // detLayer.setPlacement( detLayer_pv );

      // set the layer ID
      detLayer_pv.addPhysVolID("layer", iLayer);
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

    dd4hep::rec::Vector3D ivr1 = dd4hep::rec::Vector3D(0., 0., zmin); // defining starting vector points
    dd4hep::rec::Vector3D ivr2 = dd4hep::rec::Vector3D(0., 0., zmax); // defining end vector

    const dd4hep::rec::MaterialVec &materials = matMgr.materialsBetween(ivr1, ivr2); // calling material manager to get material info between two points
    auto mat = matMgr.createAveragedMaterial(materials); // creating average of all the material between two points to calculate X0 and lambda of averaged material
    const double nRadiationLengths = cylinderDim.dz()*2. / mat.radiationLength();
    const double nInteractionLengths = cylinderDim.dz()*2. / mat.interactionLength();

    caloLayer.distance                  = zmin;
    caloLayer.sensitive_thickness       = cylinderDim.dz()*2.;
    caloLayer.inner_thickness           = cylinderDim.dz();
    caloLayer.outer_thickness           = cylinderDim.dz();

    caloLayer.inner_nRadiationLengths   = nRadiationLengths / 2.0;
    caloLayer.inner_nInteractionLengths = nInteractionLengths / 2.0;
    caloLayer.outer_nRadiationLengths   = nRadiationLengths / 2.0;
    caloLayer.outer_nInteractionLengths = nInteractionLengths / 2.0;
  }
  else
  {
    // top volume of barrel is an assembly
    dd4hep::Assembly barrelAssembly("Barrel_assembly");

    // create the layers
    float drLayer = (cylinderDim.rmax() - cylinderDim.rmin())/nLayers;
    float rIn = cylinderDim.rmin();
    for (int iLayer=0; iLayer<nLayers; iLayer++) {
      // calculate radial extent
      float rOut = rIn + drLayer;

      // define the geometrical shape of the detector layer
      dd4hep::Tube cylinder(rIn, rOut, cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

      // define the volume (shape + material) of the detector
      dd4hep::Volume cylinderVol(
        x_det.nameStr() + _toString(iLayer, "_layer%d"), cylinder, lcdd.material(cylinderDim.materialStr()));
      if (x_det.isSensitive()) {
        dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
        cylinderVol.setSensitiveDetector(sensDet);
        sensDet.setType(sdType.typeStr());
      }
      detMaster.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);

      // create DetElement for layer, as daughter of the barrel assembly
      // dd4hep::DetElement detLayer(barrelAssembly);

      // place the layer volume inside the assembly
      auto detLayer_pv = barrelAssembly.placeVolume( cylinderVol );

      // link PV with corresponding det element (does not work, needs to be called for DetElement)
      // cylinderVol.setPlacement( detLayer_pv );

      // set the layer ID
      detLayer_pv.addPhysVolID("layer", iLayer);

      rIn += drLayer;
    }

    // place the assembly volume in the world
    auto barrelAssembly_pv = experimentalHall.placeVolume(barrelAssembly);

    // assign the system ID to the assembly volume
    barrelAssembly_pv.addPhysVolID("system", detID);
    barrelAssembly_pv.addPhysVolID("subsystem", 0);
    
    // link volume with top DetElement to be returned
    detMaster.setPlacement(barrelAssembly_pv);

    // Fill caloData object
    caloData->extent[2] = 0;
    caloData->extent[3] = cylinderDim.dz();
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;

    dd4hep::rec::Vector3D ivr1 = dd4hep::rec::Vector3D(0., cylinderDim.rmin(), 0.); // defining starting vector points 
    dd4hep::rec::Vector3D ivr2 = dd4hep::rec::Vector3D(0., cylinderDim.rmax(), 0.); // defining end vector

    const dd4hep::rec::MaterialVec &materials = matMgr.materialsBetween(ivr1, ivr2); // calling material manager to get material info between two points
    auto mat = matMgr.createAveragedMaterial(materials); // creating average of all the material between two points to calculate X0 and lambda of averaged material
    const double nRadiationLengths = (cylinderDim.rmax() - cylinderDim.rmin()) / mat.radiationLength();
    const double nInteractionLengths = (cylinderDim.rmax() - cylinderDim.rmin()) / mat.interactionLength();

    caloLayer.distance                  = cylinderDim.rmin();
    caloLayer.sensitive_thickness       = (cylinderDim.rmax() - cylinderDim.rmin());
    caloLayer.inner_thickness           = (cylinderDim.rmax() - cylinderDim.rmin()) / 2.0;
    caloLayer.outer_thickness           = (cylinderDim.rmax() - cylinderDim.rmin()) / 2.0;

    caloLayer.inner_nRadiationLengths   = nRadiationLengths / 2.0;
    caloLayer.inner_nInteractionLengths = nInteractionLengths / 2.0;
    caloLayer.outer_nRadiationLengths   = nRadiationLengths / 2.0;
    caloLayer.outer_nInteractionLengths = nInteractionLengths / 2.0;
  }

  caloLayer.cellSize0 = 20 * dd4hep::mm; // FIXME! AD: should be corrected from DDGeometryCreatorALLEGRO. GM: get it from segmentation class
  caloLayer.cellSize1 = 20 * dd4hep::mm; // FIXME! AD: should be corrected from DDGeometryCreatorALLEGRO. GM: get it from segmentation class

  // attach the layer to the caloData
  caloData->layers.push_back(caloLayer);

  // attach the layer to the cylinderDet
  detMaster.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

  // Set type flags
  dd4hep::xml::setDetectorTypeFlag(x_det, detMaster);

  return detMaster;
}
}
DECLARE_DETELEMENT(SimpleCylinder_o1_v03, det::createSimpleCylinder)


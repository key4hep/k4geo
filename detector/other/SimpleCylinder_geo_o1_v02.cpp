#include "DD4hep/DetFactoryHelper.h"
#include <DDRec/DetectorData.h>
#include "XML/Utilities.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Vector3D.h"


namespace det {
/**
  Simple cylinder using Tube to be used to define cylinder composed of 1 single material
  Based on SimpleCylinder_geo_o1_v01.cpp prepared by Clement Helsens
  When used for an endcap detector, create both endcap physical volumes and places them 
  in a single detector
  @author A. Durglishvili
  @author G. Marchiori
**/
static dd4hep::Ref_t
createSimpleCylinder(dd4hep::Detector& lcdd, xml_h e, dd4hep::SensitiveDetector sensDet) {
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  dd4hep::DetElement cylinderDet(name, x_det.id());

  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(cylinderDet);

  xml_comp_t cylinderDim(x_det.child(_U(dimensions)));

  dd4hep::Tube cylinder(
      cylinderDim.rmin(), cylinderDim.rmax(), cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

  dd4hep::Volume cylinderVol(
      x_det.nameStr() + "_SimpleCylinder", cylinder, lcdd.material(cylinderDim.materialStr()));

  if (x_det.isSensitive()) {
    dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
    cylinderVol.setSensitiveDetector(sensDet);
    sensDet.setType(sdType.typeStr());
  }

  // Create caloData object
  auto caloData = new dd4hep::rec::LayeredCalorimeterData;
  dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;
  dd4hep::rec::MaterialManager matMgr(experimentalHall);

  caloData->extent[0] = cylinderDim.rmin();
  caloData->extent[1] = cylinderDim.rmax();

  double zoff = cylinderDim.z_offset();
  double zmin = zoff - cylinderDim.dz();
  double zmax = zoff + cylinderDim.dz();
  bool isEndcap = (zmin*zmax > 0.);
  
  if (isEndcap)
  {
    dd4hep::PlacedVolume cylinderPhys1; // negative endcap
    dd4hep::PlacedVolume cylinderPhys2; // positive endcap
    dd4hep::Position trans1(0., 0., -zoff);
    dd4hep::Position trans2(0., 0., zoff);
    cylinderPhys1 = experimentalHall.placeVolume(cylinderVol, dd4hep::Transform3D(dd4hep::RotationZ(0.), trans1));
    cylinderPhys2 = experimentalHall.placeVolume(cylinderVol, dd4hep::Transform3D(dd4hep::RotationZ(0.), trans2));

    cylinderPhys1.addPhysVolID("system", x_det.id());
    cylinderPhys1.addPhysVolID("subsystem", 0); // negative endcap
    cylinderPhys1.addPhysVolID("layer", 0);

    cylinderPhys2.addPhysVolID("system", x_det.id());
    cylinderPhys1.addPhysVolID("subsystem", 1); // positive endcap
    cylinderPhys2.addPhysVolID("layer", 0);

    cylinderDet.setPlacement(cylinderPhys1);
    cylinderDet.setPlacement(cylinderPhys2);
    cylinderDet.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);

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
    dd4hep::PlacedVolume cylinderPhys;
    cylinderPhys = experimentalHall.placeVolume(cylinderVol);

    cylinderPhys.addPhysVolID("system", x_det.id());
    cylinderPhys.addPhysVolID("layer", 0);
    cylinderDet.setPlacement(cylinderPhys);
    cylinderDet.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);

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
  cylinderDet.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

  // Set type flags
  dd4hep::xml::setDetectorTypeFlag(x_det, cylinderDet);

  return cylinderDet;
}
}
DECLARE_DETELEMENT(SimpleCylinder_o1_v02, det::createSimpleCylinder)


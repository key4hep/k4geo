#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Objects.h"
#include "DDRec/DetectorData.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"

#include "DRTubesconstructor.h"
#include "DRutils.h"

using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h entities, SensitiveDetector sens) {
  xml_det_t x_det = entities;
  int det_id = x_det.id();
  std::string det_name = x_det.nameStr();

  Material air = description.air();

  sens.setType("calorimeter");

  xml_dim_t x_dim = x_det.dimensions();
  double calo_inner_r = x_dim.inner_radius();
  double calo_outer_r = x_dim.outer_radius();
  double calo_half_z = x_dim.z_length();

  // Safety check to make sure dimensions are sensible
  double tower_length = calo_outer_r - calo_inner_r;
  if (tower_length <= 0 * mm)
    throw std::runtime_error("Outer calorimeter radius needs to be larger than inner radius");

  // The barrel volume is an assembly of staves (Trap volumes)
  // A cylindrical volume would only be an approximation and lead to overlaps
  Assembly barrel_volume("calorimeter_barrel");
  barrel_volume.setMaterial(air);
  barrel_volume.setVisAttributes(description, "DRBTassembly_vis");

  DetElement s_detElement(det_name, det_id);
  Volume mother_volume = description.pickMotherVolume(s_detElement);

  // Helper class to construct the calorimeter
  DRBarrelTubes::DRTubesconstructor constructor(&description, entities, &sens);
  constructor.construct_calorimeter(barrel_volume);

  PlacedVolume barrel_placed = mother_volume.placeVolume(barrel_volume);
  barrel_placed.addPhysVolID("system", det_id);
  s_detElement.setPlacement(barrel_placed);

  // apply <type_flags> and expose r-z extent for downstream reco
  dd4hep::xml::setDetectorTypeFlag(entities, s_detElement);

  dd4hep::rec::LayeredCalorimeterData* caloData = new dd4hep::rec::LayeredCalorimeterData;
  caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;
  caloData->extent[0] = calo_inner_r; // rmin
  caloData->extent[1] = calo_outer_r; // rmax
  caloData->extent[2] = 0.;           // zmin (barrel)
  caloData->extent[3] = calo_half_z;  // zmax (half-length)
  s_detElement.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

  return s_detElement;
}

DECLARE_DETELEMENT(DRBarrelTubes, create_detector)

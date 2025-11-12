//====================================================================
//  Modified Generic Endcap Driver for the CLIC detector
//--------------------------------------------------------------------
//
//  Author     : N.Nikiforou
//  Adapted from PolyhedraBarrel Calorimeter by M. Frank
//  Note: This driver supports cutouts in the endcap controlled by z2 and rmin2
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"

using namespace std;

using dd4hep::_toString;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Cone;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Volume;

using dd4hep::rec::LayeredCalorimeterData;

static constexpr auto LOG_SOURCE = "GenericCalEndcap_o1_v01";

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  int det_id = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, e, sdet);
  dd4hep::xml::setDetectorTypeFlag(e, sdet);

  if (theDetector.buildType() == BUILD_ENVELOPE)
    return sdet;

  //-----------------------------------------------------------------------------------

  xml_dim_t dim = x_det.dimensions();
  Material air = theDetector.air();
  int nsides_inner = dim.nsides_inner();
  int nsides_outer = dim.nsides_outer();
  double rmin = dim.rmin();
  double rmax = dim.rmax(); /// FIXME: IS THIS RIGHT?
  double zmin = dim.zmin();

  const bool conicalCutout = [&dim]() {
    if (dim.hasAttr(_Unicode(conicalCutout))) {
      return dim.attr<bool>("conicalCutout");
    }
    return false;
  }();

  // rcutout defines the second inner radius if a conical cutout is used
  double rcutout = dim.hasAttr(_U(rmin2)) ? dim.rmin2() : 0.;
  double zcutout = dim.hasAttr(_U(z2)) ? dim.z2() : 0.;

  Layering layering(x_det);
  double totalThickness = layering.totalThickness();
  sens.setType("calorimeter");
  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();

  std::vector<double> cellSizeVector =
      seg.segmentation()->cellDimensions(0); // Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX = cellSizeVector[0];
  double cell_sizeY = cellSizeVector[1];

  PolyhedraRegular polyVolume(nsides_outer, conicalCutout ? 0 : rmin, rmax, totalThickness);
  Volume endcapVol("endcap", polyVolume, air);

  if (conicalCutout) {
    if (rcutout != rmin) {
      // Make the Cone slightly longer than just the total thickness (==
      // dividing by exactly 2). This makes it extend slightly outside the
      // original volume to avoid potential edge effects in the boolean shape
      // subtraction
      Cone cutoutPolyVolume(totalThickness / 1.99, 0, rmin, 0, rcutout);
      dd4hep::printout(dd4hep::INFO, LOG_SOURCE, "Conical cutout with radius %f", rcutout);
      endcapVol = Volume("endcap", SubtractionSolid(polyVolume, cutoutPolyVolume, Position(0, 0, 0)), air);
    }
  } else {
    if (zcutout > 0. || rcutout > 0.) {
      PolyhedraRegular cutoutPolyVolume(nsides_inner, 0, rmin + rcutout, zcutout);
      Position cutoutPos(0, 0, (zcutout - totalThickness) / 2.0);
      dd4hep::printout(dd4hep::INFO, LOG_SOURCE, "Coutout z width will be %f", zcutout);
      endcapVol = Volume("endcap", SubtractionSolid(polyVolume, cutoutPolyVolume, cutoutPos), air);
    }
  }

  DetElement endcapA(sdet, "endcap", det_id);
  Ref_t(endcapA)->SetName((det_name + "_A").c_str());

  int layer_num = 0;
  int layerType = 0;
  double layerZ = -totalThickness / 2;

  // Create caloData object to extend driver with data required for reconstruction
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout;
  caloData->inner_symmetry = nsides_inner;
  caloData->outer_symmetry = nsides_outer;

  /** NOTE: phi0=0 means lower face flat parallel to experimental floor
   *  This is achieved by rotating the modules with respect to the envelope
   *  which is assumed to be a Polyhedron and has its axes rotated with respect
   *  to the world by 180/nsides. In any other case (e.g. if you want to have
   *  a tip of the calorimeter touching the ground) this value needs to be computed
   */

  caloData->inner_phi0 = 0.;
  caloData->outer_phi0 = 0.;
  caloData->gap0 = 0.; // FIXME
  caloData->gap1 = 0.; // FIXME
  caloData->gap2 = 0.; // FIXME

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rmin;
  caloData->extent[1] = rmax; /// FIXME: CHECK WHAT IS NEEDED (EXSCRIBED?)
  caloData->extent[2] = zmin;
  caloData->extent[3] = zmin + totalThickness;

  endcapVol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
    xml_comp_t x_layer = c;
    double layer_thick = layering.layer(layer_num)->thickness();
    string layer_type_name = _toString(layerType, "layerType%d");
    int layer_repeat = x_layer.repeat();
    const double layer_rmin = [&]() {
      if (conicalCutout) {
        return x_layer.rmin();
      }
      return rmin + (x_layer.hasAttr(_U(gap)) ? x_layer.gap() : 0);
    }();

    dd4hep::printout(dd4hep::INFO, LOG_SOURCE, "Number of layers in group %d : %d", layerType, layer_repeat);

    Volume layer_vol(layer_type_name, PolyhedraRegular(nsides_outer, layer_rmin, rmax, layer_thick), air);

    int slice_num = 0;
    double sliceZ = -layer_thick / 2;

    // Create a caloLayer struct for thiss layer type to store copies of in the parent struct
    LayeredCalorimeterData::Layer caloLayer;
    caloLayer.cellSize0 = cell_sizeX;
    caloLayer.cellSize1 = cell_sizeY;

    double nRadiationLengths = 0.;
    double nInteractionLengths = 0.;
    double thickness_sum = 0;

    for (xml_coll_t s(x_layer, _U(slice)); s; ++s) {
      xml_comp_t x_slice = s;
      string slice_name = _toString(slice_num, "slice%d");
      double slice_thickness = x_slice.thickness();
      Material slice_material = theDetector.material(x_slice.materialStr());
      Volume slice_vol(slice_name, PolyhedraRegular(nsides_outer, layer_rmin, rmax, slice_thickness), slice_material);

      slice_vol.setVisAttributes(theDetector.visAttributes(x_slice.visStr()));
      sliceZ += slice_thickness / 2;
      layer_vol.placeVolume(slice_vol, Position(0, 0, sliceZ));

      nRadiationLengths += slice_thickness / (2. * slice_material.radLength());
      nInteractionLengths += slice_thickness / (2. * slice_material.intLength());
      thickness_sum += slice_thickness / 2;

      if (x_slice.isSensitive()) {
        slice_vol.setSensitiveDetector(sens);

        // Store "inner" quantities
        caloLayer.inner_nRadiationLengths = nRadiationLengths;
        caloLayer.inner_nInteractionLengths = nInteractionLengths;
        caloLayer.inner_thickness = thickness_sum;
        // Store scintillator thickness
        caloLayer.sensitive_thickness = slice_thickness;

        // Reset counters to measure "outside" quantitites
        nRadiationLengths = 0.;
        nInteractionLengths = 0.;
        thickness_sum = 0.;
      }

      nRadiationLengths += slice_thickness / (2. * slice_material.radLength());
      nInteractionLengths += slice_thickness / (2. * slice_material.intLength());
      thickness_sum += slice_thickness / 2;

      sliceZ += slice_thickness / 2;
      slice_num++;
    }

    // Store "outer" quantities
    caloLayer.outer_nRadiationLengths = nRadiationLengths;
    caloLayer.outer_nInteractionLengths = nInteractionLengths;
    caloLayer.outer_thickness = thickness_sum;
    layer_vol.setVisAttributes(theDetector.visAttributes(x_layer.visStr()));

    if (layer_repeat <= 0)
      throw std::runtime_error(x_det.nameStr() + "> Invalid repeat value");

    for (int j = 0; j < layer_repeat; ++j) {
      string phys_lay = _toString(layer_num, "layer%d");

      // The rest of the data is constant; only the distance needs to be updated
      // Store the position up to the inner face of the layer
      caloLayer.distance = zmin + totalThickness / 2 + layerZ;
      // Push back a copy to the caloData structure
      caloData->layers.push_back(caloLayer);

      layerZ += layer_thick / 2;
      DetElement layer_elt(endcapA, phys_lay, layer_num);
      PlacedVolume pv = endcapVol.placeVolume(layer_vol, Position(0, 0, layerZ));
      pv.addPhysVolID("layer", layer_num);
      layer_elt.setPlacement(pv);

      layerZ += layer_thick / 2;
      ++layer_num;
    }
    ++layerType;
  }

  double z_pos = zmin + totalThickness / 2;
  PlacedVolume pv;
  // Reflect it.

  DetElement endcapB = endcapA.clone(det_name + "_B", x_det.id());

  // Removed rotations to align with envelope
  // NOTE: If the envelope is not a polyhedron (eg. if you use a tube)
  // you may need to rotate so the axes match

  pv = envelope.placeVolume(endcapVol, Transform3D(RotationZYX(0, 0, 0), Position(0, 0, z_pos)));
  pv.addPhysVolID("side", 1);
  endcapA.setPlacement(pv);

  // Removed rotations
  pv = envelope.placeVolume(endcapVol, Transform3D(RotationZYX(0, M_PI, 0), Position(0, 0, -z_pos)));
  pv.addPhysVolID("side", 2);
  endcapB.setPlacement(pv);

  sdet.add(endcapB);

  sdet.addExtension<LayeredCalorimeterData>(caloData);

  return sdet;
}

DECLARE_DETELEMENT(GenericCalEndcap_o1_v01, create_detector)

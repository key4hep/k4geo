
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"
#include <cmath>
#include <string>

using namespace std;

using dd4hep::_toString;
using dd4hep::Assembly;
using dd4hep::Box;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::IntersectionSolid;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Polycone;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Solid;
using dd4hep::SubtractionSolid;
using dd4hep::Torus;
using dd4hep::Transform3D;
using dd4hep::Trapezoid;
using dd4hep::Tube;
using dd4hep::UnionSolid;
using dd4hep::Volume;

using dd4hep::rec::LayeredCalorimeterData;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::volSurfaceList;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector) {
  // static double tolerance = 0e0;

  xml_det_t x_det = e;
  int det_id = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, e, sdet);

  if (theDetector.buildType() == BUILD_ENVELOPE)
    return sdet;

  Material air = theDetector.air();
  PlacedVolume pv;
  int n = 0;

  for (xml_coll_t i(x_det, _U(layer)); i; ++i, ++n) {
    xml_comp_t x_layer = i;
    string l_name = det_name + _toString(n, "_layer%d");
    double z = x_layer.outer_z();
    double rmin = x_layer.inner_r();
    double r = rmin;
    DetElement layer(sdet, _toString(n, "layer%d"), x_layer.id());
    Tube l_tub(rmin, 2 * rmin, z);
    Volume l_vol(l_name, l_tub, air);
    int m = 0;

    for (xml_coll_t j(x_layer, _U(slice)); j; ++j, ++m) {
      xml_comp_t x_slice = j;
      Material mat = theDetector.material(x_slice.materialStr());
      string s_name = l_name + _toString(m, "_slice%d");
      double thickness = x_slice.thickness();

      Tube s_tub(r, r + thickness, z, 2 * M_PI);
      Volume s_vol(s_name, s_tub, mat);

      // Add surface to the support
      Vector3D ocyl(r + thickness / 2., 0., 0.);

      VolCylinder cylSurf1(s_vol, SurfaceType(SurfaceType::Helper), 0.5 * thickness, 0.5 * thickness, ocyl);

      volSurfaceList(sdet)->push_back(cylSurf1);

      r += thickness;

      // Set Attributes
      s_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
      pv = l_vol.placeVolume(s_vol);
      // Slices have no extra id. Take the ID of the layer!
      pv.addPhysVolID("slice", m);
    }
    l_tub.setDimensions(rmin, r, z);
    // cout << l_name << " " << rmin << " " << r << " " << z << endl;
    l_vol.setVisAttributes(theDetector, x_layer.visStr());

    pv = envelope.placeVolume(l_vol);
    pv.addPhysVolID("layer", n);
    layer.setPlacement(pv);
  }

  return sdet;
}

DECLARE_DETELEMENT(TrackerBarrelSupport_o1_v01, create_detector)

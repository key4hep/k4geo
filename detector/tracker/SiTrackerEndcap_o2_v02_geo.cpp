// $Id: $
//==========================================================================
//  AIDA Detector description implementation for LCD
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
// Tweaked by : D.Protopopescu
//==========================================================================
//
// Specialized generic detector constructor
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include <map>

using namespace std;
using dd4hep::_toString;
using dd4hep::Box;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationY;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Trapezoid;
using dd4hep::Volume;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens) {
  typedef vector<PlacedVolume> Placements;
  xml_det_t x_det = e;
  Material vacuum = theDetector.vacuum();
  int det_id = x_det.id();
  string det_name = x_det.nameStr();
  bool reflect = x_det.reflect(false);
  DetElement sdet(det_name, det_id);
  int m_id = 0, c_id = 0, n_sensor = 0;
  map<string, Volume> modules;
  map<string, Placements> sensitives;
  PlacedVolume pv;

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, e, sdet);
  dd4hep::xml::setDetectorTypeFlag(e, sdet);

  if (theDetector.buildType() == BUILD_ENVELOPE)
    return sdet;

  //-----------------------------------------------------------------------------------

  envelope.setVisAttributes(theDetector.invisible());
  sens.setType("tracker");

  // Build the sensor units
  // Loop over 'modules' as defined in the XML
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string m_nam = x_mod.nameStr();
    xml_comp_t trd = x_mod.trd();
    double posY;
    double x1 = trd.x1();
    double x2 = trd.x2();
    double z = trd.z();
    double y1, y2, total_thickness = 0.;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci)
      total_thickness += xml_comp_t(ci).thickness();

    y1 = y2 = total_thickness / 2;
    Volume m_volume(m_nam, Trapezoid(x1, x2, y1, y2, z), vacuum);
    m_volume.setVisAttributes(theDetector.visAttributes(x_mod.visStr()));

    std::cout << m_nam << ", thickness=" << total_thickness << std::endl;

    // Loop over the module_components ('slices') in the 'module'
    // The first component (top in the XML) is placed at the 'bottom'
    for (ci.reset(), n_sensor = 1, c_id = 0, posY = -y1; ci; ++ci, ++c_id) {
      xml_comp_t c = ci;
      double c_thick = c.thickness();
      Material c_mat = theDetector.material(c.materialStr());
      string c_name = _toString(c_id, "component%d");
      Volume c_vol(c_name, Trapezoid(x1, x2, c_thick / 2e0, c_thick / 2e0, z), c_mat);

      std::cout << " + sensor " << n_sensor << " " << c_name;

      c_vol.setVisAttributes(theDetector.visAttributes(c.visStr()));
      pv = m_volume.placeVolume(c_vol, Position(0, posY + c_thick / 2, 0));
      if (c.isSensitive()) {
        sdet.check(n_sensor > 2, "SiTrackerEndcap::fromCompact: " + c_name + " Max of 2 modules allowed!");
        pv.addPhysVolID("slice", n_sensor);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        std::cout << " (" << n_sensor << " is sensitive) ";
        ++n_sensor;
      }
      std::cout << std::endl;
      posY += c_thick;
    }
    modules[m_nam] = m_volume;
  }
  // done building the 2 modules, of 12 layers each

  // Build now the detector itself
  // Loop over layers as defined in the XML
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer(li);
    int l_id = x_layer.id();
    int ring_num = 0;

    std::cout << "Layer " << l_id << ":" << std::endl;

    // Loop over rings, as defined in the XML
    for (xml_coll_t ri(x_layer, _U(ring)); ri; ++ri) {
      xml_comp_t x_ring = ri;
      double r = x_ring.r();
      double phi0 = x_ring.phi0(0);
      double zstart = x_ring.zstart();
      double dz = x_ring.dz(0);
      int nmodules = x_ring.nmodules();
      string m_nam = x_ring.moduleStr();
      Volume m_vol = modules[m_nam];
      double iphi = 2 * M_PI / nmodules;
      double phi = phi0;
      Placements& sensVols = sensitives[m_nam];

      std::cout << " Ring " << ring_num << ":" << std::endl;

      // We encode side->side, layer->layer, ring->module, module->sensor;  to
      // adhere to the standard bit encoder string (side=-1 and +1 for endcaps)
      // <id>system:5,side:-2,layer:9,module:8,sensor:8</id>

      string r_base = _toString(l_id, "layer%d") + _toString(ring_num, "_module%d");

      // Loop over modules in each ring, modules are either type 1 or 2
      for (int k = 0; k < nmodules; ++k) {

        double x = -r * std::cos(phi);
        double y = -r * std::sin(phi);

        std::cout << " + " << r_base << ", ring_num=" << ring_num << ", k=" << k << ", type=" << m_nam << std::endl;

        for (int s = 1 - 2 * int(reflect); s < 2; s += 1 + int(reflect)) {

          DetElement module(sdet, _toString(s, "side%d") + r_base + _toString(k, "_sensor%d"), det_id);
          pv = envelope.placeVolume(
              m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2), Position(x, y, s * (zstart + dz))));
          pv.addPhysVolID("side", s)
              .addPhysVolID("layer", l_id)
              .addPhysVolID("module", ring_num)
              .addPhysVolID("sensor", k);
          ;
          module.setPlacement(pv);

          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            // sens_pv.setPhysVolID("sensor", 100*ic + k);  // <-- need such function
            DetElement comp_elt(module, sens_pv.volume().name(), det_id);
            comp_elt.setPlacement(sens_pv);
            std::cout << "    side" << s << "_sens_pv " << ic << ": " << sens_pv.volume().name() << std::endl;
          }
        }

        /*
        DetElement module(sdet, "pos_" + r_base + _toString(k, "_sensor%d"), det_id);
        pv = envelope.placeVolume(m_vol, Transform3D(RotationZYX(0,-M_PI/2-phi,-M_PI/2), Position(x,y,zstart+dz)));
        pv.addPhysVolID("side", 1).addPhysVolID("layer", l_id).addPhysVolID("module", ring_num);
        module.setPlacement(pv);

        for(size_t ic=0; ic<sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          // sens_pv.setPhysVolID("sensor", 100*ic + k);  // <-- need such function
          DetElement comp_elt(module, sens_pv.volume().name(), det_id);
          comp_elt.setPlacement(sens_pv);
          std::cout << "    pos_sens_pv " << ic << ": " << sens_pv.volume().name() << std::endl;
        }

        if (reflect) {
          DetElement r_module(sdet, "neg_" + r_base + _toString(k, "_sensor%d"), det_id);
          pv = envelope.placeVolume(m_vol, Transform3D(RotationZYX(M_PI,-M_PI/2-phi,-M_PI/2),Position(x,y,-zstart-dz)));
          pv.addPhysVolID("side", -1).addPhysVolID("layer", l_id).addPhysVolID("module", ring_num);
          r_module.setPlacement(pv);

          for(size_t ic=0; ic<sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            // sens_pv.setPhysVolID("sensor", 100*ic + k);  // <-- need such function
            DetElement comp_elt(r_module, sens_pv.volume().name(), det_id);
            comp_elt.setPlacement(sens_pv);
            std::cout << "    neg_sens_pv " << ic << ": " << sens_pv.volume().name() << std::endl;
          }
        }
        */

        dz = -dz;
        phi += iphi;
      }
      ++ring_num;
    }
  }

  return sdet;
}

DECLARE_DETELEMENT(SiTrackerEndcap_o2_v02, create_detector)

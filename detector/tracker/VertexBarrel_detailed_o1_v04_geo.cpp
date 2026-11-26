// $Id: $
//====================================================================
//  Based on ZPlanarTracker module from F. Gaede

//  Tracking detector to describe the FCC-ee IDEA vertex detector barrel.
//  The vertex detector is assembled of stave structures which feature
//  support and readout (flex) elements. Each stave features multiple
//  individual modules, that consist of sensitive and insensitive
//  sensor elements.
//  From o1_v01 to o1_v02: It is possible to define curved
//  sensor and support elements (using the 'isCurved' option).
//  From o1_v02 to o1_v03: One can now have a partially-
//  insensitive sensor using 'sensor_insensitive_thickness_below' and
//  'sensor_insensitive_thickness_above'. Components with zero thickness
//  are now ignored when building the detector. Curved sensors can now
//  also be approximated by trapezoidal volumes with the 'nsegments'
//  option. This allows to perform track reconstruction with conformal
//  tracking, which currently doesn't support curved sensors (=sensitive
//  surfaces) yet. Sensor components can now also be defined in an
//  external xml file to avoid repetion in the detector xml file.
//  From o1_v03 to o1_v04: Fixed the treatment of insensitive volumes.
//  They now receive a surface, which is needed for the material evaluation
//  in DDKalTest and thus conformal tracking. This surface can be skipped
//  with the 'notPlaceSurface' option, e.g. if the on-stave volumes have
//  already been included in the 'other_insensitive_thickness_below'/above
//  of the sensor definition.
//--------------------------------------------------------------------
//
//  Author     : Armin Ilg
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"
#include "XMLHandlerDB.h"
#include <exception>

using namespace std;

using dd4hep::_toString;
using dd4hep::Assembly;
using dd4hep::Box;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::DEBUG;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::getAttrOrDefault;
using dd4hep::INFO;
using dd4hep::Material;
using dd4hep::NamedObject;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationX;
using dd4hep::RotationY;
using dd4hep::RotationZ;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Solid;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Trapezoid;
using dd4hep::Tube;
using dd4hep::Volume;
using dd4hep::WARNING;

using dd4hep::Trd1;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::VolPlane;
using dd4hep::rec::volSurfaceList;

static Ref_t create_element(Detector& theDetector, xml_h e, SensitiveDetector sens) {

  xml_det_t x_det = e;
  std::string det_name = x_det.nameStr();

  DetElement sdet(det_name, x_det.id());
  PlacedVolume pv;

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, e, sdet);
  dd4hep::xml::setDetectorTypeFlag(e, sdet);

  if (theDetector.buildType() == BUILD_ENVELOPE)
    return sdet;

  envelope.setVisAttributes(theDetector, x_det.visStr());
  sens.setType("tracker");

  // Struct to support multiple readout or support layers (both using this struct)
  struct componentsStruct {
    string name;
    double r;
    double offset;
    double length;
    double z_offset;
    vector<double> thicknesses;
    vector<double> offsets;
    vector<double> z_offsets;
    vector<double> rs;
    vector<double> phis;
    vector<Volume> volumes;
    vector<bool> isCurved;
    vector<bool> notPlaceSurface;
  };

  // Struct to support end-of-stave structures
  struct endOfStaveStruct {
    string name;
    double r;
    double offset;
    double z_offset;
    vector<double> thicknesses;
    vector<double> lengths;
    vector<double> dzs;
    vector<double> offsets;
    vector<double> z_offsets;
    vector<double> rs;
    vector<Volume> volumes;
    vector<int> side;
    vector<bool> isCurved;
    vector<bool> notPlaceSurface;
  };

  // Struct for sensors
  struct sensorsStruct {
    string name;
    double r;
    double offset;
    int nx;
    Material material;
    vector<bool> sensitives;
    vector<double> thicknesses;
    vector<double> insensitive_thicknesses_below;
    vector<double> insensitive_thicknesses_above;
    vector<double> xmin;
    vector<double> xmax;
    vector<double> ymin;
    vector<double> ymax;
    vector<double> rs;
    vector<bool> isCurved;
    double width = 0.;
    double length = 0.;
    vector<Volume> volumes;
    vector<int> nsegments;
    int nSensitivePerModule = 0;
    int nGroupingModules;
  };

  // --- Module information struct ---
  struct stave_information {
    string name;
    double motherVolThickness;
    double motherVolWidth;

    vector<componentsStruct> components_vec;
    vector<endOfStaveStruct> endOfStaves_vec;
    vector<sensorsStruct> sensorsVec;

    double stave_dr;
    double stave_r;
    double stave_length;
    string stave_vis;
    string stave_layer_encoding;
  };
  list<stave_information> stave_information_list;

  // --- Collect stave(s) information
  for (xml_coll_t mi(x_det, _U(stave)); mi; ++mi) {
    xml_comp_t x_stave = mi;

    stave_information m;
    m.name = x_stave.nameStr();
    m.stave_length = x_stave.length();

    m.stave_dr = x_stave.dr(0); // Offset for every second module in r
    m.stave_r = x_stave.r(0);   // Radius of stave

    m.motherVolThickness = getAttrOrDefault(x_stave, _Unicode(motherVolThickness), double(0.0));
    m.motherVolWidth = getAttrOrDefault(x_stave, _Unicode(motherVolWidth), double(0.0));

    m.stave_vis = x_stave.visStr(x_det.visStr());

    m.stave_layer_encoding = getAttrOrDefault<string>(x_stave, _Unicode(layer_id), "");

    // Components
    xml_coll_t c_components(x_stave, _U(components));
    for (c_components.reset(); c_components; ++c_components) {
      componentsStruct components;
      components.name = xml_comp_t(c_components).nameStr();
      components.r = xml_comp_t(c_components).r(0);
      components.length = xml_comp_t(c_components).length();
      components.offset = xml_comp_t(c_components).offset(0);
      components.z_offset = xml_comp_t(c_components).z_offset(0);

      xml_coll_t c_component(c_components, _U(component));
      int iComponent = 0;
      for (c_component.reset(); c_component; ++c_component) {
        xml_comp_t component = c_component;
        double thickness =
            (component.hasAttr(_Unicode(width)))
                ? component.thickness()
                : (2. * component.outer_radius()); // Use 2 * R_out for cylindrical components in rectangular stave
        if (thickness == 0.)
          continue; // Skip components with zero thickness
        components.thicknesses.push_back(thickness);
        components.offsets.push_back(component.offset(0));
        components.z_offsets.push_back(component.z_offset(0));
        components.rs.push_back(component.r(0));
        components.phis.push_back(
            component.phi(0)); // Rotation of component around its own axis, not applicable for curved components

        bool isCurved = getAttrOrDefault(component, _Unicode(isCurved), bool(false));
        components.isCurved.push_back(isCurved);
        components.notPlaceSurface.push_back(getAttrOrDefault(
            component, _Unicode(notPlaceSurface),
            bool(false))); // Don't put a surface for this volume. Only use if volume material is included in sensor
                           // 'other_insensitive_thickness_below'/'other_insensitive_thickness_above'
        Volume ele_vol;

        if (isCurved) {
          double rmin = m.stave_r + components.r + components.rs.back();
          double half_width = component.width() / (2. * M_PI * rmin) * (2.0 * M_PI) / 2.;
          double phi_offset = getAttrOrDefault(component, _Unicode(phi_offset), double(0.0));
          Tube ele_box =
              Tube(rmin, rmin + thickness, components.length, -half_width + phi_offset, half_width + phi_offset);
          ele_vol = Volume(components.name + _toString(iComponent, "_%d"), ele_box,
                           theDetector.material(component.materialStr()));
        } else {
          Solid ele_shape;
          if (component.hasAttr(_Unicode(width))) {
            ele_shape = Box(thickness / 2., component.width() / 2., components.length);
          } else {
            ele_shape = Tube(component.inner_radius(), component.outer_radius(),
                             components.length); // Thickness ignored to follow TGeo constructor more closely
          }
          ele_vol = Volume(components.name + _toString(iComponent, "_%d"), ele_shape,
                           theDetector.material(component.materialStr()));
        }
        ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());
        components.volumes.push_back(ele_vol);

        iComponent++;
      }
      m.components_vec.push_back(components);
    }

    // End of stave structures
    xml_coll_t c_endOfStave(x_stave, _U(end_z));
    for (c_endOfStave.reset(); c_endOfStave; ++c_endOfStave) {
      endOfStaveStruct endOfStave;
      endOfStave.offset = xml_comp_t(c_endOfStave).offset(0);
      endOfStave.z_offset = xml_comp_t(c_endOfStave).z_offset(0);
      endOfStave.r = xml_comp_t(c_endOfStave).r(0);
      endOfStave.name = xml_comp_t(c_endOfStave).nameStr();

      xml_coll_t c_component = xml_coll_t(c_endOfStave, _U(component));
      int iEndOfStave = 0;
      for (c_component.reset(); c_component; ++c_component) {
        xml_comp_t component = c_component;
        double thickness = component.thickness();
        if (thickness == 0.)
          continue; // Skip components with zero thickness
        endOfStave.thicknesses.push_back(thickness);
        endOfStave.dzs.push_back(component.dz(0));
        endOfStave.offsets.push_back(component.offset(0));
        endOfStave.z_offsets.push_back(component.z_offset(0));
        endOfStave.lengths.push_back(component.length());
        endOfStave.rs.push_back(component.r(0));

        int side = getAttrOrDefault(component, _Unicode(side), int(0));

        endOfStave.side.push_back(side); // 0 for both sides (default), +1 for +z side, -1 for -z side
        bool isCurved = getAttrOrDefault(component, _Unicode(isCurved), bool(false));
        endOfStave.notPlaceSurface.push_back(getAttrOrDefault(
            component, _Unicode(notPlaceSurface),
            bool(false))); // Don't put a surface for this volume. Only use if volume material is included in sensor
                           // 'other_insensitive_thickness_below'/'other_insensitive_thickness_above'
        endOfStave.isCurved.push_back(isCurved);
        Volume ele_vol;

        if (isCurved) {
          double rmin = m.stave_r + endOfStave.r + endOfStave.rs.back();
          double half_width = component.width() / (2. * M_PI * rmin) * (2.0 * M_PI) / 2.;
          double phi_offset = getAttrOrDefault(component, _Unicode(phi_offset), double(0.0));
          Tube ele_box =
              Tube(rmin, rmin + thickness, component.length() / 2., -half_width + phi_offset, half_width + phi_offset);
          ele_vol = Volume(endOfStave.name + _toString(iEndOfStave, "_%d"), ele_box,
                           theDetector.material(component.materialStr()));
        } else {
          Box ele_box = Box(thickness / 2., component.width() / 2., component.length() / 2.);
          ele_vol = Volume(endOfStave.name + _toString(iEndOfStave, "_%d"), ele_box,
                           theDetector.material(component.materialStr()));
        }
        ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());

        endOfStave.volumes.push_back(ele_vol);
        iEndOfStave++;
      }
      m.endOfStaves_vec.push_back(endOfStave);
    }

    // Sensor
    xml_coll_t c_sensor(x_stave, _U(sensor));
    for (c_sensor.reset(); c_sensor; ++c_sensor) {
      sensorsStruct sensor;

      sensor.r = xml_comp_t(c_sensor).r(0);
      sensor.offset = xml_comp_t(c_sensor).offset(0);
      sensor.material = theDetector.material(xml_comp_t(c_sensor).materialStr());
      sensor.name = xml_comp_t(c_sensor).nameStr();
      sensor.nx = getAttrOrDefault(
          xml_comp_t(c_sensor), _Unicode(nx),
          int(1)); // Duplicate entries nx times in phi (local x) direction, default is 1 (no duplication)
      sensor.nGroupingModules =
          getAttrOrDefault(xml_comp_t(c_sensor), _Unicode(nGroupingModules),
                           int(1)); // Number of modules that are together using the same module id, but different
                                    // sensor ids (to more efficiently use GlobalTrackerReadoutID bits)
      double default_sensor_thickness = xml_comp_t(c_sensor).thickness();

      // Try to load sensor components from external include file first
      bool use_include = false;
      dd4hep::xml::Handle_t component_source;
      std::unique_ptr<dd4hep::xml::DocumentHolder> doc;
      try {
        xml_coll_t incl(c_sensor, _U(include));
        if (incl) {
          doc = std::make_unique<dd4hep::xml::DocumentHolder>(
              dd4hep::xml::DocumentHandler().load(incl, incl.attr_value(_U(ref))));
          printout(DEBUG, det_name,
                   "Using sensor description from external file: " + _toString(incl.attr_value(_U(ref))));
          component_source = doc->root();
          use_include = true;
        }
      } catch (...) {
        printout(DEBUG, det_name, "No sensor description from external file found, using inline description instead");
        use_include = false;
      }

      // Either use external include file or inline sensor description
      xml_coll_t c_component =
          use_include ? xml_coll_t(component_source, _U(component)) : xml_coll_t(c_sensor, _U(component));
      int iSensor = 0;

      for (c_component.reset(); c_component; ++c_component) {
        xml_comp_t component = c_component;

        double sensitive_thickness = getAttrOrDefault(
            component, _Unicode(thickness), default_sensor_thickness); // Take the overall sensor thickness as default

        double insensitive_sensor_below = getAttrOrDefault(component, _Unicode(sensor_insensitive_thickness_below), 0.);
        double insensitive_sensor_above = getAttrOrDefault(component, _Unicode(sensor_insensitive_thickness_above), 0.);

        // Thickness below/above the sensor (not part of sensor) that is insensitive. E.g. stave support structure
        double other_below = getAttrOrDefault(component, _Unicode(other_insensitive_thickness_below), 0.);
        double other_above = getAttrOrDefault(component, _Unicode(other_insensitive_thickness_above), 0.);
        sensor.nSensitivePerModule += component.isSensitive() ? 1 : 0;

        int nsegment = getAttrOrDefault(component, _Unicode(nsegments),
                                        int(1)); // To approximate curved sensors by trapezoid segments
        bool isCurved = getAttrOrDefault(component, _Unicode(isCurved), bool(false));
        double r_layer = component.r(0); // running radial offset within this component

        // Lambda function to add the sensitive part of the sensor and the insensitive parts below and above it
        auto addSensorPart = [&](double thickness, bool sensitive, double ins_below, double ins_above,
                                 const string& name_suffix) {
          if (thickness == 0.)
            return; // Skip parts with zero thickness
          string part_name = component.nameStr() + _toString(iSensor, "_%d") + name_suffix;
          sensor.thicknesses.push_back(thickness);
          sensor.insensitive_thicknesses_below.push_back(
              ins_below); // Volume below sensor that is insensitive, e.g. a sensor substrate
          sensor.insensitive_thicknesses_above.push_back(
              ins_above); // Volume above sensor that is insensitive, e.g. a metal layer
          sensor.sensitives.push_back(sensitive);
          sensor.xmin.push_back(component.xmin());
          sensor.xmax.push_back(component.xmax());
          sensor.ymin.push_back(component.ymin());
          sensor.ymax.push_back(component.ymax());
          sensor.rs.push_back(r_layer);
          r_layer += thickness;
          sensor.nsegments.push_back(nsegment);
          sensor.isCurved.push_back(isCurved);

          Volume ele_vol;
          if (isCurved) {
            double rmin = m.stave_r + sensor.r + sensor.rs.back();
            double half_width = abs(component.xmax() - component.xmin()) / rmin / 2.;
            if (nsegment > 1) {
              // The trapezoid has an average width of 'abs(component.xmax() - component.xmin()) / 2. / nsegment', so
              // the inner width must be a bit smaller and the outer width a bit larger depending on the radius where it
              // is placed. This is given by the scaling factor '(rmin + thickness) / rmin'
              Trd1 ele_box =
                  Trd1(abs(component.xmax() - component.xmin()) / 2. / nsegment * (1 + (1 - (rmin + thickness) / rmin)),
                       abs(component.xmax() - component.xmin()) / 2. / nsegment * (1 - (1 - (rmin + thickness) / rmin)),
                       abs(component.ymax() - component.ymin()) / 2., thickness / 2.);
              ele_vol = Volume(part_name, ele_box, sensor.material);
            } else {
              double phi_offset = getAttrOrDefault(component, _Unicode(phi_offset), double(0.0));
              Tube ele_box = Tube(rmin, rmin + thickness, abs(component.ymax() - component.ymin()) / 2.,
                                  -half_width + phi_offset, half_width + phi_offset);
              ele_vol = Volume(part_name, ele_box, sensor.material);
            }
          } else {
            Box ele_box = Box(abs(component.xmax() - component.xmin()) / 2.,
                              abs(component.ymax() - component.ymin()) / 2., thickness / 2.);
            ele_vol = Volume(part_name, ele_box, sensor.material);
          }

          if (sensor.sensitives.back())
            ele_vol.setSensitiveDetector(sens);
          ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());
          sensor.volumes.push_back(ele_vol);
          iSensor++;
        };

        // Finally add all the sensor parts. The material of the insensitive sensor parts below/above the sensitive part
        // is included in the sensitive sensor part
        addSensorPart(insensitive_sensor_below, false, 0., 0., "_insensitive_below");
        addSensorPart(sensitive_thickness, component.isSensitive(), other_below + insensitive_sensor_below,
                      other_above + insensitive_sensor_above, "");
        addSensorPart(insensitive_sensor_above, false, 0., 0., "_insensitive_above");
      }
      sensor.width =
          *max_element(sensor.xmax.begin(), sensor.xmax.end()) - *min_element(sensor.xmin.begin(), sensor.xmin.end());
      sensor.length =
          *max_element(sensor.ymax.begin(), sensor.ymax.end()) - *min_element(sensor.ymin.begin(), sensor.ymin.end());

      printout(DEBUG, det_name,
               "Module: " + sensor.name + ", sensor width: " + to_string(sensor.width) +
                   ", sensor length: " + to_string(sensor.length));
      m.sensorsVec.push_back(sensor);
    }

    stave_information_list.push_back(m);
    printout(DEBUG, det_name, "Read stave information of stave " + m.name);
  }

  //=========  loop over layer elements in xml  ======================================

  printout(INFO, det_name, "Building of detector ...");
  int nLayers = 0;
  vector<int> layer_ids;
  DetElement layerDE;
  std::map<int, int> iModuleTot_map; // Map to keep track of the total number of modules for each layer id, to give
                                     // unique IDs to modules even when there are multiple layers with the same layer ID
  for (xml_coll_t c(e, _U(layer)); c; ++c) {

    xml_comp_t x_layer(c);

    if (x_layer.attr<bool>(_Unicode(ignore), false))
      continue; // Skip layers marked to be ignored. This can be used for example to easily switch
                // off layers that are outside the envelopes of the detector
    else
      nLayers++;

    int layer_id = x_layer.id();
    double dr = x_layer.dr(0); // Spacing in r for every second stave.
    double layer_offset = x_layer.offset(0);
    double z_offset = x_layer.z_offset(0);

    string nameStr = x_layer.nameStr();
    int nmodules = x_layer.nmodules();
    int iModuleTot = iModuleTot_map[layer_id]; // Counter to give unique IDs to modules
    double step_z = getAttrOrDefault(
        x_layer, dd4hep::xml::Strng_t("step_z"),
        x_layer.step(
            0)); // Spacing of modules in z, accept both "step" (backwards-comptible) and "step_z" as attributes
    double step_rphi = getAttrOrDefault(x_layer, dd4hep::xml::Strng_t("step_rphi"),
                                        double(0.)); // Spacing of modules in rphi (Optional)

    // Use the correct stave
    auto m = *find_if(stave_information_list.cbegin(), stave_information_list.cend(),
                      [&nameStr](const stave_information& stave) { return stave.name == nameStr; });

    double motherVolThickness = getAttrOrDefault(x_layer, _Unicode(motherVolThickness), double(0.0));
    double motherVolLength = getAttrOrDefault(x_layer, _Unicode(motherVolLength), double(m.stave_length));
    double motherVolOffset =
        getAttrOrDefault(x_layer, _Unicode(motherVolOffset), double(0.0)); // In case wafer/stave is asymmetric

    std::string layer_name = _toString(layer_id, "layer%d");

    Volume whole_layer_volume_v;
    Assembly whole_layer_volume_a;
    if (motherVolThickness > 0.0) {
      double motherVolRmin = getAttrOrDefault(x_layer, _Unicode(motherVolRmin), double(x_layer.r()));
      Tube whole_layer_tube = Tube(motherVolRmin, motherVolRmin + motherVolThickness, motherVolLength / 2.);
      whole_layer_volume_v = Volume(layer_name, whole_layer_tube, theDetector.material("Air"));
      whole_layer_volume_v.setVisAttributes(theDetector, x_det.visStr());
      pv = envelope.placeVolume(whole_layer_volume_v, Position(0., 0., z_offset));
    } else {
      //// Use just assembly to avoid overlap between layers
      whole_layer_volume_a = Assembly(layer_name);
      pv = envelope.placeVolume(whole_layer_volume_a, Position(0., 0., z_offset));
    }
    if (find(layer_ids.begin(), layer_ids.end(), layer_id) != layer_ids.end()) {
      printout(
          INFO, det_name,
          "Defining multiple layers with the same layer ID: layer_id: " + _toString(layer_id) +
              ", this is okay. Just make sure you have enough bits in GlobalTrackerReadoutID for layers and modules.");
    }
    if (m.stave_layer_encoding.empty()) {
      pv.addPhysVolID("layer", layer_id).addPhysVolID("side", 0);
      layer_ids.push_back(layer_id);
    } else if (m.stave_layer_encoding != "unique") {
      throw invalid_argument("Unknown layer_id encoding: " + m.stave_layer_encoding);
    }

    layerDE = DetElement(sdet,
                         _toString(layer_id, "layer_%d") +
                             _toString(int(count(layer_ids.begin(), layer_ids.end(), layer_id)), "_%d"),
                         layer_id);
    layerDE.setPlacement(pv);

    int nLadders = x_layer.attr<int>(_Unicode(nLadders));
    double dphi = 2. * M_PI / double(nLadders);
    double phi0 = x_layer.phi0(0.);

    //--------- loop over ladders ---------------------------
    for (int iStave = 0; iStave < nLadders; ++iStave) {
      double layer_r = x_layer.r() + ((iStave % 2 == 0) ? 0.0 : dr); // Offset every second stave in r
      double phi = phi0 + iStave * dphi;
      RotationZYX rot(phi, 0, 0);

      double r = layer_r;
      double r_offset_component = layer_offset;

      double x_pos = r * cos(phi) - r_offset_component * sin(phi);
      double y_pos = r * sin(phi) + r_offset_component * cos(phi);
      double z_pos = motherVolOffset;
      Position stave_pos = Position(x_pos, y_pos, z_pos);

      string stave_name = _toString(iStave, "stave%d");

      Volume whole_stave_volume_v;
      Assembly whole_stave_volume_a;
      if (m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0) {

        r = layer_r + m.motherVolThickness / 2.;
        x_pos = r * cos(phi) - r_offset_component * sin(phi);
        y_pos = r * sin(phi) + r_offset_component * cos(phi);
        Box whole_stave_box = Box(m.motherVolThickness / 2., m.motherVolWidth / 2., m.stave_length / 2.);
        whole_stave_volume_v = Volume(stave_name, whole_stave_box, theDetector.material("Air"));
        whole_stave_volume_v.setVisAttributes(theDetector, m.stave_vis);
        if (motherVolThickness > 0.)
          pv = whole_layer_volume_v.placeVolume(whole_stave_volume_v, Transform3D(rot, Position(x_pos, y_pos, z_pos)));
        else
          pv = whole_layer_volume_a.placeVolume(whole_stave_volume_v, Transform3D(rot, Position(x_pos, y_pos, z_pos)));
      } else {
        //// Use just assembly to avoid overlap between stave volume and other volumes
        whole_stave_volume_a = Assembly(stave_name);

        if (motherVolThickness > 0.)
          pv = whole_layer_volume_v.placeVolume(whole_stave_volume_a, Transform3D(rot, stave_pos));
        else
          pv = whole_layer_volume_a.placeVolume(whole_stave_volume_a, Transform3D(rot, stave_pos));
      }

      // Place components
      for (auto& component : m.components_vec) {
        Assembly component_assembly(stave_name + "_" + component.name);
        if (m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0)
          pv = whole_stave_volume_v.placeVolume(component_assembly, Position(-m.motherVolThickness / 2., 0., 0.));
        else
          pv = whole_stave_volume_a.placeVolume(component_assembly);

        for (int i = 0; i < int(component.thicknesses.size()); i++) {
          if (component.isCurved[i]) {
            double r_component_curved = 0.0; // component.r + component.rs[i]+ component.thicknesses[i]/2.;
            r_offset_component = component.offset + component.offsets[i] +
                                 (m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0 ? 0. : layer_offset);
            x_pos = r_component_curved * cos(phi) - r_offset_component * sin(phi);
            y_pos = r_component_curved * sin(phi) + r_offset_component * cos(phi);
            z_pos = component.z_offset + component.z_offsets[i] + motherVolOffset;
            Position pos(x_pos, y_pos, z_pos);
            component_assembly.placeVolume(component.volumes[i],
                                           Transform3D(rot, stave_pos).Inverse() * Transform3D(rot, pos));

            // Helper surface for this non-sensitive component volume
            SurfaceType type = SurfaceType::Helper;
            type.setProperty(SurfaceType::Cylinder, true);
            if (!component.notPlaceSurface[i])
              volSurfaceList(sdet)->push_back(VolCylinder(component.volumes[i], type, component.thicknesses[i] / 2.,
                                                          component.thicknesses[i] / 2., Vector3D(0., 0., 0.)));
          } else {
            x_pos = component.r + component.rs[i] + component.thicknesses[i] / 2.;
            y_pos = component.offset + component.offsets[i];
            z_pos = component.z_offset + component.z_offsets[i] + motherVolOffset;
            Position pos(x_pos, y_pos, z_pos);
            component_assembly.placeVolume(component.volumes[i],
                                           Transform3D(RotationZYX(component.phis[i], 0, 0), pos));

            if (!component.notPlaceSurface[i])
              volSurfaceList(sdet)->push_back(VolPlane(
                  component.volumes[i], SurfaceType::Helper, component.thicknesses[i] / 2.,
                  component.thicknesses[i] / 2., Vector3D(0., 1., 0.), Vector3D(0., 0., 1.), Vector3D(1., 0., 0.)));
          }
        }
        component_assembly->GetShape()->ComputeBBox();
      }

      // Place end of stave structures
      for (auto& endOfStave : m.endOfStaves_vec) {
        Assembly endOfStave_assembly(stave_name + "_" + endOfStave.name);
        if (m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0)
          pv = whole_stave_volume_v.placeVolume(endOfStave_assembly, Position(-m.motherVolThickness / 2., 0., 0.));
        else
          pv = whole_stave_volume_a.placeVolume(endOfStave_assembly);

        for (int i = 0; i < int(endOfStave.volumes.size()); i++) {
          std::vector<int> endOfStave_sides = {1, -1};
          if (endOfStave.side[i] != 0)
            endOfStave_sides = {endOfStave.side[i]};
          for (auto& endOfStave_side : endOfStave_sides) { // Place it on both sides of the stave
            if (endOfStave.isCurved[i]) {
              double r_component_curved =
                  0.0; // endOfStave.r + endOfStave.rs[i] + endOfStave.thicknesses[i]/2; // Correct for the fact that a
                       // tube element's origin is offset compared to the origin of a box
              r_offset_component = endOfStave.offset + endOfStave.offsets[i] +
                                   (m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0 ? 0. : layer_offset);
              x_pos = r_component_curved * cos(phi) - r_offset_component * sin(phi);
              y_pos = r_component_curved * sin(phi) + r_offset_component * cos(phi);
              z_pos = m.stave_length / 2. + endOfStave.lengths[i] / 2. + endOfStave.dzs[i] + endOfStave.z_offset +
                      endOfStave.z_offsets[i];
              Position pos = Position(x_pos, y_pos, z_pos * endOfStave_side + motherVolOffset);
              endOfStave_assembly.placeVolume(endOfStave.volumes[i],
                                              Transform3D(rot, stave_pos).Inverse() * Transform3D(rot, pos));

              SurfaceType type = SurfaceType::Helper;
              type.setProperty(SurfaceType::Cylinder, true);
              if (!endOfStave.notPlaceSurface[i])
                volSurfaceList(sdet)->push_back(VolCylinder(endOfStave.volumes[i], type, endOfStave.thicknesses[i] / 2.,
                                                            endOfStave.thicknesses[i] / 2., Vector3D(0., 0., 0.)));
            } else {
              x_pos = endOfStave.r + endOfStave.rs[i] + endOfStave.thicknesses[i] / 2.;
              y_pos = endOfStave.offset + endOfStave.offsets[i];
              z_pos = m.stave_length / 2. + endOfStave.lengths[i] / 2. + endOfStave.dzs[i] + endOfStave.z_offset +
                      endOfStave.z_offsets[i];
              Position pos(x_pos, y_pos, z_pos * endOfStave_side + motherVolOffset);
              endOfStave_assembly.placeVolume(endOfStave.volumes[i], pos);

              if (!endOfStave.notPlaceSurface[i])
                volSurfaceList(sdet)->push_back(VolPlane(
                    endOfStave.volumes[i], SurfaceType::Helper, endOfStave.thicknesses[i] / 2.,
                    endOfStave.thicknesses[i] / 2., Vector3D(0., 1., 0.), Vector3D(0., 0., 1.), Vector3D(1., 0., 0.)));
            }
          }
        }
        endOfStave_assembly->GetShape()->ComputeBBox();
      }

      // Place sensor
      bool multiType =
          m.sensorsVec.size() > 1; // Flag used to introducing sensor "type" if there are multiple sensors in .xml
      for (size_t iSensorType = 0; iSensorType < m.sensorsVec.size(); ++iSensorType) {
        auto& sensor = m.sensorsVec[iSensorType];
        auto logical_layer_id = 2 * m.sensorsVec.size() * layer_id + 2 * iSensorType + (iStave % 2);

        for (int iModule = 0; iModule < nmodules; iModule++) {
          x_pos = sensor.r + (iModule % 2 == 0 ? 0.0 : m.stave_dr);
          y_pos = sensor.offset;
          z_pos = motherVolOffset + -(nmodules - 1) / 2. * (sensor.length) - (nmodules - 1) / 2. * step_z +
                  iModule * sensor.length + iModule * step_z;
          Position pos(x_pos, y_pos, z_pos);

          string module_name = stave_name + _toString(iModuleTot, "_module%d");
          if (multiType)
            module_name += _toString(iSensorType, "_type%d"); // Append sensor type if multiple sensor types present
          Assembly module_assembly(module_name);
          if (m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0)
            pv = whole_stave_volume_v.placeVolume(module_assembly, Position(-m.motherVolThickness / 2., 0., 0.));
          else
            pv = whole_stave_volume_a.placeVolume(module_assembly);

          int iModule_VolID =
              int(iModuleTot / sensor.nGroupingModules); // To more efficiently use GlobalTrackerReadoutID bits, use the
                                                         // same module id for every nGroupingModules modules and
                                                         // distinguish them by the sensor id instead

          if (m.stave_layer_encoding == "unique") {
            pv.addPhysVolID("layer", logical_layer_id).addPhysVolID("side", 0);
            layer_ids.push_back(layer_id);
          }

          pv.addPhysVolID("module", iModule_VolID);
          DetElement moduleDE(layerDE, module_name, iModuleTot);
          moduleDE.setPlacement(pv);

          // Place all sensor parts
          int iSensitive = iModuleTot % sensor.nGroupingModules *
                           sensor.nSensitivePerModule; // To more efficiently use GlobalTrackerReadoutID bits, use the
                                                       // same module id for every nGroupingModules modules and
                                                       // distinguish them by the sensor id instead
          for (int x_i = 0; x_i < sensor.nx; x_i++) {  // Loop over duplicates in local x (global phi) direction
            for (int i = 0; i < int(sensor.volumes.size()); i++) {
              if (sensor.isCurved[i] && sensor.nsegments[i] == 1) { // Truly curved, part of tube
                r_offset_component = m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0 ? 0. : layer_offset;
                double phi_i =
                    phi + (sensor.xmin[i] + abs(sensor.xmax[i] - sensor.xmin[i]) / 2. + x_i * sensor.width) / m.stave_r;
                double r_component_curved = (iModule % 2 == 0 ? 0.0 : m.stave_dr);
                x_pos = r_component_curved * cos(phi_i);
                y_pos = r_component_curved * sin(phi_i);
                z_pos = motherVolOffset - (nmodules - 1) / 2. * (sensor.length) - (nmodules - 1) / 2. * step_z +
                        iModule * sensor.length + iModule * step_z;

                pos = Position(x_pos, y_pos, z_pos);
                Position pos2(0., 0., sensor.ymin[i] + abs(sensor.ymax[i] - sensor.ymin[i]) / 2.);
                RotationZYX rot2(phi_i, 0, 0);

                Position pos_offset(0., r_offset_component, 0.);
                pv = module_assembly.placeVolume(
                    sensor.volumes[i],
                    Translation3D(pos_offset) * Transform3D(rot, stave_pos).Inverse() * Transform3D(rot2, pos) *
                        Translation3D(pos2)); // To do: Simplify transformations, the orientation of Trd1 is strange...

                SurfaceType type = sensor.sensitives[i] ? SurfaceType::Sensitive : SurfaceType::Helper;
                type.setProperty(SurfaceType::Cylinder, true);
                Vector3D ocyl(-(r_component_curved + m.stave_r + sensor.r), 0., 0.);
                VolCylinder surf(sensor.volumes[i], type,
                                 sensor.thicknesses[i] / 2. + sensor.insensitive_thicknesses_below[i],
                                 sensor.thicknesses[i] / 2. + sensor.insensitive_thicknesses_above[i], ocyl);

                if (sensor.sensitives[i]) { // sensitive: keep DetElement and register surface on it
                  pv.addPhysVolID("sensor", iSensitive);
                  DetElement sensorDE(moduleDE, _toString(iSensitive, "_sensor%d"), iSensitive);
                  sensorDE.setPlacement(pv);
                  volSurfaceList(sensorDE)->push_back(surf);
                  iSensitive++;
                } else if (std::string(sensor.volumes[i].name()).find("above") == std::string::npos &&
                           std::string(sensor.volumes[i].name()).find("below") ==
                               std::string::npos) { // Don't add surfaces for 'below' or 'above' sensor elements, as
                                                    // they are already accounted for using
                                                    // other_insensitive_thickness_above/below
                  volSurfaceList(sdet)->push_back(surf);
                }
              } else if (sensor.isCurved[i]) { // Approximated quasi-curved sensors, using sensor.nsegments[i] tapezoids
                r_offset_component = m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0 ? 0. : layer_offset;
                for (int iSegment = 0; iSegment < sensor.nsegments[i]; iSegment++) {
                  double phi_i =
                      phi + (sensor.xmin[i] + abs(sensor.xmax[i] - sensor.xmin[i]) * iSegment / sensor.nsegments[i] +
                             abs(sensor.xmax[i] - sensor.xmin[i]) / 2. / sensor.nsegments[i] + x_i * sensor.width) /
                                m.stave_r;
                  double r_component_curved =
                      layer_r + sensor.thicknesses[i] / 2. + (iModule % 2 == 0 ? 0.0 : m.stave_dr) + sensor.rs[i];
                  x_pos = r_component_curved * cos(phi_i);
                  y_pos = r_component_curved * sin(phi_i);
                  z_pos = motherVolOffset - (nmodules - 1) / 2. * (sensor.length) - (nmodules - 1) / 2. * step_z +
                          iModule * sensor.length + iModule * step_z + sensor.ymin[i] +
                          abs(sensor.ymax[i] - sensor.ymin[i]) / 2.;

                  Position pos_offset(0., r_offset_component, 0.);

                  pv = module_assembly.placeVolume(
                      sensor.volumes[i],
                      Translation3D(pos_offset) * Transform3D(rot, stave_pos).Inverse() * RotationY(M_PI / 2.) *
                          Transform3D(
                              RotationZYX(M_PI / 2., 0., -phi_i),
                              Position(
                                  -z_pos, y_pos,
                                  x_pos))); // To do: Simplify transformations, the orientation of Trd1 is strange...

                  SurfaceType type = sensor.sensitives[i] ? SurfaceType::Sensitive : SurfaceType::Helper;
                  VolPlane surf(sensor.volumes[i], type,
                                sensor.thicknesses[i] / 2. + sensor.insensitive_thicknesses_below[i],
                                sensor.thicknesses[i] / 2. + sensor.insensitive_thicknesses_above[i],
                                Vector3D(1., 0., 0.), Vector3D(0., 1., 0.), Vector3D(0., 0., 1.));

                  if (sensor.sensitives[i]) {
                    pv.addPhysVolID("sensor", iSensitive);
                    DetElement sensorDE(moduleDE, _toString(iSensitive, "_sensor%d"), iSensitive);
                    sensorDE.setPlacement(pv);
                    volSurfaceList(sensorDE)->push_back(surf);
                    iSensitive++;
                  } else if (std::string(sensor.volumes[i].name()).find("above") == std::string::npos &&
                             std::string(sensor.volumes[i].name()).find("below") ==
                                 std::string::npos) { // Don't add surfaces for 'below' or 'above' sensor elements, as
                                                      // they are already accounted for using
                                                      // other_insensitive_thickness_above/below
                    volSurfaceList(sdet)->push_back(surf);
                  }
                }
              } else { // not curved (rectangular), use boxes
                x_pos = sensor.rs[i] + sensor.thicknesses[i] / 2.;
                y_pos = 0. + -(sensor.nx - 1) / 2. * sensor.width - (sensor.nx - 1) / 2. * step_rphi +
                        x_i * sensor.width + x_i * step_rphi;                        // Compute module centers in rphi
                y_pos += sensor.xmin[i] + abs(sensor.xmax[i] - sensor.xmin[i]) / 2.; // Place volume within nth module
                z_pos = sensor.ymin[i] + abs(sensor.ymax[i] - sensor.ymin[i]) / 2.;
                Position pos2(x_pos, y_pos, z_pos);

                pv = module_assembly.placeVolume(sensor.volumes[i], Translation3D(pos + pos2) * RotationY(M_PI / 2.) *
                                                                        RotationZ(M_PI / 2.));

                SurfaceType type = sensor.sensitives[i] ? SurfaceType::Sensitive : SurfaceType::Helper;
                VolPlane surf(sensor.volumes[i], type,
                              sensor.thicknesses[i] / 2. + sensor.insensitive_thicknesses_below[i],
                              sensor.thicknesses[i] / 2. + sensor.insensitive_thicknesses_above[i],
                              Vector3D(1., 0., 0.), Vector3D(0, 1., 0.), Vector3D(0., 0., 1.));

                if (sensor.sensitives[i]) { // sensitive: keep DetElement, register surface on it
                  pv.addPhysVolID("sensor", iSensitive);
                  DetElement sensorDE(moduleDE, _toString(iSensitive, "_sensor%d"), iSensitive);
                  sensorDE.setPlacement(pv);
                  volSurfaceList(sensorDE)->push_back(surf);
                  iSensitive++;
                } else if (std::string(sensor.volumes[i].name()).find("above") == std::string::npos &&
                           std::string(sensor.volumes[i].name()).find("below") ==
                               std::string::npos) { // Don't add surfaces for 'below' or 'above' sensor elements, as
                                                    // they are already accounted for using
                                                    // other_insensitive_thickness_above/below
                  volSurfaceList(sdet)->push_back(surf);
                }
              }
            }
            module_assembly->GetShape()->ComputeBBox();
          }
          iModuleTot++;
        }
      }

      if (m.motherVolThickness > 0.0 && m.motherVolWidth > 0.0) {
      } else {
        whole_stave_volume_a->GetShape()->ComputeBBox();
      }
    }

    // Save layer DetElement and maximal number of module IDs for later in case multiple layers with the same layer ID
    // are called
    iModuleTot_map[layer_id] += iModuleTot;
  }

  if (nLayers > 0) {
    pv.addPhysVolID("system", x_det.id());
    sdet.setAttributes(theDetector, envelope, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  } else {
    printout(WARNING, det_name,
             "No layer found for this detector. Did you ignore all layers using 'ignore' in layer definition?");
  }

  printout(INFO, det_name, "Building of detector successfully completed");
  return sdet;
}

DECLARE_DETELEMENT(VertexBarrel_detailed_o1_v04, create_element)

//----------------------------------
//         DCH detector v2
//----------------------------------

/*!
 *  \brief     Detector constructor of DCH v2
 *  \details   This code creates full geometry of IDEA DCH subdetector
 *  \author    Alvaro Tolosa-Delgado alvaro.tolosa.delgado@cern.ch
 *  \author    Brieuc Francois       brieuc.francois@cern.ch
 *  \version   2
 *  \date      2024
 *  \pre       DD4hep compiled with Geant4+Qt
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Shapes.h"
#include "XML/Utilities.h"

#include "DDRec/DCH_info.h"

namespace DCH_v2 {

using DCH_length_t = dd4hep::rec::DCH_info_struct::DCH_length_t;
using DCH_angle_t = dd4hep::rec::DCH_info_struct::DCH_angle_t;
using DCH_layer = dd4hep::rec::DCH_info_struct::DCH_layer;

dd4hep::Solid CompositeTT(double twist_angle, double cell_rin_z0, double cell_rout_z0, double dz, double dphi,
                          const dd4hep::rec::DCH_info& DCH_i);

/// Function to build DCH
static dd4hep::Ref_t create_DCH_o1_v02(dd4hep::Detector& desc, dd4hep::xml::Handle_t handle,
                                       dd4hep::SensitiveDetector sens) {
  dd4hep::xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  dd4hep::DetElement det(detName, detID);
  dd4hep::xml::setDetectorTypeFlag(detElem, det);
  sens.setType("tracker");

  // initialize empty DCH_info object
  // data extension mechanism requires it to be a raw pointer
  dd4hep::rec::DCH_info* DCH_i = new dd4hep::rec::DCH_info();
  // fill DCH_i with information from the XML file
  {
    // DCH outer geometry dimensions
    DCH_i->Set_rin(desc.constantAsDouble("DCH_gas_inner_cyl_R"));
    DCH_i->Set_rout(desc.constantAsDouble("DCH_gas_outer_cyl_R"));
    DCH_i->Set_lhalf(desc.constantAsDouble("DCH_gas_Lhalf"));

    // guard wires position, fix position
    DCH_i->Set_guard_rin_at_z0(desc.constantAsDouble("DCH_guard_inner_r_at_z0"));
    DCH_i->Set_guard_rout_at_zL2(desc.constantAsDouble("DCH_guard_outer_r_at_zL2"));

    DCH_angle_t dch_alpha = desc.constantAsDouble("DCH_alpha");
    DCH_i->Set_twist_angle(2 * dch_alpha);

    DCH_i->Set_nsuperlayers(desc.constantAsLong("DCH_nsuperlayers"));
    DCH_i->Set_nlayersPerSuperlayer(desc.constantAsLong("DCH_nlayersPerSuperlayer"));

    DCH_i->Set_ncell0(desc.constantAsLong("DCH_ncell"));
    DCH_i->Set_ncell_increment(desc.constantAsLong("DCH_ncell_increment"));
    DCH_i->Set_ncell_per_sector(desc.constantAsLong("DCH_ncell_per_sector"));

    DCH_i->Set_first_width(desc.constantAsDouble("DCH_first_width"));
    DCH_i->Set_first_sense_r(desc.constantAsDouble("DCH_first_sense_r"));

    bool buildLayers = detElem.attr<bool>(_Unicode(buildLayers));
    if (buildLayers) {
      DCH_i->BuildLayerDatabase();
      // safety check just in case something went wrong...
      if (DCH_i->IsDatabaseEmpty())
        throw std::runtime_error("Empty database");
    }

    bool printExcelTable = detElem.attr<bool>(_Unicode(printExcelTable));
    if (printExcelTable)
      DCH_i->Show_DCH_info_database(std::cout);
  }

  bool debugGeometry = detElem.hasChild(_Unicode(debugGeometry));
  bool useG4TT = detElem.hasChild(_Unicode(useG4TT));
  auto gasElem = detElem.child("gas");
  auto gasvolMat = desc.material(gasElem.attr<std::string>(_Unicode(material)));
  auto gasvolVis = desc.visAttributes(gasElem.attr<std::string>(_Unicode(vis)));

  auto vesselElem = detElem.child("vessel");
  auto vesselSkinVis = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(visSkin)));
  auto vesselBulkVis = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(visBulk)));
  auto vessel_mainMaterial = desc.material(vesselElem.attr<std::string>(_Unicode(mainMaterial)));
  auto vessel_fillmaterial_outerR = desc.material(vesselElem.attr<std::string>(_Unicode(fillmaterial_outerR)));
  auto vessel_fillmaterial_endcap = desc.material(vesselElem.attr<std::string>(_Unicode(fillmaterial_endcap)));
  DCH_length_t vessel_fillmaterial_fraction_outerR = vesselElem.attr<double>(_Unicode(fillmaterial_fraction_outerR));
  DCH_length_t vessel_fillmaterial_fraction_endcap = vesselElem.attr<double>(_Unicode(fillmaterial_fraction_endcap));

  if (0 > vessel_fillmaterial_fraction_outerR || 1 < vessel_fillmaterial_fraction_outerR)
    throw std::runtime_error("vessel_fillmaterial_fraction_outerR must be between 0 and 1");
  if (0 > vessel_fillmaterial_fraction_endcap || 1 < vessel_fillmaterial_fraction_endcap)
    throw std::runtime_error("vessel_fillmaterial_fraction_z must be between 0 and 1");

  auto wiresElem = detElem.child("wires");
  auto wiresVis = desc.visAttributes(wiresElem.attr<std::string>(_Unicode(vis)));
  bool buildSenseWires = wiresElem.attr<bool>(_Unicode(buildSenseWires));
  bool buildFieldWires = wiresElem.attr<bool>(_Unicode(buildFieldWires));

  DCH_length_t dch_SWire_thickness = wiresElem.attr<double>(_Unicode(SWire_thickness));
  DCH_length_t dch_FSideWire_thickness = wiresElem.attr<double>(_Unicode(FSideWire_thickness));
  DCH_length_t dch_FCentralWire_thickness = wiresElem.attr<double>(_Unicode(FCentralWire_thickness));

  auto dch_SWire_material = desc.material(wiresElem.attr<std::string>(_Unicode(SWire_material)));
  auto dch_FSideWire_material = desc.material(wiresElem.attr<std::string>(_Unicode(FSideWire_material)));
  auto dch_FCentralWire_material = desc.material(wiresElem.attr<std::string>(_Unicode(FCentralWire_material)));

  /* Geometry tree:
   * Gas (tube) -> Layer_1 (hyp) -> cell_1 (twisted tube)
   *                             -> cell_... (twisted tube)
   *            -> Layer_... (hyp) -> cell_1 (twisted tube)
   *                               -> cell_... (twisted tube)
   *            -> Inner radius vessel wall
   *            -> Outer radius vessel wall -> fill made of foam
   *            -> Endcap disk  vessel wall -> fill made of kapton
   *
   * Layers represent a segmentation in radius
   * Sectors represent a segmentation in phi
   * Each cell corresponds to a Detector Element
   * Vessel wall has to be defined as 3 volumes to account for independent thickness and materials
   */

  DCH_length_t safety_r_interspace = 1 * dd4hep::nm;
  DCH_length_t safety_z_interspace = 1 * dd4hep::nm;
  DCH_length_t safety_phi_interspace = 1e-6 * dd4hep::rad;

  DCH_length_t vessel_thickness_innerR = desc.constantAsDouble("DCH_vessel_thickness_innerR");
  DCH_length_t vessel_thickness_outerR = desc.constantAsDouble("DCH_vessel_thickness_outerR");
  DCH_length_t vessel_endcapdisk_zmin = desc.constantAsDouble("DCH_vessel_disk_zmin");
  DCH_length_t vessel_endcapdisk_zmax = desc.constantAsDouble("DCH_vessel_disk_zmax");

  // if( 0 > vessel_thickness_z )
  // throw std::runtime_error("vessel_thickness_z must be positive");
  if (0 > vessel_thickness_innerR)
    throw std::runtime_error("vessel_thickness_innerR must be positive");
  if (0 > vessel_thickness_outerR)
    throw std::runtime_error("vessel_thickness_outerR must be positive");

  // // // // // // // // // // // // // // //
  // // // // // MAIN VOLUME // // // // // //
  // // // // // // // // // // // // // // //
  dd4hep::Tube gas_s(DCH_i->rin - vessel_thickness_innerR, DCH_i->rout + vessel_thickness_outerR,
                     vessel_endcapdisk_zmax);
  dd4hep::Volume gas_v(detName + "_gas", gas_s, gasvolMat);
  gas_v.setVisAttributes(gasvolVis);
  gas_v.setRegion(desc, detElem.regionStr());
  gas_v.setLimitSet(desc, detElem.limitsStr());

  DCH_length_t vessel_innerR_start = DCH_i->rin - vessel_thickness_innerR + safety_r_interspace;
  DCH_length_t vessel_innerR_end = DCH_i->rin;
  DCH_length_t vessel_outerR_start = DCH_i->rout;
  DCH_length_t vessel_outerR_end = DCH_i->rout + vessel_thickness_outerR - safety_r_interspace;
  DCH_length_t vessel_R_zhalf = vessel_endcapdisk_zmax - safety_z_interspace;
  // // // // // // // // // // // // // // //
  // // // // //  INNER R WALL  // // // // //
  // // // // // // // // // // // // // // //
  dd4hep::Tube vessel_innerR_s(vessel_innerR_start, vessel_innerR_end, vessel_R_zhalf);
  dd4hep::Volume vessel_innerR_v(detName + "_vessel_innerR", vessel_innerR_s, vessel_mainMaterial);
  vessel_innerR_v.setVisAttributes(vesselSkinVis);
  gas_v.placeVolume(vessel_innerR_v);
  // // // // // // // // // // // // // // //
  // // // // //  OUTER R WALL  // // // // //
  // // // // // // // // // // // // // // //
  dd4hep::Tube vessel_outerR_s(vessel_outerR_start, vessel_outerR_end, vessel_R_zhalf);
  dd4hep::Volume vessel_outerR_v(detName + "_vessel_outerR", vessel_outerR_s, vessel_mainMaterial);
  vessel_outerR_v.setVisAttributes(vesselSkinVis);

  // if thickness fraction of bulk material is defined, build the bulk material
  if (0 < vessel_fillmaterial_fraction_outerR) {
    double f = vessel_fillmaterial_fraction_outerR;
    DCH_length_t fillmaterial_thickness = f * (vessel_outerR_end - vessel_outerR_start);
    DCH_length_t rstart = vessel_outerR_start + 0.5 * (1 - f) * fillmaterial_thickness;
    DCH_length_t rend = vessel_outerR_end - 0.5 * (1 - f) * fillmaterial_thickness;
    dd4hep::Tube vessel_fillmat_outerR_s(rstart, rend, vessel_R_zhalf - safety_z_interspace);
    dd4hep::Volume vessel_fillmat_outerR_v(detName + "_vessel_fillmat_outerR", vessel_fillmat_outerR_s,
                                           vessel_fillmaterial_outerR);
    vessel_fillmat_outerR_v.setVisAttributes(vesselBulkVis);
    vessel_outerR_v.placeVolume(vessel_fillmat_outerR_v);
  }
  gas_v.placeVolume(vessel_outerR_v);

  // // // // // // // // // // // // // // //
  // // // // //  ENDCAP WALL   // // // // //
  // // // // // // // // // // // // // // //
  DCH_length_t vessel_endcap_thickness = vessel_endcapdisk_zmax - vessel_endcapdisk_zmin - 2 * safety_z_interspace;
  DCH_length_t vessel_endcap_zpos = 0.5 * (vessel_endcapdisk_zmax + vessel_endcapdisk_zmin);
  DCH_length_t vessel_endcap_rstart = vessel_innerR_end + safety_r_interspace;
  DCH_length_t vessel_endcap_rend = vessel_outerR_start - safety_r_interspace;
  dd4hep::Tube vessel_endcap_s(vessel_endcap_rstart, vessel_endcap_rend, 0.5 * vessel_endcap_thickness);
  dd4hep::Volume vessel_endcap_v(detName + "_vessel_endcap", vessel_endcap_s, vessel_mainMaterial);
  vessel_endcap_v.setVisAttributes(vesselSkinVis);

  // if thickness fraction of bulk material is defined, build the bulk material
  if (0 < vessel_fillmaterial_fraction_endcap) {
    double f = vessel_fillmaterial_fraction_endcap;
    DCH_length_t fillmaterial_thickness = f * vessel_endcap_thickness;
    dd4hep::Tube vessel_fillmat_endcap_s(vessel_endcap_rstart + safety_r_interspace,
                                         vessel_endcap_rend - safety_r_interspace, 0.5 * fillmaterial_thickness);
    dd4hep::Volume vessel_fillmat_endcap_v(detName + "_vessel_fillmat_endcap", vessel_fillmat_endcap_s,
                                           vessel_fillmaterial_endcap);
    vessel_fillmat_endcap_v.setVisAttributes(vesselBulkVis);
    vessel_endcap_v.placeVolume(vessel_fillmat_endcap_v);
  }
  // place endcap wall at +/- z
  gas_v.placeVolume(vessel_endcap_v, dd4hep::Position(0, 0, vessel_endcap_zpos));
  gas_v.placeVolume(vessel_endcap_v, dd4hep::Position(0, 0, -vessel_endcap_zpos));

  // // // // // // // // // // // // // // //
  // // // // //  DCH layers    // // // // //
  // // // // // // // // // // // // // // //
  for (const auto& [ilayer, l] : DCH_i->database) {

    // // // // // // // // // // // // // // // // // // // // /
    // // // // // INITIALIZATION OF THE LAYER // // // // // //
    // // // // // // // // // // // // // // // // // // // //
    // Hyperboloid parameters:
    /// inner radius at z=0
    DCH_length_t rin = l.radius_fdw_z0 + safety_r_interspace;
    /// inner stereoangle, calculated from rin(z=0)
    DCH_angle_t stin = DCH_i->stereoangle_z0(rin);
    /// outer radius at z=0
    DCH_length_t rout = l.radius_fuw_z0 - safety_r_interspace;
    /// outer stereoangle, calculated from rout(z=0)
    DCH_angle_t stout = DCH_i->stereoangle_z0(rout);
    /// half-length
    DCH_length_t dz = DCH_i->Lhalf + safety_z_interspace;

    dd4hep::Hyperboloid layer_s(rin, stin, rout, stout, dz);

    std::string layer_name = detName + "_layer" + std::to_string(ilayer);
    dd4hep::Volume layer_v(layer_name, layer_s, gasvolMat);
    layer_v.setVisAttributes(desc.visAttributes(Form("dch_layer_vis%d", ilayer % 22)));
    auto layer_pv = gas_v.placeVolume(layer_v);
    // ilayer is a counter that runs from 1 to 112 (nsuperlayers * nlayersPerSuperlayer)
    // it seems more convenient to store the layer number within the superlayer
    // ilayerWithinSuperlayer runs from 0 to 7 (nlayersPerSuperlayer-1)
    int ilayerWithinSuperlayer = (ilayer - 1) % DCH_i->nlayersPerSuperlayer;
    layer_pv.addPhysVolID("layer", ilayerWithinSuperlayer);
    // add superlayer bitfield
    int nsuperlayer_minus_1 = DCH_i->Get_nsuperlayer_minus_1(ilayer);
    layer_pv.addPhysVolID("superlayer", nsuperlayer_minus_1);

    dd4hep::DetElement layer_DE(det, layer_name + "DE", ilayer);
    layer_DE.setPlacement(layer_pv);

    // // // // // // // // // // // // // // // // // // // //
    // // // // // SEGMENTATION OF THE LAYER  // // // // // //
    // // // // // INTO CELLS (TWISTED TUBES) // // // // // //
    // // // // // // // // // // // // // // // // // // // //

    // ncells in this layer = 2x number of wires
    int ncells = l.nwires / 2;
    DCH_angle_t phi_step = (TMath::TwoPi() / ncells) * dd4hep::rad;

    // unitary cell (Twisted tube) is repeated for each layer l.nwires/2 times
    // Twisted tube parameters
    DCH_angle_t cell_twistangle = l.StereoSign() * DCH_i->twist_angle;
    DCH_length_t cell_rin_z0 = l.radius_fdw_z0 + 2 * safety_r_interspace;
    DCH_length_t cell_rout_z0 = l.radius_fuw_z0 - 2 * safety_r_interspace;
    DCH_length_t cell_rin_zLhalf = DCH_i->Radius_zLhalf(cell_rin_z0);
    DCH_length_t cell_rout_zLhalf = DCH_i->Radius_zLhalf(cell_rout_z0);
    DCH_length_t cell_dz = DCH_i->Lhalf;
    DCH_angle_t cell_phi_width = phi_step - safety_phi_interspace;
    dd4hep::Solid cell_s;
    if (useG4TT)
      cell_s = dd4hep::TwistedTube(cell_twistangle, cell_rin_zLhalf, cell_rout_zLhalf, cell_dz, 1, cell_phi_width);
    else
      cell_s = CompositeTT(cell_twistangle, cell_rin_z0, cell_rout_z0, cell_dz, cell_phi_width, *DCH_i);

    // initialize cell volume
    std::string cell_name = detName + "_layer" + std::to_string(ilayer) + "_cell";
    dd4hep::Volume cell_v(cell_name, cell_s, gasvolMat);
    cell_v.setSensitiveDetector(sens);
    cell_v.setVisAttributes(desc.visAttributes("dch_no_vis_nodaughters"));

    // // // // // // // // // // // // // // // // // // // //
    // // // // // // POSITIONING OF WIRES // // // // // // //
    // // // // // // // // // // // // // // // // // // // //
    {
      // // // // // // // // // // // // // // // // // // // //
      // // // // // // POSITIONING OF SENSE WIRES // // // // //
      // // // // // // // // // // // // // // // // // // // //
      // average radius to position sense wire
      DCH_length_t cell_rave_z0 = 0.5 * (cell_rin_z0 + cell_rout_z0);
      DCH_length_t cell_swire_radius = dch_SWire_thickness / 2;
      DCH_length_t swlength = 0.5 * DCH_i->WireLength(ilayer, cell_rave_z0) -
                              cell_swire_radius * cos(DCH_i->stereoangle_z0(cell_rave_z0)) - safety_z_interspace;
      if (buildSenseWires) {
        dd4hep::Tube swire_s(0., dch_SWire_thickness, swlength);
        dd4hep::Volume swire_v(cell_name + "_swire", swire_s, dch_SWire_material);
        swire_v.setVisAttributes(wiresVis);
        // Change sign of stereo angle to place properly the wire inside the twisted tube
        dd4hep::RotationX stereoTr((-1.) * l.StereoSign() * DCH_i->stereoangle_z0(cell_rave_z0));
        dd4hep::Transform3D swireTr(stereoTr * dd4hep::Translation3D(cell_rave_z0, 0., 0.));
        cell_v.placeVolume(swire_v, swireTr);
      }
      if (buildFieldWires) {
        // // // // // // // // // // // // // // // // // // // //
        // // // // // // POSITIONING OF FIELD WIRES // // // // //
        // // // // // // // // // // // // // // // // // // // //
        //
        //  The following sketch represents the crossection of a DCH cell, where
        //      O symbol = Field wires, the number in parenthesis is used as ID
        //      X symbol = sense wire
        //
        //   ^ radius
        //
        //   O(1)---O(4)---O(6)    radius_z0 = l.radius_fuw_z0 == (++l).radius_fdw_z0
        //
        //   O(2)   X      O(7)    radius_z0 = average(l.radius_fuw_z0, l.radius_fdw_z0)
        //
        //   O(3)---O(5)---O(8)    radius_z0 = l.radius_fdw_z0 == (--l).radius_fuw_z0
        //
        //   --> phi axis
        //
        //  In the previous sketch, the wires are shared among several cells.
        //  Since we are using an actual shape to contain each cell,
        //  it is not feasible.
        //
        //  As a workaround, we introduce an offset in phi and radially to the cell center,
        //  in such a manner that the wires are fully contained in one cell.
        //  The following code implements the following sketch:
        //
        //   O(1)---O(4)---    radius_z0 = l.radius_fuw_z0 - wire_thickness/2
        //
        //   O(2)   X          radius_z0 = average(l.radius_fuw_z0, l.radius_fdw_z0)
        //
        //   O(3)---O(5)---    radius_z0 = l.radius_fdw_z0 + wire_thickness/2
        //
        //  phi_offset(n) = atan(  wire_thickness/2 / radius_z0 )
        //
        //  notice that the field wires are offcentered with respect to the sense wire
        //  by about 20um/1cm ~ 0.1 mrad, which is not expected to have any impact

        /// encapsulate the calculation of the phi offset into a function
        /// since it will be different for each field wire
        /// it includes the safety phi distance
        auto fwire_phi_offset = [&](DCH_length_t radial_distance, DCH_length_t wire_radius) -> DCH_angle_t {
          return atan(wire_radius / radial_distance) * dd4hep::rad + safety_phi_interspace;
        };

        // // // // // // // // // // // // // // // // // // // //
        // // // // // // POSITIONING OF F WIRE 2 // // // // // //
        // // // // // // REQUIRES OFFSET OF PHI  // // // // // //
        // // // // // // // // // // // // // // // // // // // //
        {
          DCH_length_t fwire_radius = dch_FCentralWire_thickness / 2;
          DCH_length_t fwire_r_z0 = cell_rave_z0;
          DCH_angle_t fwire_stereo = (-1.) * l.StereoSign() * DCH_i->stereoangle_z0(fwire_r_z0);
          DCH_angle_t fwire_phi = -cell_phi_width / 2 + fwire_phi_offset(fwire_r_z0, fwire_radius);
          DCH_length_t fwire_length = 0.5 * DCH_i->WireLength(ilayer, fwire_r_z0) -
                                      fwire_radius * cos(DCH_i->stereoangle_z0(fwire_r_z0)) - safety_z_interspace;

          dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
          dd4hep::Volume fwire_v(cell_name + "_f2wire", fwire_s, dch_FCentralWire_material);
          fwire_v.setVisAttributes(wiresVis);
          // Change sign of stereo angle to place properly the wire inside the twisted tube
          dd4hep::RotationX fwireStereoTr(fwire_stereo);
          dd4hep::RotationZ fwirePhoTr(fwire_phi);
          dd4hep::Transform3D fwireTr(fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0, 0., 0.));
          cell_v.placeVolume(fwire_v, fwireTr);
        }

        // // // // // // // // // // // // // // // // // // // //
        // // // // // // POSITIONING OF F WIRE 1    // // // // //
        // // // // // // REQUIRES OFFSET OF PHI & R // // // // //
        // // // // // // // // // // // // // // // // // // // //
        {
          DCH_length_t fwire_radius = dch_FSideWire_thickness / 2;
          // decrease radial distance, move it closer to the sense wire
          DCH_length_t fwire_r_z0 = cell_rout_z0 - fwire_radius;
          DCH_angle_t fwire_stereo = (-1.) * l.StereoSign() * DCH_i->stereoangle_z0(fwire_r_z0);
          DCH_angle_t fwire_phi = -cell_phi_width / 2 + fwire_phi_offset(fwire_r_z0, fwire_radius);
          DCH_length_t fwire_length = 0.5 * DCH_i->WireLength(ilayer, fwire_r_z0) -
                                      fwire_radius * cos(DCH_i->stereoangle_z0(fwire_r_z0)) - safety_z_interspace;

          dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
          dd4hep::Volume fwire_v(cell_name + "_f1wire", fwire_s, dch_FSideWire_material);
          fwire_v.setVisAttributes(wiresVis);
          // Change sign of stereo angle to place properly the wire inside the twisted tube
          dd4hep::RotationX fwireStereoTr(fwire_stereo);
          dd4hep::RotationZ fwirePhoTr(fwire_phi);
          dd4hep::Transform3D fwireTr(fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0, 0., 0.));
          cell_v.placeVolume(fwire_v, fwireTr);
        }
        // // // // // // // // // // // // // // // // // // // //
        // // // // // // POSITIONING OF F WIRE 3    // // // // //
        // // // // // // REQUIRES OFFSET OF PHI & R // // // // //
        // // // // // // // // // // // // // // // // // // // //
        {
          DCH_length_t fwire_radius = dch_FSideWire_thickness / 2;
          // increase radial distance, move it closer to the sense wire
          DCH_length_t fwire_r_z0 = cell_rin_z0 + fwire_radius;
          DCH_angle_t fwire_stereo = (-1.) * l.StereoSign() * DCH_i->stereoangle_z0(fwire_r_z0);
          DCH_angle_t fwire_phi = -cell_phi_width / 2 + fwire_phi_offset(fwire_r_z0, fwire_radius);
          DCH_length_t fwire_length = 0.5 * DCH_i->WireLength(ilayer, fwire_r_z0) -
                                      fwire_radius * cos(DCH_i->stereoangle_z0(fwire_r_z0)) - safety_z_interspace;

          dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
          dd4hep::Volume fwire_v(cell_name + "_f3wire", fwire_s, dch_FSideWire_material);
          fwire_v.setVisAttributes(wiresVis);
          // Change sign of stereo angle to place properly the wire inside the twisted tube
          dd4hep::RotationX fwireStereoTr(fwire_stereo);
          dd4hep::RotationZ fwirePhoTr(fwire_phi);
          dd4hep::Transform3D fwireTr(fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0, 0., 0.));
          cell_v.placeVolume(fwire_v, fwireTr);
        }
        // // // // // // // // // // // // // // // // // // // //
        // // // // // // POSITIONING OF F WIRE 5    // // // // //
        // // // // // // REQUIRES OFFSET OF R       // // // // //
        // // // // // // // // // // // // // // // // // // // //
        {
          DCH_length_t fwire_radius = dch_FSideWire_thickness / 2;
          // increase radial distance, move it closer to the sense wire
          DCH_length_t fwire_r_z0 = cell_rin_z0 + fwire_radius;
          DCH_angle_t fwire_stereo = (-1.) * l.StereoSign() * DCH_i->stereoangle_z0(fwire_r_z0);
          DCH_length_t fwire_length = 0.5 * DCH_i->WireLength(ilayer, fwire_r_z0) -
                                      fwire_radius * cos(DCH_i->stereoangle_z0(fwire_r_z0)) - safety_z_interspace;

          dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
          dd4hep::Volume fwire_v(cell_name + "_f5wire", fwire_s, dch_FSideWire_material);
          fwire_v.setVisAttributes(wiresVis);
          // Change sign of stereo angle to place properly the wire inside the twisted tube
          dd4hep::RotationX fwireStereoTr(fwire_stereo);
          dd4hep::Transform3D fwireTr(fwireStereoTr * dd4hep::Translation3D(fwire_r_z0, 0., 0.));
          cell_v.placeVolume(fwire_v, fwireTr);
        }
        // // // // // // // // // // // // // // // // // // // //
        // // // // // // POSITIONING OF F WIRE 4    // // // // //
        // // // // // // REQUIRES OFFSET OF R       // // // // //
        // // // // // // // // // // // // // // // // // // // //
        {
          DCH_length_t fwire_radius = dch_FSideWire_thickness / 2;
          // increase radial distance, move it closer to the sense wire
          DCH_length_t fwire_r_z0 = cell_rout_z0 - fwire_radius;
          DCH_angle_t fwire_stereo = (-1.) * l.StereoSign() * DCH_i->stereoangle_z0(fwire_r_z0);
          DCH_length_t fwire_length = 0.5 * DCH_i->WireLength(ilayer, fwire_r_z0) -
                                      fwire_radius * cos(DCH_i->stereoangle_z0(fwire_r_z0)) - safety_z_interspace;

          dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
          dd4hep::Volume fwire_v(cell_name + "_f4wire", fwire_s, dch_FSideWire_material);
          fwire_v.setVisAttributes(wiresVis);
          // Change sign of stereo angle to place properly the wire inside the twisted tube
          dd4hep::RotationX fwireStereoTr(fwire_stereo);
          dd4hep::Transform3D fwireTr(fwireStereoTr * dd4hep::Translation3D(fwire_r_z0, 0., 0.));
          cell_v.placeVolume(fwire_v, fwireTr);
        }
      } // end building field wires
    }   /// end building wires
    int maxphi = ncells;
    if (debugGeometry)
      maxphi = 3;
    for (int nphi = 0; nphi < maxphi; ++nphi) {
      // TODO: check if staggering is just + 0.25*cell_phi_width*(ilayer%2);
      // phi positioning, adding offset for odd ilayers
      DCH_angle_t cell_phi_angle = phi_step * nphi + 0.25 * cell_phi_width * (ilayer % 2);
      // conversion of RotationZ into Transform3D using constructor)
      dd4hep::Transform3D cellTr{dd4hep::RotationZ(cell_phi_angle)};
      auto cell_pv = layer_v.placeVolume(cell_v, cellTr);
      cell_pv.addPhysVolID("nphi", nphi);
      cell_pv.addPhysVolID("stereosign", l.StereoSign());

      dd4hep::DetElement cell_DE(layer_DE, cell_name + std::to_string(nphi) + "DE", nphi);
      cell_DE.setPlacement(cell_pv);
    }
  }

  // Place our mother volume in the world
  dd4hep::Volume wVol = desc.pickMotherVolume(det);
  dd4hep::PlacedVolume vessel_pv = wVol.placeVolume(gas_v);
  // Associate the wall to the detector element.
  det.setPlacement(vessel_pv);
  // Assign the system ID to our mother volume
  vessel_pv.addPhysVolID("system", detID);
  det.addExtension<dd4hep::rec::DCH_info>(DCH_i);

  return det;
}

inline double Circumradius(double Apothem, double dphi) { return Apothem / cos(0.5 * dphi / dd4hep::rad); }

/// Solid equivalent to a twisted tube, resulting from the intersection of an hyperboloid and a generic trapezoid
/// the hyperboloid provides the hyperboloidal surfaces, the trapezoid provides the other two types of surfaces
/// the generic trapezoid is built in such a manner that circumscribe the twisted tube
dd4hep::Solid CompositeTT(double twist_angle, double cell_rin_z0, double cell_rout_z0, double dz, double dphi,
                          const dd4hep::rec::DCH_info& DCH_i) {

  //----------------- G. trapezoid -------------
  // make generic trapezoid bigger, later intersected with hyperboloid of proper radii
  double rmin_zLhalf = DCH_i.Radius_zLhalf(cell_rin_z0);
  double rout_zLhalf = DCH_i.Radius_zLhalf(cell_rout_z0);
  double trap_rin = 0.9 * rmin_zLhalf;
  double trap_rout = Circumradius(1.1 * rout_zLhalf, dphi);

  double poly_angle = dphi / 2;
  double twist_angle_half = twist_angle / 2.;
  // change sign, so the final shape has the same orientation as G4 twisted tube
  twist_angle_half *= -1;

  // define points of 8 genenric trapezoid
  struct point2d {
    double x = {0.0};
    double y = {0.0};
  };

  struct face {
    point2d A;
    point2d B;
    point2d C;
    point2d D;
  };

  face fZneg;
  face fZpos;

  // The generic trapezoid is built in such a manner as to circumscribe the twisted tube
  // The following points correspond to the corners of the twisted tube

  fZpos.A = {trap_rin * cos(poly_angle + twist_angle_half), trap_rin * sin(poly_angle + twist_angle_half)};
  fZpos.B = {trap_rin * cos(-poly_angle + twist_angle_half), trap_rin * sin(-poly_angle + twist_angle_half)};

  fZneg.A = {trap_rin * cos(poly_angle - twist_angle_half), trap_rin * sin(poly_angle - twist_angle_half)};
  fZneg.B = {trap_rin * cos(-poly_angle - twist_angle_half), trap_rin * sin(-poly_angle - twist_angle_half)};

  fZpos.C = {trap_rout * cos(poly_angle + twist_angle_half), trap_rout * sin(poly_angle + twist_angle_half)};
  fZpos.D = {trap_rout * cos(-poly_angle + twist_angle_half), trap_rout * sin(-poly_angle + twist_angle_half)};

  fZneg.C = {trap_rout * cos(poly_angle - twist_angle_half), trap_rout * sin(poly_angle - twist_angle_half)};
  fZneg.D = {trap_rout * cos(-poly_angle - twist_angle_half), trap_rout * sin(-poly_angle - twist_angle_half)};

  std::vector<double> vertices_array = {fZpos.B.x, fZpos.B.y, fZpos.A.x, fZpos.A.y, fZpos.C.x, fZpos.C.y,
                                        fZpos.D.x, fZpos.D.y, fZneg.B.x, fZneg.B.y, fZneg.A.x, fZneg.A.y,
                                        fZneg.C.x, fZneg.C.y, fZneg.D.x, fZneg.D.y};

  dd4hep::EightPointSolid gtrap_shape(dz, vertices_array.data());

  //----------------- Hyperboloid -------------
  // ROOT hyperboloid require stereoangles stin and stout
  //-- stereo for rmin
  double stin = DCH_i.stereoangle_z0(cell_rin_z0);
  //-- stereo for rout
  double stout = DCH_i.stereoangle_z0(cell_rout_z0);

  // make hyperboloid longer, later intersection with gtrap will lead to proper length
  double dz_safe = 1.1 * dz;
  dd4hep::Hyperboloid layer_s(cell_rin_z0, stin, cell_rout_z0, stout, dz_safe);

  // create a twisted tube as intersection of the generic trapezoid and the hyperboloid
  dd4hep::Solid mytt = dd4hep::IntersectionSolid(gtrap_shape, layer_s);
  return mytt;
}

}; // namespace DCH_v2

DECLARE_DETELEMENT(DriftChamber_o1_v02_T, DCH_v2::create_DCH_o1_v02)

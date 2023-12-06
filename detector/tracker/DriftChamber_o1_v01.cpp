/***************************************************************************************\
* DD4hep geometry code for the central drift chamber of the IDEA detector               *
* Author: Lorenzo Capriotti, Modified by Brieuc Francois to have sensitive cell volumes *
\***************************************************************************************/
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/detail/DetectorInterna.h"
#include "TClass.h"
#include "TMath.h"
#include "XML/Utilities.h"
#include <XML/Layering.h>
#include <iostream>
#include <math.h>

using namespace std;
using namespace dd4hep;

struct wire
{
  dd4hep::Volume mother_volume;
  string type; // F = field wire, S = sense wire, G = guard wire
  int num; // how many wires for the layer being considered
  double radius; // radial position of the wire at z = 0
  double delta_phi; // delta phi between two wires of this type
  double phioffset; // phi angle of the first wire of this type
  double stereo; // stereo angle
  double thickness; // thickness of the core of the wire (no coating)
  double halflength; // wire z extent before stereo angle rotation
  dd4hep::Volume volume; // actual volume of the wire
  string name;
};

namespace {

struct CDCHBuild : public dd4hep::xml::tools::VolumeBuilder {
  std::vector<dd4hep::DetElement> deSuperLayer, deLayer, deSWire;

  CDCHBuild(dd4hep::Detector& description, xml_elt_t e, dd4hep::SensitiveDetector sens);

  double diff_of_squares(double a, double b);
  void apply_wire_coating(struct wire& w, double outwrap, double halflength, string material);
  void PlaceGuardWires(struct wire& w, double outwrap, double halflength, int SL, int ilayer);
  void build_layer(DetElement parent, Volume parentVol, dd4hep::SensitiveDetector sens);
};

// ******************************************************
// Initializing constructor
// ******************************************************

CDCHBuild::CDCHBuild(dd4hep::Detector& dsc, xml_elt_t e, dd4hep::SensitiveDetector sens)
    : dd4hep::xml::tools::VolumeBuilder(dsc, e, sens) {}

double CDCHBuild::diff_of_squares(double a, double b) {

  double diff = pow(a, 2) - pow(b, 2);
  return diff;
}

void CDCHBuild::apply_wire_coating(struct wire& w, double outwrap, double halflength, string material = "Silver"){
  dd4hep::Tube WrapTube(w.thickness, w.thickness + 0.5 * outwrap, halflength);
  dd4hep::Volume WireWrapVol(w.name + "_coating", WrapTube, description.material(material));
  dd4hep::Tube TotalWire(0.0, w.thickness + 0.5 * outwrap, halflength);
  dd4hep::Volume WireVol(w.name + "_totalWire", TotalWire, description.material("Air"));
  WireVol.placeVolume(w.volume, dd4hep::Position(0.0, 0.0, 0.0));
  WireVol.placeVolume(WireWrapVol, dd4hep::Position(0.0, 0.0, 0.0));
  w.volume = WireVol;
}

void CDCHBuild::PlaceGuardWires(struct wire& w, double outwrap, double halflength, int SL = 999,
                           int ilayer = 999) {

  dd4hep::RotationZYX rot(0., 0., w.stereo);
  dd4hep::RotationX rot_stereo(w.stereo);
  dd4hep::Translation3D transl(w.radius, 0., 0.);

  dd4hep::Transform3D T(transl * rot_stereo);

  string wirewrapname = "WireWrap_SL";
  wirewrapname += std::to_string(SL);
  wirewrapname += "_layer";
  wirewrapname += std::to_string(ilayer);
  wirewrapname += "_type";
  wirewrapname += w.type;
  wirewrapname += "_stereo";
  wirewrapname += std::to_string(w.stereo);

  string wirename = "Wire_SL";
  wirename += std::to_string(SL);
  wirename += "_layer";
  wirename += std::to_string(ilayer);
  wirename += "_type";
  wirename += w.type;
  wirename += "_stereo";
  wirename += std::to_string(w.stereo);

  apply_wire_coating(w, outwrap, halflength);

  // repeat the placement of wires over phi
  for (int n = 0; n < w.num; n++) {
    dd4hep::RotationZ iRot(w.phioffset + w.delta_phi * n);
    if (n % 1 == 0) w.mother_volume.placeVolume(w.volume, dd4hep::Transform3D(iRot * T));
  }
}

void CDCHBuild::build_layer(DetElement parent, Volume parentVol, dd4hep::SensitiveDetector sens_det) {

  // ******************************************************
  // Loading parameters
  // ******************************************************

  double halfalpha = 0.5 * dd4hep::_toDouble("CDCH:alpha");
  double inner_radius = dd4hep::_toDouble("CDCH:inner_radius");
  double outer_radius = dd4hep::_toDouble("CDCH:outer_radius");
  double zHalfExtentWithServices = dd4hep::_toDouble("CDCH:zHalfExtentWithServices");
  double CarbonInnerWallThick = dd4hep::_toDouble("CDCH:CarbonInnerWallThick");
  double CopperInnerWallThick = dd4hep::_toDouble("CDCH:CopperInnerWallThick");
  double GasInnerWallThick = dd4hep::_toDouble("CDCH:GasInnerWallThick");
  double Carbon1OuterWallThick = dd4hep::_toDouble("CDCH:Carbon1OuterWallThick");
  double Carbon2OuterWallThick = dd4hep::_toDouble("CDCH:Carbon2OuterWallThick");
  double CopperOuterWallThick = dd4hep::_toDouble("CDCH:CopperOuterWallThick");
  double FoamOuterWallThick = dd4hep::_toDouble("CDCH:FoamOuterWallThick");
  double GasEndcapWallThick = dd4hep::_toDouble("CDCH:GasEndcapWallThick");
  double CopperEndcapWallThick = dd4hep::_toDouble("CDCH:CopperEndcapWallThick");
  double KaptonEndcapWallThick = dd4hep::_toDouble("CDCH:KaptonEndcapWallThick");
  double CarbonEndcapWallThick = dd4hep::_toDouble("CDCH:CarbonEndcapWallThick");
  double FWireShellThickIn = dd4hep::_toDouble("CDCH:FWireShellThickIn");
  double FWireShellThickOut = dd4hep::_toDouble("CDCH:FWireShellThickOut");
  double centerFWireShellThickIn = dd4hep::_toDouble("CDCH:centerFWireShellThickIn");
  double centerFWireShellThickOut = dd4hep::_toDouble("CDCH:centerFWireShellThickOut");
  double SWireShellThickIn = dd4hep::_toDouble("CDCH:SWireShellThickIn");
  double SWireShellThickOut = dd4hep::_toDouble("CDCH:SWireShellThickOut");
  double InGWireShellThickIn = dd4hep::_toDouble("CDCH:InGWireShellThickIn");
  double InGWireShellThickOut = dd4hep::_toDouble("CDCH:InGWireShellThickOut");
  double OutGWireShellThickIn = dd4hep::_toDouble("CDCH:OutGWireShellThickIn");
  double OutGWireShellThickOut = dd4hep::_toDouble("CDCH:OutGWireShellThickOut");
  //------------------------------------------------------------------------
  // The wireThicknessDilution parameter is used to see the wires in the rendering. NB: not propagated to the envelope volumes --> you will have volume extrusions if != 1
  //------------------------------------------------------------------------
  double wireThicknessDilution = dd4hep::_toDouble("CDCH:wireThicknessDilution");
  double secure = dd4hep::_toDouble("CDCH:secure");
  double capGasLayer = dd4hep::_toDouble("CDCH:capGasLayer");
  double extShiftFW = dd4hep::_toDouble("CDCH:extShiftFW");
  double cellDimension = dd4hep::_toDouble("CDCH:cellDimension");
  double inGuardRad = dd4hep::_toDouble("CDCH:inGuardRad");
  double outGuardRad = dd4hep::_toDouble("CDCH:outGuardRad");
  int nSDeltaWire = dd4hep::_toInt("CDCH:nSDeltaWire");
  int nSWireFirstLayer = dd4hep::_toInt("CDCH:nSWireFirstLayer");
  int nFtoSWireRatio = dd4hep::_toInt("CDCH:nFtoSWireRatio");
  int nCenterFWirePerCell = dd4hep::_toInt("CDCH:nCenterFWirePerCell");
  int nSuperLayer = dd4hep::_toInt("CDCH:nSuperLayer");
  int nLayer = dd4hep::_toInt("CDCH:nLayer");
  //bool setWireSensitive = true; // FIXME: add the possibility to have wires sensitive (parameter in the xml) which could be useful for detailed chamber behavior studies
  double halflength = zHalfExtentWithServices - (GasEndcapWallThick + CopperEndcapWallThick + KaptonEndcapWallThick + CarbonEndcapWallThick); // this will be the sensitive volume z extent

  double epsilon = 0.0;
  double phi_offset = 0.0;
  double delta_phi_top_bottom = 0.0;
  int nFWireTopAndBottom = 0;
  int nFWireTopOrBottom = 0;
  int nSWire = 0;
  int nHorizontalFtoSWireRatio = nFtoSWireRatio - nCenterFWirePerCell; // nHorizontalFtoSWireRatio is, per cell, how many wires do we have in top and bottom layer (summed), knowing that the rightmost wires are not in the current cell
  int sign_epsilon = -1;
  double sense_wire_delta_phi = 0.0;
  double scaleFactor = 0.0;
  double dropFactor = 0.0;
  double epsilonFactor = 0.0;
  double delta_radius_layer = cellDimension;
  double sense_wire_radius_at_z_0 = 0.0;
  double iradius = 0.0;
  double idelta_radius = 0.0;

  double envelop_Inner_thickness = CarbonInnerWallThick + CopperInnerWallThick + GasInnerWallThick;
  double envelop_Outer_thickness = Carbon1OuterWallThick + Carbon2OuterWallThick + CopperOuterWallThick + FoamOuterWallThick;
  double FWireDiameter = FWireShellThickIn + FWireShellThickOut;
  double FWradii = 0.5 * FWireDiameter;
  double inGWireDiameter = InGWireShellThickIn + InGWireShellThickOut;
  double inGWradii = 0.5 * inGWireDiameter;
  double fakeLayerInIWthick = -0.0001 + GasInnerWallThick;
  double guard_layer_inner_radius_at_z_0 = inner_radius + envelop_Inner_thickness - fakeLayerInIWthick;

  double radius_layer_0 = inner_radius + envelop_Inner_thickness + FWradii + secure + capGasLayer;
  double layer_outer_radius_at_z_0 = radius_layer_0 - FWradii - secure;
  double layer_inner_radius_at_z_0 = 0.0;
  double layer_inner_radius_at_z_end = 0.0;
  double layer_outer_radius_at_z_end = 0.0;

  double drop = 0.0;
  double radius_layer = 0.0;
  double epsilonIn = 0.0;
  double epsilonOut = 0.0;
  double phi_offset_cell_start = 0.0;
  double zlength = 0.0;
  double cellStaggering = 0.0;
  double epsilonInGwRing = 0.0;
  double epsilonOutGwRing = 0.0;
  double layer_inner_radius_at_z_end_whole_cell = 0.0;
  double epsilonIn_whole_cell = 0.0;

  //------------------------------------------------------------------------
  // Build the inner, outer and endcap walls first
  //------------------------------------------------------------------------

  dd4hep::Tube Endcap_Gas(inner_radius, outer_radius, 0.5 * GasEndcapWallThick);
  dd4hep::Tube Endcap_Copper(inner_radius, outer_radius, 0.5 * CopperEndcapWallThick);
  dd4hep::Tube Endcap_Kapton(inner_radius, outer_radius, 0.5 * KaptonEndcapWallThick);
  dd4hep::Tube Endcap_Carbon(inner_radius, outer_radius, 0.5 * CarbonEndcapWallThick);

  dd4hep::Volume lvEndcapWallGas =
      dd4hep::Volume("lvEndcapWallGasVol", Endcap_Gas, description.material("GasHe_90Isob_10"));
  dd4hep::Volume lvEndcapWallCopper =
      dd4hep::Volume("lvEndcapWallCopperVol", Endcap_Copper, description.material("G4_Cu"));
  dd4hep::Volume lvEndcapWallKapton =
      dd4hep::Volume("lvEndcapWallKaptonVol", Endcap_Kapton, description.material("Kapton"));
  dd4hep::Volume lvEndcapWallCarbon =
      dd4hep::Volume("lvEndcapWallCarbonVol", Endcap_Carbon, description.material("CarbonFiber"));

  dd4hep::Tube InnerWall_Carbon(inner_radius, inner_radius + CarbonInnerWallThick, halflength);
  dd4hep::Tube InnerWall_Copper(inner_radius + CarbonInnerWallThick,
                                inner_radius + CarbonInnerWallThick + CopperInnerWallThick, halflength);
  dd4hep::Tube InnerWall_Gas(inner_radius + CarbonInnerWallThick + CopperInnerWallThick,
                             inner_radius + envelop_Inner_thickness, halflength);

  dd4hep::Volume lvInnerWallCarbon =
      dd4hep::Volume("lvInnerWallCarbonVol", InnerWall_Carbon, description.material("CarbonFiber"));
  dd4hep::Volume lvInnerWallCopper =
      dd4hep::Volume("lvInnerWallCopperVol", InnerWall_Copper, description.material("G4_Cu"));
  dd4hep::Volume lvInnerWallGas =
      dd4hep::Volume("lvInnerWallGasVol", InnerWall_Gas, description.material("GasHe_90Isob_10"));

  dd4hep::Tube OuterWall_Copper(outer_radius - envelop_Outer_thickness,
                                outer_radius - Carbon1OuterWallThick - Carbon2OuterWallThick - FoamOuterWallThick,
                                halflength);
  dd4hep::Tube OuterWall_Carbon1(outer_radius - Carbon1OuterWallThick - Carbon2OuterWallThick - FoamOuterWallThick,
                                 outer_radius - Carbon2OuterWallThick - FoamOuterWallThick, halflength);
  dd4hep::Tube OuterWall_Foam(outer_radius - Carbon2OuterWallThick - FoamOuterWallThick,
                              outer_radius - Carbon2OuterWallThick, halflength);
  dd4hep::Tube OuterWall_Carbon2(outer_radius - Carbon2OuterWallThick, outer_radius, halflength);

  dd4hep::Volume lvOuterWallCarbon1 =
      dd4hep::Volume("lvOuterWallCarbon1Vol", OuterWall_Carbon1, description.material("CarbonFiber"));
  dd4hep::Volume lvOuterWallCarbon2 =
      dd4hep::Volume("lvOuterWallCarbon2Vol", OuterWall_Carbon2, description.material("CarbonFiber"));
  dd4hep::Volume lvOuterWallCopper =
      dd4hep::Volume("lvOuterWallCopperVol", OuterWall_Copper, description.material("G4_Cu"));
  dd4hep::Volume lvOuterWallFoam =
      dd4hep::Volume("lvOuterWallFoamVol", OuterWall_Foam, description.material("GasHe_90Isob_10"));


  //------------------------------------------------------------------------
  // Now we are ready to loop over the SuperLayers and fill the gas volume!
  //------------------------------------------------------------------------

  std::vector<dd4hep::Volume> lvLayerVol;
  std::vector<dd4hep::Volume> lvFwireVol, lvGwireVol;

  string wirecol, gascol, wholeHyperboloidVolumeName;
  string lvFwireName, lvSwireName;

  struct wire guard_wires{}, field_wires_bottom{}, field_wires_center{}, field_wires_top{}, sense_wires{};

  for (int SL = 0; SL < nSuperLayer; ++SL) {

    nSWire = nSWireFirstLayer + SL * nSDeltaWire;
    sense_wire_delta_phi = 2. * TMath::Pi() / nSWire;
    nFWireTopAndBottom = nHorizontalFtoSWireRatio * nSWire;
    phi_offset = 2. * TMath::Pi() / nFWireTopAndBottom;
    nFWireTopOrBottom = nFWireTopAndBottom / 2;
    if (ceilf(nFWireTopOrBottom) != nFWireTopOrBottom)
      throw std::runtime_error("Error: Failed to build CDCH. Please make sure that '(nFtoSWireRatio - nCenterFWirePerCell) * (nSWireFirstLayer + SuperLayerIndex * nSDeltaWire)' is always an even number");
    delta_phi_top_bottom = 2.0 * phi_offset;
    scaleFactor = (1.0 + TMath::Pi() / nSWire) / (1.0 - TMath::Pi() / nSWire); // used to scale the radial extent of each layer which grows towards outer radius
    dropFactor = (1.0 / cos(halfalpha) - 1.0); // used to determine the radius of the hyperboloid in z = +- halflength with r_out = r_min + r_min * dropFactor
    epsilonFactor = sin(halfalpha) / halflength;
    phi_offset_cell_start = -0.5 * sense_wire_delta_phi;

    gascol = "vCDCH:Gas1";
    if (SL % 3 == 0)
      gascol = "vCDCH:Gas1";
    else if ((SL + 1) % 3 == 0)
      gascol = "vCDCH:Gas2";
    else if ((SL + 2) % 3 == 0)
      gascol = "vCDCH:Gas3";

    if (SL % 3 == 0)
      wirecol = "vCDCH:Wire1";
    else if ((SL + 1) % 3 == 0)
      wirecol = "vCDCH:Wire2";
    else if ((SL + 2) % 3 == 0)
      wirecol = "vCDCH:Wire3";

    if (SL == 0) {// SL = 0 is special due to the guard wires and the first field wires that lie outside of the sensitive volume
      double stereoOut0 = atan(sqrt(diff_of_squares((inGuardRad - inGWradii) + (inGuardRad - inGWradii) * dropFactor, inGuardRad - inGWradii)) / halflength);

      dd4hep::Hyperboloid HypeLayer0(guard_layer_inner_radius_at_z_0, 0.0, layer_outer_radius_at_z_0 - secure, stereoOut0, halflength);
      lvLayerVol.push_back(dd4hep::Volume("hyperboloid_inner_guard_layer", HypeLayer0, description.material("GasHe_90Isob_10")));
      lvLayerVol.back().setVisAttributes(description, gascol);

      epsilonInGwRing = atan(inGuardRad * (1.0 + dropFactor) * epsilonFactor);
      zlength = halflength;
      zlength -= sin(epsilonInGwRing) * inGWradii;
      zlength /= cos(epsilonInGwRing);

      guard_wires.mother_volume = lvLayerVol.back();
      guard_wires.type = "G";
      guard_wires.num = nFWireTopOrBottom; //(#guard wires == # field wires)
      guard_wires.radius = inGuardRad - inGWradii;
      guard_wires.delta_phi = 2. * TMath::Pi() / guard_wires.num;
      guard_wires.phioffset = phi_offset_cell_start;
      guard_wires.stereo = epsilonInGwRing;
      guard_wires.thickness = 0.5 * InGWireShellThickIn * wireThicknessDilution;  // half the inner thickness as radius of tube
      guard_wires.halflength = zlength;
      guard_wires.name = string("Gwire_inner_stereoplus");

      dd4hep::Tube Gwire(0.0, guard_wires.thickness, halflength);
      lvGwireVol.push_back(dd4hep::Volume("Gwire_inner", Gwire, description.material("G4_Al")));
      lvGwireVol.back().setVisAttributes(description, wirecol);

      guard_wires.volume = lvGwireVol.back();
      CDCHBuild::PlaceGuardWires(guard_wires, FWireShellThickOut, halflength, SL, -1);

      guard_wires.volume = lvGwireVol.back(); // needed because applyWireCoating acts on it
      guard_wires.radius = inGuardRad + inGWradii + extShiftFW;
      guard_wires.phioffset = phi_offset_cell_start + phi_offset;
      guard_wires.stereo = -1.0 * epsilonInGwRing;
      guard_wires.name = string("Gwire_inner_stereominus");
      CDCHBuild::PlaceGuardWires(guard_wires, FWireShellThickOut, halflength, SL, -1);

      drop = radius_layer_0 * dropFactor;
      radius_layer = radius_layer_0 + drop;
      epsilon = atan(radius_layer * epsilonFactor);
      layer_inner_radius_at_z_0 = radius_layer_0 - FWradii - 2.0 * secure;
      layer_inner_radius_at_z_end = layer_inner_radius_at_z_0 + drop;
      layer_outer_radius_at_z_0 = radius_layer_0 + FWradii;
      layer_outer_radius_at_z_end = layer_outer_radius_at_z_0 + drop;
      epsilonIn = atan(sqrt(pow(layer_inner_radius_at_z_end, 2) - pow(layer_inner_radius_at_z_0, 2)) / halflength);
      epsilonOut = atan(sqrt(pow(layer_outer_radius_at_z_end, 2) - pow(layer_outer_radius_at_z_0, 2)) / halflength);


      zlength = halflength;
      zlength -= sin(epsilon) * FWradii;
      zlength /= cos(epsilon);

      field_wires_top.type = "F";
      field_wires_top.num = nFWireTopOrBottom;
      field_wires_top.radius = layer_inner_radius_at_z_0 - FWradii - extShiftFW;
      field_wires_top.delta_phi = 2. * TMath::Pi() /field_wires_top.num;;
      field_wires_top.phioffset = phi_offset_cell_start + cellStaggering - phi_offset;
      field_wires_top.stereo = sign_epsilon * epsilon;
      field_wires_top.thickness = 0.5 * FWireShellThickIn * wireThicknessDilution;
      field_wires_top.halflength = zlength;

      dd4hep::Hyperboloid HypeLayer1(layer_inner_radius_at_z_0, epsilonIn, layer_outer_radius_at_z_0 + field_wires_top.thickness + 0.5 * FWireShellThickOut, epsilonOut, halflength);
      lvLayerVol.push_back(dd4hep::Volume("hyperboloid_first_field_wire_ring", HypeLayer1, description.material("GasHe_90Isob_10")));
      lvLayerVol.back().setVisAttributes(description, gascol);
      field_wires_top.mother_volume = lvLayerVol.back();

      lvFwireName = dd4hep::_toString(SL, "lvFwire_%d_init");
      field_wires_top.name = string(lvFwireName);

      dd4hep::Tube Fwire(0.0, field_wires_top.thickness, halflength);
      lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material("G4_Al")));
      lvFwireVol.back().setVisAttributes(description, wirecol);

      field_wires_top.volume = lvFwireVol.back();
      CDCHBuild::PlaceGuardWires(field_wires_top, FWireShellThickOut, halflength, SL, -1);

      radius_layer_0 += FWradii;

    } else {
      delta_radius_layer = 2. * TMath::Pi() * layer_outer_radius_at_z_0 / (nSWire - TMath::Pi());
    }

    //------------------------------------------------------------------------
    // Starting the layer loop. nLayer=8
    //------------------------------------------------------------------------

    for (int ilayer = 0; ilayer < nLayer; ilayer++) {
      //------------------------------------------------------------------------
      // Fill the geometry parameters of the layer. Each layer lies
      // on top of the following one, so new layerIn = old layerOut
      //------------------------------------------------------------------------

      sense_wire_radius_at_z_0 = radius_layer_0 + 0.5 * delta_radius_layer;
      sign_epsilon *= -1;

      layer_inner_radius_at_z_0 = layer_outer_radius_at_z_0;
      layer_inner_radius_at_z_end = layer_outer_radius_at_z_end;
      epsilonIn = epsilonOut;

      layer_outer_radius_at_z_0 = layer_inner_radius_at_z_0 + FWireDiameter + 2.0 * secure;
      layer_outer_radius_at_z_end = layer_outer_radius_at_z_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(layer_outer_radius_at_z_end, layer_outer_radius_at_z_0)) / halflength);

      zlength = halflength;

      // save bottom layer inner radius and epsilon before they are modified to build the whole layer hyperboloid volume 
      layer_inner_radius_at_z_end_whole_cell = layer_inner_radius_at_z_0;
      epsilonIn_whole_cell = epsilonIn;

      //------------------------------------------------------------------------
      // Reduce zlength to avoid volume extrusions and check the staggering
      //------------------------------------------------------------------------

      zlength -= sin(epsilon) * FWradii;
      zlength /= cos(epsilon);

      if (ilayer % 2 == 1)
        cellStaggering = phi_offset;
      else
        cellStaggering = 0.0;

      //------------------------------------------------------------------------
      // Fill the field wire struct with all the relevant information
      // This is the bottom of the cell.
      //------------------------------------------------------------------------

      field_wires_bottom.type = "F";
      field_wires_bottom.num = nFWireTopOrBottom;
      field_wires_bottom.radius = layer_inner_radius_at_z_0 + FWradii + extShiftFW;
      field_wires_bottom.delta_phi = delta_phi_top_bottom;
      field_wires_bottom.phioffset = phi_offset_cell_start + cellStaggering;
      field_wires_bottom.stereo = sign_epsilon * epsilon;
      field_wires_bottom.thickness = 0.5 * FWireShellThickIn * wireThicknessDilution;
      field_wires_bottom.halflength = zlength;
      field_wires_bottom.name = dd4hep::_toString(SL, "Wire_SL%d") + dd4hep::_toString(ilayer, "_layer%d") + string("_type") + field_wires_bottom.type + dd4hep::_toString(field_wires_bottom.stereo, "_stereo%f_bottom");

      //------------------------------------------------------------------------
      // Define the field wire name and build the field wire volume
      //------------------------------------------------------------------------

      lvFwireName = dd4hep::_toString(SL, "lvFwire_%d") + dd4hep::_toString(ilayer, "_%d");

      dd4hep::Tube Fwire(0.0, field_wires_bottom.thickness, halflength);
      dd4hep::Volume FwireVol(lvFwireName, Fwire, description.material("G4_Al"));
      FwireVol.setVisAttributes(description, wirecol);

      //------------------------------------------------------------------------
      // Add the field wire volume to the struct 
      //------------------------------------------------------------------------

      field_wires_bottom.volume = FwireVol;
      apply_wire_coating(field_wires_bottom, FWireShellThickOut, halflength);

      //------------------------------------------------------------------------
      // Next, fill the geometry parameters of the central layer.
      //------------------------------------------------------------------------

      iradius = radius_layer_0;
      radius_layer_0 += delta_radius_layer;
      drop = radius_layer_0 * dropFactor;

      layer_inner_radius_at_z_0 = layer_outer_radius_at_z_0;
      layer_inner_radius_at_z_end = layer_outer_radius_at_z_end;
      epsilonIn = epsilonOut;
      layer_outer_radius_at_z_0 = radius_layer_0 - FWireDiameter - 2.0 * secure;
      layer_outer_radius_at_z_end = layer_outer_radius_at_z_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(layer_outer_radius_at_z_end, layer_outer_radius_at_z_0)) / halflength);
      zlength = halflength;

      //------------------------------------------------------------------------
      // Reduce zlength to avoid volume extrusions
      //------------------------------------------------------------------------

      zlength -= sin(epsilon) * 0.5 * FWireShellThickIn;
      zlength /= cos(epsilon);


      ////------------------------------------------------------------------------
      //// Tune the radius and epsilon of the central field wires
      ////------------------------------------------------------------------------

      idelta_radius = delta_radius_layer * 0.5;
      iradius += idelta_radius;
      //epsilon = atan(iradius * (1.0 + dropFactor) * epsilonFactor);

      //------------------------------------------------------------------------
      // Fill the sense wire struct with all the relevant information
      //------------------------------------------------------------------------

      sense_wires.type = "S";
      sense_wires.num = nSWire;
      sense_wires.radius = sense_wire_radius_at_z_0;
      sense_wires.delta_phi = sense_wire_delta_phi;
      sense_wires.phioffset = cellStaggering;
      sense_wires.stereo = sign_epsilon * epsilon;
      sense_wires.thickness = 0.5 * SWireShellThickIn * wireThicknessDilution;
      sense_wires.halflength = zlength;
      sense_wires.name = dd4hep::_toString(SL, "Wire_SL%d") + dd4hep::_toString(ilayer, "_layer%d") + string("_type") + sense_wires.type + dd4hep::_toString(sense_wires.stereo, "_stereo%f");

      //------------------------------------------------------------------------
      // Define the sense wire name and build the sense wire volume
      //------------------------------------------------------------------------

      lvSwireName = dd4hep::_toString(SL, "lvSwire_%d") + dd4hep::_toString(ilayer, "_%d");

      dd4hep::Tube Swire(0.0, sense_wires.thickness, halflength);
      dd4hep::Volume lvSwireVol(lvSwireName, Swire, description.material("G4_W"));
      lvSwireVol.setVisAttributes(description, wirecol);
      sense_wires.volume = lvSwireVol;
      apply_wire_coating(sense_wires, SWireShellThickOut, halflength, "G4_Au");

      //------------------------------------------------------------------------
      // Fill the central field wire struct with all the relevant information
      //------------------------------------------------------------------------

      field_wires_center.type = "F";
      field_wires_center.num = nSWire * nCenterFWirePerCell;
      field_wires_center.radius = sense_wire_radius_at_z_0; //iradius;
      field_wires_center.delta_phi = 2. * TMath::Pi() / field_wires_center.num;
      field_wires_center.phioffset = phi_offset_cell_start + cellStaggering;
      field_wires_center.stereo = sign_epsilon * epsilon;
      field_wires_center.thickness = 0.5 * centerFWireShellThickIn * wireThicknessDilution;
      field_wires_center.halflength = zlength;
      field_wires_center.volume = FwireVol;
      field_wires_center.name = dd4hep::_toString(SL, "Wire_SL%d") + dd4hep::_toString(ilayer, "_layer%d") + string("_type") + field_wires_center.type + dd4hep::_toString(field_wires_center.stereo, "_stereo%f_center");
      apply_wire_coating(field_wires_center, centerFWireShellThickOut, halflength);

      // derive a phi offset to englobe the 'left' (clockwise) wires inside the sensitive volume (center field wires are the thickest)
      double phi_offset_to_englobe_wires = atan((field_wires_center.thickness + 0.5 * centerFWireShellThickOut) / field_wires_center.radius);

      //------------------------------------------------------------------------
      // Next, fill the geometry parameters of the upper layer.
      //------------------------------------------------------------------------

      layer_inner_radius_at_z_0 = layer_outer_radius_at_z_0;
      layer_inner_radius_at_z_end = layer_outer_radius_at_z_end;
      epsilonIn = epsilonOut;
      layer_outer_radius_at_z_0 = layer_inner_radius_at_z_0 + FWireDiameter + 2.0 * secure;
      layer_outer_radius_at_z_end = layer_outer_radius_at_z_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(layer_outer_radius_at_z_end, layer_outer_radius_at_z_0)) / halflength);
      zlength = halflength;

      // Create hyperboloid volume of the whole layer for the cell sensitive volume definition, not needed per se but helps having a well balanced volume tree (having too many volumes inside one mother volume harms performance)
      wholeHyperboloidVolumeName = dd4hep::_toString(SL, "hyperboloid_SL_%d") + dd4hep::_toString(ilayer, "_layer_%d");
      double radial_inner_offset_to_englobe_wires = field_wires_bottom.thickness + 0.5 * FWireShellThickOut;
      double radial_outer_offset_to_englobe_wires = field_wires_top.thickness + 0.5 * FWireShellThickOut;
      dd4hep::Hyperboloid whole_layer_hyperboloid = dd4hep::Hyperboloid(layer_inner_radius_at_z_end_whole_cell - radial_inner_offset_to_englobe_wires, epsilonIn_whole_cell, layer_outer_radius_at_z_0 + radial_outer_offset_to_englobe_wires, epsilonOut, zlength);
      dd4hep::Volume whole_layer_hyperboloid_volume = dd4hep::Volume(wholeHyperboloidVolumeName, whole_layer_hyperboloid, description.material("GasHe_90Isob_10"));
      whole_layer_hyperboloid_volume.setVisAttributes(description, gascol);
      //registerVolume(wholeHyperboloidVolumeName, whole_layer_hyperboloid_volume);
      dd4hep::PlacedVolume whole_layer_hyperboloid_placedVolume;
      whole_layer_hyperboloid_placedVolume = parentVol.placeVolume(whole_layer_hyperboloid_volume);
      whole_layer_hyperboloid_placedVolume.addPhysVolID("superLayer", SL).addPhysVolID("layer", ilayer);
      dd4hep::DetElement whole_layer_hyperboloid_detElement(parent, "superLayer_" + dd4hep::_toString(SL) + "_layer_" + dd4hep::_toString(ilayer) + "_hyperboloid", SL * nLayer + ilayer);
      whole_layer_hyperboloid_detElement.setPlacement(whole_layer_hyperboloid_placedVolume);

      //------------------------------------------------------------------------
      // Reduce zlength to avoid volume extrusions
      //------------------------------------------------------------------------

      zlength -= sin(epsilon) * FWradii;
      zlength /= cos(epsilon);

      //------------------------------------------------------------------------
      // Fill the field wire struct with all the relevant information
      // This is the top of the cell.
      //------------------------------------------------------------------------

      field_wires_top.type = "F";
      field_wires_top.num = nFWireTopOrBottom;
      field_wires_top.radius = layer_inner_radius_at_z_0 - FWradii - extShiftFW;
      field_wires_top.delta_phi = delta_phi_top_bottom;
      field_wires_top.phioffset = phi_offset_cell_start + cellStaggering;
      field_wires_top.stereo = sign_epsilon * epsilon;
      field_wires_top.thickness = 0.5 * FWireShellThickIn * wireThicknessDilution;
      field_wires_top.halflength = zlength;
      field_wires_top.volume = FwireVol;
      field_wires_top.name = dd4hep::_toString(SL, "Wire_SL%d") + dd4hep::_toString(ilayer, "_layer%d") + string("_type") + field_wires_top.type + dd4hep::_toString(field_wires_top.stereo, "_stereo%f_top");
      apply_wire_coating(field_wires_top, FWireShellThickOut, halflength);

      //if(setWireSensitive){
      //  field_wires_bottom.volume.setSensitiveDetector(sens_det);
      //  field_wires_center.volume.setSensitiveDetector(sens_det);
      //  sense_wires.volume.setSensitiveDetector(sens_det);
      //  field_wires_top.volume.setSensitiveDetector(sens_det);
      //}

      // Create the tube segment volume to identify the cell sensitive regions
      // FIXME: this leads to some volume extrusion (corner of the tube going outside of the hyperboloid mother volume)
      // Using the intersection with the hyperboloid mother volume solves it but severely impacts the perfomance (memory and CPU) of the geometry building
      // Other paths to investigate are twisted tubes (caveat: no TGeo shape equivalent), extruded volumes or tessalated solids
      dd4hep::Tube cellID_tube_segment(layer_inner_radius_at_z_end_whole_cell - radial_inner_offset_to_englobe_wires, layer_outer_radius_at_z_0 + radial_outer_offset_to_englobe_wires, halflength, (- sense_wires.delta_phi / 2.0) - phi_offset_to_englobe_wires, (sense_wires.delta_phi / 2.0) - phi_offset_to_englobe_wires);
      
      // Radial translation 
      dd4hep::Translation3D radial_translation_sense_wire(sense_wires.radius, 0., 0.);
      // stereo rotation
      dd4hep::RotationX rot_stereo_sense_wire(sense_wires.stereo);
      // extract the number of wire ratio to place field wires in the loop for sense wires
      // it is not very elegant but the sense wire define the sensitive volume in which wires are placed 
      float middle_to_middle_num_wire_ratio = field_wires_center.num/float(sense_wires.num);
      float middle_to_bottom_num_wire_ratio = field_wires_bottom.num/float(sense_wires.num);
      float middle_to_top_num_wire_ratio = field_wires_top.num/float(sense_wires.num);
      if(ceilf(middle_to_middle_num_wire_ratio) != middle_to_middle_num_wire_ratio || ceilf(middle_to_bottom_num_wire_ratio) != middle_to_bottom_num_wire_ratio || ceilf(middle_to_top_num_wire_ratio) != middle_to_top_num_wire_ratio)
          throw std::runtime_error("Error: Failed to build CDCH. Please make sure that the number of wires in top/center cell rings is always a multiple of the number of wires in the middle of the cell");
      // loop to arrange the wires in phi, starting with the sense wires
      for (int phi_index = 0; phi_index < sense_wires.num; phi_index++) {
        // Place the sensitive volume inside the hyperbioloid with phi and stereo angle rotation
        string cellID_volume_name = dd4hep::_toString(SL, "cellIDvolume_SL_%d") + dd4hep::_toString(ilayer, "_layer_%d") + dd4hep::_toString(phi_index, "_phi_%d");
        dd4hep::Volume cellID_volume = dd4hep::Volume(cellID_volume_name, cellID_tube_segment, description.material("GasHe_90Isob_10"));
        cellID_volume.setSensitiveDetector(sens_det);
        // phi rotation
        double phi_angle_sense_wire_rotation = sense_wires.phioffset + sense_wires.delta_phi * phi_index;
        dd4hep::RotationZ rot_phi_sense_wire(phi_angle_sense_wire_rotation);
        dd4hep::PlacedVolume cellID_placedvolume = whole_layer_hyperboloid_volume.placeVolume(cellID_volume, dd4hep::Transform3D(rot_phi_sense_wire * rot_stereo_sense_wire));
        cellID_placedvolume.addPhysVolID("phi", phi_index).addPhysVolID("hitorigin", 0).addPhysVolID("stereo", sense_wires.stereo > 0 ? 0 : 1).addPhysVolID("layerInCell", 0);
        dd4hep::DetElement cellID_detElement(whole_layer_hyperboloid_detElement, "superLayer_" + dd4hep::_toString(SL) + "_layer_" + dd4hep::_toString(ilayer) + "_phi_" + dd4hep::_toString(phi_index) + "_cellID", phi_index);
        cellID_detElement.setPlacement(cellID_placedvolume);

        // place the wires
        // sense wires in the radial middle of the cell
        dd4hep::PlacedVolume sense_wire_placedvolume = cellID_volume.placeVolume(sense_wires.volume, dd4hep::Transform3D(radial_translation_sense_wire)); // only the radial translation is needed as the cellID volume already had the stereo and phi angle rotation
        // add the sense wire as detElement to be able to retrive its matrix with DD4hep 1.23 (with later verion we can use Volumes daugthers)
        dd4hep::DetElement senseWire_detElement(cellID_detElement, "superLayer_" + dd4hep::_toString(SL) + "_layer_" + dd4hep::_toString(ilayer) + "_phi_" + dd4hep::_toString(phi_index) + "_wire", phi_index);
        senseWire_detElement.setPlacement(sense_wire_placedvolume);

        // field wires: the cellID volume has the stereo angle rotation of the sense wires which is different than the top/bototm wires --> Rotations are defined on top of the already applied cellID volume rotation
        // bottom field wires
        for(int sub_phi_index = phi_index * middle_to_bottom_num_wire_ratio; sub_phi_index < (phi_index * middle_to_bottom_num_wire_ratio) + middle_to_bottom_num_wire_ratio; sub_phi_index++){
          // phi rotation
          dd4hep::RotationZ rot_phi_bottom_wire((field_wires_bottom.phioffset + field_wires_bottom.delta_phi * sub_phi_index) - phi_angle_sense_wire_rotation);
          cellID_volume.placeVolume(field_wires_bottom.volume, dd4hep::Transform3D(rot_phi_bottom_wire * dd4hep::Translation3D(field_wires_bottom.radius, 0., 0.) * dd4hep::RotationX(field_wires_bottom.stereo - sense_wires.stereo)));
          //dd4hep::PlacedVolume field_wire_bottom_placedvolume = cellID_volume.placeVolume(field_wires_bottom.volume, dd4hep::Transform3D(rot_phi_bottom_wire * dd4hep::Translation3D(field_wires_bottom.radius, 0., 0.) * dd4hep::RotationX(field_wires_bottom.stereo - sense_wires.stereo)));
          //if(setWireSensitive)
          //  field_wire_bottom_placedvolume.addPhysVolID("phi", sub_phi_index).addPhysVolID("hitorigin", 2).addPhysVolID("stereo", field_wires_center.stereo > 0 ? 0 : 1).addPhysVolID("layerInCell", 1);
        }

        // central field wires
        for(int sub_phi_index = phi_index * middle_to_middle_num_wire_ratio; sub_phi_index < (phi_index * middle_to_middle_num_wire_ratio) + middle_to_middle_num_wire_ratio; sub_phi_index++){
          // phi rotation
          dd4hep::RotationZ rot_phi_center_wire((field_wires_center.phioffset + field_wires_center.delta_phi * sub_phi_index) - phi_angle_sense_wire_rotation);
          cellID_volume.placeVolume(field_wires_center.volume, dd4hep::Transform3D(rot_phi_center_wire * dd4hep::Translation3D(field_wires_center.radius, 0., 0.) * dd4hep::RotationX(field_wires_center.stereo - sense_wires.stereo)));
          //dd4hep::PlacedVolume field_wire_center_placedvolume = cellID_volume.placeVolume(field_wires_center.volume, dd4hep::Transform3D(rot_phi_center_wire * dd4hep::Translation3D(field_wires_center.radius, 0., 0.) * dd4hep::RotationX(field_wires_center.stereo - sense_wires.stereo)));
          //if(setWireSensitive)
          //  field_wire_center_placedvolume.addPhysVolID("phi", sub_phi_index).addPhysVolID("hitorigin", 2).addPhysVolID("stereo", field_wires_center.stereo > 0 ? 0 : 1).addPhysVolID("layerInCell", 2);
        }

        // top field wires
        for(int sub_phi_index = phi_index * middle_to_top_num_wire_ratio; sub_phi_index < (phi_index * middle_to_top_num_wire_ratio) + middle_to_top_num_wire_ratio; sub_phi_index++){
          // phi rotation
          dd4hep::RotationZ rot_phi_top_wire((field_wires_top.phioffset + field_wires_top.delta_phi * sub_phi_index) - phi_angle_sense_wire_rotation);
          cellID_volume.placeVolume(field_wires_top.volume, dd4hep::Transform3D(rot_phi_top_wire * dd4hep::Translation3D(field_wires_top.radius, 0., 0.) * dd4hep::RotationX(field_wires_top.stereo - sense_wires.stereo)));
          //dd4hep::PlacedVolume field_wire_top_placedvolume = cellID_volume.placeVolume(field_wires_top.volume, dd4hep::Transform3D(rot_phi_top_wire * dd4hep::Translation3D(field_wires_top.radius, 0., 0.) * dd4hep::RotationX(field_wires_top.stereo - sense_wires.stereo)));
          //if(setWireSensitive)
          //  field_wire_top_placedvolume.addPhysVolID("phi", sub_phi_index).addPhysVolID("hitorigin", 2).addPhysVolID("stereo", field_wires_center.stereo > 0 ? 0 : 1).addPhysVolID("layerInCell", 3);
        }
      }

      //------------------------------------------------------------------------
      // Scale the delta radius of the layer for next iteration
      //------------------------------------------------------------------------

      delta_radius_layer *= scaleFactor;
      epsilon = atan(iradius * (1.0 + dropFactor) * epsilonFactor); // FIXME the stereo angle is here assumed to be constant for a given cell while it should be different for the bottom, middle and top rings of the cell
    }

    if (SL == (nSuperLayer - 1)) {// the last super layer is special since we need to add the field wires outside of the sensitive volume and the guard wires

      // Take care of the field wires outside of the sensitive volume
      layer_inner_radius_at_z_0 = layer_outer_radius_at_z_0;
      layer_inner_radius_at_z_end = layer_outer_radius_at_z_end;
      epsilonIn = epsilonOut;
      layer_outer_radius_at_z_0 = radius_layer_0 + FWireDiameter + 2.0 * secure;
      layer_outer_radius_at_z_end = layer_outer_radius_at_z_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(layer_outer_radius_at_z_end, layer_outer_radius_at_z_0)) / halflength);

      dd4hep::Hyperboloid HypeLayerOut(layer_inner_radius_at_z_0, epsilonIn, layer_outer_radius_at_z_0, epsilonOut, halflength);
      lvLayerVol.push_back(dd4hep::Volume("hyperboloid_last_field_wire_ring", HypeLayerOut, description.material("GasHe_90Isob_10")));
      lvLayerVol.back().setVisAttributes(description, gascol);

      zlength = halflength;
      zlength -= sin(epsilon) * FWradii;
      zlength /= cos(epsilon);

      field_wires_bottom.mother_volume = lvLayerVol.back();
      field_wires_bottom.type = "F";
      field_wires_bottom.num = nFWireTopOrBottom;
      field_wires_bottom.radius = layer_inner_radius_at_z_0 + FWradii + extShiftFW;
      field_wires_bottom.delta_phi = 2. * TMath::Pi() / field_wires_bottom.num;
      field_wires_bottom.phioffset = phi_offset_cell_start + cellStaggering + phi_offset;
      field_wires_bottom.stereo = -1. * sign_epsilon * epsilon;
      field_wires_bottom.thickness = 0.5 * FWireShellThickIn * wireThicknessDilution;
      field_wires_bottom.halflength = zlength;

      lvFwireName = dd4hep::_toString(SL, "lvFwire_%d_out");
      field_wires_bottom.name = lvFwireName;

      dd4hep::Tube Fwire(0.0, field_wires_bottom.thickness, halflength);
      lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material("G4_Al")));
      lvFwireVol.back().setVisAttributes(description, wirecol);

      field_wires_bottom.volume = lvFwireVol.back();
      CDCHBuild::PlaceGuardWires(field_wires_bottom, FWireShellThickOut, halflength, SL, -1);

      //------------------------------------------------------------------------
      // Start placing the outer layer of guard wires (#guard wires == # field wires)
      //------------------------------------------------------------------------

      layer_inner_radius_at_z_0 = layer_outer_radius_at_z_0;
      layer_inner_radius_at_z_end = layer_outer_radius_at_z_end;
      epsilonIn = epsilonOut;
      layer_outer_radius_at_z_0 = radius_layer_0 + FWireDiameter + 2.0 * secure;
      layer_outer_radius_at_z_end = layer_outer_radius_at_z_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(layer_outer_radius_at_z_end, layer_outer_radius_at_z_0)) / halflength);

      dd4hep::Hyperboloid HypeLayerOutG(layer_inner_radius_at_z_0, epsilonOut, outer_radius - envelop_Outer_thickness - 0.0001,
                                        0.0, halflength);
      lvLayerVol.push_back(dd4hep::Volume("hyperboloid_outer_guard_layer", HypeLayerOutG, description.material("GasHe_90Isob_10")));
      lvLayerVol.back().setVisAttributes(description, gascol);

      epsilonOutGwRing = atan(outGuardRad * (1.0 + dropFactor) * epsilonFactor);
      zlength = halflength;
      zlength -= sin(epsilonOutGwRing) * inGWradii;
      zlength /= cos(epsilonOutGwRing);

      guard_wires.mother_volume = lvLayerVol.back();
      guard_wires.type = "G";
      guard_wires.num = nFWireTopOrBottom;
      guard_wires.radius = outGuardRad - inGWradii;
      guard_wires.delta_phi = 2. * TMath::Pi() / guard_wires.num;
      guard_wires.phioffset = phi_offset_cell_start;
      guard_wires.stereo = epsilonOutGwRing;
      guard_wires.thickness = 0.5 * OutGWireShellThickIn * wireThicknessDilution;
      guard_wires.halflength = zlength;
      guard_wires.name = string("Gwire_outer_stereominus");

      dd4hep::Tube Gwire(0.0, guard_wires.thickness, halflength);
      lvGwireVol.push_back(dd4hep::Volume("Gwire_outer", Gwire, description.material("G4_Al")));
      lvGwireVol.back().setVisAttributes(description, wirecol);

      guard_wires.volume = lvGwireVol.back();
      CDCHBuild::PlaceGuardWires(guard_wires, OutGWireShellThickOut, halflength, SL, -1);

      guard_wires.volume = lvGwireVol.back(); // needed because applyWireCoating acts on it
      guard_wires.radius = outGuardRad + inGWradii + extShiftFW;
      guard_wires.phioffset = phi_offset_cell_start + phi_offset;
      guard_wires.stereo = -1.0 * epsilonOutGwRing;
      guard_wires.name = string("Gwire_outer_stereoplus");
      CDCHBuild::PlaceGuardWires(guard_wires, OutGWireShellThickOut, halflength, SL, -1);
    }
  }


  Int_t sizeLayer = lvLayerVol.size();


  for (Int_t i = 0; i < sizeLayer; i++) {
    registerVolume(lvLayerVol.at(i).name(), lvLayerVol.at(i));
    //cout << "Placing Volume: " << lvLayerVol.at(i).name() << endl;
    //pv = parentVol.placeVolume(volume(lvLayerVol.at(i).name()));
    //CDCHDetector.setPlacement(pv);
    parentVol.placeVolume(volume(lvLayerVol.at(i).name()));
  }

  double PosEndcapGas = halflength + 0.5 * GasEndcapWallThick;
  double PosEndcapCopper = halflength + GasEndcapWallThick + 0.5 * CopperEndcapWallThick;
  double PosEndcapKapton = halflength + GasEndcapWallThick + CopperEndcapWallThick + 0.5 * KaptonEndcapWallThick;
  double PosEndcapCarbon =
      halflength + GasEndcapWallThick + CopperEndcapWallThick + KaptonEndcapWallThick + 0.5 * CarbonEndcapWallThick;

  parentVol.placeVolume(lvInnerWallCarbon);
  parentVol.placeVolume(lvInnerWallCopper);
  parentVol.placeVolume(lvInnerWallGas);
  parentVol.placeVolume(lvOuterWallCarbon1);
  parentVol.placeVolume(lvOuterWallCarbon2);
  parentVol.placeVolume(lvOuterWallCopper);
  parentVol.placeVolume(lvOuterWallFoam);
  parentVol.placeVolume(lvEndcapWallGas, dd4hep::Position(0., 0., PosEndcapGas));
  parentVol.placeVolume(lvEndcapWallCopper, dd4hep::Position(0., 0., PosEndcapCopper));
  parentVol.placeVolume(lvEndcapWallKapton, dd4hep::Position(0., 0., PosEndcapKapton));
  parentVol.placeVolume(lvEndcapWallCarbon, dd4hep::Position(0., 0., PosEndcapCarbon));
  parentVol.placeVolume(lvEndcapWallGas, dd4hep::Position(0., 0., -PosEndcapGas));
  parentVol.placeVolume(lvEndcapWallCopper, dd4hep::Position(0., 0., -PosEndcapCopper));
  parentVol.placeVolume(lvEndcapWallKapton, dd4hep::Position(0., 0., -PosEndcapKapton));
  parentVol.placeVolume(lvEndcapWallCarbon, dd4hep::Position(0., 0., -PosEndcapCarbon));
}
} //namespace

static dd4hep::Ref_t create_element(dd4hep::Detector& description, xml_h e, dd4hep::SensitiveDetector sens_det) {

  xml_det_t x_det = e;
  CDCHBuild builder(description, x_det, sens_det);
  string det_name = x_det.nameStr();

  dd4hep::printout(dd4hep::DEBUG, "CreateCDCH", "Detector name: %s with ID: %s", det_name.c_str(), x_det.id());

  DetElement CDCH_det = builder.detector;  // ( det_name, x_det.id() );
  dd4hep::Tube CDCH_envelope(dd4hep::_toDouble("CDCH:inner_radius"), dd4hep::_toDouble("CDCH:outer_radius"), dd4hep::_toDouble("CDCH:zHalfExtentWithServices"));

  dd4hep::Volume envelope("lvCDCH", CDCH_envelope, description.air());
  envelope.setVisAttributes(description, "vCDCH:Air");

  // ******************************************************
  // Build CDCH cable
  // ******************************************************

  builder.build_layer(CDCH_det, envelope, sens_det);

  dd4hep::printout(dd4hep::DEBUG, "CreateCDCH", "MotherVolume is: %s", envelope.name());
  dd4hep::xml::Dimension sdType = x_det.child(_U(sensitive));
  sens_det.setType(sdType.typeStr());

  builder.buildVolumes(e);
  builder.placeDaughters(CDCH_det, envelope, e);

  // ******************************************************
  // Build CDCH cell and beam plug
  // ******************************************************

  //  builder.build_cell();
  //  builder.build_beamplug();

  // ******************************************************
  // Assemble CDCH
  // ******************************************************

  //  builder.build_CDCH( Ecal_det, envelope );

  // ******************************************************
  // Place the CDCH in the world
  // ******************************************************

  PlacedVolume pv;
  pv = builder.placeDetector(envelope);
  pv.addPhysVolID("system", x_det.id());

  return CDCH_det;
}

DECLARE_DETELEMENT(DriftChamber_o1_v01, create_element)

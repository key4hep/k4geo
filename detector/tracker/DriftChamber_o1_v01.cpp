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
  string type;
  int num;
  double radius;
  double phi;
  double phioffset;
  double stereo;
  double thickness;
  double halflength;
  dd4hep::Volume volume;
  string name;
};

namespace {

struct CDCHBuild : public dd4hep::xml::tools::VolumeBuilder {
  std::vector<dd4hep::DetElement> deSuperLayer, deLayer, deSWire;

  CDCHBuild(dd4hep::Detector& description, xml_elt_t e, dd4hep::SensitiveDetector sens);

  double diff_of_squares(double a, double b);
  void apply_wire_coating(struct wire& w, double outwrap, double halflength, string material);
  void PlaceWires(struct wire& w, double outwrap, double halflength, int SL, int ilayer);
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

// deprecated function, only use it for the wires before and after the sensitive zone (mainly guard wires)
void CDCHBuild::PlaceWires(struct wire& w, double outwrap, double halflength, int SL = 999,
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

  //cout << "wirewrapname: " << wirewrapname << endl;

  string wirename = "Wire_SL";
  wirename += std::to_string(SL);
  wirename += "_layer";
  wirename += std::to_string(ilayer);
  wirename += "_type";
  wirename += w.type;
  wirename += "_stereo";
  wirename += std::to_string(w.stereo);

  apply_wire_coating(w, outwrap, halflength);

  //    registerVolume(WireWrapVol.name(), WireWrapVol);
  //    registerVolume(WireVol.name(), WireVol);

  // repeat the placement of wires over phi
  for (int n = 0; n < w.num; n++) {
    dd4hep::RotationZ iRot(w.phioffset + w.phi * n);
    if (n % 1 == 0) w.mother_volume.placeVolume(w.volume, dd4hep::Transform3D(iRot * T));
  }
}

void CDCHBuild::build_layer(DetElement parent, Volume parentVol, dd4hep::SensitiveDetector sens_det) {

  // ******************************************************
  // Loading parameters
  // ******************************************************

  double halfalpha = 0.5 * dd4hep::_toDouble("CDCH:alpha");
  double inner_radius = dd4hep::_toDouble("CDCH:r0");
  double outer_radius = dd4hep::_toDouble("CDCH:rOut");
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
  double OutGWireShellThickIn = dd4hep::_toDouble("CDCH:InGWireShellThickIn");
  double OutGWireShellThickOut = dd4hep::_toDouble("CDCH:InGWireShellThickOut");
  double secure = dd4hep::_toDouble("CDCH:secure");
  double capGasLayer = dd4hep::_toDouble("CDCH:capGasLayer");
  double extShiftFW = dd4hep::_toDouble("CDCH:extShiftFW");
  double cellDimension = dd4hep::_toDouble("CDCH:cellDimension");
  double inGuardRad = dd4hep::_toDouble("CDCH:inGuardRad");
  double outGuardRad = dd4hep::_toDouble("CDCH:outGuardRad");
  int nSDeltaWire = dd4hep::_toInt("CDCH:nSDeltaWire");
  int nSWire = dd4hep::_toInt("CDCH:nSWire");
  int nStoFWireRatio = dd4hep::_toInt("CDCH:nStoFWireRatio");
  int nVerticalFWire = dd4hep::_toInt("CDCH:nVerticalFWire");
  int nSuperLayer = dd4hep::_toInt("CDCH:nSuperLayer");
  int nLayer = dd4hep::_toInt("CDCH:nLayer");
  int nFieldWireShells = dd4hep::_toInt("CDCH:nFieldWireShells");
  //bool setWireSensitive = true; // FIXME: add the possibility to have wires sensitive (parameter in the xml) which could be useful for detailed chamber behavior studies, current attempt never lead to a hit in the wire, even with enlarged wires...
  double halflength = zHalfExtentWithServices - (GasEndcapWallThick + CopperEndcapWallThick + KaptonEndcapWallThick + CarbonEndcapWallThick); // this will be the sensitive volume z extent

  double epsilon = 0.0;
  double phi_layer = 0.0;
  double phi_layer1 = 0.0;
  int nFWire = 0;
  int nFWire1 = 0;
  int num_wire = 0;
  int nHorizontalFWire = nStoFWireRatio - nVerticalFWire;
  int sign_epsilon = -1;
  double phi = 0.0;
  double scaleFactor = 0.0;
  double dropFactor = 0.0;
  double epsilonFactor = 0.0;
  double delta_radius_layer = cellDimension;
  double senseWireRing_radius_0 = 0.0;
  double iradius = 0.0;
  double idelta_radius = 0.0;

  double envelop_Inner_thickness = CarbonInnerWallThick + CopperInnerWallThick + GasInnerWallThick;
  double envelop_Outer_thickness =
      Carbon1OuterWallThick + Carbon2OuterWallThick + CopperOuterWallThick + FoamOuterWallThick;
  double FWireDiameter = FWireShellThickIn + FWireShellThickOut;
  double FWradii = 0.5 * FWireDiameter;
  double SWireDiameter = SWireShellThickIn + SWireShellThickOut;
  double SWradii = 0.5 * SWireDiameter;
  double inGWireDiameter = InGWireShellThickIn + InGWireShellThickOut;
  double inGWradii = 0.5 * inGWireDiameter;
  double fakeLayerInIWthick = -0.0001 + GasInnerWallThick;
  double inner_radius_0 = inner_radius + envelop_Inner_thickness - fakeLayerInIWthick;

  double radius_layer_0 = inner_radius + envelop_Inner_thickness + FWradii + secure + capGasLayer;
  double radius_layerOut_0 = radius_layer_0 - FWradii - secure;

  double drop = 0.0;
  double radius_layer = 0.0;
  double radius_layerIn_0 = 0.0;
  double radius_layerIn = 0.0;
  double radius_layerOut = 0.0;
  double epsilonIn = 0.0;
  double epsilonOut = 0.0;
  double layerangle = 0.0;
  double cellBase = 0.0;
  double inscribedRadius = 0.0;
  double circumscribedRadius = 0.0;
  double zlength = 0.0;
  double cellStaggering = 0.0;
  double epsilonInGwRing = 0.0;
  double epsilonOutGwRing = 0.0;
  double radius_layerIn_whole_cell = 0.0;
  double epsilonIn_whole_cell = 0.0;

  //------------------------------------------------------------------------
  // The enlarge parameter is used to see the wires in the rendering
  //------------------------------------------------------------------------

  double enlarge = 50.;

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
                                 outer_radius - Carbon2OuterWallThick - FoamOuterWallThick, halflength);// FIXME there is an overlap with OuterWall_Carbon1 and the last guard wire layer
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
  // nSuperLayer = 1;

  for (int SL = 0; SL < nSuperLayer; ++SL) {

    num_wire = nSWire + SL * nSDeltaWire;
    phi = 2. * TMath::Pi() / num_wire;
    nFWire = nHorizontalFWire * num_wire;
    phi_layer = 2. * TMath::Pi() / nFWire;
    nFWire1 = nFWire / 2;
    if (ceilf(nFWire1) != nFWire1)
      throw std::runtime_error("Error: Failed to build CDCH. Please make sure that '(nStoFWireRatio - nVerticalFWire) * (nSWire + SuperLayerIndex * nSDeltaWire)' is always an even number");
    phi_layer1 = 2.0 * phi_layer;
    scaleFactor = (1.0 + TMath::Pi() / num_wire) / (1.0 - TMath::Pi() / num_wire);
    dropFactor = (1.0 / cos(halfalpha) - 1.0); // used to determine the radius of the hyperboloid in z = +- halflength with r_out = r_min + r_min * dropFactor
    epsilonFactor = sin(halfalpha) / halflength;
    layerangle = -0.5 * phi;

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

      dd4hep::Hyperboloid HypeLayer0(inner_radius_0, 0.0, radius_layerOut_0 - secure, stereoOut0, halflength);
      lvLayerVol.push_back(dd4hep::Volume("hyperboloid_inner_guard_layer", HypeLayer0, description.material("GasHe_90Isob_10")));
      lvLayerVol.back().setVisAttributes(description, "vCDCH:Pb");

      epsilonInGwRing = atan(inGuardRad * (1.0 + dropFactor) * epsilonFactor);
      zlength = halflength;
      zlength -= sin(epsilonInGwRing) * inGWradii;
      zlength /= cos(epsilonInGwRing);

      guard_wires.mother_volume = lvLayerVol.back();
      guard_wires.type = "G";
      guard_wires.num = nFWire1; //(#guard wires == # field wires)
      guard_wires.radius = inGuardRad - inGWradii;
      guard_wires.phi = 2. * TMath::Pi() / guard_wires.num;
      guard_wires.phioffset = layerangle;
      guard_wires.stereo = epsilonInGwRing;
      guard_wires.thickness = 0.5 * InGWireShellThickIn * enlarge;  // half the inner thickness as radius of tube
      guard_wires.halflength = zlength;
      guard_wires.name = string("Gwire_inner_stereoplus");

      dd4hep::Tube Gwire(0.0, guard_wires.thickness, halflength);
      lvGwireVol.push_back(dd4hep::Volume("Gwire_inner", Gwire, description.material("G4_Al")));
      lvGwireVol.back().setVisAttributes(description, wirecol);

      guard_wires.volume = lvGwireVol.back();
      CDCHBuild::PlaceWires(guard_wires, FWireShellThickOut, halflength, SL, -1);

      guard_wires.volume = lvGwireVol.back(); // needed because applyWireCoating acts on it
      guard_wires.radius = inGuardRad + inGWradii + extShiftFW;
      guard_wires.phioffset = layerangle + phi_layer;
      guard_wires.stereo = -1.0 * epsilonInGwRing;
      guard_wires.name = string("Gwire_inner_stereominus");
      CDCHBuild::PlaceWires(guard_wires, FWireShellThickOut, halflength, SL, -1);

      drop = radius_layer_0 * dropFactor;
      radius_layer = radius_layer_0 + drop;
      epsilon = atan(radius_layer * epsilonFactor);
      radius_layerIn_0 = radius_layer_0 - FWradii - 2.0 * secure;
      radius_layerIn = radius_layerIn_0 + drop;
      radius_layerOut_0 = radius_layer_0 + FWradii;
      radius_layerOut = radius_layerOut_0 + drop;
      epsilonIn = atan(sqrt(pow(radius_layerIn, 2) - pow(radius_layerIn_0, 2)) / halflength);
      epsilonOut = atan(sqrt(pow(radius_layerOut, 2) - pow(radius_layerOut_0, 2)) / halflength);

      dd4hep::Hyperboloid HypeLayer1(radius_layerIn_0, epsilonIn, radius_layerOut_0, epsilonOut, halflength);
      lvLayerVol.push_back(dd4hep::Volume("hyperboloid_first_field_wire_ring", HypeLayer1, description.material("GasHe_90Isob_10")));
      lvLayerVol.back().setVisAttributes(description, "vCDCH:Plastic");

      zlength = halflength;
      zlength -= sin(epsilon) * FWradii;
      zlength /= cos(epsilon);

      field_wires_top.mother_volume = lvLayerVol.back();
      field_wires_top.type = "F";
      field_wires_top.num = nFWire1;
      field_wires_top.radius = radius_layerIn_0 - FWradii - extShiftFW;
      field_wires_top.phi = 2. * TMath::Pi() /field_wires_top.num;;
      field_wires_top.phioffset = layerangle + cellStaggering - phi_layer;
      field_wires_top.stereo = sign_epsilon * epsilon;
      field_wires_top.thickness = 0.5 * FWireShellThickIn * enlarge;
      field_wires_top.halflength = zlength;

      lvFwireName = dd4hep::_toString(SL, "lvFwire_%d_init");
      field_wires_top.name = string(lvFwireName);

      dd4hep::Tube Fwire(0.0, field_wires_top.thickness, halflength);
      lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material("G4_Al")));
      lvFwireVol.back().setVisAttributes(description, wirecol);

      field_wires_top.volume = lvFwireVol.back();
      CDCHBuild::PlaceWires(field_wires_top, FWireShellThickOut, halflength, SL, -1);

      radius_layer_0 += FWradii;
      
      
      // the commented part below is a trial to get totally rid of PlaceWire function and to treat wires outside of the sensitive volume the same way as the one inside

      //// deal with the first guard layer
      //double stereoOut0 = atan((radius_layerOut_0 + FWradii) * (1.0 * dropFactor * epsilonFactor));
      ////stereoOut0 = atan(sqrt(diff_of_squares(radius_layerOut_0 + FWradii,

      ////radius_layerOut_0 = radius_layerIn_0 + FWireDiameter + 2.0 * secure;
      ////radius_layerOut = radius_layerOut_0 + drop;
      ////epsilonOut = atan(sqrt(diff_of_squares(radius_layerOut, radius_layerOut_0)) / halflength);
      ////dd4hep::Hyperboloid HypeLayer0(inner_radius_0, 0.0, radius_layerOut_0 + FWradii - secure, stereoOut0, halflength);
      //dd4hep::Hyperboloid HypeLayer0(inner_radius_0, 0.0, radius_layerOut_0 + FWradii, stereoOut0, halflength);
      //dd4hep::Volume inner_guard_layer_volume("hyperboloid_inner_guard_layer", HypeLayer0, description.material("GasHe_90Isob_10"));
      //dd4hep::PlacedVolume inner_guard_layer_placedVolume = parentVol.placeVolume(inner_guard_layer_volume);
      //CDCHDetector.setPlacement(inner_guard_layer_placedVolume);

      //epsilonInGwRing = atan(inGuardRad * (1.0 + dropFactor) * epsilonFactor);
      //zlength = halflength;
      //zlength -= sin(epsilonInGwRing) * inGWradii;
      //zlength /= cos(epsilonInGwRing);

      //guard_wires.mother_volume = inner_guard_layer_volume;
      //guard_wires.type = "G";
      //guard_wires.num = nInGWire;
      //guard_wires.radius = inGuardRad - inGWradii;
      //guard_wires.phi = 2. * TMath::Pi() / guard_wires.num;
      //guard_wires.phioffset = layerangle;
      //guard_wires.stereo = epsilonInGwRing;
      //guard_wires.thickness = 0.5 * InGWireShellThickIn * enlarge;  // half the inner thickness as radius of tube
      //guard_wires.halflength = zlength;
      //guard_wires.name = string("InnerGuardWire") + dd4hep::_toString(guard_wires.stereo, "_stereo%f");
      //dd4hep::Tube Gwire(0.0, guard_wires.thickness, halflength);
      //guard_wires.volume = dd4hep::Volume("InnerGuardWire", Gwire, description.material("G4_Al"));
      //apply_wire_coating(guard_wires, InGWireShellThickIn, halflength);
      //// Radial translation
      //dd4hep::Translation3D radial_translation_guard_wire(guard_wires.radius, 0., 0.);
      //// stereo rotation
      //dd4hep::RotationX rot_stereo_guard_wire(guard_wires.stereo);
      //dd4hep::PlacedVolume inner_guard_wire_placedvolume;
      //for (int phi_index = 0; phi_index < guard_wires.num; phi_index++) {
      //  dd4hep::RotationZ iRot(guard_wires.phioffset + guard_wires.phi * phi_index);
      //  dd4hep::Transform3D total_transformation(iRot * radial_translation_guard_wire * rot_stereo_guard_wire);
      //  inner_guard_wire_placedvolume = inner_guard_layer_volume.placeVolume(guard_wires.volume, total_transformation);
      //  //if(setWireSensitive)
      //  //  inner_guard_wire_placedvolume.addPhysVolID("phi", phi_index).addPhysVolID("hitorigin", 3).addPhysVolID("stereo", guard_wires.stereo > 0 ? 0 : 1).addPhysVolID("layerInCell", 0);
      //}

      //guard_wires.radius = inGuardRad + inGWradii + extShiftFW;
      //guard_wires.phioffset = layerangle + phi_layer;
      //guard_wires.stereo = -1.0 * epsilonInGwRing;
      //guard_wires.name = string("InnerGuardWire") + dd4hep::_toString(guard_wires.stereo, "_stereo%f");

      //drop = radius_layer_0 * dropFactor;
      //radius_layer = radius_layer_0 + drop;
      //epsilon = atan(radius_layer * epsilonFactor);
      //radius_layerIn_0 = radius_layer_0 - FWradii - 2.0 * secure;
      //radius_layerIn = radius_layerIn_0 + drop;
      //radius_layerOut_0 = radius_layer_0 + FWradii;
      //radius_layerOut = radius_layerOut_0 + drop;
      //epsilonIn = atan(sqrt(pow(radius_layerIn, 2) - pow(radius_layerIn_0, 2)) / halflength);
      //epsilonOut = atan(sqrt(pow(radius_layerOut, 2) - pow(radius_layerOut_0, 2)) / halflength);

      //dd4hep::Hyperboloid HypeLayer1(radius_layerIn_0, epsilonIn, radius_layerOut_0, epsilonOut, halflength);
      //lvLayerVol.push_back(dd4hep::Volume("lvLayer_0", HypeLayer1, description.material("GasHe_90Isob_10")));
      //lvLayerVol.back().setVisAttributes(description, "vCDCH:Plastic");

      //zlength = halflength;
      //zlength -= sin(epsilon) * FWradii;
      //zlength /= cos(epsilon);

      //field_wires_top.mother_volume = lvLayerVol.back();
      //field_wires_top.type = "F";
      //field_wires_top.num = nFWire1;
      //field_wires_top.radius = radius_layerIn_0 - FWradii - extShiftFW;
      //field_wires_top.phi = phi_layer1;
      //field_wires_top.phioffset = layerangle + cellStaggering - phi_layer;
      //field_wires_top.stereo = sign_epsilon * epsilon;
      //field_wires_top.thickness = 0.5 * FWireShellThickIn * enlarge;
      //field_wires_top.halflength = zlength;
      //field_wires_top.name = dd4hep::_toString(SL, "Wire_SL%d") + dd4hep::_toString(-1, "_layer%d") + string("_type") + field_wires_top.type + dd4hep::_toString(field_wires_top.stereo, "_stereo%f_top");

      //lvFwireName = dd4hep::_toString(SL, "lvFwire_%d_init");

      //dd4hep::Tube Fwire(0.0, field_wires_top.thickness, halflength);
      //lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material("G4_Al")));
      //lvFwireVol.back().setVisAttributes(description, wirecol);

      ////field_wires_top.volume = lvFwireVol.back();
      //////CDCHBuild::PlaceWires(field_wires_top, FWireShellThickOut, halflength, SL, -1);

      //radius_layer_0 += FWradii;

    } else {

      delta_radius_layer = 2. * TMath::Pi() * radius_layerOut_0 / (num_wire - TMath::Pi());
    }

    //------------------------------------------------------------------------
    // Starting the layer loop. nLayer=8
    //------------------------------------------------------------------------

    for (int ilayer = 0; ilayer < nLayer; ilayer++) {
      //------------------------------------------------------------------------
      // Fill the geometry parameters of the layer. Each layer lies
      // on top of the following one, so new layerIn = old layerOut
      //------------------------------------------------------------------------

      inscribedRadius = 0.5 * delta_radius_layer;
      circumscribedRadius = inscribedRadius * sqrt(2.0);
      senseWireRing_radius_0 = radius_layer_0 + inscribedRadius;
      sign_epsilon *= -1;

      radius_layerIn_0 = radius_layerOut_0;
      radius_layerIn = radius_layerOut;
      epsilonIn = epsilonOut;

      radius_layerOut_0 = radius_layerIn_0 + FWireDiameter + 2.0 * secure;
      radius_layerOut = radius_layerOut_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(radius_layerOut, radius_layerOut_0)) / halflength);

      zlength = halflength;

      // save bottom layer inner radius and epsilon before they are modified to build the whole layer hyperboloid volume 
      radius_layerIn_whole_cell = radius_layerIn_0;
      epsilonIn_whole_cell = epsilonIn;

      //------------------------------------------------------------------------
      // Reduce zlength to avoid volume extrusions and check the staggering
      //------------------------------------------------------------------------

      zlength -= sin(epsilon) * FWradii;
      zlength /= cos(epsilon);

      if (ilayer % 2 == 1)
        cellStaggering = phi_layer;
      else
        cellStaggering = 0.0;

      //------------------------------------------------------------------------
      // Fill the field wire struct with all the relevant information
      // This is the bottom of the cell.
      //------------------------------------------------------------------------

      field_wires_bottom.type = "F";
      field_wires_bottom.num = nFWire1;
      field_wires_bottom.radius = radius_layerIn_0 + FWradii + extShiftFW;
      field_wires_bottom.phi = phi_layer1;
      field_wires_bottom.phioffset = layerangle + cellStaggering;
      field_wires_bottom.stereo = sign_epsilon * epsilon;
      field_wires_bottom.thickness = 0.5 * FWireShellThickIn * enlarge;
      field_wires_bottom.halflength = zlength;
      field_wires_bottom.name = dd4hep::_toString(SL, "Wire_SL%d") + dd4hep::_toString(ilayer, "_layer%d") + string("_type") + field_wires_bottom.type + dd4hep::_toString(field_wires_bottom.stereo, "_stereo%f_bottom");

      //------------------------------------------------------------------------
      // Define the field wire name and build the field wire volume
      //------------------------------------------------------------------------

      lvFwireName = dd4hep::_toString(SL, "lvFwire_%d") + dd4hep::_toString(ilayer, "_%d");

      dd4hep::Tube Fwire(0.0, field_wires_bottom.thickness, halflength);
      lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material("G4_Al")));
      lvFwireVol.back().setVisAttributes(description, wirecol);

      //------------------------------------------------------------------------
      // Add the field wire volume to the struct 
      //------------------------------------------------------------------------

      field_wires_bottom.volume = lvFwireVol.back();
      apply_wire_coating(field_wires_bottom, FWireShellThickOut, halflength);

      //------------------------------------------------------------------------
      // Next, fill the geometry parameters of the central layer.
      //------------------------------------------------------------------------

      iradius = radius_layer_0;
      radius_layer_0 += delta_radius_layer;
      drop = radius_layer_0 * dropFactor;

      radius_layerIn_0 = radius_layerOut_0;
      radius_layerIn = radius_layerOut;
      epsilonIn = epsilonOut;
      radius_layerOut_0 = radius_layer_0 - FWireDiameter - 2.0 * secure;
      radius_layerOut = radius_layerOut_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(radius_layerOut, radius_layerOut_0)) / halflength);
      zlength = halflength;

      //------------------------------------------------------------------------
      // Reduce zlength to avoid volume extrusions
      //------------------------------------------------------------------------

      zlength -= sin(epsilon) * 0.5 * FWireShellThickIn;
      zlength /= cos(epsilon);

      //------------------------------------------------------------------------
      // Fill the sense wire struct with all the relevant information
      //------------------------------------------------------------------------

      sense_wires.type = "S";
      sense_wires.num = num_wire;
      sense_wires.radius = senseWireRing_radius_0;
      sense_wires.phi = phi;
      sense_wires.phioffset = cellStaggering;
      sense_wires.stereo = sign_epsilon * epsilon;
      sense_wires.thickness = 0.5 * SWireShellThickIn * enlarge;
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
      // Tune the radius and epsilon of the central field wires
      //------------------------------------------------------------------------

      idelta_radius = delta_radius_layer * 0.5;
      iradius += idelta_radius;
      epsilon = atan(iradius * (1.0 + dropFactor) * epsilonFactor);

      //------------------------------------------------------------------------
      // Fill the central field wire struct with all the relevant information
      //------------------------------------------------------------------------

      field_wires_center.type = "F";
      field_wires_center.num = num_wire;
      field_wires_center.radius = iradius;
      field_wires_center.phi = phi;
      field_wires_center.phioffset = layerangle + cellStaggering;
      field_wires_center.stereo = sign_epsilon * epsilon;
      field_wires_center.thickness = 0.5 * centerFWireShellThickIn * enlarge;
      field_wires_center.halflength = zlength;
      field_wires_center.volume = lvFwireVol.back();
      field_wires_center.name = dd4hep::_toString(SL, "Wire_SL%d") + dd4hep::_toString(ilayer, "_layer%d") + string("_type") + field_wires_center.type + dd4hep::_toString(field_wires_center.stereo, "_stereo%f_center");
      apply_wire_coating(field_wires_center, centerFWireShellThickOut, halflength);


      //------------------------------------------------------------------------
      // Next, fill the geometry parameters of the upper layer.
      //------------------------------------------------------------------------

      radius_layerIn_0 = radius_layerOut_0;
      radius_layerIn = radius_layerOut;
      epsilonIn = epsilonOut;
      radius_layerOut_0 = radius_layerIn_0 + FWireDiameter + 2.0 * secure;
      radius_layerOut = radius_layerOut_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(radius_layerOut, radius_layerOut_0)) / halflength);
      zlength = halflength;

      // Create hyperboloid volume of the whole layer for the cell sensitive volume definition, not needed per se but helps having a well balanced volume tree (having too many volumes inside one mother volume harms performance)
      wholeHyperboloidVolumeName = dd4hep::_toString(SL, "hyperboloid_SL_%d") + dd4hep::_toString(ilayer, "_layer_%d");
      dd4hep::Hyperboloid whole_layer_hyperboloid = dd4hep::Hyperboloid(radius_layerIn_whole_cell, epsilonIn_whole_cell, radius_layerOut_0, epsilonOut, zlength);
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
      field_wires_top.num = nFWire1;
      field_wires_top.radius = radius_layerIn_0 - FWradii - extShiftFW;
      field_wires_top.phi = phi_layer1;
      field_wires_top.phioffset = layerangle + cellStaggering;
      field_wires_top.stereo = sign_epsilon * epsilon;
      field_wires_top.thickness = 0.5 * FWireShellThickIn * enlarge;
      field_wires_top.halflength = zlength;
      field_wires_top.volume = lvFwireVol.back();
      field_wires_top.name = dd4hep::_toString(SL, "Wire_SL%d") + dd4hep::_toString(ilayer, "_layer%d") + string("_type") + field_wires_top.type + dd4hep::_toString(field_wires_top.stereo, "_stereo%f_top");
      apply_wire_coating(field_wires_top, FWireShellThickOut, halflength);

      //if(setWireSensitive){
      //  field_wires_bottom.volume.setSensitiveDetector(sens_det);
      //  field_wires_center.volume.setSensitiveDetector(sens_det);
      //  sense_wires.volume.setSensitiveDetector(sens_det);
      //  field_wires_top.volume.setSensitiveDetector(sens_det);
      //}

      // Create the tube segment volume to identify the cell sensitive regions (FIXME could be improved with e.g. extruded valumes or tessalatedSolids)
      dd4hep::Tube cellID_tube_segment(radius_layerIn_whole_cell, radius_layerOut_0, halflength, - sense_wires.phi / 2.0, sense_wires.phi / 2.0);
      
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
        double phi_angle_sense_wire_rotation = sense_wires.phioffset + sense_wires.phi * phi_index;
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
          dd4hep::RotationZ rot_phi_bottom_wire((field_wires_bottom.phioffset + field_wires_bottom.phi * sub_phi_index) - phi_angle_sense_wire_rotation);
          dd4hep::PlacedVolume field_wire_bottom_placedvolume = cellID_volume.placeVolume(field_wires_bottom.volume, dd4hep::Transform3D(rot_phi_bottom_wire * dd4hep::Translation3D(field_wires_bottom.radius, 0., 0.) * dd4hep::RotationX(field_wires_bottom.stereo - sense_wires.stereo)));
          //if(setWireSensitive)
          //  field_wire_bottom_placedvolume.addPhysVolID("phi", sub_phi_index).addPhysVolID("hitorigin", 2).addPhysVolID("stereo", field_wires_center.stereo > 0 ? 0 : 1).addPhysVolID("layerInCell", 1);
        }

        // central field wires
        for(int sub_phi_index = phi_index * middle_to_middle_num_wire_ratio; sub_phi_index < (phi_index * middle_to_middle_num_wire_ratio) + middle_to_middle_num_wire_ratio; sub_phi_index++){
          // phi rotation
          dd4hep::RotationZ rot_phi_center_wire((field_wires_center.phioffset + field_wires_center.phi * sub_phi_index) - phi_angle_sense_wire_rotation);
          dd4hep::PlacedVolume field_wire_center_placedvolume = cellID_volume.placeVolume(field_wires_center.volume, dd4hep::Transform3D(rot_phi_center_wire * dd4hep::Translation3D(field_wires_center.radius, 0., 0.) * dd4hep::RotationX(field_wires_center.stereo - sense_wires.stereo)));
          //if(setWireSensitive)
          //  field_wire_center_placedvolume.addPhysVolID("phi", sub_phi_index).addPhysVolID("hitorigin", 2).addPhysVolID("stereo", field_wires_center.stereo > 0 ? 0 : 1).addPhysVolID("layerInCell", 2);
        }

        // top field wires
        for(int sub_phi_index = phi_index * middle_to_top_num_wire_ratio; sub_phi_index < (phi_index * middle_to_top_num_wire_ratio) + middle_to_top_num_wire_ratio; sub_phi_index++){
          // phi rotation
          dd4hep::RotationZ rot_phi_top_wire((field_wires_top.phioffset + field_wires_top.phi * sub_phi_index) - phi_angle_sense_wire_rotation);
          dd4hep::PlacedVolume field_wire_top_placedvolume = cellID_volume.placeVolume(field_wires_top.volume, dd4hep::Transform3D(rot_phi_top_wire * dd4hep::Translation3D(field_wires_top.radius, 0., 0.) * dd4hep::RotationX(field_wires_top.stereo - sense_wires.stereo)));
          //if(setWireSensitive)
          //  field_wire_top_placedvolume.addPhysVolID("phi", sub_phi_index).addPhysVolID("hitorigin", 2).addPhysVolID("stereo", field_wires_center.stereo > 0 ? 0 : 1).addPhysVolID("layerInCell", 3);
        }
      }

      //------------------------------------------------------------------------
      // Scale the delta radius of the layer for next iteration
      //------------------------------------------------------------------------

      delta_radius_layer *= scaleFactor;
    }

    if (SL == (nSuperLayer - 1)) {// the last super layer is special since we need to add the field wires outside of the sensitive volume and the guard wires

      // Take care of the field wires outside of the sensitive volume
      radius_layerIn_0 = radius_layerOut_0;
      radius_layerIn = radius_layerOut;
      epsilonIn = epsilonOut;
      radius_layerOut_0 = radius_layer_0 + FWireDiameter + 2.0 * secure;
      radius_layerOut = radius_layerOut_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(radius_layerOut, radius_layerOut_0)) / halflength);

      dd4hep::Hyperboloid HypeLayerOut(radius_layerIn_0, epsilonIn, radius_layerOut_0, epsilonOut, halflength);
      lvLayerVol.push_back(dd4hep::Volume("hyperboloid_last_field_wire_ring", HypeLayerOut, description.material("GasHe_90Isob_10")));
      lvLayerVol.back().setVisAttributes(description, "vCDCH:Plastic");

      zlength = halflength;
      zlength -= sin(epsilon) * FWradii;
      zlength /= cos(epsilon);

      field_wires_bottom.mother_volume = lvLayerVol.back();
      field_wires_bottom.type = "F";
      field_wires_bottom.num = nFWire1;
      field_wires_bottom.radius = radius_layerIn_0 + FWradii + extShiftFW;
      field_wires_bottom.phi = 2. * TMath::Pi() / field_wires_bottom.num;
      field_wires_bottom.phioffset = layerangle + cellStaggering + phi_layer;
      field_wires_bottom.stereo = -1. * sign_epsilon * epsilon;
      field_wires_bottom.thickness = 0.5 * FWireShellThickIn * enlarge;
      field_wires_bottom.halflength = zlength;

      lvFwireName = dd4hep::_toString(SL, "lvFwire_%d_out");
      field_wires_bottom.name = lvFwireName;

      dd4hep::Tube Fwire(0.0, field_wires_bottom.thickness, halflength);
      lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material("G4_Al")));
      lvFwireVol.back().setVisAttributes(description, wirecol);

      field_wires_bottom.volume = lvFwireVol.back();
      CDCHBuild::PlaceWires(field_wires_bottom, FWireShellThickOut, halflength, SL, -1);

      //------------------------------------------------------------------------
      // Start placing the outer layer of guard wires (#guard wires == # field wires)
      //------------------------------------------------------------------------

      radius_layerIn_0 = radius_layerOut_0;
      radius_layerIn = radius_layerOut;
      epsilonIn = epsilonOut;
      radius_layerOut_0 = radius_layer_0 + FWireDiameter + 2.0 * secure;
      radius_layerOut = radius_layerOut_0 + drop;
      epsilonOut = atan(sqrt(diff_of_squares(radius_layerOut, radius_layerOut_0)) / halflength);

      dd4hep::Hyperboloid HypeLayerOutG(radius_layerIn_0, epsilonOut, outer_radius - envelop_Outer_thickness - 0.0001,
                                        0.0, halflength);
      lvLayerVol.push_back(dd4hep::Volume("hyperboloid_outer_guard_layer", HypeLayerOutG, description.material("GasHe_90Isob_10")));
      lvLayerVol.back().setVisAttributes(description, "vCDCH:Pb");

      epsilonOutGwRing = atan(outGuardRad * (1.0 + dropFactor) * epsilonFactor);
      zlength = halflength;
      zlength -= sin(epsilonOutGwRing) * inGWradii;
      zlength /= cos(epsilonOutGwRing);

      guard_wires.mother_volume = lvLayerVol.back();
      guard_wires.type = "G";
      guard_wires.num = nFWire1;
      guard_wires.radius = outGuardRad - inGWradii;
      guard_wires.phi = 2. * TMath::Pi() / guard_wires.num;
      guard_wires.phioffset = layerangle;
      guard_wires.stereo = epsilonOutGwRing;
      guard_wires.thickness = 0.5 * OutGWireShellThickIn * enlarge;
      guard_wires.halflength = zlength;
      guard_wires.name = string("Gwire_outer_stereominus");

      dd4hep::Tube Gwire(0.0, guard_wires.thickness, halflength);
      lvGwireVol.push_back(dd4hep::Volume("Gwire_outer", Gwire, description.material("G4_Al")));
      lvGwireVol.back().setVisAttributes(description, wirecol);

      guard_wires.volume = lvGwireVol.back();
      CDCHBuild::PlaceWires(guard_wires, FWireShellThickOut, halflength, SL, -1);

      guard_wires.volume = lvGwireVol.back(); // needed because applyWireCoating acts on it
      guard_wires.radius = outGuardRad + inGWradii + extShiftFW;
      guard_wires.phioffset = layerangle + phi_layer;
      guard_wires.stereo = -1.0 * epsilonOutGwRing;
      guard_wires.name = string("Gwire_outer_stereoplus");
      CDCHBuild::PlaceWires(guard_wires, FWireShellThickOut, halflength, SL, -1);
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
  dd4hep::Tube CDCH_envelope(dd4hep::_toDouble("CDCH:r0"), dd4hep::_toDouble("CDCH:rOut"), dd4hep::_toDouble("CDCH:zHalfExtentWithServices"));

  dd4hep::Volume envelope("lvCDCH", CDCH_envelope, description.air());
  envelope.setVisAttributes(description, "vCDCH:Air");

  // ******************************************************
  // Build CDCH cable
  // ******************************************************

  builder.build_layer(CDCH_det, envelope, sens_det);

  dd4hep::printout(dd4hep::DEBUG, "CreateCDCH", "MotherVolume is: %s", envelope.name());
  sens_det.setType("tracker");

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


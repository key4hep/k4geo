//**************************************************************************
// \file DREndcapTubes.cpp
// \brief: Implementation of the endcap geometry of the IDEA dual-readout
//         calorimeter using the capillary tubes technology
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 July 2024
//**************************************************************************

// Includers from DD4hep
#include "DDRec/Vector3D.h"
#include <DD4hep/DetFactoryHelper.h>

// Includers from stl
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace dd4hep;
using namespace dd4hep::rec;  // for dd4hep::rec::Vector3D

// Includers from project files
#include "DREndcapTubes.hh"

#define COUNTTUBES  // if defined it counts the tubes in each tower
// #define DUMPTOWEREDGES // if defined it prints tower edges x and y position (cm)

// This methods unifies in 32 bits two IDs of 16 bits max
// It is used to create a G4 copynumber for volumes with two IDs
unsigned int CalcCpNo32bits(int id1, int id2)
{
  if (id1 > 0xFFF || id2 > 0xFFF){
    throw std::invalid_argument("row and col IDs should not take more than 16 bits");
  }
  unsigned int CpNo = (id1 << 16) | id2; 
  return CpNo;
}

// This methods unifies in 32 bits two IDs of 1 bit only
unsigned int CalcCpNo2bits(int id1, int id2)
{
  if (id1 > 1 || id2 > 1) {
    throw std::invalid_argument("core and cherenkov ID must be 0 or 1"); 
  } 
  unsigned int CpNo = (id1 << 1) | id2;
  return CpNo;
}

// Create the endcap calorimeter
//
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  std::cout << "--> DREndcapTubes::create_detector() start" << std::endl;

  // Get info from the xml file
  //
  // Info for the main detector element DREndcapTube
  sens.setType("calorimeter");
  xml_det_t x_det = e;
  std::string det_name = x_det.nameStr();  // DREndcapTube volume
  std::cout << "--> DREndcapTubes: going to create " << det_name << ", with ID: " << x_det.id()
            << std::endl;
  xml_dim_t x_dim = x_det.dimensions();
  const double innerR = x_dim.inner_radius();  // inner radius at theta = 0. rad
  const double tower_height = x_dim.z_length();  // tower height/length
  const double length = tower_height;  // alias for tower height/length
  const int NbOfZRot = static_cast<int>(x_dim.deltaphi());  // number or rations aroung Z axis
  const double phi_unit = 2 * M_PI / (double)NbOfZRot;  // slice delta phi
  std::cout << "--> DREndcapTubes: from XML description: innerR " << innerR / m
            << " m, tower length " << tower_height / m << " m, z-rotations " << NbOfZRot;
  // Info for subdetectors
  xml_det_t x_assembly = x_det.child(_Unicode(assembly));
  xml_det_t x_stave = x_det.child(_Unicode(stave));
  xml_det_t x_tower = x_det.child(_Unicode(tower));
  xml_det_t x_capillary_S = x_det.child(_Unicode(tube_S));
  xml_det_t x_capillary_C = x_det.child(_Unicode(tube_C));
  xml_det_t x_clad_S = x_det.child(_Unicode(clad_S));
  xml_det_t x_clad_C = x_det.child(_Unicode(clad_C));
  auto x_clad_C_sens = x_clad_C.isSensitive();
  xml_det_t x_core_S = x_det.child(_Unicode(core_S));
  auto x_core_S_sens = x_core_S.isSensitive();
  xml_det_t x_core_C = x_det.child(_Unicode(core_C));
  auto x_core_C_sens = x_core_C.isSensitive();
  const double tubeRadius = x_capillary_S.outer_radius();
  std::cout << ", tube radius " << tubeRadius / mm << " mm" << std::endl;
  const double cladRadius = x_clad_S.outer_radius();
  const double coreRadius = x_core_S.outer_radius();

  // Hardcoded parameters (not supposed to be changed)
  constexpr double thetaB{M_PI / 4};  // theta at the end of barrel (45 deg)
  constexpr int NbOfEndcap{40};  // number of eta towers for the endcap.
                                 // De facto I will stop at 35 to leave
                                 // room for the beam pipe.
  constexpr int NbOfEndcapReduced{35};
  // Length correction to be applied to capillary tubes
  // with length reduced by the intersection
  // with the tower's planes.
  // It is needed for the intersections with
  // the lateral sides for which the intersecting
  // point is not exactly the left or right edge of the tube
  constexpr double TubeLengthOffset{0.1 * mm};
  const double tubeDiameter = tubeRadius * 2.;
  const double y_pitch = tubeDiameter * sqrt(3) / 2.;  // y-distance from tube center to center
                                                       // fixed by the tubes positioning
  // The 8 vertices of a tower to be defined later
  Vector3D pt[8];
  // Create the geometry helper and set paramters
  DREndcapTubeHelper Helper;
  Helper.SetInnerR(innerR);
  Helper.SetTowerHeight(tower_height);
  Helper.SetNumZRot(NbOfZRot);
  Helper.SetTubeRadius(tubeRadius);

  // Start building the geometry
  //

  // The phi-slices (staves) will be placed inside an assembly volume
  Assembly AssemblyEndcap("DREndcapTubes");
  AssemblyEndcap.setVisAttributes(description, x_assembly.visStr());

  // Volume that contains a slice of the right endcap
  // I use the EightPointSolid/Arb8/G4Generictrap so I define directly its 8 points
  //
  // Distance between z-axis and the starting point of the Air stave containing
  // the towers, it is included to avoid overlaps with the compensating solenoid
  const double DistancetoSolenoid = 22*cm;
  // The first two points of the inner face are defined by DistancetoSolenoid
  // similarly to the second two points
  // The second two points of the inner face are at (x=tan(0.5*phi_unit, y=innerR)
  // x with plus or minus sign
  double vertices[16];
  vertices[0] = static_cast<double>(-DistancetoSolenoid * tan(0.5 * phi_unit));
  vertices[1] = static_cast<double>(DistancetoSolenoid);
  vertices[2] = static_cast<double>(-innerR * tan(0.5 * phi_unit));
  vertices[3] = static_cast<double>(innerR);
  vertices[4] = static_cast<double>(innerR * tan(0.5 * phi_unit));
  vertices[5] = static_cast<double>(innerR);
  vertices[6] = static_cast<double>(DistancetoSolenoid * tan(0.5 * phi_unit));
  vertices[7] = static_cast<double>(DistancetoSolenoid);
  // The first two points of the outer face are at the same distance to the z-axis
  // as in the inner face
  // The second two poits of the outer face are same as before with innerR+tower_height
  vertices[8] = static_cast<double>(-DistancetoSolenoid * tan(0.5 * phi_unit));
  vertices[9] = static_cast<double>(DistancetoSolenoid);
  vertices[10] = static_cast<double>(-(innerR + tower_height) * tan(0.5 * phi_unit));
  vertices[11] = static_cast<double>(innerR + tower_height);
  vertices[12] = static_cast<double>((innerR + tower_height) * tan(0.5 * phi_unit));
  vertices[13] = static_cast<double>(innerR + tower_height);
  vertices[14] = static_cast<double>(DistancetoSolenoid * tan(0.5 * phi_unit));
  vertices[15] = static_cast<double>(DistancetoSolenoid);
  // Equivalent of Geant4 GenericTrap shape constructor
  EightPointSolid phiER("phiER", tower_height / 2., vertices);
  Volume phiERLog("phiER", phiER, description.material(x_stave.attr<std::string>(_U(material))));
  phiERLog.setVisAttributes(description, x_stave.visStr());

  // Place logical volume containing the Endcap R slice multiple times
  //
  // Rotate the endcap R slice around the Z axis
  for (std::size_t j = 0; j < static_cast<std::size_t>(NbOfZRot); j++) {
    // Placement of right endcap (positive z-coordinate)
    RotationZ rotz1(M_PI / 2.);  // here I discovered that dd4hep rotations around Z are inverted
    RotationZ rotz2(j * phi_unit);  // w.r.t. Geant4 (+->-)
    RotationX rotx(M_PI);
    RotationY roty(M_PI);
    Transform3D slice_trnsform(rotz1 * rotz2 * rotx * roty,
                               Position(0, 0, (innerR)*tan(thetaB) + length / 2.));
    PlacedVolume phiERPlaced = AssemblyEndcap.placeVolume(phiERLog, j, slice_trnsform);
    // ID this volume with one number, the rotation
    phiERPlaced.addPhysVolID("stave", j);

    // Placement of left endcap (negative z-coordinate)
    RotationY rotyleft(0.);
    Transform3D slice_trnsformleft(rotz1 * rotz2 * rotx * rotyleft,
                                   Position(0, 0, -1. * ((innerR)*tan(thetaB) + length / 2.)));
    PlacedVolume phiELPlaced = AssemblyEndcap.placeVolume(phiERLog, j+NbOfZRot, slice_trnsformleft);
    phiELPlaced.addPhysVolID("stave", j+NbOfZRot);
  }  // end of slice/stave placement

  // Create an S tube with full tower length
  Tube capillary_S(0. * mm, tubeRadius, tower_height / 2., 2 * M_PI);
  Volume capillary_SLog("capillary_SLog", capillary_S,
                        description.material(x_capillary_S.attr<std::string>(_U(material))));
  capillary_SLog.setVisAttributes(description, x_capillary_S.visStr());
  Tube clad_S(0. * mm, cladRadius, tower_height / 2., 2 * M_PI);  // cladding S
  Volume clad_SLog("clad_SLog", clad_S,
                   description.material(x_clad_S.attr<std::string>(_U(material))));
  clad_SLog.setVisAttributes(description, x_clad_S.visStr());
  PlacedVolume clad_SPlaced = capillary_SLog.placeVolume(clad_SLog, 1, Position(0., 0., 0.));
  // ID this volume with clad ID (0/1)
  clad_SPlaced.addPhysVolID("clad", 1);
  Tube core_S(0. * mm, coreRadius, tower_height / 2., 2 * M_PI);  // core S
  Volume core_SLog("DRETScoreS_Log", core_S,
                   description.material(x_core_S.attr<std::string>(_U(material))));
  core_SLog.setVisAttributes(description, x_core_S.visStr());
  if (x_core_S_sens) core_SLog.setSensitiveDetector(sens);
  PlacedVolume core_SPlaced = clad_SLog.placeVolume(core_SLog, CalcCpNo2bits(1,0));
  core_SPlaced.addPhysVolID("core", 1).addPhysVolID("cherenkov", 0);

  // Create a C tube with full tower length
  Tube capillary_C(0. * mm, tubeRadius, tower_height / 2., 2 * M_PI);
  Volume capillary_CLog("capillary_CLog", capillary_C,
                        description.material(x_capillary_C.attr<std::string>(_U(material))));
  capillary_CLog.setVisAttributes(description, x_capillary_C.visStr());
  Tube clad_C(0. * mm, cladRadius, tower_height / 2., 2 * M_PI);  // cladding C
  Volume clad_CLog("DRETSclad_CLog", clad_C,
                   description.material(x_clad_C.attr<std::string>(_U(material))));
  clad_CLog.setVisAttributes(description, x_clad_C.visStr());
  if (x_clad_C_sens) clad_CLog.setSensitiveDetector(sens);
  PlacedVolume clad_CPlaced = capillary_CLog.placeVolume(clad_CLog, 1, Position(0., 0., 0.));
  clad_CPlaced.addPhysVolID("clad", 1);
  Tube core_C(0. * mm, coreRadius, tower_height / 2., 2 * M_PI);  // core C
  Volume core_CLog("DRETScoreC_Log", core_C,
                   description.material(x_core_C.attr<std::string>(_U(material))));
  core_CLog.setVisAttributes(description, x_core_C.visStr());
  if (x_core_C_sens) core_CLog.setSensitiveDetector(sens);
  PlacedVolume core_CPlaced = clad_CLog.placeVolume(core_CLog, CalcCpNo2bits(1,1));
  core_CPlaced.addPhysVolID("core", 1).addPhysVolID("cherenkov", 1);

  // Build the towers inside and endcap R slice
  //
  Helper.SetRbool(1);  // Build the right endcap (will reflect it for left part)
  double thetaofcenter = 0.;  // theta angle to center of tower being constructed
  double thetaofcenter2 = 0.;  // theta angle to center of next tower
  // I fix delta theta of every tower to be the same for every tower
  double deltatheta_endcap[40] = {0.};
  for (int i = 0; i < NbOfEndcap; i++)
    deltatheta_endcap[i] = M_PI / 4 / (NbOfEndcap);
  double thetaE = 0.;
  for (int i = 0; i < NbOfEndcap; i++)
    thetaE += deltatheta_endcap[i];
  double fulltheta = thetaE;  // 45 deg by construction

#ifdef COUNTTUBES
  unsigned int totalTubeNo{0};
  unsigned int tubeNo{0};
  double totalTubeLength{0.};
#endif

  // Loop over tower number, stop at tower 35 to leave room for beam pipe
  for (int i = 0; i < NbOfEndcapReduced; i++) {
#ifdef COUNTTUBES
    tubeNo = 0;
#endif

    // Center of first (highest) tower is 45 deg - deltatheta_endcap[1]/2
    thetaofcenter = fulltheta - deltatheta_endcap[i] / 2.;
    // Center of the next tower
    thetaofcenter2 = thetaofcenter - deltatheta_endcap[i] / 2. - deltatheta_endcap[i + 1] / 2.;
    // Update Helper class parameters accordingly
    Helper.SetDeltaTheta(deltatheta_endcap[i]);
    Helper.SetThetaOfCenter(thetaofcenter);
    Helper.SetDeltaTheta2(deltatheta_endcap[i + 1]);
    Helper.SetThetaOfCenter2(thetaofcenter2);
    Helper.CalBasic();  // Perform internal calculations
    Helper.Getpt(pt);  // Update 8 Vectors defining the tower edges

#ifdef DUMPTOWEREDGES
    std::cout << "DREndcapTubes: Tower " << i << " edges (cm):" << std::endl;
    for (std::size_t edge = 0; edge < 8; edge++) {
      std::cout << "x " << pt[edge].x() << " y " << pt[edge].y() << std::endl;
    }
#endif

    // Create now the tower as a Trap
    //
    // Problem: G4Trap has a constructors using the 8 vertices of the trapezoid
    // but DD4hep Trap does not have such constructor.
    // Therefore I calculate here the parameters to be passed to the DD4hep Trap.
    auto pDz = (pt[7]).z();
    auto pDy1 = ((pt[2]).y() - (pt[1]).y()) * 0.5;
    auto pDx1 = ((pt[1]).x() - (pt[0]).x()) * 0.5;
    auto pDx2 = ((pt[3]).x() - (pt[2]).x()) * 0.5;
    auto fTalpha1 = ((pt[2]).x() + (pt[3]).x() - (pt[1]).x() - (pt[0]).x()) * 0.25 / pDy1;
    auto pAlp1 = std::atan(fTalpha1);
    auto pDy2 = ((pt[6]).y() - (pt[5]).y()) * 0.5;
    auto pDx3 = ((pt[5]).x() - (pt[4]).x()) * 0.5;
    auto pDx4 = ((pt[7]).x() - (pt[6]).x()) * 0.5;
    auto fTalpha2 = ((pt[6]).x() + (pt[7]).x() - (pt[5]).x() - (pt[4]).x()) * 0.25 / pDy2;
    auto pAlp2 = std::atan(fTalpha2);
    auto fThetaCphi = ((pt[4]).x() + pDy2 * fTalpha2 + pDx3) / pDz;
    auto fThetaSphi = ((pt[4]).y() + pDy2) / pDz;
    double pPhi = 0.;
    double pTheta = 0.;
    if (fThetaSphi == 0. && fThetaCphi == 0.) {
    }
    else {
      pPhi = std::atan(fThetaSphi / fThetaCphi);
      pTheta = std::atan(fThetaCphi / std::cos(pPhi));
    }  // end of Trap parameters calculation

    Trap tower("phiER", pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3, pDx4, pAlp2);
    Volume towerLog("towerLog", tower,
                    description.material(x_tower.attr<std::string>(_U(material))));
    towerLog.setVisAttributes(description, x_tower.visStr());

    RotationX rotX(-thetaofcenter);  // Rotation matrix for this tower
    RotationY rotY(0.);
    RotationZ rotZ(0.);
    Vector3D c = Helper.GetOrigin(0);  // Get origin of tower
    Vector3D c_new(-c.y(), c.z(), c.x() - (innerR + 0.5 * length));
    if (i < NbOfEndcapReduced)
    {  // Place towers up to 35, this "if" is just a protection in case the loop range is changed
#ifndef COUNTTUBES
      std::cout << "----> DREndcapTubes: tower " << i << " being constructed" << std::endl;
#endif
#ifdef COUNTTUBES
      std::cout << "----> DREndcapTubes: tower " << i << " being constructed";
#endif
      Transform3D tower_trnsform(rotX * rotY * rotZ, Position(c_new.x(), c_new.y(), c_new.z()));
      PlacedVolume towerPlaced = phiERLog.placeVolume(towerLog, i, tower_trnsform);
      // ID this volume with tower ID, for the moment I leave air ID to 0 (dummy)
      towerPlaced.addPhysVolID("tower", i).addPhysVolID("air", 0);
    }
    // Or, to debug, place towers one next to each other in assembly volume
    //if(i<35) {
    //    double z = static_cast<int>(i/15)*(length+40*cm);
    //    double x = (i-static_cast<int>(i/15)*15)*100*cm - 5*m;
    //    AssemblyEndcap.placeVolume(towerLog,i,Position(-1.*x,0.,-1.*z));
    //}

    // Capillary placement inside tower (both S and C)
    //
    // Fill a tower along y (vertical direction)
    // per each row calculate x and y of the starting tube (left most)
    const double y_backplane = pt[4].y();  // y-coordinate of bottom left corner of back face
    const double x_backplane = pt[4].x();  // x-coordinate of bottom left corner of back face
    for (std::size_t k = 0; k < 200; k++) {
      double x_start = x_backplane;
      double y_tube = 0.;
      double delta_x = 0.;  // correction to tune x-coordinate
      // if it is the first row fix the y coordinate
      if (k == 0)
        y_tube = y_backplane + tubeRadius
                 + 1.0 * mm;  // add 1 mm otherwise the tube immediately intersect with a plane
      else {
        y_tube = y_backplane + tubeRadius + 1.0 * mm + k * 2. * y_pitch;

        // Adjust starting x given the opening angle of the trapezoidal back face
        double hypotenuse = sqrt(pow((-1 * y_backplane) * 2, 2) + pow(pt[6].x() - pt[4].x(), 2));
        double angle = acos(((-1 * y_backplane) * 2) / hypotenuse);
        delta_x = ((-1. * y_backplane) - (-1. * y_tube)) * tan(angle);
      }
      if (delta_x > tubeDiameter) {
        // Calculate how many fibers should be inserted in x.
        // Caveat: casting to int does not round takes only the
        // integer part (i.e. 6.9 -> 6)
        auto newTubesNo = static_cast<int>(delta_x / tubeDiameter);
        x_start = x_start - newTubesNo * tubeDiameter;
      }

      // Fill a tower row along x (lateral direction)
      for (std::size_t j = 0; j < 1000; j++) {
        auto x_tube =
          x_start + tubeRadius + j * (tubeDiameter);  // locate the center of the tube in x
        Vector3D capillaryPos(x_tube, y_tube, length / 2.);  // locate tube on tower back face
        auto capillaryLength = Helper.GetTubeLength(pt, capillaryPos);  // calculate tube length
        if (std::fabs(capillaryLength - length) < 0.0001 * mm) {
          PlacedVolume capillaryPlaced =
            towerLog.placeVolume(capillary_SLog, CalcCpNo32bits(j,k), Position(x_tube, y_tube, 0.));
          // ID this volume with row ID and column ID
          capillaryPlaced.addPhysVolID("row", k).addPhysVolID("col", j);
#ifdef COUNTTUBES
          tubeNo++;
          totalTubeLength += tower_height;
#endif
        }
        else if (capillaryLength > 5.0 * cm) {
          // Note: the root visualizer does not display tubes with capillaryLength < 4.5 cm.
          // Such tubes are usually the ones closest to the left or right side of a tower.
          // They seem to be placed correctly but not displayed.
          // I am adding a cut over tube length of 5.0 cm: tubes below it will not be placed.
          Tube capillaryShort(0. * mm, tubeRadius, (capillaryLength - TubeLengthOffset) / 2.,
                              2 * M_PI);  // reduced capillary length
                                          // by a fixed offset
          Volume capillaryShortLog(
            "capillaryShortLog", capillaryShort,
            description.material(x_capillary_S.attr<std::string>(_U(material))));
          capillaryShortLog.setVisAttributes(description, x_capillary_S.visStr());

          Tube cladShort_S(0. * mm, cladRadius, (capillaryLength - TubeLengthOffset) / 2.,
                           2 * M_PI);  // cladding S
          Volume cladShort_SLog("cladShort_SLog", cladShort_S,
                                description.material(x_clad_S.attr<std::string>(_U(material))));
          cladShort_SLog.setVisAttributes(description, x_clad_S.visStr());
          PlacedVolume cladShort_SPlaced =
            capillaryShortLog.placeVolume(cladShort_SLog, 1, Position(0., 0., 0.));
          cladShort_SPlaced.addPhysVolID("clad", 1);
          Tube coreShort_S(0. * mm, coreRadius, (capillaryLength - TubeLengthOffset) / 2.,
                           2 * M_PI);  // core S
          Volume coreShort_SLog("DRETScoreSShort_Log", coreShort_S,
                                description.material(x_core_S.attr<std::string>(_U(material))));
          coreShort_SLog.setVisAttributes(description, x_core_S.visStr());
          if (x_core_S_sens) coreShort_SLog.setSensitiveDetector(sens);
          PlacedVolume coreShort_SPlaced = cladShort_SLog.placeVolume(coreShort_SLog, CalcCpNo2bits(1,0));
          coreShort_SPlaced.addPhysVolID("core", 1).addPhysVolID("cherenkov", 0);

          PlacedVolume capillaryShortPlaced = towerLog.placeVolume(
            capillaryShortLog, CalcCpNo32bits(j,k),
            Position(x_tube, y_tube, length / 2. - capillaryLength / 2. + TubeLengthOffset / 2.));
          // ID this volume with row ID and column ID
          capillaryShortPlaced.addPhysVolID("row", k).addPhysVolID("col", j);
#ifdef COUNTTUBES
          tubeNo++;
          totalTubeLength += capillaryLength;  // neglecting the offset
#endif
        }
        else {
          ;
        }

        // Now C tubes are placed.
        // Check there is enough room in y-direction
        bool IsLastRow_C = (-1. * y_backplane) - y_tube < y_pitch + tubeRadius ? true : false;
        // Check there is enough room in x-direction to place a C tube after its S one
        bool IsLastTube_C = -1. * x_backplane + delta_x - x_tube < tubeDiameter ? true : false;

        // After the S tube placement I place the closest C tube above it
        // according to the fixed structure of the tubes placement (gluing)
        //
        // If the S tube below was not placed (too short) do not place the C either
        // && if there is not enough room in x-direction do not place the C tube
        // && if there is not enough room in y-direction do not place the C tube
        if (capillaryLength > 5.0 * cm && !IsLastTube_C && !IsLastRow_C) {
          double x_tube_C = x_tube + tubeRadius;
          double y_tube_C = y_tube + y_pitch;
          Vector3D capillaryPos_C(x_tube_C, y_tube_C, length / 2.);
          auto capillaryLength_C = Helper.GetTubeLength(pt, capillaryPos_C);
          if (std::fabs(capillaryLength_C - length) < 0.0001 * mm) {
            PlacedVolume capillaryPlaced_C = towerLog.placeVolume(capillary_CLog, CalcCpNo32bits(j,k),
                                                                  Position(x_tube_C, y_tube_C, 0.));
            capillaryPlaced_C.addPhysVolID("row", k).addPhysVolID("col", j);
#ifdef COUNTTUBES
            tubeNo++;
            totalTubeLength += tower_height;
#endif
          }
          else if (capillaryLength_C > 5.0 * cm) {
            Tube capillaryShort_C(0. * mm, tubeRadius, (capillaryLength_C - TubeLengthOffset) / 2.,
                                  2 * M_PI);  // reduced capillary length
                                              // by a fixed offset
            Volume capillaryShortLog_C(
              "capillaryShortLog_C", capillaryShort_C,
              description.material(x_capillary_C.attr<std::string>(_U(material))));
            capillaryShortLog_C.setVisAttributes(description, x_capillary_C.visStr());

            Tube cladShort_C(0. * mm, cladRadius, (capillaryLength_C - TubeLengthOffset) / 2.,
                             2 * M_PI);  // cladding C
            Volume cladShort_CLog("DRETScladShort_CLog", cladShort_C,
                                  description.material(x_clad_C.attr<std::string>(_U(material))));
            cladShort_CLog.setVisAttributes(description, x_clad_C.visStr());
            if (x_clad_C_sens) cladShort_CLog.setSensitiveDetector(sens);
            PlacedVolume cladShort_CPlaced =
              capillaryShortLog_C.placeVolume(cladShort_CLog, 1, Position(0., 0., 0.));
            cladShort_CPlaced.addPhysVolID("clad", 1);
            Tube coreShort_C(0. * mm, coreRadius, (capillaryLength_C - TubeLengthOffset) / 2.,
                             2 * M_PI);  // core C
            Volume coreShort_CLog("DRETScoreCShort_Log", coreShort_C,
                                  description.material(x_core_C.attr<std::string>(_U(material))));
            coreShort_CLog.setVisAttributes(description, x_core_C.visStr());
            if (x_core_C_sens) coreShort_CLog.setSensitiveDetector(sens);
            PlacedVolume coreShort_CPlaced = cladShort_CLog.placeVolume(coreShort_CLog, CalcCpNo2bits(1,1));
            coreShort_CPlaced.addPhysVolID("core", 1).addPhysVolID("cherenkov", 1);

            PlacedVolume capillaryShortPlaced_C = towerLog.placeVolume(
              capillaryShortLog_C, CalcCpNo32bits(j,k),
              Position(x_tube_C, y_tube_C,
                       length / 2. - capillaryLength_C / 2. + TubeLengthOffset / 2.));
            capillaryShortPlaced_C.addPhysVolID("row", k).addPhysVolID("col", j);
#ifdef COUNTTUBES
            tubeNo++;
            totalTubeLength += capillaryLength_C;  // neglecting the offset
#endif
          }
          else {
            ;
          }
        }

        // condition for stopping S capillary placement along x
        if (-1. * x_backplane + delta_x - x_tube < tubeDiameter + tubeRadius) break;
      }  // end x loop

      // Condition for stopping S capillary placement along y.
      // y_backplane is equal up and down so I can keep the same for exiting loop
      if ((-1. * y_backplane) - y_tube < (2. * y_pitch + tubeRadius)) break;
    }  // End y loop and tube placement

    // Update parameters
    Helper.Getpt(pt);
    fulltheta = fulltheta - deltatheta_endcap[i];
#ifdef COUNTTUBES
    totalTubeNo += tubeNo;
    std::cout << " with " << tubeNo << " tubes" << std::endl;
#endif
  }  // End of towers creation and placement

#ifdef COUNTTUBES
  std::cout << "--> DREndcapTubes: number of tubes per stave/phi-slice " << totalTubeNo
            << std::endl;
  std::cout << "--> DREndcapTubes: number of tubes for both endcaps " << totalTubeNo * NbOfZRot * 2
            << " and total length " << (NbOfZRot * 2 * totalTubeLength) / km << " km" << std::endl;
#endif

  // Create a (DetElement) corresponding to MyDetector.
  // From DD4hep docs
  // https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManualch2.html "Construct the
  // main detector element of this subdetector.This will be the unique entry point to access any
  // information of the subdetector."
  DetElement sdet(det_name, x_det.id());
  // Then "Place the subdetector envelope into its mother (typically the top level (world) volume)."
  Volume motherVolume = description.pickMotherVolume(sdet);
  // Place the assembly container inside the mother volume
  PlacedVolume AssemblyEndcapPV = motherVolume.placeVolume(AssemblyEndcap);
  sdet.setPlacement(AssemblyEndcapPV);

  std::cout << "--> DREndcapTubes::create_detector() end" << std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DREndcapTubes, create_detector)

//**************************************************************************

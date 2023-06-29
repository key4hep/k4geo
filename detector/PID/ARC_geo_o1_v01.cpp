//----------------------------------
//         ARC detector v01
//----------------------------------

/*!
 *  \brief     Detector constructor of barrel of ARC detector
 *  \details   This code creates full geometry of ARC detector
 *             Evolved from the pfRICH example in DD4hep.
 *  \author    Alvaro Tolosa-Delgado alvaro.tolosa.delgado@cern.ch
 *  \author    Martin Tat            martin.tat@cern.ch
 *  \version   1
 *  \date      2023
 *  \pre       DD4hep compiled with Geant4+Qt
 *  \bug       Walls do not reflect optical photons. Hard-coded values in many places.
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"

using namespace dd4hep;

// #define DUMP_SENSOR_POSITIONS

/// Function to build ARC endcaps
static Ref_t create_arc_endcap_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  double zpos_endcap = detElem.attr<double>(_Unicode(zpos));

  auto gasElem    = detElem.child(_Unicode(radiatorgas));
  auto gasvolMat  = desc.material(gasElem.attr<std::string>(_Unicode(material)));
  auto gasvolVis  = desc.visAttributes(gasElem.attr<std::string>(_Unicode(vis)));

  auto vesselElem = detElem.child(_Unicode(vessel));
  auto vesselSkinMat  = desc.material(vesselElem.attr<std::string>(_Unicode(skinMaterial)));
  auto vesselSkinVis  = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(skin_vis)));

  auto vesselBulkMat  = desc.material(vesselElem.attr<std::string>(_Unicode(bulk_material)));
  auto vesselBulkVis  = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(bulk_vis)));

    double bulk_skin_ratio = vesselElem.attr<double>(_Unicode(bulk_skin_ratio));

    if( 0 > bulk_skin_ratio || 1 < bulk_skin_ratio )
        throw std::runtime_error("ARC: bulk_skin_ratio must be a number between 0 and 1");


  // read Martin file and store parameters by name in the map
//   fill_cell_parameters_m();

  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);

  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //          VESSEL PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double vessel_outer_r = desc.constantAsDouble("ARC_ENDCAP_R_OUTER");
  double vessel_inner_r = desc.constantAsDouble("ARC_ENDCAP_R_INNER");
  double vessel_length = desc.constantAsDouble("ARC_ENDCAP_LENGTH");
  double vessel_wall_thickness = desc.constantAsDouble("ARC_VESSEL_WALL_THICKNESS");
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");
  // // //-------------------------------------------------------------// // //


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //         AEROGEL PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double aerogel_thickness = desc.constantAsDouble("ARC_AEROGEL_THICKNESS");
  auto aerogelMat = desc.material("Aerogel_PFRICH");
  // // //-------------------------------------------------------------// // //

  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //         COOLING PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double cooling_thickness = desc.constantAsDouble("ARC_COOLING_THICKNESS");
  // // //-------------------------------------------------------------// // //


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //           CELL PARAMETERS           // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  /// Cell is an hexagonal prysm
  double hexagon_side_length = 14.815 * cm;
  double hexagon_apothem = hexagon_side_length * cos(30 * deg);
  /// each cell corresponds to one object of the following class
  /// which gathers all the important parameters
  /// this info can be moved later to compact file
  struct mycell_t
  {
    /// Martin number for row
    int row = {0};
    /// Martin number for column
    int col = {0};
    /// Roger ID number
    int RID = {-1};
    /// x position of the cell
    double x = {0.};
    /// y position of the cell
    double y = {0.};
    /// if reflected
    bool isReflected = {false};
  };

// // // // // // // // // // // // // // // // // // // // // // //
//   SCHEME OF UNIQUE CELLS INSIDE A SECTOR
//   CELL NUMBERING CORRESPONDS TO ROGERS
//   MARTIN NUMBERING SPECIFIED BY ROW/COLUMN
//   THIS SECTOR MUST BE MIRRORED, AND THEN
//   BOTH (ORIGINAL AND MIRRORED) REPEATED 6 TIMES
//
//                         _____         _____
//                        /     \       /     \    .
//    7             _____/  21   \_____/  18   \   .
//                 /     \       /     \       /   .
//    7           /  20   \_____/  17   \_____/    .
//                \       /     \       /     \    .
//    6            \_____/  16   \_____/  14   \   .
//                 /     \       /     \       /   .
//    6           /  15   \_____/  13   \_____/    .
//                \       /     \       /     \    .
//    5            \_____/  12   \_____/  10   \   .
//                 /     \       /     \       /   .
//    5           /  11   \_____/   9   \_____/    .
//                \       /     \       /     \    .
//    4            \_____/   8   \_____/   7   \   .
//                       \       /     \       /   .
//    4                   \_____/   6   \_____/    .
//                        /     \       /     \    .
//    3                  /   5   \_____/   4   \   .
//                       \       /     \       /   .
//    3                   \_____/   3   \_____/    .
//                              \       /     \    .
//    2                          \_____/   2   \   .
//                               /     \       /   .
//    2                         /   1   \_____/    .
//                              \       /          .
//    COLUMN ^                   \_____/           .
//    ROW->       4        3        2     1
//
//   Y axis = column
//   X axis = row
// // // // // // // // // // // // // // // // // // // // // // //

  /// vector with cell geometric parameters
  std::vector<mycell_t> mycell_v(21);
  {
    double hx_u = hexagon_apothem;
    double hx_x = hexagon_side_length;
    mycell_v[0] = {1, 2, 2, 0, 4 * hx_u};
    mycell_v[1] = {1, 3, 4, 0, 6 * hx_u};
    mycell_v[2] = {1, 4, 7, 0, 8 * hx_u};
    mycell_v[3] = {1, 5, 10, 0, 10 * hx_u};
    mycell_v[4] = {1, 6, 14, 0, 12 * hx_u};
    mycell_v[5] = {1, 7, 18, 0, 14 * hx_u};
    mycell_v[6] = {2, 2, 1, -1.5 * hx_x, 3 * hx_u};
    mycell_v[7] = {2, 3, 3, -1.5 * hx_x, 5 * hx_u, true};
    mycell_v[8] = {2, 4, 6, -1.5 * hx_x, 7 * hx_u, true};
    mycell_v[9] = {2, 5, 9, -1.5 * hx_x, 9 * hx_u, true};
    mycell_v[10] = {2, 6, 13, -1.5 * hx_x, 11 * hx_u, true};
    mycell_v[11] = {2, 7, 17, -1.5 * hx_x, 13 * hx_u, true};
    mycell_v[12] = {3, 3, 5, -3.0 * hx_x, 6 * hx_u};
    mycell_v[13] = {3, 4, 8, -3.0 * hx_x, 8 * hx_u, true};
    mycell_v[14] = {3, 5, 12, -3.0 * hx_x, 10 * hx_u, true};
    mycell_v[15] = {3, 6, 16, -3.0 * hx_x, 12 * hx_u, true};
    mycell_v[16] = {3, 7, 21, -3.0 * hx_x, 14 * hx_u, true};
    mycell_v[17] = {4, 5, 11, -4.5 * hx_x, 9 * hx_u};
    mycell_v[18] = {4, 6, 15, -4.5 * hx_x, 11 * hx_u, true};
    mycell_v[19] = {4, 7, 20, -4.5 * hx_x, 13 * hx_u, true};
    mycell_v[20] = {5, 6, 19, -6.0 * hx_x, 12 * hx_u};
  }

  /// Distance in phi angle between complete sectors
  double phistep = 60 * deg;
  /// number of repetition of sectors
  int phinmax = 6; // 6;


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //          MIRROR PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double mirror_z_origin_Martin = vessel_length / 2. - vessel_wall_thickness - 37 * cm;
  auto mirrorElem = detElem.child(_Unicode(mirror));
  double mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
  auto mirrorSurf = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));
  auto mirrorMat = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  // // //-------------------------------------------------------------// // //


    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // //          LIGHT SENSOR PARAMETERS          // // // // // //
    // // // // // // // // // // // // // // // // // // // // // // // // // //
    //default values
    double sensor_sidex     = 8 * cm;
    double sensor_sidey     = 8 * cm;
    double sensor_thickness = 0.2 * cm;
    // empirical distance to keep the sensor inside the cell volume
    double sensor_z_origin_Martin = -vessel_length / 2. + vessel_wall_thickness + 0.5 * cooling_thickness;
    auto sensorMat = desc.material("SiliconOptical");
    auto sensorVis = desc.visAttributes("no_vis");

    // Read from xml the parameters for the sensor module
    {
        auto sensorElem  = detElem.child(_Unicode(sensors));
        sensor_sidex     = sensorElem.attr<double>(_Unicode(sensor_side_X));
        sensor_sidey     = sensorElem.attr<double>(_Unicode(sensor_side_Y));
        sensor_thickness = sensorElem.attr<double>(_Unicode(thickness));
        sensorMat        = desc.material(sensorElem.attr<std::string>(_Unicode(material)));
        sensorVis        = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
    }
    // // //-------------------------------------------------------------// // //


    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //
    // // //+++++++++++  BUILD VESSEL, CELL AND SENSOR VOLUMES ++++++++++// // //
    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

    // Build cylinder for gas, and the vessel for the gas
    Tube gasenvelopeS(  vessel_inner_r + vessel_wall_thickness,
                        vessel_outer_r - vessel_wall_thickness,
                        vessel_length/2.);
    Volume endcap_cells_gas_envelope (detName+"_gasEnvelope", gasenvelopeS, gasvolMat );
    endcap_cells_gas_envelope.setVisAttributes( desc.visAttributes("arc_envelope_vis") );

    Tube vesselEnvelopeSolid(  vessel_inner_r,
                               vessel_outer_r,
                               vessel_length/2. + vessel_wall_thickness);
    Volume endcap_cells_vessel_envelope (detName+"_vesselEnvelope", vesselEnvelopeSolid, vesselSkinMat );
    endcap_cells_vessel_envelope.setVisAttributes( vesselSkinVis );

    // if 0==bulk_skin_ratio do not create bulk at all
    if(0<bulk_skin_ratio)
    {
      // build bulk for inner wall
      double vessel_bulk_inner_r_ini = vessel_inner_r + (1 - bulk_skin_ratio)*0.5*vessel_wall_thickness;
      double vessel_bulk_inner_r_fin = vessel_inner_r + (1 + bulk_skin_ratio)*0.5*vessel_wall_thickness;

      Tube vesselInnerBulkSolid( vessel_bulk_inner_r_ini,
                            vessel_bulk_inner_r_fin,
                            vessel_length/2. + vessel_wall_thickness - (1-bulk_skin_ratio)*0.5*vessel_wall_thickness);
      Volume vessel_innerbulk_vol (detName+"_vesselInnerBulk", vesselInnerBulkSolid, vesselBulkMat );
      vessel_innerbulk_vol.setVisAttributes( vesselBulkVis );
      endcap_cells_vessel_envelope.placeVolume(vessel_innerbulk_vol);

      // build bulk for outer wall
      double vessel_bulk_outer_r_ini = vessel_outer_r - (1 + bulk_skin_ratio)*0.5*vessel_wall_thickness;
      double vessel_bulk_outer_r_fin = vessel_outer_r - (1 - bulk_skin_ratio)*0.5*vessel_wall_thickness;

      Tube vesselOuterBulkSolid( vessel_bulk_outer_r_ini,
                                vessel_bulk_outer_r_fin,
                                vessel_length/2. + vessel_wall_thickness -  (1-bulk_skin_ratio)*0.5*vessel_wall_thickness);
      Volume vessel_outerbulk_vol (detName+"_vesselOuterBulk", vesselOuterBulkSolid, vesselBulkMat );
      vessel_outerbulk_vol.setVisAttributes( vesselBulkVis );
      endcap_cells_vessel_envelope.placeVolume(vessel_outerbulk_vol);

      Tube vesselBaseBulkSolid(  vessel_bulk_inner_r_fin,
                                vessel_bulk_outer_r_ini,
                                bulk_skin_ratio*0.5*vessel_wall_thickness);
      Volume vessel_base_bulk_vol (detName+"_vesselBaseBulk", vesselBaseBulkSolid, vesselBulkMat );
      vessel_base_bulk_vol.setVisAttributes( vesselBulkVis );
      auto posZPositive = Position(0, 0, vessel_length/2. + 0.5*vessel_wall_thickness);
      endcap_cells_vessel_envelope.placeVolume(vessel_base_bulk_vol,posZPositive);

      auto posZNegative = Position(0, 0, -vessel_length/2. - 0.5*vessel_wall_thickness);
      endcap_cells_vessel_envelope.placeVolume(vessel_base_bulk_vol,posZNegative);
    }


    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

  // Use regular polyhedra for endcaps cells
  PolyhedraRegular cellS(6, 0, 0., hexagon_apothem, vessel_length);

  // Build sensor shape
  Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
  Volume sensorVol(detName + "_sensor", sensorSol, sensorMat);
  sensorVol.setSensitiveDetector(sens);
  sensorVol.setVisAttributes( sensorVis );

  // Build cooling plate
  double cooling_z_offset =   sensor_thickness  + cooling_thickness + 0.5*mm;
  Tube coolingSol_tube(0, 1.5*hexagon_side_length, cooling_thickness);

  // Build aerogel plate
  double aerogel_z_offset =   sensor_thickness  + aerogel_thickness + 0.5*mm;
  Tube aerogelSol_tube(0, 1.5*hexagon_side_length, aerogel_thickness);

  // Build cells of a sector
//   mycell_v = {mycell_v[16], mycell_v[19]};
//   phinmax = 1;
  int cellCounter = 0;
  int physicalVolumeCounter = 0;
  auto createPhysVolID = [&](){return physicalVolumeCounter++;};

#ifdef DUMP_SENSOR_POSITIONS
  std::ofstream ofile_sensor_pos("ofile_sensor_pos_endcap.txt");
#endif
  for (auto &ncell : mycell_v)
  {
    // sanity check, skip non initialized cells
    if (-1 == ncell.RID)
      continue;

    // The following line skips even number cells
    // if ( 1 != ncell.row )
    //   continue;

    /// repeat the sector 6 times
    for (int phin = 0; phin < phinmax; phin++, cellCounter++)
    {

      /// function to create names in a systematic way
      /// final name = detName + part + cell parameters
      auto create_part_name_ff = [ncell,detName,phin](std::string  partName){
          std::string fullName = detName + "_" + partName;
          fullName += std::to_string(ncell.RID);
          fullName += "_phi" +  std::to_string(phin);
          dd4hep::printout(dd4hep::DEBUG,"ARCENDCAP_T", "+++ New name:%s",fullName.c_str());
          return fullName;
        };

      /// cell volume, hex prism
      /// the elements must be placed inside
      std::string cellName = create_part_name_ff("cell");
      Volume cellV(cellName, cellS, gasvolMat);
      cellV.setVisAttributes( gasvolVis );
      /// Detector element that will contain cellVol later
      /// there are 3 elements with ID:
      /// the cell, ID= 6 * cellCounter
      /// its mirror. ID = 6 * cellCounter +1
      /// and its sensor, ID = 6 * cellCounter +2
      DetElement cellDE(det, cellName+"DE", 6 * cellCounter + 0);



      // // initialize parameters for creating the mirror
      double center_of_sphere_x(-999.);
      double center_of_sphere_z(-999.);
      double radius_of_sphere(-999.);

      double center_of_sensor_x(-999.);
      double angle_of_sensor(-999.);
      double zoffset_of_sensor(0);


      // retrieve cell parameters
      // if parameter not present, exception is thrown and not catched
      {
        // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
        std::string name_col_s = std::to_string(ncell.col);
        std::string name_row_s = std::to_string(ncell.row);
        std::string MartinCellName = "EndCapRadiator_c" + name_col_s + "_r" + name_row_s;

        radius_of_sphere = desc.constantAsDouble(MartinCellName + "_Curvature");

        center_of_sphere_x = desc.constantAsDouble(MartinCellName + "_XPosition");

        double zposition = desc.constantAsDouble(MartinCellName + "_ZPosition");

        center_of_sphere_z = mirror_z_origin_Martin + zposition;

        center_of_sensor_x = desc.constantAsDouble(MartinCellName + "_DetPosition");

        angle_of_sensor = desc.constantAsDouble(MartinCellName + "_DetTilt");

        std::string ZOffsetSensorParName = MartinCellName + "_DetZOffset";


        if( desc.constants().count(MartinCellName + "_DetPositionZ") )
            zoffset_of_sensor = desc.constantAsDouble(MartinCellName + "_DetPositionZ");
        else
          dd4hep::printout(dd4hep::WARNING,"ARCENDCAP_T", "+++ Constant %s is missing in xml file, default is 0",ZOffsetSensorParName.c_str());

        if (radius_of_sphere <= mirrorThickness)
          throw std::runtime_error(Form("Ilegal parameters cell %d: %g <= %g", ncell.RID, radius_of_sphere, mirrorThickness));

      }
      double sensor_z_pos = zoffset_of_sensor + sensor_z_origin_Martin;

      // create the semi-sphere that will result in the mirror
      Sphere mirrorShapeFull(radius_of_sphere - mirrorThickness,
                             radius_of_sphere,
                             0.,
                             3.14 / 3);
      /// alpha: angle of position vector of first sector n-cell with respect to x-axis
      double alpha = atan(ncell.y / ncell.x) * rad;
      if (0 > alpha)
        alpha += 180 * deg;
      double dx = center_of_sphere_x * cos(alpha);
      double dy = center_of_sphere_x * sin(alpha);

      /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
      Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(dx, dy, center_of_sphere_z));

      /// Define the actual mirror as intersection of the hex cell volume and the hollow sphere just defined
      Solid mirrorSol = IntersectionSolid(cellS, mirrorShapeFull, mirrorTr);
      std::string mirrorVolName = create_part_name_ff("mirror");
      Volume mirrorVol(mirrorVolName, mirrorSol, mirrorMat);
      mirrorVol.setVisAttributes(desc.visAttributes(Form("arc_mirror_vis%d", ncell.RID)));
      PlacedVolume mirrorPV = cellV.placeVolume(mirrorVol);

      DetElement mirrorDE(cellDE, mirrorVolName + "DE", 6 * cellCounter+1 );
      mirrorDE.setPlacement(mirrorPV);
      SkinSurface mirrorSkin(desc, mirrorDE, Form("mirror_optical_surface%d", cellCounter), mirrorSurf, mirrorVol); // FIXME: 3rd arg needs `imod`?
      mirrorSkin.isValid();


      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  COOLING PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      auto coolingTrCell = RotationZYX(0, 0, angle_of_sensor ) *
                           Translation3D(0, center_of_sensor_x, sensor_z_pos-cooling_z_offset);

      Solid coolingSol = IntersectionSolid(cellS, coolingSol_tube, coolingTrCell);
      std::string coolingName = create_part_name_ff("cooling");
      /// TODO: change material
      Volume coolingVol( coolingName , coolingSol, mirrorMat );
      coolingVol.setVisAttributes( desc.visAttributes("arc_cooling_vis") );
      cellV.placeVolume(coolingVol);

      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  AEROGEL PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      auto aerogelTrCell = RotationZYX(0, 0, angle_of_sensor ) *
                           Translation3D(0, center_of_sensor_x, sensor_z_pos+aerogel_z_offset);

      Solid aerogelSol = IntersectionSolid(cellS, aerogelSol_tube, aerogelTrCell);
      std::string aerogelName = create_part_name_ff("aerogel");
      Volume aerogelVol( aerogelName , aerogelSol, aerogelMat );
      aerogelVol.setVisAttributes( desc.visAttributes("arc_aerogel_vis") );
      cellV.placeVolume(aerogelVol);

      auto sensorTr = RotationZYX(alpha - 90 * deg, 0 , angle_of_sensor )*
                           Translation3D(0, center_of_sensor_x, sensor_z_pos );


      PlacedVolume sensorPV = cellV.placeVolume(sensorVol, sensorTr);
//       sensorPV.addPhysVolID("cellnumber", 6 * cellCounter+2);
#ifdef DUMP_SENSOR_POSITIONS
      ofile_sensor_pos  << 6 * cellCounter+2 << '\t'
                         << ncell.RID << '\t'
                         << ncell.isReflected << '\t'
                         << ncell.x << '\t'
                         << ncell.y << '\t'
                         << phin << '\n';
#endif
      DetElement sensorDE(cellDE, create_part_name_ff("sensor") + "DE", 6 * cellCounter+2 );
      sensorDE.setType("tracker");
      sensorDE.setPlacement(sensorPV);

      PlacedVolume cellPV = endcap_cells_gas_envelope.placeVolume(cellV, RotationZ(phistep * phin) * Translation3D(ncell.x, ncell.y, 0));
      cellPV.addPhysVolID("cellnumber", createPhysVolID() );//6*cellCounter + 0);
      cellDE.setPlacement( cellPV );

      if ( ncell.isReflected)
      {
        std::string cellRefName = create_part_name_ff("cell_ref");
        Volume cellV_reflected(cellRefName, cellS, gasvolMat);
        cellV_reflected.setVisAttributes( gasvolVis );
        DetElement cell_reflected_DE(det, cellRefName+"DE", 6 * cellCounter + 3);
        Transform3D mirrorTr_reflected(RotationZYX(0, 0, 0), Translation3D(-dx, dy, center_of_sphere_z));

        /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
        Solid mirrorSol_reflected = IntersectionSolid(cellS, mirrorShapeFull, mirrorTr_reflected);
        Volume mirrorVol_reflected(mirrorVolName + "_ref1", mirrorSol_reflected, mirrorMat);
        mirrorVol_reflected.setVisAttributes(desc.visAttributes(Form("arc_mirror_vis%d", ncell.RID)));
        PlacedVolume mirror_ref_PV = cellV_reflected.placeVolume(mirrorVol_reflected);

        DetElement mirror_ref_DE(cell_reflected_DE, mirrorVolName + "_ref1" + "DE", 6 * cellCounter+4 );
        mirror_ref_DE.setPlacement(mirror_ref_PV);
        SkinSurface mirror_ref_Skin(desc, mirror_ref_DE, Form("mirror_ref_optical_surface%d", cellCounter), mirrorSurf, mirrorVol_reflected); // FIXME: 3rd arg needs `imod`?
        mirror_ref_Skin.isValid();

        auto sensorTr_reflected = RotationZYX(-alpha + 90 * deg, 0 /*90*deg-angle_of_sensor*/, angle_of_sensor)*
                                       Translation3D(0, center_of_sensor_x, sensor_z_origin_Martin);
        PlacedVolume sensor_ref_PV = cellV_reflected.placeVolume(sensorVol, sensorTr_reflected);
//         sensor_ref_PV.addPhysVolID("cellnumber", 6 * cellCounter+5);
        DetElement sensor_ref_DE(cell_reflected_DE, create_part_name_ff("sensor") + "_ref_DE", 6 * cellCounter+5 );
        sensor_ref_DE.setPlacement(sensor_ref_PV);

#ifdef DUMP_SENSOR_POSITIONS
        ofile_sensor_pos  << 6 * cellCounter+5 << '\t'
                    << ncell.RID << '\t'
                    << ncell.isReflected << '\t'
                    << ncell.x << '\t'
                    << ncell.y << '\t'
                    << phin << '\n';
#endif

        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  COOLING PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  AEROGEL PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
          cellV_reflected.placeVolume(coolingVol);
          cellV_reflected.placeVolume(aerogelVol);


        PlacedVolume cell_ref_PV = endcap_cells_gas_envelope.placeVolume(cellV_reflected, RotationZ(phistep * phin) * Translation3D(-ncell.x, ncell.y, 0));
        cell_ref_PV.addPhysVolID("cellnumber", createPhysVolID() );//6*cellCounter + 3);
        cell_reflected_DE.setPlacement( cell_ref_PV );
      }
    } //-- end loop for sector
  }   //-- end loop for endcap

  endcap_cells_vessel_envelope.placeVolume(endcap_cells_gas_envelope);



  Assembly endcaps_assemblyV("endcaps_assemblyV");

  Transform3D endcapZPos_Tr(RotationZYX(0,0,0), Translation3D(0, 0, zpos_endcap));
  PlacedVolume endcapZPos_PV = endcaps_assemblyV.placeVolume(endcap_cells_vessel_envelope, endcapZPos_Tr);
  endcapZPos_PV.addPhysVolID("barrel", 1);

  DetElement endcapZPos_DE(det, "endcapZPos_DE", 0 );
  endcapZPos_DE.setPlacement(endcapZPos_PV);


  Transform3D envelope_zreflected_Tr(RotationZYX( 0 ,0,180*deg), Translation3D(0, 0, -zpos_endcap));
  PlacedVolume endcapZNeg_PV = endcaps_assemblyV.placeVolume(endcap_cells_vessel_envelope, envelope_zreflected_Tr);
  endcapZNeg_PV.addPhysVolID("barrel", 2);

  DetElement endcapZNeg_DE(det, "endcapZNeg_DE", 2 );
  endcapZNeg_DE.setPlacement(endcapZNeg_PV);


  PlacedVolume endcaps_PV = motherVol.placeVolume(endcaps_assemblyV);
  endcaps_PV.addPhysVolID("system", detID);
  det.setPlacement(endcaps_PV);

  return det;
}
DECLARE_DETELEMENT(ARCENDCAP_o1_v01_T, create_arc_endcap_cell)


/// Function to build ARC barrel
/**
 * The geometry tree is the following:
 * vessel (CarbFib)-> vessel bulk (foam)
 *                 -> gas envelope (gas)-> gas cell 1 (gas) -> elements (mirror, sensor, aeogel, cooling)
 *                 -> gas envelope (gas)-> gas cell 2 (gas) -> elements (mirror, sensor, aeogel, cooling)
 *                 -> gas envelope (gas)-> gas cell ...
 */
static Ref_t create_arc_barrel_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    OpticalSurfaceManager surfMgr = desc.surfaceManager();
    DetElement det(detName, detID);
    sens.setType("tracker");

    auto gasElem    = detElem.child(_Unicode(radiatorgas));
    auto gasvolMat  = desc.material(gasElem.attr<std::string>(_Unicode(material)));
    auto gasvolVis  = desc.visAttributes(gasElem.attr<std::string>(_Unicode(vis)));

    auto vesselElem = detElem.child(_Unicode(vessel));
    auto vesselSkinMat  = desc.material(vesselElem.attr<std::string>(_Unicode(skinMaterial)));
    auto vesselSkinVis  = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(skin_vis)));

    auto vesselBulkMat  = desc.material(vesselElem.attr<std::string>(_Unicode(bulk_material)));
    auto vesselBulkVis  = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(bulk_vis)));

    double bulk_skin_ratio = vesselElem.attr<double>(_Unicode(bulk_skin_ratio));

    if( 0 > bulk_skin_ratio || 1 < bulk_skin_ratio )
        throw std::runtime_error("ARC: bulk_skin_ratio must be a number between 0 and 1");

    // mother volume corresponds to the world
    Volume motherVol = desc.pickMotherVolume(det);


    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // //          VESSEL PARAMETERS          // // // // // //
    // // // // // // // // // // // // // // // // // // // // // // // // // //
    double vessel_outer_r = desc.constantAsDouble("ARC_BARREL_R_OUTER");
    double vessel_inner_r = desc.constantAsDouble("ARC_BARREL_R_INNER");
    double vessel_length = desc.constantAsDouble("ARC_BARREL_LENGTH");
    double vessel_wall_thickness = desc.constantAsDouble("ARC_VESSEL_WALL_THICKNESS");
    if (vessel_outer_r <= vessel_inner_r)
        throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");
    // // //-------------------------------------------------------------// // //


    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // //         AEROGEL PARAMETERS          // // // // // //
    // // // // // // // // // // // // // // // // // // // // // // // // // //
    double aerogel_radial_thickness = desc.constantAsDouble("ARC_AEROGEL_THICKNESS");
    auto aerogelMat = desc.material("Aerogel_PFRICH");

    // // //-------------------------------------------------------------// // //

    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // //         COOLING PARAMETERS          // // // // // //
    // // // // // // // // // // // // // // // // // // // // // // // // // //
    double cooling_radial_thickness = desc.constantAsDouble("ARC_COOLING_THICKNESS");
    // // //-------------------------------------------------------------// // //


    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // //           CELL PARAMETERS           // // // // // //
    // // // // // // // // // // // // // // // // // // // // // // // // // //
    /// Distance in x-direction
    double hexagon_side_length = 14.815 * cm;
    double hex_apothem_length = hexagon_side_length*cos( M_PI / 6. ); //
    double zstep = 2 * hex_apothem_length;
    /// Distance in phi angle between cells
    /// since cells are regular hexagons, this distance matches
    /// the angle that one cell covers
    double phistep = 13.333 * deg;
    /// number of repetition of unique cells around the barrel
    int phinmax = 27; // 27;
    // // //-------------------------------------------------------------// // //


    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // // // //          MIRROR PARAMETERS          // // // // // //
    // // // // // // // // // // // // // // // // // // // // // // // // // //
    auto mirrorElem = detElem.child(_Unicode(mirror));
    auto mirrorSurf = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));
    auto mirrorMat = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
    double mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
    // if this z shrink is not applied, the upper tip of the mirrors are cut
    // TODO: crosscheck with Martin distances between mirrors and sensors
    double mirror_z_safe_shrink = 6*mm;
    double mirror_z_origin_Martin = vessel_outer_r - vessel_wall_thickness - 37 * cm  - mirrorThickness - mirror_z_safe_shrink;
    // // //-------------------------------------------------------------// // //



    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // //          LIGHT SENSOR PARAMETERS          // // // // // //
    // // // // // // // // // // // // // // // // // // // // // // // // // //
    //default values
    double sensor_sidex     = 8 * cm;
    double sensor_sidey     = 8 * cm;
    double sensor_thickness = 0.2 * cm;
    double sensor_z_origin_Martin = vessel_inner_r + vessel_wall_thickness + 5*mm;
    auto sensorMat = desc.material("SiliconOptical");
    auto sensorVis = desc.visAttributes("arc_no_vis");
    // auto sensorSurf = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));

    // Read from xml the parameters for the sensor module
    {
        auto sensorElem  = detElem.child(_Unicode(sensors));
        sensor_sidex     = sensorElem.attr<double>(_Unicode(sensor_side_Phi));
        sensor_sidey     = sensorElem.attr<double>(_Unicode(sensor_side_Z));
        sensor_thickness = sensorElem.attr<double>(_Unicode(thickness));
        sensorMat        = desc.material(sensorElem.attr<std::string>(_Unicode(material)));
        sensorVis        = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
    }
    // // //-------------------------------------------------------------// // //



    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //
    // // //+++++++++++  BUILD VESSEL, CELL AND SENSOR VOLUMES ++++++++++// // //
    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

    // Build cylinder for gas, and the vessel for the gas
    Tube gasenvelopeS(  vessel_inner_r + vessel_wall_thickness,
                        vessel_outer_r - vessel_wall_thickness,
                        vessel_length/2.);
    Volume barrel_cells_gas_envelope (detName+"_gasEnvelope", gasenvelopeS, gasvolMat );
    barrel_cells_gas_envelope.setVisAttributes( desc.visAttributes("arc_envelope_vis") );


    Tube vesselEnvelopeSolid(  vessel_inner_r,
                               vessel_outer_r,
                               vessel_length/2. + vessel_wall_thickness);
    Volume barrel_cells_vessel_envelope (detName+"_vesselSkin", vesselEnvelopeSolid, vesselSkinMat );
    barrel_cells_vessel_envelope.setVisAttributes( vesselSkinVis );

    // if 0==bulk_skin_ratio do not create bulk at all
    if(0<bulk_skin_ratio)
    {
        // build bulk for inner wall
        double vessel_bulk_inner_r_ini = vessel_inner_r + (1 - bulk_skin_ratio)*0.5*vessel_wall_thickness;
        double vessel_bulk_inner_r_fin = vessel_inner_r + (1 + bulk_skin_ratio)*0.5*vessel_wall_thickness;

        Tube vesselInnerBulkSolid( vessel_bulk_inner_r_ini,
                            vessel_bulk_inner_r_fin,
                            vessel_length/2. + vessel_wall_thickness - (1-bulk_skin_ratio)*0.5*vessel_wall_thickness);
        Volume vessel_innerbulk_vol (detName+"_vesselInnerBulk", vesselInnerBulkSolid, vesselBulkMat );
        vessel_innerbulk_vol.setVisAttributes( vesselBulkVis );
        barrel_cells_vessel_envelope.placeVolume(vessel_innerbulk_vol);

        // build bulk for outer wall
        double vessel_bulk_outer_r_ini = vessel_outer_r - (1 + bulk_skin_ratio)*0.5*vessel_wall_thickness;
        double vessel_bulk_outer_r_fin = vessel_outer_r - (1 - bulk_skin_ratio)*0.5*vessel_wall_thickness;

        Tube vesselOuterBulkSolid( vessel_bulk_outer_r_ini,
                                vessel_bulk_outer_r_fin,
                                vessel_length/2. + vessel_wall_thickness - (1-bulk_skin_ratio)*0.5*vessel_wall_thickness);
        Volume vessel_outerbulk_vol (detName+"_vesselOuterBulk", vesselOuterBulkSolid, vesselBulkMat );
        vessel_outerbulk_vol.setVisAttributes( vesselBulkVis );
        barrel_cells_vessel_envelope.placeVolume(vessel_outerbulk_vol);

        Tube vesselBaseBulkSolid(  vessel_bulk_inner_r_fin,
                                vessel_bulk_outer_r_ini,
                                bulk_skin_ratio*0.5*vessel_wall_thickness);
        Volume vessel_base_bulk_vol (detName+"_vesselBaseBulk", vesselBaseBulkSolid, vesselBulkMat );
        vessel_base_bulk_vol.setVisAttributes( vesselBulkVis );
        auto posZPositive = Position(0, 0, vessel_length/2. + 0.5*vessel_wall_thickness);
        barrel_cells_vessel_envelope.placeVolume(vessel_base_bulk_vol,posZPositive);

        auto posZNegative = Position(0, 0, -vessel_length/2. - 0.5*vessel_wall_thickness);
        barrel_cells_vessel_envelope.placeVolume(vessel_base_bulk_vol,posZNegative);
    }
    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

    // Define the cell shape and volume
    // while developping, I used Polyhedra
    // TODO: use regular Polyhedra instead?
    double angle_hex = 2*asin( hex_apothem_length / vessel_outer_r );
    std::vector<double> zplanes = { vessel_inner_r + vessel_wall_thickness,
                                    (vessel_outer_r - vessel_wall_thickness)*cos(angle_hex)
                                  };
    std::vector<double> rs = { hex_apothem_length -1*mm, hex_apothem_length-1*mm };
    /// Hexagonal truncated pyramid
    Polyhedra cell_shape("cellShape", 6, 30 * deg, 360 * deg, zplanes, rs);
    /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
    Transform3D pyramidTr(RotationZYX(0, -90. * deg, 0. * deg), Translation3D(0, 0, 0));

    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //


    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // Build the mirror for ncell=1..18
    // negative values correspond the cells that placed for z<0
    std::vector<int> ncell_vector = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17, /*, -18,*/
                                     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 , 17 /*, 18 */
                                    };

    // dummy counter to identify the cell number
    // TODO: improve ID bitfield
    int cellCounter(0);

    // WARNING for developping purposes
//     ncell_vector = {16,17};
//     phinmax = 1;

    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // loop to build each cell, repeated 27 times around phi    // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    for (auto ncell_original : ncell_vector)
    {
        // Activate flag to reflect parameters later
        int ncell = ncell_original;
        bool reflect_parameters = false;
        if (0 > ncell)
        {
            ncell *= -1;
            reflect_parameters = true;
        }
        for (int phin = 0; phin < phinmax; ++phin)
        {

            // WARNING for developping purposes
            // The following line skips even number cells
            // if (!(ncell % 2))
            //   continue;

            // The following line skips odd number cells
            // if ((ncell % 2))
            // continue;
            // there is no cell number 0, and cell number 1 do not need to be reflected
            if (0 == ncell || -1 == ncell)
                continue;


            auto create_part_name_ff = [ncell,detName,phin,reflect_parameters](std::string  partName) {
                std::string fullName = detName + "_" + partName;
                fullName += std::to_string(ncell);
                fullName += "_ref" + std::to_string(reflect_parameters);
                fullName += "_phi" +  std::to_string(phin);
                dd4hep::printout(dd4hep::DEBUG,"ARCBARREL_T", "+++ New name:%s",fullName.c_str());
                return fullName;
            };






            std::string cellName = create_part_name_ff( "cell");

            /// Volume that contains gas and other stuff
            Volume cellVol(cellName, cell_shape, gasvolMat);
            cellVol.setVisAttributes( gasvolVis );
            /// Detector element that will contain cellVol later
            /// there are 3 elements with ID:
            /// the cell, ID= 3 * cellCounter
            /// its mirror. ID = 3 * cellCounter +1
            /// and its sensor, ID = 3 * cellCounter +2
            DetElement cellDE(det, cellName+"DE", 3 * cellCounter);



            // initialize parameters for creating the mirror
            double center_of_sphere_x(-999.);
            double center_of_sphere_z(-999.);
            double radius_of_sphere(-999.);

            double center_of_sensor_x(-999.);
            double center_of_sensor_z_offset(0);
            double angle_of_sensor(-999.);

            // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
            int name_col = ncell / 2;
            int name_row = ncell % 2 ? 1 : 2;
            // retrieve cell parameters
            // if parameter not present, exception is thrown and not catched
            {
                std::string name_col_s = std::to_string(name_col);
                std::string name_row_s = std::to_string(name_row);
                std::string MartinCellName = "Radiator_c" + name_col_s + "_r" + name_row_s;

                radius_of_sphere = desc.constantAsDouble(MartinCellName + "_Curvature");

                center_of_sphere_x = desc.constantAsDouble(MartinCellName + "_XPosition");

                double zposition = desc.constantAsDouble(MartinCellName + "_ZPosition");

                center_of_sensor_x = desc.constantAsDouble(MartinCellName + "_DetPosition");

                if( desc.constants().count(MartinCellName + "_DetPositionZ") )
                    center_of_sensor_z_offset = desc.constantAsDouble(MartinCellName + "_DetPositionZ");

                angle_of_sensor = desc.constantAsDouble(MartinCellName + "_DetTilt");

                center_of_sphere_z = mirror_z_origin_Martin + zposition;

                // check if parameters are ok
                if (radius_of_sphere <= mirrorThickness)
                    throw std::runtime_error("Ilegal parameters: radius_of_sphere <= mirrorThickness");
            }

            // reflect parameters for cells with z<0,
            // ie, ncell < 0
            if (reflect_parameters)
            {
                center_of_sphere_x *= -1.0;
                center_of_sensor_x *= -1.0;
                angle_of_sensor *= -1.0;
            }

            // create the semi-sphere that will result in the mirror
            Sphere mirrorShapeFull(radius_of_sphere - mirrorThickness,
                                   radius_of_sphere,
                                   0.,
                                   3.14 / 2);
            /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
            Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z - mirror_z_safe_shrink));

            // TODO: cell 18 corresponds to half a pyramid, currently is full pyramid
            Solid mirrorSol = IntersectionSolid(cell_shape, mirrorShapeFull, mirrorTr);

            std::string mirrorName = create_part_name_ff("mirror"); // detName + "_mirror" + std::to_string(ncell) + "z" + std::to_string(reflect_parameters)

            Volume mirrorVol( mirrorName, mirrorSol, mirrorMat );
            mirrorVol.setVisAttributes(desc.visAttributes(Form("arc_mirror_vis%d", ncell)));
            PlacedVolume mirrorPV = cellVol.placeVolume(mirrorVol);

            DetElement mirrorDE(cellDE, mirrorName + "DE", 3 * cellCounter+1 );
            mirrorDE.setPlacement(mirrorPV);
            SkinSurface mirrorSkin(desc, mirrorDE, Form("mirror_optical_surface%d", cellCounter), mirrorSurf, mirrorVol); // FIXME: 3rd arg needs `imod`?
            mirrorSkin.isValid();

            // Place detector in cell
            // Build sensor shape
            Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
            std::string sensorName = create_part_name_ff("sensor");
            Volume sensorVol( sensorName, sensorSol, sensorMat );
            sensorVol.setVisAttributes( sensorVis );
            sensorVol.setSensitiveDetector(sens);

            double center_of_sensor_z = center_of_sensor_z_offset + sensor_z_origin_Martin;

            Transform3D sensorTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-center_of_sensor_z, 0, center_of_sensor_x));
            PlacedVolume sensorPV = cellVol.placeVolume(sensorVol, RotationZYX(0, 90. * deg, 0. * deg)*sensorTr);
            sensorPV.addPhysVolID("cellnumber", 3 * cellCounter+2);
            DetElement sensorDE(cellDE, sensorName + "DE", 3 * cellCounter+2 );
            sensorDE.setType("tracker");
            sensorDE.setPlacement(sensorPV);

            // this is an empirical parameter in order to pass the overlaps
            double safe_distance_from_sensor = 262*um;

            // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
            // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  COOLING PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
            // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
            {
                double cooling_z_offset =   sensor_thickness  + cooling_radial_thickness/2 + safe_distance_from_sensor;
                Tube coolingSol_tube(0, 1.5*hexagon_side_length, cooling_radial_thickness);
                Transform3D coolingTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-center_of_sensor_z+cooling_z_offset, 0, center_of_sensor_x));
                auto coolingTrCell = RotationZYX(0, 90. * deg, 0. * deg)*coolingTr;
                Solid coolingSol = IntersectionSolid(cell_shape, coolingSol_tube, coolingTrCell);
                std::string coolingName = create_part_name_ff("cooling");
                /// TODO: change material
                Volume coolingVol( coolingName, coolingSol, mirrorMat );
                coolingVol.setVisAttributes( desc.visAttributes("arc_cooling_vis") );
                cellVol.placeVolume(coolingVol );
            }
            // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //

            // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
            // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  AEROGEL PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
            // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
            double aerogel_z_offset =   sensor_thickness - aerogel_radial_thickness/2 - safe_distance_from_sensor;
            Tube aerogelSol_tube(0, 1.5*hexagon_side_length, cooling_radial_thickness);
            Transform3D aerogelTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-center_of_sensor_z + aerogel_z_offset, 0, center_of_sensor_x));
            auto aerogelTrCell = RotationZYX(0, 90. * deg, 0. * deg)*aerogelTr;
            Solid aerogelSol = IntersectionSolid(cell_shape, aerogelSol_tube, aerogelTrCell);
            std::string aerogelName = create_part_name_ff("aerogel");
            /// TODO: change material
            Volume aerogelVol( aerogelName, aerogelSol, aerogelMat );
            aerogelVol.setVisAttributes( desc.visAttributes("arc_aerogel_vis") );
            cellVol.placeVolume(aerogelVol );
            // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //


            // position of mirror in cylinder coordinate system
            double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
            if (reflect_parameters)
                mirror_abs_pos_z *= -1.0;

            // row 2 is shifted half step size
            double phi_offset = 0 + 0.5 * phistep * (2 == name_row);


            auto cellTr = RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z)*pyramidTr ;
            PlacedVolume cellPV = barrel_cells_gas_envelope.placeVolume(cellVol, cellTr);
            cellDE.setPlacement(cellPV);


            // increase counter
            cellCounter++;

        }
    }
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    barrel_cells_vessel_envelope.placeVolume(barrel_cells_gas_envelope);
    PlacedVolume assemblyPV = motherVol.placeVolume(barrel_cells_vessel_envelope);
    assemblyPV.addPhysVolID("system", detID);
    det.setPlacement(assemblyPV);
    return det;
}


DECLARE_DETELEMENT(ARCBARREL_o1_v01_T, create_arc_barrel_cell)

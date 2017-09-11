//====================================================================
//  DD4hep Geometry driver for HcalBarrel
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:  $
//====================================================================
#include "DD4hep/Printout.h"
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/TiledLayerGridXY.h"
#include "LcgeoExceptions.h"

#include <iostream>
#include <vector>

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::IntersectionSolid;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Trapezoid;
using dd4hep::Tube;
using dd4hep::UnionSolid;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {

  double boundarySafety = 0.0001;

  xml_det_t   x_det       = element;
  string      det_name    = x_det.nameStr();
  DetElement  sdet( det_name,x_det.id() );

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;
  
  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  xml_comp_t    x_staves          = x_det.staves();
  Material      stavesMaterial    = theDetector.material(x_staves.materialStr());
  Material      air               = theDetector.air();

  PlacedVolume pv;

  sens.setType("calorimeter");

//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build HcalBarrel
//
//====================================================================
  double      Hcal_inner_radius   = theDetector.constant<double>("Hcal_inner_radius");
  double      Hcal_outer_radius   = theDetector.constant<double>("Hcal_outer_radius");
  double      Hcal_half_length    = theDetector.constant<double>("Hcal_half_length");
  int         Hcal_inner_symmetry = theDetector.constant<int>("Hcal_inner_symmetry");
  int         Hcal_outer_symmetry = 0; // Fixed shape for Tube, and not allow to modify from compact xml.

  double      Hcal_radiator_thickness          = theDetector.constant<double>("Hcal_radiator_thickness");
  double      Hcal_chamber_thickness           = theDetector.constant<double>("Hcal_chamber_thickness");
  double      Hcal_back_plate_thickness        = theDetector.constant<double>("Hcal_back_plate_thickness");
  double      Hcal_lateral_plate_thickness     = theDetector.constant<double>("Hcal_lateral_structure_thickness");
  double      Hcal_stave_gaps                  = theDetector.constant<double>("Hcal_stave_gaps");
  double      Hcal_modules_gap                 = theDetector.constant<double>("Hcal_modules_gap"); 
  double      Hcal_middle_stave_gaps           = theDetector.constant<double>("Hcal_middle_stave_gaps");
  double      Hcal_layer_air_gap               = theDetector.constant<double>("Hcal_layer_air_gap");

  int         Hcal_nlayers                     = theDetector.constant<int>("Hcal_nlayers");
  int         Hcal_stair_ntiles                = theDetector.constant<int>("Hcal_stair_ntiles");

  double      TPC_outer_radius                = theDetector.constant<double>("TPC_outer_radius");

  double      Ecal_outer_radius               = theDetector.constant<double>("Ecal_outer_radius");

  printout( dd4hep::DEBUG,  "SHcalSc04_Barrel_v03", "TPC_outer_radius : %e   - Ecal_outer_radius: %e ", TPC_outer_radius , Ecal_outer_radius) ;

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  // check if we have a TiledLayerGridXY segmentation :
  dd4hep::DDSegmentation::TiledLayerGridXY* tileSeg = 
    dynamic_cast< dd4hep::DDSegmentation::TiledLayerGridXY*>( seg.segmentation() ) ;



  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];

  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::BarrelLayout ;
  caloData->inner_symmetry = Hcal_inner_symmetry  ;
  caloData->outer_symmetry = Hcal_outer_symmetry  ;
  caloData->phi0 = 0 ; // fg: also hardcoded below 

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = Hcal_inner_radius ;
  caloData->extent[1] = Hcal_outer_radius ;
  caloData->extent[2] = 0. ; // Barrel zmin is "0" by default.
  caloData->extent[3] = Hcal_half_length ;

//====================================================================
//
// general calculated parameters
//
//====================================================================

  double Hcal_total_dim_y   = Hcal_nlayers * (Hcal_radiator_thickness + Hcal_chamber_thickness) 
                            + Hcal_back_plate_thickness;
  // Hcal New layout design,
  int Hcal_nstairs = 3; //  not allow to be changed for this layout.
  double Hcal_stair_dim_y   = Hcal_total_dim_y / Hcal_nstairs ;
  double Hcal_stair_dim_z   = cell_sizeY * Hcal_stair_ntiles;

  double Hcal_y_dim1_for_x  = Hcal_outer_radius*cos(M_PI/Hcal_inner_symmetry) - Hcal_inner_radius;
  double Hcal_bottom_dim_x  = 2.*(Hcal_inner_radius)*tan(M_PI/Hcal_inner_symmetry)- Hcal_stave_gaps;
  double Hcal_normal_dim_z  = (2 * Hcal_half_length - Hcal_modules_gap)/2.;

 //only the middle has the steel plate.
  double Hcal_regular_chamber_dim_z = Hcal_normal_dim_z - Hcal_lateral_plate_thickness;

// ========= Create Hcal Barrel stave for module 1 ==========================
//  It will be the volume for palcing the Hcal Barrel Chamber(i.e. Layers).
//  The outside will be shaped by tube with maxium Hcal radius.
// ==========================================================================

  double THX1_DDx1   = (Hcal_stair_dim_y + Hcal_inner_radius)*tan(M_PI/Hcal_inner_symmetry);
  double BHX1_DDx2   = Hcal_inner_radius*tan(M_PI/Hcal_inner_symmetry);
  double DHZ1_DDy12  = (Hcal_normal_dim_z+Hcal_stair_dim_z - Hcal_lateral_plate_thickness) / 2.;
  double YXH1_DDz    = Hcal_stair_dim_y / 2.;
  Trapezoid m1stave_bottom(  THX1_DDx1, BHX1_DDx2, DHZ1_DDy12, DHZ1_DDy12, YXH1_DDz);

  double THX2_DDx1   = (Hcal_stair_dim_y*2 + Hcal_inner_radius)*tan(M_PI/Hcal_inner_symmetry);
  double BHX2_DDx2   = (Hcal_inner_radius+Hcal_stair_dim_y)*tan(M_PI/Hcal_inner_symmetry);
  double DHZ2_DDy12  = (Hcal_normal_dim_z - Hcal_lateral_plate_thickness) / 2.;
  double YXH2_DDz    = Hcal_stair_dim_y / 2.;
  Trapezoid m1stave_middle(  THX2_DDx1, BHX2_DDx2, DHZ2_DDy12, DHZ2_DDy12, YXH2_DDz);

  double THX3_DDx1   = (Hcal_stair_dim_y*3 + Hcal_inner_radius)*tan(M_PI/Hcal_inner_symmetry);
  double BHX3_DDx2   = (Hcal_inner_radius+Hcal_stair_dim_y*2)*tan(M_PI/Hcal_inner_symmetry);
  double DHZ3_DDy12  = (Hcal_normal_dim_z-Hcal_stair_dim_z - Hcal_lateral_plate_thickness) / 2.;
  double YXH3_DDz    = Hcal_stair_dim_y/ 2.;
  Trapezoid m1stave_top(  THX3_DDx1, BHX3_DDx2, DHZ3_DDy12, DHZ3_DDy12, YXH3_DDz);

  Position pos1(0,Hcal_stair_dim_z/2.0,Hcal_stair_dim_y-boundarySafety);
  Position pos2(0,Hcal_stair_dim_z/2.0,Hcal_stair_dim_y-boundarySafety*2.0);

  UnionSolid m1stave_bottom_middle(m1stave_middle,m1stave_bottom,pos1);
  UnionSolid m1stave_BMT(m1stave_top,m1stave_bottom_middle,pos2);

// ========= Create Hcal Barrel stave for module 2 ==========================
//  It will be the volume for palcing the Hcal Barrel Chamber(i.e. Layers).
//  The outside will be shaped by tube with maxium Hcal radius.
// ==========================================================================

  DHZ1_DDy12 = (Hcal_normal_dim_z-Hcal_stair_dim_z - Hcal_lateral_plate_thickness) / 2.;
  Trapezoid m2stave_bottom(  THX1_DDx1, BHX1_DDx2, DHZ1_DDy12, DHZ1_DDy12, YXH1_DDz);

  Trapezoid m2stave_middle(  THX2_DDx1, BHX2_DDx2, DHZ2_DDy12, DHZ2_DDy12, YXH2_DDz);

  DHZ3_DDy12  = (Hcal_normal_dim_z+Hcal_stair_dim_z - Hcal_lateral_plate_thickness) / 2.;
  Trapezoid m2stave_top(  THX3_DDx1, BHX3_DDx2, DHZ3_DDy12, DHZ3_DDy12, YXH3_DDz);

  UnionSolid m2stave_bottom_middle(m2stave_middle,m2stave_bottom,pos1);
  UnionSolid m2stave_BMT(m2stave_top,m2stave_bottom_middle,pos2);
 
// ========= Create Hcal Barrel stave   ====================================
//  It will be the volume for palcing the Hcal Barrel Chamber(i.e. Layers).
//  Itself will be placed into the world volume.
// ==========================================================================
 
  double chambers_y_off_correction = 0.;
  double chambers_z_off_correction = 0.;
 
  chambers_y_off_correction = Hcal_stair_dim_y;

  Tube solidCaloTube(0, Hcal_outer_radius, Hcal_normal_dim_z);
  
  RotationZYX mrot(0,0,M_PI/2.);

  Rotation3D mrot3D(mrot);
  Position mxyzVec(0,0,(Hcal_inner_radius + Hcal_total_dim_y / 2. + Hcal_stair_dim_y));
  Transform3D mtran3D(mrot3D,mxyzVec);

  IntersectionSolid barrelModule1Solid(m1stave_BMT, solidCaloTube, mtran3D);
  Volume  EnvLogHcalModule1Barrel(det_name+"_module1",
				  barrelModule1Solid,
				  stavesMaterial);

  IntersectionSolid barrelModule2Solid(m2stave_BMT, solidCaloTube, mtran3D);
  Volume  EnvLogHcalModule2Barrel(det_name+"_module2",
				  barrelModule2Solid,
				  stavesMaterial);

  EnvLogHcalModule1Barrel.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  EnvLogHcalModule2Barrel.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());





  //stave modules lateral plate  parameters
  double rmin_LP  = Hcal_inner_radius;
  double rmax_LP  = Hcal_inner_radius + Hcal_stair_dim_y;
  double zlen_LP  = Hcal_lateral_plate_thickness/2.0; 

  PolyhedraRegular lateral_plate_bottom(Hcal_inner_symmetry, rmin_LP, rmax_LP, zlen_LP);

  rmin_LP = rmin_LP + Hcal_stair_dim_y;
  rmax_LP = rmax_LP + Hcal_stair_dim_y;
  PolyhedraRegular lateral_plate_middle(Hcal_inner_symmetry, rmin_LP, rmax_LP, zlen_LP);

  rmin_LP = rmin_LP + Hcal_stair_dim_y;
  rmax_LP = rmax_LP + Hcal_stair_dim_y;
  PolyhedraRegular lateral_plate_top(Hcal_inner_symmetry, rmin_LP, rmax_LP, zlen_LP);

  RotationZYX Tube_LP_rot(M_PI/2., 0, 0);
  Tube solidCaloTube_LP(0, Hcal_outer_radius, zlen_LP+boundarySafety);

  IntersectionSolid Module_lateral_plate_top(lateral_plate_top, solidCaloTube_LP, Tube_LP_rot);
  Volume  EnvLogHcalModuleBarrelTop_LP(det_name+"_Module_lateral_plate_top",
				       Module_lateral_plate_top,
				       stavesMaterial);

  IntersectionSolid Module_lateral_plate_middle(lateral_plate_middle, solidCaloTube_LP, Tube_LP_rot);
  Volume  EnvLogHcalModuleBarrelMiddle_LP(det_name+"_Module_lateral_plate_middle",
					  Module_lateral_plate_middle,
					  stavesMaterial);

  IntersectionSolid Module_lateral_plate_bottom(lateral_plate_bottom, solidCaloTube_LP, Tube_LP_rot);
  Volume  EnvLogHcalModuleBarrelBottom_LP(det_name+"_Module_lateral_plate_bottom",
					  Module_lateral_plate_bottom,
					  stavesMaterial);

  EnvLogHcalModuleBarrelTop_LP.setAttributes(theDetector,x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  EnvLogHcalModuleBarrelMiddle_LP.setAttributes(theDetector,x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  EnvLogHcalModuleBarrelBottom_LP.setAttributes(theDetector,x_det.regionStr(), x_det.limitsStr(), x_det.visStr());



#ifdef SHCALSC04_DEBUG
  std::cout<< " ==> Hcal_outer_radius: "<<Hcal_outer_radius <<std::endl;
#endif





//====================================================================
//
// Chambers in the HCAL BARREL
//
//====================================================================
  // Build Layer Chamber fill with air, which include the tolarance space at the two side
  // place the slice into the Layer Chamber
  // Place the Layer Chamber into the Stave module
  // place the Stave module into the asembly Volume
  // place the module middle lateral plate into asembly Volume

  double x_length; //dimension of an Hcal barrel layer on the x-axis
  double y_height; //dimension of an Hcal barrel layer on the y-axis
  double m1_z_width;  //dimension of an Hcal barrel layer on the z-axis
  double m2_z_width;  //dimension of an Hcal barrel layer on the z-axis
  x_length = 0.; // Each Layer the x_length has new value.
  y_height = Hcal_chamber_thickness / 2.;

  double xOffset = 0.;//the x_length of a barrel layer is calculated as a
  //barrel x-dimension plus (bottom barrel) or minus
  //(top barrel) an x-offset, which depends on the angle M_PI/Hcal_inner_symmetry

  double xShift = 0.;//Geant4 draws everything in the barrel related to the 
  //center of the bottom barrel, so we need to shift the layers to
  //the left (or to the right) with the quantity xShift


  //-------------------- start loop over HCAL layers ----------------------

  for (int layer_id = 1; layer_id <= (2*Hcal_nlayers); layer_id++)
    {
 
      double TanPiDiv8 = tan(M_PI/Hcal_inner_symmetry);
      double x_total   = 0.;
      double x_halfLength;
      x_length  = 0.;

      int logical_layer_id = 0;

      if ( (layer_id < Hcal_nlayers)
	   || (layer_id > Hcal_nlayers && layer_id < (2*Hcal_nlayers)) )
	logical_layer_id = layer_id % Hcal_nlayers;
      else if ( (layer_id == Hcal_nlayers) 
		|| (layer_id == 2*Hcal_nlayers) ) logical_layer_id = Hcal_nlayers;

      if( logical_layer_id <= Hcal_nlayers/Hcal_nstairs) {
	m1_z_width  = (Hcal_regular_chamber_dim_z + Hcal_stair_dim_z)/2.;
	m2_z_width  = (Hcal_regular_chamber_dim_z - Hcal_stair_dim_z)/2.;
        chambers_z_off_correction = Hcal_stair_dim_z;
      }
      else if ( logical_layer_id > Hcal_nlayers/Hcal_nstairs*2) {
	m1_z_width  = (Hcal_regular_chamber_dim_z - Hcal_stair_dim_z)/2.;
	m2_z_width  = (Hcal_regular_chamber_dim_z + Hcal_stair_dim_z)/2.;
	chambers_z_off_correction = 0;
      }
      else { 
	m1_z_width  = Hcal_regular_chamber_dim_z/2.;
	m2_z_width  = Hcal_regular_chamber_dim_z/2.;
	chambers_z_off_correction = Hcal_stair_dim_z/2.0;
      }
      //---- bottom barrel------------------------------------------------------------
      if( logical_layer_id *(Hcal_radiator_thickness + Hcal_chamber_thickness)
	  < (Hcal_outer_radius * cos(M_PI/Hcal_inner_symmetry) - Hcal_inner_radius ) ) {
	xOffset = (logical_layer_id * Hcal_radiator_thickness 
		   + (logical_layer_id -1) * Hcal_chamber_thickness) * TanPiDiv8;

	x_total  = Hcal_bottom_dim_x/2 - Hcal_middle_stave_gaps/2 + xOffset;
	x_length = x_total - 2*Hcal_layer_air_gap;
	x_halfLength = x_length/2.;

      } else {//----- top barrel -------------------------------------------------
	double y_layerID = logical_layer_id * (Hcal_radiator_thickness + Hcal_chamber_thickness) + Hcal_inner_radius;
	double ro_layer = Hcal_outer_radius - Hcal_radiator_thickness;
	
	x_total = sqrt( ro_layer * ro_layer - y_layerID * y_layerID);
	
	x_length = x_total - Hcal_middle_stave_gaps;
	
	x_halfLength = x_length/2.;
	
	xOffset = (logical_layer_id * Hcal_radiator_thickness 
		   + (logical_layer_id - 1) * Hcal_chamber_thickness - Hcal_y_dim1_for_x) / TanPiDiv8
	  + Hcal_chamber_thickness / TanPiDiv8;
	
      }

      double xAbsShift = (Hcal_middle_stave_gaps/2 + Hcal_layer_air_gap + x_halfLength);
      
      if (layer_id <= Hcal_nlayers)     xShift = - xAbsShift;
      else if (layer_id > Hcal_nlayers) xShift = xAbsShift;

      x_length = x_length/2.;

      LayeredCalorimeterData::Layer caloLayer ;
      caloLayer.cellSize0 = cell_sizeX;
      caloLayer.cellSize1 = cell_sizeY;


      //--------------------------------------------------------------------------------
      //  build chamber box, with the calculated dimensions 
      //-------------------------------------------------------------------------------
      printout( dd4hep::DEBUG,  "SHcalSc04_Barrel_v03", " \n Start to build layer chamber - layer_id: %d", layer_id ) ;
      printout( dd4hep::DEBUG,  "SHcalSc04_Barrel_v03"," m1 chamber x:y:z:  %e:%e:%e", x_length*2., m1_z_width*2. , y_height*2. );
      printout( dd4hep::DEBUG,  "SHcalSc04_Barrel_v03"," m2 chamber x:y:z:  %e:%e:%e", x_length*2., m2_z_width*2. , y_height*2. );

      //check if x_length (scintillator length) is divisible with x_integerTileSize
      if( layer_id <= Hcal_nlayers) {
	double fracPart, intPart;
	double temp = x_length*2./cell_sizeX;
	fracPart = modf(temp, &intPart);
	int noOfIntCells = int(temp);

	tileSeg->setBoundaryLayerX(x_length);

 	if (fracPart == 0){ //divisible
	  if ( noOfIntCells%2 ) {
	    if( tileSeg !=0 ) tileSeg->setLayerOffsetX(0);
	  }
	  else {
	    if( tileSeg !=0 ) tileSeg->setLayerOffsetX(1);
	  }
	  tileSeg->setFractCellSizeXPerLayer(0);
	}
	else if (fracPart>0){
	  if ( noOfIntCells%2 ) {
	    if( tileSeg !=0 ) tileSeg->setLayerOffsetX(1);
	  }
	  else {
	    if( tileSeg !=0 ) tileSeg->setLayerOffsetX(0);
	  }
	  tileSeg->setFractCellSizeXPerLayer( (fracPart+1.0)/2.0*cell_sizeX );
	}

	if ( (int)( (m1_z_width*2.) / cell_sizeX)%2 ){
	  if( tileSeg !=0 ) tileSeg->setLayerOffsetY(0);
	}
	else {
	  if( tileSeg !=0 ) tileSeg->setLayerOffsetY(1);
	}
      }

      //x + air gaps at two side, do not need to build air gaps individualy.
      Box m1ChamberSolid((x_length + Hcal_layer_air_gap), 
			 m1_z_width,
			 y_height);

      string m1ChamberLogical_name      = det_name+_toString(layer_id,"m1_layer%d");
      Volume m1ChamberLogical(m1ChamberLogical_name, m1ChamberSolid, air);   

      Box m2ChamberSolid((x_length + Hcal_layer_air_gap), 
			 m2_z_width,
			 y_height);

      string m2ChamberLogical_name      = det_name+_toString(layer_id,"m2_layer%d");
      Volume m2ChamberLogical(m2ChamberLogical_name, m2ChamberSolid, air);   



      double layer_thickness = y_height*2.;

      double nRadiationLengths=0.;
      double nInteractionLengths=0.;
      double thickness_sum=0;

      nRadiationLengths   = Hcal_radiator_thickness/(stavesMaterial.radLength());
      nInteractionLengths = Hcal_radiator_thickness/(stavesMaterial.intLength());
      thickness_sum       = Hcal_radiator_thickness;

//====================================================================
// Create Hcal Barrel Chamber without radiator
// Place into the Hcal Barrel stave, after each radiator 
//====================================================================
      xml_coll_t c(x_det,_U(layer));
      xml_comp_t   x_layer = c;
      string layer_name      = det_name+_toString(layer_id,"_layer%d");
      
      // Create the slices (sublayers) within the Hcal Barrel Chamber.
      double slice_pos_z = layer_thickness/2. ;
      int slice_number = 0;

      for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
	xml_comp_t x_slice = k;
	string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
	double   slice_thickness = x_slice.thickness();
	Material slice_material  = theDetector.material(x_slice.materialStr());
	DetElement m1_slice(layer_name,_toString(slice_number,"m1_slice%d"),x_det.id());
	DetElement m2_slice(layer_name,_toString(slice_number,"m2_slice%d"),x_det.id());
	
	slice_pos_z -= slice_thickness/2.;
	
	// Slice volume & box
	Volume m1_slice_vol(slice_name,Box(x_length,m1_z_width,slice_thickness/2.),slice_material);
	Volume m2_slice_vol(slice_name,Box(x_length,m2_z_width,slice_thickness/2.),slice_material);

	nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
	nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
	thickness_sum       += slice_thickness/2;
	
	if ( x_slice.isSensitive() ) {
	  m1_slice_vol.setSensitiveDetector(sens);
	  m2_slice_vol.setSensitiveDetector(sens);

#if DD4HEP_VERSION_GE( 0, 15 )
	  //Store "inner" quantities
	  caloLayer.inner_nRadiationLengths   = nRadiationLengths;
	  caloLayer.inner_nInteractionLengths = nInteractionLengths;
	  caloLayer.inner_thickness = thickness_sum;
	  //Store scintillator thickness
	  caloLayer.sensitive_thickness = slice_thickness;
#endif
	  //Reset counters to measure "outside" quantitites
	  nRadiationLengths=0.;
	  nInteractionLengths=0.;
	  thickness_sum = 0.;
	}

	nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
	nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
	thickness_sum       += slice_thickness/2;


	// Set region, limitset, and vis.
	m1_slice_vol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	m2_slice_vol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	// slice PlacedVolume
	PlacedVolume m1_slice_phv = m1ChamberLogical.placeVolume(m1_slice_vol,Position(0.,0.,slice_pos_z));
	PlacedVolume m2_slice_phv = m2ChamberLogical.placeVolume(m2_slice_vol,Position(0.,0.,slice_pos_z));
	if ( x_slice.isSensitive() ) {
	  int slice_id  = (layer_id > Hcal_nlayers)? 1:-1;
	  m1_slice_phv.addPhysVolID("layer",logical_layer_id).addPhysVolID("tower",slice_id);
	  m2_slice_phv.addPhysVolID("layer",logical_layer_id).addPhysVolID("tower",slice_id);
	  printout( dd4hep::DEBUG,  "SHcalSc04_Barrel_v03", "  logical_layer_id:  %d  slice_id:  %d", logical_layer_id, slice_id  ) ;
	}
	
	m1_slice.setPlacement(m1_slice_phv);
	m2_slice.setPlacement(m2_slice_phv);
	// Increment x position for next slice.
	slice_pos_z -= slice_thickness/2.;
	// Increment slice number.
	++slice_number;             
      }

#if DD4HEP_VERSION_GE( 0, 15 )
      //Store "outer" quantities
      caloLayer.outer_nRadiationLengths = nRadiationLengths;
      caloLayer.outer_nInteractionLengths = nInteractionLengths;
      caloLayer.outer_thickness = thickness_sum;
#endif  

//---------------------------  Chamber Placements -----------------------------------------

      double chamber_x_offset, chamber_y_offset, chamber_z_offset;
      chamber_x_offset = xShift;

      chamber_z_offset = 0;

      chamber_y_offset = -(-Hcal_total_dim_y/2. 
			   + (logical_layer_id-1) *(Hcal_chamber_thickness + Hcal_radiator_thickness)
			   + Hcal_radiator_thickness + Hcal_chamber_thickness/2.);

      
      pv =  EnvLogHcalModule1Barrel.placeVolume(m1ChamberLogical,
					       Position(chamber_x_offset,
							chamber_z_offset + chambers_z_off_correction,
							chamber_y_offset + chambers_y_off_correction));
      
      pv =  EnvLogHcalModule2Barrel.placeVolume(m2ChamberLogical,
					       Position(chamber_x_offset,
							chamber_z_offset + chambers_z_off_correction,
							chamber_y_offset + chambers_y_off_correction));
      


      //-----------------------------------------------------------------------------------------
      if( layer_id <= Hcal_nlayers ){ // add the first set of layers to the reconstruction data
	
	caloLayer.distance = Hcal_inner_radius + Hcal_total_dim_y/2.0 - (chamber_y_offset + chambers_y_off_correction)
	  - caloLayer.inner_thickness ;
	caloLayer.absorberThickness = Hcal_radiator_thickness ;
	
	caloData->layers.push_back( caloLayer ) ;
      }
      //-----------------------------------------------------------------------------------------

      
    }//end loop over HCAL nlayers;

  if( tileSeg !=0 ){
    // check the offsets directly in the TileSeg ...
    std::vector<double> LOX = tileSeg->layerOffsetX();
    std::vector<double> LOY = tileSeg->layerOffsetY();

    std::stringstream sts ;
    sts <<" layerOffsetX(): ";
    for (std::vector<double>::const_iterator ix = LOX.begin(); ix != LOX.end(); ++ix)
      sts << *ix << ' ';
    printout( dd4hep::DEBUG,  "SHcalSc04_Barrel_v03", "%s" , sts.str().c_str() ) ;
    sts.clear() ; sts.str("") ;
    sts <<" layerOffsetY(): ";
    for (std::vector<double>::const_iterator iy = LOY.begin(); iy != LOY.end(); ++iy)
      sts << *iy << ' ';
    printout( dd4hep::DEBUG,  "SHcalSc04_Barrel_v03", "%s" , sts.str().c_str() ) ;

  }



//====================================================================
// Place HCAL Barrel stave module into the envelope
//====================================================================
  double stave_phi_offset,  module_z_offset,  lateral_plate_z_offset;

  double Y = Hcal_inner_radius + Hcal_total_dim_y / 2. + Hcal_stair_dim_y;

  stave_phi_offset = M_PI/Hcal_inner_symmetry -M_PI/2.;


  //-------- start loop over HCAL BARREL staves ----------------------------

  for (int stave_id = 1;
       stave_id <=Hcal_inner_symmetry;
       stave_id++)
    {
      module_z_offset = - (Hcal_normal_dim_z + Hcal_modules_gap + Hcal_lateral_plate_thickness)/2.;

      double phirot = stave_phi_offset;
      RotationZYX srot(0,phirot,M_PI*0.5);
      Rotation3D srot3D(srot);

      int module_id = 1;
      Position m1sxyzVec(-Y*sin(phirot), Y*cos(phirot), module_z_offset - Hcal_stair_dim_z/2.0);
      Transform3D m1stran3D(srot3D,m1sxyzVec);
	  
      // Place Hcal Barrel volume into the envelope volume for module 1
      pv = envelope.placeVolume(EnvLogHcalModule1Barrel,m1stran3D);
      pv.addPhysVolID("stave",stave_id);
      pv.addPhysVolID("module",module_id);

      module_id = 2;
      Position m2sxyzVec(-Y*sin(phirot), Y*cos(phirot), -module_z_offset - Hcal_stair_dim_z/2.0);
      Transform3D m2stran3D(srot3D,m2sxyzVec);
	  
      // Place Hcal Barrel volume into the envelope volume for module 2
      pv = envelope.placeVolume(EnvLogHcalModule2Barrel,m2stran3D);
      pv.addPhysVolID("stave",stave_id);
      pv.addPhysVolID("module",module_id);

      stave_phi_offset -=  M_PI*2.0/Hcal_inner_symmetry;

    }  //-------- end loop over HCAL BARREL staves ----------------------------

  lateral_plate_z_offset = - (Hcal_lateral_plate_thickness + Hcal_modules_gap)/2.;

  // place lateral plate for module 1 and 2.
  for (int m_id = 1; m_id <=2;m_id++) {
    lateral_plate_z_offset = lateral_plate_z_offset - Hcal_stair_dim_z;
    pv = envelope.placeVolume(EnvLogHcalModuleBarrelTop_LP,
			      Position(0,0,lateral_plate_z_offset) );
    
    lateral_plate_z_offset = lateral_plate_z_offset + Hcal_stair_dim_z;
    pv = envelope.placeVolume(EnvLogHcalModuleBarrelMiddle_LP,
			      Position(0,0,lateral_plate_z_offset));
    
    lateral_plate_z_offset = lateral_plate_z_offset + Hcal_stair_dim_z;
    pv = envelope.placeVolume(EnvLogHcalModuleBarrelBottom_LP,
			      Position(0,0,lateral_plate_z_offset));
    
    lateral_plate_z_offset = (Hcal_lateral_plate_thickness + Hcal_modules_gap)/2.;
  }

  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;

 
  return sdet;

}

DECLARE_DETELEMENT(SHcalSc04_Barrel_v03, create_detector)

//====================================================================
//  DD4hep Geometry driver for HcalBarrel
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:  $
//====================================================================
#include "DD4hep/Printout.h"
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/TiledLayerGridXY.h"
#include "LcgeoExceptions.h"

#include <iostream>
#include <vector>

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DDSegmentation::TiledLayerGridXY;
using dd4hep::DEBUG;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::ERROR;
using dd4hep::IntersectionSolid;
using dd4hep::Material;
using dd4hep::PlacedVolume;
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
using dd4hep::Volume;
using dd4hep::WARNING;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

using lcgeo::GeometryException;

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

// After reading in all the necessary parameters.
// To check the radius range and the space for placing the total layers
static bool validateEnvelope(double rInner, double rOuter, double radiatorThickness, double layerThickness, int layerNumber){
  
  bool Error = false;
  bool Warning = false;
  double spaceAllowed = rOuter*cos(M_PI/16.) - rInner;
  double spaceNeeded  = (radiatorThickness + layerThickness)* layerNumber;
  double spaceToleranted  = (radiatorThickness + layerThickness)* (layerNumber+1);
  double rOuterRecommaned = ( rInner + spaceNeeded )/cos(M_PI/16.);
  int layerNumberRecommaned = floor( ( spaceAllowed ) / (radiatorThickness + layerThickness) );
  
  
  if( spaceNeeded > spaceAllowed )
    {
      printout( ERROR,  "SHcalSc04_Barrel_v01", " Layer number is more than it can be built! "  ) ;
      Error = true;
    }
  else if ( spaceToleranted < spaceAllowed )
    {
      printout( WARNING,  "SHcalSc04_Barrel_v01", " Layer number is less than it is able to build!" ) ;
      Warning = true;
    }
  else
    {
      printout( DEBUG,  "SHcalSc04_Barrel_v01"," has been validated and start to build it." ) ;
      Error = false;
      Warning = false;
    }

  if( Error )
    {
      cout<<"\n ============> First Help Documentation <=============== \n"
	  <<" When you see this message, that means you are crashing the module. \n"
	  <<" Please take a cup of cafe, and think about what you want to do! \n"
	  <<" Here are few FirstAid# for you. Please read them carefully. \n"
	  <<" \n"
	  <<" ###  FirstAid 1: ###\n"
	  <<" If you want to build HCAL within the rInner and rOuter range, \n"
	  <<" please reduce the layer number to " << layerNumberRecommaned <<" \n"
	  <<" with the current layer thickness structure. \n"
	  <<" \n"
	  <<" You may redisgn the layer structure and thickness, too. \n"
	  <<" \n"
	  <<" ###  FirstAid 2: ###\n"
	  <<" If you want to build HCAL with this layer number and the layer thickness, \n"
	  <<" you have to update rOuter to "<< rOuterRecommaned*10. <<"*mm \n"
	  <<" and to inform other subdetector, you need this space for building your design. \n"
	  <<" \n"
	  <<" ###  FirstAid 3: ###\n"
	  <<" Do you think that you are looking for another type of HCAL driver? \n"
	  <<" \n"
	  <<endl; 
       throw GeometryException(  "SHcalSc04_Barrel: Error: Layer number is more than it can be built!"   ) ;
    } 
  else if( Warning )
    {
      cout<<"\n ============> First Help Documentation <=============== \n"
	  <<" When you see this warning message, that means you are changing the module. \n"
	  <<" Please take a cup of cafe, and think about what you want to do! \n"
	  <<" Here are few FirstAid# for you. Please read them carefully. \n"
	  <<" \n"
	  <<" ###  FirstAid 1: ###\n"
	  <<" If you want to build HCAL within the rInner and rOuter range, \n"
	  <<" You could build the layer number up to " << layerNumberRecommaned <<" \n"
	  <<" with the current layer thickness structure. \n"
	  <<" \n"
	  <<" You may redisgn the layer structure and thickness, too. \n"
	  <<" \n"
	  <<" ###  FirstAid 2: ###\n"
	  <<" If you want to build HCAL with this layer number and the layer thickness, \n"
	  <<" you could reduce rOuter to "<< rOuterRecommaned*10. <<"*mm \n"
	  <<" and to reduce the back plate thickness, which you may not need for placing layer. \n"
	  <<" \n"
	  <<" ###  FirstAid 3: ###\n"
	  <<" Do you think that you are looking for another type of HCAL driver? \n"
	  <<" \n"
	  <<endl; 
      return Warning;
    }
  else { return true; }

}

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
  //double      Hcal_cells_size                  = theDetector.constant<double>("Hcal_cells_size");

  int         Hcal_nlayers                     = theDetector.constant<int>("Hcal_nlayers");

  double      TPC_outer_radius               = theDetector.constant<double>("TPC_outer_radius");

  double      Ecal_outer_radius               = theDetector.constant<double>("Ecal_outer_radius");

  printout( DEBUG,  "SHcalSc04_Barrel_v01", "TPC_outer_radius : %e   - Ecal_outer_radius: %e ", TPC_outer_radius , Ecal_outer_radius) ;

  validateEnvelope(Hcal_inner_radius, Hcal_outer_radius, Hcal_radiator_thickness, Hcal_chamber_thickness, Hcal_nlayers);

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  // check if we have a TiledLayerGridXY segmentation :
  TiledLayerGridXY* tileSeg = 
    dynamic_cast< TiledLayerGridXY*>( seg.segmentation() ) ;



  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
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

  double Hcal_y_dim1_for_x  = Hcal_outer_radius*cos(M_PI/Hcal_inner_symmetry) - Hcal_inner_radius;
  double Hcal_bottom_dim_x  = 2.*Hcal_inner_radius*tan(M_PI/Hcal_inner_symmetry)- Hcal_stave_gaps;
  double Hcal_normal_dim_z  = (2 * Hcal_half_length - Hcal_modules_gap)/2.;

 //only the middle has the steel plate.
  double Hcal_regular_chamber_dim_z = Hcal_normal_dim_z - Hcal_lateral_plate_thickness;

  //double Hcal_cell_dim_x            = Hcal_cells_size;
  //double Hcal_cell_dim_z            = Hcal_regular_chamber_dim_z / floor (Hcal_regular_chamber_dim_z/Hcal_cell_dim_x);

 
// ========= Create Hcal Barrel stave   ====================================
//  It will be the volume for palcing the Hcal Barrel Chamber(i.e. Layers).
//  Itself will be placed into the world volume.
// ==========================================================================
 
  double chambers_y_off_correction = 0.;
 
  // stave modules shaper parameters
  double BHX  = (Hcal_bottom_dim_x + Hcal_stave_gaps)/2.;
  double THX  = (Hcal_total_dim_y + Hcal_inner_radius)*tan(M_PI/Hcal_inner_symmetry);
  double YXH  = Hcal_total_dim_y / 2.;
  double DHZ  = (Hcal_normal_dim_z - Hcal_lateral_plate_thickness) / 2.;

  Trapezoid stave_shaper(  THX, BHX, DHZ, DHZ, YXH);

  Tube solidCaloTube(0, Hcal_outer_radius, DHZ+boundarySafety);
  
  RotationZYX mrot(0,0,M_PI/2.);

  Rotation3D mrot3D(mrot);
  Position mxyzVec(0,0,(Hcal_inner_radius + Hcal_total_dim_y / 2.));
  Transform3D mtran3D(mrot3D,mxyzVec);

  IntersectionSolid barrelModuleSolid(stave_shaper, solidCaloTube, mtran3D);

  Volume  EnvLogHcalModuleBarrel(det_name+"_module",barrelModuleSolid,stavesMaterial);

  EnvLogHcalModuleBarrel.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());





  //stave modules lateral plate shaper parameters
  double BHX_LP  = BHX;
  double THX_LP  = THX;
  double YXH_LP  = YXH;

  //build lateral palte here to simulate lateral plate in the middle of barrel.
  double DHZ_LP  = Hcal_lateral_plate_thickness/2.0; 

  Trapezoid stave_shaper_LP(THX_LP, BHX_LP, DHZ_LP, DHZ_LP, YXH_LP);

  Tube solidCaloTube_LP(0, Hcal_outer_radius, DHZ_LP+boundarySafety);

  IntersectionSolid Module_lateral_plate(stave_shaper_LP, solidCaloTube_LP, mtran3D);

  Volume  EnvLogHcalModuleBarrel_LP(det_name+"_Module_lateral_plate",Module_lateral_plate,stavesMaterial);

  EnvLogHcalModuleBarrel_LP.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());



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
  double z_width;  //dimension of an Hcal barrel layer on the z-axis
  x_length = 0.; // Each Layer the x_length has new value.
  y_height = Hcal_chamber_thickness / 2.;
  z_width  = Hcal_regular_chamber_dim_z/2.;

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

      
      //calculate the size of a fractional tile
      //-> this sets fract_cell_dim_x
      
      //double fract_cell_dim_x = 0.;
      //this->CalculateFractTileSize(2*x_length, Hcal_cell_dim_x, fract_cell_dim_x);
      
      //Vector newFractCellDim(fract_cell_dim_x, Hcal_chamber_thickness, Hcal_cell_dim_z);
      //theBarrilRegSD->SetFractCellDimPerLayer(layer_id, newFractCellDim);

      LayeredCalorimeterData::Layer caloLayer ;
      caloLayer.cellSize0 = cell_sizeX;
      caloLayer.cellSize1 = cell_sizeY;


      //--------------------------------------------------------------------------------
      //  build chamber box, with the calculated dimensions 
      //-------------------------------------------------------------------------------
      printout( DEBUG,  "SHcalSc04_Barrel_v01", " \n Start to build layer chamber - layer_id: %d", layer_id ) ;
      printout( DEBUG,  "SHcalSc04_Barrel_v01"," chamber x:y:z:  %e:%e:%e", x_length*2., z_width*2. , y_height*2. );

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

	if ( (int)( (z_width*2.) / cell_sizeX)%2 ){
	  if( tileSeg !=0 ) tileSeg->setLayerOffsetY(0);
	}
	else {
	  if( tileSeg !=0 ) tileSeg->setLayerOffsetY(1);
	}
      }

      Box ChamberSolid((x_length + Hcal_layer_air_gap),  //x + air gaps at two side, do not need to build air gaps individualy.
			 z_width,   //z attention!
			 y_height); //y attention!

      string ChamberLogical_name      = det_name+_toString(layer_id,"_layer%d");

      Volume ChamberLogical(ChamberLogical_name, ChamberSolid, air);   



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
	DetElement slice(layer_name,_toString(slice_number,"slice%d"),x_det.id());
	
	slice_pos_z -= slice_thickness/2.;
	
	// Slice volume & box
	Volume slice_vol(slice_name,Box(x_length,z_width,slice_thickness/2.),slice_material);

	nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
	nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
	thickness_sum       += slice_thickness/2;
	
	if ( x_slice.isSensitive() ) {
	  slice_vol.setSensitiveDetector(sens);

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
	slice_vol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	// slice PlacedVolume
	PlacedVolume slice_phv = ChamberLogical.placeVolume(slice_vol,Position(0.,0.,slice_pos_z));
	if ( x_slice.isSensitive() ) {
	  int slice_id  = (layer_id > Hcal_nlayers)? 1:-1;
	  slice_phv.addPhysVolID("layer",logical_layer_id).addPhysVolID("tower",slice_id);
	  printout( DEBUG,  "SHcalSc04_Barrel_v01", "  logical_layer_id:  %d  slice_id:  %d", logical_layer_id, slice_id  ) ;
	}
	
	slice.setPlacement(slice_phv);
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


      pv =  EnvLogHcalModuleBarrel.placeVolume(ChamberLogical,
					       Position(chamber_x_offset,
							chamber_z_offset,
							chamber_y_offset + chambers_y_off_correction));
      


      //-----------------------------------------------------------------------------------------
      if( layer_id <= Hcal_nlayers ){ // add the first set of layers to the reconstruction data
	
	caloLayer.distance = Hcal_inner_radius + Hcal_total_dim_y/2.0 - (chamber_y_offset + chambers_y_off_correction)
	  - caloLayer.inner_thickness ; // Will be added later at "DDMarlinPandora/DDGeometryCreator.cc:226" to get center of sensitive element

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
    printout( DEBUG,  "SHcalSc04_Barrel_v01", "%s" , sts.str().c_str() ) ;
    sts.clear() ; sts.str("") ;
    sts <<" layerOffsetY(): ";
    for (std::vector<double>::const_iterator iy = LOY.begin(); iy != LOY.end(); ++iy)
      sts << *iy << ' ';
    printout( DEBUG,  "SHcalSc04_Barrel_v01", "%s" , sts.str().c_str() ) ;

  }



//====================================================================
// Place HCAL Barrel stave module into the envelope
//====================================================================
  double stave_phi_offset,  module_z_offset,  lateral_plate_z_offset;

  double Y = Hcal_inner_radius + Hcal_total_dim_y / 2.;

  stave_phi_offset = M_PI/Hcal_inner_symmetry -M_PI/2.;


  //-------- start loop over HCAL BARREL staves ----------------------------

  for (int stave_id = 1;
       stave_id <=Hcal_inner_symmetry;
       stave_id++)
    {
      module_z_offset = - (Hcal_normal_dim_z + Hcal_modules_gap + Hcal_lateral_plate_thickness)/2.;
      lateral_plate_z_offset = - (Hcal_lateral_plate_thickness + Hcal_modules_gap)/2.;

      double phirot = stave_phi_offset;
      RotationZYX srot(0,phirot,M_PI*0.5);
      Rotation3D srot3D(srot);

      for (int module_id = 1;
         module_id <=2;
         module_id++)
	{
	  Position sxyzVec(-Y*sin(phirot), Y*cos(phirot), module_z_offset);
	  Transform3D stran3D(srot3D,sxyzVec);
	  
	  // Place Hcal Barrel volume into the envelope volume
	  pv = envelope.placeVolume(EnvLogHcalModuleBarrel,stran3D);
	  pv.addPhysVolID("stave",stave_id);
	  pv.addPhysVolID("module",module_id);

	  const int staveCounter = (stave_id-1)*2+module_id-1;
	  DetElement stave(sdet, _toString(staveCounter,"stave%d"), staveCounter );
	  stave.setPlacement(pv);

	  Position xyzVec_LP(-Y*sin(phirot), Y*cos(phirot),lateral_plate_z_offset);
	  Transform3D tran3D_LP(srot3D,xyzVec_LP);
	  pv = envelope.placeVolume(EnvLogHcalModuleBarrel_LP,tran3D_LP);

	  module_z_offset = - module_z_offset;
	  lateral_plate_z_offset = - lateral_plate_z_offset;
	}


      stave_phi_offset -=  M_PI*2.0/Hcal_inner_symmetry;
    }  //-------- end loop over HCAL BARREL staves ----------------------------


  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;

 
  return sdet;

}

DECLARE_DETELEMENT(SHcalSc04_Barrel_v01, create_detector)

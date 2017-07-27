//====================================================================
//  DD4hep Geometry driver for HcalBarrel
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:  $
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "LcgeoExceptions.h"

using namespace std;

using dd4hep::Assembly;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::Cone;
using dd4hep::Trapezoid;
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
using dd4hep::Tube;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

// workaround for DD4hep v00-14 (and older)
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0
#endif
//

// SemiDigital RPC Hcal_Barrel_SD in TESLA geometry

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {
  
  double boundarySafety = 0.0001;

  xml_det_t   x_det       = element;
  string      det_name    = x_det.nameStr();
  int         det_id      = x_det.id();
  DetElement  sdet( det_name,det_id );

  
  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;

  dd4hep::xml::setDetectorTypeFlag(element, sdet);
  
  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  xml_comp_t    x_staves          = x_det.staves();
  Material      stavesMaterial    = theDetector.material(x_staves.materialStr());
  Material      air               = theDetector.air();

  PlacedVolume pv;

  sens.setType("calorimeter");

//  DetElement    module_det("module0",det_id);

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();

  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeZ      = cellSizeVector[1];


  // Some verbose output
  cout << " \n\n\n CREATE DETECTOR: Hcal_Barrel_SD_v02  - TESLA geometry" << endl;




  //====================================================================
  //
  // Read all the constant from ILD_o2_v01.xml
  // Use them to build HcalBarrel
  //
  //====================================================================
  double      Hcal_inner_radius   = theDetector.constant<double>("Hcal_inner_radius");
  double      Hcal_outer_radius   = theDetector.constant<double>("Hcal_outer_radius");
  double      Hcal_half_length    = theDetector.constant<double>("Hcal_half_length");
  int         Hcal_inner_symmetry = theDetector.constant<int>("Hcal_inner_symmetry");
  int         Hcal_outer_symmetry = 0; // Fixed shape for Tube, and not allow to modify from compact xml.

  double      Hcal_radiator_thickness          = theDetector.constant<double>("HcalSD_radiator_thickness");
  double      Hcal_module_wall_thickness        = theDetector.constant<double>("HcalSD_module_wall_thickness");
  double      Hcal_modules_gap                 = theDetector.constant<double>("HcalSD_modules_gap"); 
  double      Hcal_layer_air_gap               = theDetector.constant<double>("HcalSD_layer_air_gap");
  double      Hcal_airgap_thickness             = theDetector.constant<double>("HcalSD_airgap_thickness");
//  double      Hcal_cells_size                  = theDetector.constant<double>("HcalSD_cells_size");
  double      Hcal_barrel_thickness            = theDetector.constant<double>("Hcal_barrel_thickness");
  int         Hcal_MinNumCellsInTransvPlane    = theDetector.constant<int>("HcalSD_MinNumCellsInTransvPlane");
  int         Hcal_barrel_number_modules       = theDetector.constant<int>("HcalBarrelSD_number_modules");

  int         Hcal_nlayers                     = theDetector.constant<int>("HcalBarrelSD_nlayers");

  double      Hcal_stave_gaps                  = theDetector.constant<double>("Hcal_stave_gaps");
  double      Hcal_middle_stave_gaps           = theDetector.constant<double>("Hcal_middle_stave_gaps");

  double      TPC_outer_radius               = theDetector.constant<double>("TPC_outer_radius");
  std::cout << " ***********TPC_outer_radius " << TPC_outer_radius  << std::endl ;

  double      Ecal_outer_radius               = theDetector.constant<double>("Ecal_outer_radius");
  std::cout << " ***********Ecal_outer_radius " <<  Ecal_outer_radius << std::endl ;
  std::cout << " ***********Hcal_inner_radius " <<  Hcal_inner_radius << std::endl ;
  std::cout << " ***********Hcal_outer_radius " <<  Hcal_outer_radius << std::endl ;
  



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

// First get the layer (chamber) thickness caluculated from the Hcal_Endcaps_SD_v0x.xml file
// not from the fixed, defined value "Hcal_chamber_thickness"
  xml_coll_t c(x_det,_U(layer));
  xml_comp_t   x_layer = c;

  double layer_thickness = 0.0;

  for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
     xml_comp_t x_slice = k;
     layer_thickness += x_slice.thickness();
  }
  cout<<" layer_thickness (from slices) = "<<layer_thickness<<endl;

  double Hcal_total_dim_y   = Hcal_nlayers*layer_thickness ;

// trapezoid dimensions
  double Hcal_trap_height_dim_y  = Hcal_outer_radius*cos(M_PI/Hcal_inner_symmetry) - Hcal_inner_radius;
  double Hcal_trap_top_dim_x     = 2.*Hcal_outer_radius*sin(M_PI/Hcal_inner_symmetry)- Hcal_stave_gaps;
  double Hcal_trap_bottom_dim_x  = 2.*Hcal_inner_radius*tan(M_PI/Hcal_inner_symmetry)- Hcal_stave_gaps;
//  remaining space between top of the trapezoid and the tube (outer radius)
  double Hcal_remaining_space = Hcal_outer_radius * (1.0 - cos(M_PI/Hcal_inner_symmetry));

  double Hcal_normal_dim_z  = 2 * Hcal_half_length/Hcal_barrel_number_modules;
  double Hcal_regular_chamber_dim_z = Hcal_normal_dim_z - 2.*Hcal_module_wall_thickness;
// each barrel (wheel) in a "box" with wall thickness of 10mm (in z-direction)
  
  std::cout<< "Hcal_Barrel number of z-modules: "<<Hcal_barrel_number_modules<<std::endl;
  std::cout<< " ==> Hcal_cell_dim_x  (z): "<<cell_sizeX << " " << cell_sizeZ <<std::endl;
  std::cout <<"Hcal_half_length:    "<< Hcal_half_length<< " Hcal_modules_gap: "<<Hcal_modules_gap << endl;
  std::cout <<"Hcal_barrel_thickness:"<< Hcal_barrel_thickness << " Hcal_requested_thickness: "<<Hcal_total_dim_y << endl;
  std::cout <<"Hcal_regular_chamber_dim_z, Hcal_normal_dim_z:  "<< Hcal_regular_chamber_dim_z<<" "<<Hcal_normal_dim_z<<endl;
  std::cout <<"Hcal_MinNumCellsInTransvPlane:  "<< Hcal_MinNumCellsInTransvPlane<<endl;
 
  int Nlayers_trap_height = floor( Hcal_trap_height_dim_y / layer_thickness );
  int Nlayers_remaining   = floor( Hcal_remaining_space / layer_thickness );
  std::cout <<"Hcal_barrel_trapezoid_height: "<< Hcal_trap_height_dim_y <<" # of layers fitting into it: "<<Nlayers_trap_height << endl;
  std::cout <<"Hcal_barrel remainning space: "<<Hcal_remaining_space<<" # of layers fitting into it: "<<Nlayers_remaining<< endl;
  // ========= Create Hcal Barrel stave   ====================================
  //  It will be the volume for placing the Hcal Barrel Chamber(i.e. Layers).
  //  Itself will be placed into the world volume.
  // ==========================================================================
 
  Assembly EnvLogHcalModuleBarrel(det_name+"_module");

  EnvLogHcalModuleBarrel.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());


  double DHZ_LP  = Hcal_module_wall_thickness; 

  // stave modules shaper parameters
  double BHX  = (Hcal_trap_bottom_dim_x + Hcal_stave_gaps)/2.;
  double THX  = (Hcal_total_dim_y + Hcal_inner_radius)*tan(M_PI/Hcal_inner_symmetry);
  double YXH  = Hcal_total_dim_y / 2.;
  double DHZ  = Hcal_regular_chamber_dim_z / 2.;
  
  Trapezoid stave_shaper(  THX, BHX, DHZ, DHZ, YXH);
 
  Tube solidCaloTube(0, Hcal_outer_radius, DHZ+boundarySafety);
                 
//TK  RotationZYX rot(0,0,M_PI/Hcal_inner_symmetry);
  
  RotationZYX mrot(0,0,M_PI/2.0);
  Rotation3D mrot3D(mrot);
  Position mxyzVec(0,(Hcal_inner_radius + Hcal_total_dim_y / 2.),0);
//TK  Position mxyzVec(0,0,(Hcal_inner_radius + Hcal_total_dim_y / 2.));
  Transform3D mtran3D(mrot3D,mxyzVec);
 
  IntersectionSolid barrelModuleSolid(stave_shaper, solidCaloTube,mtran3D);
//  Volume  EnvLogHcalModuleBarrel(det_name+"_module",barrelModuleSolid,stavesMaterial);

//  EnvLogHcalModuleBarrel.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());


  //stave modules lateral plate shaper parameters
  double BHX_LP  = BHX;
  double THX_LP  = THX;
  double YXH_LP  = YXH;
 
  //build lateral plate here to simulate lateral plate in the middle of barrel.
//  double DHZ_LP  = Hcal_lateral_plate_thickness/2.0; 

  
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
  // Build Layer Chamber fill with air, which include the tolerance space at the two sides
  // place the slice into the Layer Chamber
  // Place the Layer Chamber into the Stave module
  // place the Stave module into the asembly Volume
  // place the module middle lateral plate into asembly Volume


  LayeredCalorimeterData::Layer caloLayer ;
  caloLayer.cellSize0 = cell_sizeX;
  caloLayer.cellSize1 = cell_sizeZ;

  //-------------------- start loop over HCAL layers ----------------------
  double AngleRatio=tan(M_PI/Hcal_inner_symmetry);//"k", updated for flexible symmetry
  double d_InnerOctoSize=2*AngleRatio*Hcal_inner_radius - Hcal_stave_gaps;//"d/2"
  cout <<" \n radiator_thickness, airgap_thickness  "<<Hcal_radiator_thickness<<" "<<Hcal_airgap_thickness <<endl ;
  cout <<" \n AngleRatio, d_InnerOctoSize, trap top size  "<< AngleRatio<<" "<<d_InnerOctoSize <<" "<<Hcal_trap_top_dim_x<<endl ;
  cout <<" \n trapezoid bottom & top size:  "<<Hcal_trap_bottom_dim_x <<" "<<Hcal_trap_top_dim_x<<endl ;
  cout <<" \n boundarySafety "<< boundarySafety <<endl ;

  double x_length = 0.; //dimension of an Hcal barrel layer on the x-axis; each layer has a nex x_length value
  double y_height = layer_thickness/2.; //dimension of an Hcal barrel layer on the y-axis
  double z_width  = Hcal_regular_chamber_dim_z/2.;  //dimension of an Hcal barrel layer on the z-axis
  double Hcal_chamber_thickness = layer_thickness - Hcal_radiator_thickness;

  double xOffset = 0.;//the x_length of a barrel layer is calculated as a
  //barrel x-dimension plus (bottom barrel) or minus
  //  //(top barrel) an x-offset, which depends on the angle M_PI/Hcal_inner_symmetry
  //
  double xShift = 0.;//Geant4 draws everything in the barrel related to the 
  //      //center of the bottom barrel, so we need to shift the layers to
  //      //the left (or to the right) with the quantity xShift


  for (int layer_id = 1; layer_id <= 2*Hcal_nlayers; layer_id++)
    {
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
      if( logical_layer_id *layer_thickness < Hcal_trap_height_dim_y ) {
	xOffset = (logical_layer_id -1) *layer_thickness * AngleRatio; 

	x_total  = Hcal_trap_bottom_dim_x/2 - Hcal_middle_stave_gaps/2 + xOffset;
	x_length = x_total - 2*Hcal_layer_air_gap;
	x_halfLength = x_length/2.;
        cout<<"Bottom : layer_id, ncells, x_length, z_width, y_height, x_offset :"<<layer_id<<" "<<floor(x_length/cell_sizeX)<<" "<<x_length<<" "<<z_width<<" "<<y_height<<" "<<xOffset<<endl;
       

      }else {//----- top barrel -------------------------------------------------
        double y_layerID = logical_layer_id * layer_thickness + Hcal_inner_radius;
	double ro_layer = Hcal_outer_radius - Hcal_radiator_thickness;
	
	x_total = sqrt( ro_layer * ro_layer - y_layerID * y_layerID);
	
	xOffset =  (logical_layer_id * Hcal_chamber_thickness - Hcal_trap_height_dim_y) /AngleRatio 
	  + Hcal_chamber_thickness / AngleRatio;
	x_length = x_total - Hcal_middle_stave_gaps/2 - 2*Hcal_layer_air_gap ;
	
	x_halfLength = x_length/2.;
/***	
	xOffset = (logical_layer_id * Hcal_radiator_thickness 
		   + (logical_layer_id - 1) * Hcal_chamber_thickness - Hcal_y_dim1_for_x) / TanPiDiv8
	  + Hcal_chamber_thickness / TanPiDiv8;
***/
        cout<<"Top barrel: layer_id, ncells, x_length, z_width, y_height, xOffset :"<<layer_id<<" "<<floor(x_length/cell_sizeX)<<" "<<x_length<<" "<<z_width<<" "<<y_height<<" "<<xOffset<<endl;
	
      }

      double xAbsShift = (Hcal_middle_stave_gaps/2 + Hcal_layer_air_gap + x_halfLength);
//      if (logical_layer_id <= Hcal_nlayers)     xShift = - xAbsShift;
//      else if (logical_layer_id > Hcal_nlayers) xShift = xAbsShift;
      if (layer_id <= Hcal_nlayers) {
         xShift = - xAbsShift;
      }else if (layer_id > Hcal_nlayers){
         xShift = xAbsShift;
      }
      //--------------------------------------------------------------------------------
      //  build chamber box, with the calculated dimensions 
      //-------------------------------------------------------------------------------
      //x + air gaps at two side, do not need to build air gaps individualy.
      Box ChamberSolid((x_halfLength + Hcal_layer_air_gap),
			 z_width,   //z attention!
                         y_height); //y attention!

      string ChamberLogical_name      = det_name+_toString(layer_id,"_layer%d");

      Volume ChamberLogical(ChamberLogical_name, ChamberSolid, air);   


      //====================================================================
      // Create Hcal Barrel Chamber with radiator
      // Place into the Hcal Barrel stave
      //====================================================================

      double nRadiationLengths=0.;
      double nInteractionLengths=0.;
      double thickness_sum=0;
    
      nRadiationLengths   = Hcal_radiator_thickness/(stavesMaterial.radLength());
      nInteractionLengths = Hcal_radiator_thickness/(stavesMaterial.intLength());

      string layer_name      = det_name +_toString(layer_id,"_layer%d");
      // Create the slices (sublayers) within the Hcal Barrel Chamber.
      double slice_pos_z = layer_thickness/2.;
      int slice_number = 0;
      for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
	xml_comp_t x_slice = k;
	string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
	double   slice_thickness = x_slice.thickness();
	Material slice_material  = theDetector.material(x_slice.materialStr());

	DetElement slice(layer_name,_toString(slice_number,"slice%d"),x_det.id());
	
	slice_pos_z -= slice_thickness/2.;
        nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
        nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
        thickness_sum       += slice_thickness/2;
	
        if(logical_layer_id==1)
           cout<<"  Layer_slice: "<<slice_name<<" slice_thickness: "<< slice_thickness<<" slice_pos_z :"<<slice_pos_z<<  endl;
	// Slice volume & box
        Volume slice_vol(slice_name,Box(x_halfLength,z_width,slice_thickness/2.),slice_material);	
	if ( x_slice.isSensitive() ) {
	  slice_vol.setSensitiveDetector(sens);
#if DD4HEP_VERSION_GE( 0, 15 )
        //Store "inner" quantities
          caloLayer.inner_nRadiationLengths = nRadiationLengths;
          caloLayer.inner_nInteractionLengths = nInteractionLengths;
          caloLayer.inner_thickness = thickness_sum;
          if(layer_id==1)
            cout<<"Hcal_Barrel:  inner_thickness= "<<thickness_sum<<endl;
        //Store readout gasgap thickness
          caloLayer.sensitive_thickness = slice_thickness;
#endif
        //Reset counters to measure "outside" quantitites
          nRadiationLengths=0.;
          nInteractionLengths=0.;
          thickness_sum = 0.;
	}
         nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
         nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
         thickness_sum += slice_thickness/2;


	// Set region, limitset, and vis.
	slice_vol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	// slice PlacedVolume
	PlacedVolume slice_phv = ChamberLogical.placeVolume(slice_vol,Position(0.,0.,slice_pos_z));
        slice_phv.addPhysVolID("layer",logical_layer_id).addPhysVolID("slice",slice_number);
	if ( x_slice.isSensitive() ) {
	  int tower_id  = (layer_id > Hcal_nlayers)? 1:0;
 	  slice_phv.addPhysVolID("tower",tower_id);
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
      if(layer_id==1)
         cout<<"Hcal_Barrel:  outer_thickness= "<<thickness_sum<<endl;
#endif


      //---------------------------  Chamber Placements -----------------------------------------

      double chamber_x_offset, chamber_y_offset, chamber_z_offset;
      chamber_x_offset = xShift;
      chamber_z_offset = 0;

      chamber_y_offset = -(-Hcal_total_dim_y/2. 
			   + (logical_layer_id-1) * layer_thickness
			   + layer_thickness/2.);


      if(logical_layer_id==1 )
        cout<<"logical_layer_id, chamber_ x_offset, z_offset, y_offset :"<<logical_layer_id<<" "<<chamber_x_offset<<" "<<chamber_z_offset<<" "<<chamber_y_offset<<endl;
      pv =  EnvLogHcalModuleBarrel.placeVolume(ChamberLogical,
					       Position(chamber_x_offset,
							chamber_z_offset,
							chamber_y_offset ));
      //-----------------------------------------------------------------------------------------
      if( layer_id <= Hcal_nlayers) { // add the first set of layers to the reconstruction data
      
        caloLayer.distance = Hcal_inner_radius + Hcal_total_dim_y/2.0 + chamber_y_offset ;
        caloLayer.absorberThickness = Hcal_radiator_thickness ;
      
        caloData->layers.push_back( caloLayer ) ;
      }
      //-----------------------------------------------------------------------------------------

      
    }//end loop over HCAL nlayers;
  





  //====================================================================
  // Place HCAL Barrel stave module into the envelope
  //====================================================================
  double stave_phi_offset,  module_z_offset,  lateral_plate_z_offset;
  
  double Y = Hcal_inner_radius + Hcal_total_dim_y / 2.;

  stave_phi_offset = M_PI/Hcal_inner_symmetry;

  //-------- start loop over HCAL BARREL staves ----------------------------

  for (int stave_id = 1;
       stave_id <=Hcal_inner_symmetry;
       stave_id++)
    {
      lateral_plate_z_offset = - Hcal_module_wall_thickness/2.;

      double phirot = stave_phi_offset;
      RotationZYX rota(0,phirot,M_PI*0.5);
      Rotation3D rot3D(rota);

      for (int module_id = 1;
         module_id <=Hcal_barrel_number_modules;
         module_id++)
	{
          module_z_offset = -Hcal_half_length+ Hcal_normal_dim_z/2. + (module_id-1)*Hcal_normal_dim_z;
          if (stave_id ==1) {cout <<"\n module_id:    "<< module_id <<"   module_z_offset:    "<< module_z_offset
                                  <<" Hcal_normal_dim_z: "<< Hcal_normal_dim_z<< " Hcal_module_wall_thickness: "<<Hcal_module_wall_thickness<<endl;}

	  Position xyzVec(-Y*sin(phirot), Y*cos(phirot), module_z_offset);
	  Transform3D tran3D(rot3D,xyzVec);
	  
	  // Place Hcal Barrel volume into the envelope volume
	  pv = envelope.placeVolume(EnvLogHcalModuleBarrel,tran3D);
	  pv.addPhysVolID("stave",stave_id);
	  pv.addPhysVolID("module",module_id);
	  pv.addPhysVolID("system",det_id);

          lateral_plate_z_offset = - Hcal_half_length + Hcal_module_wall_thickness/2. + (module_id-1)*Hcal_normal_dim_z;
             
	  if (stave_id == 1) { //place only once for the whole lateral_plate
	    Position xyzVec_LP(0, 0,lateral_plate_z_offset);
	    pv = envelope.placeVolume(EnvLogHcalModuleBarrel_LP,xyzVec_LP);
//put the barrel wall also at the other end
            lateral_plate_z_offset = lateral_plate_z_offset + Hcal_regular_chamber_dim_z + Hcal_module_wall_thickness;
	    Position xyzVec_LP2(0, 0,lateral_plate_z_offset);
	    pv = envelope.placeVolume(EnvLogHcalModuleBarrel_LP,xyzVec_LP2);
//            cout <<" lateral_plate_z_offset LP2,  module_id, stave_id:   "<< lateral_plate_z_offset <<" , "<< module_id<<" , "<<stave_id <<endl;
	  }

	}
      stave_phi_offset +=  M_PI*2.0/Hcal_inner_symmetry;
    }  //-------- end loop over HCAL BARREL staves ----------------------------


  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;

 
  return sdet;

}

DECLARE_DETELEMENT(Hcal_Barrel_SD_v02, create_detector)

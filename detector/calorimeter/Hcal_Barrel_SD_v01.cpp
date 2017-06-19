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

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens) {
  
  double boundarySafety = 0.0001;

  xml_det_t   x_det       = element;
  string      det_name    = x_det.nameStr();
  DetElement  sdet( det_name,x_det.id() );

  
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
  cout << " \n\n\n CREATE DETECTOR: Hcal_Barrel_SD_v01" << endl;




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
//

  double Hcal_total_dim_y   = Hcal_nlayers*layer_thickness ;

  double Hcal_normal_dim_z  = 2 * Hcal_half_length/Hcal_barrel_number_modules;
  double Hcal_regular_chamber_dim_z = Hcal_normal_dim_z - 2.*Hcal_module_wall_thickness;
// each barrel (wheel) in a "box" with wall thickness of 10mm (in z-direction)
  
//  double Hcal_cell_dim_x            = cell_sizeX + Hcal_pad_separation;
//  double Hcal_cell_dim_z            = cell_sizeZ + Hcal_pad_separation;

  std::cout<< "Hcal_Barrel number of z-modules: "<<Hcal_barrel_number_modules<<std::endl;
  std::cout<< " ==> Hcal_cell_dim_x  (z): "<<cell_sizeX << " " << cell_sizeZ <<std::endl;
  std::cout <<"Hcal_half_length:    "<< Hcal_half_length<< " Hcal_modules_gap: "<<Hcal_modules_gap << endl;
  std::cout <<"Hcal_barrel_thickness:"<< Hcal_barrel_thickness << " Hcal_requested_thickness: "<<Hcal_total_dim_y << endl;
  std::cout <<"Hcal_regular_chamber_dim_z, Hcal_normal_dim_z:  "<< Hcal_regular_chamber_dim_z<<" "<<Hcal_normal_dim_z<<endl;
  std::cout <<" Hcal_MinNumCellsInTransvPlane:  "<< Hcal_MinNumCellsInTransvPlane<<endl;
 
  // ========= Create Hcal Barrel stave   ====================================
  //  It will be the volume for placing the Hcal Barrel Chamber(i.e. Layers).
  //  Itself will be placed into the world volume.
  // ==========================================================================
 
  Assembly EnvLogHcalModuleBarrel(det_name+"_module");

  EnvLogHcalModuleBarrel.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());



  double DHZ_LP  = Hcal_module_wall_thickness; 

  PolyhedraRegular stave_shaper_LP(Hcal_inner_symmetry, 
				   Hcal_inner_radius, Hcal_outer_radius, DHZ_LP);

  Tube solidCaloTube_LP(0, Hcal_outer_radius, DHZ_LP+boundarySafety);
  Cone solidCaloCone_LP( DHZ_LP+boundarySafety, 0, Hcal_outer_radius,0, Hcal_outer_radius);

  RotationZYX rot(0,0, M_PI/Hcal_inner_symmetry);

  IntersectionSolid Module_lateral_plate(stave_shaper_LP, solidCaloTube_LP,rot);
//  IntersectionSolid Module_lateral_plate(stave_shaper_LP, solidCaloCone_LP,rot);

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

  double xShift = 0.;

  LayeredCalorimeterData::Layer caloLayer ;
  caloLayer.cellSize0 = cell_sizeX;
  caloLayer.cellSize1 = cell_sizeZ;

  //-------------------- start loop over HCAL layers ----------------------
//  double layer_thickness =  Hcal_radiator_thickness +Hcal_chamber_thickness + 1.*Hcal_airgap_thickness;
  double AngleRatio=tan(M_PI/Hcal_inner_symmetry);//"k", updated for flexible symmetry
  double d_InnerOctoSize=AngleRatio*Hcal_inner_radius;//"d/2"
  cout <<" \n radiator_thickness, airgap_thickness  "<<Hcal_radiator_thickness<<" "<<Hcal_airgap_thickness <<endl ;
  cout <<" \n AngleRatio, d_InnerOctoSize  "<< AngleRatio<<" "<<d_InnerOctoSize <<endl ;
  cout <<" \n boundarySafety "<< boundarySafety <<endl ;

  for (int layer_id = 1; layer_id <= Hcal_nlayers; layer_id++)
    {
 
      double yn = sqrt(Hcal_outer_radius*Hcal_outer_radius -
		       (Hcal_inner_radius + layer_id*layer_thickness)*(Hcal_inner_radius + layer_id*layer_thickness));

      //updated for flexible symmetry
      double Yn = d_InnerOctoSize - layer_id*layer_thickness/tan(M_PI*2.0/Hcal_inner_symmetry);

      double halfX = layer_thickness/2.;
      double halfY = (yn+Yn)/2. - boundarySafety;
      double halfZ = Hcal_regular_chamber_dim_z / 2.;

      xShift = (yn-Yn)/2.;
      int  CellsInTransvPlane =  floor(2*halfY/cell_sizeX);

      // Check the space and number of possible readout cells for the next layer in the envelope 
      if( (Yn <= 0 && yn <= xShift) || (CellsInTransvPlane < Hcal_MinNumCellsInTransvPlane)) continue; 

      cout <<" \n layer_id; CellsInTransvPlane  "<<layer_id<<" "<<CellsInTransvPlane <<endl;
      cout <<" yn, xShift, Yn:   "<<yn<<" "<<xShift <<" "<<Yn<<endl;

      //--------------------------------------------------------------------------------
      //  build chamber box, with the calculated dimensions 
      //-------------------------------------------------------------------------------

      Box ChamberSolid(halfY,  //Add gap at two sides, do not need to build air gaps individualy.
		       halfZ,  //z attention!
		       halfX); //y attention!

      string ChamberLogical_name      = det_name+_toString(layer_id,"_layer%d");

      Volume ChamberLogical(ChamberLogical_name, ChamberSolid, air);   


      //====================================================================
      // Create Hcal Barrel Chamber with radiator
      // Place into the Hcal Barrel stave
      //====================================================================
      string layer_name      = det_name+_toString(layer_id,"_layer%d");

      double nRadiationLengths=0.;
      double nInteractionLengths=0.;
      double thickness_sum=0;
    
      nRadiationLengths   = Hcal_radiator_thickness/(stavesMaterial.radLength());
      nInteractionLengths = Hcal_radiator_thickness/(stavesMaterial.intLength());

      // Create the slices (sublayers) within the Hcal Barrel Chamber.
      double slice_pos_z = layer_thickness/2.;
      int slice_number = 0;
      for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
	xml_comp_t x_slice = k;
	string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
	double   slice_thickness = x_slice.thickness();
	Material slice_material  = theDetector.material(x_slice.materialStr());
        if(layer_id==1)
           cout<<"  Layer_slice:  "<<  slice_name<<" slice_thickness:  "<< slice_thickness<< endl;

	DetElement slice(layer_name,_toString(slice_number,"slice%d"),x_det.id());
	
	slice_pos_z -= slice_thickness/2.;
        nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
        nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
        thickness_sum       += slice_thickness/2;
	
	// Slice volume & box
	Volume slice_vol(slice_name,Box((halfY - Hcal_layer_air_gap),halfZ,slice_thickness/2.),slice_material);
	
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
	if ( x_slice.isSensitive() ) {
	  int slice_id  = (layer_id > Hcal_nlayers)? 1:-1;
	  slice_phv.addPhysVolID("layer",layer_id).addPhysVolID("slice",slice_id);
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
			   + (layer_id-1) * layer_thickness
			   + layer_thickness/2.);


      pv =  EnvLogHcalModuleBarrel.placeVolume(ChamberLogical,
					       Position(chamber_x_offset,
							chamber_z_offset,
							chamber_y_offset ));
      
      //-----------------------------------------------------------------------------------------
      
      caloLayer.distance = Hcal_inner_radius + Hcal_total_dim_y/2.0 + chamber_y_offset ;
      caloLayer.absorberThickness = Hcal_radiator_thickness ;
      
      caloData->layers.push_back( caloLayer ) ;
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

          lateral_plate_z_offset = - Hcal_half_length + Hcal_module_wall_thickness/2. + (module_id-1)*Hcal_normal_dim_z;
             
	  if (stave_id == 1) { //place only once for the whole lateral_plate
	    Position xyzVec_LP(0, 0,lateral_plate_z_offset);
	    pv = envelope.placeVolume(EnvLogHcalModuleBarrel_LP,xyzVec_LP);
//put the barrel wall also at the other end
            lateral_plate_z_offset = lateral_plate_z_offset + Hcal_regular_chamber_dim_z + Hcal_module_wall_thickness;
	    Position xyzVec_LP2(0, 0,lateral_plate_z_offset);
	    pv = envelope.placeVolume(EnvLogHcalModuleBarrel_LP,xyzVec_LP2);
//	    lateral_plate_z_offset = - lateral_plate_z_offset;
            cout <<" lateral_plate_z_offset LP2,  module_id, stave_id:   "<< lateral_plate_z_offset <<" , "<< module_id<<" , "<<stave_id <<endl;
	  }

	}
      stave_phi_offset +=  M_PI*2.0/Hcal_inner_symmetry;
    }  //-------- end loop over HCAL BARREL staves ----------------------------


  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;

 
  return sdet;

}

DECLARE_DETELEMENT(Hcal_Barrel_SD_v01, create_detector)

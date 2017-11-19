//====================================================================
//  AIDA Detector description implementation
//  for LDC AHCAL Endcap
//--------------------------------------------------------------------
//
//  Author     : S.Lu
//
// Basic idea:
// 1. Create the Hcal Endcap module envelope (16 modules).
//    Note: with default material Steel235.
//    
// 2. Create the Hcal Endcap Chamber(i.e. Layer) for each module.
//    Create the Layer with slices (Polystyrene,Cu,FR4,air).
//    Place each slice into the chamber with the right position,
//    And registry the IDs for slice
//
// 3. Place the same Layer into the endcap module envelope.
//    It will be repeated repeat 48 times.
//    And registry the IDs for layer, and towerID.
//
// 4. Place the endcap module into the world volume,
//    with the right position and rotation.
//    And registry the IDs for stave,module and towerID.
//
// 5. Customer material FR4 and Steel235 defined in materials.xml
//
// Tibor Kurca: Modify for SemiDigital Hcal_Endcaps &
//              - calculate chamber thickness from slice-thicknesses values in xml
//              - check if possible to build requested Nb of chambers
//              - the same geometry as SHcalSC04 (boxes), but divided into 4 staves
//                not corresponding to engineering design of SDHcal
//              - inner, outer symmetries should be also different ????!!!
//              - number of staves should be defined also by parameter in xml file????
//                                                                ²
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::RotationX;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {
  xml_det_t   x_det     = element;
  Layering    layering(x_det);

  Material    air         = theDetector.air();
  Material    stavesMaterial    = theDetector.material(x_det.materialStr());

  int           det_id    = x_det.id();
  string      det_name    = x_det.nameStr();


  DetElement   sdet(det_name,det_id);
//  Volume      motherVol = theDetector.pickMotherVolume(sdet); 

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  dd4hep::xml::setDetectorTypeFlag(element, sdet);

  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
  //-----------------------------------------------------------------------------------

  sens.setType("calorimeter");

  DetElement    stave_det("module0stave0",det_id);
 
  // The way to read constant from XML/THEDETECTOR file.
  double      Hcal_radiator_thickness          = theDetector.constant<double>("HcalSD_radiator_thickness");
  int         HcalEndcap_nlayers               = theDetector.constant<int>("HcalEndcapSD_nlayers");

  double      Hcal_stave_gaps                  = theDetector.constant<double>("HcalSD_stave_gaps");
  double      Hcal_lateral_plate_thickness     = theDetector.constant<double>("HcalSD_endcap_lateral_structure_thickness");
  double      Hcal_back_plate_thickness        = theDetector.constant<double>("HcalSD_back_plate_thickness");

  double      HcalEndcap_inner_radius          = theDetector.constant<double>("HcalEndcap_inner_radius");
  double      HcalEndcap_outer_radius          = theDetector.constant<double>("HcalEndcap_outer_radius");
  double      HcalEndcap_min_z                 = theDetector.constant<double>("HcalEndcap_min_z");
  double      HcalEndcap_max_z                 = theDetector.constant<double>("HcalEndcap_max_z");

  double      Hcal_endcap_layer_air_gap        = theDetector.constant<double>("Hcal_endcap_layer_air_gap");
// Some verbose output
  cout << " \n\n\n CREATE DETECTOR: Hcal_Endcaps_SD_v02" << endl;
  cout<<"    HCAL endcaps: HcalEndcap_inner_radius: "<< HcalEndcap_inner_radius <<endl;
  cout<<"    HCAL endcaps: HcalEndcap_outer_radius: "<< HcalEndcap_outer_radius <<endl;
  cout<<"    HCAL endcaps: HcalEndcap_min_z:  "<< HcalEndcap_min_z <<endl;
  cout<<"    HCAL endcaps: HcalEndcap_max_z:  "<< HcalEndcap_max_z <<endl;


  double SpaceForLayers = HcalEndcap_max_z -HcalEndcap_min_z
                        - Hcal_back_plate_thickness;

// First get the layer (chamber) thickness caluculated from the Hcal_Endcaps_SD_v02.xml file
// not from the fixed, defined value "Hcal_chamber_thickness"
  xml_coll_t c1(x_det,_U(layer));
  xml_comp_t   x1_layer = c1;

  double layer_thickness_z = 0.0;

  for(xml_coll_t k(x1_layer,_U(slice)); k; ++k)  {
     xml_comp_t x_slice = k;
     layer_thickness_z += x_slice.thickness();
  }
  cout<<"    layer_thickness (from slices) = "<<layer_thickness_z<<endl;
// chamber placements
   int number_of_chambers = HcalEndcap_nlayers;
   int possible_number_of_chambers = (int) floor( SpaceForLayers/ layer_thickness_z );
   if(possible_number_of_chambers < number_of_chambers)
        number_of_chambers = possible_number_of_chambers;

   cout<<"    Possible/ Requested #of chambers= "<<possible_number_of_chambers<<" / "<<HcalEndcap_nlayers<<endl;

//====================================================================
//
  
  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];

  cout<<"    HCAL endcaps: cell_sizeX, cell_sizeY:  "<<cell_sizeX <<" "<<cell_sizeY <<endl;
  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = 4  ; // hard code cernter box hole
  caloData->outer_symmetry = 0  ; // outer tube, or 8 for Octagon
  caloData->phi0 = 0 ;

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = HcalEndcap_inner_radius ;
  caloData->extent[1] = HcalEndcap_outer_radius ;
  caloData->extent[2] = HcalEndcap_min_z ;
  caloData->extent[3] = HcalEndcap_max_z ;
  

  int towerID = 0;
  for(xml_coll_t c(x_det.child(_U(dimensions)),_U(dimensions)); c; ++c) 
    {
      xml_comp_t l(c);
      
      double dim_x = l.attr<double>(_Unicode(dim_x));
      double dim_y = l.attr<double>(_Unicode(dim_y));
      double dim_z = l.attr<double>(_Unicode(dim_z));
      double pos_y = l.attr<double>(_Unicode(y_offset));
    
      // Hcal Endcap module shape
      double box_half_x= dim_x/2.0; // module width, all the same
      double box_half_y= dim_y/2.0; // module length, changing, 
      double box_half_z= dim_z/2.0; // total thickness, all the same
      
//      double x_offset = box_half_x*16-box_half_x*towerID*2.0-box_half_x;
      double x_offset = dim_x*(0.5 + towerID);
//      double y_offset = pos_y + Hcal_lateral_plate_thickness/2.0;
      double y_offset = pos_y + Hcal_stave_gaps/2.0;
  
      cout << "Hcal_Endcaps: towerID, dim_x, dim_y, dim_z = "<<towerID<<" "<<dim_x<<" "<<dim_y<<" "<<dim_z <<endl;
      cout << "Offsets: x_offset, y_offset =" << x_offset<<" "<<y_offset<<endl;
   
      Box    EndcapModule(box_half_x,box_half_y,box_half_z);
      
      // define the name of each endcap Module
      string envelopeVol_name   = det_name+_toString(towerID,"_EndcapModule%d");
      
      Volume envelopeVol(envelopeVol_name,EndcapModule,stavesMaterial);
      
      // Set envelope volume attributes.
      envelopeVol.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
      
      
      // ========= Create Hcal Chamber (i.e. Layers) ==============================
      // It will be the sub volume for placing the slices.
      // Itself will be placed into the Hcal Endcap modules envelope.
      // ==========================================================================
      
      // create Layer (air) and place the slices into it. 
      // place the Layer into the Hcal Endcap Modules envelope (stavesMaterial).
      
      double layer_pos_z     = - box_half_z ;                      
      
      // Create Hcal Endcap Chamber without radiator
      // Place into the Hcal Encap module envelope, after each radiator 
      int layer_num = 1;
      for(xml_coll_t k1(x_det,_U(layer)); k1; ++k1)  {
	xml_comp_t   x_layer = k1;
        //	int          repeat = x_layer.repeat();          // Get number of layers.
	//const Layer* lay    = layering.layer(layer_num); // Get the layer from the layering engine.
	//	double layer_thickness = lay->thickness();

	double layer_thickness = layering.layer(layer_num)->thickness();
	string layer_name      = envelopeVol_name+"_layer";
	DetElement  layer(stave_det,layer_name,det_id);
	
	// Active Layer box & volume
	double active_layer_dim_x = box_half_x - Hcal_lateral_plate_thickness - Hcal_endcap_layer_air_gap;
	double active_layer_dim_z = layer_thickness/2.0;
	double active_layer_dim_y = box_half_y;// - Hcal_lateral_platee_thickness - Hcal_endcap_layer_air_gap;
	
	// Build chamber including air gap
	// The Layer will be filled with slices, 
	Volume layer_vol(layer_name, Box((active_layer_dim_x + Hcal_endcap_layer_air_gap),
					 active_layer_dim_y,active_layer_dim_z), air);

	LayeredCalorimeterData::Layer caloLayer ;
	caloLayer.cellSize0 = cell_sizeX;
	caloLayer.cellSize1 = cell_sizeY;
	
	// ========= Create sublayer slices =========================================
	// Create and place the slices into Layer
	// ==========================================================================
	
	// Create the slices (sublayers) within the Hcal Chamber.
	double slice_pos_z = -(layer_thickness / 2.0);
	int slice_number = 0;

	double nRadiationLengths=0.;
	double nInteractionLengths=0.;
	double thickness_sum=0;

	nRadiationLengths   = Hcal_radiator_thickness/(stavesMaterial.radLength());
	nInteractionLengths = Hcal_radiator_thickness/(stavesMaterial.intLength());

	for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
	  xml_comp_t x_slice = k;
	  string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
	  double   slice_thickness = x_slice.thickness();
	  Material slice_material  = theDetector.material(x_slice.materialStr());
	  DetElement slice(layer,_toString(slice_number,"slice%d"),det_id);
	  
	  slice_pos_z += slice_thickness / 2.0;
	  
	  // Slice volume & box
	  Volume slice_vol(slice_name,Box(active_layer_dim_x,active_layer_dim_y,slice_thickness/2.0),slice_material);
	  
	  nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
	  nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
	  thickness_sum       += slice_thickness/2;

	  if ( x_slice.isSensitive() ) {
	    sens.setType("calorimeter");
	    slice_vol.setSensitiveDetector(sens);

#if DD4HEP_VERSION_GE( 0, 15 )
	    //Store "inner" quantities
	    caloLayer.inner_nRadiationLengths = nRadiationLengths;
	    caloLayer.inner_nInteractionLengths = nInteractionLengths;
	    caloLayer.inner_thickness = thickness_sum;
	    //Store sensitive slice (gas) thickness
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
	  PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol,Position(0,0,slice_pos_z));
	  //slice_phv.addPhysVolID("slice",slice_number);
	  
	  slice.setPlacement(slice_phv);
	  // Increment Z position for next slice.
	  slice_pos_z += slice_thickness / 2.0;
	  // Increment slice number.
	  ++slice_number;             
	}
	// Set region, limitset, and vis.
	layer_vol.setAttributes(theDetector,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());

#if DD4HEP_VERSION_GE( 0, 15 )
	//Store "outer" quantities
	caloLayer.outer_nRadiationLengths = nRadiationLengths;
	caloLayer.outer_nInteractionLengths = nInteractionLengths;
	caloLayer.outer_thickness = thickness_sum;
#endif 	
	
	// ========= Place the Layer (i.e. Chamber) =================================
	// Place the Layer into the Hcal Endcap module envelope.
	// with the right position and rotation.
	// Registry the IDs (layer, stave, module).
	// Place the same layer 48 times into Endcap module
	// ==========================================================================
	
	for (int j = 0; j < number_of_chambers; j++)    {
	  
	  // Layer position in z within the Endcap Modules.
	  layer_pos_z += layer_thickness / 2.0;
	  
	  PlacedVolume layer_phv = envelopeVol.placeVolume(layer_vol,
							   Position(0,0,layer_pos_z));
	  // registry the ID of Layer, stave and module
	  layer_phv.addPhysVolID("layer",layer_num);

	  // then setPlacement for it.
	  layer.setPlacement(layer_phv);
	  
	  
	  //-----------------------------------------------------------------------------------------
	  if ( caloData->layers.size() < (unsigned int)number_of_chambers ) {

	    caloLayer.distance = HcalEndcap_min_z + box_half_z + layer_pos_z
	      - caloLayer.inner_thickness ; // Will be added later at "DDMarlinPandora/DDGeometryCreator.cc:226" to get center of sensitive element
	    caloLayer.absorberThickness = Hcal_radiator_thickness ;
	    
	    caloData->layers.push_back( caloLayer ) ;
	  }
	  //-----------------------------------------------------------------------------------------
	  
	  
	  // ===== Prepare for next layer (i.e. next Chamber) =========================
	  // Prepare the chamber placement position and the chamber dimension
	  // ==========================================================================
	  
	  // Increment the layer_pos_y
	  // Place Hcal Chamber after each radiator 
	  layer_pos_z += layer_thickness / 2.0;
	  ++layer_num;         
	}
	
	
      }
      
      
      // =========== Place Hcal Endcap envelope ===================================
      // Finally place the Hcal Endcap envelope into the world volume.
      // Registry the stave(up/down), module(left/right) and towerID.
      // ==========================================================================
      
      // Acording to the number of staves and modules,
      // Place the same Hcal Endcap module volume into the world volume
      // with the right position and rotation.
      for(int stave_num=0;stave_num<4;stave_num++){
	
	double EndcapModule_pos_x = 0;
	double EndcapModule_pos_y = 0;
	double EndcapModule_pos_z = 0;
	double rot_EM = 0;
//	double rot_EMZ = 0;

	double EndcapModule_center_pos_z = HcalEndcap_min_z + box_half_z;
          if(stave_num==0){ 
            EndcapModule_pos_x = -x_offset;
	    EndcapModule_pos_y = y_offset;
          }
          if(stave_num==1){ 
            EndcapModule_pos_x = x_offset;
	    EndcapModule_pos_y = y_offset;
          }
          if(stave_num==2){ 
            EndcapModule_pos_x = x_offset;
	    EndcapModule_pos_y = -y_offset;
          }
          if(stave_num==3){ 
            EndcapModule_pos_x = -x_offset;
	    EndcapModule_pos_y = -y_offset;
          }
	
	
	for(int module_num=0;module_num<2;module_num++) {

	  int module_id = (module_num==0)? 0:6;
	  
//	  rot_EMZ = (module_id==0)?-M_PI:0.0;
	  rot_EM = (module_id==0)?M_PI:0.0;
//          cout << " Hcal_Endcaps: module_id = " << module_id << " rot_EM= "<<rot_EM  << " rot_EMZ= "<<rot_EMZ <<endl;	  
	  EndcapModule_pos_z = (module_id==0)? -EndcapModule_center_pos_z:EndcapModule_center_pos_z;
          cout << " Hcal_Endcaps: module,towerID,staveID, module_pos_x,y,z = "<<module_id<<" "<<towerID<<" "<<stave_num<<" " <<EndcapModule_pos_x <<" "<<EndcapModule_pos_y<<" "<<EndcapModule_pos_z <<endl;	  

          RotationZYX rotECM(0.0,rot_EM,0.0);
//						      Transform3D(rotECM,
	  PlacedVolume env_phv = envelope.placeVolume(envelopeVol,
                                                      Transform3D(RotationX(rot_EM),
								  Translation3D(EndcapModule_pos_x,
										EndcapModule_pos_y,
										EndcapModule_pos_z)));
	  env_phv.addPhysVolID("tower",towerID);	  
	  env_phv.addPhysVolID("stave",stave_num);   // y: up /down
	  env_phv.addPhysVolID("module",module_id); // z: -/+ 0/6
	  env_phv.addPhysVolID("system",det_id);

	  DetElement sd = (module_num==0&&stave_num==0) ? stave_det : stave_det.clone(_toString(module_id,"module%d")+_toString(stave_num,"stave%d"));	  
	  sd.setPlacement(env_phv);	  

	}
	
      }

    towerID++;
      
    }
  
  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;  
  
  return sdet;
}




DECLARE_DETELEMENT(Hcal_Endcaps_SD_v02, create_detector)

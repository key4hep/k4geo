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
//    And registry the IDs for layer, and endcapID.
//
// 4. Place the endcap module into the world volume,
//    with the right position and rotation.
//    And registry the IDs for stave,module and endcapID.
//
// 5. Customer material FR4 and Steel235 defined in materials.xml
//
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
using dd4hep::DetType;
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::RotationX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Tube;
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
  xml_dim_t   dim         = x_det.dimensions();
  string      det_name    = x_det.nameStr();
  //unused: string      det_type    = x_det.typeStr();
  Material    air         = theDetector.air();
  Material    stavesMaterial    = theDetector.material(x_det.materialStr());
  int         numSides    = dim.numsides();

  int           det_id    = x_det.id();

  DetElement   sdet(det_name,det_id);

  PlacedVolume pVol;

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  sdet.setTypeFlag( DetType::CALORIMETER |  DetType::ENDCAP  | DetType::HADRONIC ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
  //-----------------------------------------------------------------------------------

  sens.setType("calorimeter");

  DetElement    stave_det("module0stave0",det_id);
 
  // The way to reaad constant from XML/Detector file.
  double      Hcal_radiator_thickness          = theDetector.constant<double>("Hcal_radiator_thickness");
  double      Hcal_endcap_lateral_structure_thickness = theDetector.constant<double>("Hcal_endcap_lateral_structure_thickness");
  double      Hcal_endcap_layer_air_gap        = theDetector.constant<double>("Hcal_endcap_layer_air_gap");

  //double      Hcal_cells_size                  = theDetector.constant<double>("Hcal_cells_size");
  double      HcalEndcap_inner_radius          = theDetector.constant<double>("HcalEndcap_inner_radius");
  double      HcalEndcap_outer_radius          = theDetector.constant<double>("HcalEndcap_outer_radius");
  double      HcalEndcap_min_z                 = theDetector.constant<double>("HcalEndcap_min_z");
  double      HcalEndcap_max_z                 = theDetector.constant<double>("HcalEndcap_max_z");

  double   Hcal_steel_cassette_thickness       = theDetector.constant<double>("Hcal_steel_cassette_thickness");
  double   HcalServices_outer_FR4_thickness    = theDetector.constant<double>("HcalServices_outer_FR4_thickness");
  double   HcalServices_outer_Cu_thickness     = theDetector.constant<double>("HcalServices_outer_Cu_thickness");
  double   Hcal_endcap_services_module_width   = theDetector.constant<double>("Hcal_endcap_services_module_width");

  Material  stainless_steel =  theDetector.material("stainless_steel");
  Material  PCB             =  theDetector.material("PCB");
  Material  copper          =  theDetector.material("Cu");

  std::cout <<"\n HcalEndcap_inner_radius = "
	    <<HcalEndcap_inner_radius/dd4hep::mm <<" mm"
	    <<"\n HcalEndcap_outer_radius = "
	    <<HcalEndcap_outer_radius/dd4hep::mm <<" mm"
	    <<"\n HcalEndcap_min_z = "
	    <<HcalEndcap_min_z/dd4hep::mm <<" mm"
	    <<"\n HcalEndcap_max_z = "
	    <<HcalEndcap_max_z/dd4hep::mm <<" mm"
	    <<std::endl;


  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];

  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = 4  ; // hard code cernter box hole
  caloData->outer_symmetry = 0  ; // outer tube, or 8 for Octagun
  caloData->phi0 = 0 ;

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = HcalEndcap_inner_radius ;
  caloData->extent[1] = HcalEndcap_outer_radius ;
  caloData->extent[2] = HcalEndcap_min_z ;
  caloData->extent[3] = HcalEndcap_max_z ;
  

  int endcapID = 0;
  for(xml_coll_t c(x_det.child(_U(dimensions)),_U(dimensions)); c; ++c) 
    {
      xml_comp_t l(c);
      
      double dim_x = l.attr<double>(_Unicode(dim_x));
      double dim_y = l.attr<double>(_Unicode(dim_y));
      double dim_z = l.attr<double>(_Unicode(dim_z));
      double pos_y = l.attr<double>(_Unicode(y_offset));
    
      // Hcal Endcap module shape
      double box_half_x= dim_x/2.0; // module width, all are same
      double box_half_y= dim_y/2.0; // total thickness, all are same
      double box_half_z= dim_z/2.0; // module length, changed, 
      
      double x_offset = box_half_x*numSides-box_half_x*endcapID*2.0-box_half_x;
      double y_offset = pos_y;
      
      Box    EndcapModule(box_half_x,box_half_y,box_half_z);
      
      // define the name of each endcap Module
      string envelopeVol_name   = det_name+_toString(endcapID,"_EndcapModule%d");
      
      Volume envelopeVol(envelopeVol_name,EndcapModule,stavesMaterial);
      
      // Set envelope volume attributes.
      envelopeVol.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
      

      double FEE_half_x = box_half_x-Hcal_endcap_services_module_width/2.0;
      double FEE_half_y = box_half_y;
      double FEE_half_Z = Hcal_endcap_services_module_width/2.0;

      Box    FEEBox(FEE_half_x,FEE_half_y,FEE_half_Z);
      Volume FEEModule("Hcal_endcap_FEE",FEEBox,air);

      double FEELayer_thickness = Hcal_steel_cassette_thickness + HcalServices_outer_FR4_thickness + HcalServices_outer_Cu_thickness;
      Box    FEELayerBox(FEE_half_x,FEELayer_thickness/2.0,FEE_half_Z);
      Volume FEELayer("FEELayer",FEELayerBox,air);

      Box    FEELayerSteelBox(FEE_half_x,Hcal_steel_cassette_thickness/2.0,FEE_half_Z);
      Volume FEELayerSteel("FEELayerSteel",FEELayerSteelBox,stainless_steel);
      pVol = FEELayer.placeVolume(FEELayerSteel,
				  Position(0,
					   (-FEELayer_thickness/2.0
					    +Hcal_steel_cassette_thickness/2.0),
					   0));

      Box    FEELayerFR4Box(FEE_half_x,HcalServices_outer_FR4_thickness/2.0,FEE_half_Z);
      Volume FEELayerFR4("FEELayerFR4",FEELayerFR4Box,PCB);
      pVol = FEELayer.placeVolume(FEELayerFR4,
				  Position(0,
					   (-FEELayer_thickness/2.0+Hcal_steel_cassette_thickness
					    +HcalServices_outer_FR4_thickness/2.0),
					   0));

      Box    FEELayerCuBox(FEE_half_x,HcalServices_outer_Cu_thickness/2.0,FEE_half_Z);
      Volume FEELayerCu("FEELayerCu",FEELayerCuBox,copper);
      pVol = FEELayer.placeVolume(FEELayerCu,
				  Position(0,
					   (-FEELayer_thickness/2.0+Hcal_steel_cassette_thickness+HcalServices_outer_FR4_thickness
					    +HcalServices_outer_Cu_thickness/2.0),
					   0));


      // ========= Create Hcal Chamber (i.e. Layers) ==============================
      // It will be the sub volume for placing the slices.
      // Itself will be placed into the Hcal Endcap modules envelope.
      // ==========================================================================
      
      // create Layer (air) and place the slices (Polystyrene,Cu,FR4,air) into it. 
      // place the Layer into the Hcal Endcap Modules envelope (stavesMaterial).
      
      // First Hcal Chamber position, start after first radiator
      double layer_pos_y     = - box_half_y + Hcal_radiator_thickness;                      
      
      // Create Hcal Endcap Chamber without radiator
      // Place into the Hcal Encap module envelope, after each radiator 
      int layer_num = 1;
      for(xml_coll_t m(x_det,_U(layer)); m; ++m)  {
	xml_comp_t   x_layer = m;
	int          repeat = x_layer.repeat();          // Get number of layers.

	double layer_thickness = layering.layer(layer_num)->thickness();
	string layer_name      = envelopeVol_name+"_layer";
	DetElement  layer(stave_det,layer_name,det_id);
	
	// Active Layer box & volume
	double active_layer_dim_x = box_half_x - Hcal_endcap_lateral_structure_thickness - Hcal_endcap_layer_air_gap;
	double active_layer_dim_y = layer_thickness/2.0;
	double active_layer_dim_z = box_half_z;
	
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
	double slice_pos_y = -(layer_thickness / 2.0);
	int slice_number = 0;

	double nRadiationLengths=0.;
	double nInteractionLengths=0.;
	double thickness_sum=0;

	nRadiationLengths   = Hcal_radiator_thickness/(stavesMaterial.radLength());
	nInteractionLengths = Hcal_radiator_thickness/(stavesMaterial.intLength());
	thickness_sum       = Hcal_radiator_thickness;

	for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
	  xml_comp_t x_slice = k;
	  string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
	  double   slice_thickness = x_slice.thickness();
	  Material slice_material  = theDetector.material(x_slice.materialStr());
	  DetElement slice(layer,_toString(slice_number,"slice%d"),det_id);
	  
	  slice_pos_y += slice_thickness / 2.0;
	  
	  // Slice volume & box
	  Volume slice_vol(slice_name,Box(active_layer_dim_x,slice_thickness/2.0,active_layer_dim_z),slice_material);
	  
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
	    //Store scintillator thickness
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
	  PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol,Position(0,slice_pos_y,0));
	  //slice_phv.addPhysVolID("slice",slice_number);
	  
	  slice.setPlacement(slice_phv);
	  // Increment Z position for next slice.
	  slice_pos_y += slice_thickness / 2.0;
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
	
	for (int j = 0; j < repeat; j++)    {
	  
	  // Layer position in y within the Endcap Modules.
	  layer_pos_y += layer_thickness / 2.0;
	  
	  PlacedVolume layer_phv = envelopeVol.placeVolume(layer_vol,
							   Position(0,layer_pos_y,0));
	  // registry the ID of Layer, stave and module
	  layer_phv.addPhysVolID("layer",layer_num);

	  // then setPlacement for it.
	  layer.setPlacement(layer_phv);
	  
	  pVol = FEEModule.placeVolume(FEELayer,
				       Position(0,layer_pos_y,0));
	  //-----------------------------------------------------------------------------------------
	  if ( caloData->layers.size() < (unsigned int)repeat ) {

	    caloLayer.distance = HcalEndcap_min_z + box_half_y + layer_pos_y
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
	  layer_pos_y += layer_thickness / 2.0;
	  layer_pos_y += Hcal_radiator_thickness;
	  ++layer_num;         
	}
	
	
      }
      
      
      // =========== Place Hcal Endcap envelope ===================================
      // Finally place the Hcal Endcap envelope into the world volume.
      // Registry the stave(up/down), module(left/right) and endcapID.
      // ==========================================================================
      
      // Acording to the number of staves and modules,
      // Place the same Hcal Endcap module volume into the world volume
      // with the right position and rotation.
      for(int stave_num=0;stave_num<2;stave_num++){
	
	double EndcapModule_pos_x = 0;
	double EndcapModule_pos_y = 0;
	double EndcapModule_pos_z = 0;
	double rot_EM = 0;

	double EndcapModule_center_pos_z = HcalEndcap_min_z + box_half_y;

	double FEEModule_pos_x = 0;
	double FEEModule_pos_y = 0;
	double FEEModule_pos_z = 0;
	double FEEModule_center_pos_z = HcalEndcap_min_z + box_half_y;
	
	switch (stave_num)
	  {
	  case 0 : 
	    EndcapModule_pos_x = x_offset;
	    EndcapModule_pos_y = y_offset;
	    FEEModule_pos_x = x_offset;
	    FEEModule_pos_y = y_offset + box_half_z + Hcal_endcap_services_module_width/2.0;
	    break;
	  case 1 : 
	    EndcapModule_pos_x = -x_offset;
	    EndcapModule_pos_y = -y_offset;
	    FEEModule_pos_x = -x_offset;
	    FEEModule_pos_y = -y_offset - box_half_z - Hcal_endcap_services_module_width/2.0;
	    break;
	  }
	
	for(int module_num=0;module_num<2;module_num++) {

	  int module_id = (module_num==0)? 0:6;
	  
	  rot_EM = (module_id==0)?(-M_PI/2.0):(M_PI/2.0);
	  
	  EndcapModule_pos_z = (module_id==0)? -EndcapModule_center_pos_z:EndcapModule_center_pos_z;

	  PlacedVolume env_phv = envelope.placeVolume(envelopeVol,
						      Transform3D(RotationX(rot_EM),
								  Translation3D(EndcapModule_pos_x,
										EndcapModule_pos_y,
										EndcapModule_pos_z)));
	  env_phv.addPhysVolID("tower",endcapID);	  
	  env_phv.addPhysVolID("stave",stave_num);   // y: up /down
	  env_phv.addPhysVolID("module",module_id); // z: -/+ 0/6
	  env_phv.addPhysVolID("system",det_id);


	  FEEModule_pos_z = (module_id==0)? -FEEModule_center_pos_z:FEEModule_center_pos_z;

	  if (!(endcapID==0))
	    env_phv = envelope.placeVolume(FEEModule,
					   Transform3D(RotationX(rot_EM),
						       Translation3D(FEEModule_pos_x,
								     FEEModule_pos_y,
								     FEEModule_pos_z)));


	  DetElement sd = (module_num==0&&stave_num==0) ? stave_det : stave_det.clone(_toString(module_id,"module%d")+_toString(stave_num,"stave%d"));	  
	  sd.setPlacement(env_phv);	  

	}
	
      }

    endcapID++;
      
    }
  
  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;  
  
  return sdet;
}




DECLARE_DETELEMENT(SHcalSc04_Endcaps, create_detector)

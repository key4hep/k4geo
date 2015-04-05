//====================================================================
//  DDSim - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for SiWEcalEndcaps
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id$
//====================================================================

 /* History:  
  
// *******************************************************
// *                                                     *
// *                      Mokka                          * 
// *   - the detailed geant4 simulation for ILC   -      *
// *                                                     *
// * For more information about Mokka, visit the         *
// *                                                     *
// *  Mokka.in2p3.fr  Mokka home page.                   *
// *                                                     *
// *******************************************************
//
// $Id$
// $Name: mokka-07-00 $
//
// 
//
// SEcal04.cc
//
   Shaojun Lu:  Ported from Mokka SEcal04 Endcaps part. Read the constants from XML
                instead of the DB. Then build the Endcap in the same way with DD4hep
		construct.
		Inside SEcal04, some parameters, which used by Ecal Endcaps, come from
		Ecal Barrel. They can be seen here again.
		Start ECRing ...
 */

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

//#define VERBOSE 1

static Ref_t create_detector(LCDD& lcdd, xml_h element, SensitiveDetector sens)  {
  static double tolerance = 0e0;

  xml_det_t     x_det     = element;
  string        det_name  = x_det.nameStr();
  Layering      layering (element);

  Material      air       = lcdd.air();
  Material      vacuum    = lcdd.vacuum();

  int           det_id    = x_det.id();
  xml_comp_t    x_staves  = x_det.staves();
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = XML::createPlacedEnvelope( lcdd,  element , sdet ) ;
  
  if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  sens.setType("calorimeter");

  DetElement    module_det("module0",det_id);


//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build Ecal ECRing
//
//====================================================================

  int N_FIBERS_W_STRUCTURE = 2; 
  int N_FIBERS_ALVOULUS = 3;

  //  read parametere from compact.xml file
  double Ecal_Alveolus_Air_Gap              = lcdd.constant<double>("Ecal_Alveolus_Air_Gap");
  double Ecal_Slab_shielding                = lcdd.constant<double>("Ecal_Slab_shielding");
  double Ecal_Slab_copper_thickness         = lcdd.constant<double>("Ecal_Slab_copper_thickness");
  double Ecal_Slab_PCB_thickness            = lcdd.constant<double>("Ecal_Slab_PCB_thickness");
  double Ecal_Slab_glue_gap                 = lcdd.constant<double>("Ecal_Slab_glue_gap");
  double Ecal_Slab_ground_thickness         = lcdd.constant<double>("Ecal_Slab_ground_thickness");
  double Ecal_fiber_thickness               = lcdd.constant<double>("Ecal_fiber_thickness");
  double Ecal_Si_thickness                  = lcdd.constant<double>("Ecal_Si_thickness");
  
  double Ecal_radiator_thickness1           = lcdd.constant<double>("Ecal_radiator_layers_set1_thickness");
  double Ecal_radiator_thickness2           = lcdd.constant<double>("Ecal_radiator_layers_set2_thickness");
  double Ecal_radiator_thickness3           = lcdd.constant<double>("Ecal_radiator_layers_set3_thickness");
  double Ecal_Barrel_halfZ                  = lcdd.constant<double>("Ecal_Barrel_halfZ");
  
  double Ecal_support_thickness             = lcdd.constant<double>("Ecal_support_thickness");
  double Ecal_front_face_thickness          = lcdd.constant<double>("Ecal_front_face_thickness");
  double Ecal_lateral_face_thickness        = lcdd.constant<double>("Ecal_lateral_face_thickness");
 
  double Ecal_Slab_Sc_PCB_thickness         = lcdd.constant<double>("Ecal_Slab_Sc_PCB_thickness");
  double Ecal_Sc_thickness                  = lcdd.constant<double>("Ecal_Sc_thickness");
  double Ecal_Sc_reflector_thickness        = lcdd.constant<double>("Ecal_Sc_reflector_thickness");

  double Ecal_EC_Ring_gap                   = lcdd.constant<double>("Ecal_EC_Ring_gap");

  double Ecal_cables_gap                    = lcdd.constant<double>("Ecal_cables_gap");
  double Lcal_outer_radius                  = lcdd.constant<double>("Lcal_outer_radius");
  double Ecal_Lcal_ring_gap                 = lcdd.constant<double>("Ecal_Lcal_ring_gap");
  double Ecal_endcap_center_box_size        = lcdd.constant<double>("Ecal_endcap_center_box_size");

  int    Ecal_nlayers1                      = lcdd.constant<int>("Ecal_nlayers1");
  int    Ecal_nlayers2                      = lcdd.constant<int>("Ecal_nlayers2");
  int    Ecal_nlayers3                      = lcdd.constant<int>("Ecal_nlayers3");
  

  double      Ecal_cells_size               = lcdd.constant<double>("Ecal_cells_size");
  double      EcalEndcapRing_inner_radius   = lcdd.constant<double>("EcalEndcapRing_inner_radius");
  double      EcalEndcapRing_outer_radius   = lcdd.constant<double>("EcalEndcapRing_outer_radius");
  double      EcalEndcapRing_min_z          = lcdd.constant<double>("EcalEndcapRing_min_z");
  double      EcalEndcapRing_max_z          = lcdd.constant<double>("EcalEndcapRing_max_z");

  //========== fill data for reconstruction ============================
  DDRec::LayeredCalorimeterData* caloData = new DDRec::LayeredCalorimeterData ;
  caloData->layoutType = DDRec::LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = 0  ; // hard code cernter pipe hole
  caloData->outer_symmetry = 4  ; // outer box
  caloData->phi0 = 0 ;

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = EcalEndcapRing_inner_radius ;
  caloData->extent[1] = EcalEndcapRing_outer_radius ;
  caloData->extent[2] = EcalEndcapRing_min_z ;
  caloData->extent[3] = EcalEndcapRing_max_z ;


//====================================================================
//
// general calculated parameters
//
//====================================================================
  
  double Ecal_total_SiSlab_thickness = 
    Ecal_Slab_shielding + 
    Ecal_Slab_copper_thickness + 
    Ecal_Slab_PCB_thickness +
    Ecal_Slab_glue_gap + 
    Ecal_Si_thickness + 
    Ecal_Slab_ground_thickness +
    Ecal_Alveolus_Air_Gap / 2;
#ifdef VERBOSE
  std::cout << " Ecal_total_SiSlab_thickness = " << Ecal_total_SiSlab_thickness  << std::endl;
#endif
  
  

  double Ecal_total_ScSlab_thickness = 
    Ecal_Slab_shielding + 
    Ecal_Slab_copper_thickness + 
    Ecal_Slab_Sc_PCB_thickness +
    Ecal_Sc_thickness + 
    Ecal_Sc_reflector_thickness * 2 +
    Ecal_Alveolus_Air_Gap / 2;
#ifdef VERBOSE
  std::cout << " Ecal_total_ScSlab_thickness = " << Ecal_total_ScSlab_thickness  << std::endl;
#endif
  

    int Number_of_Si_Layers_in_Barrel = 0;
    int Number_of_Sc_Layers_in_Barrel = 0;


#ifdef VERBOSE
  std::cout << " Ecal total number of Silicon layers = " << Number_of_Si_Layers_in_Barrel  << std::endl;
  std::cout << " Ecal total number of Scintillator layers = " << Number_of_Sc_Layers_in_Barrel  << std::endl;
#endif
  
  // In this release the number of modules is fixed to 5
  double Ecal_Barrel_module_dim_z = 2 * Ecal_Barrel_halfZ / 5. ;
#ifdef VERBOSE
  std::cout << "Ecal_Barrel_module_dim_z  = " << Ecal_Barrel_module_dim_z  << std::endl;
#endif

  int n_total_layers = Ecal_nlayers1 + Ecal_nlayers2 + Ecal_nlayers3;
  Number_of_Si_Layers_in_Barrel = n_total_layers+1;

  //int total_number_of_layers = Ecal_nlayers1 + Ecal_nlayers2 + Ecal_nlayers3;

  double module_thickness = 
    Ecal_nlayers1 * Ecal_radiator_thickness1 +
    Ecal_nlayers2 * Ecal_radiator_thickness2 +
    Ecal_nlayers3 * Ecal_radiator_thickness3 +
    
    int(n_total_layers/2) * // fiber around W struct layers
    (N_FIBERS_W_STRUCTURE * 2 *  Ecal_fiber_thickness) +
    
    Number_of_Si_Layers_in_Barrel * // Silicon slabs plus fiber around and inside
    (Ecal_total_SiSlab_thickness +
     (N_FIBERS_ALVOULUS + 1 ) * Ecal_fiber_thickness) +
    
    Number_of_Sc_Layers_in_Barrel * // Scintillator slabs plus fiber around and inside
    (Ecal_total_ScSlab_thickness +
     (N_FIBERS_ALVOULUS + 1 ) * Ecal_fiber_thickness) +
    
    Ecal_support_thickness + Ecal_front_face_thickness;
  
#ifdef VERBOSE
  std::cout << " module_thickness = " << module_thickness  << std::endl;
#endif
  

  double ECRingSiplateSize = Ecal_endcap_center_box_size 
    - 2 * Ecal_EC_Ring_gap
    - 2 * Ecal_lateral_face_thickness;


// ========= Create Ecal end cap ring   ====================================
//  It will be the volume for palcing the Ecal endcaps alveolus(i.e. Layers).
//  And the structure W plate.
//  Itself will be placed into the world volume.
// ==========================================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~              EndcapRing                           ~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  double Ecal_endcap_Tube_rmax = Lcal_outer_radius + Ecal_Lcal_ring_gap;
  
  Tube CenterECTub(0., Ecal_endcap_Tube_rmax, module_thickness);
  Box  CenterECBox((Ecal_endcap_center_box_size/2. - Ecal_EC_Ring_gap),
		   (Ecal_endcap_center_box_size/2. - Ecal_EC_Ring_gap),
		   module_thickness/2.);

  SubtractionSolid ECRingSolid( CenterECBox, CenterECTub);

  Volume EnvLogECRing("ECRing",ECRingSolid,air);
 


  //==============================================================
  // build the layer and place into the ECRing
  //==============================================================


  //-------------------------------------------------------
  // Radiator and towers placements inside the ECRing module
  //-------------------------------------------------------

  // We count the layers starting from IP and from 1,
  // so odd layers should be inside slabs and
  // even ones on the structure.
  //
  double z_floor = 
    - module_thickness/2 +
    Ecal_front_face_thickness + 
    N_FIBERS_ALVOULUS * Ecal_fiber_thickness;

  double l_pos_z = z_floor;
 
  //-------------------- start loop over ECAL layers ----------------------
  // Loop over the sets of layer elements in the detector.
  int l_num = 1;
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    xml_comp_t x_layer = li;
    int repeat = x_layer.repeat();
    // Loop over number of repeats for this layer.
    for (int j=0; j<repeat; j++)    {
      
      string l_name = _toString(l_num,"layer%d");
      double l_thickness = layering.layer(l_num-1)->thickness();  // Layer's thickness.
      l_pos_z  += l_thickness/2.;
      
      int EC_Number_of_towers = 0; 

      double radiator_dim_y = -1.0; //to be updated with slice radiator thickness 
      

	  
      // We use the same method able to create the Barrel
      // Slabs, so we have to rotate it later when placing 
      // into the EC modules.
      //
      // We use the same method able to create the Barrel
      // radiator plates between slabs, but with the good
      // dimensions to avoid to rotate it later when placing 
      // into the EC modules.
      
      // While the towers have the same shape use the same 
      // logical volumes and parameters.
      
      
      string tower_name = _toString(EC_Number_of_towers,"tower%d");

      
      // build the ECRing layer, a box with hole in the middle 
      
      Tube CenterECTubForSi(0.,
			    Lcal_outer_radius + Ecal_Lcal_ring_gap + 0.001, // + tolerance,
			    module_thickness,
			    0.,
			    2 * M_PI);

      Box  ECRingSiBox( ECRingSiplateSize/ 2. - tolerance, ECRingSiplateSize/ 2. - tolerance, l_thickness/2.0-tolerance);

      SubtractionSolid ECRingSiSolid( ECRingSiBox, CenterECTubForSi);
	  

      Volume     l_vol(det_name+"_"+l_name+"_"+tower_name,ECRingSiSolid,air);
      DetElement layer(module_det, l_name+tower_name, det_id);
	  
      l_vol.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));
      
      
      
      // Loop over the sublayers or slices for this layer.
      int s_num = 1;
      double s_pos_z = -(l_thickness / 2);
      
      //--------------------------------------------------------------------------------
      // BuildECRing keep the Alveolus structure, the same as the Barrel and Endcap
      //--------------------------------------------------------------------------------
	
      for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	xml_comp_t x_slice = si;
	string     s_name  =  _toString(s_num,"slice%d");
	double     s_thick = x_slice.thickness();
	    
	double slab_dim_x = ECRingSiplateSize/2.-tolerance;
	double slab_dim_y = s_thick/2.-tolerance;
	double slab_dim_z = ECRingSiplateSize/2.-tolerance;
	    
	Box        s_box(slab_dim_x,slab_dim_z,slab_dim_y);

	SubtractionSolid ECRingSiSliceSolid( s_box, CenterECTubForSi);

	Volume     s_vol(det_name+"_"+l_name+"_"+s_name,ECRingSiSliceSolid,lcdd.material(x_slice.materialStr()));
	DetElement slice(layer,s_name,det_id);
	    
	s_vol.setVisAttributes(lcdd.visAttributes(x_slice.visStr()));

#ifdef VERBOSE	    
	std::cout<<"x_slice.materialStr(): "<< x_slice.materialStr() <<std::endl;
#endif
	if (x_slice.materialStr().compare(x_staves.materialStr()) == 0)
	  radiator_dim_y = s_thick;
	// W StructureLayer has the same thickness as W radiator layer in the Alveolus layer
	    
	if ( x_slice.isSensitive() ) {
	  s_vol.setSensitiveDetector(sens);
	}
	slice.setAttributes(lcdd,s_vol,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	
	// Slice placement.
	PlacedVolume slice_phv = l_vol.placeVolume(s_vol,Position(0,0,s_pos_z+s_thick/2));
	
	if ( x_slice.isSensitive() ) {
	  slice_phv.addPhysVolID("slice",s_num);
	}
	slice.setPlacement(slice_phv);
	// Increment Z position of slice.
	s_pos_z += s_thick;
	
	// Increment slice number.
	++s_num;
      }        
      
      
      if(radiator_dim_y <= 0) {
	stringstream err;
	err << " \n ERROR: The subdetector " << x_det.nameStr() << " geometry parameter -- radiator_dim_y = " << radiator_dim_y ;
	err << " \n Please check the radiator material name in the subdetector xml file";
	throw runtime_error(err.str());
      }
      
      // #########################
      // BuildECRingStructureLayer
      // #########################
      
	  
      string bs_name="bs";
      
      //Box        EndcapStructureLayer_box(radiator_dim_x/2.,radiator_dim_z/2.,radiator_dim_y/2.);
      Box        EndcapStructureLayer_box(ECRingSiplateSize/ 2. - tolerance, ECRingSiplateSize/ 2. - tolerance, radiator_dim_y/2.);
      SubtractionSolid  EndcapStructureLayerSolid( EndcapStructureLayer_box, CenterECTubForSi);
      
      Volume     EndcapStructureLayer_vol(det_name+"_"+l_name+"_"+bs_name,EndcapStructureLayerSolid,lcdd.material(x_staves.materialStr()));
      
      EndcapStructureLayer_vol.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));	
      
  
      int i_stave = 1;
      int i_tower = 1;
      
      Position l_pos(0,0,l_pos_z);
      RotationZYX rot(0,0,0);
      Transform3D tran3D(rot,l_pos);
	      
      PlacedVolume layer_phv = EnvLogECRing.placeVolume(l_vol,tran3D);
      layer_phv.addPhysVolID("layer", l_num);
      layer_phv.addPhysVolID("tower", i_tower);
      layer_phv.addPhysVolID("stave", i_stave);
      layer.setPlacement(layer_phv);
	      
      
      // Without last W StructureLayer, the last part is Si SD even layer.
      // the last number of  Ecal_nlayers1, Ecal_nlayers2 and  Ecal_nlayers3 is odd.
      int even_layer = l_num*2;
      if(even_layer > Ecal_nlayers1 + Ecal_nlayers2 + Ecal_nlayers3) continue;
      
      double bsl_pos_z = l_pos_z + l_thickness/2. + (radiator_dim_y/2. + Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));
      
      Position   bsl_pos(0,0,bsl_pos_z);
      Transform3D bsl_tran3D(rot,bsl_pos);
      PlacedVolume  EndcapStructureLayer_phv = EnvLogECRing.placeVolume(EndcapStructureLayer_vol,bsl_tran3D);
	      
      //-----------------------------------------------------------------------------------------
      DDRec::LayeredCalorimeterData::Layer caloLayer ;
      
      caloLayer.distance = EcalEndcapRing_min_z + l_pos_z ;
      caloLayer.thickness = l_thickness + radiator_dim_y ;
      caloLayer.absorberThickness = radiator_dim_y ;
      caloLayer.cellSize0 = Ecal_cells_size ;
      caloLayer.cellSize1 = Ecal_cells_size ;
      
      caloData->layers.push_back( caloLayer ) ;
      //-----------------------------------------------------------------------------------------
     
      // Increment to next layer Z position.
      l_pos_z +=   (l_thickness/2. +(radiator_dim_y/2. + Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE))*2.);
      
      ++l_num;
      
    }
  }
 
  
  // Set stave visualization.
  if (x_staves)   {
    EnvLogECRing.setVisAttributes(lcdd.visAttributes(x_staves.visStr()));
  }
  

  //====================================================================
  // Place Ecal Endcap module into the assembly envelope volume
  //====================================================================
  

  double EC_module_z_offset = Ecal_Barrel_module_dim_z * 2.5 + Ecal_cables_gap + module_thickness /2;
  
  for(int module_num=0;module_num<2;module_num++) {

    int module_id = ( module_num == 0 ) ? 0:6;
    double this_module_z_offset = ( module_id == 0 ) ? - EC_module_z_offset : EC_module_z_offset; 
    double this_module_rotY = ( module_id == 0 ) ? M_PI:0; 
  
    Position xyzVec(0,0,this_module_z_offset);
    RotationZYX rot(0,this_module_rotY,0);
    Rotation3D rot3D(rot);
    Transform3D tran3D(rot3D,xyzVec);

    PlacedVolume pv = envelope.placeVolume(EnvLogECRing,tran3D);
    pv.addPhysVolID("module",module_id); // z: -/+ 0/6

    DetElement sd = (module_num==0) ? module_det : module_det.clone(_toString(module_num,"module%d"));
    sd.setPlacement(pv);

  }
  
  sdet.addExtension< DDRec::LayeredCalorimeterData >( caloData ) ; 

  return sdet;
  
}



DECLARE_DETELEMENT(SEcal04_ECRing, create_detector)


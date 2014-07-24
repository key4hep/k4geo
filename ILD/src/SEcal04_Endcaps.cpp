//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for SiWEcalEndcaps
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
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
// $Id: SEcal04.cc,v 1.1 2009/06/04 11:40:42 musat Exp $
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
		Ecal Barrel. They can be ssen here again.
 */

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "DD4hep/Shapes.h"

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

  Assembly envelope_assembly( det_name + "assembly"  ) ;  
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope_assembly);

  env_phv.addPhysVolID("system",det_id);
  sdet.setPlacement(env_phv);

  sens.setType("calorimeter");

  DetElement    module_det("module0",det_id);


//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build Ecal endcaps
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
  
  double Ecal_inner_radius                  = lcdd.constant<double>("TPC_outer_radius") +lcdd.constant<double>("Ecal_Tpc_gap");
  double Ecal_radiator_thickness1           = lcdd.constant<double>("Ecal_radiator_layers_set1_thickness");
  double Ecal_radiator_thickness2           = lcdd.constant<double>("Ecal_radiator_layers_set2_thickness");
  double Ecal_radiator_thickness3           = lcdd.constant<double>("Ecal_radiator_layers_set3_thickness");
  double Ecal_Barrel_halfZ                  = lcdd.constant<double>("Ecal_Barrel_halfZ");
  
  double Ecal_support_thickness             = lcdd.constant<double>("Ecal_support_thickness");
  double Ecal_front_face_thickness          = lcdd.constant<double>("Ecal_front_face_thickness");
  double Ecal_lateral_face_thickness        = lcdd.constant<double>("Ecal_lateral_face_thickness");
  double Ecal_Slab_H_fiber_thickness        = lcdd.constant<double>("Ecal_Slab_H_fiber_thickness");

  double Ecal_Slab_Sc_PCB_thickness         = lcdd.constant<double>("Ecal_Slab_Sc_PCB_thickness");
  double Ecal_Sc_thickness                  = lcdd.constant<double>("Ecal_Sc_thickness");
  double Ecal_Sc_reflector_thickness        = lcdd.constant<double>("Ecal_Sc_reflector_thickness");

  double Ecal_endcap_extra_size             = lcdd.constant<double>("Ecal_endcap_extra_size");
  double Ecal_cables_gap                    = lcdd.constant<double>("Ecal_cables_gap");
  double Lcal_outer_radius                  = lcdd.constant<double>("Lcal_outer_radius");
  double Ecal_Lcal_ring_gap                 = lcdd.constant<double>("Ecal_Lcal_ring_gap");
  double Ecal_endcap_center_box_size        = lcdd.constant<double>("Ecal_endcap_center_box_size");

  int    Ecal_nlayers1                      = lcdd.constant<int>("Ecal_nlayers1");
  int    Ecal_nlayers2                      = lcdd.constant<int>("Ecal_nlayers2");
  int    Ecal_nlayers3                      = lcdd.constant<int>("Ecal_nlayers3");
  int    Ecal_barrel_number_of_towers       = lcdd.constant<int>("Ecal_barrel_number_of_towers");
  

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

  // The alveolus size takes in account the module Z size
  // but also 4 fiber layers for the alveoulus wall, the all
  // divided by the number of towers
  double alveolus_dim_z = 
    (Ecal_Barrel_module_dim_z - 2. * Ecal_lateral_face_thickness) /
    Ecal_barrel_number_of_towers - 
    2 * N_FIBERS_ALVOULUS  * Ecal_fiber_thickness  - 
    2 * Ecal_Slab_H_fiber_thickness -
    2 * Ecal_Slab_shielding;


#ifdef VERBOSE
  std::cout << "alveolus_dim_z = " <<  alveolus_dim_z << std::endl;
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
  

// ========= Create Ecal endcaps   ====================================
//  It will be the volume for palcing the Ecal endcaps alveolus(i.e. Layers).
//  And the structure W plate.
//  Itself will be placed into the world volume.
// ==========================================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~              EndcapStandardModule                 ~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  double Ecal_endcap_rmax  = Ecal_inner_radius 
    + module_thickness
    + Ecal_endcap_extra_size;

  double rInner = 0;
  double rOuter = Ecal_endcap_rmax;
  double zPlane = module_thickness;

  double Ecal_endcap_Tube_rmax = Lcal_outer_radius + Ecal_Lcal_ring_gap;
  
  PolyhedraRegular ECPolyHedra(8, M_PI/8.,rInner, rOuter, zPlane);
  //Tube CenterECTub(0., Ecal_endcap_Tube_rmax, module_thickness+0.0001);
  //SubtractionSolid EndCapSolid( ECPolyHedra, CenterECTub);

  Box  CenterECBox(Ecal_endcap_center_box_size/2.,Ecal_endcap_center_box_size/2.,module_thickness+0.0001);
  SubtractionSolid EndCapSolid( ECPolyHedra, CenterECBox);

  Volume EnvLogEndCap("endcap",EndCapSolid,air);
 



  //-------------------------------------------------------
  // Radiator and towers placements inside the Endcap module
  //-------------------------------------------------------

  // We count the layers starting from IP and from 1,
  // so odd layers should be inside slabs and
  // even ones on the structure.
  //
  double z_floor = 
    - module_thickness/2 +
    Ecal_front_face_thickness + 
    N_FIBERS_ALVOULUS * Ecal_fiber_thickness;

  double EC_y_bottom =
    + Ecal_endcap_center_box_size / 2
    + Ecal_lateral_face_thickness;

  double y_middle =
    (Ecal_inner_radius + module_thickness + Ecal_endcap_extra_size) 
    * tan(M_PI/8)
    - Ecal_lateral_face_thickness;

  double EC_y_top =
    Ecal_inner_radius + module_thickness + Ecal_endcap_extra_size
    - Ecal_lateral_face_thickness;

  double x_left =
    -(Ecal_inner_radius + module_thickness + Ecal_endcap_extra_size)
    + Ecal_lateral_face_thickness; 
 
  double x_right =
    + Ecal_endcap_center_box_size / 2
    - Ecal_lateral_face_thickness;

  double fiber_inter_alveolus =
    2 * (Ecal_Slab_H_fiber_thickness 
	 + Ecal_Slab_shielding
	 + N_FIBERS_ALVOULUS * Ecal_fiber_thickness);
  
  double EC_alveolus_dim_y = alveolus_dim_z;
  double EC_alveolus_x_left = 0;
  double EC_alveolus_dim_x = 0 ;
  double alv_upper_y = 0;
  double inc = (-y_middle -x_left ) / (EC_y_top - y_middle);

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
      double y_floor = EC_y_bottom;
      double radiator_dim_y = -1.0; //to be updated with slice radiator thickness 
      
      while ( ( y_floor + EC_alveolus_dim_y) < EC_y_top )
	{
	  alv_upper_y = y_floor + EC_alveolus_dim_y;
	  
	  EC_alveolus_dim_x = x_right- x_left;
	  
	  if( alv_upper_y <= y_middle )
	    {
	      EC_alveolus_dim_x = x_right- x_left;
	    }
	  else
	    {
	      EC_alveolus_x_left = 
		(alv_upper_y - y_middle) * inc + x_left;
	      EC_alveolus_dim_x = x_right
		- EC_alveolus_x_left;
	    }
	  
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
	  
	  Box        l_box(alveolus_dim_z/2.-tolerance,EC_alveolus_dim_x/2.-tolerance,l_thickness/2.0-tolerance);
	  Volume     l_vol(det_name+"_"+l_name+"_"+tower_name,l_box,air);
	  DetElement layer(module_det, l_name+tower_name, det_id);
	  
	  l_vol.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));
	  
	  
	  // Loop over the sublayers or slices for this layer.
	  int s_num = 1;
	  double s_pos_z = -(l_thickness / 2);
	  
	  //--------------------------------------------------------------------------------
	  // BuildEndcapAlveolus: BuildSiliconSlab:
	  //--------------------------------------------------------------------------------
	
	  for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	    xml_comp_t x_slice = si;
	    string     s_name  =  _toString(s_num,"slice%d");
	    double     s_thick = x_slice.thickness();
	    
	    double slab_dim_x = alveolus_dim_z/2.-tolerance;
	    double slab_dim_y = s_thick/2.-tolerance;
	    double slab_dim_z = EC_alveolus_dim_x/2.-tolerance;
	    
	    Box        s_box(slab_dim_x,slab_dim_z,slab_dim_y);
	    Volume     s_vol(det_name+"_"+l_name+"_"+s_name,s_box,lcdd.material(x_slice.materialStr()));
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
	  // BuildEndcapStructureLayer
	  // #########################
	  
	  double radiator_dim_x = alveolus_dim_z;
	  
#ifdef VERBOSE
	  std::cout << "radiator_dim_x = " << radiator_dim_x << std::endl;
#endif  
	  
	  double radiator_dim_z = EC_alveolus_dim_x;
	  
	  string bs_name="bs";
	  
	  Box        EndcapStructureLayer_box(radiator_dim_x/2.,radiator_dim_z/2.,radiator_dim_y/2.);
	  Volume     EndcapStructureLayer_vol(det_name+"_"+l_name+"_"+bs_name,EndcapStructureLayer_box,lcdd.material(x_staves.materialStr()));
	  
	  EndcapStructureLayer_vol.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));	
	  
	  
	  int limit_staves;
	  limit_staves = 4;
	  // **********************
	  for (int i_stave = 1;
	       i_stave <= limit_staves;
	       i_stave ++)
	    {
	      double angle_module = M_PI/2. * ( i_stave - 1 );
	      
	      Position l_pos(-(y_floor + EC_alveolus_dim_y/2),-(-EC_alveolus_dim_x/2 + x_right),l_pos_z);
	      RotationZ rotz(angle_module);
	      Position l_new = rotz*l_pos;   // The position of EC_alveolus in each stave, place same EC_alveolus into 4 stave.
	      
	      RotationZYX rot(angle_module,0,0);
	      Transform3D tran3D(rot,l_new); // The rotation of EC_alveolus in each stave, place same EC_alveolus into 4 stave.
	      
	      PlacedVolume layer_phv = EnvLogEndCap.placeVolume(l_vol,tran3D);
	      layer_phv.addPhysVolID("layer", l_num);
	      layer_phv.addPhysVolID("tower", EC_Number_of_towers);
	      layer_phv.addPhysVolID("stave", i_stave);
	      layer.setPlacement(layer_phv);
	      
#ifdef VERBOSE
	      std::cout<< "Ecal_endcap_towers: "<< EC_Number_of_towers
		       << " has been placed into EnvLogEndCap stave "<<i_stave<<"."
		       << std::endl;
#endif	    
	      
	      // Without last W StructureLayer, the last part is Si SD even layer.
	      // the last number of  Ecal_nlayers1, Ecal_nlayers2 and  Ecal_nlayers3 is odd.
	      int even_layer = l_num*2;
	      if(even_layer > Ecal_nlayers1 + Ecal_nlayers2 + Ecal_nlayers3) continue;
	      
	      double bsl_pos_z = l_pos_z + l_thickness/2. + (radiator_dim_y/2. + Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));
	      
	      Position   bsl_pos(-(y_floor + EC_alveolus_dim_y/2),-(-EC_alveolus_dim_x/2 + x_right),bsl_pos_z);      // Position of the layer.
	      Position bsl_new = rotz*bsl_pos;
	      Transform3D bsl_tran3D(rot,bsl_new);
	      PlacedVolume  EndcapStructureLayer_phv = EnvLogEndCap.placeVolume(EndcapStructureLayer_vol,bsl_tran3D);
	      
	    }
	  
	  y_floor += EC_alveolus_dim_y + fiber_inter_alveolus;
	  
	  EC_Number_of_towers++;
	  
	}
      
      // Increment to next layer Z position.
      l_pos_z +=   (l_thickness/2. +(radiator_dim_y/2. + Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE))*2.);
      
      ++l_num;
      
    }
  }
    
  
    // Set stave visualization.
    if (x_staves)   {
      EnvLogEndCap.setVisAttributes(lcdd.visAttributes(x_staves.visStr()));
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

    PlacedVolume pv = envelope_assembly.placeVolume(EnvLogEndCap,tran3D);
    pv.addPhysVolID("module",module_id); // z: -/+ 0/6

    DetElement sd = (module_num==0) ? module_det : module_det.clone(_toString(module_num,"module%d"));
    sd.setPlacement(pv);

  }
  

  return sdet;
  
}



DECLARE_DETELEMENT(SEcal04_Endcaps, create_detector)


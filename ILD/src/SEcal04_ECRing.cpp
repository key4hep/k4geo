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
		Ecal Barrel. They can be seen here again.
		Start ECRing ...
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
  
  //double Ecal_inner_radius                  = lcdd.constant<double>("TPC_outer_radius") +lcdd.constant<double>("Ecal_Tpc_gap");
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

  double Ecal_EC_Ring_gap                   = lcdd.constant<double>("Ecal_EC_Ring_gap");

  //double Ecal_endcap_extra_size             = lcdd.constant<double>("Ecal_endcap_extra_size");
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
  Box  CenterECBox(Ecal_endcap_center_box_size/2.,Ecal_endcap_center_box_size/2.,module_thickness/2.);

  SubtractionSolid ECRingSolid( CenterECBox, CenterECTub);

  Volume EnvLogECRing("ECRing",ECRingSolid,air);
 


  //==============================================================
  // TODO: to build the layer and place into the ECRing
  //==============================================================
















  
    // Set stave visualization.
    if (x_staves)   {
      EnvLogECRing.setVisAttributes(lcdd.visAttributes(x_staves.visStr()));
    }
 

  //====================================================================
  // Place Ecal Endcap module into the assembly envelope volume
  //====================================================================
  

  double EC_module_z_offset = Ecal_Barrel_module_dim_z * 2.5 + Ecal_cables_gap + module_thickness /2;
  
  for(int module_num=0;module_num<2;module_num++) {

    int module_id = module_num;
    double this_module_z_offset = ( module_id == 0 ) ? EC_module_z_offset : - EC_module_z_offset; 
    double this_module_rotY = ( module_id == 0 ) ? 0:M_PI; 
  
    Position xyzVec(0,0,this_module_z_offset);
    RotationZYX rot(0,this_module_rotY,0);
    Rotation3D rot3D(rot);
    Transform3D tran3D(rot3D,xyzVec);

    PlacedVolume pv = envelope_assembly.placeVolume(EnvLogECRing,tran3D);
    pv.addPhysVolID("module",module_id); // z: +/-

    DetElement sd = (module_num==0) ? module_det : module_det.clone(_toString(module_num,"module%d"));
    sd.setPlacement(pv);

  }
  

  return sdet;
  
}



DECLARE_DETELEMENT(SEcal04_ECRing, create_detector)


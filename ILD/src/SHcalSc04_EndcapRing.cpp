//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for HcalEndcapRing
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================

 /* History:  
  
   initial version: 
   F.Gaede: identical to  Hcal03 driver except that an additional gap
            for the fibres is introduced between scintillator and steel
   PMoradeFreitas: Super driver without tmp database and able
                   to build Hcal barrel with just two modules in stave.
   Angela Lucaci: similar to SHcal03, just that the drift chambers option is 
                  not considered anymore, only the scintillator one. In addition,
                  a gap in the middle of the stave is build, and fractional cells
                  at the edges are introduced (see SDHcalBarrel.cc)
   Ralf Diener: Corrected use of GEAR interface
   Andre Sailer: Added Tungsten
   Andre Sailer: Added possibility for different endcap/barrel material
                 Which needs three additional parameters: 
		 endcap_radiator_thickness, endcap_radiator_material, endcap_layers
                 Also added parameters to the database for this driver
   Shaojun Lu:  Barrel driver has been changed for the new engineering design shape.
                Endcaps driver has been changed for the new engineering design shape.
		Updated for flexibility.
   Shaojun Lu:  Ported from Mokka SHcalSc04 EndcapRing part. Read the constants from XML
                instead of the DB. Then build the EndcapRing in the same way with DD4hep
		construct.
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

  xml_det_t   x_det     = element;
  string      det_name    = x_det.nameStr();
  Layering    layering(x_det);

  Material    air         = lcdd.air();
  Material    Steel235    = lcdd.material(x_det.materialStr());

  int           det_id    = x_det.id();
  xml_comp_t    x_staves  = x_det.staves();
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);


  Assembly envelope_assembly( det_name + "assembly"  ) ;
  PlacedVolume env_phv   = motherVol.placeVolume(envelope_assembly);

  env_phv.addPhysVolID("system",det_id);
  sdet.setPlacement(env_phv);

  sens.setType("calorimeter");

  DetElement    module_det("module0",det_id);

  // Some verbose output
  cout << " \n\n\n CREATE DETECTOR: SHcalSC04_EndcapRing" << endl;




//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build HcalEndcapRing
//
//====================================================================
  // The way to read constant from XML/LCDD file.
  double      Hcal_radiator_thickness          = lcdd.constant<double>("Hcal_radiator_thickness");
  double      Hcal_chamber_thickness           = lcdd.constant<double>("Hcal_chamber_thickness");
  double      Hcal_inner_radius                = lcdd.constant<double>("Hcal_inner_radius");
  double      Hcal_back_plate_thickness        = lcdd.constant<double>("Hcal_back_plate_thickness");
  double      Hcal_lateral_plate_thickness     = lcdd.constant<double>("Hcal_lateral_structure_thickness");
  double      Hcal_stave_gaps                  = lcdd.constant<double>("Hcal_stave_gaps");
  double      Hcal_modules_gap                 = lcdd.constant<double>("Hcal_modules_gap"); 
  double      TPC_Ecal_Hcal_barrel_halfZ       = lcdd.constant<double>("TPC_Ecal_Hcal_barrel_halfZ"); 
  double      Hcal_middle_stave_gaps           = lcdd.constant<double>("Hcal_middle_stave_gaps");
  double      Hcal_layer_air_gap               = lcdd.constant<double>("Hcal_layer_air_gap");

  double      Hcal_module_radius               = lcdd.constant<double>("Hcal_outer_radius");

  int         Hcal_nlayers                     = lcdd.constant<int>("Hcal_nlayers");

  //double      Hcal_inner_radius                = lcdd.constant<double>("HcalBarrel_rmin");
  //double      Hcal_outer_radius                = lcdd.constant<double>("Hcal_outer_radius");
  //double      Hcal_Barrel_rotation             = lcdd.constant<double>("Hcal_Barrel_rotation");

  double      Ecal_endcap_zmin                 = lcdd.constant<double>("Ecal_endcap_zmin");
  double      Hcal_endcap_radiator_thickness   = lcdd.constant<double>("Hcal_endcap_radiator_thickness");
  double      Hcal_radial_ring_inner_gap       = lcdd.constant<double>("Hcal_radial_ring_inner_gap");
  double      Hcal_endcap_cables_gap           = lcdd.constant<double>("Hcal_endcap_cables_gap");
  double      Hcal_endcap_ecal_gap             = lcdd.constant<double>("Hcal_endcap_ecal_gap");

  int         Hcal_ring                        = lcdd.constant<int>("Hcal_ring");

  //TODO: thinking about how to pass the updated value at runtime from other inner drivers? 
  double      Ecal_endcap_zmax                 = 2635*mm;  //= lcdd.constant<double>("Ecal_endcap_zmax");
  double      Ecal_endcap_outer_radius         = 2088.8*mm;//= lcdd.constant<double>("Ecal_endcap_outer_radius");


//====================================================================
//
// general calculated parameters
//
//====================================================================

  double Hcal_total_dim_y   = Hcal_nlayers * (Hcal_radiator_thickness + Hcal_chamber_thickness) 
                            + Hcal_back_plate_thickness;

  
  // Moved the calculation into ILD_o1_v05.xml
  // double Hcal_module_radius = Hcal_inner_radius /  cos(M_PI/8.)
  //   + ( Hcal_radiator_thickness + Hcal_chamber_thickness ) * Hcal_nlayers
  //   + Hcal_back_plate_thickness;
  // std::cout << " *********** hcal outer radius : " << Hcal_module_radius  << std::endl ;


  double Hcal_y_dim1_for_x  = Hcal_module_radius*cos(M_PI/8.) - Hcal_inner_radius;
  double Hcal_bottom_dim_x  = 2.*Hcal_inner_radius*tan(M_PI/8.)- Hcal_stave_gaps;
  double Hcal_normal_dim_z  = (2 * TPC_Ecal_Hcal_barrel_halfZ - Hcal_modules_gap)/2.;

 //only the middle has the steel plate.
  double Hcal_regular_chamber_dim_z = Hcal_normal_dim_z - Hcal_lateral_plate_thickness;


  //double totalThickness_Hcal_Barrel       = (Hcal_radiator_thickness + Hcal_chamber_thickness) 
  //                                                 * HcalBarrel_layers + Hcal_back_plate_thickness;
  //double Hcal_module_radius               = Hcal_inner_radius + totalThickness_Hcal_Barrel;

  //double Hcal_cell_dim_x            = Hcal_cells_size;
  //double Hcal_cell_dim_z            = Hcal_regular_chamber_dim_z / floor (Hcal_regular_chamber_dim_z/Hcal_cell_dim_x);
  //double Hcal_digi_cell_dim_x       = Hcal_digi_cells_size;
  //double Hcal_digi_cell_dim_z       = Hcal_regular_chamber_dim_z / floor (Hcal_regular_chamber_dim_z/Hcal_digi_cell_dim_x);

 

  // just two modules per stave
  //double Hcal_normal_dim_z = (2 * TPC_Ecal_Hcal_barrel_halfZ - Hcal_modules_gap)/2.;
  double Hcal_start_z =  Hcal_normal_dim_z + Hcal_modules_gap / 2. + Hcal_endcap_cables_gap;

  cout<<"    HCAL endcap rings: Hcal_start_z:  "<< Hcal_start_z <<endl<<endl;

  // Hcal_start_z is the Hcal Endcap boundary coming from the IP
  // Test Hcal_start_z against Ecal_endcap_zmax + Hcal_endcap_ecal_gap
  // to avoid overlap problems with Ecal if scaled.
  //
  double Hcal_endcap_rmax   = Hcal_inner_radius + Hcal_y_dim1_for_x;

  if( Hcal_start_z < Ecal_endcap_zmax + Hcal_endcap_ecal_gap )
    Hcal_start_z = Ecal_endcap_zmax + Hcal_endcap_ecal_gap;

  cout<<"    HCAL endcap rings: Hcal_start_z:  "<< Hcal_start_z <<endl<<endl;






// ========= Create Hcal end cap ring   ====================================
//  It will be the volume for palcing the Hcal endcaps Ring Layers.
//  And the absorber plate.
//  Itself will be placed into the world volume.
// ==========================================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~              EndcapRings                          ~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // old parameters from database
  double pRMax, pDz, pRMin;
  pRMax = Hcal_endcap_rmax;

  // The rings start from inner Ecal endcap boundary
  // and finish at inner Hcal endcap one.
  double start_z, stop_z;
  start_z = Ecal_endcap_zmin;
  double SpaceForLayers = Hcal_start_z - Hcal_endcap_ecal_gap
    - Ecal_endcap_zmin - Hcal_back_plate_thickness;

  int MaxNumberOfLayers = (int) (SpaceForLayers /
				     (Hcal_chamber_thickness + Hcal_endcap_radiator_thickness));


  cout<<"    HCAL endcap rings: Hcal_start_z:  "<< Hcal_start_z <<endl;
  cout<<"    HCAL endcap rings: Hcal_endcap_rmax: "<< Hcal_endcap_rmax <<endl;
  cout<<"    HCAL endcap rings: SpaceForLayers:  "<< SpaceForLayers <<endl;
  cout<<"    HCAL endcap rings will have "<< MaxNumberOfLayers << " layers."<<endl;

  stop_z = start_z + MaxNumberOfLayers * (Hcal_chamber_thickness + Hcal_endcap_radiator_thickness)
    + Hcal_back_plate_thickness;

  pDz = (stop_z - start_z) / 2.;

  pRMin = Ecal_endcap_outer_radius
    + Hcal_radial_ring_inner_gap;

  double zPlane[2];
  zPlane[0]=-pDz;
  zPlane[1]=-zPlane[0];

  double zlen = pDz*2.;

  double rInner[2],rOuter[2];
  rInner[0]=rInner[1]=pRMin;
  rOuter[0]=rOuter[1]=pRMax;

  double rmin = pRMin;
  double rmax = pRMax;
  

  PolyhedraRegular HcalEndCapRingSolid( 8, M_PI/8., rmin, rmax,  zlen);

  Volume  HcalEndCapRingLogical("HcalEndCapRingLogical",HcalEndCapRingSolid, Steel235);



  //==============================================================
  // build the layer and place into the Hcal EndcapRing
  //==============================================================























  // Set stave visualization.
  if (x_staves)   {
    HcalEndCapRingLogical.setVisAttributes(lcdd.visAttributes(x_staves.visStr()));
   }
  

  //====================================================================
  // Place Hcal Endcap Ring module into the assembly envelope volume
  //====================================================================
  
  double endcap_z_offset = Ecal_endcap_zmin + pDz;
  
  for(int module_num=0;module_num<2;module_num++) {

    int module_id = module_num;
    double this_module_z_offset = ( module_id == 0 ) ? endcap_z_offset : - endcap_z_offset; 
    double this_module_rotY = ( module_id == 0 ) ? 0:M_PI; 
  
    Position xyzVec(0,0,this_module_z_offset);
    RotationZYX rot(0,this_module_rotY,0);
    Rotation3D rot3D(rot);
    Transform3D tran3D(rot3D,xyzVec);

    PlacedVolume pv = envelope_assembly.placeVolume(HcalEndCapRingLogical,tran3D);
    pv.addPhysVolID("module",module_id); // z: +/-

    DetElement sd = (module_num==0) ? module_det : module_det.clone(_toString(module_num,"module%d"));
    sd.setPlacement(pv);

  }
  

  return sdet;
  
}



DECLARE_DETELEMENT(SHcalSc04_EndcapRing, create_detector)

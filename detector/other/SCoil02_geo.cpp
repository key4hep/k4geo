//====================================================================
//  DDSim - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/Surface.h"
#include "XMLHandlerDB.h"
#include "XML/Utilities.h"
#include <cmath>
#include "DDRec/DetectorData.h"

//#include "GearWrapper.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Tube;
using dd4hep::Volume;

using dd4hep::xml::_toString;

using dd4hep::rec::LayeredCalorimeterData;


/** Construction of the coil, ported from Mokka drivers SCoil02.cc and Coil01.cc
 *
 * Mokka History:  
 * SCoil02.cc
 * F.Gaede, DESY:  based on SCoil01, with added parameters for the 
 *                 clearance to the yoke:
 *     Hcal_Coil_additional_gap  : adjust the actual gap in r (Hcal_R_max allready defines a gap)
 *     Coil_Yoke_radial_clearance :  -> defines inner r of yoke
 *     Coil_Yoke_lateral_clearance : -> defines zEndcap of yoke
 * Coil00.cc
 * - first implementation P. Mora de Freitas (may 01)
 * - F.Gaede:  write out parameters to GEAR  (Oct 08)
 * Coil00.cc:
 * - V. Saveliev: replaced simple Al-tube with more detailed structure (cryostat, etc )
 * 
 *  @author: F.Gaede, DESY, Aug 2014
 */

static Ref_t create_element(Detector& theDetector, xml_h element, SensitiveDetector /*sens*/)  {

  //------------------------------------------
  //  See comments starting with '//**' for
  //     hints on porting issues
  //------------------------------------------

  xml_det_t    x_det = element;
  string       name  = x_det.nameStr();
    
  DetElement   coil(  name, x_det.id()  ) ;

 // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , coil ) ;

  dd4hep::xml::setDetectorTypeFlag( element, coil ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return coil ;

  //-----------------------------------------------------------------------------------

  PlacedVolume pv;



  //######################################################################################################################################################################
  //  code ported from Coil02::construct() :
  //##################################
  
  cout << "\nBuilding Coil..." << endl;

  xml_comp_t  x_tube (x_det.child(_U(tube)));

  double inner_radius = x_tube.rmin() ;
  double outer_radius = x_tube.rmax() ;
  double half_z       = x_tube.dz() ;

  cout << "\n... cryostat inner_radius " << inner_radius 
       << "\n... cryostat outer_radius " << outer_radius
       << "\n... cryostat half_z       " << half_z
       << endl;
  

  Material coilMaterial = theDetector.material( x_tube.materialStr() ) ;


  //FG: for now fall back to a simple tube filled with Al
  //    the code below has to many hard coded numbers
  //    that need verification and then need to be converted
  //    to xml parameters ... 

#if !code_is_cleaned_up 

  Tube   coil_tube( x_tube.rmin(), x_tube.rmax(), x_tube.dz() );

  Volume coil_vol( "coil_vol", coil_tube , coilMaterial );
  pv  =  envelope.placeVolume( coil_vol ) ;
  coil.setVisAttributes( theDetector, "BlueVis" , coil_vol );

  cout << " ... for the time being simply use a tube of aluminum ..." << endl ;

  //=========================================================================================================
#else

  //... Coil Cryostat (Al, inside vacuum)
  //... inner cylinder
  Tube CoilEnvelopeSolid_1( inner_radius, 40.* dd4hep::mm + inner_radius , half_z ) ;

  Volume CoilLogical_1( "CoilEnvelope_1", CoilEnvelopeSolid_1, coilMaterial ) ;

  coil.setVisAttributes( theDetector, "BlueVis" , CoilLogical_1 );

  pv  = envelope.placeVolume( CoilLogical_1 ) ;

  //... outer cylinder
  Tube CoilEnvelopeSolid_2 ( -30*dd4hep::mm + outer_radius, outer_radius , half_z ) ;

  Volume CoilLogical_2( "CoilEnvelope_2", CoilEnvelopeSolid_2, coilMaterial ) ;

  coil.setVisAttributes( theDetector, "BlueVis" , CoilLogical_2 );

  pv  = envelope.placeVolume( CoilLogical_2 ) ;


  //... side wall left
  Tube CoilEnvelopeSolid_3( 40*dd4hep::mm + inner_radius, -30*dd4hep::mm + outer_radius, 25.* dd4hep::mm ) ;
    
  Volume CoilLogical_3( "CoilEnvelope_3", CoilEnvelopeSolid_3, coilMaterial ) ;
  
  coil.setVisAttributes( theDetector, "BlueVis" , CoilLogical_3 );

  pv  = envelope.placeVolume( CoilLogical_3 , Position( 0., 0., -25.*dd4hep::mm + half_z )  ) ;


  //... side wall right
  // Tube CoilEnvelopeSolid_4( 40*dd4hep::mm + inner_radius, -30*dd4hep::mm + outer_radius, 25.* dd4hep::mm ) ;
  // Volume CoilLogical_4( "CoilEnvelope_4", CoilEnvelopeSolid_4, coilMaterial ) ;
  // coil.setVisAttributes( theDetector, "BlueVis" , CoilLogical_4 );
  // pv  = envelope.placeVolume( CoilLogical_4 , Position( 0., 0., 25.*dd4hep::mm - half_z )  ) ;
  //simply place the same volume again

  pv  = envelope.placeVolume( CoilLogical_3 , Position( 0., 0., 25.*dd4hep::mm - half_z )  ) ;


  //... Coil modules
  //... main coll module 1,2,3

  Tube CoilMainSolid_1( 175*dd4hep::mm + inner_radius, 560*dd4hep::mm + inner_radius,  1632.*dd4hep::mm/2-0.*dd4hep::mm ) ;

  Volume CoilMainLogical_1( "CoilMain_1" , CoilMainSolid_1, coilMaterial ) ;
  
  coil.setVisAttributes( theDetector, "GreyVis" , CoilMainLogical_1 );

  pv  = envelope.placeVolume( CoilMainLogical_1 , Position( 0., 0., -1632.*dd4hep::mm ) ) ;
  pv  = envelope.placeVolume( CoilMainLogical_1 , Position( 0., 0.,     0.            ) ) ;
  pv  = envelope.placeVolume( CoilMainLogical_1 , Position( 0., 0.,  1632.*dd4hep::mm ) ) ;




  // ######################################################################################################
  // ======================================================================================================
  //   => continue porting below .... 
  //
  // ######################################################################################################
  // ======================================================================================================


  //... corrected coil module 1,2
  Tubs *CoilCorrectSolid_1
    Tube ("CoilCorrect_1", 
		 175*dd4hep::mm + inner_radius,
		 560*dd4hep::mm + inner_radius, 
		 1224.*dd4hep::mm/2,
		 0*deg, 
		 360.*deg);
  
  Volume *CoilCorrectLogical_1=
    new Volume(CoilCorrectSolid_1,
			theDetector.material("aluminium"),
			"CoilCorrect_1", 
			0,
			0,
			0);
  
  VisAttributes * VisAtt_c1 = new VisAttributes(Colour(1.,0.5,0.5));
  VisAtt_c1->SetForceWireframe(true);
  //VisAtt_c1->SetForceSolid(true);
  CoilCorrectLogical_1  coil.setVisAttributes( theDetector, "VisAtt_c1");

  .setPlacement(0,
		    ThreeVector(0., 0., -1632.*dd4hep::mm*3/2-612.),
		    CoilCorrectLogical_1,
		    "CoilCorrect_1",
		    WorldLog,
		    false,0);

  .setPlacement(0,
		    ThreeVector(0., 0., 1632.*dd4hep::mm*3/2+612.),
		    CoilCorrectLogical_1,
		    "CoilCorrect_1",
		    WorldLog,
		    false,1);

//... Coil mandrel
  Tubs *CoilMandrelSolid
    Tube ("CoilMandrel", 
		 560*dd4hep::mm + inner_radius,
		 635*dd4hep::mm + inner_radius, 
		 3675.*dd4hep::mm,
		 0*deg, 
		 360.*deg);
  
  Volume *CoilMandrelLogical=
    new Volume(CoilMandrelSolid,
			theDetector.material("aluminium"),
			"CoilMandrel", 
			0,
			0,
			0);
  
  VisAttributes * VisAtt_md = new VisAttributes(Colour(0.7,0.7,0.7));
  VisAtt_md->SetForceWireframe(true);
  //VisAtt_md->SetForceSolid(true);
  CoilMandrelLogical  coil.setVisAttributes( theDetector, "VisAtt_md");

  .setPlacement(0,
		    ThreeVector(0., 0., 0.),
		    CoilMandrelLogical,
		    "CoilMandrel",
		    WorldLog,
		    false,0);

//... Coil sensitive detectors
  int layer_id=1; //Put in the database!!!
// Threshold is 20%. mip = 200 keV/dd4hep::mm 
  const double sensitive_thickness = 10. *dd4hep::mm;
  theCoilSD = 
    new TRKSD00("COIL", 
		sensitive_thickness * dd4hep::mm 
		* 200 * keV
		* 0.01);
  RegisterSensitiveDetector(theCoilSD);


//... Scintillator Detector layer 1

  const double zPozBarrelArray_1[2]   = {-100.*dd4hep::mm+half_z, 
				           +100.*dd4hep::mm-half_z};
  const double rInnerBarrelArray_1[2] = 
                                      {90.*dd4hep::mm+inner_radius, 
				       90.*dd4hep::mm+inner_radius};
  const double rOuterBarrelArray_1[2] = 
                                     {100.*dd4hep::mm+inner_radius, 
			              100.*dd4hep::mm+inner_radius};
//... Coil detector geometry


  Polyhedra *CoilScintSolid_1
    = new Polyhedra("CoilScint_1", 
		      0.*deg,
		      360.*deg,
		      24,
		      2,
		      zPozBarrelArray_1,
		      rInnerBarrelArray_1,
		      rOuterBarrelArray_1);
 
  Volume *CoilScintLogical_1=
    new Volume(CoilScintSolid_1,
			theDetector.material("polystyrene"),
			"CoilScint_1", 
			0,
			0,
			0);
  
  VisAttributes * VisAtt_s1 = new VisAttributes(Colour(0.5,0.5,1.));
  VisAtt_s1->SetForceWireframe(true);
  //VisAtt_s1->SetForceSolid(true);
  CoilScintLogical_1  coil.setVisAttributes( theDetector, "VisAtt_s1");
  CoilScintLogical_1->SetSensitiveDetector(theCoilSD);

  .setPlacement(0,
		    ThreeVector(0., 0., 0. ),
		    CoilScintLogical_1,
		    "CoilScint_1",
		    WorldLog,
		    false,layer_id);

  layer_id++;
  //... Scintillation Detector layer 2
  const double zPozBarrelArray_2[2]   = {-100.*dd4hep::mm+half_z, 
				          +100.*dd4hep::mm-half_z};
  const double rInnerBarrelArray_2[2] = 
                                     {105.*dd4hep::mm+inner_radius, 
				     105.*dd4hep::mm+inner_radius};
  const double rOuterBarrelArray_2[2] = 
                                     {115.*dd4hep::mm+inner_radius, 
				      115.*dd4hep::mm+inner_radius};

  Polyhedra *CoilScintSolid_2
    = new Polyhedra("CoilScint_2", 
		      0.*deg,
		      360.*deg,
		      24,
		      2,
		      zPozBarrelArray_2,
		      rInnerBarrelArray_2,
		      rOuterBarrelArray_2);
 
  Volume *CoilScintLogical_2=
    new Volume(CoilScintSolid_2,
			theDetector.material("polystyrene"),
			"CoilScint_2", 
			0,
			0,
			0);
  
  VisAttributes * VisAtt_s2 = new VisAttributes(Colour(0.75,0.75,1.));
  VisAtt_s2->SetForceWireframe(true);
  //VisAtt_s2->SetForceSolid(true);
  CoilScintLogical_2  coil.setVisAttributes( theDetector, "VisAtt_s2");
  CoilScintLogical_2->SetSensitiveDetector(theCoilSD);

  .setPlacement(0,
		    ThreeVector(0., 0., 0. ),
		    CoilScintLogical_2,
		    "CoilScint_2",
		    WorldLog,
		    false,layer_id);

  layer_id++;
//... Scint detector layer 3 

  const double zPozBarrelArray_3[2]   = {-100.*dd4hep::mm+half_z, 
				           +100.*dd4hep::mm-half_z};
  const double rInnerBarrelArray_3[2] = 
                                     {-90.*dd4hep::mm+outer_radius, 
				      -90.*dd4hep::mm+outer_radius};
  const double rOuterBarrelArray_3[2] = 
                                     {-80.*dd4hep::mm+outer_radius, 
				      -80.*dd4hep::mm+outer_radius};

  Polyhedra *CoilScintSolid_3
    = new Polyhedra("CoilScint_3", 
		      0.*deg,
		      360.*deg,
		      24,
		      2,
		      zPozBarrelArray_3,
		      rInnerBarrelArray_3,
		      rOuterBarrelArray_3);
 
  Volume *CoilScintLogical_3=
    new Volume(CoilScintSolid_3,
			theDetector.material("polystyrene"),
			"CoilScint_3", 
			0,
			0,
			0);
  
  VisAttributes * VisAtt_s3 = new VisAttributes(Colour(0.5,0.5,1.));
  VisAtt_s3->SetForceWireframe(true);
  //VisAtt_s3->SetForceSolid(true);
  CoilScintLogical_3  coil.setVisAttributes( theDetector, "VisAtt_s3");
  CoilScintLogical_3->SetSensitiveDetector(theCoilSD);

  .setPlacement(0,
		    ThreeVector(0., 0., 0. ),
		    CoilScintLogical_3,
		    "CoilScint_3",
		    WorldLog,
		    false,layer_id);
  layer_id++;
  //... Scintillation Detector layer 4
  const double zPozBarrelArray_4[2]   = 
                                  {-100.*dd4hep::mm+half_z, 
				   +100.*dd4hep::mm-half_z};
  const double rInnerBarrelArray_4[2] = 
                                  {-105.*dd4hep::mm+outer_radius, 
				   -105.*dd4hep::mm+outer_radius};
  const double rOuterBarrelArray_4[2] = 
                                  {-95.*dd4hep::mm+outer_radius, 
				   -95.*dd4hep::mm+outer_radius};

  Polyhedra *CoilScintSolid_4
    = new Polyhedra("CoilScint_4", 
		      0.*deg,
		      360.*deg,
		      24,
		      2,
		      zPozBarrelArray_4,
		      rInnerBarrelArray_4,
		      rOuterBarrelArray_4);
 
  Volume *CoilScintLogical_4=
    new Volume(CoilScintSolid_4,
			theDetector.material("polystyrene"),
			"CoilScint_4", 
			0,
			0,
			0);
  
  VisAttributes * VisAtt_s4 = new VisAttributes(Colour(0.75,0.75,1.));
  VisAtt_s4->SetForceWireframe(true);
  //VisAtt_s4->SetForceSolid(true);
  CoilScintLogical_4  coil.setVisAttributes( theDetector, "VisAtt_s4");
  CoilScintLogical_4->SetSensitiveDetector(theCoilSD);

  .setPlacement(0,
		    ThreeVector(0., 0., 0. ),
		    CoilScintLogical_4,
		    "CoilScint_4",
		    WorldLog,
		    false,layer_id);

#ifdef MOKKA_GEAR
  //----------------------------------------------------
  // MokkaGear
  //----------------------------------------------------

  MokkaGear* gearMgr = MokkaGear::getMgr() ;
  
  gear::GearParametersImpl* gp = new gear::GearParametersImpl ;

  //Inner Cylinder
  gp->setDoubleVal("Coil_cryostat_inner_cyl_inner_radius", 
                                        inner_radius);
  gp->setDoubleVal("Coil_cryostat_inner_cyl_outer_radius", 
                                        40.*dd4hep::mm+outer_radius);
  gp->setDoubleVal("Coil_cryostat_inner_cyl_half_z", half_z);
  gp->setStringVal("Coil_material_inner_cyl", "aluminium");


 //Outer Cylinder
  gp->setDoubleVal("Coil_cryostat_outer_cyl_inner_radius", 
                                        -30*dd4hep::mm+outer_radius);
  gp->setDoubleVal("Coil_cryostat_outer_cyl_outer_radius", 
                                        outer_radius);
  gp->setDoubleVal("Coil_cryostat_outer_cyl_half_z", half_z);
  gp->setStringVal("Coil_material_outer_cyl", "aluminium");

  //FG: add the parameters under the 'old' names as expected by the reconstruction:
  gp->setDoubleVal("Coil_cryostat_inner_radius", inner_radius);
  gp->setDoubleVal("Coil_cryostat_outer_radius", outer_radius);
  gp->setDoubleVal("Coil_cryostat_half_z", half_z);

  //Side wall left

  gp->setDoubleVal("Coil_cryostat_side_l_inner_radius", 
                                        40*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_side_l_outer_radius", 
                                        -30*dd4hep::mm+outer_radius);
  gp->setDoubleVal("Coil_cryostat_side_l_half_z", 25.*dd4hep::mm);
  gp->setStringVal("Coil_material_side_l", "aluminium");

  //Side wall right

  gp->setDoubleVal("Coil_cryostat_side_r_inner_radius", 
                                        40*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_side_r_outer_radius", 
                                        -30*dd4hep::mm+outer_radius);
  gp->setDoubleVal("Coil_cryostat_side_r_half_z", 25.*dd4hep::mm);
  gp->setStringVal("Coil_material_side_r", "aluminium");

 // Coil modules

  gp->setDoubleVal("Coil_cryostat_modules_inner_radius", 
                                        175*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_modules_outer_radius", 
                                        560*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_modules_half_z", 1632.*dd4hep::mm/2-20.*dd4hep::mm);
  gp->setStringVal("Coil_material_modules", "aluminium");

  gp->setDoubleVal("Coil_cryostat_c_modules_inner_radius", 
                                        175*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_c_modules_outer_radius", 
                                        560*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_c_modules_half_z", 1224.*dd4hep::mm);
  gp->setStringVal("Coil_material_c_modules", "aluminium");


  //Coil mandrel


  gp->setDoubleVal("Coil_cryostat_mandrel_inner_radius", 
                                        560*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_mandrel_outer_radius", 
                                        635*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_mandrel_half_z", 3675.*dd4hep::mm);
  gp->setStringVal("Coil_material_mandrel", "aluminium");

  //Sensitive detectors

  gp->setDoubleVal("Coil_cryostat_scint1_inner_radius", 
                                        90*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_scint1_outer_radius", 
                                        100*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_scint1_zposin", 
                                        -100*dd4hep::mm+half_z);
  gp->setDoubleVal("Coil_cryostat_scint1_zposend", 
                                        +100*dd4hep::mm+half_z);
  gp->setStringVal("Coil_material_scint1", "polystyrene");

  gp->setDoubleVal("Coil_cryostat_scint2_inner_radius", 
                                        105*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_scint2_outer_radius", 
                                        115*dd4hep::mm+inner_radius);
  gp->setDoubleVal("Coil_cryostat_scint2_zposin", 
                                        -100*dd4hep::mm+half_z);
  gp->setDoubleVal("Coil_cryostat_scint2_zposend", 
                                        +100*dd4hep::mm+half_z);
  gp->setStringVal("Coil_material_scint2", "polystyrene");

  gp->setDoubleVal("Coil_cryostat_scint3_inner_radius", 
                                        -90*dd4hep::mm+outer_radius);
  gp->setDoubleVal("Coil_cryostat_scint3_outer_radius", 
                                        -80*dd4hep::mm+outer_radius);
  gp->setDoubleVal("Coil_cryostat_scint3_zposin", 
                                        -100*dd4hep::mm+half_z);
  gp->setDoubleVal("Coil_cryostat_scint3_zposend", 
                                        +100*dd4hep::mm+half_z);
  gp->setStringVal("Coil_material_scint3", "polystyrene");

  gp->setDoubleVal("Coil_cryostat_scint4_inner_radius", 
                                        -105*dd4hep::mm+outer_radius);
  gp->setDoubleVal("Coil_cryostat_scint4_outer_radius", 
                                        -95*dd4hep::mm+outer_radius);
  gp->setDoubleVal("Coil_cryostat_scint4_zposin", 
                                        -100*dd4hep::mm+half_z);
  gp->setDoubleVal("Coil_cryostat_scint4_zposend", 
                                        +100*dd4hep::mm+half_z);
  gp->setStringVal("Coil_material_scint4", "polystyrene");
  gearMgr->setGearParameters("CoilParameters", gp);
#endif


  cout << "Coil done.\n" << endl;








    // //====== create the meassurement surface ===================
    // Vector3D u,v,n ;

    // if( faces_IP == 0 ){
    //   // will be rotated around z-axis later
    //   u.fill(  0. ,  -1. , 0. ) ;
    //   v.fill(  0. ,   0. , 1. ) ;
    //   n.fill( -1. ,   0. , 0. ) ;

    //   // implement 7 deg stereo angle 
    //   u.fill( 0. , -cos( 3.5 * dd4hep::deg  ) , -sin( 3.5 * dd4hep::deg  ) ) ;
    //   v.fill( 0. , -sin( 3.5 * dd4hep::deg  ) ,  cos( 3.5 * dd4hep::deg  ) ) ;

    // } else {

    //   u.fill( 0. , 1. , 0. ) ;
    //   v.fill( 0. , 0. , 1. ) ;
    //   n.fill( 1. , 0. , 0. ) ;

    //   // implement 7 deg stereo angle 
    //   u.fill( 0. ,  cos( 3.5 * dd4hep::deg  ) ,  sin( 3.5 * dd4hep::deg  ) ) ;
    //   v.fill( 0. , -sin( 3.5 * dd4hep::deg  ) ,  cos( 3.5 * dd4hep::deg  ) ) ;
    // }


    // double inner_thick =  sensitive_thickness / 2.0 ;
    // double outer_thick =  sensitive_thickness / 2.0 + support_thickness ;  // support is on top
   
    // VolPlane surf( sitSenLogical , SurfaceType(SurfaceType::Sensitive,SurfaceType::Measurement1D) ,inner_thick, outer_thick , u,v,n ) ; //,o ) ;

    // // vector of sensor placements - needed for DetElements in ladder loop below
    // std::vector<PlacedVolume> pvV(  layer_geom.n_sensors_per_ladder ) ;

    // //============================================================



#endif //code_is_cleaned_up
  //=========================================================================================================


  cout << "SCoil02 done.\n" << endl;
  //######################################################################################################################################################################
  
  
  //--------------------------------------
  
  //coil.setVisAttributes( theDetector, x_det.visStr(), envelope );
  //added coded by Thorben Quast
  //the coil is modelled as a calorimeter layer to be consistent with the 
  //implementation of the solenoid (layers) of CLIC
  LayeredCalorimeterData* coilData = new LayeredCalorimeterData;

  //NN: Adding the rest of the data  
  coilData->inner_symmetry = 0;
  coilData->outer_symmetry = 0;
  coilData->layoutType = LayeredCalorimeterData::BarrelLayout ;
  
  coilData->extent[0] = inner_radius  ;
  coilData->extent[1] = outer_radius;
  coilData->extent[2] = 0. ;
  coilData->extent[3] = half_z;
  
  //NN: These probably need to be fixed and ced modified to read the extent, rather than the layer
  LayeredCalorimeterData::Layer coilLayer;
  coilLayer.distance = inner_radius;
  coilLayer.inner_thickness = ( outer_radius - inner_radius ) / 2. ;
  coilLayer.outer_thickness = coilLayer.inner_thickness  ;
  coilLayer.cellSize0 = 0;        //equivalent to 
  coilLayer.cellSize1 = half_z;    //half extension along z-axis
  coilData->layers.push_back(coilLayer);
  coil.addExtension< LayeredCalorimeterData >( coilData ) ;
  return coil;
}
DECLARE_DETELEMENT(SCoil02,create_element)

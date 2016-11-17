//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for SiWEcalEndcaps
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id: SEcal04_Endcaps.cpp 1060 2016-09-05 07:48:41Z /C=JP/O=KEK/OU=CRC/CN=JEANS Daniel Thomelin Dietrich $
//====================================================================

 /* History:  

/** SEcal05.cc
 *
 * new SEcal05 endcap driver: allows removal of preshower layer. DJeans UTokyo, sep/2016
    cleaned up, modularised code, added SEcal05_Helper, used by both barrel and endcap constructors
    add possibility to remove preshower layer
    can be used for either square cells or strips or a mixture
    knows nothing about the active material (ie can be silicon or scintillator or anything else)
    some changes to overall geometry in barrel (details of barrel overlap)
    split endcap quadrant into modules
    use MegatileLayerGridXY segmentation class
    added checks on parameters in compact description
    moved some hard-coded values to model parameters
    material of module fixed to be CarbonFiber

//
// $Id: SEcal04_Endcaps.cpp 1060 2016-09-05 07:48:41Z /C=JP/O=KEK/OU=CRC/CN=JEANS Daniel Thomelin Dietrich $
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
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/WaferGridXY.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

#include "SEcal05_Helpers.h"

#undef NDEBUG
#include <assert.h>


//#define VERBOSE 1

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(LCDD& lcdd, xml_h element, SensitiveDetector sens)  {

  cout << "------------------------" << endl;
  cout << "creating SEcal05_Endcaps" << endl;
  cout << "------------------------" << endl;

  xml_det_t     x_det     = element;
  string        det_name  = x_det.nameStr();
  Layering      layering (element);

  int           det_id    = x_det.id();
  //  xml_comp_t    x_staves  = x_det.staves();
  DetElement    sdet      (det_name,det_id);
  //  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = XML::createPlacedEnvelope( lcdd,  element , sdet ) ;

  XML::setDetectorTypeFlag( element, sdet ) ;

  if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;
  //-----------------------------------------------------------------------------------

  sens.setType("calorimeter");

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();

  //====================================================================
  //
  // Read all the constant from ILD_o1_v05.xml
  // Use them to build Ecal endcaps
  //
  //====================================================================
  
  // hardcoded numbers....better to add to xmk description
  const int N_FIBERS_W_STRUCTURE = 2; 
  const int N_FIBERS_ALVEOLUS = 3;

  //  read parameters from compact.xml file
  double Ecal_Slab_shielding                = lcdd.constant<double>("Ecal_Slab_shielding");
  double Ecal_fiber_thickness               = lcdd.constant<double>("Ecal_fiber_thickness");
  
  double EcalBarrel_inner_radius            = lcdd.constant<double>("TPC_outer_radius") +lcdd.constant<double>("Ecal_Tpc_gap");
  int    Ecal_barrel_z_modules              = lcdd.constant<int>("Ecal_barrel_z_modules");

  double Ecal_radiator_thickness1           = lcdd.constant<double>("Ecal_radiator_layers_set1_thickness");
  double Ecal_radiator_thickness2           = lcdd.constant<double>("Ecal_radiator_layers_set2_thickness");
  double Ecal_radiator_thickness3           = lcdd.constant<double>("Ecal_radiator_layers_set3_thickness");
  double Ecal_Barrel_halfZ                  = lcdd.constant<double>("Ecal_Barrel_halfZ");

  int    Ecal_barrel_number_of_towers       = lcdd.constant<int>("Ecal_barrel_number_of_towers");
  
  double Ecal_support_thickness             = lcdd.constant<double>("Ecal_support_thickness");
  double Ecal_front_face_thickness          = lcdd.constant<double>("Ecal_front_face_thickness");
  double Ecal_lateral_face_thickness        = lcdd.constant<double>("Ecal_lateral_face_thickness");
  double Ecal_Slab_H_fiber_thickness        = lcdd.constant<double>("Ecal_Slab_H_fiber_thickness");

  double Ecal_endcap_extra_size             = lcdd.constant<double>("Ecal_endcap_extra_size");
  double Ecal_cables_gap                    = lcdd.constant<double>("Ecal_cables_gap");

  int    Ecal_nlayers1                      = lcdd.constant<int>("Ecal_nlayers1");
  int    Ecal_nlayers2                      = lcdd.constant<int>("Ecal_nlayers2");
  int    Ecal_nlayers3                      = lcdd.constant<int>("Ecal_nlayers3");

  // first layer is preshower?
  bool   Ecal_Endcap_PreshowerLayer         = lcdd.constant<int>("Ecal_Endcap_Preshower") > 0;

  // this is a string like "341": first (central, longest) module has 3 towers, second one 4, third one 1
  // more usual would be e.g. "333" (for ILD large model)
  std::string Ecal_endcap_number_of_towers  = lcdd.constant<string>("Ecal_endcap_number_of_towers"); 

  int Ecal_end_of_slab_strategy             = lcdd.constant<int>("Ecal_end_of_slab_strategy");


  int    Ecal_n_wafers_per_tower            = lcdd.constant<int>("Ecal_n_wafers_per_tower");
  
  double Ecal_guard_ring_size               = lcdd.constant<double>("Ecal_guard_ring_size");

  double      EcalEndcap_inner_radius       = lcdd.constant<double>("EcalEndcap_inner_radius");
  double      EcalEndcap_outer_radius       = lcdd.constant<double>("EcalEndcap_outer_radius");
  double      EcalEndcap_min_z              = lcdd.constant<double>("EcalEndcap_min_z");
  double      EcalEndcap_max_z              = lcdd.constant<double>("EcalEndcap_max_z");

  std::string Ecal_layerConfig              = lcdd.constant<string>("Ecal_layer_pattern");

  int Ecal_cells_across_megatile            = lcdd.constant <int> ("Ecal_cells_across_megatile" );
  int Ecal_strips_across_megatile           = lcdd.constant <int> ("Ecal_strips_across_megatile");
  int Ecal_strips_along_megatile            = lcdd.constant <int> ("Ecal_strips_along_megatile" );

  double Ecal_barrel_thickness              = lcdd.constant<double>("Ecal_barrel_thickness"); // what's assumed in the compact description
 

  //====================================================================
  //
  // set up the helper
  //
  //====================================================================

  SEcal05_Helpers helper;

  helper.setDet( & x_det );

  // layer configuration
  helper.setPreshower( Ecal_Endcap_PreshowerLayer );
  if ( !Ecal_Endcap_PreshowerLayer && Ecal_front_face_thickness>0 ) {
    cout << "WARNING, using front CF plate in a non-preshower setup: do you really want this?" << endl;
  }

  cout << "Preshower ? " << Ecal_Endcap_PreshowerLayer << endl;

  helper.setAbsLayers( Ecal_nlayers1, Ecal_radiator_thickness1,
                       Ecal_nlayers2, Ecal_radiator_thickness2,
                       Ecal_nlayers3, Ecal_radiator_thickness3 );

  cout << "absorber layers " <<
    Ecal_nlayers1 << "*" << Ecal_radiator_thickness1 << "mm + " <<
    Ecal_nlayers2 << "*" << Ecal_radiator_thickness2 << "mm + " <<
    Ecal_nlayers3 << "*" << Ecal_radiator_thickness3 << "mm " << endl;

  // check this setup is self-consistent
  helper.checkLayerConsistency();

  // set structural CF thicknesses
  helper.setCFthickness( N_FIBERS_W_STRUCTURE*Ecal_fiber_thickness,
                         N_FIBERS_ALVEOLUS*Ecal_fiber_thickness,
                         Ecal_front_face_thickness,
                         Ecal_support_thickness);

  // check that resulting total thickness is consistent with what's in the compact file
  //   n.b. this assumes same thickness in barrel and endcap!
  float module_thickness = helper.getTotalThickness();
  if ( fabs( Ecal_barrel_thickness - module_thickness ) > 0.1*dd4hep::mm ) {
    cout << "ERROR : thickness in comapct decription not consistent with what I calculate!" << endl;
    cout << "    calculated = " << module_thickness << "    compact description: " << Ecal_barrel_thickness << endl;
    assert( 0 ); // exit
  }

  cout << "module thickness = " << module_thickness << endl;

  // set up the sensitive layer segmentation
  helper.setSegmentation( &seg );

  helper.setNCells( Ecal_cells_across_megatile, Ecal_strips_across_megatile, Ecal_strips_along_megatile);
  helper.setMagicMegatileStrategy ( Ecal_end_of_slab_strategy );



  std::vector < int > layerConfig;
  for(std::string::size_type i = 0; i < Ecal_layerConfig.size(); ++i) {
    char c = Ecal_layerConfig[i];
    int itype = atoi( &c );
    layerConfig.push_back( itype );
  }
  helper.setLayerConfig( layerConfig );

  cout << "layer config: ";
  for (size_t i=0; i<layerConfig.size(); i++) cout << layerConfig[i] << " ";
  cout << endl;

  // set the number of towers/modules


  std::vector < int > ntowers;
  for(std::string::size_type i = 0; i < Ecal_endcap_number_of_towers.size(); ++i) {
    char c = Ecal_endcap_number_of_towers[i];
    int itype = atoi( &c );
    ntowers.push_back( itype );
  }

  cout << "ECAL endcap tower configuration " << Ecal_endcap_number_of_towers << endl;
  for (size_t i=0; i<ntowers.size(); i++) {
    cout << ntowers[i] << " ";
  }
  cout << endl;

  // check that resulting quadrant size is consistent with compact description
  // here we assume that the width of an alveolus in the endcaps is the same as that in the barrel.
  double barrel_alv_width    = ( Ecal_Barrel_halfZ*2./Ecal_barrel_z_modules - 2.*Ecal_lateral_face_thickness ) / Ecal_barrel_number_of_towers;
  double calc_endcap_rout(EcalEndcap_inner_radius);
  for (size_t i=0; i<ntowers.size(); i++) {
    calc_endcap_rout+=ntowers[i]*barrel_alv_width + 2.*Ecal_lateral_face_thickness;
  }

  cout << "alveolus width = " << barrel_alv_width << endl;

  //  compare calculated size vs that in compact description
  if ( fabs( calc_endcap_rout - EcalEndcap_outer_radius ) > 0.1*dd4hep::mm ) {
    cout << "WARNING, inconsistent ECAL endcap radial extent! calculated: " << calc_endcap_rout << 
      " cm , nominal (from compact descrition): " << EcalEndcap_outer_radius << endl;
    cout << "consider changing Ecal_endcap_extra_size from " << Ecal_endcap_extra_size << 
      " to " << calc_endcap_rout - EcalBarrel_inner_radius - Ecal_barrel_thickness << endl;
    assert(0);
  }


  helper.setTowersUnits( ntowers,
			 barrel_alv_width,
			 Ecal_n_wafers_per_tower,
			 Ecal_lateral_face_thickness,
                         N_FIBERS_ALVEOLUS * Ecal_fiber_thickness + Ecal_Slab_H_fiber_thickness + Ecal_Slab_shielding,
			 Ecal_guard_ring_size );



  // ========= Create Ecal endcaps   ====================================
  //  It will be the volume for palcing the Ecal endcaps alveolus(i.e. Layers).
  //  And the structure W plate.
  //  Itself will be placed into the world volume.
  // ==========================================================================

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~              EndcapStandardModule                 ~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // octagon
  PolyhedraRegular ECPolyHedra(8, M_PI/8., 0., calc_endcap_rout, module_thickness);

  // take just one quadrant of the octagon, and remove the centre box
  double quadr = calc_endcap_rout; // anything which is big enough
  Box quadrant( quadr, quadr , module_thickness);
  IntersectionSolid EndCapSolid( ECPolyHedra, 
				 quadrant,
				 Position( quadr - EcalEndcap_inner_radius, 
					   quadr + EcalEndcap_inner_radius, 0 ) );

  Volume EnvLogEndCap("EcalEndcapQuadrant",EndCapSolid,lcdd.material("CarbonFiber"));

  // kink position wrt bottom of quadrant (inside the lateral support)
  // Y==0 is defined as outer edge of lateral face of inner module of quadrant
  double kink_y = calc_endcap_rout*tan(M_PI/8.) - EcalEndcap_inner_radius;

  helper.setModuleDimensions(1, // xytype
			     0, // xztype
			     calc_endcap_rout + EcalEndcap_inner_radius, // dxMax
			     kink_y
			     );
  
  helper.setTranslation ( DD4hep::Geometry::Position ( -EcalEndcap_inner_radius , EcalEndcap_inner_radius, -module_thickness/2. ) );

  // make the module

  DDRec::LayeredCalorimeterData* caloData = new DDRec::LayeredCalorimeterData ;

  DetElement mod_det ("quad0",det_id);

  helper.makeModule( EnvLogEndCap, 
		     mod_det,
		     *caloData,
		     lcdd, 
		     sens );

  for (size_t i=0; i<caloData->layers.size(); i++) {
    caloData->layers[i].distance += Ecal_Barrel_halfZ + Ecal_cables_gap; // add IP->front face distance
  }

  cout << "cell sizes: " << endl;
  for (size_t i=0; i<caloData->layers.size(); i++) {
    cout << "sensitive layer " << i << " : x,y = " << caloData->layers[i].cellSize0 << " " << caloData->layers[i].cellSize1 << endl;
  }

  //====================================================================
  // Place Ecal Endcap modules into the assembly envelope volume
  //====================================================================
   
  for (int iend=0; iend<2; iend++) { // the 2 endcaps
    int module_id = ( iend == 0 ) ? 0:6;

    double this_module_z_offset = (EcalEndcap_min_z + EcalEndcap_max_z)/2.;
    if ( iend == 0 ) this_module_z_offset*=-1;

    double this_module_rotY = iend==0 ? M_PI : 0;
    double rotZ_offset = iend==0 ? M_PI/8. - M_PI/2. : M_PI/8. + 3.*M_PI/4.; // magic rotation to get modules in right place
    for (int iquad=0; iquad<4; iquad++) { // the 4 quadrants per endcap
      int stave_id=iquad+1;
      double this_module_rotZ(0);
      if ( iend==0 ) {
	this_module_rotZ = rotZ_offset - (iquad-2) * M_PI/2.;
      } else {
	this_module_rotZ = rotZ_offset + (iquad+1) * M_PI/2.;
      }
      Position xyzVec(0,0,this_module_z_offset);
      RotationZYX rot( this_module_rotZ ,this_module_rotY, 0);
      Rotation3D rot3D(rot);
      Transform3D tran3D(rot3D,xyzVec);
      PlacedVolume pv = envelope.placeVolume(EnvLogEndCap,tran3D);
      pv.addPhysVolID("module",module_id); // z: -/+ 0/6
      pv.addPhysVolID("stave",stave_id);
      DetElement sd = (iend==0 && iquad==0) ? mod_det : mod_det.clone(_toString(module_id,"module%d")+_toString(stave_id,"stave%d"));
      sd.setPlacement(pv);
    }
  }
  
  sdet.addExtension< DDRec::LayeredCalorimeterData >( caloData ) ; 

  //  cout << "finished SEcal05_Endcaps" << endl;

  return sdet;
  
}



DECLARE_DETELEMENT(SEcal05_Endcaps, create_detector)


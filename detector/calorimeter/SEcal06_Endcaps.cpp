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
//#include "DDSegmentation/WaferGridXY.h"
#include <sstream>

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::IntersectionSolid;
using dd4hep::Layering;
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
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

#include "SEcal06_Helpers.h"

#undef NDEBUG
#include <assert.h>


//#define VERBOSE 1

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {

  cout << "------------------------" << endl;
  cout << "creating SEcal06_Endcaps" << endl;
  cout << "------------------------" << endl;

  xml_det_t     x_det     = element;
  string        det_name  = x_det.nameStr();
  Layering      layering (element);

  int           det_id    = x_det.id();
  //  xml_comp_t    x_staves  = x_det.staves();
  DetElement    sdet      (det_name,det_id);
  //  Volume        motherVol = theDetector.pickMotherVolume(sdet);

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;

  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
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
  
  //  read parameters from compact.xml file
  double Ecal_Slab_shielding                = theDetector.constant<double>("Ecal_Slab_shielding");
  double Ecal_fiber_thickness_structure = theDetector.constant<double>("Ecal_fiber_thickness_structure"); // absorber wrapping thickness
  double Ecal_fiber_thickness_alveolus  = theDetector.constant<double>("Ecal_fiber_thickness_alveolus"); // alveolar wall thickness
  
  double EcalBarrel_inner_radius                  = theDetector.constant<double>("TPC_outer_radius") +theDetector.constant<double>("Ecal_Tpc_gap");
  int    Ecal_barrel_z_modules              = theDetector.constant<int>("Ecal_barrel_z_modules");

  double Ecal_radiator_thickness1           = theDetector.constant<double>("Ecal_radiator_layers_set1_thickness");
  double Ecal_radiator_thickness2           = theDetector.constant<double>("Ecal_radiator_layers_set2_thickness");
  double Ecal_radiator_thickness3           = theDetector.constant<double>("Ecal_radiator_layers_set3_thickness");
  double Ecal_Barrel_halfZ                  = theDetector.constant<double>("Ecal_Barrel_halfZ");

  int    Ecal_barrel_number_of_towers       = theDetector.constant<int>("Ecal_barrel_number_of_towers");
  
  double Ecal_support_thickness             = theDetector.constant<double>("Ecal_support_thickness");
  double Ecal_front_face_thickness          = theDetector.constant<double>("Ecal_front_face_thickness");
  double Ecal_lateral_face_thickness        = theDetector.constant<double>("Ecal_lateral_face_thickness");
  double Ecal_Slab_H_fiber_thickness        = theDetector.constant<double>("Ecal_Slab_H_fiber_thickness");

  double Ecal_endcap_extra_size             = theDetector.constant<double>("Ecal_endcap_extra_size");
  double Ecal_cables_gap                    = theDetector.constant<double>("Ecal_cables_gap");

  int    Ecal_nlayers1                      = theDetector.constant<int>("Ecal_nlayers1");
  int    Ecal_nlayers2                      = theDetector.constant<int>("Ecal_nlayers2");
  int    Ecal_nlayers3                      = theDetector.constant<int>("Ecal_nlayers3");

  // first layer is preshower?
  bool   Ecal_Endcap_PreshowerLayer         = theDetector.constant<int>("Ecal_Endcap_Preshower") > 0;

  std::string Ecal_endcap_number_of_towers       = theDetector.constant<string>("Ecal_endcap_number_of_towers");

  int Ecal_end_of_slab_strategy             = theDetector.constant<int>("Ecal_end_of_slab_strategy");


  int    Ecal_n_wafers_per_tower            = theDetector.constant<int>("Ecal_n_wafers_per_tower");
  
  double Ecal_guard_ring_size               = theDetector.constant<double>("Ecal_guard_ring_size");

  double      EcalEndcap_inner_radius          = theDetector.constant<double>("EcalEndcap_inner_radius");
  double      EcalEndcap_outer_radius          = theDetector.constant<double>("EcalEndcap_outer_radius");
  double      EcalEndcap_min_z                 = theDetector.constant<double>("EcalEndcap_min_z");
  double      EcalEndcap_max_z                 = theDetector.constant<double>("EcalEndcap_max_z");

  std::string Ecal_layerConfig              = theDetector.constant<string>("Ecal_layer_pattern");

  int Ecal_cells_across_megatile = theDetector.constant <int> ("Ecal_cells_across_megatile" );
  int Ecal_strips_across_megatile= theDetector.constant <int> ("Ecal_strips_across_megatile");
  int Ecal_strips_along_megatile = theDetector.constant <int> ("Ecal_strips_along_megatile" );

  double Ecal_barrel_thickness              = theDetector.constant<double>("Ecal_barrel_thickness"); // what's assumed in the compact description
 

  //====================================================================
  //
  // set up the helper
  //
  //====================================================================

  SEcal06_Helpers helper;

  helper.setDet( & x_det );

  // layer configuration
  helper.setPreshower( Ecal_Endcap_PreshowerLayer );

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
  helper.setCFthickness( Ecal_fiber_thickness_structure,
                         Ecal_fiber_thickness_alveolus,
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


  int ntemp;
  std::stringstream stream(Ecal_layerConfig);
  std::vector < int > layerConfig;
  while ( stream >> ntemp ) {
    assert( ntemp>=0 );
    layerConfig.push_back( ntemp );
  }
  cout << "layer config: ";
  for (size_t i=0; i<layerConfig.size(); i++) cout << layerConfig[i] << " ";
  cout << endl;

  helper.setLayerConfig( layerConfig );

  // set the number of towers/modules
  std::vector < int > ntowers;
  std::stringstream stream2( Ecal_endcap_number_of_towers );
  while ( stream2 >> ntemp ) {
    ntowers.push_back( ntemp );
  }

  cout << "ECAL endcap tower configuration " << Ecal_endcap_number_of_towers << " : ";
  for (size_t i=0; i<ntowers.size(); i++) 
    cout << ntowers[i] << " ";
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
                         Ecal_fiber_thickness_alveolus + Ecal_Slab_H_fiber_thickness + Ecal_Slab_shielding,
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

  Volume EnvLogEndCap("EcalEndcapQuadrant",EndCapSolid,theDetector.material("CarbonFiber"));

  // kink position wrt bottom of quadrant (inside the lateral support)
  // Y==0 is defined as outer edge of lateral face of inner module of quadrant
  double kink_y = calc_endcap_rout*tan(M_PI/8.) - EcalEndcap_inner_radius;

  helper.setModuleDimensions(1, // xytype
			     0, // xztype
			     calc_endcap_rout + EcalEndcap_inner_radius, // dxMax
			     kink_y
			     );
  
  helper.setTranslation ( Position ( -EcalEndcap_inner_radius , EcalEndcap_inner_radius, -module_thickness/2. ) );

  // make the module

  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;

  DetElement mod_det ("quad0",det_id);

  helper.makeModule( EnvLogEndCap, 
		     mod_det,
		     *caloData,
		     theDetector, 
		     sens );

  for (size_t i=0; i<caloData->layers.size(); i++) {
    caloData->layers[i].distance += Ecal_Barrel_halfZ + Ecal_cables_gap; // add IP->front face distance
  }

  //  cout << "cell sizes: " << endl;
  //  for (size_t i=0; i<caloData->layers.size(); i++) {
  //    cout << "sensitive layer " << i << " : x,y = " << caloData->layers[i].cellSize0 << " " << caloData->layers[i].cellSize1 << endl;
  //  }
  
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = 4  ; // hard code cernter box hole
  caloData->outer_symmetry = 8  ; // outer Octagun
  caloData->phi0 = 0 ;

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = EcalEndcap_inner_radius ;
  caloData->extent[1] = EcalEndcap_outer_radius ;
  caloData->extent[2] = EcalEndcap_min_z ;
  caloData->extent[3] = EcalEndcap_max_z ;



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
  
  sdet.addExtension< LayeredCalorimeterData >( caloData ) ; 

  //  cout << "finished SEcal06_Endcaps" << endl;

  return sdet;
  
}



DECLARE_DETELEMENT(SEcal06_Endcaps, create_detector)


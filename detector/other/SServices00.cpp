//====================================================================
//  DDSim - LC detector models in DD4hep 
//--------------------------------------------------------------------
// DD4hep Geometry driver for SService00
// Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id$
//====================================================================
//
//*******************************************************
//*                                                     *
//*                      Mokka                          * 
//*   - the detailed geant4 simulation for ILC   -      *
//*                                                     *
//* For more information about Mokka, visit the         *
//*                                                     *
//*  Mokka.in2p3.fr  Mokka home page.                   *
//*                                                     *
//*******************************************************
//
// $Id$
// $Name: mokka-07-00 $
//
// 
//
// SServices00.cc
//

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "XML/Utilities.h"
#include "DD4hep/DetType.h"
#include "XMLHandlerDB.h"

#include "SServices00.h"
 
using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Assembly;
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


static Ref_t create_element(Detector& theDetector, xml_h element, Ref_t)  {

  //unused:  static double tolerance = 0e0;

  xml_det_t x_det = element;
  string det_name = x_det.nameStr();
  
  int det_id = x_det.id();
  DetElement sdet (det_name,det_id);
  Volume motherVol = theDetector.pickMotherVolume(sdet);
  
  Assembly envelope_assembly( det_name + "assembly" ) ;
  PlacedVolume pv;

  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;

//====================================================================
// build all services
//====================================================================
  cout << "\nBuilding SServices00"<< endl;

  //==================================================
  //           BuildTPCEndplateServices
  //==================================================
  
  // start to prepare the Material and geometry as Mokka
  double TPC_Ecal_Hcal_barrel_halfZ =theDetector.constant<double>("TPC_Ecal_Hcal_barrel_halfZ");

  XMLHandlerDB db = XMLHandlerDB(  x_det.child( _Unicode( TPC_Cooling ) ) ) ;

  const int N_TPC_RINGS =  db->fetchInt("number_of_rings");;
  //fg: keeep gcc on SL5 happy:
  static const int MAX_TPC_RINGS = 1024 ;

  double tpcEndplateServices_R[MAX_TPC_RINGS];
  double tpcEndplateServices_r[MAX_TPC_RINGS];

  for(xml_coll_t c( x_det ,_U(ring)); c; ++c)  {
    
    xml_comp_t  x_ring( c );
    db = XMLHandlerDB( x_ring )  ;
    
    int ring_id = db->fetchInt("ring_id");
    cout <<" ring_id : " << ring_id <<endl;

    tpcEndplateServices_R[ring_id] = db->fetchDouble("tpcEndplateServicesRing_R");
    tpcEndplateServices_r[ring_id] = db->fetchDouble("tpcEndplateServicesRing_ro");

  }
 
  // start to build TPC cooling rings with the class BuildTPCEndplateServices
  BuildTPCEndplateServices TPCEndplateServices;
  //TPCEndplateServices.setMaterial(theDetector.material("copper"));
  TPCEndplateServices.setMaterial(theDetector.material("Cu"));
  TPCEndplateServices.sethalfZ(TPC_Ecal_Hcal_barrel_halfZ);
  for(int i=0;i<N_TPC_RINGS;i++)
    TPCEndplateServices.settpcEndplateServicesRing_R_ro(tpcEndplateServices_R[i],tpcEndplateServices_r[i]);
  
  TPCEndplateServices.DoBuildTPCEndplateServices(pv,envelope_assembly);
 


 
  //==================================================
  //           BuildEcalBarrelServices
  //==================================================

  // start to prepare the Material and geometry as Mokka
  double Ecal_outer_radius = theDetector.constant<double>("Ecal_outer_radius");
  double Ecal_inner_radius                  = theDetector.constant<double>("TPC_outer_radius") +theDetector.constant<double>("Ecal_Tpc_gap");
  double module_thickness = Ecal_outer_radius - Ecal_inner_radius;
 
  // module barrel key parameters
  double bottom_dim_x = 2. * tan(M_PI/8.) * Ecal_inner_radius +
    module_thickness/sin(M_PI/4.);
  
  double top_dim_x = bottom_dim_x - 2 * module_thickness;
  
  double RailHeight = theDetector.constant<double>("EcalBarrelServices_RailHeight");

  cout<<"\n Ecal_inner_radius: " << Ecal_inner_radius
      <<"\n Ecal_outer_radius: " << Ecal_outer_radius
      <<"\n  module_thickness: " << module_thickness 
      <<"\n      bottom_dim_x: " << bottom_dim_x
      <<"\n         top_dim_x: " << top_dim_x
      <<"\n        RailHeight: " << RailHeight
      <<"\n" <<endl;


  double RailSeparation = theDetector.constant<double>("EcalBarrelServices_RailSeparation");
  double RailWidth =  theDetector.constant<double>("EcalBarrelServices_RailWidth");

  double ZMinus_FirstInterrail_PE_Thickness = theDetector.constant<double>("EcalBarrelServices_ZMinus_FirstInterrail_PE_Thickness");
  double ZMinus_FirstInterrail_Cu_Thickness = theDetector.constant<double>("EcalBarrelServices_ZMinus_FirstInterrail_Cu_Thickness");

  double ZPlus_FirstInterrail_PE_Thickness = theDetector.constant<double>("EcalBarrelServices_ZPlus_FirstInterrail_PE_Thickness");
  double ZPlus_FirstInterrail_Cu_Thickness = theDetector.constant<double>("EcalBarrelServices_ZPlus_FirstInterrail_Cu_Thickness");

  // start to build Ecal Barrel services with the class BuildEcalBarrelServices
  BuildEcalBarrelServices EcalBarrelServices;

  EcalBarrelServices.setMaterialAir(theDetector.air());
  EcalBarrelServices.setMaterialAluminium(theDetector.material("Al"));
  EcalBarrelServices.setMaterialPolyethylene(theDetector.material("G4_POLYSTYRENE"));
  EcalBarrelServices.setMaterialCopper(theDetector.material("Cu"));
  EcalBarrelServices.setenv_safety( theDetector.constant<double>("env_safety"));

  EcalBarrelServices.sethalfZ(TPC_Ecal_Hcal_barrel_halfZ);
  EcalBarrelServices.setTopDimX(top_dim_x);
  EcalBarrelServices.setRailHeight(RailHeight);
  EcalBarrelServices.setOutRadius(Ecal_outer_radius);
  EcalBarrelServices.setModuleThickness(module_thickness);

  EcalBarrelServices.setRailSeparation(RailSeparation);
  EcalBarrelServices.setRailWidth(RailWidth);

  EcalBarrelServices.setZMinus_FirstInterrail_PE_Thickness(ZMinus_FirstInterrail_PE_Thickness);
  EcalBarrelServices.setZMinus_FirstInterrail_Cu_Thicknesst(ZMinus_FirstInterrail_Cu_Thickness);

  EcalBarrelServices.setZPlus_FirstInterrail_PE_Thickness(ZPlus_FirstInterrail_PE_Thickness);
  EcalBarrelServices.setZPlus_FirstInterrail_Cu_Thickness(ZPlus_FirstInterrail_Cu_Thickness);

  EcalBarrelServices.DoBuildEcalBarrelServices(pv,envelope_assembly);




  //==================================================
  //          BuildEcalBarrel_EndCapServices
  //==================================================

  // start to prepare the Material and geometry as Mokka
  double ZMinus_PE_Thickness = theDetector.constant<double>("EcalBarrel_EndCapServices_ZMinus_PE_Thickness");
  double ZMinus_Cu_Thickness = theDetector.constant<double>("EcalBarrel_EndCapServices_ZMinus_Cu_Thickness");
  double ZPlus_PE_Thickness = theDetector.constant<double>("EcalBarrel_EndCapServices_ZPlus_PE_Thickness");
  double ZPlus_Cu_Thickness = theDetector.constant<double>("EcalBarrel_EndCapServices_ZPlus_Cu_Thickness");
  double Ecal_cables_gap = theDetector.constant<double>("Ecal_cables_gap");
  double InnerServicesWidth = theDetector.constant<double>("HcalServicesModule_InnerServicesWidth");

  BuildEcalBarrel_EndCapServices EcalBarrel_EndCapServices;
  EcalBarrel_EndCapServices.setMaterialAir(theDetector.air());
  EcalBarrel_EndCapServices.setMaterialPolyethylene(theDetector.material("G4_POLYSTYRENE"));
  EcalBarrel_EndCapServices.setMaterialCopper(theDetector.material("Cu"));
  EcalBarrel_EndCapServices.sethalfZ(TPC_Ecal_Hcal_barrel_halfZ);
  EcalBarrel_EndCapServices.setTopDimX(top_dim_x);
  EcalBarrel_EndCapServices.setOutRadius(Ecal_outer_radius);
  EcalBarrel_EndCapServices.setModuleThickness(module_thickness);
  EcalBarrel_EndCapServices.setInnerServicesWidth(InnerServicesWidth);
  EcalBarrel_EndCapServices.setEcal_cables_gap(Ecal_cables_gap);
  EcalBarrel_EndCapServices.setZMinus_PE_Thickness(ZMinus_PE_Thickness);
  EcalBarrel_EndCapServices.setZMinus_Cu_Thickness(ZMinus_Cu_Thickness);
  EcalBarrel_EndCapServices.setZPlus_PE_Thickness(ZPlus_PE_Thickness);
  EcalBarrel_EndCapServices.setZPlus_Cu_Thickness(ZPlus_Cu_Thickness);
  EcalBarrel_EndCapServices.setenv_safety( theDetector.constant<double>("env_safety"));

  EcalBarrel_EndCapServices.DoBuildEcalBarrel_EndCapServices(pv,envelope_assembly);





  //==================================================
  //         BuildHcalBarrel_EndCapServices
  //==================================================

  // start to prepare the Material and geometry as Mokka
  BuildHcalBarrel_EndCapServices HcalBarrel_EndCapServices;
  HcalBarrel_EndCapServices.setMaterialAir( theDetector.air() );
  HcalBarrel_EndCapServices.setMaterialPolyethylene(theDetector.material("G4_POLYSTYRENE"));
  HcalBarrel_EndCapServices.setMaterialCopper(theDetector.material("Cu"));
  HcalBarrel_EndCapServices.setMaterialstainless_steel(theDetector.material("stainless_steel"));
  HcalBarrel_EndCapServices.setMaterialS235(theDetector.material("Steel235"));
  HcalBarrel_EndCapServices.setMaterialPCB(theDetector.material("PCB"));

  HcalBarrel_EndCapServices.setHcal_stave_gaps( theDetector.constant<double>("Hcal_stave_gaps") );
  HcalBarrel_EndCapServices.setHcal_inner_radius( theDetector.constant<double>("Hcal_inner_radius") );
  HcalBarrel_EndCapServices.setHcal_R_max( theDetector.constant<double>("Hcal_outer_radius") );
  HcalBarrel_EndCapServices.setBuildHcalElectronicsInterface( theDetector.constant<int>("BuildHcalElectronicsInterface") );
  HcalBarrel_EndCapServices.sethalfZ( theDetector.constant<double>("TPC_Ecal_Hcal_barrel_halfZ") );
  HcalBarrel_EndCapServices.setEcal_cables_gap( theDetector.constant<double>("Ecal_cables_gap") );

//-Z thicknesses
  HcalBarrel_EndCapServices.setZMinus_StainlessSteel_Thickness( theDetector.constant<double>("HcalServicesModule_ZMinus_StainlessSteel_Thickness") );
  HcalBarrel_EndCapServices.setZMinus_PE_Thickness( theDetector.constant<double>("HcalServicesModule_ZMinus_PE_Thickness") );
  HcalBarrel_EndCapServices.setZMinus_Cu_Thickness( theDetector.constant<double>("HcalServicesModule_ZMinus_Cu_Thickness") );
//+Z thicknesses
  HcalBarrel_EndCapServices.setZPlus_StainlessSteel_Thickness( theDetector.constant<double>("HcalServicesModule_ZPlus_StainlessSteel_Thickness") );
  HcalBarrel_EndCapServices.setZPlus_PE_Thickness( theDetector.constant<double>("HcalServicesModule_ZPlus_PE_Thickness") );
  HcalBarrel_EndCapServices.setZPlus_Cu_Thickness( theDetector.constant<double>("HcalServicesModule_ZPlus_Cu_Thickness") );

  HcalBarrel_EndCapServices.setInnerServicesWidth( theDetector.constant<double>("HcalServicesModule_InnerServicesWidth") );
  HcalBarrel_EndCapServices.setHcal_back_plate_thickness( theDetector.constant<double>("Hcal_back_plate_thickness") );
  HcalBarrel_EndCapServices.setHcal_nlayers( theDetector.constant<int>("Hcal_nlayers") );
  HcalBarrel_EndCapServices.setHcal_radiator_thickness( theDetector.constant<double>("Hcal_radiator_thickness") );

  HcalBarrel_EndCapServices.setHcal_steel_cassette_thickness( theDetector.constant<double>("Hcal_steel_cassette_thickness") );
  HcalBarrel_EndCapServices.setHcalServices_outer_FR4_thickness( theDetector.constant<double>("HcalServices_outer_FR4_thickness") );
  HcalBarrel_EndCapServices.setHcalServices_outer_Cu_thickness( theDetector.constant<double>("HcalServices_outer_Cu_thickness") );
  HcalBarrel_EndCapServices.setenv_safety( theDetector.constant<double>("env_safety"));

  HcalBarrel_EndCapServices.DoBuildHcalBarrel_EndCapServices(pv,envelope_assembly);





  //==================================================
  //          BuildSitCables
  //==================================================

  // start to prepare the Material and geometry as Mokka

  BuildSitCables SitCables;
  SitCables.setMaterialAluminium( theDetector.material("Al") );
  SitCables.setSit_cables_cylinder_thickness( theDetector.constant<double>("SIT12_cable_thickness") );
  SitCables.setftd4to7_tpc_radial_gap(  theDetector.constant<double>("ftd4to7_tpc_radial_gap") );
  SitCables.sethalfZ( TPC_Ecal_Hcal_barrel_halfZ );
  SitCables.setTPC_Inner_Radius( theDetector.constant<double>("TPC_inner_radius") );
  SitCables.setz2_position_ReltoTPCLength( theDetector.constant<double>("FTD_disk2_zPosRelToTpcLength") );
  SitCables.setz3_position_ReltoTPCLength( theDetector.constant<double>("FTD_disk3_zPosRelToTpcLength") );
  SitCables.setpetal_cp_support_thickness( theDetector.constant<double>("petal_cp_support_thickness") );
  SitCables.setdisk_si_thickness( theDetector.constant<double>("disk_si_thickness") );
  SitCables.setpetal_support_zoffset( theDetector.constant<double>("petal_support_zoffset") );
  SitCables.setftd2_sit1_radial_diff( theDetector.constant<double>("ftd2_sit1_radial_diff") );
  SitCables.setftd3_sit2_radial_diff( theDetector.constant<double>("ftd3_sit2_radial_diff") );
  SitCables.setSIT1_Radius( theDetector.constant<double>("SIT1_Radius") );
  SitCables.setSIT2_Radius( theDetector.constant<double>("SIT2_Radius") );
  SitCables.setFTD2_cone_thickness( theDetector.constant<double>("SServices_FTD2_cone_thickness") );
  SitCables.setFTD3_cone_thickness( theDetector.constant<double>("SServices_FTD3_cone_thickness") );
  SitCables.setTUBE_IPOuterBulge_end_radius( theDetector.constant<double>("TUBE_IPOuterBulge_end_radius") );
  SitCables.setTUBE_IPOuterBulge_end_envradius( theDetector.constant<double>("TUBE_IPOuterBulge_end_envradius") );
  SitCables.setTUBE_IPOuterBulge_end_z( theDetector.constant<double>("TUBE_IPOuterBulge_end_z") );
  SitCables.setSServices_FTD7_cables_thickness( theDetector.constant<double>("SServices_FTD7_cables_thickness") );

  SitCables.DoBuildSitCables(pv,envelope_assembly);


  //==================================================
  //          BuildVXDCables
  //==================================================

  BuildVXDCables VXDCables;
  VXDCables.setMaterialCopper ( theDetector.material("Cu") );
  VXDCables.setVXD_cable_cross_section_area( theDetector.constant<double>("VXD_cable_cross_section_area") );
  VXDCables.setVXD_cable_z_start( theDetector.constant<double>("FTD_min_z_0") );
  VXDCables.setVXD_cable_z_end( theDetector.constant<double>("TUBE_IPOuterTube_end_z") ) ;
  VXDCables.setVXD_cable_inner1_radius( theDetector.constant<double>("VXD_cable_inner1_radius") );
  VXDCables.setVXD_cable_inner2_radius( theDetector.constant<double>("TUBE_IPOuterBulge_start_envradius") );

  VXDCables.DoBuildVXDCables(pv,envelope_assembly);


//====================================================================
// Place services into the world volume
//====================================================================

  pv = motherVol.placeVolume(envelope_assembly);
  pv.addPhysVolID("system", det_id);
  sdet.setVisAttributes( theDetector, x_det.visStr(),  envelope_assembly);
  sdet.setPlacement(pv);

  return sdet;
}

DECLARE_DETELEMENT(SServices00,create_element)

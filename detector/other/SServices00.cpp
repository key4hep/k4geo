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
#include "DDRec/Surface.h"
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

  EcalBarrelServices.setenv_safety( theDetector.constant<double>("env_safety"));

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

  BuildSitCables SitCables( sdet ) ;
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

  BuildVXDCables VXDCables( sdet );
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






// to build TPC cooling rings into the service assembly 
bool BuildTPCEndplateServices::DoBuildTPCEndplateServices(dd4hep::PlacedVolume &pVol,dd4hep::Assembly &envelope){
  for( std::vector< std::pair<double,double> >::const_iterator it = tpcEndplateServicesRing_R_ro.begin() ;
       it != tpcEndplateServicesRing_R_ro.end(); ++it){
    dd4hep::Torus solidTube(it->first, 0, it->second, 0, 2*M_PI);
    double z_position = TPC_barrel_halfZ + it->second;
    dd4hep::Volume EnvTube("service_tube",solidTube,cooling_Material);
    //EnvTube.setVisAttributes("RedVis");
    dd4hep::Position pos1(0,0,z_position);
    pVol = envelope.placeVolume(EnvTube,pos1);
    dd4hep::Position pos2(0,0,-z_position);
    pVol = envelope.placeVolume(EnvTube,pos2);
  } 
  return true;
}


/// to fill the detail service layers into the container
bool BuildEcalBarrelServices::FillEcalBarrelServicesContainer(dd4hep::PlacedVolume &pVol,
							      dd4hep::Volume &pContainerLogical){
  
  const int NRAILS=2;
  if ( RailSeparation*(NRAILS-1) > top_dim_x - RailWidth ) {
    std::cout << "WARNING, requested rail separation too large! " << std::endl;
    std::cout << "Must be smaller than " << (top_dim_x - RailWidth)/(NRAILS-1) << " cm" << std::endl;
    assert(0);
  }
  
  float railPositions[NRAILS]; // positions wrt module centre: equally spaced, centered on module
  for (int i=0; i<NRAILS; i++) {
    railPositions[i] = i*RailSeparation - RailSeparation*(NRAILS-1)/2.;
  }
  
  dd4hep::Box railSolid(RailWidth/2., RailHeight/2., Ecal_barrel_halfZ); 
  
  dd4hep::Volume railLogical("railLogical", railSolid, aluminium);
  
  for(int i=0; i<NRAILS; i++)
  {
    dd4hep::Position pos(railPositions[i],0,0);
    pVol = pContainerLogical.placeVolume(railLogical,pos);
  }
  
  double moduleLength = Ecal_barrel_halfZ * 2. / 5.;
  
  double internalZoneWidth = RailSeparation - RailWidth - 0.1; // space between rails (in x-dir) (minus a small gap)
  
  // now fit fit service materials into the zone between rails
  
  const int NMODULES_IN_Z=5;
  for(int i=0; i<NMODULES_IN_Z; i++){ // these are the 5 modules in z
    
    float pe_thick = ZMinus_FirstInterrail_PE_Thickness;
    float cu_thick = ZMinus_FirstInterrail_Cu_Thickness;
    if ( i>2 ) {
      pe_thick = ZPlus_FirstInterrail_PE_Thickness;
      cu_thick = ZPlus_FirstInterrail_Cu_Thickness;
    }
    
    // polyethylene
    dd4hep::Box PESolid( internalZoneWidth/2. , pe_thick / 2., moduleLength / 2.); 
    
    std::string PELogical_name  = "PELogical"+dd4hep::_toString(i,"_%d");
    dd4hep::Volume PELogical(PELogical_name,PESolid,polyethylene);
    //PELogical.setVisAttributes("GreenVis");
    
    dd4hep::Position posPE( 0 , 
                            -RailHeight/2. + pe_thick/2.,
                            -Ecal_barrel_halfZ + moduleLength*(i+1/2.));
    
    pVol = pContainerLogical.placeVolume(PELogical,posPE);
    
    // Cu_1
    dd4hep::Box Cu_1_Solid( internalZoneWidth/2., cu_thick / 2., moduleLength / 2.); 
    
    std::string Cu_1_Logical_name  = "Cu_1_Logical"+dd4hep::_toString(i,"_%d");
    dd4hep::Volume Cu_1_Logical(Cu_1_Logical_name, Cu_1_Solid,copper);
    //Cu_1_Logical.setVisAttributes("BlueVis");
    
    dd4hep::Position posCu1( 0 , 
                             -RailHeight/2. + pe_thick + cu_thick/2.,
                             -Ecal_barrel_halfZ + moduleLength*(i+1/2.));
    
    pVol = pContainerLogical.placeVolume(Cu_1_Logical,posCu1);
    
  }

  return true;
}


// to build Ecal Barrel service into the service assembly 
bool BuildEcalBarrelServices::DoBuildEcalBarrelServices(dd4hep::PlacedVolume &pVol,dd4hep::Assembly &envelope){
  
  dd4hep::Box ContainerSolid(top_dim_x/2., RailHeight/2., Ecal_barrel_halfZ); 
  dd4hep::Volume EBScontainerLogical("EcalBarrelServicesContainerLogical",ContainerSolid,air);
  //containerLogical.setVisAttributes("SeeThrough");
  //containerLogical.setVisAttributes("MagentaVis");
  
  if(!FillEcalBarrelServicesContainer(pVol,EBScontainerLogical))
    return false;
  
  for (int stave_id = 1; stave_id < 9 ; stave_id++)
  {
    double phirot = (stave_id-1) * M_PI/4.;
    dd4hep::RotationZYX rot(-phirot,0,0);
    dd4hep::Position stavePosition((Ecal_outer_radius + env_safety*2.0 + RailHeight/2.)*sin(phirot) - module_thickness * sin(M_PI/4.)*cos(-phirot),
				   (Ecal_outer_radius + env_safety*2.0 + RailHeight/2.)*cos(phirot) - module_thickness * sin(M_PI/4.)*sin(-phirot),
				   0 );
    dd4hep::Transform3D tran3D(rot,stavePosition);
    pVol = envelope.placeVolume(EBScontainerLogical,tran3D);
  }
  
  return true;
}


// to build Ecal Barrel service into the service assembly 
bool BuildEcalBarrel_EndCapServices::DoBuildEcalBarrel_EndCapServices(dd4hep::PlacedVolume &pVol,
								      dd4hep::Assembly &envelope){

  double containerThickness = ZMinus_PE_Thickness + ZMinus_Cu_Thickness;
  double z_position = -Ecal_barrel_halfZ -containerThickness/2. - env_safety*2.0;
  double container_x_dim = top_dim_x/2. + module_thickness*sin(M_PI/4.) - InnerServicesWidth/2.;
  double Cu_Thickness = ZMinus_Cu_Thickness;
  double PE_Thickness = ZMinus_PE_Thickness;

  for(int i=0; i<=1; i++){

    if(Ecal_cables_gap < containerThickness) return false;
    //Control::Abort(
    //	       "EcalBarrel_EndCap services: Cu+PE thicknesses exceed Ecal_cables_gap",
    //	       MOKKA_ERROR_BAD_DATABASE_PARAMETERS);

    dd4hep::Box ContainerSolid(container_x_dim/2., module_thickness/2., containerThickness/2.); 

    std::string containerLogical_name  = "containerLogical"+dd4hep::_toString(i,"_%d");
    dd4hep::Volume containerLogical(containerLogical_name,ContainerSolid,air);


    dd4hep::Box PESolid(container_x_dim/2., module_thickness/2., PE_Thickness/2.);
 
    std::string PELogical_name  = "PELogical"+dd4hep::_toString(i,"_%d");
    dd4hep::Volume PELogical(PELogical_name,PESolid,polyethylene);

    dd4hep::Position  PosPE(0,0,containerThickness/2. - PE_Thickness/2.);
      
    pVol = containerLogical.placeVolume(PELogical,PosPE);
      
    dd4hep::Box Cu_Solid(container_x_dim/2., module_thickness/2., Cu_Thickness/2.); 
      
    std::string Cu_Logical_name  = "Cu_Logical"+dd4hep::_toString(i,"_%d");
    dd4hep::Volume Cu_Logical(Cu_Logical_name,Cu_Solid,copper);
      
    dd4hep::Position  PosCu(0,0,-containerThickness/2. + Cu_Thickness/2.);
      
    pVol = containerLogical.placeVolume(Cu_Logical,PosCu);

      
    for (int stave_id = 1; stave_id < 9 ; stave_id++){
	
      double phirot = (stave_id-1) * M_PI/4;
	
      dd4hep::RotationZYX rot;
      dd4hep::Position stavePosition;
	
      if(z_position > 0) {
	dd4hep::RotationZYX rot_P(-phirot,0,0);
	rot = rot_P;
	dd4hep::Position stavePosition_P((Ecal_outer_radius - module_thickness/2.)*sin(phirot)+(InnerServicesWidth/2. + container_x_dim/2.)*cos(-phirot),
					 (Ecal_outer_radius - module_thickness/2.)*cos(phirot)+(InnerServicesWidth/2. + container_x_dim/2.)*sin(-phirot), 
					 z_position);

	stavePosition = stavePosition_P;
      }
      else
      {
	dd4hep::RotationZYX rot_M(phirot,M_PI,0);
	rot = rot_M;
	dd4hep::Position stavePosition_M((Ecal_outer_radius - module_thickness/2.)*sin(phirot)+(-InnerServicesWidth/2. - container_x_dim/2.)*cos(-phirot),
					 (Ecal_outer_radius - module_thickness/2.)*cos(phirot)+(-InnerServicesWidth/2. - container_x_dim/2.)*sin(-phirot),
					 z_position);
	    
	stavePosition = stavePosition_M;
      }
	
      dd4hep::Transform3D tran3D(rot,stavePosition);
      pVol = envelope.placeVolume(containerLogical,tran3D);
	
    }
      
    containerThickness = ZPlus_PE_Thickness + ZPlus_Cu_Thickness;
    z_position = Ecal_barrel_halfZ + containerThickness/2. + env_safety*2.0;
      
    Cu_Thickness = ZPlus_Cu_Thickness;
    PE_Thickness = ZPlus_PE_Thickness;
      
  }


  return true;
}




// to build Ecal Barrel service into the service assembly
bool BuildSitCables::DoBuildSitCables(dd4hep::PlacedVolume &pVol, dd4hep::Assembly &envelope){

  //First place the Cylinder

  if(ftd4to7_tpc_radial_gap < Sit_cables_cylinder_thickness)
    throw std::runtime_error("SServices00: the ftd-tpc radial gap is less than Sit cable thickness");

  double SitTube_inner_radius = TPC_inner_radius -
    ftd4to7_tpc_radial_gap/2.;

  double z_start_3 = TPC_Ecal_Hcal_barrel_halfZ * z3_position_ReltoTPCLength;

  double petalairthickness_half = 0.5 * ( petal_cp_support_thickness + 2.0*disk_si_thickness);

  double max_half_thickness_disk_3 = petal_support_zoffset + petalairthickness_half ;

  double z_half_len = (TPC_Ecal_Hcal_barrel_halfZ -
		       z_start_3) / 2.;

  dd4hep::Tube SitTubeSolid (SitTube_inner_radius,
			     SitTube_inner_radius + Sit_cables_cylinder_thickness,
			     z_half_len,
			     0., 2 * M_PI);

  dd4hep::Volume SitTubeLog0("SitTubeLog0",SitTubeSolid, aluminium);
  dd4hep::Volume SitTubeLog1("SitTubeLog1",SitTubeSolid, aluminium);

  dd4hep::rec::Vector3D ocyl(  SitTube_inner_radius + Sit_cables_cylinder_thickness/2.  , 0. , 0. ) ;

  dd4hep::rec::VolCylinder cylSurf0( SitTubeLog0 , dd4hep::rec::SurfaceType( dd4hep::rec::SurfaceType::Helper ) ,
				     0.5*Sit_cables_cylinder_thickness  ,
				     0.5*Sit_cables_cylinder_thickness , ocyl );

  dd4hep::rec::VolCylinder cylSurf1( SitTubeLog1 , dd4hep::rec::SurfaceType( dd4hep::rec::SurfaceType::Helper ) ,
				     0.5*Sit_cables_cylinder_thickness  ,
				     0.5*Sit_cables_cylinder_thickness , ocyl );

  dd4hep::rec::volSurfaceList( detElt )->push_back( cylSurf0 );
  dd4hep::rec::volSurfaceList( detElt )->push_back( cylSurf1 );

  dd4hep::Position Cylinder_pos1(0, 0, z_start_3 + z_half_len);
  pVol = envelope.placeVolume(SitTubeLog0,Cylinder_pos1);

  dd4hep::Position Cylinder_pos2(0, 0, -z_start_3 - z_half_len);
  pVol = envelope.placeVolume(SitTubeLog1,Cylinder_pos2);


  //Then place the cone

  double z_start_2 = TPC_Ecal_Hcal_barrel_halfZ * z2_position_ReltoTPCLength;

  petalairthickness_half = 0.5 * ( petal_cp_support_thickness + 2.0*disk_si_thickness);

  double max_half_thickness_disk_2 = petal_support_zoffset + petalairthickness_half ;

  double FTD2_outer_radius = SIT1_Radius + ftd2_sit1_radial_diff;

  double FTD3_outer_radius = SIT2_Radius + ftd3_sit2_radial_diff;

  double zPlane[2];
  zPlane[0] = z_start_2 + max_half_thickness_disk_2 + SurfaceTolerance;

  zPlane[1] = z_start_3 - 2.0*max_half_thickness_disk_3 - SurfaceTolerance;

  double rInner[2];
  rInner[0] = FTD2_outer_radius + SurfaceTolerance;
  rInner[1] = FTD3_outer_radius + SurfaceTolerance;

  double rOuter[2];
  rOuter[0] = FTD2_outer_radius + FTD2_cone_thickness;
  rOuter[1] = FTD3_outer_radius + FTD3_cone_thickness;

  std::vector<double> rmin;
  std::vector<double> rmax;
  std::vector<double> z;
  for(int i=0;i<=1;i++) {
    rmin.push_back(rInner[i]);
    rmax.push_back(rOuter[i]);
    z.push_back(zPlane[i]);
  }

  dd4hep::ConeSegment SitConeSolid( (z[1]-z[0])/2., rmin[0], rmax[0], rmin[1], rmax[1],0, 2.0 * M_PI);

  dd4hep::Volume SitConeLog0("SitConeLog0", SitConeSolid, aluminium);
  dd4hep::Volume SitConeLog1("SitConeLog1", SitConeSolid, aluminium);

  // Print out some parameters
  std::cout <<"\n   - SIT Cone Cable services: "
           <<"\n   - zPlane[0] = "
           <<zPlane[0]
           <<"\n   - zPlane[1] = "
           <<zPlane[1]
           <<"\n   - rInner[0] = "
           << rInner[0]
           <<"\n   - rInner[1] = "
           << rInner[1]
           <<"\n   - rOuter[0] = "
           <<rOuter[0]
           <<"\n   - rOuter[1] = "
           <<rOuter[1]
           <<std::endl;

  const double dr    = rmin[1] - rmin[0] ;
  const double theta = atan2( dr , z[1] - z[0] ) ;

  double coneThickness = ( ( rmax[0] - rmin[0] ) +  ( rmax[1] - rmin[1] ) )   / 2. ;

  dd4hep::rec::Vector3D ocon( rmin[0] + 0.5 * ( dr + coneThickness ), 0. , 0. );

  dd4hep::rec::Vector3D v( 1. , 0. , theta, dd4hep::rec::Vector3D::spherical ) ;

  dd4hep::rec::VolCone conSurf0( SitConeLog0 , dd4hep::rec::SurfaceType( dd4hep::rec::SurfaceType::Helper ) ,
                                 0.5*coneThickness  , 0.5*coneThickness , v, ocon );

  dd4hep::rec::VolCone conSurf1( SitConeLog1, dd4hep::rec::SurfaceType( dd4hep::rec::SurfaceType::Helper ) ,
                                 0.5*coneThickness  , 0.5*coneThickness , v, ocon );

  dd4hep::rec::volSurfaceList( detElt )->push_back( conSurf0 );
  dd4hep::rec::volSurfaceList( detElt )->push_back( conSurf1 );


  double coneZPos = (z[1]+z[0]) /2. ;
  pVol = envelope.placeVolume(SitConeLog0, Transform3D( RotationZYX() , Position(0., 0., coneZPos ) ) );

  dd4hep::RotationZYX rot(0,0, M_PI);
  pVol = envelope.placeVolume(SitConeLog1, Transform3D( rot , Position(0., 0., -coneZPos ) ));


  //Then place the disk
  double Sit_cables_disk_thickness =((1+TUBE_IPOuterBulge_end_radius/TPC_inner_radius)*0.5*
				     SServices_FTD7_cables_thickness );

  dd4hep::Tube SitDiskSolid(TUBE_IPOuterBulge_end_envradius + SurfaceTolerance,
			    TPC_inner_radius,
			    Sit_cables_disk_thickness/2. - SurfaceTolerance,
			    0., 2 * M_PI);

  dd4hep::Volume SitDiskLog("SitDiskLog", SitDiskSolid, aluminium);

  dd4hep::Position disk_pos1(0, 0, TUBE_IPOuterBulge_end_z + Sit_cables_disk_thickness/2.);
  pVol = envelope.placeVolume(SitDiskLog,disk_pos1);

  dd4hep::Position disk_pos2(0, 0, -TUBE_IPOuterBulge_end_z - Sit_cables_disk_thickness/2.);
  pVol = envelope.placeVolume(SitDiskLog,disk_pos2);

  // Print out some parameters
  std::cout <<"\n   - SIT12_cable_thickness = "
	    <<Sit_cables_cylinder_thickness
	    <<"\n   - SIT1_Radius = "
	    <<SIT1_Radius
	    <<"\n   - SIT2_Radius = "
	    <<SIT2_Radius
	    <<"\n   - SServices_FTD2_cone_thickness = "
	    <<FTD2_cone_thickness
	    <<"\n   - SServices_FTD3_cone_thickness = "
	    <<FTD3_cone_thickness
	    <<"\n   - TUBE_IPOuterBulge_end_z = "
	    <<TUBE_IPOuterBulge_end_z
	    <<"\n   - SServices_FTD7_cables_thickness = "
	    <<SServices_FTD7_cables_thickness
	    <<std::endl;


  return true;
}


bool BuildHcalBarrel_EndCapServices::FillHcalServicesModuleWithInnerServices(dd4hep::PlacedVolume &pVol,
									     dd4hep::Volume &ModuleLogicalZMinus,
									     dd4hep::Volume &ModuleLogicalZPlus){

  double StainlessSteel_Thickness = ZMinus_StainlessSteel_Thickness;
  double Cu_Thickness = ZMinus_Cu_Thickness;
  double PE_Thickness = ZMinus_PE_Thickness;

  dd4hep::Volume motherLogical = ModuleLogicalZMinus;

  double yShift = Hcal_y_dim2_for_x/2.;

  for(int i=0; i<=1; i++)
  {
    dd4hep::Position layerPositionSteel(0,
					Ecal_cables_gap/2. - StainlessSteel_Thickness/2.,
					yShift);

    if(!PlaceHcalInnerServicesLayer(pVol,motherLogical,
				    stainless_steel,
				    StainlessSteel_Thickness, layerPositionSteel))

      return false;

    dd4hep::Position layerPositionPoly(0,
				       Ecal_cables_gap/2. -StainlessSteel_Thickness - PE_Thickness/2.,
				       yShift);

    if(!PlaceHcalInnerServicesLayer(pVol,motherLogical,
				    polyethylene,
				    PE_Thickness, layerPositionPoly))

      return false;

    dd4hep::Position layerPositionCopper(0,
					 Ecal_cables_gap/2.-StainlessSteel_Thickness-PE_Thickness - Cu_Thickness/2.,
					 yShift);

    if(!PlaceHcalInnerServicesLayer(pVol,motherLogical,
				    copper,
				    Cu_Thickness, layerPositionCopper))

      return false;

    StainlessSteel_Thickness = ZPlus_StainlessSteel_Thickness;
    Cu_Thickness = ZPlus_Cu_Thickness;
    PE_Thickness = ZPlus_PE_Thickness;


    motherLogical = ModuleLogicalZPlus;
  }

  // Print out some parameters
  std::cout <<"\n   - HcalServicesModule_InnerServicesWidth = "
	    <<InnerServicesWidth
	    <<"\n   - HcalServicesModule_ZMinus_Cu_Thickness = "
	    <<ZMinus_Cu_Thickness
	    <<"\n   - HcalServicesModule_ZMinus_PE_Thickness = "
	    <<ZMinus_PE_Thickness
	    <<"\n   - HcalServicesModule_ZMinus_StainlessSteel_Thickness = "
	    <<ZMinus_StainlessSteel_Thickness
	    <<"\n   - HcalServicesModule_ZPlus_Cu_Thickness = "
	    <<ZPlus_Cu_Thickness
	    <<"\n   - HcalServicesModule_ZPlus_PE_Thickness = "
	    <<ZPlus_PE_Thickness
	    <<"\n   - HcalServicesModule_ZPlus_StainlessSteel_Thickness = "
	    <<ZPlus_StainlessSteel_Thickness
	    <<"\n   - HcalServices_outer_Cu_thickness = "
	    <<HcalServices_outer_Cu_thickness
	    <<"\n   - HcalServices_outer_FR4_thickness = "
	    <<HcalServices_outer_FR4_thickness
	    <<"\n   - Hcal_R_max = "
	    <<Hcal_R_max
	    <<"\n   - Hcal_back_plate_thickness = "
	    <<Hcal_back_plate_thickness
	    <<"\n   - Hcal_nlayers = "
	    <<Hcal_nlayers
	    <<"\n   - Hcal_radiator_thickness = "
	    <<Hcal_radiator_thickness
	    <<"\n   - Hcal_stave_gaps = "
	    <<Hcal_stave_gaps
	    <<"\n   - Hcal_steel_cassette_thickness = "
	    <<Hcal_steel_cassette_thickness
	    <<std::endl;

  return true;
}



bool BuildHcalBarrel_EndCapServices::PlaceHcalInnerServicesLayer(dd4hep::PlacedVolume &pVol,
								 dd4hep::Volume &motherLogical,
								 dd4hep::Material layerMaterial,
								 double layerThickness,
								 dd4hep::Position &layerPosition) {
  
  dd4hep::RotationZYX rot(0,0,M_PI*0.5);

  dd4hep::Box Center_Solid(InnerServicesWidth/2. - 2*(SurfaceTolerance),
			   Hcal_total_dim_y/2. - 2*(SurfaceTolerance),
			   layerThickness/2. - 2*(SurfaceTolerance));

  dd4hep::Volume Center_Logical("HcalBarrel_EndCap_Center_Logical",Center_Solid,
				layerMaterial);

  dd4hep::Transform3D tran3D(rot,layerPosition);
  pVol = motherLogical.placeVolume(Center_Logical,tran3D);

  dd4hep::Box Box_Side(InnerServicesWidth/2. - 2*(SurfaceTolerance),
		       Hcal_total_dim_y*cos(M_PI/16)/2. - 4*(SurfaceTolerance),
		       layerThickness/2. - 2*(SurfaceTolerance));

  dd4hep::Solid Solid_Side(Box_Side);

  dd4hep::Solid motherSolid = motherLogical.solid();

  double xShift = (Hcal_bottom_dim_x + Hcal_total_dim_y * tan(M_PI/8.)) / 2.;

  dd4hep::Position sideLayerPosition = layerPosition;
  sideLayerPosition.SetX(xShift);

  dd4hep::RotationZYX  rotRight(-M_PI/8.,0,M_PI*0.5);
  dd4hep::Transform3D tran3DRight(rotRight,sideLayerPosition);

  dd4hep::IntersectionSolid intersectionSolidRight(motherSolid,Solid_Side,tran3DRight);

  dd4hep::Volume Logical_Right("HcalBarrel_EndCap_Logical_Right",intersectionSolidRight,
			       layerMaterial);

  dd4hep::Position pos(0,0,0);
  pVol = motherLogical.placeVolume(Logical_Right,pos);

  sideLayerPosition = layerPosition;
  sideLayerPosition.SetX(-xShift);

  dd4hep::RotationZYX  rotLeft(M_PI/8.,0,M_PI*0.5);
  dd4hep::Transform3D tran3DLeft(rotLeft,sideLayerPosition);

  dd4hep::IntersectionSolid intersectionSolidLeft(motherSolid,Solid_Side,tran3DLeft);

  dd4hep::Volume Logical_Left("HcalBarrel_EndCap_Logical_Left",intersectionSolidLeft,
			      layerMaterial);

  pVol = motherLogical.placeVolume(Logical_Left,pos);


  return true;

}


bool BuildHcalBarrel_EndCapServices::FillHcalServicesModuleWithHcalElectronicsInterface(
  dd4hep::PlacedVolume &pVol,
  dd4hep::Volume &ModuleLogicalZMinus,
  dd4hep::Volume &ModuleLogicalZPlus) {

  std::cout<<"   - BuildHcalElectronicsInterface = true " <<std::endl;

  double Hcal_layer_thickenss =
    (Hcal_total_dim_y - Hcal_back_plate_thickness) / Hcal_nlayers;

  double Hcal_chamber_thickness =
    Hcal_layer_thickenss - Hcal_radiator_thickness;

  if(Hcal_chamber_thickness <= 0)
    throw std::runtime_error("Hcal Barrel-EndCap services: Hcal chamber thicknesses  <= 0");

#ifdef VERBOSE
  std::cout << "Hcal Barrel-EndCap services: Hcal_chamber_thickness = " <<
    Hcal_chamber_thickness << std::endl;
#endif

  double Hcal_y_dim1_for_x  = Hcal_total_dim_y - Hcal_y_dim2_for_x;

  double layer_x_dim = 0;
  double layer_y_offset = 0;

  for(int layer_id=1; layer_id<=Hcal_nlayers; layer_id++) {

    layer_y_offset = (layer_id - 1)*Hcal_layer_thickenss +
      Hcal_radiator_thickness;

    //---- bottom barrel----
    if(layer_id * Hcal_layer_thickenss  < Hcal_y_dim1_for_x )
      layer_x_dim = Hcal_bottom_dim_x + 2 *
	layer_y_offset * tan(M_PI/8.);
    else   //----- top barrel ---
      layer_x_dim = Hcal_midle_dim_x - 2*
	( layer_y_offset + Hcal_chamber_thickness -
	  Hcal_y_dim1_for_x ) / tan(M_PI/8.);

    if(!FillHcalElectronicsInterfaceLayer(pVol,
					  ModuleLogicalZMinus, ModuleLogicalZPlus,
					  layer_y_offset, layer_x_dim))

      return false;

  }

  return true;
}


bool BuildHcalBarrel_EndCapServices::FillHcalElectronicsInterfaceLayer(dd4hep::PlacedVolume &pVol,
								       dd4hep::Volume &ModuleLogicalZMinus,
								       dd4hep::Volume &ModuleLogicalZPlus,
								       double layer_y_offset, double layer_x_dim) {

  double Hcal_y_dim1_for_x  = Hcal_total_dim_y - Hcal_y_dim2_for_x;

  double y_position = -Hcal_y_dim1_for_x/2. + layer_y_offset +
    Hcal_steel_cassette_thickness/2.;

  if(!PlaceHcalElectronicsInterfaceComponent(pVol,
					     ModuleLogicalZMinus,
					     ModuleLogicalZPlus,
					     S235,
					     Hcal_steel_cassette_thickness,y_position,layer_x_dim)
    )
    return false;

  y_position += (Hcal_steel_cassette_thickness/2. +
		 HcalServices_outer_FR4_thickness/2.);

  if(!PlaceHcalElectronicsInterfaceComponent(pVol,
					     ModuleLogicalZMinus,
					     ModuleLogicalZPlus,
					     PCB,
					     HcalServices_outer_FR4_thickness,y_position,layer_x_dim)
    )
    return false;


  y_position += (HcalServices_outer_FR4_thickness/2. + HcalServices_outer_Cu_thickness/2.);

  if(!PlaceHcalElectronicsInterfaceComponent(pVol,
					     ModuleLogicalZMinus,
					     ModuleLogicalZPlus,
					     copper,
					     HcalServices_outer_Cu_thickness,y_position,layer_x_dim)
    )
    return false;

  return true;
}


dd4hep::Solid BuildHcalBarrel_EndCapServices::CutLayer(dd4hep::Solid &layerSolid, double y_position) {

  double tolerance = 2 * dd4hep::mm;

  dd4hep::Box Solid_Side(InnerServicesWidth/2. + tolerance,
			 Hcal_total_dim_y + 8*(SurfaceTolerance),
			 Ecal_cables_gap + 8*(SurfaceTolerance));

  double xShift = (Hcal_bottom_dim_x + Hcal_total_dim_y * tan(M_PI/8.)) / 2.;
  double yShift = Hcal_y_dim2_for_x/2. - y_position;

  dd4hep::Position rightSideLayerPosition(xShift,yShift,0);

  dd4hep::RotationZYX rightRotationSide(-M_PI/8.,0,0);
  dd4hep::Transform3D tran3Dright(rightRotationSide,rightSideLayerPosition);

  dd4hep::SubtractionSolid subtractionRight(layerSolid,Solid_Side,tran3Dright);

  dd4hep::Position leftSideLayerPosition(-xShift,yShift,0);
  dd4hep::RotationZYX leftRotationSide(+M_PI/8.,0,0);
  dd4hep::Transform3D tran3Dleft(leftRotationSide,leftSideLayerPosition);

  dd4hep::SubtractionSolid subtractionLeft(subtractionRight,Solid_Side,tran3Dleft);

  dd4hep::Position centerLayerPosition(0,yShift,0);

  dd4hep::SubtractionSolid theFinalCut(subtractionLeft,Solid_Side,centerLayerPosition);

  return theFinalCut;

}


bool BuildHcalBarrel_EndCapServices::DoBuildHcalBarrel_EndCapServices(dd4hep::PlacedVolume &pVol,
								      dd4hep::Assembly &envelope) {

  Hcal_total_dim_y = Hcal_R_max * cos(M_PI/16) - Hcal_inner_radius;
  double Hcal_module_radius = Hcal_inner_radius + Hcal_total_dim_y;
  Hcal_y_dim2_for_x  =
    (Hcal_module_radius - Hcal_module_radius*cos(M_PI/8));
  double Hcal_y_dim1_for_x  = Hcal_total_dim_y - Hcal_y_dim2_for_x;
  Hcal_bottom_dim_x  =
    2.*Hcal_inner_radius*tan(M_PI/8.)- Hcal_stave_gaps;
  Hcal_midle_dim_x   =
    Hcal_bottom_dim_x + 2* Hcal_y_dim1_for_x*tan(M_PI/8.);
  Hcal_top_dim_x     =
    Hcal_midle_dim_x - 2 * Hcal_y_dim2_for_x/tan(M_PI/8.);

#ifdef VERBOSE
  double Hcal_outer_radius = Hcal_inner_radius + Hcal_total_dim_y;

  G4cout << "BuildHcalBarrel_EndCapServices information: "
	 << "\n                       Hcal_outer_radius = "
	 << Hcal_outer_radius
	 << "\n                       module thickness = "
	 << Hcal_total_dim_y
	 << "\n                       Hcal_R_max = "
	 << Hcal_R_max
	 << "\n                       Hcal_bottom_dim_x = "
	 << Hcal_bottom_dim_x
	 << G4endl;
#endif

  double BHX  = Hcal_bottom_dim_x /2. - SurfaceTolerance;
  double MHX  = Hcal_midle_dim_x / 2. - SurfaceTolerance;
  double THX  = Hcal_top_dim_x / 2. - SurfaceTolerance;
  double YX1H = Hcal_y_dim1_for_x / 2.;
  double YX2H = Hcal_y_dim2_for_x / 2.;
  double DHZ  = Ecal_cables_gap/ 2. - SurfaceTolerance;

  // Attention: on bâtit le module dans la verticale
  // à cause du G4Trd et on le tourne avant de le positioner
  dd4hep::Trapezoid Bottom(BHX, MHX, DHZ, DHZ, YX1H);

  dd4hep::Trapezoid Top(MHX, THX, DHZ, DHZ, YX2H);

  dd4hep::Position pos(0, 0, YX1H + YX2H);
  dd4hep::UnionSolid ModuleSolid( Bottom, Top, pos);


  dd4hep::Volume ModuleLogicalZMinus("ServicesHcalModuleZMinus", ModuleSolid, air);

  dd4hep::Volume  ModuleLogicalZPlus("ServicesHcalModuleZPlus", ModuleSolid, air);

  //First place layer models of services coming from Ecal and TPC:
  if(!FillHcalServicesModuleWithInnerServices(pVol,ModuleLogicalZMinus,ModuleLogicalZPlus))
    return false;

  //Then place layer models of HCAL electronics interface:
  if(!(BuildHcalElectronicsInterface == 0))
    if(!FillHcalServicesModuleWithHcalElectronicsInterface(pVol,
							   ModuleLogicalZMinus,ModuleLogicalZPlus))
      return false;

  double Y = Hcal_inner_radius + YX1H;
  double stave_phi_offset = 0;

  for (int stave_id = 1;
       stave_id <= 8;
       stave_id++)
  {
    double module_z_offset =
      - TPC_Ecal_Hcal_barrel_halfZ - Ecal_cables_gap/2. - env_safety*2.0;

    double phirot = stave_phi_offset;

    dd4hep::RotationZYX rot(0,phirot,-M_PI*0.5);
    dd4hep::Position stave_pos1(Y*sin(phirot),
				Y*cos(phirot),
				module_z_offset);
    dd4hep::Transform3D tran3D1(rot,stave_pos1);
    pVol = envelope.placeVolume(ModuleLogicalZMinus,tran3D1);


    module_z_offset = - module_z_offset;

    dd4hep::Position stave_pos2(Y*sin(phirot),
				Y*cos(phirot),
				module_z_offset);
    dd4hep::Transform3D tran3D2(rot,stave_pos2);
    pVol = envelope.placeVolume(ModuleLogicalZPlus,tran3D2);

    stave_phi_offset -=  M_PI/4;
  }

  return true;
}



// to build VXD calbe cone into the service assembly
bool BuildVXDCables::DoBuildVXDCables(dd4hep::PlacedVolume &pVol, dd4hep::Assembly &envelope){
  
  double zPlane[2];
  zPlane[0] = VXD_cable_z_start;

  zPlane[1] = VXD_cable_z_end;

  double rInner[2];
  rInner[0] = VXD_cable_inner1_radius + SurfaceTolerance ;
  rInner[1] = VXD_cable_inner2_radius + SurfaceTolerance ;

  double rOuter[2];

  VXD_cable_cone_middle_thickness = std::sqrt( VXD_cable_cross_section_area/M_PI
					       + (VXD_cable_inner1_radius+VXD_cable_inner2_radius)/2.0
					       * (VXD_cable_inner1_radius+VXD_cable_inner2_radius)/2.0 )
    - (VXD_cable_inner1_radius+VXD_cable_inner2_radius)/2.0 ;

  VXD_cable_outer1_radius = VXD_cable_inner1_radius + VXD_cable_cone_middle_thickness ;
  VXD_cable_outer2_radius = VXD_cable_inner2_radius + VXD_cable_cone_middle_thickness ;

  rOuter[0] = VXD_cable_outer1_radius ;
  rOuter[1] = VXD_cable_outer2_radius ;

  std::vector<double> rmin;
  std::vector<double> rmax;
  std::vector<double> z;
  for(int i=0;i<=1;i++) {
    rmin.push_back(rInner[i]);
    rmax.push_back(rOuter[i]);
    z.push_back(zPlane[i]);
  }

#if 0  // fixme: intersection with 'inverse' cones currently not working -> use a cylinder for now ....

  dd4hep::ConeSegment VXDConeSolid( (z[1]-z[0])/2., rmin[0], rmax[0], rmin[1], rmax[1],0, 2.0 * M_PI);

  dd4hep::Volume VXDConeLog0("VXDConeLog1", VXDConeSolid, copper);
  dd4hep::Volume VXDConeLog1("VXDConeLog1", VXDConeSolid, copper);

  const double dr    = rmin[1] - rmin[0] ;
  const double theta = atan2( dr , z[1] - z[0] ) ;

  if( theta < 0 )   // reverse direction of v ...
    theta += M_PI ;

  double coneThickness = ( ( rmax[0] - rmin[0] ) +  ( rmax[1] - rmin[1] ) )   / 2. ;

  dd4hep::rec::Vector3D ocon( rmin[0] + 0.5 * ( dr + coneThickness ) , 0. , 0. );

  dd4hep::rec::Vector3D v( 1. , 0. , theta , dd4hep::rec::Vector3D::spherical ) ;

  dd4hep::rec::VolCone conSurf0( VXDConeLog0 , dd4hep::rec::SurfaceType( dd4hep::rec::SurfaceType::Helper ) ,
				 0.5*coneThickness  , 0.5*coneThickness , v, ocon );

  dd4hep::rec::VolCone conSurf1( VXDConeLog1, dd4hep::rec::SurfaceType( dd4hep::rec::SurfaceType::Helper ) ,
				 0.5*coneThickness  , 0.5*coneThickness , v, ocon );

  dd4hep::rec::volSurfaceList( detElt )->push_back( conSurf0 );
  dd4hep::rec::volSurfaceList( detElt )->push_back( conSurf1 );


  double coneZPos = (z[1]+z[0]) /2. ;
  pVol = envelope.placeVolume(VXDConeLog0, Transform3D( RotationZYX() , Position(0., 0., coneZPos ) ) );

  dd4hep::RotationZYX rot(0,0, M_PI);
  pVol = envelope.placeVolume(VXDConeLog1, Transform3D( rot , Position(0., 0., -coneZPos ) ));

#else

  double VXDTube_inner_radius = ( rmin[0]+rmin[1]) / 2. ;
  double VXD_cables_cylinder_thickness = VXD_cable_cone_middle_thickness ;
  double z_half_len = std::fabs( (z[1]-z[0])/2.) ;

  dd4hep::Tube VXDTubeSolid (VXDTube_inner_radius,
			     VXDTube_inner_radius + VXD_cables_cylinder_thickness,
			     z_half_len,
			     0., 2 * M_PI);

  dd4hep::Volume VXDTubeLog0("VXDTubeLog0",VXDTubeSolid, copper);
  dd4hep::Volume VXDTubeLog1("VXDTubeLog1",VXDTubeSolid, copper);

  dd4hep::rec::Vector3D ocyl(  VXDTube_inner_radius + VXD_cables_cylinder_thickness/2.  , 0. , 0. ) ;

  dd4hep::rec::VolCylinder cylSurf0( VXDTubeLog0 , dd4hep::rec::SurfaceType( dd4hep::rec::SurfaceType::Helper ) ,
				     0.5*VXD_cables_cylinder_thickness  ,
				     0.5*VXD_cables_cylinder_thickness , ocyl );

  dd4hep::rec::VolCylinder cylSurf1( VXDTubeLog1 , dd4hep::rec::SurfaceType( dd4hep::rec::SurfaceType::Helper ) ,
				     0.5*VXD_cables_cylinder_thickness  ,
				     0.5*VXD_cables_cylinder_thickness , ocyl );

  dd4hep::rec::volSurfaceList( detElt )->push_back( cylSurf0 );
  dd4hep::rec::volSurfaceList( detElt )->push_back( cylSurf1 );

  dd4hep::Position Cylinder_pos1(0, 0, z[0] + z_half_len);
  pVol = envelope.placeVolume(VXDTubeLog0,Cylinder_pos1);

  dd4hep::Position Cylinder_pos2(0, 0, -z[0] - z_half_len);
  pVol = envelope.placeVolume(VXDTubeLog1,Cylinder_pos2);


#endif




  // Print out some parameters
  std::cout <<"\n   - VXD_cable_cross_section_area = "
	    <<VXD_cable_cross_section_area/dd4hep::mm/dd4hep::mm  <<" mm^2"
	    <<"\n   - VXD_cable_inner1_radius = "
	    <<VXD_cable_inner1_radius/dd4hep::mm  <<" mm"
	    <<"\n   - VXD_cable_inner2_radius = "
	    <<VXD_cable_inner2_radius/dd4hep::mm  <<" mm"
	    <<"\n   - VXD_cable_cone_middle_thickness = "
	    <<VXD_cable_cone_middle_thickness/dd4hep::mm  <<" mm"
	    <<"\n   - VXD_cable_outer1_radius = "
	    <<VXD_cable_outer1_radius/dd4hep::mm  <<" mm"
	    <<"\n   - VXD_cable_outer2_radius = "
	    <<VXD_cable_outer2_radius/dd4hep::mm  <<" mm"
	    <<"\n   - VXD_cable_z_start = "
	    <<VXD_cable_z_start/dd4hep::mm  <<" mm"
	    <<"\n   - VXD_cable_z_end = "
	    <<VXD_cable_z_end/dd4hep::mm  <<" mm"
	    <<std::endl;


  return true;
}



DECLARE_DETELEMENT(SServices00,create_element)

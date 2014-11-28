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
#include "XMLHandlerDB.h"

#include "SServices00.h"
 
using namespace std;
using namespace DD4hep;
using namespace dd4hep;
using namespace DD4hep::Geometry;


static Ref_t create_element(LCDD& lcdd, xml_h element, Ref_t)  {

  static double tolerance = 0e0;

  xml_det_t x_det = element;
  string det_name = x_det.nameStr();
  
  int det_id = x_det.id();
  DetElement sdet (det_name,det_id);
  Volume motherVol = lcdd.pickMotherVolume(sdet);
  
  Assembly envelope_assembly( det_name + "assembly" ) ;
  PlacedVolume pv;


//====================================================================
// build all services
//====================================================================
  cout << "\nBuilding SServices00"<< endl;

  //==================================================
  //           BuildTPCEndplateServices
  //==================================================
  
  // start to prepare the Material and geometry as Mokka
  double TPC_Ecal_Hcal_barrel_halfZ =lcdd.constant<double>("TPC_Ecal_Hcal_barrel_halfZ");

  XMLHandlerDB db = XMLHandlerDB(  x_det.child( _Unicode( TPC_Cooling ) ) ) ;

  const int MAX_TPC_RINGS =  db->fetchInt("number_of_rings");;
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
  //TPCEndplateServices.setMaterial(lcdd.material("copper"));
  TPCEndplateServices.setMaterial(lcdd.material("Cu"));
  TPCEndplateServices.sethalfZ(TPC_Ecal_Hcal_barrel_halfZ);
  for(int i=0;i<MAX_TPC_RINGS;i++)
    TPCEndplateServices.settpcEndplateServicesRing_R_ro(tpcEndplateServices_R[i],tpcEndplateServices_r[i]);
  
  TPCEndplateServices.DoBuildTPCEndplateServices(&pv,&envelope_assembly);
 


 
  //==================================================
  //           BuildEcalBarrelServices
  //==================================================

  // start to prepare the Material and geometry as Mokka
  double Ecal_outer_radius = lcdd.constant<double>("Ecal_outer_radius");
  double Ecal_inner_radius = lcdd.constant<double>("Ecal_inner_radius");
  double module_thickness = Ecal_outer_radius - Ecal_inner_radius;
 
  // module barrel key parameters
  double bottom_dim_x = 2. * tan(M_PI/8.) * Ecal_inner_radius +
    module_thickness/sin(M_PI/4.);
  
  double top_dim_x = bottom_dim_x - 2 * module_thickness;
  
  double RailHeight = lcdd.constant<double>("EcalBarrelServices_RailHeight");

  cout<<"\n Ecal_inner_radius: " << Ecal_inner_radius
      <<"\n Ecal_outer_radius: " << Ecal_outer_radius
      <<"\n  module_thickness: " << module_thickness 
      <<"\n      bottom_dim_x: " << bottom_dim_x
      <<"\n         top_dim_x: " << top_dim_x
      <<"\n        RailHeight: " << RailHeight
      <<"\n" <<endl;


  double RailDistanceToRight = lcdd.constant<double>("EcalBarrelServices_RailDistanceToRight");
  double RailSeparation = lcdd.constant<double>("EcalBarrelServices_RailSeparation");
  double RailWidth =  lcdd.constant<double>("EcalBarrelServices_RailWidth");

  double ZMinus_FirstInterrail_PE_Thickness = lcdd.constant<double>("EcalBarrelServices_ZMinus_FirstInterrail_PE_Thickness");
  double ZMinus_FirstInterrail_Cu_Thickness = lcdd.constant<double>("EcalBarrelServices_ZMinus_FirstInterrail_Cu_Thickness");
  double ZMinus_SecondInterrail_Cu_Thickness = lcdd.constant<double>("EcalBarrelServices_ZMinus_SecondInterrail_Cu_Thickness");

  double ZPlus_FirstInterrail_PE_Thickness = lcdd.constant<double>("EcalBarrelServices_ZPlus_FirstInterrail_PE_Thickness");
  double ZPlus_FirstInterrail_Cu_Thickness = lcdd.constant<double>("EcalBarrelServices_ZPlus_FirstInterrail_Cu_Thickness");
  double ZPlus_SecondInterrail_Cu_Thickness = lcdd.constant<double>("EcalBarrelServices_ZPlus_SecondInterrail_Cu_Thickness");

  double OldRailSeparation = lcdd.constant<double>("EcalBarrelServices_RailSeparation");

  // start to build Ecal Barrel services with the class BuildEcalBarrelServices
  BuildEcalBarrelServices EcalBarrelServices;

  EcalBarrelServices.setMaterialAir(lcdd.air());
  //EcalBarrelServices.setMaterialAluminium(lcdd.material("aluminium"));
  EcalBarrelServices.setMaterialAluminium(lcdd.material("Al"));
  //EcalBarrelServices.setMaterialPolyethylene(lcdd.material("polyethylene"));
  EcalBarrelServices.setMaterialPolyethylene(lcdd.material("G4_POLYSTYRENE"));
  //EcalBarrelServices.setMaterialCopper(lcdd.material("copper"));
  EcalBarrelServices.setMaterialCopper(lcdd.material("Cu"));

  EcalBarrelServices.sethalfZ(TPC_Ecal_Hcal_barrel_halfZ);
  EcalBarrelServices.setTopDimX(top_dim_x);
  EcalBarrelServices.setRailHeight(RailHeight);
  EcalBarrelServices.setOutRadius(Ecal_outer_radius);
  EcalBarrelServices.setModuleThickness(module_thickness);

  EcalBarrelServices.setRailDistanceToRight(RailDistanceToRight);
  EcalBarrelServices.setRailSeparation(RailSeparation);
  EcalBarrelServices.setRailWidth(RailWidth);

  EcalBarrelServices.setZMinus_FirstInterrail_PE_Thickness(ZMinus_FirstInterrail_PE_Thickness);
  EcalBarrelServices.setZMinus_FirstInterrail_Cu_Thicknesst(ZMinus_FirstInterrail_Cu_Thickness);
  EcalBarrelServices.setZMinus_SecondInterrail_Cu_Thickness(ZMinus_SecondInterrail_Cu_Thickness);

  EcalBarrelServices.setZPlus_FirstInterrail_PE_Thickness(ZPlus_FirstInterrail_PE_Thickness);
  EcalBarrelServices.setZPlus_FirstInterrail_Cu_Thickness(ZPlus_FirstInterrail_Cu_Thickness);
  EcalBarrelServices.setZPlus_SecondInterrail_Cu_Thickness(ZPlus_SecondInterrail_Cu_Thickness);
  EcalBarrelServices.setOldRailSeparation(OldRailSeparation);

  EcalBarrelServices.DoBuildEcalBarrelServices(&pv,&envelope_assembly);




  //==================================================
  //          BuildEcalBarrel_EndCapServices
  //==================================================

  // start to prepare the Material and geometry as Mokka
  double ZMinus_PE_Thickness = lcdd.constant<double>("EcalBarrel_EndCapServices_ZMinus_PE_Thickness");
  double ZMinus_Cu_Thickness = lcdd.constant<double>("EcalBarrel_EndCapServices_ZMinus_Cu_Thickness");
  double ZPlus_PE_Thickness = lcdd.constant<double>("EcalBarrel_EndCapServices_ZPlus_PE_Thickness");
  double ZPlus_Cu_Thickness = lcdd.constant<double>("EcalBarrel_EndCapServices_ZPlus_Cu_Thickness");
  double Ecal_cables_gap = lcdd.constant<double>("Ecal_cables_gap");
  double InnerServicesWidth = lcdd.constant<double>("HcalServicesModule_InnerServicesWidth");

  BuildEcalBarrel_EndCapServices EcalBarrel_EndCapServices;
  EcalBarrel_EndCapServices.setMaterialAir(lcdd.air());
  EcalBarrel_EndCapServices.setMaterialPolyethylene(lcdd.material("G4_POLYSTYRENE"));
  EcalBarrel_EndCapServices.setMaterialCopper(lcdd.material("Cu"));
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

  EcalBarrel_EndCapServices.DoBuildEcalBarrel_EndCapServices(&pv,&envelope_assembly);






  //BuildHcalBarrel_EndCapServices(&pv,&envelope_assembly);

  //BuildSitCables(&pv,&envelope_assembly);


//====================================================================
// Place services into the world volume
//====================================================================

  pv = motherVol.placeVolume(envelope_assembly);
  pv.addPhysVolID("system", det_id);
  sdet.setVisAttributes( lcdd, x_det.visStr(),  envelope_assembly);
  sdet.setPlacement(pv);

  return sdet;
}

DECLARE_DETELEMENT(SServices00,create_element)

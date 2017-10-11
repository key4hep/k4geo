//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "LcgeoExceptions.h"
#include "lcgeo.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"
#include "XMLHandlerDB.h"

#include <math.h>

using namespace std;
using namespace dd4hep;
using namespace lcgeo;

using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::volSurfaceList;
using dd4hep::rec::VolPlane;
using dd4hep::rec::FixedPadSizeTPCData;

/** Construction of TPC detector, ported from Mokka driver TPC10.cc
 * Mokka History:
 * - modified version of TPC driver by Ties Behnke
 * - modified version of TPC02 as TPC03 with selectable chamber gas -- Adrian Vogel, 2005-06-09
 * - modified version of TPC03 as TPC04 with limit of step length   -- Adrian Vogel, 2006-02-01
 * - introduced self-scalability, no superdriver needed anymore     -- Adrian Vogel, 2006-03-11
 * - modified version of TPC04 as TPC05 in order to have full MC
 *   information both at entry and exit hits in the TPC ,
 *   more realistic central electrode and endplate             -- Predrag Krstonosic, 2006-07-12
 * - implemented new GEAR interface -- K. Harder, T. Pinto Jayawardena                2007-07-31
 * - TPC10 implemented readout within the Gas volume and layered inner and outer wall -- SJA -- 2010-11-19
 *
 *  @author: F.Gaede, DESY, Nov 2013
 *
 */
static Ref_t create_element(Detector& theDetector, xml_h e, SensitiveDetector sens)  {

  //------------------------------------------
  //  See comments starting with '//**' for
  //     hints on porting issues
  //------------------------------------------

  
  xml_det_t    x_det = e;
  string       name  = x_det.nameStr();

  DetElement   tpc(  name, x_det.id()  ) ;
  

 // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , tpc ) ;

  dd4hep::xml::setDetectorTypeFlag( e, tpc ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return tpc ;

  //-----------------------------------------------------------------------------------

  PlacedVolume pv;  

  sens.setType("tracker");

  std::cout << " ** building TPC10_geo in lcgeo " << lcgeo::versionString() << std::endl ;

  //   //######################################################################################################################################################################
  //   //  code ported from TPC10::construct() :
  //   //##################################
  
  // const double phi1 =   0.0 * deg ;
  // const double phi2 = 360.0 * deg;
  const double phi1 =   0.0 ;
  const double phi2 =  2*M_PI ;
  //  
  const double dzTotal           = theDetector.constant<double>("TPC_Ecal_Hcal_barrel_halfZ") * 2. ; 
  const double rInner            = theDetector.constant<double>("TPC_inner_radius") ;
  const double rOuter            = theDetector.constant<double>("TPC_outer_radius") ;
    
    
  // Geometry parameters from the geometry environment and from the database
  // Database *db = new Database(env.GetDBName());
    
  XMLHandlerDB db(  x_det.child( _Unicode( global ) ) );
    
  const double TPCMaxStepLength    = db->fetchDouble("TPC_max_step_length") ;
  const double padHeight           = db->fetchDouble("TPC_pad_height") ;
  const double padWidth            = db->fetchDouble("TPC_pad_width") ;
  const double dr_InnerWall        = db->fetchDouble("dr_InnerWall") ;
  const double dr_InnerServiceArea = db->fetchDouble("dr_InnerServiceArea") ;
  const double dr_OuterServiceArea = db->fetchDouble("dr_OuterServiceArea") ;
  const double dr_OuterWall        = db->fetchDouble("dr_OuterWall") ;
  const double dz_Readout          = db->fetchDouble("dz_Readout") ;
  const double dz_Endplate         = db->fetchDouble("dz_Endplate") ;

  //    Material* const material_TPC_Gas = CGAGeometryManager::GetMaterial(db->fetchString("chamber_Gas"));
  Material material_TPC_Gas =  theDetector.material(db->fetchString("chamber_Gas") ) ;

  // #ifdef MOKKA_GEAR
  //   _gear_gas_material = material_TPC_Gas; 
  // #endif
    
  //unused:  const double sensitive_threshold_eV  = db->fetchDouble("sensitive_threshold_eV") ;
    

  //   db->exec("SELECT * FROM `cathode`;");
  //   db->getTuple();
  db = XMLHandlerDB(  x_det.child( _Unicode( cathode ) ) );
    
  const double dz_Cathode_Insulator   = db->fetchDouble("dz_Cathode_Insulator") ;
  const double dz_Cathode_Conductor   = db->fetchDouble("dz_Cathode_Conductor") ;
  Material material_Cathode_Insulator = theDetector.material(db->fetchString("material_Cathode_Insulator"));
  Material material_Cathode_Conductor = theDetector.material(db->fetchString("material_Cathode_Conductor"));
    
  const double dr_Cathode_Grip = db->fetchDouble("dr_Cathode_Grip") ;
  const double dz_Cathode_Grip = db->fetchDouble("dz_Cathode_Grip") ;
  Material material_Cathode_Grip = theDetector.material(db->fetchString("material_Cathode_Grip"));
    
  cout << " Cathode Grip Ring Material: " << material_Cathode_Grip->GetName() << " : Rad length = " << material_Cathode_Grip->GetMaterial()->GetRadLen() / mm << " mm." << std::endl;
    
  const double dz_Cathode = 2*(dz_Cathode_Insulator+dz_Cathode_Conductor);
    
  double tracking_tpc_ecal_gap = theDetector.constant<double>("Ecal_Tpc_gap") ;;
  double tracking_region_rmax = rOuter + tracking_tpc_ecal_gap - 0.1*mm;  // give 100 micron clearance 
    
  std::stringstream tracking_region_rmax_as_string;
  tracking_region_rmax_as_string <<  tracking_region_rmax;
    
  //fg needed ???
  // (*Control::globalModelParameters)["tracker_region_rmax"] = tracking_region_rmax_as_string.str();
  // (*Control::globalModelParameters)["tracker_region_zmax"] = env.GetParameterAsString("TPC_Ecal_Hcal_barrel_halfZ");
    
    
  // Calculate Dimentions needed for later. Note gas volume and endplate will be mirror placed ...
  const double rMin_GasVolume      = rInner + dr_InnerWall;
  const double rMax_GasVolume      = rOuter - dr_OuterWall;
    
  // note the will be two gas volumes one in each z-half. The cathode and readout are considered to be placed inside the Gas volume
  const double dz_GasVolume        = ( dzTotal/2.0 ) - dz_Endplate; 
    
  const double rMin_Sensitive      = rMin_GasVolume + dr_InnerServiceArea;
  const double rMax_Sensitive      = rMax_GasVolume - dr_OuterServiceArea;
  const double dz_Sensitive        = dz_GasVolume - ( dz_Cathode/2.0 + dz_Readout ); // the d_Cathode spans both halfs of the TPC 
  const double dz_Wall             = dzTotal - 2.0 * dz_Endplate ; // note field cage spans the complete length of the TPC Gas volume
    
    
  const int numberPadRows = (int)((rMax_Sensitive-rMin_Sensitive)/padHeight) ;


  // Materials to be used
  Material materialAir     = theDetector.material("G4_AIR");

  // Material mixture for end of endplate zone 
  //  db->exec("SELECT * FROM `endplate_mixture`;");

  //unused:  double endplate_mixture_total = 0.0;
  //unused:  double densityTotal = 0.0 ;

  //fg: don't create the material on the fly but define it in the xml file
  // // std::map<Material* const, double> material_fractions;
  // // //  while (db->getTuple()) 
  // // for(xml_coll_t c( x_det ,_U(endplate_mixture)); c; ++c)  {
  // //   xml_comp_t  x_row( c );
  // //   db = XMLHandlerDB( x_row )  ;
  // //   Material material = theDetector.material(db->fetchString("material"));
  // //   material_fractions[material]    = db->fetchDouble("percentage") * perCent; // fraction of material mix
  // //   endplate_mixture_total         += material_fractions[material];
  // //   densityTotal                   += material->GetDensity() * material_fractions[material]; 
  // // }
  // // if (fabs( endplate_mixture_total - 1) > 1E-06 || (fabs( 1 - endplate_mixture_total) > 1E-06 )) 
  // //   {
  // //     cout << "endplate_mixture_total = " << endplate_mixture_total << endl;
  // //     Control::Abort("TPC endplate material fractions do not add up to 100%",MOKKA_ERROR_BAD_DATABASE_PARAMETERS);
  // //   }
  // // Material *endplate_MaterialMix = new Material("TPC_endplate_mix", densityTotal, material_fractions.size());
  // // for(  std::map<Material* const, double>::iterator it=material_fractions.begin(); it!=material_fractions.end();++it )
  // //   {
  // //     endplate_MaterialMix->AddMaterial(it->first, it->second);
  // //   }
  
Material endplate_MaterialMix = theDetector.material( "TPC_endplate_mix" ) ;
  
  cout << "Endplate material mix Density   = " << endplate_MaterialMix->GetMaterial()->GetDensity() / g * cm3 << " g/cm3" << endl;
  cout << "Endplate material mix Radlength = " << endplate_MaterialMix->GetMaterial()->GetRadLen() / mm << " mm" << endl;
  
  
  //   // Visualisation attributes
  
  //   VisAttributes *wallVisAttributes = new VisAttributes(Colour(0.0, 0.5, 0.5)); // dull cyan
  //   wallVisAttributes->SetForceWireframe(false);
  //   wallVisAttributes->SetDaughtersInvisible(true);
  
  //   VisAttributes *cathodeVisAttributes = new VisAttributes(Colour(0.9, 0.3, 0.1)); // coppery brown
  //   cathodeVisAttributes->SetForceWireframe(false);
  
  
  // Some verbose output
  cout << "TPC10: Inner radius of the gas volume is " << std::setw(4) << rMin_GasVolume / mm << " mm." << endl;
  cout << "TPC10: Outer radius of the gas volume is " << std::setw(4) << rMax_GasVolume / mm << " mm." << endl;
  cout << "TPC10: Inner wall thickness is " << std::setw(4) << dr_InnerWall / mm << " mm." << endl;
  cout << "TPC10: Outer wall thickness is " << std::setw(4) << dr_OuterWall / mm << " mm." << endl;
  cout << "TPC10: Outer wall thickness is " << std::setw(4) << dr_OuterWall / mm << " mm." << endl;
  
  cout << "TPC10: Inner radius of the sensitive volume is " << std::setw(4) << rMin_Sensitive / mm << " mm." << endl;
  cout << "TPC10: Outer radius of the sensitive volume is " << std::setw(4) << rMax_Sensitive / mm << " mm." << endl;
  cout << "TPC10: Number of Pad Rows in the TPC  " << std::setw(4) << numberPadRows << endl;
  cout << "TPC10: Limiting the step length in the TPC to  " << std::setw(4) << TPCMaxStepLength / mm << " mm." << endl;
  
  //-------------------------------------------------------------------------------------------------------//
  
  //-------------------------------- TPC mother volume ----------------------------------------------------//
  //------------ Volume for the whole TPC, Field Cage, Cathode, and Endplate and Sensitive ----------------//
  
  Tube tpc_motherSolid(rInner ,rOuter ,dzTotal/2.0 , phi1 , phi1+phi2 ); 
  Volume tpc_motherLog(  "TPCLog", tpc_motherSolid, material_TPC_Gas );
  pv = envelope.placeVolume( tpc_motherLog ) ;
  tpc.setVisAttributes(theDetector,  "TPCMotherVis" ,  tpc_motherLog ) ;

  // VisAttributes* motherVisAttributes = new VisAttributes(Colour(0.0, 0.5, 0.5)); // dull cyan
  // motherVisAttributes->SetVisibility(false);
  // motherVisAttributes->SetDaughtersInvisible(true);
  // motherLog->SetVisAttributes(motherVisAttributes);
  
  
  cout << "TPC10: Total Gas material corresponds to " << ( ( (rOuter-dr_OuterWall) - (rInner + dr_InnerWall) ) / (material_TPC_Gas->GetMaterial()->GetRadLen() / mm ) * 100.0 ) 
       << "% of a radiation length." << endl;

  //-------------------------------------------------------------------------------------------------------//

  //-------------------------------- inner wall construction ----------------------------------------//

  Tube innerWallSolid(rInner ,rInner + dr_InnerWall ,dz_Wall / 2.0 , phi1 ,phi1+phi2 ); 
  Volume innerWallLog( "TPCInnerWallLog", innerWallSolid, materialAir ) ; 
  pv = tpc_motherLog.placeVolume( innerWallLog ) ;
  tpc.setVisAttributes(theDetector,  "CyanVis" ,  innerWallLog ) ;


  Vector3D ocyl(  rInner + 0.5*dr_InnerWall , 0. , 0. ) ;
  VolCylinder surfI( innerWallLog , SurfaceType( SurfaceType::Helper ) ,0.5*dr_InnerWall  , 0.5*dr_InnerWall , ocyl ) ;
  volSurfaceList( tpc )->push_back(  surfI ) ;

  int layerCounter = 0;
  double fracRadLengthInnerWall = 0;
  double rCursor = rInner ;
  //  double gear_inner_wall_material_total_density = 0.0;

//   db->exec("SELECT * FROM `innerWall`;");
//   while (db->getTuple()) {

  xml_comp_t  x_iWall=  x_det.child( _Unicode( innerWall) );
  for(xml_coll_t c( x_iWall , _Unicode( row) ); c; ++c)  {

    xml_comp_t  x_row( c );
    db = XMLHandlerDB( x_row )  ;

    const double dr = db->fetchDouble("dr") ;
    Material layerMaterial = theDetector.material(db->fetchString("material"));    
    Tube  layerSolid( rCursor, rCursor + dr , dz_Wall / 2.0, phi1, phi1+phi2);
    Volume layerLog( _toString( layerCounter ,"TPCInnerWallLayerLog_%02d") , layerSolid, layerMaterial );
    

    //    layerLog->SetVisAttributes(VisAttributes::Invisible);
    pv = innerWallLog.placeVolume( layerLog ) ;
    ++layerCounter;
    rCursor += dr ;
    fracRadLengthInnerWall += dr / layerMaterial->GetMaterial()->GetRadLen();
    
    cout << "TPC10: Add Material to Inner Wall: dr =  " << std::setw(4) << dr / mm << " mm. Material = " 
	 << layerMaterial->GetName() << " X0 = " << layerMaterial->GetMaterial()->GetRadLen() << "  " <<  dr / layerMaterial->GetMaterial()->GetRadLen() << "% X0" << endl;
    
// #ifdef MOKKA_GEAR
//     gear_inner_wall_material_total_density += layerMaterial->GetDensity() * (dr / dr_InnerWall);    
//     std::string material_name = db->fetchString("material");
//     cout << "TPC10: gear_inner_wall_material_total_density = " << std::setw(4) << gear_inner_wall_material_total_density << endl;   
//     if( gear_inner_wall_material_thicknesses.find( material_name ) ==  gear_inner_wall_material_thicknesses.end() ) {
//       gear_inner_wall_material_thicknesses[ material_name ] = dr;
//     }
//     else {
//       gear_inner_wall_material_thicknesses[ material_name ] += dr;
//     }
// #endif
    
  }
  
  cout << "TPC10: Inner wall material corresponds to " << int( fracRadLengthInnerWall * 1000) / 10. << "% of a radiation length." << endl;
  cout << "TPC10: Inner wall effective X0 = " << std::setw(4) << dr_InnerWall / fracRadLengthInnerWall<< endl;  


  //-------------------------------------------------------------------------------------------------------//

  //-------------------------------- outer wall construction ---------------------------------------------//

  Tube outerWallSolid( rOuter - dr_OuterWall ,rOuter ,dz_Wall / 2.0 ,phi1 , phi1+phi2) ;
  Volume outerWallLog( "TPCOuterWallLog" , outerWallSolid, materialAir);
  pv = tpc_motherLog.placeVolume( outerWallLog ) ;
  tpc.setVisAttributes(theDetector,  "CyanVis" ,  outerWallLog );

  ocyl.fill(   rOuter - 0.5*dr_OuterWall  , 0., 0. ) ;
  VolCylinder surfO( outerWallLog , SurfaceType( SurfaceType::Helper ) ,0.5*dr_OuterWall  , 0.5*dr_OuterWall , ocyl ) ;
  volSurfaceList( tpc )->push_back(  surfO ) ;

  layerCounter = 0;
  double fracRadLengthOuterWall = 0;
  rCursor = rOuter - dr_OuterWall ;
  //unused:  double gear_outer_wall_material_total_density = 0.0;
 
  // db->exec("SELECT * FROM `outerWall`;");
  // while (db->getTuple()) {
  xml_comp_t  x_oWall=  x_det.child( _Unicode( outerWall) );
  for(xml_coll_t c( x_oWall , _Unicode( row) ); c; ++c)  {
    xml_comp_t  x_row( c );
    db = XMLHandlerDB( x_row )  ;
    
    const double dr = db->fetchDouble("dr") ;
    Material layerMaterial = theDetector.material(db->fetchString("material"));    
    Tube  layerSolid( rCursor, rCursor + dr , dz_Wall / 2.0, phi1, phi1+phi2 );
    Volume layerLog(  _toString( layerCounter ,"TPCOuterWallLayerLog_%02d") , layerSolid, layerMaterial );
    //    layerLog->SetVisAttributes(VisAttributes::Invisible);
    pv = outerWallLog.placeVolume( layerLog ) ;
    ++layerCounter;
    rCursor += dr ;
    fracRadLengthOuterWall += dr / layerMaterial->GetMaterial()->GetRadLen();
    
    cout << "TPC10: Add Material to Outer Wall: dr =  " << std::setw(4) << dr / mm << " mm. Material = " << layerMaterial->GetName() << " X0 = "
	 << layerMaterial->GetMaterial()->GetRadLen() << "  " <<  dr / layerMaterial->GetMaterial()->GetRadLen() << "% X0" << endl;
    
// #ifdef MOKKA_GEAR
//     gear_outer_wall_material_total_density += layerMaterial->GetDensity() * (dr / dr_OuterWall);    
//     std::string material_name = db->fetchString("material");
//     cout << "TPC10: gear_outer_wall_material_total_density = " << std::setw(4) << gear_outer_wall_material_total_density << endl;   
//     if( gear_outer_wall_material_thicknesses.find( material_name ) ==  gear_outer_wall_material_thicknesses.end() )
//       {
// 	gear_outer_wall_material_thicknesses[ material_name ] = dr;
//       }
//     else
//       {
//       	gear_outer_wall_material_thicknesses[ material_name ] += dr;
//       }
//#endif  
    
  }
  cout << "TPC10: Outer wall material corresponds to " << int( fracRadLengthOuterWall * 1000) / 10.0 << "% of a radiation length." << endl;
  cout << "TPC10: Outer wall effective X0 = " << std::setw(4) << dr_OuterWall / fracRadLengthOuterWall   << endl;  
  //-----------------------------------------------------------------------------------------------//  


  //-------------------------------- cathode grip ring construction ----------------------------------------//
  // inner grip ring
  Tube cathodeInnerGripSolid( rMin_GasVolume , rMin_GasVolume + dr_Cathode_Grip, dz_Cathode_Grip / 2.0 , phi1, phi2);
  Volume cathodeInnerGripLog( "TPCCathodeInnerGripLog", cathodeInnerGripSolid, material_Cathode_Grip );
  pv = tpc_motherLog.placeVolume( cathodeInnerGripLog ) ;
  tpc.setVisAttributes(theDetector,  "GrayVis" ,  cathodeInnerGripLog );
  
  // outer grip ring
  Tube cathodeOuterGripSolid( rMax_GasVolume - dr_Cathode_Grip, rMax_GasVolume, dz_Cathode_Grip / 2.0 , phi1, phi2);
  Volume cathodeOuterGripLog("TPCCathodeOuterGripLog", cathodeOuterGripSolid, material_Cathode_Grip );
  pv = tpc_motherLog.placeVolume(cathodeOuterGripLog ) ;
  tpc.setVisAttributes(theDetector,  "GrayVis" ,  cathodeOuterGripLog );

  //-----------------------------------------------------------------------------------------------//  

  //-------------------------------- cathode construction ----------------------------------------//
  Tube cathodeSolid(  rMin_Sensitive, rMax_Sensitive, dz_Cathode / 2.0, phi1, phi2);
  Volume cathodeLog( "TPCCathodeLog", cathodeSolid, materialAir ) ;
  pv = tpc_motherLog.placeVolume( cathodeLog );
  tpc.setVisAttributes(theDetector,  "GrayVis" ,  cathodeLog );

  // insulator 
  Tube cathodeInsulatorSolid( rMin_Sensitive, rMax_Sensitive, (dz_Cathode_Insulator / 2.0)-0.00000001*mm, phi1, phi2);
  Volume cathodeInsulatorLog( "TPCcathodeInsulatorLog",cathodeInsulatorSolid, material_Cathode_Insulator);
  
  // place plus and minus z 
  pv = cathodeLog.placeVolume( cathodeInsulatorLog , Position( 0.0, 0.0, + dz_Cathode_Insulator / 2.0) );
  pv = cathodeLog.placeVolume( cathodeInsulatorLog , Position( 0.0, 0.0, - dz_Cathode_Insulator / 2.0) );

  tpc.setVisAttributes(theDetector,  "GrayVis" ,  outerWallLog );

  cout << "Cathode dz = " <<  dz_Cathode_Insulator << endl; 
  cout << "Place cathode +z at " <<  dz_Cathode_Insulator / 2.0 << endl; 
  cout << "Place cathode -z at " << -dz_Cathode_Insulator / 2.0 << endl;
  
  // conductor 
  Tube cathodeConductorSolid( rMin_Sensitive, rMax_Sensitive, dz_Cathode_Conductor / 2.0, phi1, phi2);
  Volume cathodeConductorLog("TPCCathodeConductorLog",cathodeConductorSolid, material_Cathode_Conductor );
  
  // place plus and minus z 
  pv = cathodeLog.placeVolume( cathodeConductorLog , Position( 0.0, 0.0, +(dz_Cathode_Insulator + (dz_Cathode_Conductor / 2.0) ) ) );
  pv = cathodeLog.placeVolume( cathodeConductorLog , Position( 0.0, 0.0, -(dz_Cathode_Insulator + (dz_Cathode_Conductor / 2.0) ) ) );
  //-----------------------------------------------------------------------------------------------//


  //----------------------------------------------- TPC Sensitive Detector (Pad Rings) ---------------------------------------------------------------//
  //fg: fixme: put this to SD and corresponding xml
  //  TPCSD04 *sensitiveDetector = new TPCSD04("TPC", sensitive_threshold_eV);
  //  RegisterSensitiveDetector(sensitiveDetector);
  //  UserLimits *userLimits = new UserLimits(TPCMaxStepLength);

  Tube senstiveGasSolid( rMin_Sensitive, rMax_Sensitive, dz_Sensitive / 2.0, phi1, phi2);

  // ThreeVector translation(0,0,0);
  // RotationMatrix rot;
  // Transform3D transform(rot,translation);

  Volume sensitiveGasLog( "TPCSensitiveLog", senstiveGasSolid, material_TPC_Gas );

  //  sensitiveGasLog->SetVisAttributes(VisAttributes::Invisible);

  // new PVPlacement(Transform3D(RotationMatrix().rotateY(   0 * deg), ThreeVector(0, 0, +( dz_Cathode/2.0 + dz_Sensitive/2.0 ) )), sensitiveGasLog, "TPCSensitiveLog_+z", motherLog, false, 0);
  // new PVPlacement(Transform3D(RotationMatrix().rotateY( 180 * deg), ThreeVector(0, 0, -( dz_Cathode/2.0 + dz_Sensitive/2.0 ) )), sensitiveGasLog, "TPCSensitiveLog_-z", motherLog, false, 1);

  DetElement   sensGasDEfwd( tpc ,  "tpc_sensGas_fwd", x_det.id() );
  DetElement   sensGasDEbwd( tpc ,  "tpc_sensGas_bwd", x_det.id() );

  pv = tpc_motherLog.placeVolume( sensitiveGasLog , Transform3D( RotationY( 0.) , Position(0, 0, +( dz_Cathode/2.0 + dz_Sensitive/2.0 ) ) ) ) ;
  pv.addPhysVolID("side", +1 ) ; 
  sensGasDEfwd.setPlacement( pv ) ;
 
  pv = tpc_motherLog.placeVolume( sensitiveGasLog , Transform3D( RotationY( pi ) , Position(0, 0, -( dz_Cathode/2.0 + dz_Sensitive/2.0 ) ) ) ) ;
  pv.addPhysVolID("side", -1 ) ; 
  sensGasDEbwd.setPlacement( pv ) ;
  

  //debug:  tpc.setVisAttributes(theDetector,  "RedVis" ,  sensitiveGasLog) ;
  tpc.setVisAttributes(theDetector,  "Invisible" ,  sensitiveGasLog) ;

  //---------------------------------------------------- Pad row doublets -------------------------------------------------------------------------------//

  for (int layer = 0; layer < numberPadRows; layer++) {
    
#if 1
    // create twice the number of rings as there are pads, producing an lower and upper part of the pad with the boundry between them the pad-ring centre
    
    const double inner_lowerlayer_radius = rMin_Sensitive + (layer * (padHeight));
    const double outer_lowerlayer_radius = inner_lowerlayer_radius + (padHeight/2.0);
    
    const double inner_upperlayer_radius = outer_lowerlayer_radius ;
    const double outer_upperlayer_radius = inner_upperlayer_radius + (padHeight/2.0);
    
    Tube lowerlayerSolid( inner_lowerlayer_radius, outer_lowerlayer_radius, dz_Sensitive / 2.0, phi1, phi2);
    Tube upperlayerSolid( inner_upperlayer_radius, outer_upperlayer_radius, dz_Sensitive / 2.0, phi1, phi2);

    //fixme: layerstring
    Volume lowerlayerLog( _toString( layer ,"TPC_lowerlayer_log_%02d") ,lowerlayerSolid, material_TPC_Gas );
    Volume upperlayerLog( _toString( layer ,"TPC_upperlayer_log_%02d") ,upperlayerSolid, material_TPC_Gas );

    tpc.setVisAttributes(theDetector,  "Invisible" ,  lowerlayerLog) ;
    tpc.setVisAttributes(theDetector,  "Invisible" ,  upperlayerLog) ;
    

    DetElement   layerDEfwd( sensGasDEfwd ,   _toString( layer, "tpc_row_fwd_%03d") , x_det.id() );
    DetElement   layerDEbwd( sensGasDEbwd ,   _toString( layer, "tpc_row_bwd_%03d") , x_det.id() );
 
    Vector3D o(  inner_upperlayer_radius + 1e-10  , 0. , 0. ) ;
    // create an unbounded surface (i.e. an infinite cylinder) and assign it to the forward gaseous volume only
    VolCylinder surf( upperlayerLog , SurfaceType(SurfaceType::Sensitive, SurfaceType::Invisible, SurfaceType::Unbounded ) ,  (padHeight/2.0) ,  (padHeight/2.0) ,o ) ;

    volSurfaceList( layerDEfwd )->push_back( surf ) ;
//    volSurfaceList( layerDEbwd )->push_back( surf ) ;


    pv = sensitiveGasLog.placeVolume( lowerlayerLog ) ;
    pv.addPhysVolID("layer", layer ).addPhysVolID( "module", 0 ).addPhysVolID("sensor", 1 ) ;

    pv = sensitiveGasLog.placeVolume( upperlayerLog ) ;
    pv.addPhysVolID("layer", layer ).addPhysVolID( "module", 0 ).addPhysVolID("sensor", 0 ) ;
    layerDEfwd.setPlacement( pv ) ;
    layerDEbwd.setPlacement( pv ) ;

    lowerlayerLog.setSensitiveDetector(sens);
    upperlayerLog.setSensitiveDetector(sens);

#else
    // create just one volume per pad ring
    
    const double inner_radius = rMin_Sensitive + (layer * (padHeight) );
    const double outer_radius = inner_radius +  padHeight ;
    
    Tube layerSolid( inner_radius, outer_radius, dz_Sensitive / 2.0, phi1, phi2);

    Volume layerLog( _toString( layer ,"TPC_layer_log_%02d") , layerSolid, material_TPC_Gas );

    tpc.setVisAttributes(theDetector,  "Invisible" ,  layerLog) ;
    
    DetElement   layerDEfwd( sensGasDEfwd ,   _toString( layer, "tpc_row_fwd_%03d") , x_det.id() );
    DetElement   layerDEbwd( sensGasDEbwd ,   _toString( layer, "tpc_row_bwd_%03d") , x_det.id() );
 
    Vector3D o(  inner_radius + (padHeight/2.0)  , 0. , 0. ) ;

    VolCylinder surf( layerLog , SurfaceType(SurfaceType::Sensitive, SurfaceType::Invisible ) ,  (padHeight/2.0) ,  (padHeight/2.0) ,o ) ;

    volSurfaceList( layerDEfwd )->push_back( surf ) ;
    volSurfaceList( layerDEbwd )->push_back( surf ) ;

    pv = sensitiveGasLog.placeVolume( layerLog ) ;
    pv.addPhysVolID("layer", layer  ).addPhysVolID( "module", 0 ) ;

    layerDEfwd.setPlacement( pv ) ;
    layerDEbwd.setPlacement( pv ) ;

    layerLog.setSensitiveDetector(sens);

#endif
  }

  // Assembly of the TPC Readout
  Tube readoutSolid( rMin_GasVolume, rMax_GasVolume, dz_Readout / 2.0, phi1, phi2);
  Volume readoutLog( "TPCReadoutLog", readoutSolid, material_TPC_Gas );
  tpc.setVisAttributes(theDetector,  "CyanVis" ,  readoutLog );


  pv = tpc_motherLog.placeVolume( readoutLog , Transform3D( RotationY( 0.) , Position(0, 0, +( (dz_GasVolume - (dz_Readout/2.0) ) ))) ) ;
  pv = tpc_motherLog.placeVolume( readoutLog , Transform3D( RotationY( pi ) , Position(0, 0, -( (dz_GasVolume - (dz_Readout/2.0) ) ))) ) ;
 
  // new PVPlacement(Transform3D(RotationMatrix().rotateY(  0 * deg), ThreeVector( 0, 0, +(dz_GasVolume - (dz_Readout / 2.0) ))), readoutLog, "TPCReadout_+z", motherLog, false, 0);
  // new PVPlacement(Transform3D(RotationMatrix().rotateY(180 * deg), ThreeVector( 0, 0, -(dz_GasVolume - (dz_Readout / 2.0) ))), readoutLog, "TPCReadout_-z", motherLog, false, 1);

  int pieceCounter = 0;
  double fracRadLengthReadout = 0;
  double zCursor = -dz_Readout / 2;
  
  xml_comp_t  x_ro=  x_det.child( _Unicode( readout ) );
  for(xml_coll_t c( x_ro , _Unicode( row) ); c; ++c)  {
    
    xml_comp_t  x_row( c );
    db = XMLHandlerDB( x_row )  ;

    const double dzPiece = db->fetchDouble("dz") ;
    Material pieceMaterial = theDetector.material( db->fetchString("material") );
    
    Tube pieceSolid(  rMin_GasVolume, rMax_GasVolume, dzPiece / 2, phi1, phi2);
    //fixme namestring
    Volume pieceLog (  _toString( pieceCounter ,"TPCReadoutPieceLog_%02d"), pieceSolid, pieceMaterial ) ;

    //    pieceLog->SetVisAttributes(VisAttributes::Invisible);
    pv = readoutLog.placeVolume( pieceLog  , Position(0, 0,  zCursor + dzPiece/2. ) ) ;
    
    ++pieceCounter;
    fracRadLengthReadout += dzPiece / pieceMaterial->GetMaterial()->GetRadLen();
    zCursor += dzPiece;

    if (zCursor > +dz_Readout / 2) {
      throw GeometryException(  "TPC10: Overfull TPC readout - check your xml file - section <readout>."   ) ;
    }
  }

  // Some verbose output
  cout << "TPC10: Readout material corresponds to " << int(fracRadLengthReadout * 1000) / 10.0 << "% of a radiation length." << endl;

  
  Tube endplateSolid( rInner, rOuter, dz_Endplate / 2.0, phi1, phi2);
  Volume endplateLog( "TPCEndplateLog",endplateSolid, endplate_MaterialMix );
  tpc.setVisAttributes(theDetector,  "CyanVis" ,  endplateLog );
 
  // add a plane to the endcap volume 
  // note: u and v are exchanged: normal is along z ...      
  DetElement   endcapDEfwd( tpc ,  "tpc_endcap_fwd", x_det.id() );
  DetElement   endcapDEbwd( tpc ,  "tpc_endcap_bwd", x_det.id() );

  // vectors for endplate plane
  Vector3D u( 0. , 1. , 0. ) ;
  Vector3D v( 1. , 0. , 0. ) ;
  Vector3D n( 0. , 0. , 1. ) ;

  // need to set the origin of this helper plane to be inside the material ( otherwise it would pick up the vacuum at the origin)
  double mid_r = 0.5 * ( rOuter + rInner ) ;
  Vector3D o( 0. , mid_r , 0. ) ;
  
  VolPlane surf( endplateLog , SurfaceType( SurfaceType::Helper ) ,  dz_Endplate / 2.0 + dz_Readout / 2 ,  dz_Endplate / 2.0 , u , v, n , o ) ;
  volSurfaceList( endcapDEfwd )->push_back( surf ) ;
  volSurfaceList( endcapDEbwd )->push_back( surf ) ;

  // note: as opposed to the readout, the endpate is not placed inside the gas volume
  pv = tpc_motherLog.placeVolume( endplateLog , Transform3D( RotationY( 0. ) , Position(0, 0, +(dz_GasVolume + (dz_Endplate / 2.0))) ) );
  endcapDEfwd.setPlacement( pv ) ;
  pv = tpc_motherLog.placeVolume( endplateLog , Transform3D( RotationY( pi ) , Position(0, 0, -(dz_GasVolume + (dz_Endplate / 2.0))) ) );
  endcapDEbwd.setPlacement( pv ) ;
  
  // new PVPlacement(Transform3D(RotationMatrix().rotateY(  0 * deg), ThreeVector( 0, 0, +(dz_GasVolume + (dz_Endplate / 2.0) ))), endplateLog, "TPCEndplate_+z", motherLog, false, 0);
  // new PVPlacement(Transform3D(RotationMatrix().rotateY(180 * deg), ThreeVector( 0, 0, -(dz_GasVolume + (dz_Endplate / 2.0) ))), endplateLog, "TPCEndplate_-z", motherLog, false, 1);

  cout << "TPC10: Total Endplate material corresponds to " << (fracRadLengthReadout * 100.0) + (dz_Endplate / (endplate_MaterialMix->GetMaterial()->GetRadLen() / mm) * 100.0 ) << "% of a radiation length." << endl;
  

// #ifdef MOKKA_GEAR
//   // save the parameters needed to write the gear file ...
//   _gear_r_min = rInner;
//   _gear_r_max = rOuter;
//   _gear_inner_wall_thickness = dr_InnerWall;
//   _gear_outer_wall_thickness = dr_OuterWall;
//   _gear_r_min_readout = rMin_Sensitive;
//   _gear_r_max_readout = rMin_Sensitive + numberPadRows*padHeight;
//   _gear_n_rows_readout = numberPadRows;
//   _gear_pad_height = padHeight;
//   _gear_pad_width = padWidth;
//   _gear_max_drift_length = dz_Sensitive + dz_Cathode/2.0; // SJA: cathode has to be added as the sensitive region does not start at 0.00    
//   _gear_z_anode = dzTotal - dz_Endplate; // the edge of the readout terminating the drift volume
// #endif
  

  FixedPadSizeTPCData* tpcData = new FixedPadSizeTPCData ;

  tpcData->zHalf = dzTotal/2.0 ; 
  tpcData->rMin = rInner;
  tpcData->rMax = rOuter;
  tpcData->innerWallThickness = dr_InnerWall;
  tpcData->outerWallThickness = dr_OuterWall;
  tpcData->rMinReadout = rMin_Sensitive;
  tpcData->rMaxReadout = rMin_Sensitive + numberPadRows*padHeight;
  tpcData->maxRow = numberPadRows;
  tpcData->padHeight = padHeight;
  tpcData->padWidth = padWidth;
  tpcData->driftLength = dz_Sensitive + dz_Cathode/2.0; // SJA: cathode has to be added as the sensitive region does not start at 0.00    
  tpcData->zMinReadout = dz_Cathode/2.0; 
  //  tpcData.z_anode = dzTotal - dz_Endplate; // the edge of the readout terminating the drift volume
  
  tpc.addExtension< FixedPadSizeTPCData >( tpcData )   ; 

  //######################################################################################################################################################################
  
  
  //--------------------------------------
  
  // Volume mother =  theDetector.pickMotherVolume( tpc ) ;
  // pv = mother.placeVolume(envelope);
  // pv.addPhysVolID( "system", x_det.id() ) ; //.addPhysVolID("side", 0 ) ;
  
  tpc.setVisAttributes( theDetector, x_det.visStr(), envelope );
  //  if( tpc.isValid() ) 
  // tpc.setPlacement(pv);
  
  return tpc;
}
DECLARE_DETELEMENT(TPC10,create_element)

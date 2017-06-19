//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "XMLHandlerDB.h"
#include "XML/Utilities.h"
#include <cmath>

//#include "GearWrapper.h"

using namespace std;
using namespace dd4hep::rec;

/** helper struct */
struct SET_Layer {
  int     n_ladders;
  int     n_sensors_per_ladder;
  double  sensor_length;
  double  half_z;
  double  sensitive_inner_radius ;
  double  support_inner_radius ;
  double  ladder_width ;
  double  ladder_dphi ;
};    

//std::vector<SET_Layer> _SET_Layers;
  
/** helper struct */
struct extended_reconstruction_parameters {
  double sensor_length_mm;
  double strip_width_mm;
  double strip_length_mm;
  double strip_pitch_mm;
  double strip_angle_deg;
};

//extended_reconstruction_parameters _e_r_p;


/** Construction of the SET detector, ported from Mokka driver SET_Simple_Planar.cc
 *
 *  Mokka History:
 *  Feb 7th 2011, Steve Aplin - original version
 *
 *  @author: F.Gaede, DESY, Jan 2014
 */
static dd4hep::Ref_t create_element(dd4hep::Detector& theDetector, xml_h e, dd4hep::SensitiveDetector sens)  {

  //------------------------------------------
  //  See comments starting with '//**' for
  //     hints on porting issues
  //------------------------------------------

  
  xml_det_t    x_det = e;
  string       name  = x_det.nameStr();

  dd4hep::DetElement   set(  name, x_det.id()  ) ;
  
  // --- create an envelope volume and position it into the world ---------------------
  
  dd4hep::Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , set ) ;
  
  dd4hep::xml::setDetectorTypeFlag( e, set ) ;

  if( theDetector.buildType() == dd4hep::BUILD_ENVELOPE ) return set ;
  
  //-----------------------------------------------------------------------------------

  
  //PlacedVolume pv;
  
  
  sens.setType("tracker");

  
  dd4hep::rec::ZPlanarData*  zPlanarData = new ZPlanarData ;

  //######################################################################################################################################################################
  //  code ported from SET_Simple_Planar::construct() :
  //##################################

  // extended_reconstruction_parameters _e_r_p;
  
  // *********************
  //  Read and Store the Extended Reconstruction Parameters which are passed directly through to gear. Note others may be added below
  // db->exec("select * from extended_reconstruction_parameters;");
  // db->getTuple();
  XMLHandlerDB db = XMLHandlerDB(  x_det.child( _Unicode( reconstruction ) ) ) ;
  
  zPlanarData->widthStrip  = db->fetchDouble("strip_width")  ;
  zPlanarData->lengthStrip = db->fetchDouble("strip_length") ;
  zPlanarData->pitchStrip  = db->fetchDouble("strip_pitch")  ;
  zPlanarData->angleStrip  = db->fetchDouble("strip_angle") ;
  double strip_angle = zPlanarData->angleStrip ;
  // *********************
  
  
  //... db common_parameters
  // // db->exec("select * from global;");
  // // db->getTuple();
  db = XMLHandlerDB(  x_det.child( _Unicode( global ) ) ) ;

  // Sensitive Thickness  
  double sensitive_thickness = db->fetchDouble("sensitive_thickness") ;
  // Support Thickness
  double support_thickness = db->fetchDouble("support_thickness") ;
  // Sensor Length
  double sensor_length = db->fetchDouble("sensor_length") ;
  
  // _e_r_p.sensor_length_mm  =sensor_length;

  dd4hep::Material air = theDetector.air()  ;
  dd4hep::Material sensitiveMat = theDetector.material(db->fetchString("sensitive_mat"));  
  dd4hep::Material supportMat   = theDetector.material(db->fetchString("support_mat"));  
  
  
  // // // setup the encoder 
  // // UTIL::BitField64 encoder( LCTrackerCellID::encoding_string() ) ; 
  
  // // encoder.reset() ;  // reset to 0
  
  // // encoder[LCTrackerCellID::subdet()] = ILDDetID::NOTUSED ;
  // // encoder[LCTrackerCellID::side()] = 0 ;
  // // encoder[LCTrackerCellID::layer()]  = 0 ;
  // // encoder[LCTrackerCellID::module()] = 0 ;
  // // encoder[LCTrackerCellID::sensor()] = 0 ;
  // // int cellID0 = encoder.lowWord() ;
  
  //... The SET Sensitive detector
  //unused:  double sensitive_threshold_KeV = db->fetchDouble("sensitive_threshold_KeV")  ;
  
  //FIXME: the SD  ...
  // // _theSETSD = 
  // // new TRKSD02("SET",
  // //             _sensitive_thickness * mm 
  // //             * sensitive_threshold_KeV ,
  // //             10.0 * MeV);
  
  // // RegisterSensitiveDetector(_theSETSD);
  

  const double TPC_outer_radius = theDetector.constant<double>("TPC_outer_radius");
  const double TPC_Ecal_Hcal_barrel_halfZ = theDetector.constant<double>("TPC_Ecal_Hcal_barrel_halfZ");

  for(xml_coll_t c( x_det ,_U(layer)); c; ++c)  {
    
    xml_comp_t  x_layer( c );
    db = XMLHandlerDB( x_layer )  ;
    
    int layer_id = db->fetchInt("layer_id");
    
    double sensitive_distance_from_tpc = db->fetchDouble("sensitive_distance_from_tpc") ;
    double coverage_of_TPC_Ecal_Hcal_barrel = db->fetchDouble("coverage_of_TPC_Ecal_Hcal_barrel") ;
    
    double max_half_z = coverage_of_TPC_Ecal_Hcal_barrel * TPC_Ecal_Hcal_barrel_halfZ;
        
    int number_of_sensors_per_half = floor(max_half_z/sensor_length);
        
    double half_z = sensor_length * number_of_sensors_per_half;
        
    double sensitive_radius = TPC_outer_radius + sensitive_distance_from_tpc;
        
    int n_ladders        = db->fetchInt("n_ladders") ;
    int faces_IP         = db->fetchInt("faces_IP") ;
    double ladder_clearance = db->fetchDouble("ladder_clearance") ;

        
    // create assembly and DetElement for the layer
    std::string layerName = dd4hep::_toString( layer_id , "layer_%d" );
    dd4hep::Assembly layer_assembly( layerName ) ;
    dd4hep::PlacedVolume pv = envelope.placeVolume( layer_assembly ) ;
    dd4hep::DetElement layerDE( set , layerName  , x_det.id() );
    layerDE.setPlacement( pv ) ;


    const double ladder_dphi = ( dd4hep::twopi / n_ladders ) ;

    double sensitive_inner_radius = sensitive_radius - 0.5 * sensitive_thickness;
    double ladder_width = 2*(tan(ladder_dphi*0.5)*sensitive_inner_radius - ladder_clearance) ;
        
    double support_radius(0) ;

    if( faces_IP == 1 ){ // support is on the outside 
      support_radius = sensitive_radius + (0.5 * sensitive_thickness) ;
      ladder_width = 2*(tan(ladder_dphi*0.5)*sensitive_inner_radius - ladder_clearance) ;
    }
    else{ // support is on the inside
      support_radius = sensitive_radius - (0.5 * sensitive_thickness) - support_thickness;
      ladder_width = 2*(tan(ladder_dphi*0.5)*support_radius - ladder_clearance) ;
    }


    //FIXME: GEAR....
    // std::ostringstream ossradius;
    // std::ostringstream osshalfz;
    // ossradius << inner_most_radius / mm;
    // osshalfz << half_z / mm;

    // if(is_SET1 == 1){
    //   (*Control::globalModelParameters)["SET1_Radius"] = ossradius.str();
    //   (*Control::globalModelParameters)["SET1_Half_Length_Z"] = osshalfz.str();
    // }
    // if(is_SET2 == 1){
    //   (*Control::globalModelParameters)["SET2_Radius"] = ossradius.str();
    //   (*Control::globalModelParameters)["SET2_Half_Length_Z"] = osshalfz.str();
    // }
    
    dd4hep::rec::ZPlanarData::LayerLayout thisLayer ;
    thisLayer.sensorsPerLadder = number_of_sensors_per_half * 2.0 ;
    thisLayer.lengthSensor     = sensor_length ;
    
    thisLayer.distanceSupport  = support_radius;
    thisLayer.offsetSupport    = 0. ;
    thisLayer.thicknessSupport = support_thickness ;
    thisLayer.zHalfSupport     = half_z ;
    thisLayer.widthSupport     = ladder_width ;
    
    thisLayer.distanceSensitive  = sensitive_radius - 0.5 *sensitive_thickness;
    thisLayer.offsetSensitive    = 0. ;
    thisLayer.thicknessSensitive = sensitive_thickness ;
    thisLayer.zHalfSensitive     = half_z ;
    thisLayer.widthSensitive     = ladder_width ;
    
    thisLayer.ladderNumber =  n_ladders;
    thisLayer.phi0         =  0. ;
    
    zPlanarData->layers.push_back( thisLayer ) ;
    
    SET_Layer layer_geom ;
    std::vector<SET_Layer> SET_Layers;
    
    layer_geom.n_ladders = n_ladders;
    layer_geom.half_z = half_z ;
    layer_geom.sensor_length = sensor_length;
    layer_geom.n_sensors_per_ladder = number_of_sensors_per_half * 2.0;
    layer_geom.sensitive_inner_radius = sensitive_radius - 0.5 * sensitive_thickness;
    layer_geom.support_inner_radius = support_radius;
    layer_geom.ladder_width = ladder_width ;
    layer_geom.ladder_dphi = ladder_dphi;
    
    SET_Layers.push_back(layer_geom) ;
    
    
    std::cout << "SET_Simple_Planar: Layer:" << layer_id
	      << "\t half length = " << layer_geom.half_z
	      << "\t sensor length = " << layer_geom.sensor_length
	      << "\t n sensors per ladder = " << layer_geom.n_sensors_per_ladder
	      << "\t min r sensitive = " << layer_geom.sensitive_inner_radius
	      << "\t min r support = " << layer_geom.support_inner_radius
	      << "\t n ladders = " << layer_geom.n_ladders
	      << "\t ladder width = " << layer_geom.ladder_width
	      << "\t ladder clearance = " << ladder_clearance
	      << "\t ladder dphi = " << ladder_dphi
	      << "\t sensitive mat = " <<sensitiveMat->GetName()
	      << "\t support mat = " <<supportMat->GetName()
	      << "\t faces_IP = " << faces_IP
	      << std::endl;
    
        
    
    //************************************************************************************************
    // Geometric Values Established. Start Creating Volumes.
        

    // create an enclosing ladder volume that will be placed in the world volume for every ladder
        
    dd4hep::Box setLadderSolid( (sensitive_thickness +support_thickness ) / 2.0 ,
                                layer_geom.ladder_width / 2.0,
                                layer_geom.half_z);

    dd4hep::Volume setLadderLogical (dd4hep::_toString( layer_id,"SET_LadderLogical_%02d"), setLadderSolid, air ) ; 
        
    // now create an envelope volume to represent the sensitive area, which will be divided up into individual sensors         
        
    dd4hep::Box setSenEnvelopeSolid( (sensitive_thickness ) / 2.0 ,
                                     layer_geom.ladder_width  / 2.0,
                                     layer_geom.half_z);
    
    //fixme: material ???    Volume setSenEnvelopeLogical( _toString( layer_id,"SET_SenEnvelopeLogical_%02d"), setSenEnvelopeSolid, sensitiveMat )  ;
    dd4hep::Volume setSenEnvelopeLogical( dd4hep::_toString( layer_id,"SET_SenEnvelopeLogical_%02d"),
                                          setSenEnvelopeSolid, air )  ;
    
    // create the sensor volumes and place them in the senstive envelope volume 
    
    dd4hep::Box setSenSolid( (sensitive_thickness ) / 2.0 ,
                             layer_geom.ladder_width  / 2.0,
                             (layer_geom.sensor_length / 2.0 ) - 1.e-06*dd4hep::mm ); // added tolerance to avoid false overlap detection
    
    dd4hep::Volume setSenLogical( dd4hep:: _toString( layer_id,"SET_SenLogical_%02d"), setSenSolid,sensitiveMat ) ; 
    
    setSenLogical.setSensitiveDetector(sens);
    
    
    
    //====== create the meassurement surface ===================
    Vector3D u,v,n ;
    
    if( faces_IP == 0 ){

      n.fill( -1. ,   0. , 0. ) ;

      // implement stereo angle 
      u.fill( 0. , -cos( strip_angle  ) , -sin( strip_angle  ) ) ;
      v.fill( 0. , -sin( strip_angle  ) ,  cos( strip_angle  ) ) ;

    } else {

      n.fill( 1. , 0. , 0. ) ;

      // implement stereo angle 
      u.fill( 0. ,  cos( strip_angle  ) ,  sin( strip_angle  ) ) ;
      v.fill( 0. , -sin( strip_angle  ) ,  cos( strip_angle  ) ) ;
    }

    double inner_thick =  sensitive_thickness / 2.0 ;
    double outer_thick =  sensitive_thickness / 2.0 + support_thickness ;  // support is on top
   
    VolPlane surf( setSenLogical , SurfaceType(SurfaceType::Sensitive,SurfaceType::Measurement1D) ,inner_thick, outer_thick , u,v,n ) ; //,o ) ;
 
    // vector of sensor placements - needed for DetElements in ladder loop below
    std::vector<dd4hep::PlacedVolume> pvV(  layer_geom.n_sensors_per_ladder ) ;

  //============================================================

    for (int isensor=0; isensor < layer_geom.n_sensors_per_ladder ; ++isensor) {
      
      // encoder.reset() ;  // reset to 0
      // encoder[LCTrackerCellID::subdet()] = ILDDetID::NOTUSED ;
      // encoder[LCTrackerCellID::sensor()] =  isensor+1;
      // cellID0 = encoder.lowWord() ;
      
      double xpos = 0.0;
      double ypos = 0.0;
      double zpos = -layer_geom.half_z + (0.5*layer_geom.sensor_length) + (isensor*layer_geom.sensor_length) ;
      
      pv = setSenEnvelopeLogical.placeVolume( setSenLogical,
                                              dd4hep::Transform3D( dd4hep::RotationY(0.) ,
                                                                   dd4hep::Position( xpos, ypos, zpos)  ) );
      
      pv.addPhysVolID("sensor",  isensor ) ; 
      //fixme: what is the correct numbering convention ?
      // pv.addPhysVolID("sensor",  isensor + 1 ) ; 
      pvV[isensor] = pv ;
   }					      
    
    set.setVisAttributes(theDetector, "SeeThrough",  setLadderLogical ) ;
    set.setVisAttributes(theDetector, "SeeThrough",  setSenEnvelopeLogical ) ;

    set.setVisAttributes(theDetector, "BlueVis",       setSenLogical ) ;
    
    
    // encoder.reset() ;  // reset to 0
    // encoder[LCTrackerCellID::subdet()] = ILDDetID::NOTUSED ;
    // encoder[LCTrackerCellID::layer()]  = layer_id ;
    // cellID0 = encoder.lowWord() ;
        

    pv = setLadderLogical.placeVolume( setSenEnvelopeLogical ,
                                       dd4hep::Transform3D( dd4hep::RotationY( 0.), 
                                                            dd4hep::Position( (-(sensitive_thickness +support_thickness ) / 2.0
                                                                               + ( sensitive_thickness / 2.0) ), 0.,0.) ) );
    // pv = setSenEnvelopeLogical.placeVolume( setLadderLogical, Transform3D( RotationY( 0.), 
    // 									   Position( (-(sensitive_thickness +support_thickness ) / 2.0 + ( sensitive_thickness / 2.0) ), 0.,0.) ) );

    //fixme: needed ??    pv.addPhysVolID("layer", layer_id ) ; 
    
    

    // create support volume which will be placed in the enclosing ladder volume together with the senstive envelope volume
    
    dd4hep::Box setSupSolid( (support_thickness ) / 2.0 ,
                             layer_geom.ladder_width / 2.0,
                             layer_geom.half_z);
    
    dd4hep::Volume setSupLogical( dd4hep::_toString( layer_id,"SET_SupLogical_%02d"),  setSupSolid, supportMat ) ;
    
    
    set.setVisAttributes(theDetector, "RedVis",  setSupLogical ) ;
    
    
    pv = setLadderLogical.placeVolume( setSupLogical,
                                       dd4hep::Transform3D( dd4hep::RotationY( 0.), 
                                                            dd4hep::Position( (-(sensitive_thickness +support_thickness ) / 2.0
                                                                               +sensitive_thickness + ( support_thickness / 2.0)   ), 0.,0.) ) );
    
    for( int i = 0 ; i < n_ladders ; ++i ){
      
      std::stringstream ladder_enum; ladder_enum << "set_ladder_" << layer_id << "_" << i;
      
      dd4hep::DetElement ladderDE( layerDE ,  ladder_enum.str() , x_det.id() );

      for (int isensor=0; isensor < layer_geom.n_sensors_per_ladder ; ++isensor) {

	std::stringstream sensor_ss ;  sensor_ss << ladder_enum.str() << "_" << isensor ;
	
	dd4hep::DetElement sensorDE( ladderDE, sensor_ss.str() ,  x_det.id() );
	sensorDE.setPlacement( pvV[isensor] ) ;

	volSurfaceList( sensorDE )->push_back(  surf ) ;
      }					      
    

     // RotationMatrix *rot = new RotationMatrix();
      // rot->rotateZ( i * -ladder_dphi );
      
      // // rotate by 180 degrees around z if facing away from the IP
      // if( faces_IP == 0 ) rot->rotateZ( 180 * deg );
      
      // encoder[LCTrackerCellID::subdet()] = ILDDetID::SET ;
      // encoder[LCTrackerCellID::layer()]  = layer_id ;
      // encoder[LCTrackerCellID::module()] = i + 1 ;
      // cellID0 = encoder.lowWord() ;  
      
      float dr = ( (sensitive_thickness +support_thickness ) / 2.0 ) - (sensitive_thickness / 2.0 ) ;
      
      //      double phi_rot =  i * -ladder_dphi ;
      double phi_rot =  i * ladder_dphi ;

      if( faces_IP == 0 ) { 

	dr = -dr;

	phi_rot += M_PI ;
      }

      pv = layer_assembly.placeVolume( setLadderLogical,
                                       dd4hep::Transform3D( dd4hep::RotationZYX(  phi_rot, 0. , 0. ), 
                                                            dd4hep::Position( (sensitive_radius+dr) * cos(i * ladder_dphi), 
                                                                              (sensitive_radius+dr) * sin(i * ladder_dphi), 
                                                                              0. ) ) ) ;
      
      pv.addPhysVolID("layer", layer_id ).addPhysVolID("module", i ) ; 
      //fixme: what is the correct numbering convention ?
      //pv.addPhysVolID("layer", layer_id ).addPhysVolID("module", i+1 ) ; 
      
      ladderDE.setPlacement( pv ) ;
    }
    
    
  }

  cout << "SET_Simple_Planar done.\n" << endl;
  //######################################################################################################################################################################
  
  set.addExtension< ZPlanarData >( zPlanarData ) ;

  //--------------------------------------
  
  
  set.setVisAttributes( theDetector, x_det.visStr(), envelope );
  
  return set;
}
DECLARE_DETELEMENT(SET_Simple_Planar,create_element)


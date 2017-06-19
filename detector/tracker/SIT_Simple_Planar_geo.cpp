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
#include "XML/Utilities.h"
#include "XMLHandlerDB.h"
#include <cmath>

using namespace std;

using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Position;
using dd4hep::RotationY;
using dd4hep::RotationZYX;
using dd4hep::Transform3D;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::rec::volSurfaceList;
using dd4hep::rec::ZPlanarData;

/** helper struct */
struct SIT_Layer {
  int     n_ladders;
  int     n_sensors_per_ladder;
  double  sensor_length;
  double  half_z;
  double  sensitive_inner_radius ;
  double  support_inner_radius ;
  double  ladder_width ;
  double  ladder_dphi ;
};    

//std::vector<SIT_Layer> _SIT_Layers;
  
// /** helper struct */
// struct extended_reconstruction_parameters {
//   double sensor_length_mm;
//   double strip_width_mm;
//   double strip_length_mm;
//   double strip_pitch_mm;
//   double strip_angle_deg;
// };

//extended_reconstruction_parameters _e_r_p;


/** Construction of the SIT detector, ported from Mokka driver SIT_Simple_Planar.cc
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
  
  dd4hep::DetElement   sit(  name, x_det.id()  ) ;
  
  // --- create an envelope volume and position it into the world ---------------------
  
  dd4hep::Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sit ) ;
  
  dd4hep::xml::setDetectorTypeFlag( e, sit ) ;

  if( theDetector.buildType() == dd4hep::BUILD_ENVELOPE ) return sit ;
  
  //-----------------------------------------------------------------------------------
  
  dd4hep::PlacedVolume pv;
  
  
  sens.setType("tracker");


  dd4hep::rec::ZPlanarData*  zPlanarData = new dd4hep::rec::ZPlanarData ;

  //######################################################################################################################################################################
  //  code ported from SIT_Simple_Planar::construct() :
  //##################################

  //  extended_reconstruction_parameters _e_r_p;

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
  
  //... The SIT Sensitive detector
  //unused:  double sensitive_threshold_KeV = db->fetchDouble("sensitive_threshold_KeV")  ;
  
  //FIXME: the SD  ...
  // // _theSITSD = 
  // // new TRKSD02("SIT",
  // //             _sensitive_thickness * mm 
  // //             * sensitive_threshold_KeV ,
  // //             10.0 * MeV);
  
  // // RegisterSensitiveDetector(_theSITSD);
  

  for(xml_coll_t c( x_det ,_U(layer)); c; ++c)  {
    
    xml_comp_t  x_layer( c );
    db = XMLHandlerDB( x_layer )  ;
    
    int layer_id = db->fetchInt("layer_id");
    
    double half_z(0);
    double sensitive_radius(0);
    double sensitive_inner_radius(0);
    double support_radius(0);
    int    n_sensors_per_ladder(0) ;
    int    n_ladders(0) ;
    double ladder_width(0) ;
    double ladder_clearance(0) ;
    int    faces_IP(0) ;
    int    is_SIT1(0) ;
    int    is_SIT2(0) ;

      
    sensitive_radius     = db->fetchDouble("sensitive_radius") ;
    n_sensors_per_ladder = db->fetchInt("n_sensors_per_ladder") ;
    half_z               = (n_sensors_per_ladder *sensor_length) / 2.0 ;
    n_ladders            = db->fetchInt("n_ladders") ;
    faces_IP             = db->fetchInt("faces_IP") ;
    is_SIT1              = db->fetchInt("is_SIT1") ;
    is_SIT2              = db->fetchInt("is_SIT2") ;
    ladder_clearance     = db->fetchDouble("ladder_clearance") ;


    // create assembly and DetElement for the layer
    std::string layerName = dd4hep::_toString( layer_id , "layer_%d"  );
    dd4hep::Assembly layer_assembly( layerName ) ;
    pv = envelope.placeVolume( layer_assembly ) ;
    dd4hep::DetElement layerDE( sit , layerName  , x_det.id() );
    layerDE.setPlacement( pv ) ;


    const double ladder_dphi = ( dd4hep::twopi / n_ladders ) ;

    sensitive_inner_radius = sensitive_radius - 0.5 *sensitive_thickness;
    ladder_width = 2*(tan(ladder_dphi*0.5)*sensitive_inner_radius - ladder_clearance) ;
                    
    // double inner_most_radius = 0.0;
    
    if( faces_IP == 1 ){ // support is on the outside 
      support_radius = sensitive_radius + (0.5 *sensitive_thickness) ;
      ladder_width = 2*(tan(ladder_dphi*0.5)*sensitive_inner_radius - ladder_clearance) ;
      // inner_most_radius = sensitive_inner_radius;
    }
    else{ // support is on the inside
      support_radius = sensitive_radius - (0.5 *sensitive_thickness) -support_thickness;
      ladder_width = 2*(tan(ladder_dphi*0.5)*support_radius - ladder_clearance) ;
      // inner_most_radius = support_radius;
    }
    
    //FIXME: GEAR....
    // std::ostringstream ossradius;
    // std::ostringstream osshalfz;
    // ossradius << inner_most_radius / mm;
    // osshalfz << half_z / mm;

    // if(is_SIT1 == 1){
    //   (*Control::globalModelParameters)["SIT1_Radius"] = ossradius.str();
    //   (*Control::globalModelParameters)["SIT1_Half_Length_Z"] = osshalfz.str();
    // }
    // if(is_SIT2 == 1){
    //   (*Control::globalModelParameters)["SIT2_Radius"] = ossradius.str();
    //   (*Control::globalModelParameters)["SIT2_Half_Length_Z"] = osshalfz.str();
    // }
    
    dd4hep::rec::ZPlanarData::LayerLayout thisLayer ;
    thisLayer.sensorsPerLadder = n_sensors_per_ladder ;
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

    SIT_Layer layer_geom ;
    layer_geom.n_ladders = n_ladders;
    layer_geom.sensor_length =sensor_length;
    layer_geom.n_sensors_per_ladder = n_sensors_per_ladder;
    layer_geom.half_z = half_z ;
    layer_geom.sensitive_inner_radius = sensitive_radius - 0.5 *sensitive_thickness;
    layer_geom.support_inner_radius = support_radius;
    layer_geom.ladder_width = ladder_width ;
    layer_geom.ladder_dphi = ladder_dphi;
    std::vector<SIT_Layer>SIT_Layers;
    SIT_Layers.push_back(layer_geom) ;
    
    
    std::cout << "SIT_Simple_Planar: Layer:" << layer_id
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
    	      << "\t is_SIT1 = " << is_SIT1
    	      << "\t is_SIT2 = " << is_SIT2
    	      << std::endl;
    
        
    // std::stringstream name_base;
    // name_base << "SIT";
    // std::stringstream name_enum;
    // name_enum << layer_id;
        
    // create an enclosing ladder volume that will be placed in the world volume for every ladder
        
    dd4hep::Box sitLadderSolid( (sensitive_thickness +support_thickness ) / 2.0 ,
                                layer_geom.ladder_width / 2.0,
                                layer_geom.half_z);

    dd4hep::Volume sitLadderLogical( dd4hep::_toString( layer_id,"SIT_LadderLogical_%02d"), sitLadderSolid, air ) ; 
        
    // now create an envelope volume to represent the sensitive area, which will be divided up into individual sensors         
        
    dd4hep::Box sitSenEnvelopeSolid( (sensitive_thickness ) / 2.0 ,
                                     layer_geom.ladder_width  / 2.0,
                                     layer_geom.half_z);
    
    //fixme: material ???    Volume sitSenEnvelopeLogical( _toString( layer_id,"SIT_SenEnvelopeLogical_%02d"), sitSenEnvelopeSolid, sensitiveMat )  ;
    dd4hep::Volume sitSenEnvelopeLogical( dd4hep::_toString( layer_id,"SIT_SenEnvelopeLogical_%02d"),
                                          sitSenEnvelopeSolid, air )  ;
    
    // create the sensor volumes and place them in the senstive envelope volume 
    
    dd4hep::Box sitSenSolid( (sensitive_thickness ) / 2.0 ,
                             layer_geom.ladder_width  / 2.0,
                             (layer_geom.sensor_length / 2.0 ) - 1.e-06 * dd4hep::mm ); // added tolerance to avoid false overlap detection
    
    dd4hep::Volume sitSenLogical( dd4hep::_toString( layer_id,"SIT_SenLogical_%02d"), sitSenSolid,sensitiveMat ) ; 
    
    sitSenLogical.setSensitiveDetector(sens);


    //====== create the meassurement surface ===================
    dd4hep::rec::Vector3D u,v,n ;

    if( faces_IP == 0 ){

      n.fill( -1. ,   0. , 0. ) ;

      // implement stereo angle 
      u.fill( 0. ,  -cos( strip_angle ) , -sin( strip_angle  ) ) ;
      v.fill( 0. ,  -sin( strip_angle ) ,  cos( strip_angle  ) ) ;

    } else {

      n.fill( 1. , 0. , 0. ) ;

      // implement stereo angle 
      u.fill( 0. ,  cos( strip_angle  ) ,  sin( strip_angle  ) ) ;
      v.fill( 0. , -sin( strip_angle  ) ,  cos( strip_angle  ) ) ;
    }


    double inner_thick =  sensitive_thickness / 2.0 ;
    double outer_thick =  sensitive_thickness / 2.0 ;
   
    if( faces_IP ){
      outer_thick += support_thickness ;
    } else { 
      inner_thick += support_thickness ;
    }


    dd4hep::rec::VolPlane surf( sitSenLogical ,
                                dd4hep::rec::SurfaceType(dd4hep::rec::SurfaceType::Sensitive,
                                                         dd4hep::rec::SurfaceType::Measurement1D),
                                inner_thick, outer_thick , u,v,n ) ; //,o ) ;

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
      
      pv = sitSenEnvelopeLogical.placeVolume( sitSenLogical,
                                              Transform3D( RotationY(0.) ,
                                                           Position( xpos, ypos, zpos)  ) );
      
      pv.addPhysVolID("sensor",  isensor ) ; 
      //fixme: what is the correct numbering convention ?
      // pv.addPhysVolID("sensor",  isensor + 1 ) ; 
      pvV[isensor] = pv ;
    }					      
    
    sit.setVisAttributes(theDetector, "SeeThrough",  sitLadderLogical ) ;
    sit.setVisAttributes(theDetector, "SeeThrough",  sitSenEnvelopeLogical ) ;

    sit.setVisAttributes(theDetector, "BlueVis",       sitSenLogical ) ;
    
    
    // encoder.reset() ;  // reset to 0
    // encoder[LCTrackerCellID::subdet()] = ILDDetID::NOTUSED ;
    // encoder[LCTrackerCellID::layer()]  = layer_id ;
    // cellID0 = encoder.lowWord() ;
        

    pv = sitLadderLogical.placeVolume( sitSenEnvelopeLogical , Transform3D( RotationY( 0.), 
									   Position( (-(sensitive_thickness +support_thickness ) / 2.0 + ( sensitive_thickness / 2.0) ), 0.,0.) ) );
    // pv = sitSenEnvelopeLogical.placeVolume( sitLadderLogical, Transform3D( RotationY( 0.), 
    // 									   Position( (-(sensitive_thickness +support_thickness ) / 2.0 + ( sensitive_thickness / 2.0) ), 0.,0.) ) );

    //fixme: needed ??    pv.addPhysVolID("layer", layer_id ) ; 
    
    

    // create support volume which will be placed in the enclosing ladder volume together with the senstive envelope volume
    
    Box sitSupSolid( (support_thickness ) / 2.0 ,
		     layer_geom.ladder_width / 2.0,
		     layer_geom.half_z);
    
    Volume sitSupLogical(   _toString( layer_id,"SIT_SupLogical_%02d"),  sitSupSolid, supportMat ) ;
    
    
    sit.setVisAttributes(theDetector, "RedVis",  sitSupLogical ) ;
    
    
    pv = sitLadderLogical.placeVolume( sitSupLogical, Transform3D( RotationY( 0.), 
								   Position( (-(sensitive_thickness +support_thickness ) / 2.0 +sensitive_thickness + ( support_thickness / 2.0)   ), 0.,0.) ) );
    
    for( int i = 0 ; i < n_ladders ; ++i ){
      
      std::stringstream ladder_enum; ladder_enum << "sit_ladder_" << layer_id << "_" << i;
      

      DetElement   ladderDE( layerDE ,  ladder_enum.str() , x_det.id() );

      for (int isensor=0; isensor < layer_geom.n_sensors_per_ladder ; ++isensor) {

	std::stringstream sensor_ss ;  sensor_ss << ladder_enum.str() << "_" << isensor ;
	
	DetElement sensorDE( ladderDE, sensor_ss.str() ,  x_det.id() );
	sensorDE.setPlacement( pvV[isensor] ) ;

	volSurfaceList( sensorDE )->push_back(  surf ) ;
      }					      
    

      // RotationMatrix *rot = new RotationMatrix();
      // rot->rotateZ( i * -ladder_dphi );
      
      // // rotate by 180 degrees around z if facing away from the IP
      // if( faces_IP == 0 ) rot->rotateZ( 180 * deg );
      
      // encoder[LCTrackerCellID::subdet()] = ILDDetID::SIT ;
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

      pv = layer_assembly.placeVolume( sitLadderLogical, Transform3D( RotationZYX(  phi_rot, 0. , 0. ), 
								      Position( (sensitive_radius+dr) * cos(i * ladder_dphi), 
										(sensitive_radius+dr) * sin(i * ladder_dphi), 
										0. ) ) ) ;
      
      pv.addPhysVolID("layer", layer_id ).addPhysVolID("module", i ) ; 
      //fixme: what is the correct numbering convention ?
      //pv.addPhysVolID("layer", layer_id ).addPhysVolID("module", i+1 ) ; 
 

      ladderDE.setPlacement( pv ) ;
    }
    
    
  }

  cout << "SIT_Simple_Planar done.\n" << endl;
  //######################################################################################################################################################################
  
  sit.addExtension< ZPlanarData >( zPlanarData ) ;
  
  //--------------------------------------
  
  sit.setVisAttributes( theDetector, x_det.visStr(), envelope );
  
  return sit;
}
DECLARE_DETELEMENT(SIT_Simple_Planar,create_element)

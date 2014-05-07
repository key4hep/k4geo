//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id:$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/TGeoUnits.h"
#include "DDRec/Surface.h"
#include <cmath>

//#include "GearWrapper.h"

using namespace std;
using namespace DD4hep;
using namespace tgeo ;
using namespace DD4hep::Geometry;
using namespace DDRec ;


/** Wrapper class to replace the Database class used in Mokka to read the parameters.
 *  Assumes parameters are stored as attributes of the corresponding xml element.
 */
struct XMLHandlerDB{
  xml_comp_t x_det ;
  /** C'tor initializes the handle */
  XMLHandlerDB(xml_comp_t det) : x_det(det) {}
  
  double fetchDouble( const char* _name){ return  x_det.attr<double>( _name )  ; }

  int    fetchInt( const char* _name){ return  x_det.attr<int>( _name ) ; }

  std::string fetchString( const char* _name){ return  x_det.attr<string>( _name ) ;}

  /** allow this to be used as a 'pointer' ( as was used for Mokka Database object)*/
  XMLHandlerDB* operator->() { return this ; }
};


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
  
/** helper struct */
struct extended_reconstruction_parameters {
  double sensor_length_mm;
  double strip_width_mm;
  double strip_length_mm;
  double strip_pitch_mm;
  double strip_angle_deg;
};

//extended_reconstruction_parameters _e_r_p;


/** Construction of the SIT detector, ported from Mokka driver SIT_Simple_Planar.cc
 *
 *  Mokka History:
 *  Feb 7th 2011, Steve Aplin - original version
 *
 *  @author: F.Gaede, DESY, Jan 2014
 */
static Ref_t create_element(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {

  //------------------------------------------
  //  See comments starting with '//**' for
  //     hints on porting issues
  //------------------------------------------

  
  xml_det_t    x_det = e;
  string       name  = x_det.nameStr();
  
  //---- envelope: cylinder of air:
  //xml_comp_t  x_tube (x_det.child(_U(tubs)));
  //Tube        envelope_cylinder( x_tube.rmin(), x_tube.rmax(), x_tube.zhalf() );
  //Volume      envelope( "sit_envelope_cyl", envelope_cylinder , lcdd.air() );
  //--------------------------------
  Assembly envelope( name + "_assembly"  ) ;
  //--------------------------------
  
  PlacedVolume pv;
  
  DetElement   sit(  name, x_det.id()  ) ;
  
  sens.setType("tracker");

  //######################################################################################################################################################################
  //  code ported from SIT_Simple_Planar::construct() :
  //##################################

  extended_reconstruction_parameters _e_r_p;

  // *********************
  //  Read and Store the Extended Reconstruction Parameters which are passed directly through to gear. Note others may be added below
  // db->exec("select * from extended_reconstruction_parameters;");
  // db->getTuple();
  
  // _e_r_p.strip_width_mm  = db->fetchDouble("strip_width_mm")  *mm;
  // _e_r_p.strip_length_mm = db->fetchDouble("strip_length_mm") *mm;
  // _e_r_p.strip_pitch_mm  = db->fetchDouble("strip_pitch_mm")  *mm;
  // _e_r_p.strip_angle_deg = db->fetchDouble("strip_angle_deg") *deg;
  
  // *********************
  
  
  //... db common_parameters
  // // db->exec("select * from global;");
  // // db->getTuple();
  XMLHandlerDB db = XMLHandlerDB(  x_det.child( _Unicode( global ) ) ) ;

  // Sensitive Thickness  
  double sensitive_thickness = db->fetchDouble("sensitive_thickness") ;
  // Support Thickness
  double support_thickness = db->fetchDouble("support_thickness") ;
  // Sensor Length
  double sensor_length = db->fetchDouble("sensor_length") ;
  
  Material air = lcdd.air()  ;
  Material sensitiveMat = lcdd.material(db->fetchString("sensitive_mat"));  
  Material supportMat   = lcdd.material(db->fetchString("support_mat"));  
  
  
  // // // setup the encoder 
  // // UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  
  // // encoder.reset() ;  // reset to 0
  
  // // encoder[ILDCellID0::subdet] = ILDDetID::NOTUSED ;
  // // encoder[ILDCellID0::side] = 0 ;
  // // encoder[ILDCellID0::layer]  = 0 ;
  // // encoder[ILDCellID0::module] = 0 ;
  // // encoder[ILDCellID0::sensor] = 0 ;
  // // int cellID0 = encoder.lowWord() ;
  
  //... The SIT Sensitive detector
  double sensitive_threshold_KeV = db->fetchDouble("sensitive_threshold_KeV")  ;
  
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
    std::string layerName = _toString( layer_id , "layer_%d"  );
    Assembly layer_assembly( layerName ) ;
    PlacedVolume pv = envelope.placeVolume( layer_assembly ) ;
    DetElement layerDE( sit , layerName  , x_det.id() );
    layerDE.setPlacement( pv ) ;


    const double ladder_dphi = ( twopi / n_ladders ) ;

    sensitive_inner_radius = sensitive_radius - 0.5 *sensitive_thickness;
    ladder_width = 2*(tan(ladder_dphi*0.5)*sensitive_inner_radius - ladder_clearance) ;
                    
    double inner_most_radius = 0.0;
    
    if( faces_IP == 1 ){ // support is on the outside 
      support_radius = sensitive_radius + (0.5 *sensitive_thickness) ;
      ladder_width = 2*(tan(ladder_dphi*0.5)*sensitive_inner_radius - ladder_clearance) ;
      inner_most_radius = sensitive_inner_radius;
    }
    else{ // support is on the inside
      support_radius = sensitive_radius - (0.5 *sensitive_thickness) -support_thickness;
      ladder_width = 2*(tan(ladder_dphi*0.5)*support_radius - ladder_clearance) ;
      inner_most_radius = support_radius;
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
    
    _e_r_p.sensor_length_mm  =sensor_length;
    
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
        
    Box sitLadderSolid( (sensitive_thickness +support_thickness ) / 2.0 ,
			layer_geom.ladder_width / 2.0,
			layer_geom.half_z);

    Volume sitLadderLogical (_toString( layer_id,"SIT_LadderLogical_%02d"), sitLadderSolid, air ) ; 
        
    // now create an envelope volume to represent the sensitive area, which will be divided up into individual sensors         
        
    Box sitSenEnvelopeSolid( (sensitive_thickness ) / 2.0 ,
			     layer_geom.ladder_width  / 2.0,
			     layer_geom.half_z);
    
    //fixme: material ???    Volume sitSenEnvelopeLogical( _toString( layer_id,"SIT_SenEnvelopeLogical_%02d"), sitSenEnvelopeSolid, sensitiveMat )  ;
    Volume sitSenEnvelopeLogical( _toString( layer_id,"SIT_SenEnvelopeLogical_%02d"), sitSenEnvelopeSolid, air )  ;
    
    // create the sensor volumes and place them in the senstive envelope volume 
    
    Box sitSenSolid( (sensitive_thickness ) / 2.0 ,
		     layer_geom.ladder_width  / 2.0,
		     (layer_geom.sensor_length / 2.0 ) - 1.e-06*mm ); // added tolerance to avoid false overlap detection
    
    Volume sitSenLogical(  _toString( layer_id,"SIT_SenLogical_%02d"), sitSenSolid,sensitiveMat ) ; 
    
    sitSenLogical.setSensitiveDetector(sens);


    //====== create the meassurement surface ===================
    Vector3D u,v,n ;

    if( faces_IP == 0 ){
      // will be rotated around z-axis later
      u.fill(  0. ,  -1. , 0. ) ;
      v.fill(  0. ,   0. , 1. ) ;
      n.fill( -1. ,   0. , 0. ) ;

      // implement 7 deg stereo angle 
      u.fill( 0. , -cos( 3.5/180.*M_PI  ) , -sin( 3.5/180.*M_PI  ) ) ;
      v.fill( 0. , -sin( 3.5/180.*M_PI  ) ,  cos( 3.5/180.*M_PI  ) ) ;

    } else {

      u.fill( 0. , 1. , 0. ) ;
      v.fill( 0. , 0. , 1. ) ;
      n.fill( 1. , 0. , 0. ) ;

      // implement 7 deg stereo angle 
      u.fill( 0. ,  cos( 3.5/180.*M_PI  ) ,  sin( 3.5/180.*M_PI  ) ) ;
      v.fill( 0. , -sin( 3.5/180.*M_PI  ) ,  cos( 3.5/180.*M_PI  ) ) ;
    }


    double inner_thick =  sensitive_thickness / 2.0 ;
    double outer_thick =  sensitive_thickness / 2.0 + support_thickness ;  // support is on top
   
    VolPlane surf( sitSenLogical , SurfaceType(SurfaceType::Sensitive) ,inner_thick, outer_thick , u,v,n ) ; //,o ) ;

    // vector of sensor placements - needed for DetElements in ladder loop below
    std::vector<PlacedVolume> pvV(  layer_geom.n_sensors_per_ladder ) ;

    //============================================================


    for (int isensor=0; isensor < layer_geom.n_sensors_per_ladder ; ++isensor) {
      
      // encoder.reset() ;  // reset to 0
      // encoder[ILDCellID0::subdet] = ILDDetID::NOTUSED ;
      // encoder[ILDCellID0::sensor] =  isensor+1;
      // cellID0 = encoder.lowWord() ;
      
      double xpos = 0.0;
      double ypos = 0.0;
      double zpos = -layer_geom.half_z + (0.5*layer_geom.sensor_length) + (isensor*layer_geom.sensor_length) ;
      
      pv = sitSenEnvelopeLogical.placeVolume( sitSenLogical, Transform3D( RotationY(0.) , Position( xpos, ypos, zpos)  ) );
      
      pv.addPhysVolID("sensor",  isensor + 1 ) ; 
      pvV[isensor] = pv ;
    }					      
    
    sit.setVisAttributes(lcdd, "SeeThrough",  sitLadderLogical ) ;
    sit.setVisAttributes(lcdd, "SeeThrough",  sitSenEnvelopeLogical ) ;

    sit.setVisAttributes(lcdd, "BlueVis",       sitSenLogical ) ;
    
    
    // encoder.reset() ;  // reset to 0
    // encoder[ILDCellID0::subdet] = ILDDetID::NOTUSED ;
    // encoder[ILDCellID0::layer]  = layer_id ;
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
    
    
    sit.setVisAttributes(lcdd, "RedVis",  sitSupLogical ) ;
    
    
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
      
      // encoder[ILDCellID0::subdet] = ILDDetID::SIT ;
      // encoder[ILDCellID0::layer]  = layer_id ;
      // encoder[ILDCellID0::module] = i + 1 ;
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
      
      pv.addPhysVolID("layer", layer_id ).addPhysVolID("module", i+1 ) ; 


      ladderDE.setPlacement( pv ) ;
    }
    
    
  }

  cout << "SIT_Simple_Planar done.\n" << endl;
  //######################################################################################################################################################################
  
  
  //--------------------------------------
  
  Volume mother =  lcdd.pickMotherVolume( sit ) ;
  
  pv = mother.placeVolume(envelope);

  pv.addPhysVolID( "system", x_det.id() ) ; //.addPhysVolID("side", 0 ) ;
  
  sit.setVisAttributes( lcdd, x_det.visStr(), envelope );
  //  if( sit.isValid() ) 
  sit.setPlacement(pv);
  
  return sit;
}
DECLARE_DETELEMENT(SIT_Simple_Planar,create_element);

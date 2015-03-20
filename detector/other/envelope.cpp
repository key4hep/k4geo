//====================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------
//  DD4hep Geometry generic envelope
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================
#include "envelope.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;


DD4hep::Geometry::Ref_t create_placed_envelope( DD4hep::Geometry::LCDD& lcdd, DD4hep::XML::Handle_t e , 
						DD4hep::Geometry::DetElement sdet ){
  
  xml_det_t     x_det     = e;
  
  int           det_id    = x_det.id();
  
  string        det_name  = x_det.nameStr();
  
  xml_comp_t    x_env     =  x_det.child( DD4hep::XML::Strng_t("envelope") ) ;
  xml_comp_t    x_shape   =  x_env.child( _U(shape) ); 
  
  Material      env_mat   = lcdd.material( x_shape.materialStr() );
  

  bool useRot = false ;
  bool usePos = false ; 
  Position    pos ;
  RotationZYX rot ;

  if( x_env.hasChild( _U(position) ) ) {
    usePos = true ;
    xml_comp_t env_pos = x_env.position();
    pos = Position( env_pos.x(),env_pos.y(),env_pos.z() );  
  }
  if( x_env.hasChild( _U(rotation) ) ) {
    useRot = true ;
    xml_comp_t    env_rot     = x_env.rotation();
    rot = RotationZYX( env_rot.z(),env_rot.y(),env_rot.x() ) ;
  }


  // ---- create a shape from the specified xml element --------
  Box  env_solid = xml_comp_t( x_shape ).createShape();
  
  if( !env_solid.isValid() ){
    
    throw std::runtime_error( std::string(" Cannot create envelope volume : ") + x_shape.typeStr() + 
			      std::string(" for detector " ) + det_name ) ;
  }
  
  Volume        envelope  ( det_name+"_envelope", env_solid, env_mat );
  
  PlacedVolume  env_pv  ; 

  Volume        mother = lcdd.pickMotherVolume(sdet);


  // ---- place the envelope into the mother volume 
  //      only specify transformations given in xml
  //      to allow for optimization 

  if( useRot && usePos ){
    env_pv =  mother.placeVolume( envelope , Transform3D( rot, pos )  ) ;

  } else if( useRot ){
    env_pv =  mother.placeVolume( envelope , rot  ) ;

  } else if( usePos ){
    env_pv =  mother.placeVolume( envelope , pos  ) ;

  } else {
    env_pv = mother.placeVolume( envelope );
  }

  // ----------------------------------------------

  env_pv.addPhysVolID("system", sdet.id());

  sdet.setPlacement( env_pv ) ;

  envelope.setAttributes( lcdd,x_det.regionStr(),x_det.limitsStr(),x_env.visStr());

  return envelope ;
}



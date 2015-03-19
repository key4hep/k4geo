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

/** Create an envelope volume that is placed into the world volume.
 *  to be called by subdetector drivers.
 *  @author S.Lu DESY, F. Gaede CERN/DESY 
 *  @version $Id:$
 */

DD4hep::Geometry::Ref_t create_envelope( DD4hep::Geometry::LCDD& lcdd, DD4hep::XML::Handle_t e , 
					 DD4hep::Geometry::DetElement sdet ){
  
  xml_det_t     x_det     = e;
  
  int           det_id    = x_det.id();
  
  string        det_name  = x_det.nameStr();

  xml_comp_t    x_shape   =  x_det.child( _U(shape) ); 

  Material      env_mat   = lcdd.material( x_shape.materialStr() );

  xml_comp_t    env_pos     = x_det.position();
  xml_comp_t    env_rot     = x_det.rotation();

  Position      pos( env_pos.x(),env_pos.y(),env_pos.z() );
  RotationZYX   rotZYX( env_rot.z(),env_rot.y(),env_rot.x() );

  Transform3D   tr( rotZYX, pos );

  //  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  Box  env_solid = xml_comp_t(x_shape).createShape();

  if( !env_solid.isValid() ){

    throw std::runtime_error( std::string(" Cannot create envelope volume : ") + x_shape.typeStr() + 
			      std::string(" for detector " ) + det_name ) ;
  }

  Volume        envelope  ( det_name+"_envelope", env_solid, env_mat );

  PlacedVolume  env_phv   =  motherVol.placeVolume(envelope,tr);

  env_phv.addPhysVolID("system", sdet.id());

  sdet.setPlacement( env_phv ) ;

  envelope.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

  return envelope ;

  //  return sdet;
}



//====================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------
//  DD4hep Geometry tube envelope
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
  static double tolerance = 0e0;
  Layering      layering (e);
  xml_det_t     x_det     = e;
  Material      air       = lcdd.air();

  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();

  xml_comp_t    x_dim     = x_det.dimensions();

  double        thickness = layering.totalThickness();

  double        inner_r   = x_dim.rmin();
  double        outer_r   = x_dim.rmax();

  //assert(thickness < (outer_r - inner_r));
  
  Material      env_mat   = lcdd.material(x_dim.materialStr());

  xml_comp_t    env_pos     = x_det.position();
  xml_comp_t    env_rot     = x_det.rotation();

  Position      pos(env_pos.x(),env_pos.y(),env_pos.z());
  RotationZYX   rotZYX(env_rot.z(),env_rot.y(),env_rot.x());

  Transform3D   tr(rotZYX,pos);

  xml_comp_t    x_staves  = x_det.staves();


  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  // The shape Tube for envelope. You may define your shape.
  Tube          tubeSolid (inner_r, outer_r, x_dim.zhalf());
  Volume        envelope  (det_name+"_envelope",tubeSolid,env_mat);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,tr);


  // Then build and place the layers into the envelope in your way here.


  envelope.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  return sdet;
}

DECLARE_DETELEMENT(tube_envelope,create_detector)

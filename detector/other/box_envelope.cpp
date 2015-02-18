//====================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------
//  DD4hep Geometry box envelope
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector /*sens*/)  {
  //unused:  static double tolerance = 0e0;
  Layering      layering (e);
  xml_det_t     x_det     = e;
  Material      air       = lcdd.air();

  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();

  xml_comp_t    x_dim     = x_det.dimensions();

  //unused:  double        thickness = layering.totalThickness();

  double        dim_x   = x_dim.X();
  double        dim_y   = x_dim.Y();
  double        dim_z   = x_dim.Z();

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

  // The shape Box for envelope. You may define your shape.
  Box           boxSolid  (dim_x, dim_y, dim_z);
  Volume        envelope  (det_name+"_envelope",boxSolid,env_mat);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,tr);


  // Then build and place the layers into the envelope in your way here.


  envelope.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  return sdet;
}

DECLARE_DETELEMENT(box_envelope,create_detector)

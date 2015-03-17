//====================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------
//  DD4hep Geometry generic envelope
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector /*sens*/)  {

  xml_det_t     x_det     = e;

  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();

  xml_comp_t    x_dim     = x_det.dimensions();

  Material      env_mat   = lcdd.material(x_dim.materialStr());

  xml_comp_t    env_pos     = x_det.position();
  xml_comp_t    env_rot     = x_det.rotation();

  Position      pos(env_pos.x(),env_pos.y(),env_pos.z());
  RotationZYX   rotZYX(env_rot.z(),env_rot.y(),env_rot.x());

  Transform3D   tr(rotZYX,pos);

  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  // https://root.cern.ch/root/html/TGeoShape.html
  // TNamed <- TGeoShape <- TGeoBBox <- TGeoCompositeShape
  // Use "Box" for all type shapes of solid
  Box  env_solid = xml_comp_t(x_dim).createShape();
  if( !env_solid.isValid()) cout <<" createShape() has problem! "<<endl;

  Volume        envelope  (det_name+"_envelope",env_solid,env_mat);

  PlacedVolume  env_phv   =  motherVol.placeVolume(envelope,tr);
  sdet.setPlacement( env_phv ) ;

  envelope.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  return sdet;
}

DECLARE_DETELEMENT(Envelope,create_detector)

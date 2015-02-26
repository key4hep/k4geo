//====================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------
//  DD4hep Geometry polyhedra outer and circle hole inner envelope
//                  for symmetry endcap
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
  static double tolerance = 0e0;
  Layering      layering (e);
  xml_det_t     x_det     = e;
 //unused:  Material      air       = lcdd.air();

  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();

  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();

  //unused: double        dphi      = (2*M_PI/nsides);
  //unused: double        hphi      = dphi/2;

  //unused: double        thickness = layering.totalThickness();

  double        inner_r   = x_dim.rmin();
  double        outer_r   = x_dim.rmax();
  double        endcap_zhalf   = x_dim.zhalf();
  double        endcap_zstart   = abs(x_dim.zstart());

  //assert(thickness < (outer_r - inner_r));
  
  Material      env_mat   = lcdd.material(x_dim.materialStr());

  xml_comp_t    env_pos     = x_det.position();
  xml_comp_t    env_rot     = x_det.rotation();

  Position      pos(env_pos.x(),env_pos.y(),env_pos.z());
  RotationZYX   rotZYX(env_rot.z(),env_rot.y(),env_rot.x());

  Transform3D   tr(rotZYX,pos);

  //unused:  xml_comp_t    x_staves  = x_det.staves();


  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  // The shape PolyhedraRegular for envelope. You may define your shape.
  PolyhedraRegular hedra_solid (nsides,0,outer_r+tolerance*2e0, (endcap_zstart + endcap_zhalf*2.0)*2.0 );
  //Tube             tube_hole(0, inner_r+tolerance*2e0,endcap_zstart+endcap_zhalf*2.0+2.0e-5);
  Cone             tube_hole(endcap_zstart+endcap_zhalf*2.0+2.0e-5, 0, inner_r+tolerance*2e0, 0, inner_r+tolerance*2e0);
  SubtractionSolid hedratube(hedra_solid,tube_hole);

  if( endcap_zstart > 0.0 ){
    PolyhedraRegular hedra_cut (nsides,0,outer_r+tolerance*2e0+1.0e-5, endcap_zstart*2.0 );
    SubtractionSolid hedratube_env(hedratube,hedra_cut);

    Volume        envelope  (det_name+"_envelope",hedratube_env,env_mat);
    envelope.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

    // Then build and place the layers into the envelope in your way here.

    //unused:  PlacedVolume  env_phv   = 
    motherVol.placeVolume(envelope,tr);
  }else {
    cout<<" WARNING: hedra_tube_sys_envelope: Nothing has been built!"<<endl;
    cout<<" WARNING: Please check \"zstart\" value larger than 0.0*mm"<<endl;
  }


  return sdet;
}

DECLARE_DETELEMENT(hedra_tube_symmetry_envelope,create_detector)

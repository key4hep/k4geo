//====================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------
//  DD4hep Geometry polyhedra envelope
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t createPolyhedraEnvelope(LCDD& lcdd, xml_h e ) {

  static double tolerance = 0e0;

  //xml_det_t     x_det     = e;
  xml_comp_t x_det     = e;

  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();

  double        inner_r   = x_dim.rmin();
  double        outer_r   = x_dim.rmax();

  PolyhedraRegular hedra  ( nsides, 
			    inner_r, 
			    outer_r+tolerance*2e0, 
			    x_dim.zhalf()*2.0
			    );

  Material   mat   = lcdd.material( x_det.materialStr() );
  string     det_name = x_det.nameStr();

  return hedra ;
}

DECLARE_XMLELEMENT(PolyhedraEnvelope__shape_constructor,createPolyhedraEnvelope)

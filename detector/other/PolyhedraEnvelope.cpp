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

  xml_det_t     x_det     = e;

  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();

  double        inner_r   = x_dim.rmin();
  double        outer_r   = x_dim.rmax();

  // The shape PolyhedraRegular for envelope. You may define your shape.
  PolyhedraRegular hedra  ( nsides, 
			    inner_r, 
			    outer_r+tolerance*2e0, 
			    x_dim.zhalf()*2.0
			    );

  return hedra ;
}


//DECLARE_XMLELEMENT(PolyhedraEnvelope,createPolyhedraEnvelope)
DECLARE_XMLELEMENT(PolyhedraEnvelope,createPolyhedraEnvelope)

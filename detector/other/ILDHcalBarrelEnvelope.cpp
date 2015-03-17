//====================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------
//  DD4hep Geometry user defined envelope shape
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t createILDHcalBarrelEnvelope(LCDD&, xml_h element ) {

  xml_dim_t e(element);

  // a composite shape
  Cone             tube_solid(e.zhalf(), 0, e.rmax(), 0, e.rmax());
  PolyhedraRegular hedra_hole (e.numsides(),0,e.rmin(),e.zhalf()*2.0+1.0e-5);
  SubtractionSolid tubehedra(tube_solid,hedra_hole);
  return tubehedra;

  // Or a simple shape
  //return PolyhedraRegular(e.numsides(),e.rmin(),e.rmax(),e.zhalf()*2.0);
}

DECLARE_XMLELEMENT(ILDHcalBarrelEnvelope__shape_constructor, createILDHcalBarrelEnvelope)

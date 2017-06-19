#include "OtherDetectorHelpers.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"
#include "XMLHandlerDB.h"
#include <cmath>
#include <map>
#include <string>

using namespace std;

using dd4hep::Assembly;
using dd4hep::ConeSegment;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotateY;
using dd4hep::RotationY;
using dd4hep::Solid;
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Tube;
using dd4hep::Volume;
using dd4hep::Detector;
using dd4hep::SensitiveDetector;
using dd4hep::DetElement;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::_toString;

using dd4hep::rec::Vector3D;
using dd4hep::rec::VolSurfaceHandle;
using dd4hep::rec::VolPlaneImpl;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::volSurfaceList;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector)
{
  xml_det_t     x_det     = e;
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  DetElement    sdet(det_name, det_id);

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector,  e , sdet) ;

  if (theDetector.buildType() == BUILD_ENVELOPE) return sdet ;

  Material   air       = theDetector.air();
  bool       reflect   = x_det.reflect();

  PlacedVolume pv;
  int l_num = 0;

  for (xml_coll_t i(x_det, _U(layer)); i; ++i, ++l_num)  {
    xml_comp_t x_layer = i;
    string l_nam = det_name + _toString(l_num, "_layer%d");
    double  zmin = x_layer.inner_z();
    double  rmin = x_layer.inner_r();
    double  rmax = x_layer.outer_r();
    double  z    = zmin, layerWidth = 0.;
    int     s_num = 0;

    for (xml_coll_t j(x_layer, _U(slice)); j; ++j)  {
      double thickness = xml_comp_t(j).thickness();
      layerWidth += thickness;
    }
    Tube    l_tub(rmin, rmax, layerWidth, 2 * M_PI);
    Volume  l_vol(l_nam, l_tub, air);
    l_vol.setVisAttributes(theDetector, x_layer.visStr());
    for (xml_coll_t j(x_layer, _U(slice)); j; ++j, ++s_num)  {
      xml_comp_t x_slice = j;
      double thick = x_slice.thickness();
      Material mat = theDetector.material(x_slice.materialStr());
      string s_nam = l_nam + _toString(s_num, "_slice%d");
      Volume s_vol(s_nam, Tube(rmin, rmax, thick/2.), mat); //NN: Tube is a MyConeSeg, takes half thickness

      //Add surface to the support
      double mid_r = 0.5 * ( rmin + rmax ) ;
      //      Vector3D ocyl(  0., mid_r , z - zmin - layerWidth / 2 + thick / 2 );
      Vector3D ocyl(  0., mid_r ,  0. );
      Vector3D u(1.,0.,0.), v(0.,1.,0.), n(0.,0.,1.);      

	  VolSurfaceHandle<VolPlaneImpl> cylSurf1( s_vol , SurfaceType( SurfaceType::Helper ) , 0.5*thick  , 0.5*thick , u, v, n, ocyl );

	  volSurfaceList( sdet )->push_back( cylSurf1 );
	   

      s_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
      pv = l_vol.placeVolume(s_vol, Position(0, 0, z - zmin - layerWidth / 2 + thick / 2));
      pv.addPhysVolID("sensor", s_num);
    }

    DetElement layer_pos(sdet, l_nam + "_pos", l_num);
    pv = envelope.placeVolume(l_vol, Position(0, 0, zmin + layerWidth / 2.));
    pv.addPhysVolID("layer", l_num);
    pv.addPhysVolID("barrel", 1);
    layer_pos.setPlacement(pv);
    
    if (reflect)  {
	  Tube    l_tub2(rmin, rmax, layerWidth, 2 * M_PI);
      Volume  l_vol2(l_nam, l_tub2, air);
      l_vol2.setVisAttributes(theDetector, x_layer.visStr());
      for (xml_coll_t j(x_layer, _U(slice)); j; ++j, ++s_num)  {
        xml_comp_t x_slice = j;
        double thick = x_slice.thickness();
        Material mat = theDetector.material(x_slice.materialStr());
        string s_nam = l_nam + _toString(s_num, "_slice%d");
        Volume s_vol2(s_nam, Tube(rmin, rmax, thick), mat);

        //Add surface to the support
	double mid_r = 0.5 * ( rmin + rmax ) ;
	//        Vector3D ocyl(  0., mid_r , z - zmin - layerWidth / 2 + thick / 2 );
        Vector3D ocyl(  0., mid_r , 0. );
        Vector3D u(1.,0.,0.), v(0.,1.,0.), n(0.,0.,1.);      

	    VolSurfaceHandle<VolPlaneImpl> cylSurf2( s_vol2 , SurfaceType( SurfaceType::Helper ) , 0.5*thick  , 0.5*thick , u, v, -1.*n, -1.*ocyl );

	    volSurfaceList( sdet )->push_back( cylSurf2 );
	   

        s_vol2.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
        pv = l_vol2.placeVolume(s_vol2, Position(0, 0, -1.*(z - zmin - layerWidth / 2 + thick / 2)));
        pv.addPhysVolID("sensor", s_num);
      }

    DetElement layer_neg(sdet, l_nam + "_neg", l_num);
    pv = envelope.placeVolume(l_vol2, Position(0, 0, -1.*(zmin + layerWidth / 2.)));
    pv.addPhysVolID("layer", l_num);
    pv.addPhysVolID("barrel", 2);
    layer_neg.setPlacement(pv);
    }
  }


  return sdet;
}

DECLARE_DETELEMENT(TrackerEndcapSupport_o1_v01, create_detector)

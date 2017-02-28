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
using namespace DD4hep;
using namespace DD4hep::Geometry;
using namespace DD4hep::DDRec ;
using namespace DDSurfaces ;

using DD4hep::Geometry::Transform3D;
using DD4hep::Geometry::Position;
using DD4hep::Geometry::RotationY;
using DD4hep::Geometry::RotateY;
using DD4hep::Geometry::ConeSegment;
using DD4hep::Geometry::SubtractionSolid;
using DD4hep::Geometry::Material;
using DD4hep::Geometry::Volume;
using DD4hep::Geometry::Solid;
using DD4hep::Geometry::Tube;
using DD4hep::Geometry::PlacedVolume;
using DD4hep::Geometry::Assembly;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector /*sens*/)
{

  xml_det_t     x_det     = e;
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  DetElement    sdet(det_name, det_id);
  bool          reflect   = x_det.reflect();

  Volume envelope = XML::createPlacedEnvelope(lcdd,  e , sdet) ;

  if (lcdd.buildType() == BUILD_ENVELOPE) return sdet ;

  for (xml_coll_t c(x_det, Unicode("section")); c; ++c) {

    xml_comp_t xmlSection(c);

    const double zStart       = xmlSection.attr< double > (_Unicode(start));
    const double zEnd         = xmlSection.attr< double > (_Unicode(end));
    const double rInner       = xmlSection.attr< double > (_Unicode(rMin));
    const double rOuter       = xmlSection.attr< double > (_Unicode(rMax));
    const double thickness    = rOuter - rInner;
    Material sectionMat       = lcdd.material(xmlSection.materialStr());
    const std::string volName      = "section_" + xmlSection.nameStr();

    const double zHalf       = fabs(zEnd - zStart) * 0.5; // half z length of the cone
    const double zPosition   = fabs(zEnd + zStart) * 0.5; // middle z position

    // solid for the tube (including vacuum and wall): a solid cone
    Tube tubeSolid(rInner, rOuter, zHalf);

    Volume tubeVol(volName + "_pos", tubeSolid, sectionMat);
    tubeVol.setAttributes(lcdd, xmlSection.regionStr(), xmlSection.limitsStr(), xmlSection.visStr());
    envelope.placeVolume(tubeVol, Position(0, 0, zPosition));

    //Add surface to the support
    Vector3D ocyl(  rInner + thickness/2.  , 0. , 0. );
    VolCylinder cylSurf1( tubeVol , SurfaceType( SurfaceType::Helper ) , 0.5*thickness  , 0.5*thickness , ocyl );
    volSurfaceList( sdet )->push_back( cylSurf1 );


    if (reflect) {

      Volume tubeVol2(volName + "_neg", tubeSolid, sectionMat);
      tubeVol2.setAttributes(lcdd, xmlSection.regionStr(), xmlSection.limitsStr(), xmlSection.visStr());
      Transform3D Vol2Place(RotationY(-180.0 * dd4hep::degree), Position(0, 0, -1.*zPosition));
      envelope.placeVolume(tubeVol2, Vol2Place);

      VolCylinder cylSurf2( tubeVol2 , SurfaceType( SurfaceType::Helper ) , 0.5*thickness  , 0.5*thickness , ocyl );
      volSurfaceList( sdet )->push_back( cylSurf2 );


    }


  }

  return sdet;

}

DECLARE_DETELEMENT(TubeSupport_o1_v01, create_detector)

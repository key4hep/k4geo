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
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotateY;
using dd4hep::RotationY;
using dd4hep::SensitiveDetector;
using dd4hep::Solid;
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Tube;
using dd4hep::Volume;
using dd4hep::BUILD_ENVELOPE;

using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::volSurfaceList;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector /*sens*/)
{

  xml_det_t     x_det     = e;
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  DetElement    sdet(det_name, det_id);
  bool          reflect   = x_det.reflect();

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector,  e , sdet) ;

  if (theDetector.buildType() == BUILD_ENVELOPE) return sdet ;

  for (xml_coll_t c(x_det, Unicode("section")); c; ++c) {

    xml_comp_t xmlSection(c);

    const double zStart       = xmlSection.attr< double > (_Unicode(start));
    const double zEnd         = xmlSection.attr< double > (_Unicode(end));
    const double rInner       = xmlSection.attr< double > (_Unicode(rMin));
    const double rOuter       = xmlSection.attr< double > (_Unicode(rMax));
    const double thickness    = rOuter - rInner;
    Material sectionMat       = theDetector.material(xmlSection.materialStr());
    const std::string volName      = "section_" + xmlSection.nameStr();

    const double zHalf       = fabs(zEnd - zStart) * 0.5; // half z length of the cone
    const double zPosition   = fabs(zEnd + zStart) * 0.5; // middle z position

    // solid for the tube (including vacuum and wall): a solid cone
    Tube tubeSolid(rInner, rOuter, zHalf);

    Volume tubeVol(volName + "_pos", tubeSolid, sectionMat);
    tubeVol.setAttributes(theDetector, xmlSection.regionStr(), xmlSection.limitsStr(), xmlSection.visStr());
    envelope.placeVolume(tubeVol, Position(0, 0, zPosition));

    //Add surface to the support
    Vector3D ocyl(  rInner + thickness/2.  , 0. , 0. );
    VolCylinder cylSurf1( tubeVol , SurfaceType( SurfaceType::Helper ) , 0.5*thickness  , 0.5*thickness , ocyl );
    volSurfaceList( sdet )->push_back( cylSurf1 );


    if (reflect) {

      Volume tubeVol2(volName + "_neg", tubeSolid, sectionMat);
      tubeVol2.setAttributes(theDetector, xmlSection.regionStr(), xmlSection.limitsStr(), xmlSection.visStr());
      Transform3D Vol2Place(RotationY(-180.0 * dd4hep::degree), Position(0, 0, -1.*zPosition));
      envelope.placeVolume(tubeVol2, Vol2Place);

      VolCylinder cylSurf2( tubeVol2 , SurfaceType( SurfaceType::Helper ) , 0.5*thickness  , 0.5*thickness , ocyl );
      volSurfaceList( sdet )->push_back( cylSurf2 );


    }


  }

  return sdet;

}

DECLARE_DETELEMENT(TubeSupport_o1_v01, create_detector)

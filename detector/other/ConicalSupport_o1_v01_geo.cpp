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
using dd4hep::BUILD_ENVELOPE;
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

using dd4hep::rec::ConicalSupportData;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCone;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::VolCylinderImpl;
using dd4hep::rec::VolSurface;
using dd4hep::rec::volSurfaceList;

static Ref_t create_detector(Detector& lcdd, xml_h e, SensitiveDetector /*sens*/)
{

  xml_det_t     x_det     = e;
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  DetElement    sdet(det_name, det_id);
  bool          reflect   = x_det.reflect();

  Volume envelope = dd4hep::xml::createPlacedEnvelope(lcdd,  e , sdet) ;

  if (lcdd.buildType() == BUILD_ENVELOPE) return sdet ;

  const double phi1 = 0 ;
  const double phi2 = 360.0 * dd4hep::degree;


  for (xml_coll_t c(x_det, Unicode("section")); c; ++c) {

    xml_comp_t xmlSection(c);

    const double zStart       = xmlSection.attr< double > (_Unicode(start));
    const double zEnd         = xmlSection.attr< double > (_Unicode(end));
    const double rInnerStart  = xmlSection.attr< double > (_Unicode(rMin1));
    const double rInnerEnd    = xmlSection.attr< double > (_Unicode(rMin2));
    const double rOuterStart  = xmlSection.attr< double > (_Unicode(rMax1));
    const double rOuterEnd    = xmlSection.attr< double > (_Unicode(rMax2));
    const double thickness    = rOuterStart - rInnerStart;
    Material sectionMat       = lcdd.material(xmlSection.materialStr());
    const std::string volName      = "section_" + xmlSection.nameStr();

    const double zHalf       = fabs(zEnd - zStart) * 0.5; // half z length of the cone
    const double zPosition   = fabs(zEnd + zStart) * 0.5; // middle z position

    // solid for the tube (including vacuum and wall): a solid cone
    ConeSegment tubeSolid(zHalf, rInnerStart, rOuterStart, rInnerEnd, rOuterEnd , phi1, phi2);

    Volume tubeVol(volName + "_pos", tubeSolid, sectionMat);
    tubeVol.setAttributes(lcdd, xmlSection.regionStr(), xmlSection.limitsStr(), xmlSection.visStr());
    envelope.placeVolume(tubeVol, Position(0, 0, zPosition));

    const double dr    = rInnerEnd - rInnerStart ;
    const double theta = atan2(dr , 2.* zHalf) ;

    Vector3D ocon(rInnerStart + 0.5 * (dr + thickness), 0. , 0.);
    Vector3D v(1. , 0. , theta, Vector3D::spherical) ;
    VolCone conSurf1(tubeVol , SurfaceType(SurfaceType::Helper) , 0.5 * thickness  , 0.5 * thickness , v, ocon);
    volSurfaceList(sdet)->push_back(conSurf1);


    if (reflect) {

      Volume tubeVol2(volName + "_neg", tubeSolid, sectionMat);
      tubeVol2.setAttributes(lcdd, xmlSection.regionStr(), xmlSection.limitsStr(), xmlSection.visStr());
      Transform3D Vol2Place(RotationY(-180.0 * dd4hep::degree), Position(0, 0, -1.*zPosition));
      envelope.placeVolume(tubeVol2, Vol2Place);

      VolCone conSurf2(tubeVol2, SurfaceType(SurfaceType::Helper) , 0.5 * thickness  , 0.5 * thickness , v, ocon);
      volSurfaceList(sdet)->push_back(conSurf2);


    }


  }

  return sdet;

}

DECLARE_DETELEMENT(ConicalSupport_o1_v01, create_detector)

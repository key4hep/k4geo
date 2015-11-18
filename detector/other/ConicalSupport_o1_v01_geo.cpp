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

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)
{
  static double tolerance = 0e0;

  xml_det_t     x_det     = e;
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  DetElement    sdet(det_name, det_id);
  bool          reflect   = x_det.reflect();
  
  Volume envelope = XML::createPlacedEnvelope(lcdd,  e , sdet) ;

  if (lcdd.buildType() == BUILD_ENVELOPE) return sdet ;

  const double phi1 = 0 ;
  const double phi2 = 360.0*dd4hep::degree;


  for(xml_coll_t c(x_det, Unicode("section")); c; ++c) {

    xml_comp_t xmlSection( c );
    
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
    ConeSegment tubeSolid1( zHalf, rInnerStart, rOuterStart, rInnerEnd, rOuterEnd , phi1, phi2);
    ConeSegment tubeSolid2( zHalf*1.01, 0. , rInnerStart, 0. , rInnerEnd , phi1, phi2);
	Solid tubeSolid = SubtractionSolid(tubeSolid1,tubeSolid2);
	
	Volume tubeVol( volName + "_pos", tubeSolid, sectionMat );
    tubeVol.setAttributes(lcdd, xmlSection.regionStr(), xmlSection.limitsStr(), xmlSection.visStr());
	envelope.placeVolume( tubeVol, Position(0, 0, zPosition) );
	
	const double dr    = rInnerEnd - rInnerStart ;
	const double theta = atan2( dr , 2.* zHalf ) ;

	Vector3D ocon( rInnerStart + 0.5 * ( dr + thickness ), 0. , 0. );
	Vector3D v( 1. , 0. , theta, Vector3D::spherical ) ;
	VolCone conSurf1( tubeVol , SurfaceType( SurfaceType::Helper ) , 0.5*thickness  , 0.5*thickness , v, ocon );
	volSurfaceList( sdet )->push_back( conSurf1 );

	
	if(reflect){
		
	  Volume tubeVol2( volName + "_neg", tubeSolid, sectionMat );
      tubeVol2.setAttributes(lcdd, xmlSection.regionStr(), xmlSection.limitsStr(), xmlSection.visStr());
	  Transform3D Vol2Place(RotationY(-180.0*dd4hep::degree), Position(0, 0, -1.*zPosition));
	  envelope.placeVolume( tubeVol2, Vol2Place );
	  
	  VolCone conSurf2( tubeVol2, SurfaceType( SurfaceType::Helper ) , 0.5*thickness  , 0.5*thickness , v, ocon );
	  volSurfaceList( sdet )->push_back( conSurf2 );
	
	  
	}
	
	
  }

  return sdet;

}

DECLARE_DETELEMENT(ConicalSupport_o1_v01, create_detector)

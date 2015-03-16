//====================================================================
//  LCGeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
// Beampipe without the use of global constants, just reads one xml Section after the other
//  A.Sailer, CERN
//  $Id$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "XMLHandlerDB.h"
#include <cmath>
#include <map>
#include <string>

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

namespace {
  typedef enum { // These constants are also used in the MySQL database:
    kCenter                     = 0, // centered on the z-axis
    kUpstream                   = 1, // on the upstream branch, rotated by half the crossing angle
    kDnstream                   = 2, // on the downstream branch, rotated by half the crossing angle
    kPunchedCenter              = 3, // centered, with one or two inner holes
    kPunchedUpstream            = 4, // on the upstream branch, with two inner holes
    kPunchedDnstream            = 5, // on the downstrem branch, with two inner holes
    kUpstreamClippedFront       = 6, // upstream, with the front face parallel to the xy-plane
    kDnstreamClippedFront       = 7, // downstream, with the front face parallel to the xy-plane
    kUpstreamClippedRear        = 8, // upstream, with the rear face parallel to the xy-plane
    kDnstreamClippedRear        = 9, // downstream, with the rear face parallel to the xy-plane
    kUpstreamClippedBoth        = 10, // upstream, with both faces parallel to the xy-plane
    kDnstreamClippedBoth        = 11, // downstream, with both faces parallel to the xy-plane
    kUpstreamSlicedFront        = 12, // upstream, with the front face parallel to a tilted piece
    kDnstreamSlicedFront        = 13, // downstream, with the front face parallel to a tilted piece
    kUpstreamSlicedRear         = 14, // upstream, with the rear face parallel to a tilted piece
    kDnstreamSlicedRear         = 15, // downstream, with the rear face parallel to a tilted piece
    kUpstreamSlicedBoth         = 16, // upstream, with both faces parallel to a tilted piece
    kDnstreamSlicedBoth         = 17 // downstream, with both faces parallel to a tilted piece
  } ECrossType;
}//namespace



bool checkForSensibleGeometry(double crossingAngle, ECrossType crossType) {
  if (crossingAngle == 0 && crossType != kCenter) {
    std::cout << "TubeX01: You are trying to build a crossing geometry without a crossing angle.\n"
      "This is probably not what you want - better check your geometry data!" << std::endl ;
    return false; // premature exit, dd4hep will abort now
  }
  return true;
}


double getCurrentAngle( double crossingAngle, ECrossType crossType ) {
  double tmpAngle;
  switch (crossType) {
  case kUpstream:
  case kPunchedUpstream:
  case kUpstreamClippedFront:
  case kUpstreamClippedRear:
  case kUpstreamClippedBoth:
  case kUpstreamSlicedFront:
  case kUpstreamSlicedRear:
  case kUpstreamSlicedBoth:
    tmpAngle = -crossingAngle; break;
  case kDnstream:
  case kPunchedDnstream:
  case kDnstreamClippedFront:
  case kDnstreamClippedRear:
  case kDnstreamClippedBoth:
  case kDnstreamSlicedFront:
  case kDnstreamSlicedRear:
  case kDnstreamSlicedBoth:
    tmpAngle = +crossingAngle; break;
  default:
    tmpAngle = 0; break;
  }

  return tmpAngle;
}

/** Construction of VTX detector, ported from Mokka driver TubeX01.cc
 *
 *  Mokka History:
 * - first implementation as Tube00: Paulo Mora de Freitas, Sep 2002
 * - modified from Tube00 to Tube01DT: Ties Behnke, 2003-02-11
 * - modified for a crossing angle as TubeX00: Adrian Vogel, 2005-05-18
 * - modified for fancier geometries as TubeX01: Adrian Vogel, 2006-04-20
 *
 *  @author: F.Gaede, DESY, Jan 2014
 *
 */
static DD4hep::Geometry::Ref_t create_element(DD4hep::Geometry::LCDD& lcdd,
					      DD4hep::XML::Handle_t xmlHandle,
					      DD4hep::Geometry::SensitiveDetector /*sens*/) {


  std::map< std::string, ECrossType > CrossTypes;
  CrossTypes["Center"]                = kCenter               ;
  CrossTypes["Upstream"]              = kUpstream             ;
  CrossTypes["Dnstream"]              = kDnstream             ;
  CrossTypes["PunchedCenter"]         = kPunchedCenter        ;
  CrossTypes["PunchedUpstream"]       = kPunchedUpstream      ;
  CrossTypes["PunchedDnstream"]       = kPunchedDnstream      ;
  CrossTypes["UpstreamClippedFront"]  = kUpstreamClippedFront ;
  CrossTypes["DnstreamClippedFront"]  = kDnstreamClippedFront ;
  CrossTypes["UpstreamClippedRear"]   = kUpstreamClippedRear  ;
  CrossTypes["DnstreamClippedRear"]   = kDnstreamClippedRear  ;
  CrossTypes["UpstreamClippedBoth"]   = kUpstreamClippedBoth  ;
  CrossTypes["DnstreamClippedBoth"]   = kDnstreamClippedBoth  ;
  CrossTypes["UpstreamSlicedFront"]   = kUpstreamSlicedFront  ;
  CrossTypes["DnstreamSlicedFront"]   = kDnstreamSlicedFront  ;
  CrossTypes["UpstreamSlicedRear"]    = kUpstreamSlicedRear   ;
  CrossTypes["DnstreamSlicedRear"]    = kDnstreamSlicedRear   ;
  CrossTypes["UpstreamSlicedBoth"]    = kUpstreamSlicedBoth   ;
  CrossTypes["DnstreamSlicedBoth"]    = kDnstreamSlicedBoth   ;


  //------------------------------------------
  //  See comments starting with '//**' for
  //     hints on porting issues
  //------------------------------------------

  std::cout << "This is the Beampipe:"  << std::endl;

  //Access to the XML File
  DD4hep::XML::DetElement xmlBeampipe = xmlHandle;
  const std::string name = xmlBeampipe.nameStr();

  //--------------------------------
  Assembly envelope( name + "_assembly"  ) ;
  //--------------------------------

  DD4hep::Geometry::DetElement tube(  name, xmlBeampipe.id()  ) ;

  DD4hep::DDRec::ConicalSupportData* beampipeData = new DD4hep::DDRec::ConicalSupportData ;

  //######################################################################################################################################################################
  //  code ported from TubeX01::construct() :
  //##################################
  
  //** DD4hep/TGeo seems to need rad (as opposed to the manual)
  const double phi1 = 0 ;
  const double phi2 = 360.0*dd4hep::degree;
  
  //Parameters we have to know about
  DD4hep::XML::Component xmlParameter = xmlBeampipe.child(_Unicode(parameter));
  const double crossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle))*0.5; //  only half the angle

  for(xml_coll_t c( xmlBeampipe ,Unicode("section")); c; ++c) {

    xml_comp_t xmlSection( c );
    
    ECrossType crossType = CrossTypes[xmlSection.attr< std::string >(_Unicode(type))];
    const double zStart       = xmlSection.attr< double > (_Unicode(start));
    const double zEnd         = xmlSection.attr< double > (_Unicode(end));
    const double rInnerStart  = xmlSection.attr< double > (_Unicode(rMin1));
    const double rInnerEnd    = xmlSection.attr< double > (_Unicode(rMin2));
    const double rOuterStart  = xmlSection.attr< double > (_Unicode(rMax1));
    const double rOuterEnd    = xmlSection.attr< double > (_Unicode(rMax2));
    const double thickness    = rOuterStart - rInnerStart;
    Material sectionMat  = lcdd.material(xmlSection.materialStr());
    const std::string volName      = "tube_" + xmlSection.nameStr();

    std::cout << std::setw(8) << zStart
	      << std::setw(8) << zEnd
	      << std::setw(8) << rInnerStart
	      << std::setw(8) << rInnerEnd
	      << std::setw(8) << rOuterStart
	      << std::setw(8) << rOuterEnd
	      << std::setw(8) << thickness
	      << std::setw(8) << crossType
	      << std::setw(15) << volName
	      << std::setw(15) << sectionMat.name()  << std::endl;    

    if( crossType == kCenter ) { // store only the central sections !
      DD4hep::DDRec::ConicalSupportData::Section section ;
      section.rInner = rInnerStart ;
      section.rOuter = rOuterStart ;
      section.zPos   = zStart ;
      beampipeData->sections.push_back( section ) ;
    }

    // things which can be calculated immediately
    const double zHalf       = fabs(zEnd - zStart) * 0.5; // half z length of the cone
    const double zPosition   = fabs(zEnd + zStart) * 0.5; // middle z position
    Material coreMaterial    = lcdd.material("beam"); // always the same
    Material wallMaterial    = sectionMat;

    // this could mess up your geometry, so better check it
    if (not checkForSensibleGeometry(crossingAngle, crossType)){

      throw std::runtime_error( " Beampipe_o1_v01_geo.cpp : checkForSensibleGeometry() failed " ) ;
      //      return false;
    }
    const double rotateAngle = getCurrentAngle(crossingAngle, crossType); // for the placement at +z (better make it const now)
    const double mirrorAngle = M_PI - rotateAngle; // for the "mirrored" placement at -z
    // the "mirroring" in fact is done by a rotation of (almost) 180 degrees around the y-axis

    switch (crossType) {
    case kCenter:
    case kUpstream:
    case kDnstream: {
      // a volume on the z-axis, on the upstream branch, or on the downstream branch
      
      // absolute transformations for the placement in the world
      Transform3D transformer(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition), rotateAngle) );
      Transform3D transmirror(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition), mirrorAngle) );
      
      // solid for the tube (including vacuum and wall): a solid cone
      ConeSegment tubeSolid( zHalf, 0, rOuterStart, 0, rOuterEnd , phi1, phi2);
        
      // tube consists of vacuum
      Volume tubeLog( volName, tubeSolid, coreMaterial ) ;
      
      // placement of the tube in the world, both at +z and -z
      envelope.placeVolume( tubeLog,  transformer );
      envelope.placeVolume( tubeLog,  transmirror );

      // if inner and outer radii are equal, then omit the tube wall
      if (rInnerStart != rOuterStart || rInnerEnd != rOuterEnd) {
	
	// the wall solid: a tubular cone
	ConeSegment wallSolid( zHalf, rInnerStart, rOuterStart, rInnerEnd, rOuterEnd, phi1, phi2);
        
	// the wall consists of the material given in the XML
	Volume wallLog ( volName + "_wall", wallSolid, wallMaterial);
	
	wallLog.setVisAttributes(lcdd, "TubeVis");
	tubeLog.setVisAttributes(lcdd, "VacVis");
	
	// placement as a daughter volume of the tube, will appear in both placements of the tube
	tubeLog.placeVolume( wallLog,  Transform3D() );
      }
    }  
      break;
      
    case kPunchedCenter: {
      // a volume on the z-axis with one or two inner holes
      // (implemented as a cone from which tubes are punched out)
      
      const double rUpstreamPunch = rInnerStart; // just alias names denoting what is meant here
      const double rDnstreamPunch = rInnerEnd; // (the database entries are "abused" in this case)
      
      // relative transformations for the composition of the SubtractionVolumes
      Transform3D upstreamTransformer(RotationY(-crossingAngle), Position(zPosition * tan(-crossingAngle), 0, 0));
      Transform3D dnstreamTransformer(RotationY(+crossingAngle), Position(zPosition * tan(+crossingAngle), 0, 0));
  
      // absolute transformations for the final placement in the world (angles always equal zero and 180 deg)
      Transform3D placementTransformer(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition) , rotateAngle) );
      Transform3D placementTransmirror(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition) , mirrorAngle) );
  
      // solid for the tube (including vacuum and wall): a solid cone
      ConeSegment tubeSolid( zHalf, 0, rOuterStart, 0, rOuterEnd, phi1, phi2);
        
      // tube consists of vacuum (will later have two different daughters)
      Volume tubeLog0( volName + "_0", tubeSolid, coreMaterial );
      Volume tubeLog1( volName + "_1", tubeSolid, coreMaterial );
        
      // placement of the tube in the world, both at +z and -z
      envelope.placeVolume( tubeLog0, placementTransformer );
      envelope.placeVolume( tubeLog1, placementTransmirror );

      // the wall solid and placeholders for possible SubtractionSolids
      ConeSegment wholeSolid(  zHalf, 0, rOuterStart, 0, rOuterEnd, phi1, phi2);
      
      Solid tmpSolid0, tmpSolid1, wallSolid0, wallSolid1;
      
      // the punched subtraction solids can be asymmetric and therefore have to be created twice:
      // one time in the "right" way, another time in the "reverse" way, because the "mirroring"
      // rotation around the y-axis will not only exchange +z and -z, but also +x and -x

      if ( rUpstreamPunch > 1e-6 ) { // do we need a hole on the upstream branch?
       	Tube upstreamPunch( 0, rUpstreamPunch, 5 * zHalf, phi1, phi2); // a bit longer
       	tmpSolid0 = SubtractionSolid( wholeSolid, upstreamPunch, upstreamTransformer);
       	tmpSolid1 = SubtractionSolid( wholeSolid, upstreamPunch, dnstreamTransformer); // [sic]
      } else { // dont't do anything, just pass on the unmodified shape
       	tmpSolid0 = wholeSolid;
       	tmpSolid1 = wholeSolid;
      }
      
      if (rDnstreamPunch > 1e-6 ) { // do we need a hole on the downstream branch?
       	Tube dnstreamPunch( 0, rDnstreamPunch, 5 * zHalf, phi1, phi2); // a bit longer
       	wallSolid0 = SubtractionSolid( tmpSolid0, dnstreamPunch, dnstreamTransformer);
       	wallSolid1 = SubtractionSolid( tmpSolid1, dnstreamPunch, upstreamTransformer); // [sic]
      } else { // dont't do anything, just pass on the unmodified shape
       	wallSolid0 = tmpSolid0;
       	wallSolid1 = tmpSolid1;
      }

      // the wall consists of the material given in the XML
      Volume wallLog0( volName + "_wall_0", wallSolid0, wallMaterial );
      Volume wallLog1( volName + "_wall_1", wallSolid1, wallMaterial );

      wallLog0.setVisAttributes(lcdd, "TubeVis");
      wallLog1.setVisAttributes(lcdd, "TubeVis");
      tubeLog0.setVisAttributes(lcdd, "VacVis");
      tubeLog1.setVisAttributes(lcdd, "VacVis");

      // placement as a daughter volumes of the tube
      tubeLog0.placeVolume( wallLog0, Position() );
      tubeLog1.placeVolume( wallLog1, Position() );
      
      break;
    }

    case kPunchedUpstream:
    case kPunchedDnstream: {
      // a volume on the upstream or downstream branch with two inner holes
      // (implemented as a cone from which another tube is punched out)
      
      const double rCenterPunch = (crossType == kPunchedUpstream) ? (rInnerStart) : (rInnerEnd); // just alias names denoting what is meant here
      const double rOffsetPunch = (crossType == kPunchedDnstream) ? (rInnerStart) : (rInnerEnd); // (the database entries are "abused" in this case)
      
      // relative transformations for the composition of the SubtractionVolumes
      Transform3D punchTransformer(RotationY(-2 * rotateAngle), Position(zPosition * tan(-2 * rotateAngle), 0, 0));
      Transform3D punchTransmirror(RotationY(+2 * rotateAngle), Position(zPosition * tan(+2 * rotateAngle), 0, 0));
      
      // absolute transformations for the final placement in the world
      Transform3D placementTransformer(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition) , rotateAngle) );
      Transform3D placementTransmirror(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition) , mirrorAngle) );
      
      // solid for the tube (including vacuum and wall): a solid cone
      ConeSegment tubeSolid( zHalf, 0, rOuterStart, 0, rOuterEnd, phi1, phi2);
      
      // tube consists of vacuum (will later have two different daughters)
      Volume tubeLog0( volName + "_0", tubeSolid, coreMaterial );
      Volume tubeLog1( volName + "_1", tubeSolid, coreMaterial );
      
      // placement of the tube in the world, both at +z and -z
      envelope.placeVolume( tubeLog0, placementTransformer );
      envelope.placeVolume( tubeLog1, placementTransmirror );

      // the wall solid and the piece (only a tube, for the moment) which will be punched out
      ConeSegment wholeSolid( zHalf, rCenterPunch , rOuterStart, rCenterPunch, rOuterEnd, phi1, phi2);
      
      Tube punchSolid( 0, rOffsetPunch, 5 * zHalf, phi1, phi2); // a bit longer
      
      // the punched subtraction solids can be asymmetric and therefore have to be created twice:
      // one time in the "right" way, another time in the "reverse" way, because the "mirroring"
      // rotation around the y-axis will not only exchange +z and -z, but also +x and -x
      SubtractionSolid wallSolid0( wholeSolid, punchSolid, punchTransformer);
      SubtractionSolid wallSolid1( wholeSolid, punchSolid, punchTransmirror);

      // the wall consists of the material given in the database
      Volume wallLog0( volName + "_wall_0", wallSolid0, wallMaterial );
      Volume wallLog1( volName + "_wall_1", wallSolid1, wallMaterial );
      
      wallLog0.setVisAttributes(lcdd, "TubeVis");
      wallLog1.setVisAttributes(lcdd, "TubeVis");

      tubeLog0.setVisAttributes(lcdd, "VacVis");
      tubeLog1.setVisAttributes(lcdd, "VacVis");
      
      // placement as a daughter volumes of the tube
      tubeLog0.placeVolume( wallLog0 , Position() );
      tubeLog1.placeVolume( wallLog1 , Position() );
      
      break;
    }
      
    case kUpstreamClippedFront:
    case kDnstreamClippedFront:
    case kUpstreamSlicedFront:
    case kDnstreamSlicedFront: {
      // a volume on the upstream or donwstream branch, but with the front face parallel to the xy-plane
      // or to a piece tilted in the other direction ("sliced" like a salami with 2 * rotateAngle)
      // (implemented as a slightly longer cone from which the end is clipped off)
      
      // the volume which will be used for clipping: a solid tube
      const double clipSize = rOuterStart; // the right order of magnitude for the clipping volume (alias name)
      Tube clipSolid( 0, 2 * clipSize, clipSize, phi1, phi2); // should be large enough
        
      // relative transformations for the composition of the SubtractionVolumes
      const double clipAngle = (crossType == kUpstreamClippedFront || crossType == kDnstreamClippedFront) ? (rotateAngle) : (2 * rotateAngle);
      const double clipShift = (zStart - clipSize) / cos(clipAngle) - (zPosition - clipSize / 2); // question: why is this correct?
      Transform3D clipTransformer(RotationY(-clipAngle), Position(0, 0, clipShift));
      Transform3D clipTransmirror(RotationY(+clipAngle), Position(0, 0, clipShift));
  
      // absolute transformations for the final placement in the world
      Transform3D placementTransformer(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition - clipSize / 2) , rotateAngle) );
      Transform3D placementTransmirror(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition - clipSize / 2) , mirrorAngle) );
  
      // solid for the tube (including vacuum and wall): a solid cone

      ConeSegment wholeSolid(  zHalf + clipSize / 2, 0, rOuterStart, 0, rOuterEnd,  phi1, phi2); // a bit longer
  
      // clip away the protruding end
      SubtractionSolid tubeSolid0( wholeSolid, clipSolid, clipTransformer);
      SubtractionSolid tubeSolid1( wholeSolid, clipSolid, clipTransmirror);

      // tube consists of vacuum (will later have two different daughters)
      Volume tubeLog0( volName + "_0", tubeSolid0, coreMaterial );
      Volume tubeLog1( volName + "_1", tubeSolid1, coreMaterial );
        
      // placement of the tube in the world, both at +z and -z
      envelope.placeVolume( tubeLog0, placementTransformer );
      envelope.placeVolume( tubeLog1, placementTransmirror );
       
      if (rInnerStart != rOuterStart || rInnerEnd != rOuterEnd) {
	// the wall solid: a tubular cone
	ConeSegment wallWholeSolid(  zHalf + clipSize / 2, rInnerStart, rOuterStart, rInnerEnd, rOuterEnd, phi1, phi2); // a bit longer
        
	// clip away the protruding end
	SubtractionSolid wallSolid0( wallWholeSolid, clipSolid, clipTransformer);
	SubtractionSolid wallSolid1( wallWholeSolid, clipSolid, clipTransmirror);
        
	// the wall consists of the material given in the database
	Volume wallLog0( volName + "_wall_0", wallSolid0, wallMaterial );
	Volume wallLog1( volName + "_wall_1", wallSolid1, wallMaterial );

	wallLog0.setVisAttributes(lcdd, "TubeVis");
	wallLog1.setVisAttributes(lcdd, "TubeVis");

	tubeLog0.setVisAttributes(lcdd, "VacVis");
	tubeLog1.setVisAttributes(lcdd, "VacVis");

	// placement as a daughter volumes of the tube
	tubeLog0.placeVolume( wallLog0, Position() );
	tubeLog1.placeVolume( wallLog1, Position() );
      }
    }
      break;

    case kUpstreamClippedRear:
    case kDnstreamClippedRear:
    case kUpstreamSlicedRear:
    case kDnstreamSlicedRear: {
      // a volume on the upstream or donwstream branch, but with the rear face parallel to the xy-plane
      // or to a piece tilted in the other direction ("sliced" like a salami with 2 * rotateAngle)
      // (implemented as a slightly longer cone from which the end is clipped off)
      
      // the volume which will be used for clipping: a solid tube
      const double clipSize = rOuterEnd; // the right order of magnitude for the clipping volume (alias name)
      Tube clipSolid( 0, 2 * clipSize, clipSize, phi1, phi2); // should be large enough
      
      // relative transformations for the composition of the SubtractionVolumes
      const double clipAngle = (crossType == kUpstreamClippedRear || crossType == kDnstreamClippedRear) ? (rotateAngle) : (2 * rotateAngle);
      const double clipShift = (zEnd + clipSize) / cos(clipAngle) - (zPosition + clipSize / 2); // question: why is this correct?
      Transform3D clipTransformer(RotationY(-clipAngle), Position(0, 0, clipShift));
      Transform3D clipTransmirror(RotationY(+clipAngle), Position(0, 0, clipShift));
      
      // absolute transformations for the final placement in the world
      Transform3D placementTransformer(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition + clipSize / 2) , rotateAngle) );
      Transform3D placementTransmirror(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition + clipSize / 2) , mirrorAngle) );
      
      // solid for the tube (including vacuum and wall): a solid cone
      ConeSegment wholeSolid( 0, rOuterStart, 0, rOuterEnd, zHalf + clipSize / 2, phi1, phi2); // a bit longer
      
      // clip away the protruding end
      SubtractionSolid tubeSolid0( wholeSolid, clipSolid, clipTransformer);
      SubtractionSolid tubeSolid1( wholeSolid, clipSolid, clipTransmirror);
      
      // tube consists of vacuum (will later have two different daughters)
      Volume tubeLog0( volName + "_0", tubeSolid0, coreMaterial );
      Volume tubeLog1( volName + "_1", tubeSolid1, coreMaterial );
      
      // placement of the tube in the world, both at +z and -z
      envelope.placeVolume( tubeLog0, placementTransformer );
      envelope.placeVolume( tubeLog1, placementTransmirror );
      
      if (rInnerStart != rOuterStart || rInnerEnd != rOuterEnd) {
	// the wall solid: a tubular cone
	ConeSegment wallWholeSolid( rInnerStart, rOuterStart, rInnerEnd, rOuterEnd, zHalf + clipSize / 2, phi1, phi2); // a bit longer
        
	// clip away the protruding end
	SubtractionSolid wallSolid0( wallWholeSolid, clipSolid, clipTransformer);
	SubtractionSolid wallSolid1( wallWholeSolid, clipSolid, clipTransmirror);
        
	// the wall consists of the material given in the database
	Volume wallLog0( volName + "_wall_0", wallSolid0, wallMaterial );
	Volume wallLog1( volName + "_wall_1", wallSolid1, wallMaterial );

	wallLog0.setVisAttributes(lcdd, "TubeVis");
	wallLog1.setVisAttributes(lcdd, "TubeVis");

	tubeLog0.setVisAttributes(lcdd, "VacVis");
	tubeLog1.setVisAttributes(lcdd, "VacVis");
        
	// placement as a daughter volumes of the tube
	tubeLog0.placeVolume( wallLog0, Transform3D() );
	tubeLog1.placeVolume( wallLog1, Transform3D() );
      }
      break;
    }

    case kUpstreamClippedBoth:
    case kDnstreamClippedBoth:
    case kUpstreamSlicedBoth:
    case kDnstreamSlicedBoth: {
      // a volume on the upstream or donwstream branch, but with both faces parallel to the xy-plane
      // or to a piece tilted in the other direction ("sliced" like a salami with 2 * rotateAngle)
      // (implemented as a slightly longer cone from which the end is clipped off)
      
      // the volume which will be used for clipping: a solid tube
      const double clipSize = rOuterEnd; // the right order of magnitude for the clipping volume (alias name)
      Tube clipSolid( 0, 2 * clipSize, clipSize, phi1, phi2); // should be large enough
        
      // relative transformations for the composition of the SubtractionVolumes
      const double clipAngle = (crossType == kUpstreamClippedBoth || crossType == kDnstreamClippedBoth) ? (rotateAngle) : (2 * rotateAngle);
      const double clipShiftFrnt = (zStart - clipSize) / cos(clipAngle) - zPosition;
      const double clipShiftRear = (zEnd   + clipSize) / cos(clipAngle) - zPosition;
      Transform3D clipTransformerFrnt(RotationY(-clipAngle), Position(0, 0, clipShiftFrnt));
      Transform3D clipTransformerRear(RotationY(-clipAngle), Position(0, 0, clipShiftRear));
      Transform3D clipTransmirrorFrnt(RotationY(+clipAngle), Position(0, 0, clipShiftFrnt));
      Transform3D clipTransmirrorRear(RotationY(+clipAngle), Position(0, 0, clipShiftRear));
  
      // absolute transformations for the final placement in the world
      Transform3D placementTransformer(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition) , rotateAngle) );
      Transform3D placementTransmirror(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition) , mirrorAngle) );
      
      // solid for the tube (including vacuum and wall): a solid cone
      ConeSegment wholeSolid( 0, rOuterStart, 0, rOuterEnd, zHalf + clipSize, phi1, phi2); // a bit longer
      
      // clip away the protruding ends
      SubtractionSolid tmpSolid0 ( wholeSolid, clipSolid, clipTransformerFrnt);
      SubtractionSolid tmpSolid1 ( wholeSolid, clipSolid, clipTransmirrorFrnt);
      SubtractionSolid tubeSolid0( tmpSolid0,  clipSolid, clipTransformerRear);
      SubtractionSolid tubeSolid1( tmpSolid1,  clipSolid, clipTransmirrorRear);
        
      // tube consists of vacuum (will later have two different daughters)
      Volume tubeLog0( volName + "_0", tubeSolid0, coreMaterial );
      Volume tubeLog1( volName + "_1", tubeSolid1, coreMaterial );
        
      // placement of the tube in the world, both at +z and -z
      envelope.placeVolume( tubeLog0, placementTransformer );
      envelope.placeVolume( tubeLog1, placementTransmirror );
       
      if (rInnerStart != rOuterStart || rInnerEnd != rOuterEnd) {
	// the wall solid: a tubular cone
	ConeSegment wallWholeSolid( rInnerStart, rOuterStart, rInnerEnd, rOuterEnd, zHalf + clipSize, phi1, phi2); // a bit longer
        
	// clip away the protruding ends
	SubtractionSolid wallTmpSolid0( wallWholeSolid, clipSolid, clipTransformerFrnt);
	SubtractionSolid wallTmpSolid1( wallWholeSolid, clipSolid, clipTransmirrorFrnt);
	SubtractionSolid wallSolid0   ( wallTmpSolid0,  clipSolid, clipTransformerRear);
	SubtractionSolid wallSolid1   ( wallTmpSolid1,  clipSolid, clipTransmirrorRear);
        
	// the wall consists of the material given in the database
	Volume wallLog0(volName + "_wall_0", wallSolid0, wallMaterial );
	Volume wallLog1(volName + "_wall_1", wallSolid1, wallMaterial );

	wallLog0.setVisAttributes(lcdd, "TubeVis");
	wallLog1.setVisAttributes(lcdd, "TubeVis");

	tubeLog0.setVisAttributes(lcdd, "VacVis");
	tubeLog1.setVisAttributes(lcdd, "VacVis");
        
	// placement as a daughter volumes of the tube
	tubeLog0.placeVolume( wallLog0, Transform3D() );
	tubeLog1.placeVolume( wallLog1, Transform3D() );
      }
      break;
    }
    default: {
      throw std::runtime_error( " Beampipe_o1_v01_geo.cpp : fatal failure !! ??  " ) ;

      //      return false; // fatal failure
    }

    }//end switch
  }//for all xmlSections

  //######################################################################################################################################################################
  
  tube.addExtension< DD4hep::DDRec::ConicalSupportData >( beampipeData ) ;

  //--------------------------------------
  
  Volume mother =  lcdd.pickMotherVolume( tube ) ;
  PlacedVolume pv(mother.placeVolume(envelope));
  pv.addPhysVolID( "system", xmlBeampipe.id() ) ; //.addPhysVolID("side", 0 ) ;

  tube.setVisAttributes( lcdd, xmlBeampipe.visStr(), envelope );

  tube.setPlacement(pv);
  
  return tube;
}
DECLARE_DETELEMENT(Beampipe_o1_v01,create_element)

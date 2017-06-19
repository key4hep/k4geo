//====================================================================
//  DDSim - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "XMLHandlerDB.h"
#include <cmath>
#include <map>

typedef std::map<std::string, double> TReferenceMap;
 
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

using dd4hep::rec::ConicalSupportData;

namespace {
  typedef enum { // These constants are also used in the MySQL database:
    kCenter = 0, // centered on the z-axis
    kUpstream = 1, // on the upstream branch, rotated by half the crossing angle
    kDnstream = 2, // on the downstream branch, rotated by half the crossing angle
 	
    kPunchedCenter = 3, // centered, with one or two inner holes
    kPunchedUpstream = 4, // on the upstream branch, with two inner holes
    kPunchedDnstream = 5, // on the downstrem branch, with two inner holes
 	
    kUpstreamClippedFront = 6, // upstream, with the front face parallel to the xy-plane
    kDnstreamClippedFront = 7, // downstream, with the front face parallel to the xy-plane
    kUpstreamClippedRear = 8, // upstream, with the rear face parallel to the xy-plane
    kDnstreamClippedRear = 9, // downstream, with the rear face parallel to the xy-plane
    kUpstreamClippedBoth = 10, // upstream, with both faces parallel to the xy-plane
    kDnstreamClippedBoth = 11, // downstream, with both faces parallel to the xy-plane
 	
    kUpstreamSlicedFront = 12, // upstream, with the front face parallel to a tilted piece
    kDnstreamSlicedFront = 13, // downstream, with the front face parallel to a tilted piece
    kUpstreamSlicedRear = 14, // upstream, with the rear face parallel to a tilted piece
    kDnstreamSlicedRear = 15, // downstream, with the rear face parallel to a tilted piece
    kUpstreamSlicedBoth = 16, // upstream, with both faces parallel to a tilted piece
    kDnstreamSlicedBoth = 17 // downstream, with both faces parallel to a tilted piece
  } ECrossType;
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
static Ref_t create_element(Detector& theDetector, xml_h e, SensitiveDetector /*sens*/)  {

  //------------------------------------------
  //  See comments starting with '//**' for
  //     hints on porting issues
  //------------------------------------------

  
  xml_det_t    x_det = e;
  string       name  = x_det.nameStr();
  
  //--------------------------------
  Assembly envelope( name + "_assembly"  ) ;
  //--------------------------------
  
  PlacedVolume pv;
  
  DetElement   tube(  name, x_det.id()  ) ;
  
  ConicalSupportData* beampipeData = new ConicalSupportData ;

  //######################################################################################################################################################################
  //  code ported from TubeX01::construct() :
  //##################################
  
  // double phi1 =   0 * deg;
  // double phi2 = 360 * deg;
  //** DD4hep/TGeo seems to need rad (as opposed to the manual)
  double phi1 = 0 ;
  double phi2 = 2*M_PI;
  
  //  // some visualization attributes for the tube wall and the vacuum inside
  // G4VisAttributes *wallVisAttrib = new G4VisAttributes(G4Colour(1.0, 0.7, 0.5)); // light brown
  // //wallVisAttrib->SetForceSolid(true);
  // G4VisAttributes *vacuumVisAttrib = new G4VisAttributes(G4Colour(0.0, 0.0, 0.5)); // dark blue
  // vacuumVisAttrib->SetVisibility(false); // there isn't anything, so what do you expect?
  
  const double crossingAngle = theDetector.constant<double>("ILC_Main_Crossing_Angle") / 2 ; //  only half the angle
  
  // const String dbName = env.GetDBName() + "_" + env.GetParameterAsString("ILC_Main_Crossing_Angle");
  // Database *db = new Database(dbName.c_str());
  
  bool usingOffsets = false;
  TReferenceMap referenceOffsets;

  for(xml_coll_t c( x_det ,Unicode("reference")); c; ++c)  {
    
    xml_comp_t  x_ref( c );
    XMLHandlerDB db = XMLHandlerDB( x_ref )  ;
    
    const std::string globalName   = db->fetchString("globalName");
    const std::string localName    = db->fetchString("localName");
    const double assumedValue = db->fetchDouble("assumption") ;
    const double currentValue = theDetector.constant<double>( globalName );
    const double offsetValue  = currentValue - assumedValue;
    
    referenceOffsets[localName] = offsetValue;
    
    if (offsetValue != 0) {
      cout
	<< "TubeX01: Using " << globalName << " = "
	<< currentValue / dd4hep::mm << " mm instead of "
	<< assumedValue / dd4hep::mm << " mm" << endl;
      usingOffsets = true;
    }
  }
  if (usingOffsets) std::cout << "TubeX01: Be sure you know what you're doing!" << std::endl ;
  
  
  //unused: bool saveToGear = true; // as long as this is true the vectors gearValZ, gearValRInner and gearValROuter will be filled in the coming loop
  
  // db->exec("SELECT * FROM `tube`;");
  // while (db->getTuple()) {
  for(xml_coll_t c( x_det ,Unicode("section")); c; ++c)  {
    
    xml_comp_t  x_sec( c );
    XMLHandlerDB db = XMLHandlerDB( x_sec )  ;
    
    // reference values for r- and z-values
    const std::string zStartRef         = db->fetchString("zStartRef");
    const std::string zEndRef           = db->fetchString("zEndRef");
    const std::string rInnerStartRef    = db->fetchString("rInnerStartRef");
    const std::string rInnerEndRef      = db->fetchString("rInnerEndRef");
    const std::string rOuterStartRef    = db->fetchString("rOuterStartRef");
    const std::string rOuterEndRef      = db->fetchString("rOuterEndRef");
    
    const double zStartOffset      = (zStartRef      == "") ? (0) : (referenceOffsets[zStartRef]);
    const double zEndOffset        = (zEndRef        == "") ? (0) : (referenceOffsets[zEndRef]);
    const double rInnerStartOffset = (rInnerStartRef == "") ? (0) : (referenceOffsets[rInnerStartRef]);
    const double rInnerEndOffset   = (rInnerEndRef   == "") ? (0) : (referenceOffsets[rInnerEndRef]);
    const double rOuterStartOffset = (rOuterStartRef == "") ? (0) : (referenceOffsets[rOuterStartRef]);
    const double rOuterEndOffset   = (rOuterEndRef   == "") ? (0) : (referenceOffsets[rOuterEndRef]);
  
    // fields in the data tuple
    const ECrossType crossType  = ECrossType(db->fetchInt("crossType")); // positioning of the volume
    const double zStart       = db->fetchDouble("zStart")      + zStartOffset;
    const double zEnd         = db->fetchDouble("zEnd")        + zEndOffset;
    const double rInnerStart  = db->fetchDouble("rInnerStart") + rInnerStartOffset;
    const double rInnerEnd    = db->fetchDouble("rInnerEnd")   + rInnerEndOffset;
    const double rOuterStart  = db->fetchDouble("rOuterStart") + rOuterStartOffset;
    const double rOuterEnd    = db->fetchDouble("rOuterEnd")   + rOuterEndOffset;
    const std::string materialName = db->fetchString("material");
    const std::string volName      = "tube_" + db->fetchString("name");
    
    

    if( crossType == kCenter ) { // store only the central sections !

      ConicalSupportData::Section section ;
      section.rInner = rInnerStart ;
      section.rOuter = rOuterStart ;
      section.zPos   = zStart ;
      beampipeData->sections.push_back( section ) ;
    }

    //fixme
    // if( saveToGear ){
    
    //   gearValZ.push_back( zStart );
    //   gearValRInner.push_back( rInnerStart );
    //   gearValROuter.push_back( rOuterStart );
      
    // }
    
 
    // fixme:  how to handle this export of parameters ????
    //   if(volName == "tube_IPOuterTube"){
    //     std::ostringstream oss1;
    //     oss1 << zStart;
    //     (*Control::globalModelParameters)["TUBE_IPOuterTube_start_z"] = oss1.str();
    //     std::ostringstream oss2;
    //     oss2 << zEnd;
    //     (*Control::globalModelParameters)["TUBE_IPOuterTube_end_z"] = oss2.str();
    //     std::ostringstream oss3;
    //     oss3 << rOuterStart;
    //     (*Control::globalModelParameters)["TUBE_IPOuterTube_start_radius"] = oss3.str();
    //     std::ostringstream oss4;
    //     oss4 << rOuterEnd;
    //     (*Control::globalModelParameters)["TUBE_IPOuterTube_end_radius"] = oss4.str();
   	
    //   }
    
    //   if(volName == "tube_IPOuterBulge"){

    //     std::ostringstream oss1;
    //     oss1 << zEnd;
    //     (*Control::globalModelParameters)["TUBE_IPOuterBulge_end_z"] = oss1.str();
    //     std::ostringstream oss2;
    //     oss2 << rOuterEnd;
    //     (*Control::globalModelParameters)["TUBE_IPOuterBulge_end_radius"] = oss2.str();
      
      
    //     saveToGear = false;
    //     gearValZ.push_back( zEnd );
    //     gearValRInner.push_back( rInnerEnd );
    //     gearValROuter.push_back( rOuterEnd ); 
      
    //   }
    
    // things which can be calculated immediately
    double zHalf        = fabs(zEnd - zStart) / 2; // half z length of the cone
    const double zPosition    = fabs(zEnd + zStart) / 2; // middle z position
    Material coreMaterial    = theDetector.material("beam"); // always the same
    Material wallMaterial    = theDetector.material(materialName);

    // this could mess up your geometry, so better check it
    if (crossingAngle == 0 && crossType != kCenter) {

      std::cout << "TubeX01: You are trying to build a crossing geometry without a crossing angle.\n"
	"This is probably not what you want - better check your geometry data!" << std::endl ;

      return 0 ;//false; // premature exit, Mokka will abort now
    }

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
    const double rotateAngle = tmpAngle; // for the placement at +z (better make it const now)
    const double mirrorAngle = M_PI - rotateAngle; // for the "mirrored" placement at -z
    // const double mirrorAngle = 180 * deg - rotateAngle; // for the "mirrored" placement at -z
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
      //      tubeLog->SetVisAttributes(vacuumVisAttrib);
      
      // placement of the tube in the world, both at +z and -z
      pv = envelope.placeVolume( tubeLog,  transformer );
      pv = envelope.placeVolume( tubeLog,  transmirror );
      
      // std::cout << " *** placing tubes for volume " << volName << ": "  << transformer  << transmirror 
      // 		<< " rotation : " << RotationY(rotateAngle) 	<< std::endl 
      // 		<< " angle : " << rotateAngle	<< std::endl 
      // 		<< " zPos: " << zPosition	<< std::endl ;

      // if inner and outer radii are equal, then omit the tube wall
      if (rInnerStart != rOuterStart || rInnerEnd != rOuterEnd) {
	
	// the wall solid: a tubular cone
	ConeSegment wallSolid( zHalf, rInnerStart, rOuterStart, rInnerEnd, rOuterEnd, phi1, phi2);
        
	// the wall consists of the material given in the database
	Volume wallLog ( volName + "_wall", wallSolid, wallMaterial);
	
	tube.setVisAttributes(theDetector, "TubeVis"  , wallLog );
	
	// placement as a daughter volume of the tube, will appear in both placements of the tube
	pv = tubeLog.placeVolume( wallLog,  Transform3D() );
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

      // tubeLog0->SetVisAttributes(vacuumVisAttrib);
      // tubeLog1->SetVisAttributes(vacuumVisAttrib);
      // //FIXME debug
      // tube.setVisAttributes(theDetector, "TubeVis"  , tubeLog0 );
      // tube.setVisAttributes(theDetector, "TubeVis"  , tubeLog1 );
        
      // placement of the tube in the world, both at +z and -z
      pv = envelope.placeVolume( tubeLog0, placementTransformer );
      pv = envelope.placeVolume( tubeLog1, placementTransmirror );
      // new PVPlacement(placementTransformer, tubeLog0, volName, worldLog, false, 0);
      // new PVPlacement(placementTransmirror, tubeLog1, volName, worldLog, false, 1);

      // std::cout << " *** placing tubes for volume " << volName << ": "  << placementTransformer  << placementTransmirror 
      // 		<< " rotation : " << RotationY(rotateAngle) 	<< std::endl 
      // 		<< " angle : " << rotateAngle	<< std::endl 
      // 		<< " zPos: " << zPosition	<< std::endl ;


        
      // the wall solid and placeholders for possible SubtractionSolids
      ConeSegment wholeSolid(  zHalf, 0, rOuterStart, 0, rOuterEnd, phi1, phi2);
      
      Solid tmpSolid0, tmpSolid1, wallSolid0, wallSolid1;
      // VSolid *tmpSolid0, *tmpSolid1, *wallSolid0, *wallSolid1;
      
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


      // the wall consists of the material given in the database
      Volume wallLog0( volName + "_wall_0", wallSolid0, wallMaterial );
      Volume wallLog1( volName + "_wall_1", wallSolid1, wallMaterial );
      
      tube.setVisAttributes(theDetector, "TubeVis"  , wallLog0);
      tube.setVisAttributes(theDetector, "TubeVis"  , wallLog1);
      
      // placement as a daughter volumes of the tube
      pv = tubeLog0.placeVolume( wallLog0, Position() );
      pv = tubeLog1.placeVolume( wallLog1, Position() );
      // new PVPlacement(0, Position(), wallLog0, volName + "_wall", tubeLog0, false, 0);
      // new PVPlacement(0, Position(), wallLog1, volName + "_wall", tubeLog1, false, 1);
      
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

      // tubeLog0->SetVisAttributes(vacuumVisAttrib);
      // tubeLog1->SetVisAttributes(vacuumVisAttrib);
      
      // placement of the tube in the world, both at +z and -z
      pv = envelope.placeVolume( tubeLog0, placementTransformer );
      pv = envelope.placeVolume( tubeLog1, placementTransmirror );
      // new PVPlacement(placementTransformer, tubeLog0, volName, worldLog, false, 0);
      // new PVPlacement(placementTransmirror, tubeLog1, volName, worldLog, false, 1);
      // std::cout << " *** placing tubes for volume " << volName << ": "  << placementTransformer  << placementTransmirror 
      // 		<< " rotation : " << RotationY(rotateAngle) 	<< std::endl 
      // 		<< " angle : " << rotateAngle	<< std::endl 
      // 		<< " zPos: " << zPosition	<< std::endl ;

      
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
      
      tube.setVisAttributes(theDetector, "TubeVis"  , wallLog0 );
      tube.setVisAttributes(theDetector, "TubeVis"  , wallLog1 );
      
      
      // placement as a daughter volumes of the tube
      pv = tubeLog0.placeVolume( wallLog0 , Position() );
      pv = tubeLog1.placeVolume( wallLog1 , Position() );
      // new PVPlacement(0, Position(), wallLog0, volName + "_wall", tubeLog0, false, 0);
      // new PVPlacement(0, Position(), wallLog1, volName + "_wall", tubeLog1, false, 1);
      
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
        // tubeLog0->SetVisAttributes(vacuumVisAttrib);
        // tubeLog1->SetVisAttributes(vacuumVisAttrib);
        
        // placement of the tube in the world, both at +z and -z
	pv = envelope.placeVolume( tubeLog0, placementTransformer );
	pv = envelope.placeVolume( tubeLog1, placementTransmirror );
        // new PVPlacement(placementTransformer, tubeLog0, volName, worldLog, false, 0);
        // new PVPlacement(placementTransmirror, tubeLog1, volName, worldLog, false, 1);
      // std::cout << " *** placing tubes for volume " << volName << ": "  << placementTransformer  << placementTransmirror 
      // 		<< " rotation : " << RotationY(rotateAngle) 	<< std::endl 
      // 		<< " angle : " << rotateAngle	<< std::endl 
      // 		<< " zPos: " << zPosition	<< std::endl ;
       
        if (rInnerStart != rOuterStart || rInnerEnd != rOuterEnd) {
          // the wall solid: a tubular cone
          ConeSegment wallWholeSolid(  zHalf + clipSize / 2, rInnerStart, rOuterStart, rInnerEnd, rOuterEnd, phi1, phi2); // a bit longer
        
          // clip away the protruding end
          SubtractionSolid wallSolid0( wallWholeSolid, clipSolid, clipTransformer);
          SubtractionSolid wallSolid1( wallWholeSolid, clipSolid, clipTransmirror);
        
          // the wall consists of the material given in the database
          Volume wallLog0( volName + "_wall_0", wallSolid0, wallMaterial );
          Volume wallLog1( volName + "_wall_1", wallSolid1, wallMaterial );

	  tube.setVisAttributes(theDetector, "TubeVis"  , wallLog0 );
	  tube.setVisAttributes(theDetector, "TubeVis"  , wallLog1 );
        
          // placement as a daughter volumes of the tube
	  pv = tubeLog0.placeVolume( wallLog0, Position() );
	  pv = tubeLog1.placeVolume( wallLog1, Position() );
          // new PVPlacement(0, Position(), wallLog0, volName + "_wall", tubeLog0, false, 0);
          // new PVPlacement(0, Position(), wallLog1, volName + "_wall", tubeLog1, false, 1);
        }
        break;
      }
	
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
      // tubeLog0->SetVisAttributes(vacuumVisAttrib);
      // tubeLog1->SetVisAttributes(vacuumVisAttrib);
      
      // placement of the tube in the world, both at +z and -z
      pv = envelope.placeVolume( tubeLog0, placementTransformer );
      pv = envelope.placeVolume( tubeLog1, placementTransmirror );
      // new PVPlacement(placementTransformer, tubeLog0, volName, worldLog, false, 0);
      // new PVPlacement(placementTransmirror, tubeLog1, volName, worldLog, false, 1);
      // std::cout << " *** placing tubes for volume " << volName << ": "  << placementTransformer  << placementTransmirror 
      // 		<< " rotation : " << RotationY(rotateAngle) 	<< std::endl 
      // 		<< " angle : " << rotateAngle	<< std::endl 
      // 		<< " zPos: " << zPosition	<< std::endl ;
      
      if (rInnerStart != rOuterStart || rInnerEnd != rOuterEnd) {
	// the wall solid: a tubular cone
	ConeSegment wallWholeSolid( rInnerStart, rOuterStart, rInnerEnd, rOuterEnd, zHalf + clipSize / 2, phi1, phi2); // a bit longer
        
	// clip away the protruding end
	SubtractionSolid wallSolid0( wallWholeSolid, clipSolid, clipTransformer);
	SubtractionSolid wallSolid1( wallWholeSolid, clipSolid, clipTransmirror);
        
	// the wall consists of the material given in the database
	Volume wallLog0( volName + "_wall_0", wallSolid0, wallMaterial );
	Volume wallLog1( volName + "_wall_1", wallSolid1, wallMaterial );

	tube.setVisAttributes(theDetector, "TubeVis"  , wallLog0 );
	tube.setVisAttributes(theDetector, "TubeVis"  , wallLog1 );
        
	// placement as a daughter volumes of the tube
	pv = tubeLog0.placeVolume( wallLog0, Transform3D() );
	pv = tubeLog1.placeVolume( wallLog1, Transform3D() );
	// new PVPlacement(0, Position(), wallLog0, volName + "_wall", tubeLog0, false, 0);
	// new PVPlacement(0, Position(), wallLog1, volName + "_wall", tubeLog1, false, 1);
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
      // tubeLog0->SetVisAttributes(vacuumVisAttrib);
      // tubeLog1->SetVisAttributes(vacuumVisAttrib);
        
      // placement of the tube in the world, both at +z and -z
      pv = envelope.placeVolume( tubeLog0, placementTransformer );
      pv = envelope.placeVolume( tubeLog1, placementTransmirror );
      // new PVPlacement(placementTransformer, tubeLog0, volName, worldLog, false, 0);
      // new PVPlacement(placementTransmirror, tubeLog1, volName, worldLog, false, 1);
      // std::cout << " *** placing tubes for volume " << volName << ": "  << placementTransformer  << placementTransmirror 
      // 		<< " rotation : " << RotationY(rotateAngle) 	<< std::endl 
      // 		<< " angle : " << rotateAngle	<< std::endl 
      // 		<< " zPos: " << zPosition	<< std::endl ;
       
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

	tube.setVisAttributes(theDetector, "TubeVis"  , wallLog0 );
	tube.setVisAttributes(theDetector, "TubeVis"  , wallLog1 );
        
	// placement as a daughter volumes of the tube
	pv = tubeLog0.placeVolume( wallLog0, Transform3D() );
	pv = tubeLog1.placeVolume( wallLog1, Transform3D() );
	// new PVPlacement(0, Position(), wallLog0, volName + "_wall", tubeLog0, false, 0);
	// new PVPlacement(0, Position(), wallLog1, volName + "_wall", tubeLog1, false, 1);
      }
      break;
    }
    default: {
      //        Control::Log("TubeX01: Unimplemented \"crossType\" code.");
      return 0 ; // fatal failure
    }
    } // switch (crossType)
  } // while (db->getTuple())
  
  
  
  //######################################################################################################################################################################
  
  tube.addExtension< ConicalSupportData >( beampipeData ) ;

  //--------------------------------------
  
  Volume mother =  theDetector.pickMotherVolume( tube ) ;
  
  pv = mother.placeVolume(envelope);

  pv.addPhysVolID( "system", x_det.id() ) ; //.addPhysVolID("side", 0 ) ;
  
  tube.setVisAttributes( theDetector, x_det.visStr(), envelope );
  //  if( tube.isValid() ) 
  tube.setPlacement(pv);
  
  return tube;
}
DECLARE_DETELEMENT(TubeX01,create_element)

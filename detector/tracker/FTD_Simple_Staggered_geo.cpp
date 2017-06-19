//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "XMLHandlerDB.h"
#include "XML/Utilities.h"

#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "FTD_Simple_Staggered.h"

// #define DEBUG_VALUES
// #define DEBUG_PETAL 4

//#include "DDRec/DDGear.h"
//#define MOKKA_GEAR

#include <cmath>
#include <string>
#include <map>
#include <assert.h>

using namespace std;
using namespace dd4hep;
using namespace rec;

/// helper to wrap access to global constants
struct EnvDetector{
  Detector* _theDetector;
  EnvDetector() : _theDetector(0) {} ;
  EnvDetector(Detector& theDetector) : _theDetector( &theDetector) {} ;
  inline double GetParameterAsDouble(const std::string& name) const {  return _theDetector->constant<double>( name ) ; } 
} ;


/// Helper vector for creation of DetElements holding Volumes and their placement 
typedef std::vector< std::pair< Volume, PlacedVolume > > VolVec ;

///Helper vector for DetElements (e.g. petals in disk)
typedef std::vector< DetElement > DEVec ;


/// Storing all the gear parameters
std::map<int,std::vector<double> > _ftdparameters;


/// Global environment variables
glEnviron _glEnv;
/// common DB Parameters
dbInfoCommon _dbParCommon;
/// disk DB Parameters
dbInfoDisk _dbParDisk;
/// reco DB Parameters
dbExtended_reconstruction_parameters _dbParExReco;
double _inner_radius = 0.0;
double _outer_radius = 0.0;
double _z_position = 0.0;
double _beamTubeRadius = 0.0;
double _zEnd = 0.0;
Material _SiMat ;
Material _KaptonMat ;
Material _CuMat ;
Material _AirMat ;
Material _CarbonFiberMat ;


// function prototpyes
double Getdy(const double & innerRadius );
double Getdx( const double & innerRadius );
//void DoAndPlaceDisk( Detector& theDetector,DetElement det,SensitiveDetector sens, std::map<std::string,double> valuesDict, Volume  mother ) ;
void petalSupport( Detector& theDetector,DetElement det,  std::map<std::string,double> valuesDict, Volume  FTDPetalAirLogical ) ;
VolVec  petalSensor(  Detector& theDetector, DetElement ftd, SensitiveDetector sens, std::map<std::string,double> valuesDict, Volume  FTDPetalAirLogical ) ;
Trap SemiPetalSolid(const double& petal_cp_support_dy, const std::string& whereItgoes, const bool isSilicon ) ;


//debug print function
void printVolume( Volume v ){
  std::cout << " ++++ create Volume " << v.name() << std::endl ;
}
//=========================== PARAMETERS SETTERS FUNCTIONS ====================================/
//*********************************************************************************************
// Set Environment variables (dependent of other subdetectors)
void SetEnvironPar(const EnvDetector& env)
{

  _glEnv.TPC_Ecal_Hcal_barrel_halfZ = env.GetParameterAsDouble("TPC_Ecal_Hcal_barrel_halfZ") ;
  _glEnv.Ecal_endcap_zmin = env.GetParameterAsDouble("Ecal_endcap_zmin")  ;
  _glEnv.TPC_inner_radius = env.GetParameterAsDouble("TPC_inner_radius")  ;
  
  _glEnv.SIT1_Half_Length_Z = env.GetParameterAsDouble("SIT1_Half_Length_Z") ;
  _glEnv.SIT2_Half_Length_Z = env.GetParameterAsDouble("SIT2_Half_Length_Z") ;
  _glEnv.SIT1_Radius = env.GetParameterAsDouble("SIT1_Radius") ;
  _glEnv.SIT2_Radius = env.GetParameterAsDouble("SIT2_Radius") ;
  _glEnv.VXD_layer3_maxZ = env.GetParameterAsDouble("VXD_length_r3") ;
  
  _glEnv.zEnd_IPOuterTube = (env.GetParameterAsDouble("TUBE_IPOuterTube_end_z"))  ;  // ---> A lo mejor no hacen falta
  _glEnv.rEnd_IPOuterTube = (env.GetParameterAsDouble("TUBE_IPOuterTube_end_radius")) ;
  _glEnv.zEnd_IPOuterBulge = (env.GetParameterAsDouble("TUBE_IPOuterBulge_end_z")) ;
  _glEnv.rEnd_IPOuterBulge = (env.GetParameterAsDouble("TUBE_IPOuterBulge_end_radius")) ;
  
  _glEnv.beamTubeTangent = ( _glEnv.rEnd_IPOuterBulge - _glEnv.rEnd_IPOuterTube ) / (_glEnv.zEnd_IPOuterBulge - _glEnv.zEnd_IPOuterTube);
}


//*********************************************************************************************
// Set variables common to all disk, dumping 'common_parameters' table from 'ftd08' database
void SetdbParCommon(xml_comp_t x_det)
{
  // Getting common_parameters table
  //db->exec("select * from common_parameters;");
  { 
    XMLHandlerDB db(  x_det.child( _Unicode( common_parameters ) ) );
    
    //  db->getTuple();
    
    _dbParCommon.beamTubeClearance = db->fetchDouble("beamtube_clearance") ; 
    _dbParCommon.outer_cylinder_total_thickness = db->fetchDouble("outer_cylinder_total_thickness") ;
    _dbParCommon.inner_cylinder_total_thickness = _dbParCommon.outer_cylinder_total_thickness;
    _dbParCommon.cable_shield_thickness = db->fetchDouble("cable_shield_thickness") ;
    _dbParCommon.cables_thickness = db->fetchDouble("cables_thickness") ;
    
    // check that there is enough space for the cables and support
    if( _dbParCommon.beamTubeClearance < (_dbParCommon.cables_thickness + (2.0*_dbParCommon.cable_shield_thickness) + 0.5 *mm) )  
      {
	cout << "FTD_Simple_Staggered:Stop: Not enough space for inner support structure and cables: increase beamTubeClearance" << endl;
	exit(1);
      }
    
    _dbParCommon.ftd1_vtx3_distance_z =  db->fetchDouble("ftd1_vtx3_distance_z") ; 
    _dbParCommon.ftd7_ecal_distance_z =  db->fetchDouble("ftd7_ecal_distance_z") ; 
    _dbParCommon.ftd1_sit1_radial_diff =  db->fetchDouble("ftd1_sit1_radial_diff") ; 
    _dbParCommon.ftd2_sit1_radial_diff =  db->fetchDouble("ftd2_sit1_radial_diff") ; 
    _dbParCommon.ftd3_sit2_radial_diff =  db->fetchDouble("ftd3_sit2_radial_diff") ; 
    _dbParCommon.ftd4to7_tpc_radial_gap =  db->fetchDouble("ftd4to7_tpc_radial_gap") ; 
    // Petal Central Part Support: X-Y dimensions, thickness and angles. Same constant values
    // for all the micro-strips disks
    _dbParCommon.petal_half_angle_support = db->fetchDouble("petal_half_angle_support") ;
    _dbParCommon.petal_y_ratio = db->fetchDouble("petal_y_ratio") ;
    //fg: add additional parameter:
    _dbParCommon.support_spaceframe_width= db->fetchDouble("support_spaceframe_width") ;

  } 

  // db->exec("select * from extended_reconstruction_parameters;");
  // db->getTuple();
  {
    XMLHandlerDB db(  x_det.child( _Unicode( extended_reconstruction_parameters ) ) );
    
    _dbParExReco.strip_width  = db->fetchDouble("strip_width")   ;
    _dbParExReco.strip_length = db->fetchDouble("strip_length")  ;
    _dbParExReco.strip_pitch  = db->fetchDouble("strip_pitch")   ;
    _dbParExReco.strip_angle  = db->fetchDouble("strip_angle")   ;
  }

#ifdef DEBUG_VALUES
  cout << "FTD_Simple_Staggered:SetdbParCommon:\n" 
       << "beamTubeClearance = " << _dbParCommon.beamTubeClearance << " \n"   
       << "outer_cylinder_total_thickness = " << _dbParCommon.outer_cylinder_total_thickness << " \n"   
       << "inner_cylinder_total_thickness = " << _dbParCommon.inner_cylinder_total_thickness << " \n"   
       << "cable_shield_thickness = " << _dbParCommon.cable_shield_thickness << " \n"   
       << "cables_thickness = " << _dbParCommon.cables_thickness << " \n"   
       << "ftd1_vtx3_distance_z = " << _dbParCommon.ftd1_vtx3_distance_z << " \n"   
       << "ftd7_ecal_distance_z = " << _dbParCommon.ftd7_ecal_distance_z << " \n"   
       << "ftd1_sit1_radial_diff = " << _dbParCommon.ftd1_sit1_radial_diff << " \n"   
       << "ftd2_sit1_radial_diff = " << _dbParCommon.ftd2_sit1_radial_diff << " \n"   
       << "ftd3_sit2_radial_diff = " << _dbParCommon.ftd3_sit2_radial_diff << " \n"   
       << "ftd4to7_tpc_radial_gap = " << _dbParCommon.ftd4to7_tpc_radial_gap << " \n"   
       << "petal_half_angle_support = " << _dbParCommon.petal_half_angle_support << " \n"   
       << "petal_y_ratio = " << _dbParCommon.petal_y_ratio << " \n"   
       << "strip_width = " << _dbParExReco.strip_width << " \n"   
       << "strip_length = " << _dbParExReco.strip_length << " \n"   
       << "strip_pitch = " << _dbParExReco.strip_pitch << " \n"   
       << "strip_angle = " << _dbParExReco.strip_angle << " \n"   
       << endl;

#endif
  
}

//*********************************************************************************************
// Set variables disk number specific, dumping 'disk' table from 'ftd08' database
void SetParDisk( XMLHandlerDB db)
{

  _dbParDisk.disk_number = db->fetchInt( "disk_number" );
  _dbParDisk.disks_Si_thickness = db->fetchDouble("disk_si_thickness")  ;
  _dbParDisk.petal_cp_support_thickness = db->fetchDouble("petal_cp_support_thickness")  ;
  _dbParDisk.petal_cp_support_dxMax = db->fetchDouble("petal_cp_support_dxMax") ; 
  _dbParDisk.petal_support_zoffset = db->fetchDouble("petal_support_zoffset") ; //NEW
  _dbParDisk.sensor_is_pixel = db->fetchInt("sensor_is_pixel"); //NEW
  _dbParDisk.double_sided = db->fetchInt("double_sided"); //NEW
  
#ifdef DEBUG_VALUES
  cout << "FTD_Simple_Staggered:SetParDisk:\n" 
       << "disk_number = " << _dbParDisk.disk_number << " \n"   
       << "sensor_is_pixel = " << _dbParDisk.sensor_is_pixel << " \n"   
       << "double_sided = " << _dbParDisk.double_sided << " \n"   
       << "petal_support_zoffset = " << _dbParDisk.petal_support_zoffset << " \n"   
       << "disks_Si_thickness = " << _dbParDisk.disks_Si_thickness << " \n"   
       << "petal_cp_support_thickness = " << _dbParDisk.petal_cp_support_thickness << " \n"   
       << "petal_cp_support_dxMax = " << _dbParDisk.petal_cp_support_dxMax << " \n"   
       << "petal_support_zoffset = " << _dbParDisk.petal_support_zoffset << " \n"   
       << endl;
#endif
  
}
//=END======================= PARAMETERS SETTERS FUNCTIONS ================================END=/



/** Construction of FTD detector, ported from Mokka driver FTD_simple_Staggered.cc
 *
 * FTD_Simple_Staggered.cc
 *
 * Simplified Implementation of a self scaling 7 disk FTD
 * Based on SFtd06 but using simple trapezoid sensitive and support structures
 * All disks have the same structure
 * Sensitive material Silicon 
 * Support material Carbon Fiber, Foam 
 *
 * All disks' envelop:
 * Dimensions and coordinates are specified for the sensitive layer, support disks are built on to these
 * _inner_radius = (  _beamTubeRadius + beamTubeClearance)
 *
 * First Disk:
 * z defined by distance from end of VTX layer 3
 * outer r defined by radial difference to SIT layer 1
 *
 * Second Disk:
 * z defined relative to TPC half-length: to ensure positioning with SIT set these numbers to the same value in DB
 * outer r defined by radial difference to SIT layer 1
 *
 * Third Disk:
 * z defined relative to TPC half-length: to ensure positioning with SIT set these numbers to the same value in DB
 * outer r defined by radial difference to SIT layer 1
 *
 * Fourth, Fifth and Sixth Disk:
 * z defined relative to TPC half-length
 * outer r defined by gap between TPC inner radius and FTD disks
 *
 * Last Disk:
 * z defined by distance from front of ECal endcap
 * outer r defined by gap between TPC inner radius and FTD disks
 *
 * Parameters Set in Model Parameter DB Table:
 * TPC_Ecal_Hcal_barrel_halfZ
 * _glEnv.Ecal_endcap_zmin
 * _glEnv.TPC_inner_radius
 * VXD_length_r3
 *
 * Parameters shared with other drivers:
 * SSit03:_glEnv.SIT1_Half_Length_Z
 * SSit03:_glEnv.SIT2_Half_Length_Z 
 * SSit03:_glEnv.SIT1_Radius 
 * SSit03:_glEnv.SIT2_Radius 
 * TubeX01:TUBE_IPOuterTube_end_z
 * TubeX01:TUBE_IPOuterTube_end_radius
 * TubeX01:TUBE_IPOuterBulge_end_z
 * TubeX01:TUBE_IPOuterBulge_end_radius
 *
 *
 * History:
 *  May 2014: FG: original port from Mokka
 *  Oct 2014  FG: - added class ZDiskPetalsData as interface to reconstruction (initially to 
 *                  instantiate gear::FTDParameters )  
 *                - changed rotations and translations to a more natural one (positive sense of rotation, etc.)
 *                - 'replaced' usage of reflection in original code with individual placements
 *                  in distinct disks on either side of the origin (see comments in code for placements
                    of 
 * 
 * Mokka History:  
 * - first implementation P. Mora de Freitas (sept 02)
 * - fixed geometry overlap -- Adrian Vogel, 2005-12-05
 * - implemented new GEAR interface -- K. Harder, T. Pinto Jayawardena  2007-07-31
 * - SFtd03: Modified version of SFtd02: Rewritten as a self scaling driver which does not
 *   make use of a seperate super driver. Steve Aplin (May 2008)
 *   October 15th 2008, Steve Aplin using description from SiLC Collaboration
 * - Fixes a bug in SFtd04.cc which meant that the copper cables on thin 
 *   inside of the cylinder were far to thick (SJA 28/05/09)
 *   September 7th 2010, Jordi Duarte using mechanical design from IFCA group
 * - SFtd06: Modified version of SFtd05 implementing realistic details of the disks 
 *           4,5,6,7 structure. -- J. Duarte Campderros (Sept. 2010)
 *           Added realistic description to disks 1,2,3. Changed disk 3 to micro-strips 
 *           technology --- J. Duarte Campderros (Oct. 2010)
 *           Included the alternative z-offset between petals --|
 *           Included the use of the GEAR class FTDParameter  --| J. Duarte Campderros (July, 2011)
 *           Modified the placement of the Volumes using the 
 *           G4ReflectionFactory Place methods. Now the volumes in Z
 *           negatives are specular images from the positives -- J. Duarte Campderros (Sept, 2011)
 * - FTD_Simple_Staggered: created to enable development of tracking code while SFtd06 is finalised.
 *   S.J. Aplin (Nov 2011)
 *
 *  @author: F.Gaede, DESY, May 2014
 *  @version $Id$
 */
static Ref_t create_element(Detector& theDetector, xml_h e, SensitiveDetector sens)  {


  //  std::cout << " FTD04 - Detector.BuildType = " << theDetector.buildType() << std::endl ;

  //------------------------------------------
  //  See comments starting with '//**' for
  //     hints on porting issues
  //------------------------------------------

  
  xml_det_t    x_det = e;
  string       name  = x_det.nameStr();
  
  DetElement   ftd(  name, x_det.id()  ) ;

 // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , ftd ) ;
  
  dd4hep::xml::setDetectorTypeFlag( e, ftd ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return ftd ;
  
  //-----------------------------------------------------------------------------------
  
  PlacedVolume pv;

  sens.setType("tracker");

  // --- create assembly and DetElement for support and service volumes 

  Assembly supp_assembly( name + "_support_assembly"  ) ;

  pv = envelope.placeVolume( supp_assembly ) ;

  DetElement suppDE( ftd , name+"_support" , x_det.id() )  ;
  suppDE.setPlacement( pv ) ;
  //--------------------------------


  dd4hep::rec::ZDiskPetalsData*  zDiskPetalsData = new ZDiskPetalsData ;

 //######################################################################################################################################################################
  //  code ported from FTD_Simple_Staggered::construct() :
  //##################################
  
  double phi1 = 0 ;
  double phi2 = 2*M_PI;
  
  //fg: vis attributes defined in xml now
  //     // Globals and cosmetics definitions
  // 	_VisAttSensitive = new G4VisAttributes(G4Colour(1.,1.,.45));
  //   //G4VisAttributes *
  // 	_VisAttSupport = new G4VisAttributes(G4Colour(1,.5,.5));
  // 	_VisAttHolePetal = new G4VisAttributes(G4Color(1,.5,.5,1.0));
  // 	G4VisAttributes *VisAttAirDisk = new G4VisAttributes(G4Colour(.5,.3,.8,0.98));
  // 	VisAttAirDisk->SetVisibility(0);
  // 	G4VisAttributes *VisAttAirPetal = new G4VisAttributes(G4Colour(.5,.3,.8,0.98));
  // 	VisAttAirPetal->SetVisibility(0);
  // 	G4VisAttributes *VisAttCyl = new G4VisAttributes(G4Colour(0.45,.2,0.9,.98));
  //   G4VisAttributes *VisAttCables = new G4VisAttributes(G4Colour(0.,0.9,0.));
  //   VisAttCables->SetForceWireframe(false);
  
  // 	PhysicalVolumesPair Phys;
  
  // Get and set the Globals from the surrounding environment TPC ECAL SIT VTX and Beam-Pipe
  SetEnvironPar(   EnvDetector( theDetector ) );
	
  // Get and set the variables global to the FTD cables_thickness, ftd1_vtx3_distance_z, etc
  //  Database * db = new Database(env.GetDBName());
  SetdbParCommon( x_det );
  
  // Materials definitions
  _SiMat     = theDetector.material("G4_Si") ; // silicon_2.33gccm");
  _KaptonMat = theDetector.material("G4_KAPTON"); //kapton");
  _CuMat     = theDetector.material("G4_Cu"); //copper"); 
  _AirMat    = theDetector.material("G4_AIR" ); //air");
  _CarbonFiberMat = theDetector.material("CarbonFiber");
 

  //fg: replace with standard from materials file
  // -- PROVISIONAL -- Carbon Fiber definition from database ?? 
  // double density;
  // std::string matname, symbol;
  // int nel;
  // double fractionmass, volumefraction;
  // volumefraction = 0.5;
  // density = (1.3 + volumefraction / 3 ) * g/cm3;
  // fractionmass = 1 - 1.3 * (1 - volumefraction) / (density / (g/cm3));
  // _CarbonFiberMat = new Material(matname = "CarbonFiber", density, nel=2);
  // _CarbonFiberMat->AddElement(CGAGeometryManager::GetElement("C"), fractionmass);
  // _CarbonFiberMat->AddMaterial(CGAGeometryManager::GetMaterial("epoxy"),1.0-fractionmass);
  // cout << "CarbonFiber->GetRadlen() = " << _CarbonFiberMat->GetRadlen() / cm << " cm" << endl;
  // <-------------------------------------------------------------------------------------------
  
  
	
  cout << "FTD_Simple_Staggered:"  
       << "\t inner support thickness = " << _dbParCommon.inner_cylinder_total_thickness  
       << "\t cables thickness = " << _dbParCommon.cables_thickness
       << "\t 2 x cable shield thickness = " << 2 * _dbParCommon.cable_shield_thickness
       << "\t beamTubeClearance = " << _dbParCommon.beamTubeClearance
       << endl;
  
  // Now we can start to build the disks -------------------------------------------
  const double theta = _dbParCommon.petal_half_angle_support;
  
  // Disk parameters
  _dbParDisk.ZStartOuterCylinder=0;
  _dbParDisk.ZStopOuterCylinder=0;
  double OuterCylinderInnerRadius=0;
  
  _dbParDisk.ZStartInnerCylinder=0;
  _dbParDisk.ZStopInnerCylinder=0;
  double InnerCylinderOuterRadius1=0;
  double InnerCylinderOuterRadius2=0;
  
  //   db->exec("select * from disks;");
  //   db->getTuple();
  //   //... assembling detector
  //   do {

  for(xml_coll_t c( x_det ,_U(disk)); c; ++c)  {
    
    xml_comp_t  x_disk( c );
    XMLHandlerDB db( x_disk )  ;

    //... 
    int disk_number(-1);
    // double _inner_radius = 0.0;
    // double _outer_radius = 0.0;
    // double _z_position = 0.0;
    // double _beamTubeRadius = 0.0;
    // double _zEnd = 0.0;

    // Get and set the parameters disk specific
    SetParDisk( db );
    
    disk_number = _dbParDisk.disk_number;
#ifdef ONE_DISK
    if (disk_number != ONE_DISK )
      {
	continue;
      }
#endif
		
    switch (disk_number) 
      {
      case 1:
	// z defined by distance from end of VTX layer 3
	_z_position = ( _glEnv.VXD_layer3_maxZ + _dbParCommon.ftd1_vtx3_distance_z );
          
	//          _z_position = disk_number * 100.0 * mm;
          
	// outer r defined by radial difference to SIT layer 1
	_outer_radius = ( _glEnv.SIT1_Radius + _dbParCommon.ftd1_sit1_radial_diff ); 
          
          
	// beam tube radius at backside of disk 
	_zEnd = _z_position +  _dbParDisk.petal_support_zoffset + 0.5 * _dbParDisk.petal_cp_support_thickness + ( _dbParDisk.double_sided * _dbParDisk.disks_Si_thickness ) ;
          
	// check which part of the beam tube this disk lies above
	_beamTubeRadius = (_zEnd < _glEnv.zEnd_IPOuterTube ) ? _glEnv.rEnd_IPOuterTube : _glEnv.rEnd_IPOuterTube + ( (_zEnd - _glEnv.zEnd_IPOuterTube ) * _glEnv.beamTubeTangent );
          
	_inner_radius = (  _beamTubeRadius + _dbParCommon.beamTubeClearance);
          
	// check that there is no overlap with SIT1
	if( _z_position <= _glEnv.SIT1_Half_Length_Z && _outer_radius>=_glEnv.SIT1_Radius) 
	  {
            cout << "FTD_Simple_Staggered:Stop: Overlap between FTD1 and SIT1" << endl;
            cout << "FTD_Simple_Staggered:FTD1 Radius = " << _outer_radius << "SIT1 Radius = " << _glEnv.SIT1_Radius << endl;
            exit(1);
	  }
	if( db->fetchDouble("z_position_ReltoTPCLength") != 0.0) 
	  {
            cout << "FTD_Simple_Staggered:Stop: The z position of FTD1 is not relative. The relative value will not be used. It should be set to 0.0 in the DB." << endl;
            cout << "FTD_Simple_Staggered:Stop: The z position of FTD1 is set by the distance between the centre of the sensitive layer and the max z of VTX layer 3." << endl;
            exit(1);
	  }
	break;
          
      case 2:
	// z defined relative to TPC half-length: to ensure positioning with SIT set these numbers to the same value in DB
	_z_position = (_glEnv.TPC_Ecal_Hcal_barrel_halfZ * db->fetchDouble("z_position_ReltoTPCLength")) ;
          
	//          _z_position = disk_number * 100.0 * mm;
          
	// outer r defined by radial difference to SIT layer 1
	_outer_radius = _glEnv.SIT1_Radius + _dbParCommon.ftd2_sit1_radial_diff; 
          
	// beam tube radius at backside of disk 
	_zEnd = _z_position +  _dbParDisk.petal_support_zoffset + 0.5 * _dbParDisk.petal_cp_support_thickness + ( _dbParDisk.double_sided * _dbParDisk.disks_Si_thickness ) ;
          
	// check which part of the beam tube this disk lies above
	_beamTubeRadius = (_zEnd < _glEnv.zEnd_IPOuterTube ) ? _glEnv.rEnd_IPOuterTube : _glEnv.rEnd_IPOuterTube + ( (_zEnd - _glEnv.zEnd_IPOuterTube ) * _glEnv.beamTubeTangent );
          
	_inner_radius = (  _beamTubeRadius + _dbParCommon.beamTubeClearance) ;
          
	//... keep information for inner support cylinder with 0.5mm saftey clearance from inner radius of disks
	_dbParDisk.ZStartInnerCylinder = _glEnv.zEnd_IPOuterTube;
          
	InnerCylinderOuterRadius1 = _inner_radius - ( ( _zEnd - _glEnv.zEnd_IPOuterTube ) * _glEnv.beamTubeTangent ) - 0.5 * mm; 
          
	// check that there is no overlap with SIT1
	if( _z_position <= _glEnv.SIT1_Half_Length_Z && _outer_radius>=_glEnv.SIT1_Radius) 
	  {
            cout << "FTD_Simple_Staggered:Stop:Overlap between FTD2 and SIT1" << endl;
            cout << "FTD_Simple_Staggered:FTD2 Radius = " << _outer_radius << "SIT1 Radius = " << _glEnv.SIT1_Radius << endl;
            exit(1);
	  }
	break;
	
      case 3:
	// z defined relative to TPC half-length: to ensure positioning with SIT set these numbers to the same value in DB
	_z_position = (_glEnv.TPC_Ecal_Hcal_barrel_halfZ * db->fetchDouble("z_position_ReltoTPCLength")) ;
          
	//          _z_position = disk_number * 100.0 * mm;
          
	// outer r defined by radial difference to SIT layer 2
	_outer_radius = _glEnv.SIT2_Radius + _dbParCommon.ftd3_sit2_radial_diff; 
          
	// beam tube radius at backside of disk 
	_zEnd = _z_position +  _dbParDisk.petal_support_zoffset + 0.5 * _dbParDisk.petal_cp_support_thickness + ( _dbParDisk.double_sided * _dbParDisk.disks_Si_thickness ) ;
          
	// check which part of the beam tube this disk lies above
	_beamTubeRadius = (_zEnd < _glEnv.zEnd_IPOuterTube ) ? _glEnv.rEnd_IPOuterTube : _glEnv.rEnd_IPOuterTube + ( (_zEnd - _glEnv.zEnd_IPOuterTube ) * _glEnv.beamTubeTangent );
          
	_inner_radius = _beamTubeRadius + _dbParCommon.beamTubeClearance ;
          
	// check that there is no overlap with SIT1
	if( _z_position <= _glEnv.SIT2_Half_Length_Z && _outer_radius>=_glEnv.SIT2_Radius) 
	  {
            cout << "FTD_Simple_Staggered:Stop:Overlap between FTD3 and SIT2" <<  endl;
            cout << "FTD_Simple_Staggered:FTD3 Radius = " << _outer_radius << "SIT2 Radius = " << _glEnv.SIT2_Radius << endl;
            exit(1);
	  }
	break;
	
      case 4:
      case 5:
      case 6:
	// z defined relative to TPC half-length
	_z_position = (_glEnv.TPC_Ecal_Hcal_barrel_halfZ * db->fetchDouble("z_position_ReltoTPCLength")) ;
          
	//          _z_position = disk_number * 100.0 * mm;
          
	// outer r defined by gap between TPC inner radius and FTD disks
	_outer_radius = _glEnv.TPC_inner_radius - _dbParCommon.ftd4to7_tpc_radial_gap; 
          
	// beam tube radius at backside of disk 
	_zEnd = _z_position +  _dbParDisk.petal_support_zoffset + 0.5 * _dbParDisk.petal_cp_support_thickness + ( _dbParDisk.double_sided * _dbParDisk.disks_Si_thickness ) ;

          
	// check which part of the beam tube this disk lies above
	_beamTubeRadius = (_zEnd < _glEnv.zEnd_IPOuterTube ) ? _glEnv.rEnd_IPOuterTube : _glEnv.rEnd_IPOuterTube + ( (_zEnd - _glEnv.zEnd_IPOuterTube ) * _glEnv.beamTubeTangent );
          
	_inner_radius = _beamTubeRadius + _dbParCommon.beamTubeClearance ;
          
	// keep the information for outer cylinder
	if(disk_number==4)
	  {
            _dbParDisk.ZStartOuterCylinder = _z_position;
	  }
	break;
          
      case 7:
	// z defined by distance from front of ECal endcap
	_z_position = _glEnv.Ecal_endcap_zmin - _dbParCommon.ftd7_ecal_distance_z;
          
	//          _z_position = disk_number * 100.0 * mm;
          
	// outer r defined by gap between TPC inner radius and FTD disks
	_outer_radius = _glEnv.TPC_inner_radius - _dbParCommon.ftd4to7_tpc_radial_gap; 
          
	// beam tube radius at backside of disk 
	_zEnd = _z_position +  _dbParDisk.petal_support_zoffset + 0.5 * _dbParDisk.petal_cp_support_thickness + ( _dbParDisk.double_sided * _dbParDisk.disks_Si_thickness ) ;
          
	// check which part of the beam tube this disk lies above
	_beamTubeRadius = (_zEnd < _glEnv.zEnd_IPOuterTube ) ? _glEnv.rEnd_IPOuterTube : _glEnv.rEnd_IPOuterTube + ( (_zEnd - _glEnv.zEnd_IPOuterTube ) * _glEnv.beamTubeTangent );
          
	_inner_radius = _beamTubeRadius + _dbParCommon.beamTubeClearance ;
          
	// End of Support Structure: 0.5mm clearance from disks
	_dbParDisk.ZStopOuterCylinder = _zEnd;
	_dbParDisk.ZStopInnerCylinder = _zEnd;
          
	OuterCylinderInnerRadius = _outer_radius + 0.5 * mm;
	InnerCylinderOuterRadius2 = _inner_radius - 0.5 * mm; 
          
	if( db->fetchDouble("z_position_ReltoTPCLength") != 0.0) 
	  {
            cout << "FTD_Simple_Staggered:Stop: The z position of FTD7 is not relative. The relative value will not be used. It should be set to 0.0 in the DB." << endl;
            cout << "FTD_Simple_Staggered:Stop: The z position of FTD7 is set by the distance between the centre of the sensitive layer and the min z of the ECal Endcap." << endl;
            exit(1);
	  }
	break;
          
      default:
	cout << "FTD_Simple_Staggered: Error disk number must be between 1-7: disk number = " << disk_number << endl;
	exit(1);
      }
    
    cout << "FTD_Simple_Staggered: Disk:" << disk_number
	 << "\t z = " << _z_position
	 << "\t inner rad = " << _inner_radius
	 << "\t outer rad = " << _outer_radius
	 << "\t beamtube rad = " << _beamTubeRadius
	 << "\t free space = " << (_inner_radius - 0.5 * mm - _dbParCommon.inner_cylinder_total_thickness - (2*_dbParCommon.cable_shield_thickness) - _dbParCommon.cables_thickness) - _beamTubeRadius 
	 << endl;
		
    
    /**************************************************************************************
     ** Begin construction of disks with appropiate parameters    **
     **************************************************************************************/
    
    //================================== AIR DISK =======================================//
    //  The air-disk is the container and mother volume of the petals. There will be
    //  7x2 air disks copies placed in the world volume.
    //  
    //  Check the comments at the beginning of this file for the description of 
    //  each # disk parameter.
    //  
    //  Input parameters:
    //       _inner_radius: inner radius of the whole structure
    //       _outer_radius: outer radius of the whole structure
    //       max_half_thickness_disk: 
    //           the maximum thickness of the disk = 2.0 * ( sensitive thickness + support thickness + Zoffset ) 
    
    //                                 Zoffset=the displacement of the disks in z-direction
    
    // 
    //		// The thickness of the air petal (containing the support and sensors)
    //		double petalairthickness_half = 0.5 * ( _dbParDisk.petal_cp_support_thickness +
    //			_dbParDisk.disks_Si_thickness )  
    
    
    
    // need enough space for double sided 
    double petalairthickness_half = 0.5 * ( _dbParDisk.petal_cp_support_thickness
					    + 2.0*_dbParDisk.disks_Si_thickness ) ;
    
    double max_half_thickness_disk = _dbParDisk.petal_support_zoffset + petalairthickness_half ;

		
    // Tubs *FTDDiskSolid = new Tubs("FTDAirDiskSolid",
    // 				  _inner_radius,
    // 				  _outer_radius,
    // 				  max_half_thickness_disk,
    // 				  phi1,
    // 				  phi2
    // 				  );
    
    Tube FTDDiskSolid( _inner_radius, _outer_radius, max_half_thickness_disk, phi1, phi2 );

    // LogicalVolume *FTDDiskLogical = new LogicalVolume(FTDDiskSolid,
    // 						      _AirMat,
    // 						      "FTDAirDiskLogical", 
    // 						      0, 
    // 						      0, 
    // 						      0);

    //fg: Volume FTDDiskLogical(  _toString(  _dbParDisk.disk_number, "FTDAirDiskLogical_%d" ), FTDDiskSolid, _AirMat ) ;
    //fg: replace the logical volume for the disks with two individual ones for the pos. and neg. z axis respectively
    //fg: this way we do not need a reflection and can position the petals with different transforms on either side...
    Volume FTDDiskLogicalPZ(  _toString(  _dbParDisk.disk_number, "FTDAirDiskLogicalPZ_%d" ), FTDDiskSolid, _AirMat ) ;
    Volume FTDDiskLogicalNZ(  _toString(  _dbParDisk.disk_number, "FTDAirDiskLogicalNZ_%d" ), FTDDiskSolid, _AirMat ) ;


    //    FTDDiskLogical->SetVisAttributes(VisAttAirDisk);
    //    ftd.setVisAttributes(theDetector,  "SeeThrough", FTDDiskLogical ) ;
    ftd.setVisAttributes(theDetector,  "SeeThrough", FTDDiskLogicalPZ ) ;
    ftd.setVisAttributes(theDetector,  "SeeThrough", FTDDiskLogicalNZ ) ;

		
    
    // RotationMatrix *rotDiskPositive = new RotationMatrix();
    // // Sensors facing the IP)
    // rotDiskPositive->rotateY(pi);
    // // Re-allocating the local disk frame to the global frame
    // rotDiskPositive->rotateZ(-pi/2.0);
    
    // Transform3D transPositive( *rotDiskPositive, Position( 0.,0.,_z_position) );

    // // Place the positive copy in the world
    // Phys = ReflectionFactory::Instance()->Place( transPositive,
    // 						 "FTDAirDisk",
    // 						 FTDDiskLogical,
    // 						 worldLog,
    // 						 false,
    // 						 disk_number);
    // registerPV( Phys );

    ///fg RotationZYX rotDiskPositive( -pi/2.0, pi , 0. ) ; 
    //fg use unrotated air disks...
    RotationZYX rotDiskPositive(0,0,0) ; 
    Transform3D transPositive( rotDiskPositive,  Position( 0.,0.,_z_position) );      

    pv = envelope.placeVolume( FTDDiskLogicalPZ, transPositive ) ;

    DetElement   diskDEposZ( ftd ,   _toString(  _dbParDisk.disk_number, "FTDDisk_%d_posZ" ) , x_det.id() );
    diskDEposZ.setPlacement( pv ) ;

    pv.addPhysVolID("layer", disk_number - 1  ).addPhysVolID("side", 1 )   ;

    
#ifdef DEBUG_VALUES
    cout << "===================================================================== " << "\n" <<
      "FTDAirDisk:\n" << 
      " Inner Radius= " << _inner_radius <<  "\n" <<
      " Outer Radius= " << _outer_radius <<  "\n" <<
      " thickness =   " << max_half_thickness_disk*2.0 << "\n" <<
      " placed at \n" << 
      " x =   " <<  transPositive.Translation().Vect().X() << "\n" <<
      " y =   " <<  transPositive.Translation().Vect().Y() << "\n" <<
      " z =   " <<  transPositive.Translation().Vect().Z() << "\n" <<
      endl;
#endif
    
    
//     // Place negative copy
#ifndef DEBUG_POSITIVE
    // 		RotationMatrix *rotDiskNegative = new RotationMatrix();
    // 		rotDiskNegative->rotateZ(-pi/2.0);
    
    // 		Transform3D transNegative( *rotDiskNegative, Position( 0.,0.,-_z_position ) );
    //     //Specular image
    // 		transNegative = transNegative*ReflectX3D();
    
    // 		Phys = ReflectionFactory::Instance()->Place( transNegative,
    //                                                   "FTDAirDisk",
    //                                                   FTDDiskLogical,
    //                                                   worldLog,
    //                                                   false,
    //                                                   -disk_number);
    // 		registerPV( Phys );
    
    
    //fg RotationZYX rotDiskNegative( -pi/2.0, 0 , 0. ) ; 
    //fg use unrotated air disks...
    RotationZYX rotDiskNegative(0,0,0) ; 
    Transform3D transNegative( rotDiskNegative,  Position( 0.,0., -_z_position) );      
    pv = envelope.placeVolume( FTDDiskLogicalNZ, transNegative ) ;

    DetElement   diskDEnegZ( ftd ,   _toString(  _dbParDisk.disk_number, "FTDDisk_%d_negZ" ) , x_det.id() );
    diskDEnegZ.setPlacement( pv ) ;


    pv.addPhysVolID("layer", disk_number -1  ).addPhysVolID("side", -1 )   ;

    
    
    
#ifdef DEBUG_VALUES
    cout << "===================================================================== " << "\n" <<
      "FTDAirDisk:\n" << 
      " Inner Radius= " << _inner_radius <<  "\n" <<
      " Outer Radius= " << _outer_radius <<  "\n" <<
      " thickness =   " << max_half_thickness_disk*2.0 << "\n" <<
      " placed at \n" << 
      " x =   " <<  transNegative.Translation().Vect().X() << "\n" <<
      " y =   " <<  transNegative.Translation().Vect().Y() << "\n" <<
      " z =   " <<  transNegative.Translation().Vect().Z() << "\n" <<
      endl;
#endif
    
#endif
    //=END=============================== AIR DISK  =================================END=/
    
    //=================================== AIR PETAL =====================================/
    // Air container for the petal: the mother of the real support petal and the silicon 
    // sensors. This air petal will be placed inside the Air Disk,
    // generating N rotated copies along the z-axis.               
    //  Input parameters:     dxMax                                _
    //                      --------                              | |    
    //                      \      /   |                          | |              
    //       XY-Plane        \    /    | dy          YZ-Plane     | |    
    //                        \__/     |                          |_|     
    //                        dxMin                                dz
    // 
    //                     dxMax: given by the database
    //                     dxMin: depends of the _inner_radius of each disk
    //                     dy:    heigth, depends of each disk
    //                     dz:    thickness of the supports + thickness of Si
    //                     theta: given by the db, semi-angle which defines the trapezoid
		
    // Dimensions for the disk
 
    const double petal_cp_supp_half_dxMin = Getdx( _inner_radius )/2.0;
    const double petal_cp_support_dy = Getdy(_inner_radius);
 
    // ------------------------------------------------------------------------
 
#ifdef DEBUG_VALUES
    std::cout << "*** Petal parameters : petal_cp_supp_half_dxMin=" << petal_cp_supp_half_dxMin
	      << " _dbParDisk.petal_cp_support_dxMax/2.0 =" << _dbParDisk.petal_cp_support_dxMax/2.0
	      << " petal_cp_support_dy/2.0 =" << petal_cp_support_dy/2.0
	      << " petalairthickness_half =" << petalairthickness_half << std::endl ;
#endif

    Trap FTDPetalAirSolid( petalairthickness_half, //thickness (calculated in the disk zone)
    			   0.0,
    			   0.0,
    			   petal_cp_support_dy/2.0,  // dy
    			   petal_cp_supp_half_dxMin, //dxMin 
    			   _dbParDisk.petal_cp_support_dxMax/2.0, //dxMax
    			   0.0,
    			   petal_cp_support_dy/2.0,  // dy
    			   petal_cp_supp_half_dxMin,  // dxMin
    			   _dbParDisk.petal_cp_support_dxMax/2.0, //dxMax
    			   0.0);
 
    Volume FTDPetalAirLogical( _toString(  _dbParDisk.disk_number, "FTDPetalAirLogical_%d" ) , FTDPetalAirSolid, _AirMat ) ;
    //printVolume( FTDPetalAirLogical ) ; 

    ftd.setVisAttributes(theDetector,  "SeeThrough" , FTDPetalAirLogical ) ;

		
    // Placing N-copies of the air petal inside the air disk. The copies are built using the z-axis as
    // the axis of rotation
    const int petal_max_number = (int)(360.0*deg/(2.0*theta)) ; 

    DEVec petVecposZ(petal_max_number) ;
    DEVec petVecnegZ(petal_max_number) ;

    double zSignPetal0 = 1. ;

    for (int i = 0; i < petal_max_number; i++){ 

#ifdef DEBUG_PETAL
      if(i != DEBUG_PETAL ) {
	continue;
      }
#endif

      // Put the petal in the position inside the disk
      double petalCdtheta = i*2.0*theta;
      //fg: changed this to positve rotation
            
      //fg: the petals at negative z are numbered in positive sense of rotation
      //    for the ones at positve z, we have to flip around the petal (rotY(pi)) 
      //    and then rotate in the negative direction around the new (inverted) z-axis
      //    -> this really mimicks the reflection of the petal disk:
      //       changed sense of rotation and flipped petals !!
      RotationZ   rotPetalNZ(    petalCdtheta - pi/2.             );
      RotationZYX rotPetalPZ( -( petalCdtheta - pi/2. ) , pi , 0. );

      int zsign = pow((double)-1,i);
     
      // Petal i=0 parameters for gear
      if( i == 0 ) {
	_ftdparameters[gearpar::PHI0].push_back(petalCdtheta );
	_ftdparameters[gearpar::PETAL0SIGNOFFSET].push_back(zsign);
	zSignPetal0 = zsign ;
      }
#ifdef DEBUG_PETAL
      _ftdparameters[gearpar::PHI0].push_back(0);
      _ftdparameters[gearpar::PETAL0SIGNOFFSET].push_back(1);
#endif
      
      //fg: exchanged sin() and cos() in order to have a normal positve sense 
      //    of rotation around z-axis
      double dx = (petal_cp_support_dy/2.0 + _inner_radius)*cos( petalCdtheta );
      double dy = (petal_cp_support_dy/2.0 + _inner_radius)*sin( petalCdtheta ); 
      double dz = zsign*( _dbParDisk.petal_support_zoffset) ;
      
      Transform3D transPetalPZ( rotPetalPZ, Position( dx, dy, -dz) );

      //fg: at negative z we just exchange the sign of the z-offset
      dx = (petal_cp_support_dy/2.0 + _inner_radius)*cos( petalCdtheta );
      dy = (petal_cp_support_dy/2.0 + _inner_radius)*sin( petalCdtheta ); 
      Transform3D transPetalNZ( rotPetalNZ, Position( dx, dy,  dz) );
      
      
      // Phys = ReflectionFactory::Instance()->Place(
      // 						  transPetal,
      // 						  "FTDPetalAir",
      // 						  FTDPetalAirLogical,
      // 						  FTDDiskLogical,
      // 						  false,
      // 						  i+1);
      // registerPV( Phys );

      //      pv = FTDDiskLogical.placeVolume( FTDPetalAirLogical, transPetal ) ;

      // create DetElements for every petal
      std::stringstream sspz ;  sspz << "ftd_petal_posZ_" << disk_number << "_"  << i  ;
      std::stringstream ssnz ;  ssnz << "ftd_petal_negZ_" << disk_number << "_"  << i  ;
      DetElement petalDEposZ( diskDEposZ, sspz.str() ,  x_det.id() );
      DetElement petalDEnegZ( diskDEnegZ, ssnz.str() ,  x_det.id() );

      pv = FTDDiskLogicalPZ.placeVolume( FTDPetalAirLogical, transPetalPZ ) ;
      pv.addPhysVolID("module", i ) ;
      petalDEposZ.setPlacement( pv ) ;


      pv = FTDDiskLogicalNZ.placeVolume( FTDPetalAirLogical, transPetalNZ ) ;
      pv.addPhysVolID("module", i ) ;
      petalDEnegZ.setPlacement( pv ) ;

      petVecposZ[i] =  petalDEposZ ;
      petVecnegZ[i] =  petalDEnegZ ;


#ifdef DEBUG_VALUES
      cout << "===================================================================== " << "\n" <<
	"FTDPetalAir:\n" << 
	" Petal Offset = " << zsign*_dbParDisk.petal_support_zoffset << 
	" Inner Radius= " << _inner_radius <<  "\n" <<
	" Outer Radius= " << _outer_radius <<  "\n" <<
	" xMax = " << _dbParDisk.petal_cp_support_dxMax <<  "\n" <<
	" xMin = " << 2.0*petal_cp_supp_half_dxMin << "\n" <<
	" dy =   " << petal_cp_support_dy << "\n" <<
	" thickness =   " << petalairthickness_half*2.0 << "\n" <<
	" placed at (pos z)  \n" << 
	" x =   " <<  transPetalPZ.Translation().Vect().X() << "\n" <<
	" y =   " <<  transPetalPZ.Translation().Vect().Y() << "\n" <<
	" z =   " <<  transPetalPZ.Translation().Vect().Z() << "\n" <<
	endl;
#endif
      
    } // end petal loop ...  
    
    
    // -------- reconstruction parameters  ----------------
    dd4hep::rec::ZDiskPetalsData::LayerLayout thisLayer ;
    
    int isDoubleSided = false;
    int nSensors = 1;

    if( _dbParDisk.sensor_is_pixel != 1 ) {
      isDoubleSided = true;
      nSensors = 2;
    }

    thisLayer.typeFlags[ dd4hep::rec::ZDiskPetalsData::SensorType::DoubleSided ] = isDoubleSided ;
    thisLayer.typeFlags[ dd4hep::rec::ZDiskPetalsData::SensorType::Pixel ]	   = _dbParDisk.sensor_is_pixel ;

    thisLayer.petalHalfAngle	  = _dbParCommon.petal_half_angle_support ;
    thisLayer.alphaPetal	  = 0. ;	// petals are othogonal to z-axis
    thisLayer.zPosition		  = _z_position ;
    thisLayer.petalNumber	  = petal_max_number ;
    thisLayer.sensorsPerPetal	  = nSensors ; 
    thisLayer.phi0		  = 0.  ;
    thisLayer.zOffsetSupport	  = - zSignPetal0 *  fabs( _dbParDisk.petal_support_zoffset ) ; // sign of offset is negative (!?)
    thisLayer.distanceSupport	  = _inner_radius ;
    thisLayer.thicknessSupport	  = _dbParDisk.petal_cp_support_thickness ;
    thisLayer.widthInnerSupport	  = 2. * petal_cp_supp_half_dxMin ;
    thisLayer.widthOuterSupport	  = _dbParDisk.petal_cp_support_dxMax ;
    thisLayer.lengthSupport	  = petal_cp_support_dy ;
    thisLayer.zOffsetSensitive	  = zSignPetal0 * ( fabs( _dbParDisk.petal_support_zoffset ) +  0.5 * (_dbParDisk.disks_Si_thickness+_dbParDisk.petal_cp_support_thickness)  )  ;
    thisLayer.distanceSensitive	  = _inner_radius ; 
    thisLayer.thicknessSensitive  = _dbParDisk.disks_Si_thickness ;
    thisLayer.widthInnerSensitive =  2. * petal_cp_supp_half_dxMin ;
    thisLayer.widthOuterSensitive = _dbParDisk.petal_cp_support_dxMax ;
    thisLayer.lengthSensitive	  = petal_cp_support_dy ;
    
    zDiskPetalsData->layers.push_back( thisLayer ) ;

    // -------- end reconstruction parameters  ----------------

  
#ifdef DD4HEP_WITH_GEAR // ------------------------ Gear disk parameters
    int sensorType = gear::FTDParameters::PIXEL;  
    int isDoubleSided = false;
    int nSensors = 1;
    if( _dbParDisk.sensor_is_pixel != 1 ) {
      
      sensorType = gear::FTDParameters::STRIP;
      isDoubleSided = true;
      nSensors = 2;
    }
#ifdef DEBUG_PETAL
    _ftdparameters[gearpar::NPETALS].push_back(1);
#else
    _ftdparameters[gearpar::NPETALS].push_back(petal_max_number);
#endif
    _ftdparameters[gearpar::SENSORTYPE].push_back(sensorType);
    _ftdparameters[gearpar::ISDOUBLESIDED].push_back(isDoubleSided);
    _ftdparameters[gearpar::NSENSORS].push_back(nSensors);
    _ftdparameters[gearpar::ZPOSITION].push_back(_z_position);
    _ftdparameters[gearpar::ZOFFSET].push_back(_dbParDisk.petal_support_zoffset);		
    _ftdparameters[gearpar::ALPHA].push_back(0.0); // staggered design has no tilt
    _ftdparameters[gearpar::HALFANGLEPETAL].push_back(_dbParCommon.petal_half_angle_support);
    
#endif
    //=END=============================== AIR PETAL =================================END=/
    
    //=========================== PETALS & SENSORS ==============================/ 
    
    /******************************************************
     ** Support, sensors and electronics are built via   **
     ** DoAnPlaceDisk, see the appropiate functions:     **
     **                                                  **  
     **   +---------------------++--------------------+  **
     **   |    Petal Supports   ||       sensors      |  **
     **   +---------------------++--------------------+  **
     **   | petalSupportPixels  || pixelSensors       |  **
     **   +---------------------++--------------------+  **
     **                                                  **
     **                                                  **
     **                                                  **
     ******************************************************/
    
    std::map<std::string,double> valuesDict;
    
    valuesDict["petal_cp_supp_half_dxMin"] = petal_cp_supp_half_dxMin;
    valuesDict["petal_cp_support_dy"] = petal_cp_support_dy;
    valuesDict["_inner_radius"] = _inner_radius;
    

      
    //    DoAndPlaceDisk( theDetector, ftd, sens, valuesDict, FTDPetalAirLogical );		

    petalSupport(theDetector, ftd, valuesDict, FTDPetalAirLogical ) ; 

    VolVec volV = petalSensor( theDetector, ftd, sens, valuesDict, FTDPetalAirLogical );


    //---- meassurement surface vectors 

    Vector3D u0( -1. , 0. ,  0. ) ;
    Vector3D v0(  0. , 1. ,  0. ) ;
    Vector3D n0(  0. , 0. , -1. ) ;

    Vector3D u1( -1. , 0. ,  0. ) ;
    Vector3D v1(  0. , 1. ,  0. ) ;
    Vector3D n1(  0. , 0. , -1. ) ;

    
    double supp_thick = _dbParDisk.petal_cp_support_thickness ;
    double active_silicon_thickness =  _dbParDisk.disks_Si_thickness  ;

    SurfaceType surfType(SurfaceType::Sensitive) ;


    if( ! _dbParDisk.sensor_is_pixel ){  // strip sensor
      
      surfType.setProperty( SurfaceType::Measurement1D , true ) ;
      
      // implement stereo angle 
      double strip_angle  = _dbParExReco.strip_angle  ;
      
      // choose the rotation here such that u x v = n
      
      u0.fill( -cos( strip_angle ) ,  sin( strip_angle  ) , 0. ) ;
      v0.fill(  sin( strip_angle ) ,  cos( strip_angle  ) , 0. ) ;

      u1.fill( -cos( strip_angle ) , -sin( strip_angle  ) , 0. ) ;
      v1.fill( -sin( strip_angle ) ,  cos( strip_angle  ) , 0. ) ;

    }

    // surf0 is used for the first sensor and includes the complete support material - surf1 is used for the second sensor and has only the silicon
    VolPlane surf0( volV[0].first , surfType , active_silicon_thickness/2 , active_silicon_thickness/2 + supp_thick,  u0,v0,n0 ) ;
    VolPlane surf1( volV[1].first , surfType , active_silicon_thickness/2 , active_silicon_thickness/2             ,  u1,v1,n1 ) ; ;

    //----


    // create DetElements for every sensor and assign to the petal DEs
    // one or two (for double layers ) for positve and negative z each 
    for (int i = 0; i < petal_max_number; i++){ 

#ifdef DEBUG_PETAL
      if(i != DEBUG_PETAL ) {
	continue;
      }
#endif
      //create DetElements for sensors - one per sensitive petal
      std::stringstream sspz ;  sspz << "ftd_sensor_posZ_" << disk_number << "_"  << i << "_0"  ;
      std::stringstream ssnz ;  ssnz << "ftd_sensor_negZ_" << disk_number << "_"  << i << "_0"  ;
      
      DetElement sensorDEposZ( petVecposZ[i], sspz.str() ,  x_det.id() );
      DetElement sensorDEnegZ( petVecnegZ[i], ssnz.str() ,  x_det.id() );

      sensorDEposZ.setPlacement( volV[0].second ) ;
      sensorDEnegZ.setPlacement( volV[0].second ) ;
      
      volSurfaceList( sensorDEposZ )->push_back( surf0 ) ;
      volSurfaceList( sensorDEnegZ )->push_back( surf0 ) ;

      if(_dbParDisk.double_sided == 1 ) { // first two disks are single sided pixel

	std::stringstream sspz1 ;  sspz1 << "ftd_sensor_posZ_" << disk_number << "_"  << i << "_1"  ;
	std::stringstream ssnz1 ;  ssnz1 << "ftd_sensor_negZ_" << disk_number << "_"  << i << "_1"  ;

	DetElement sensorDEposZ_2( petVecposZ[i], sspz1.str() ,  x_det.id() );
	DetElement sensorDEnegZ_2( petVecnegZ[i], ssnz1.str() ,  x_det.id() );

	sensorDEposZ_2.setPlacement( volV[1].second ) ;
	sensorDEnegZ_2.setPlacement( volV[1].second ) ;

	volSurfaceList( sensorDEposZ_2 )->push_back( surf1 ) ;
	volSurfaceList( sensorDEnegZ_2 )->push_back( surf1 ) ;
      }
    }

    //=END======================= PETALS, SENSORS & ELECT. ==========================END=/ 
    
    
  } //**************** LOOP over disks *******************************


	
  //================================ OUTER CYLINDER ==================================/

#ifndef DEBUG_PETAL
#ifndef ONE_DISK

  assert(_dbParDisk.ZStartOuterCylinder>0);
  assert(_dbParDisk.ZStopOuterCylinder>0);
  
  double OuterCylinder_half_z = (_dbParDisk.ZStopOuterCylinder-_dbParDisk.ZStartOuterCylinder)/2.;
  assert(OuterCylinder_half_z>0);
  
  double OuterCylinder_position = _dbParDisk.ZStartOuterCylinder + OuterCylinder_half_z;
  
  Tube FTDOuterCylinderSolid(OuterCylinderInnerRadius,
			     OuterCylinderInnerRadius+_dbParCommon.outer_cylinder_total_thickness,
			     OuterCylinder_half_z,
			     phi1, 
			     phi2);
  
  Volume FTDOuterCylinderLogical("FTDOuterCylinder", FTDOuterCylinderSolid, _KaptonMat ) ;

  ftd.setVisAttributes( theDetector, "FTDCylVis", FTDOuterCylinderLogical ) ;
	
  Transform3D transCylPlus(  RotationZYX() , Position(0.,0.,OuterCylinder_position));
  Transform3D transCylMinus( RotationZYX(), Position(0.,0.,-OuterCylinder_position));

  //fixme: do we need a special transform here ?
  //       nothing is placed anyways - see below...
  // transCylMinus = transCylMinus*ReflectZ3D();
  
  //	Phys= ReflectionFactory::Instance()->Place(
  //                                               transCylPlus,
  //                                               "FTDOuterCylinder",
  //                                               FTDOuterCylinderLogical,
  //                                               worldLog,
  //                                               false,
  //                                               0);      
  //	registerPV( Phys );
  //  
  //	Phys= ReflectionFactory::Instance()->Place(
  //                                               transCylMinus,
  //                                               "FTDOuterCylinder",
  //                                               FTDOuterCylinderLogical,
  //                                               worldLog,
  //                                               false,
  //                                               0);      
  //	registerPV( Phys );
  //=END============================ OUTER CYLINDER =============================END==/
	
  //================================ INNER CYLINDER ==================================/
  //... Inner cylinder (cone)
  assert(_dbParDisk.ZStartInnerCylinder>0);
  assert(_dbParDisk.ZStopInnerCylinder>0);
  
  double InnerCylinder_half_z =  (_dbParDisk.ZStopInnerCylinder-_dbParDisk.ZStartInnerCylinder)/2.;
  assert(InnerCylinder_half_z>0);
  
  //double InnerCylinder_position = _dbParDisk.ZStartInnerCylinder + InnerCylinder_half_z; NOT USED
  
  double InnerCylinderRmin1 = InnerCylinderOuterRadius1 - _dbParCommon.inner_cylinder_total_thickness - (2.0*_dbParCommon.cable_shield_thickness) - _dbParCommon.cables_thickness ;
  double InnerCylinderRmax1 = InnerCylinderOuterRadius1;
  double InnerCylinderRmin2 = InnerCylinderOuterRadius2 - _dbParCommon.inner_cylinder_total_thickness - (2.0*_dbParCommon.cable_shield_thickness) - _dbParCommon.cables_thickness ;
  double InnerCylinderRmax2 = InnerCylinderOuterRadius2;
	
  double cableShieldRmin1 = InnerCylinderRmin1;  double cableShieldRmax1 = cableShieldRmin1 + (2.0*_dbParCommon.cable_shield_thickness) + _dbParCommon.cables_thickness ;
  double cableShieldRmin2 = InnerCylinderRmin2;
  double cableShieldRmax2 = cableShieldRmin2 + (2.0*_dbParCommon.cable_shield_thickness) + _dbParCommon.cables_thickness;
	
  double cablesRmin1 = cableShieldRmin1 + _dbParCommon.cable_shield_thickness; 
  double cablesRmax1 = cablesRmin1 + _dbParCommon.cables_thickness;
  double cablesRmin2 = cableShieldRmin2 + _dbParCommon.cable_shield_thickness; 
  double cablesRmax2 = cablesRmin2 + _dbParCommon.cables_thickness;
  
  ConeSegment FTDInnerCylinderSolid( InnerCylinder_half_z, 
				     InnerCylinderRmin1,
				     InnerCylinderRmax1,
				     InnerCylinderRmin2,
				     InnerCylinderRmax2,		 
				     phi1, 
				     phi2);
  
  Volume FTDInnerCylinderLogical("FTDInnerCylinder", FTDInnerCylinderSolid, _KaptonMat ) ;

  ftd.setVisAttributes( theDetector, "FTDCylVis", FTDInnerCylinderLogical ) ;
  
  ConeSegment FTDCableShieldSolid( InnerCylinder_half_z,
				   cableShieldRmin1,
				   cableShieldRmax1,
				   cableShieldRmin2,
				   cableShieldRmax2,
				   phi1, 
				   phi2);
  
  Volume FTDCableShieldLogical( "FTDInnerCableShield", FTDCableShieldSolid, _KaptonMat ) ;
			        
  ftd.setVisAttributes( theDetector, "FTDCylVis",   FTDCableShieldLogical );
  
  ConeSegment FTDCablesSolid( InnerCylinder_half_z,
			      cablesRmin1,
			      cablesRmax1,
			      cablesRmin2,
			      cablesRmax2,
			      phi1, 
			      phi2);
  
  Volume FTDCablesLogical("FTDInnerCables", FTDCablesSolid, _CuMat ) ;
			   
  ftd.setVisAttributes( theDetector, "FTDCylVis",  FTDCablesLogical ) ;
  
  // fg: the placements are all commented out in the Mokka original
  //     so no outer cylinder is created - do we need one ???

  //  //... the cables are placed inside the cylinder
  //	Phys = ReflectionFactory::Instance()->Place(
  //                                                Transform3D(),
  //                                                "FTDInnerCables",
  //                                                FTDCablesLogical,
  //                                                FTDCableShieldLogical,
  //                                                false,
  //                                                0);      
  //	registerPV( Phys );
  //	
  //	Phys = ReflectionFactory::Instance()->Place(
  //                                                Transform3D(),
  //                                                "FTDInnerCableShield",
  //                                                FTDCableShieldLogical,
  //                                                FTDInnerCylinderLogical,
  //                                                false,
  //                                                0);      
  //	registerPV( Phys );
  //  
  //	Phys = ReflectionFactory::Instance()->Place(
  //                                                Transform3D(RotationMatrix(), 
  //                                                              Position(0., 0., InnerCylinder_position) ),
  //                                                "FTDInnerCylinder",
  //                                                FTDInnerCylinderLogical,
  //                                                worldLog,
  //                                                false,
  //                                                0);
  //	registerPV( Phys );
  //  
  //  Transform3D Tcyl( RotationMatrix(), Position(0.,0.,-InnerCylinder_position) );
  //	Tcyl = Tcyl*ReflectZ3D();
  //	Phys = ReflectionFactory::Instance()->Place(
  //                                                Tcyl,
  //                                                "FTDInnerCylinder",
  //                                                FTDInnerCylinderLogical,
  //                                                worldLog,
  //                                                false,
  //                                                0);
//fg  registerPV( Phys );
  //=END============================ INNER CYLINDER =============================END==/
#endif
#endif
	

  
  //######################################################################################################################################################################
  

  zDiskPetalsData->widthStrip  = _dbParExReco.strip_width  ;
  zDiskPetalsData->lengthStrip = _dbParExReco.strip_length ;
  zDiskPetalsData->pitchStrip  = _dbParExReco.strip_pitch  ;
  zDiskPetalsData->angleStrip  = _dbParExReco.strip_angle  ;
  


  ftd.addExtension< ZDiskPetalsData >( zDiskPetalsData ) ;
  
  //--------------------------------------
  

  return ftd;
}

DECLARE_DETELEMENT(FTD_Simple_Staggered ,create_element)





//===================================================================================================================================================================================

//================================ PETAL BUILD FUNCTIONS ======================================/

//fg: obsolete - call the two functions directly 
// //***********************************************************************************************
// // Build the support, sensors and electronics, choosing what technology have the current disk    
// //
// // Input Parameters: 
// //                   valuesDict: map (name variable, its value)  containing some dimension-disk 
// //                               parameters, to be passed to the real building functions
// //                   mother:     Volume where will be placed the volumes are going to
// //                               build.
// //
// void DoAndPlaceDisk( Detector& theDetector,DetElement det,  SensitiveDetector sens, std::map<std::string,double> valuesDict, Volume  mother )
// {
//   petalSupport(theDetector, det, valuesDict, mother ) ; // support is placed at 0,0,0 withing the petal
//   petalSensor( theDetector, det, sens, valuesDict, mother );
//}

//***********************************************************************************************
// Build the petal  support. The support is a trapezoid made of foam.
//
// Input Parameters: 
//                   valuesDict: map (name variable, its value)  containing some dimension-disk 
//                               parameters, to be passed to the real building functions
//                   mother:     Volume the volumes built are to be placed.


void petalSupport( Detector& theDetector, DetElement ftd,  std::map<std::string,double> valuesDict, Volume  FTDPetalAirLogical )
{
  double petal_cp_supp_half_dxMin = valuesDict["petal_cp_supp_half_dxMin"];
  double petal_cp_support_dy = valuesDict["petal_cp_support_dy"];
  double _inner_radius_petal = valuesDict["_inner_radius"];

  if( _dbParDisk.sensor_is_pixel == 1) {

    Trap FTDPetalSupportSolid( _dbParDisk.petal_cp_support_thickness/2.0, //thickness
			       0.0,
			       0.0,
			       petal_cp_support_dy/2.0,  // dy
			       petal_cp_supp_half_dxMin, //dxMin 
			       _dbParDisk.petal_cp_support_dxMax/2.0, //dxMax
			       0.0,
			       petal_cp_support_dy/2.0,  // dy
			       petal_cp_supp_half_dxMin,  // dxMin
			       _dbParDisk.petal_cp_support_dxMax/2.0, //dxMax
			       0.0);
    
    Volume FTDPetalSupportLogical( _toString(  _dbParDisk.disk_number, "FTDPetalSupportLogical_%d" ), FTDPetalSupportSolid, _CarbonFiberMat );
    //printVolume( FTDPetalSupportLogical ) ;

    ftd.setVisAttributes(theDetector,  "FTDSupportVis" , FTDPetalSupportLogical ) ; 

    // Position Ta;
    // Ta.setX(0.0); 
    // Ta.setY(0.0); 
    // Ta.setZ(0.0); 

    // PhysicalVolumesPair Phys = ReflectionFactory::Instance()->Place(
    // 								    Transform3D(RotationMatrix(),Ta),
    // 								    "FTDPetalSupport",
    // 								    FTDPetalSupportLogical,
    // 								    FTDPetalAirLogical,
    // 								    false,
    // 								    0);
    //   registerPV(Phys);

    // PlacedVolume pv =
    FTDPetalAirLogical.placeVolume( FTDPetalSupportLogical , Transform3D() ) ;



  }  else {

    Trap FTDPetalSupportCPSolid( _dbParDisk.petal_cp_support_thickness/2.0,//thickness
				 0.0,
				 0.0,
				 petal_cp_support_dy/2.0,  // height
				 petal_cp_supp_half_dxMin, 
				 _dbParDisk.petal_cp_support_dxMax/2.0,
				 0.0,
				 petal_cp_support_dy/2.0,  // height
				 petal_cp_supp_half_dxMin, 
				 _dbParDisk.petal_cp_support_dxMax/2.0,
				 0.0);
    
	      
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Holes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
    //fg   Trap FTDPetalSupportHoleDownSolid = SemiPetalSolid( petal_cp_support_dy, "DOWN", false );
    //fg   Trap FTDPetalSupportHoleUpSolid   = SemiPetalSolid( petal_cp_support_dy, "UP", false );
    //fg	      
    //fg   // some transformation needed 
    //fg   const double petal_cp_holes_separation = 10.0*mm; // WARNING! HARDCODED...
    //fg   
    //fg   // RotationMatrix * idRot = new RotationMatrix;
    //fg   Position movDown( 0.0, -petal_cp_support_dy/4.0+petal_cp_holes_separation/4.0, 0.0 );
    //fg   Position movUp( 0.0, petal_cp_support_dy/4.0-petal_cp_holes_separation/4.0, 0.0);
    
    //fg----- SemiPetalSolid does not seem to work - compute hole parameters here instead:

    // the space frame width 
    double spfw = _dbParCommon.support_spaceframe_width ; 

    double dxmin = petal_cp_supp_half_dxMin*2. ;
    double dxmax = _dbParDisk.petal_cp_support_dxMax  ;

    double ybase = petal_cp_support_dy / ( dxmax / dxmin  - 1. )  ;

    // compute y_base to use intercept theorem for computing
    //  the x-widths of the wholes at various y-values
    //                        dxMax                                  _
    //                     -----------   |                          | |    
    //                     \         /   |                          | |              
    //       XY-Plane       \ dxMin /    | dy          YZ-Plane     | |    
    //                       \_____/     |                          |_|     
    //                        \   /    |                            dz
    //                         \ /     | ybase
    //                          v      |

    double petal_hole_dy =  ( petal_cp_support_dy - 3 * spfw ) / 2. ;

    double petal_hole_down_dxMin = ( ( ybase + spfw                        ) / ybase  ) * dxmin  - 2. * spfw ; 
    double petal_hole_down_dxMax = ( ( ybase + spfw    + petal_hole_dy     ) / ybase  ) * dxmin  - 2. * spfw ; 
    double petal_hole_up_dxMin   = ( ( ybase + spfw *2 + petal_hole_dy     ) / ybase  ) * dxmin  - 2. * spfw ; 
    double petal_hole_up_dxMax   = ( ( ybase + spfw *2 + petal_hole_dy * 2 ) / ybase  ) * dxmin  - 2. * spfw ; 
    
    Trap FTDPetalSupportHoleDownSolid(  2.*_dbParDisk.petal_cp_support_thickness/2.0,//thickness
    					0.0,
    					0.0,
    					petal_hole_dy/2.0,  // height
    					petal_hole_down_dxMin/2., 
    					petal_hole_down_dxMax/2., 
    					0.0,
    					petal_hole_dy/2.0,  // height
    					petal_hole_down_dxMin/2., 
    					petal_hole_down_dxMax/2., 
    					0.0) ;
    
    Trap FTDPetalSupportHoleUpSolid(   2.*_dbParDisk.petal_cp_support_thickness/2.0,//thickness
    				      0.0,
    				      0.0,
    				      petal_hole_dy/2.0,  // height
    				      petal_hole_up_dxMin/2., 
    				      petal_hole_up_dxMax/2., 
    				      0.0,
    				      petal_hole_dy/2.0,  // height
    				      petal_hole_up_dxMin/2., 
    				      petal_hole_up_dxMax/2.,
    				      0.0 ) ;
    
    Position movDown( 0.0,  -( petal_hole_dy /2.  + spfw / 2. ) , 0.0 );
    Position movUp(   0.0,   ( petal_hole_dy /2.  + spfw / 2. ) , 0.0 );


#ifdef DEBUG_VALUES
    std::cout << " *** support holes parameters : " 
    	      << " petal_hole_dy "  <<   petal_hole_dy 
    	      << " petal_hole_down_dxMin "  <<  petal_hole_down_dxMin 
    	      << " petal_hole_down_dxMax "  <<   petal_hole_down_dxMax
    	      << " petal_hole_up_dxMin "  <<   petal_hole_up_dxMin 
    	      << " petal_hole_up_dxMax "  <<   petal_hole_up_dxMax 
    	      << " y shift: " << ( petal_hole_dy /4.  + spfw / 2. )  ;
#endif    
    
    
    //fg----- END: SemiPetalSolid does not seem to work - compute hole parameters here instead:
    
    
    SubtractionSolid FTDPetalSupportSolid_Prov( FTDPetalSupportCPSolid, FTDPetalSupportHoleDownSolid, movDown ) ;
    SubtractionSolid FTDPetalSupportSolid( FTDPetalSupportSolid_Prov, FTDPetalSupportHoleUpSolid, movUp ) ;
    
    //fg: debug stuff - cutting out tubes and boxes...
    //    Tube tubehole( 0, 20*mm , 2.*_dbParDisk.petal_cp_support_thickness ); 
    //    Box tubehole( 20*mm , 20*mm , 2.*_dbParDisk.petal_cp_support_thickness ); 
    // SubtractionSolid FTDPetalSupportSolid_Prov( FTDPetalSupportCPSolid, tubehole , movDown ) ;
    // SubtractionSolid FTDPetalSupportSolid( FTDPetalSupportSolid_Prov, tubehole , movUp ) ;


    //%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Holes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%/
    
    // Petal support with two holes substracted
    Volume FTDPetalSupportLogical (_toString(  _dbParDisk.disk_number, "FTDPetalSupportLogical_%d" ), FTDPetalSupportSolid,  _CarbonFiberMat ) ;
    //printVolume( FTDPetalSupportLogical ) ;

    ftd.setVisAttributes(theDetector,  "FTDSupportVis" , FTDPetalSupportLogical ) ; 
    //ftd.setVisAttributes(theDetector,  "Green Vis" , FTDPetalSupportLogical ) ; 


      // // Placing all together ( support petal boolean  ) inside the air petal container
      // PhysicalVolumesPair Phys = ReflectionFactory::Instance()->Place(
      // 								      Transform3D(),
      // 								      "FTDPetalSupport",
      // 								      FTDPetalSupportLogical,
      // 								      FTDPetalAirLogical,
      // 								      false,
      // 								      0);
      // registerPV(Phys);

      // PlacedVolume pv =
	   FTDPetalAirLogical.placeVolume( FTDPetalSupportLogical , Transform3D() ) ;

    }
  //-END--------------------------- Central Part ----------------------------------END-/
	

#ifdef DEBUG_VALUES
  cout << "===================================================================== " << "\n" <<
    "FTDPetalSupport:\n" << 
    " Inner Radius= " << _inner_radius_petal <<  "\n" <<
    " Outer Radius= " << _outer_radius <<  "\n" <<
    " xMax = " << _dbParDisk.petal_cp_support_dxMax <<  "\n" <<
    " xMin = " << 2.0*petal_cp_supp_half_dxMin << "\n" <<
    " dy =   " << petal_cp_support_dy << "\n" <<
    " thickness =   " << _dbParDisk.petal_cp_support_thickness << "\n" <<
    // " placed at \n" << 
    // " x =   " <<  Ta.getX() << "\n" <<
    // " y =   " <<  Ta.getY() << "\n" <<
    // " z =   " <<  Ta.getZ() << "\n" <<
    endl;
#endif




  // Gear Ladder 
  _ftdparameters[gearpar::SUPPORTRINNER].push_back(_inner_radius_petal);
  _ftdparameters[gearpar::SUPPORTLENGTHMIN].push_back(2.0*petal_cp_supp_half_dxMin);
  _ftdparameters[gearpar::SUPPORTLENGTHMAX].push_back(_dbParDisk.petal_cp_support_dxMax);
  _ftdparameters[gearpar::SUPPORTWIDTH].push_back(petal_cp_support_dy);
  _ftdparameters[gearpar::SUPPORTTHICKNESS].push_back(_dbParDisk.petal_cp_support_thickness);
}

//***********************************************************************************************
// Build the petal sensitive. The sensitive volume is a trapezoid made of silicon.
//
// Input Parameters: 
//                   valuesDict: map (name variable, its value)  containing some dimension-disk 
//                               parameters, to be passed to the real building functions
//                   mother:     Volume the volumes built are to be placed.

VolVec petalSensor(  Detector& theDetector, DetElement ftd, SensitiveDetector sens, std::map<std::string,double> valuesDict, Volume  FTDPetalAirLogical ) {
  
  VolVec volV ;

  double petal_half_dxMin = valuesDict["petal_cp_supp_half_dxMin"];
  double petal_dy = valuesDict["petal_cp_support_dy"];
  double _inner_radius_petalSensor = valuesDict["_inner_radius"];
  
  Trap FTDPetalSensitiveSolid( _dbParDisk.disks_Si_thickness/2.0, //thickness
				0.0,
				0.0,
				petal_dy/2.0,  // dy
				petal_half_dxMin, //dxMin 
				_dbParDisk.petal_cp_support_dxMax/2.0, //dxMax
				0.0,
				petal_dy/2.0,  // dy
				petal_half_dxMin,  // dxMin
				_dbParDisk.petal_cp_support_dxMax/2.0, //dxMax
				0.0);
  
  // Now check 
  // FIXME: sensitive detectors
  // TRKSD_FTD01* sensitive_det = 0 ;
  // if ( _dbParDisk.sensor_is_pixel ) {
  //   sensitive_det = _theFTDSD_pixel;
  // } 
  // else {
  //   sensitive_det = _theFTDSD_strip;
  // }
  
  Volume FTDPetalSensitiveLogical (_toString( _dbParDisk.disk_number , "FTDPetalSensitiveLogical_%d" ) , FTDPetalSensitiveSolid, _SiMat ) ; 
  //printVolume( FTDPetalSensitiveLogical ) ;

  FTDPetalSensitiveLogical.setSensitiveDetector( sens ) ;

  //ftd.setVisAttributes(  theDetector, "SeeThrough" , FTDPetalSensitiveLogical );
  ftd.setVisAttributes(  theDetector, "FTDSensitiveVis" , FTDPetalSensitiveLogical );

  
  // front sensor
  Position Ta( 0. , 0. , (_dbParDisk.petal_cp_support_thickness + _dbParDisk.disks_Si_thickness)/2.0 ) ;    
  
  // PhysicalVolumesPair Phys_front = ReflectionFactory::Instance()->Place(
  // 									Transform3D(RotationMatrix(),Ta),
  // 									"FTDPetalSensitive",
  // 									FTDPetalSensitiveLogical,
  // 									FTDPetalAirLogical,
  // 									false,
  // 									1);
  
  PlacedVolume pv = FTDPetalAirLogical.placeVolume( FTDPetalSensitiveLogical, Ta)  ;
  pv.addPhysVolID("sensor", 1 )   ;
  

  volV.push_back(  std::make_pair( FTDPetalSensitiveLogical , pv )  ) ;

#ifdef DEBUG_VALUES
  cout << "===================================================================== " << "\n" <<
    "FTDPetalSensitive:\n" << 
    " Inner Radius= " << _inner_radius_petalSensor <<  "\n" <<
    " Outer Radius= " << _outer_radius <<  "\n" <<
    " xMax = " <<  _dbParDisk.petal_cp_support_dxMax <<  "\n" <<
    " xMin = " << 2.0*petal_half_dxMin << "\n" <<
    " dy =   " << petal_dy << "\n" <<
    " thickness =   " << _dbParDisk.disks_Si_thickness << "\n" <<
    " placed at\n " << 
    " x =   " <<  Ta.X() << "\n" <<
    " y =   " <<  Ta.Y() << "\n" <<
    " z =   " <<  Ta.Z() << "\n" <<
    endl;
#endif
  
  if(_dbParDisk.double_sided == 1 ) { // first two disks are single sided pixel
    
    // rear sensor
    Ta.SetZ( -(_dbParDisk.petal_cp_support_thickness + _dbParDisk.disks_Si_thickness)/2.0 );    

    // PhysicalVolumesPair Phys_rear = ReflectionFactory::Instance()->Place(
    //                                                                          Transform3D(RotationMatrix(),Ta),
    //                                                                          "FTDPetalSensitive",
    //                                                                          FTDPetalSensitiveLogical,
    //                                                                          FTDPetalAirLogical,
    //                                                                          false,
    //                                                                          2);
    
    pv = FTDPetalAirLogical.placeVolume( FTDPetalSensitiveLogical, Ta )  ;
    pv.addPhysVolID("sensor", 2 )   ;

    volV.push_back(  std::make_pair( FTDPetalSensitiveLogical , pv )  ) ;

    // registerPV(Phys_front);
    // registerPV(Phys_rear);
    
    
#ifdef DEBUG_VALUES
    cout << "===================================================================== " << "\n" <<
    "FTDPetalSensitive:\n" << 
    " Inner Radius= " << _inner_radius_petalSensor <<  "\n" <<
    " Outer Radius= " << _outer_radius <<  "\n" <<
    " xMax = " <<  _dbParDisk.petal_cp_support_dxMax <<  "\n" <<
    " xMin = " << 2.0*petal_half_dxMin << "\n" <<
    " dy =   " << petal_dy << "\n" <<
    " thickness =   " << _dbParDisk.disks_Si_thickness << "\n" <<
    " placed at\n " << 
    " x =   " <<  Ta.X() << "\n" <<
    " y =   " <<  Ta.Y() << "\n" <<
    " z =   " <<  Ta.Z() << "\n" <<
    endl;
#endif
  }
  
  //        //================================ SILICON SENSORS ==================================/
  //        // Sensors build as boxes of silicon. Assemblyblblaa
  //        // 
  //        //
  //	double pixel_si_width = 9.9 * mm; //FIXME: From DB 
  //	double pixel_si_length = 7.0 * mm; //FIXME: From DB 
  //	double pixel_si_interspace = 5.0 * um; //FIXME: From DB
  //    
  //	Box * FTDPixelSolid = new Box( "FTDPixelSensor",
  //                                      pixel_si_length/2.0,
  //                                      pixel_si_width/2.0,
  //                                      _dbParDisk.disks_Si_thickness/2.0
  //                                      );
  //    
  //	Volume  FTDPixelLogical ( FTDPixelSolid,
  //                                                            _SiMat,
  //                                                            "FTDPixelSensor",
  //                                                            0,
  //                                                            0,
  //                                                            0);
  //    
  //	FTDPixelLogical->SetVisAttributes( _VisAttSensitive );
  //	FTDPixelLogical->SetSensitiveDetector(_theFTDSD); //Sensitive
  //    
  //        // Defining two (one) rows of pixels per assembly
  //	std::vector<AssemblyVolume*> pixelsAssemblyRows;
  //	
  //	int howManyTotalRows = (int)(petal_cp_support_dy/pixel_si_width);
  //        // How many rows have in the first assembly (disk 1 = 2, disk 2 =2 );
  //        // except the first assembly, all the others have two rows per assembly
  //	int numberRows = 2;
  //	if( howManyTotalRows % 2 != 0 )
  //        {
  //		numberRows = 1;
  //        }
  //	int howManyMinPixels = (int)(2.0*petal_cp_supp_half_dxMin/pixel_si_length);
  //#ifdef DEBUG_VALUES
  //	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" <<
  //    "Total Pixel rows: " << howManyTotalRows << " ( really done " << howManyTotalRows/2.0 << " assemblies ) \n";
  //#endif
  //	
  //	Position rowTrans( 0.0, 0.0, 0.0);
  //	for(int row = 0; row < (howManyTotalRows/2 + howManyTotalRows%2); ++row)
  //        {
  //            //Instantiating the assembly
  //		pixelsAssemblyRows.push_back( new AssemblyVolume() );
  //            // Number of pixels (lowest row should have 4)
  //		int nPixel = howManyMinPixels+row;
  //            // Positioning: note that if there are a odd number of
  //            // pixels, the center pixel is centered in x=0
  //		double x_offset = 0.0;
  //		if( nPixel % 2 == 1 )
  //            {
  //			x_offset = pixel_si_length/2.0;
  //                // Placing the central pixel
  //			rowTrans.setX( 0.0 );
  //			for(int j = 0; j < numberRows; ++j)
  //                {
  //				rowTrans.setY( pixel_si_width*(1.0/2.0 + j) +pixel_si_interspace );
  //				pixelsAssemblyRows.back()->AddPlacedVolume( FTDPixelLogical, 
  //                                                           rowTrans, (RotationMatrix*)0 );
  //                }
  //            }
  //            // The others pixels except the central one, if there is
  //		for(int pixelId = 0 ; pixelId < nPixel/2; ++pixelId)
  //            {
  //			rowTrans.setX( x_offset + pixel_si_length/2.0 + pixel_si_interspace + pixelId*pixel_si_length); //pixel_si_offset
  //			for(int j = 0; j < numberRows; ++j)
  //                {
  //				rowTrans.setY( pixel_si_width*(1.0/2.0 + j) + pixel_si_interspace );
  //                    // Assembly built as two (or one) rows of pixels
  //				pixelsAssemblyRows.back()->AddPlacedVolume( FTDPixelLogical, 
  //                                                           rowTrans, (RotationMatrix*)0 );
  //				rowTrans.setX( -rowTrans.getX() ) ;
  //				pixelsAssemblyRows.back()->AddPlacedVolume( FTDPixelLogical, 
  //                                                           rowTrans, (RotationMatrix*)0 );
  //                }
  //            }
  //            // All the others assemblies have two rows
  //		numberRows = 2;
  //        }
  //        //Placing the assemblies inside the air petal, begining from the bottom	
  //        //	double dz = _dbParDisk.petal_cp_support_thickness/2.0 + _dbParDisk.kapton_petal_thickness + _dbParDisk.disks_Si_thickness/2.0; 
  //	double dz = _dbParDisk.disks_Si_thickness /2.0; 
  //    
  //	double dx = 0.0;
  //	double dy = -petal_cp_support_dy/2.0;
  //	Position inMotherTrans(dx, dy, dz); 
  //	if( howManyTotalRows % 2 != 0 )
  //        {	
  //            numberRows = 1;
  //        }
  //	for( std::vector<AssemblyVolume*>::iterator assemblyRow = pixelsAssemblyRows.begin();
  //        assemblyRow != pixelsAssemblyRows.end(); ++assemblyRow )
  //        {
  //#ifdef DEBUG_VALUES
  //		static int i = 0;
  //        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" <<
  //        " Placement pixels sensors: \n" 
  //        << "   Row " << i << ": dx=" << dx << ", dy=" << dy << ", dz=" << dz 
  //        << " -- # pixel: " << (*assemblyRow)->GetInstanceCount() << "\n";
  //        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" << endl;
  //		i++;
  //#endif
  //        (*assemblyRow)->MakeImprint( FTDPetalSupportAirLogical, inMotherTrans, (RotationMatrix*)0 );
  //            // Ready to put the next one
  //		dy += 2.0*pixel_si_width;
  //		inMotherTrans.setY( dy );
  //		if( numberRows == 1 )
  //            {
  //			dy -= pixel_si_width;
  //			inMotherTrans.setY( dy );
  //			numberRows = 2;
  //            }
  //        }
  //        // Gear
  //	int howManyPixelsUp = (int)(_dbParDisk.petal_cp_support_dxMax/pixel_si_length);
	
  _ftdparameters[gearpar::SENSITIVERINNER].push_back(_inner_radius_petalSensor);
  _ftdparameters[gearpar::SENSITIVELENGTHMIN].push_back(2.0*petal_half_dxMin);
  _ftdparameters[gearpar::SENSITIVELENGTHMAX].push_back(_dbParDisk.petal_cp_support_dxMax);
  _ftdparameters[gearpar::SENSITIVEWIDTH].push_back(petal_dy);
  _ftdparameters[gearpar::SENSITIVETHICKNESS].push_back(_dbParDisk.disks_Si_thickness);
  
  
  return volV ;
}



//=END============================ PETAL BUILD FUNCTIONS ==================================END=/


//============================= PETAL DIMENSION FUNCTIONS =====================================/
// Functions to extract relative dimensions of the petals. The mechanical design for the petals
// is made taking as reference one disk (disk 1 for pixels, disk 4 for strips), so the petal
// is parametrized using variables given for the database (or calculates from them). All the
// dimensions of the petal are extracted having the inner radius and outer radius of the disk,
// the top length of the petal and the angle defined by the petal using trigonometry.
//
//   (*)outer radius
//                    dxMax/2                     dxMax
//	           ============               =============
//	          |          //               \           /
//      dy        |  (*)<--/ /                 \  PETAL  /
//	          |      /  /                   \       /
//	          -====/===/--> dxMin/2          =======   dxMin
//	          |  /    /   
//  inner radius  |/     /
//                -·····/
//	          |    /             dy = _outer_radius*cos( arcsin(dxMax/(2*_outer_radius)) ) - _inner_radius
//	          |   /              ( a = dxMax/(2*tag(theta)) - _inner_radius - dy )
//	 a        |  /               dxMin = 2*tan(theta)*( _inner_radius + a )
//	          |-/ theta
//	          |/
//    
//
// So, changing the inner radius it'll change the dy and dxMin as well, providing the diferents widths
// of the petal, needed to build the sensors, holes, etc...

//------------------------------------------------------------------------------------------
// Get the dy of the petal which corresponds to a given radius
//
// Input Parameters:   inner radius
// Output Parameters:  dy
  double Getdy(const double & innerRadius ){
    
    return _outer_radius*cos( asin(_dbParDisk.petal_cp_support_dxMax/(2.0*_outer_radius)) ) - innerRadius;
  }
//------------------------------------------------------------------------------------------
// Get the dxMin of the petal which corresponds to a given radius
//
// Input Parameters:   inner radius
// Output Parameters:  dxMin
 double Getdx( const double & innerRadius ) {
   
   double a = _dbParDisk.petal_cp_support_dxMax/(2.0*tan(_dbParCommon.petal_half_angle_support)) - innerRadius - 
     Getdy( innerRadius );
   
   return 2.0*(innerRadius + a)*tan(_dbParCommon.petal_half_angle_support);
 }
 
//------------------------------------------------------------------------------------------
// Get the dimensions of the up and down holes or the silicon up or down sensors:
// 
// Input parameters:      
//                   petal_cp_support_dy: height of the air petal container
//                   whereItgoes:         "UP" or "DOWN", where is placed
//                   isSilicon:           True of False, define if is the sensor or the holes
//
// Output:         
//                  std::vector<double>* =  [ xMin_half, xMax_half, dy_half, thickness_half ]
std::vector<double> GetPetalDimensions( const double& petal_cp_support_dy, const std::string & whereItgoes, const bool isSilicon )
{
  const double theta = _dbParCommon.petal_half_angle_support;
	
  double central_separation_y = 10.0*mm/2.0; //HARDCODED _dbParDisk.petal_cp_holes_separation/2.0;
  double x_dim = 6.0*mm/cos(theta);          //HARDCODED _dbParDisk.petal_cp_holes_width_support/cos(theta);
  double y_dimension = 10.0*mm;        // HARDCODED _dbParDisk.petal_cp_holes_separation;
  double half_thickness = _dbParDisk.petal_cp_support_thickness/2.0;
    
  //Silicon detector or Hole?
  if(isSilicon)
    {
      central_separation_y = 0.0;
      const double padUp_Si_dxMax = 118.46*mm; 
      x_dim = (_dbParDisk.petal_cp_support_dxMax-padUp_Si_dxMax)/2.0; // HARDCODED
      //x_dim = (_dbParDisk.petal_cp_support_dxMax-_dbParDisk.padUp_Si_dxMax)/2.0;
      y_dimension = y_dimension/2.0;
      half_thickness = _dbParDisk.disks_Si_thickness/2.0;
    }
	
  double pseudo_radius_up;
  double pseudo_radius_down;
  // Up or down?
  if( whereItgoes == "UP" )
    {
        
      pseudo_radius_up   = _inner_radius+petal_cp_support_dy - y_dimension;
      pseudo_radius_down = _inner_radius+petal_cp_support_dy*_dbParCommon.petal_y_ratio + central_separation_y ;
    }
  else if( whereItgoes == "DOWN" )
    {
      pseudo_radius_up = _inner_radius + petal_cp_support_dy*_dbParCommon.petal_y_ratio - central_separation_y;
      pseudo_radius_down = _inner_radius + y_dimension;
    }
  else
    {
      cout << "FTD_Simple_Staggered: Internal Error: The function SemiPetalSolid is not well called, the " <<
        "4th argument must be \"UP\" or \"DOWN\".\n Check the code!!" << endl;
      exit(-1);
    }
  const double xMin_half = (Getdx( pseudo_radius_down ) - x_dim)/2.0;
  const double xMax_half = (Getdx( pseudo_radius_up )- x_dim)/2.0;
  const double dy = pseudo_radius_up - pseudo_radius_down;
#ifdef DEBUG_VALUES
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" <<
    "    IsSilicon?: " << isSilicon <<  
    " " << whereItgoes <<
    " xMin=" << 2.*xMin_half << 
    " xMax=" << 2.0*xMax_half << 
    " dy= " << dy << std::endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" << endl;
#endif
  std::vector<double>  dimensions ;
  dimensions.push_back( xMin_half );
  dimensions.push_back( xMax_half );
  dimensions.push_back( dy/2.0 );
  dimensions.push_back( half_thickness );
    
  return dimensions;
}

//------------------------------------------------------------------------------------------
// Construction of the up and down holes or the silicon up or down sensors
// 
// Input parameters:      
//                   petal_cp_support_dy: height of the air petal container
//                   whereItgoes:         "UP" or "DOWN", where is placed
//                   isSilicon:           True of False, define if is the sensor or the holes

//-------------------------------------------------------------------------------------------------------------------
Trap SemiPetalSolid( const double& petal_cp_support_dy, const std::string& whereItgoes, const bool isSilicon ){
  // Get dimensions
  const std::vector<double>&  dimensions = GetPetalDimensions( petal_cp_support_dy, whereItgoes, isSilicon );
  double xMin_half      = dimensions[0];
  double xMax_half      = dimensions[1];
  double dy_half        = dimensions[2];
  double half_thickness = dimensions[3];
  
  return Trap( half_thickness,//thickness
	       0.0,
	       0.0,
	       dy_half,  // height
	       xMin_half, 
	       xMax_half,
	       0.0,
	       dy_half,  // height
	       xMin_half, 
	       xMax_half,
	       0.0 ) ;
  
  //  return FTDSemiPetalSolid ;
 }
 //=END========================= PETAL DIMENSION FUNCTIONS =================================END=/

//fixme: registering sensitive ...
// //*********************************************************************************************
// // Register two ftd sensitive detectors, one for the pixel disks and one for the strip disks
// void RegisterSDs( Database * db )
// {
//   // Getting parameters disk-specific
//   db->exec("select * from disks;");
//   db->getTuple();
  
//   double minDiskThickness (MAXFLOAT);
//   do
//     	{
//     _dbParDisk.disks_Si_thickness = db->fetchDouble("disk_si_thickness")  ;
//     if(minDiskThickness>_dbParDisk.disks_Si_thickness)
//         {
// 			minDiskThickness = _dbParDisk.disks_Si_thickness;
//         }
//      	} while(db->getTuple()!=NULL);
  
//   this->_theFTDSD_pixel =  new TRKSD_FTD01( "FTD_PIXEL", minDiskThickness * 340 * keV/mm * 0.2,true);
//   RegisterSensitiveDetector(_theFTDSD_pixel);

//   this->_theFTDSD_strip =  new TRKSD_FTD01( "FTD_STRIP", minDiskThickness * 340 * keV/mm * 0.2,true);
//   RegisterSensitiveDetector(_theFTDSD_strip);

// }


// void registerPV(const PhysicalVolumesPair & pvPair )
// {
// 	if( pvPair.first != 0 )
//       {
// 		_registerPV.push_back( pvPair.first );
//       }
// 	if( pvPair.second != 0 )
//       {
// 		_registerPV.push_back( pvPair.second );
//       }
// }

// //================================ GEAR STUFF FUNCTIONS ====================================/
#ifdef DD4HEP_WITH_GEAR
//fixme: seems to be never called in Mokka class (maybe outside ? )

void GearSetup()
{	
  //--- Added carbon fiber.  
  //    TODO: It is needed some other changes??  
  //    October, 2010, J.Duarte
  double Si_RadLen, Si_dEdx;
  double Kapton_RadLen, Kapton_dEdx;
  double CarbonFiber_RadLen, CarbonFiber_dEdx;
  double Cu_RadLen, Cu_dEdx;
	
  Si_RadLen = _SiMat->GetRadlen();
  Kapton_RadLen = _KaptonMat->GetRadlen();
  Cu_RadLen = _CuMat->GetRadlen();
  CarbonFiber_RadLen = _CarbonFiberMat->GetRadlen();
	
  //... Looping over bins in the DEDX table to obtain the mip DEDX 
  //... From energy 0.0001MeV to 1000MeV in steps of 10 (See GetdEdx function)
  Si_dEdx=GetdEdx( _SiMat );
  Kapton_dEdx=GetdEdx(_KaptonMat);
  CarbonFiber_dEdx = GetdEdx(_CarbonFiberMat);
  Cu_dEdx=GetdEdx( _CuMat );
  
  // Parameters for FTD
  gear::GearMgr* gearMgr = MokkaGear::getMgr() ;
  
  gear::FTDParametersImpl* ftdParam = new gear::FTDParametersImpl();
  
  
  // Write gearParameters to GearMgr
  for(unsigned int layer = 0; layer < _ftdparameters[gearpar::NPETALS].size(); layer++)
    {
      // Extract all the param
      int	nPetals       =   (int)_ftdparameters[gearpar::NPETALS].at(layer);
      int	nSensors      =   (int) _ftdparameters[gearpar::NSENSORS].at(layer);
      bool	isDoubleSided =   (bool)_ftdparameters[gearpar::ISDOUBLESIDED].at(layer);
      int	sensorType    =   (int)_ftdparameters[gearpar::SENSORTYPE].at(layer);
      double	phalfangle    =   _ftdparameters[gearpar::HALFANGLEPETAL].at(layer);
      double	phi0	      =   _ftdparameters[gearpar::PHI0].at(layer);
      // Correct the sign: axis of the trapezoids built inside a ref. with Z --> -Z
      double	signoffset    = - _ftdparameters[gearpar::PETAL0SIGNOFFSET].at(layer);
      double	alpha	      =   _ftdparameters[gearpar::ALPHA].at(layer);
      double	zposition     =   _ftdparameters[gearpar::ZPOSITION].at(layer);
      double	zoffset	      =   _ftdparameters[gearpar::ZOFFSET].at(layer);
      double	suprtRin      =   _ftdparameters[gearpar::SUPPORTRINNER].at(layer);
      double	suprtThic     =   _ftdparameters[gearpar::SUPPORTTHICKNESS].at(layer);
      double	suprtLMin     =   _ftdparameters[gearpar::SUPPORTLENGTHMIN].at(layer);
      double	suprtLMax     =   _ftdparameters[gearpar::SUPPORTLENGTHMAX].at(layer);
      double	suprtW	      =   _ftdparameters[gearpar::SUPPORTWIDTH].at(layer);
      //double suprtRL	      =   _ftdparameters[gearpar::SUPPORTRADLENGTH].at(layer); FIXME
      double	suprtRL	      =    Si_RadLen;
      double	sensitRin     =   _ftdparameters[gearpar::SENSITIVERINNER].at(layer);
      double	sensitThic    =   _ftdparameters[gearpar::SENSITIVETHICKNESS].at(layer);
      double	sensitLMin    =   _ftdparameters[gearpar::SENSITIVELENGTHMIN].at(layer);
      double	sensitLMax    =   _ftdparameters[gearpar::SENSITIVELENGTHMAX].at(layer);
      double	sensitW	      =   _ftdparameters[gearpar::SENSITIVEWIDTH].at(layer);
      //double sensitRL	      =   _ftdparameters[gearpar::SENSITIVERADLENGTH].at(layer); //FIXME
      double	sensitRL      = Si_RadLen;
    
      ftdParam->addLayer( nPetals, nSensors, isDoubleSided, sensorType, phalfangle, phi0, alpha,zposition, zoffset, signoffset,
			  suprtRin, suprtThic, 
			  suprtLMin, suprtLMax,
			  suprtW, suprtRL,
			  sensitRin, sensitThic,
			  sensitLMin, sensitLMax,
			  sensitW, sensitRL ) ;
    }
  
  
  // Add the extended_reconstruction_parameters
  
  ftdParam->setDoubleVal("strip_width", _dbParExReco.strip_width );
  ftdParam->setDoubleVal("strip_length",_dbParExReco.strip_length);
  ftdParam->setDoubleVal("strip_pitch", _dbParExReco.strip_pitch );
  ftdParam->setDoubleVal("strip_angle", _dbParExReco.strip_angle );
  
  
  gearMgr->setFTDParameters(ftdParam);
}
#endif

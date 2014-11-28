/*
 * Helper classes for FTD Self-Scaling Driver from Mokka
 * copied from original FTD_Simple_Staggered.hh. 
 * 
 * @author F.Gaede, DESY 
 * @date  May 2014
 * @version $Id$
 */
#ifndef FTD_Simple_Staggered_h
#define FTD_Simple_Staggered_h 1

// class LogicalVolume;
// class Database;
// class TRKSD_FTD01;
// class VisAttributes;

// #include "VSubDetectorDriver.hh"
// #include "Material.hh"
// #include "Trap.hh"
// #include "VisAttributes.hh"
// #include "ReflectionFactory.hh"

#include <map>

class gearpar {

public:
  static const short NPETALS             = 1; 
  static const short ZOFFSET             = 2;
  static const short ALPHA               = 3;
  static const short PHI0                = 4;
  static const short PETAL0SIGNOFFSET    = 5;
  static const short HALFANGLEPETAL      = 6;
  static const short ZPOSITION           = 7;
  static const short SUPPORTRINNER       = 8;
  static const short SUPPORTLENGTHMIN    = 9;
  static const short SUPPORTLENGTHMAX    =10; 
  static const short SUPPORTWIDTH        =11; 
  static const short SUPPORTTHICKNESS    =12;
  static const short SUPPORTRANDLENGTH   =13; 
  static const short SENSITIVERINNER     =14; 
  static const short SENSITIVELENGTHMIN  =15; 
  static const short SENSITIVELENGTHMAX  =16;
  static const short SENSITIVEWIDTH      =17; 
  static const short SENSIVIVERANLENGHT  =18;
  static const short SENSITIVETHICKNESS  =19; 
  static const short SENSORTYPE          =20; 
  static const short NSENSORS            =21;
  static const short ISDOUBLESIDED       =22;
    
};

/** Structures to store common parameters (database, etc.. )
 * Possible improvement: map<string,double>, then no need
 * of struct, string = parameter name in database
 */
struct glEnviron
{
  double TPC_Ecal_Hcal_barrel_halfZ;
  double Ecal_endcap_zmin;
  double TPC_inner_radius;
	
  double SIT1_Half_Length_Z;
  double SIT2_Half_Length_Z;
  double SIT1_Radius;
  double SIT2_Radius;
  double VXD_layer3_maxZ;
	
  double zEnd_IPOuterTube;
  double rEnd_IPOuterTube;
  double zEnd_IPOuterBulge;
  double rEnd_IPOuterBulge;
	
  double beamTubeTangent;
};

/// Helper struct
struct dbInfoCommon
{
  double ftd1_vtx3_distance_z;
  double ftd7_ecal_distance_z;
  double ftd1_sit1_radial_diff;
  double ftd2_sit1_radial_diff;
  double ftd3_sit2_radial_diff;
  double ftd4to7_tpc_radial_gap;
    
  double beamTubeClearance;
  double cables_thickness;
  double cable_shield_thickness;
    
  double outer_cylinder_total_thickness;
  double inner_cylinder_total_thickness;
    
  // Petal specific
  double petal_half_angle_support;  // theta
  double petal_y_ratio;

  //fg: additional parameter:
  double support_spaceframe_width ;
};

/// helper struct
struct dbExtended_reconstruction_parameters
{

  double strip_width;    
  double strip_length;   
  double strip_pitch;    
  double strip_angle;   
  
};

struct dbInfoDisk
{
  int disk_number;
  
  int sensor_is_pixel;
  int double_sided;
  
  double disks_Si_thickness;  
  // Support Cylinders
  double ZStartOuterCylinder;
  double ZStopOuterCylinder;
  double ZStartInnerCylinder;
  double ZStopInnerCylinder;
  // Petal	
  double petal_cp_support_thickness;
  double petal_cp_support_dxMax;
  double petal_support_zoffset;
  
};


// class FTD_Simple_Staggered : public VSubDetectorDriver
// {
// public:
//     FTD_Simple_Staggered(void) : VSubDetectorDriver("FTD_Simple_Staggered","ftd")  {} 
//     ~FTD_Simple_Staggered(void) {};
    
//     bool ContextualConstruct(const CGAGeometryEnvironment &env, LogicalVolume *worldLog);
    
// #ifdef MOKKA_GEAR    
//     void GearSetup();
//     double GetdEdx( const Material* material );
// #endif
    
// private:
// 		// register PV
//     void registerPV( const PhysicalVolumesPair & pvPair );
    
// 		// Setters for parameters -----
//     void SetEnvironPar( const CGAGeometryEnvironment & env );
//     void SetdbParCommon( Database * db );
//     void SetParDisk( Database * db );
//     // ----------------------------
//     void RegisterSDs( Database * db ); // registers two sensitive detectors one for pixel disks and one for strips
    
// 		// Construction of parametrized values for the petals (dy, dxMin,..)
//     double Getdy( const double & the_inner_radius );
//     double Getdx( const double & the_inner_radius );
// 		// -----------------------------------------------------------------
    
// 		// Petal support construction ------------------------------------
//     void DoAndPlaceDisk( std::map<std::string,double> valuesDict, LogicalVolume * mother );
    
//     void petalSupport( std::map<std::string,double>  valuesDict, LogicalVolume * mother ); 
//     void petalSensor( std::map<std::string,double>  valuesDict, LogicalVolume * mother ); 
    

//     std::vector<double> * GetPetalDimensions( const double& petal_cp_support_dy, 
//                                                const String& whereItgoes, const bool isSilicon );

//     Trap * SemiPetalSolid(const double& petal_cp_support_dy, const String& whereItgoes,
//                             const bool isSilicon ); 
// 		// -----------------------------------------------------------------
    
// 		// DB Parameters
//     dbInfoCommon _dbParCommon;
//     dbInfoDisk _dbParDisk;
// 	  dbExtended_reconstruction_parameters _dbParExReco;

//     // Global environment variables
//     glEnviron _glEnv;
  
  
//     Material* _SiMat;
//     Material* _KaptonMat;
//     Material* _CuMat;
//     Material* _AirMat;
//     Material* _CarbonFiberMat;
    
//     VisAttributes * _VisAttSupport;
//     VisAttributes * _VisAttSensitive;
//     VisAttributes * _VisAttHolePetal;
    
// 		// Disk positioning parameters
//     double _z_position;
//     double _zEnd;
//     double _inner_radius;
//     double _outer_radius;
//     double _beamTubeRadius;
    
//     TRKSD_FTD01* _theFTDSD_pixel;
//     TRKSD_FTD01* _theFTDSD_strip;
    
// 		// Storing all the gear parameters
//     std::map<int,std::vector<double> > _ftdparameters;
    
// 		// Registring the Physical Volumes to be deleted
//     std::vector<VPhysicalVolume*> _registerPV;
// };

#endif



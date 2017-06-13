#ifndef SServices00_v01_h
#define SServices00_v01_h 1
//
//*******************************************************
//*                                                     *
//*                      Mokka                          * 
//*   - the detailed geant4 simulation for ILC   -      *
//*                                                     *
//* For more information about Mokka, visit the         *
//*                                                     *
//*  Mokka.in2p3.fr  Mokka home page.                   *
//*                                                     *
//*******************************************************
//
// $Id: SServices00.hh,v 1.0 2009/06/04 11:40:42 musat Exp $
// $Name: mokka-07-00 $
//

// #define VERBOSE 1

#define MAX_TPC_RINGS 8

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"
//#include "XMLHandlerDB.h"
#include <cmath>
#include <map>
#include <string>

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;
using namespace DD4hep::DDRec;
using namespace DDSurfaces;

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

class SServices00_v01
{
 public:
  SServices00_v01();
  
  ~SServices00_v01();

   bool
   BuildServices(PlacedVolume &pVol,
		 Assembly &envelope,
		 LCDD &lcdd);
private:

   bool
   BuildTPCEndplateServices(PlacedVolume &pVol,
			    Assembly &envelope,
			    LCDD &lcdd);
   
   bool
   BuildEcalBarrelServices(PlacedVolume &pVol,
			   Assembly &envelope,
			   LCDD &lcdd);
   
   bool
   BuildEcalBarrel_EndCapServices(PlacedVolume &pVol,
				  Assembly &envelope,
				  LCDD &lcdd);

   bool
   BuildHcalBarrel_EndCapServices(PlacedVolume &pVol,
				  Assembly &envelope,
				  LCDD &lcdd);

   bool
   FillEcalBarrelServicesContainer(PlacedVolume &pVol,
				   Volume &pContainerLogical,
				   LCDD &lcdd);
   
   bool
   FillHcalServicesModuleWithInnerServices(PlacedVolume &pVol,
					   Volume &ModuleLogicalZMinus,
					   Volume &ModuleLogicalZPlus,
					   LCDD &lcdd);

   
   bool
   PlaceHcalInnerServicesLayer(PlacedVolume &pVol,Volume &,
			       Material, double, Position &);
   
   bool
   FillHcalServicesModuleWithHcalElectronicsInterface(PlacedVolume &pVol,
						      Volume &ModuleLogicalZMinus,
						      Volume &ModuleLogicalZPlus,
						      LCDD &lcdd);

   

   bool
   FillHcalElectronicsInterfaceLayer(PlacedVolume &pVol,
				     Volume &ModuleLogicalZMinus,
				     Volume &ModuleLogicalZPlus,
				     double, double,
				     LCDD &lcdd);
   
   bool
   PlaceHcalElectronicsInterfaceComponent(PlacedVolume &pVol,
					  Volume &ModuleLogicalZMinus,
					  Volume &ModuleLogicalZPlus,
					  Material layerMaterial,
					  double, double, double);

   Solid CutLayer(Solid &, double);
      
   void
   BuildSitCables(PlacedVolume &pVol,
		  Assembly &envelope,
		  LCDD &lcdd);

   string FTD_db_name;
   double TPC_inner_radius;
   double Sit_cables_cylinder_thickness;
   double TUBE_IPOuterBulge_end_z;
   double TUBE_IPOuterBulge_end_radius;
   double Sit_cables_disk_thickness;
   double SIT1_Radius, SIT2_Radius;
   double FTD2_cone_thickness, FTD3_cone_thickness;

   double TPC_Ecal_Hcal_barrel_halfZ;
   double Ecal_cables_gap;
   double TPC_outer_radius;
   double InnerServicesWidth;
   double RailHeight;
   double Ecal_outer_radius;
   double module_thickness;
   double bottom_dim_x;
   double top_dim_x;
   double Hcal_total_dim_y;
   double Hcal_y_dim2_for_x;
   double Hcal_bottom_dim_x;
   double Hcal_midle_dim_x;
   double Hcal_top_dim_x;

   const double SurfaceTolerance = 0.0001;
   /*
   G4VisAttributes * VisAttAir;
   G4VisAttributes * VisAttPE;
   G4VisAttributes * VisAttCu;
   G4VisAttributes * VisAttStainlessSteel;

#ifdef MOKKA_GEAR
  // MokkaGear

  struct helpParameters{
    double innerRadius;
    double outerRadius;
    double zMax;
    double phi0;
    std::vector<double> layerPos;
    std::vector<double> radiThickness;
    int count;
    double leastZ;
    double mostZ;
  };
  
  helpParameters helpBarrel;
  helpParameters helpEndcap;
  helpParameters helpPlug;
#endif
   */  
};

#endif


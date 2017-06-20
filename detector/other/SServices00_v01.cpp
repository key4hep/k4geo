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
// $Id: SServices00.cc,v 1.1 2009/06/04 11:40:42 musat Exp $
// $Name: mokka-07-00 $
//
// 
//
// SServices00_v01.cpp
//

#include "SServices00_v01.h"

using dd4hep::Assembly;
using dd4hep::Box;
using dd4hep::Detector;
using dd4hep::IntersectionSolid;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Polycone;
using dd4hep::Position;
using dd4hep::RotationZYX;
using dd4hep::Solid;
using dd4hep::SubtractionSolid;
using dd4hep::Torus;
using dd4hep::Transform3D;
using dd4hep::Trapezoid;
using dd4hep::Tube;
using dd4hep::UnionSolid;
using dd4hep::Volume;


SServices00_v01::SServices00_v01() 
{
  std::cout<<"SServices00_v01 is being created "<<std::endl;
}  

SServices00_v01::~SServices00_v01() 
{
}  

bool SServices00_v01::BuildServices(PlacedVolume &pVol,
				    Assembly &envelope,
				    Detector &theDetector) 
{

  std::cout << "\nBuilding SServices00"<< std::endl;

  TPC_inner_radius =
    theDetector.constant<double>("TPC_inner_radius");

  Sit_cables_cylinder_thickness =
    theDetector.constant<double>("SIT12_cable_thickness");

  TUBE_IPOuterBulge_end_z  =
    theDetector.constant<double>("TUBE_IPOuterBulge_end_z");

  TUBE_IPOuterBulge_end_radius  =
    theDetector.constant<double>("TUBE_IPOuterBulge_end_radius");

  Sit_cables_disk_thickness =((1+TUBE_IPOuterBulge_end_radius/TPC_inner_radius)*0.5*
			      (theDetector.constant<double>("SServices_FTD7_cables_thickness")));

  SIT1_Radius = theDetector.constant<double>("SIT1_Radius");
  SIT2_Radius = theDetector.constant<double>("SIT2_Radius");

  FTD2_cone_thickness =
    theDetector.constant<double>("SServices_FTD2_cone_thickness");

  FTD3_cone_thickness =
    theDetector.constant<double>("SServices_FTD3_cone_thickness");

  TPC_Ecal_Hcal_barrel_halfZ = 
    theDetector.constant<double>("TPC_Ecal_Hcal_barrel_halfZ");

  TPC_outer_radius = 
    theDetector.constant<double>("TPC_outer_radius");

  Ecal_cables_gap =
    theDetector.constant<double>("Ecal_cables_gap");

  InnerServicesWidth = 
    theDetector.constant<double>("HcalServicesModule_InnerServicesWidth");

  RailHeight = 
    theDetector.constant<double>("EcalBarrelServices_RailHeight");

  if(RailHeight > theDetector.constant<double>("Hcal_Ecal_gap"))
	throw std::runtime_error("EcalBarrel services: RailHeight exceeds Hcal_Ecal_gap");

  double Ecal_inner_radius = TPC_outer_radius +
    theDetector.constant<double>("Ecal_Tpc_gap");

  Ecal_outer_radius =
    theDetector.constant<double>("Ecal_outer_radius");

  module_thickness = Ecal_outer_radius - Ecal_inner_radius;

  // module barrel key parameters
  bottom_dim_x = 2. * tan(M_PI/8.) * Ecal_inner_radius +
    module_thickness/sin(M_PI/4.);

  top_dim_x = bottom_dim_x - 2 * module_thickness;

#ifdef VERBOSE
  G4cout << "top_dim_x = " << top_dim_x << G4endl;
#endif

  if (!BuildTPCEndplateServices(pVol,envelope,theDetector))  return false;

  if (!BuildEcalBarrelServices(pVol,envelope,theDetector))  return false;

  if (!BuildEcalBarrel_EndCapServices(pVol,envelope,theDetector))  
  		return false;

  if (!BuildHcalBarrel_EndCapServices(pVol,envelope,theDetector))  
  		return false;

  BuildSitCables(pVol,envelope,theDetector);

  return true;

}

void SServices00_v01::BuildSitCables(PlacedVolume &pVol,
				     Assembly &envelope,
				     Detector &theDetector) 
{

//First place the Cylinder

  double ftd4to7_tpc_radial_gap  = 
        theDetector.constant<double>("ftd4to7_tpc_radial_gap");

  if(ftd4to7_tpc_radial_gap < Sit_cables_cylinder_thickness)
	throw std::runtime_error("SServices_02_v00: the ftd-tpc radial gap is less than Sit cable thickness");

  double SitTube_inner_radius = TPC_inner_radius -
	ftd4to7_tpc_radial_gap/2.;

  double z_start_3 = TPC_Ecal_Hcal_barrel_halfZ * 
	theDetector.constant<double>("z_position_ReltoTPCLength");

  double petalairthickness_half = 0.5 * ( theDetector.constant<double>("petal_cp_support_thickness")
                                             + 2.0*(theDetector.constant<double>("disk_si_thickness")));

  double max_half_thickness_disk_3 = theDetector.constant<double>("petal_support_zoffset") + petalairthickness_half ;

  double z_half_len = (TPC_Ecal_Hcal_barrel_halfZ -
	z_start_3) / 2.;

  Tube SitTubeSolid (SitTube_inner_radius,
		     SitTube_inner_radius + Sit_cables_cylinder_thickness,
		     z_half_len,
		     0., 2 * M_PI);

  Volume SitTubeLog("SitTubeLog",SitTubeSolid,
	       theDetector.material("aluminium"));
  
  Position Cylinder_pos1(0, 0, z_start_3 + z_half_len);
  pVol = envelope.placeVolume(SitTubeLog,Cylinder_pos1);

  Position Cylinder_pos2(0, 0, -z_start_3 - z_half_len);
  pVol = envelope.placeVolume(SitTubeLog,Cylinder_pos2);


//Then place the cone

  double z_start_2 = TPC_Ecal_Hcal_barrel_halfZ * 
	theDetector.constant<double>("z_position_ReltoTPCLength");

  petalairthickness_half = 0.5 * ( theDetector.constant<double>("petal_cp_support_thickness")
                                             + 2.0*(theDetector.constant<double>("disk_si_thickness")));

  double max_half_thickness_disk_2 = theDetector.constant<double>("petal_support_zoffset") + petalairthickness_half ;

  double FTD2_outer_radius = SIT1_Radius + 
	theDetector.constant<double>("ftd2_sit1_radial_diff");

  double FTD3_outer_radius = SIT2_Radius +
	theDetector.constant<double>("ftd3_sit2_radial_diff");

  double zPlane[2];
  zPlane[0] = z_start_2 + max_half_thickness_disk_2 + SurfaceTolerance;

  zPlane[1] = z_start_3 - max_half_thickness_disk_3 - SurfaceTolerance;

  double rInner[2];
  rInner[0] = FTD2_outer_radius + SurfaceTolerance;
  rInner[1] = FTD3_outer_radius + SurfaceTolerance;

  double rOuter[2];
  rOuter[0] = FTD2_outer_radius + FTD2_cone_thickness;
  rOuter[1] = FTD3_outer_radius + FTD3_cone_thickness;

  std::vector<double> rmin;
  std::vector<double> rmax;
  std::vector<double> z;
  for(int i=0;i<=1;i++) {
    rmin.push_back(rInner[i]);
    rmax.push_back(rOuter[i]);
    z.push_back(zPlane[i]);
  }

  Polycone SitConeSolid (0., 2.0 * M_PI, rmin, rmax, z);
		
  Volume SitConeLog("SitConeLog",SitConeSolid,
		    theDetector.material("aluminium"));
  
  Position cone_pos(0, 0, 0);
  pVol = envelope.placeVolume(SitConeLog,cone_pos);

  RotationZYX rot(0,M_PI,0);
  Transform3D tran3D(rot,cone_pos);
  pVol = envelope.placeVolume(SitConeLog,tran3D);


  //Then place the disk

  Tube SitDiskSolid(TUBE_IPOuterBulge_end_radius + SurfaceTolerance,
		    TPC_inner_radius,
		    Sit_cables_disk_thickness/2. - SurfaceTolerance,
		    0., 2 * M_PI);

  Volume SitDiskLog("SitDiskLog",SitDiskSolid,
		    theDetector.material("aluminium"));
  
  Position disk_pos1(0, 0, TUBE_IPOuterBulge_end_z + Sit_cables_disk_thickness/2.);
  pVol = envelope.placeVolume(SitDiskLog,disk_pos1);

  Position disk_pos2(0, 0, -TUBE_IPOuterBulge_end_z - Sit_cables_disk_thickness/2.);
  pVol = envelope.placeVolume(SitDiskLog,disk_pos2);

  // Print out some parameters
  std::cout <<"\n   - SIT12_cable_thickness = "
	    <<Sit_cables_cylinder_thickness
	    <<"\n   - SIT1_Radius = "
	    <<SIT1_Radius
	    <<"\n   - SIT2_Radius = "
	    <<SIT2_Radius
	    <<"\n   - SServices_FTD2_cone_thickness = "
	    <<FTD2_cone_thickness
	    <<"\n   - SServices_FTD3_cone_thickness = "
	    <<FTD3_cone_thickness
	    <<"\n   - SServices_FTD7_cables_thickness = "
	    <<theDetector.constant<double>("SServices_FTD7_cables_thickness")
	    <<std::endl;

}

bool
SServices00_v01::
BuildTPCEndplateServices(PlacedVolume &pVol,
			 Assembly &envelope,
			 Detector &theDetector)
{
   double tpcEndplateServices_R[MAX_TPC_RINGS], 
	tpcEndplateServices_r[MAX_TPC_RINGS];

   tpcEndplateServices_R[0] = 
	theDetector.constant<double>("tpcEndplateServicesRing1_R");
   tpcEndplateServices_r[0] = 
	theDetector.constant<double>("tpcEndplateServicesRing1_ro");

   tpcEndplateServices_R[1] = 
	theDetector.constant<double>("tpcEndplateServicesRing2_R");
   tpcEndplateServices_r[1] = 
	theDetector.constant<double>("tpcEndplateServicesRing2_ro");

   tpcEndplateServices_R[2] = 
	theDetector.constant<double>("tpcEndplateServicesRing3_R");
   tpcEndplateServices_r[2] = 
	theDetector.constant<double>("tpcEndplateServicesRing3_ro");

   tpcEndplateServices_R[3] = 
	theDetector.constant<double>("tpcEndplateServicesRing4_R");
   tpcEndplateServices_r[3] = 
	theDetector.constant<double>("tpcEndplateServicesRing4_ro");

   tpcEndplateServices_R[4] = 
	theDetector.constant<double>("tpcEndplateServicesRing5_R");
   tpcEndplateServices_r[4] = 
	theDetector.constant<double>("tpcEndplateServicesRing5_ro");

   tpcEndplateServices_R[5] = 
	theDetector.constant<double>("tpcEndplateServicesRing6_R");
   tpcEndplateServices_r[5] = 
	theDetector.constant<double>("tpcEndplateServicesRing6_ro");

   tpcEndplateServices_R[6] = 
	theDetector.constant<double>("tpcEndplateServicesRing7_R");
   tpcEndplateServices_r[6] = 
	theDetector.constant<double>("tpcEndplateServicesRing7_ro");

   tpcEndplateServices_R[7] = 
	theDetector.constant<double>("tpcEndplateServicesRing8_R");
   tpcEndplateServices_r[7] = 
	theDetector.constant<double>("tpcEndplateServicesRing8_ro");

   for(int i=0; i<MAX_TPC_RINGS; i++) {

      if(Ecal_cables_gap < 2*tpcEndplateServices_r[i])
	throw std::runtime_error("TPC endplate services: ring thickness exceeds Ecal_cables_gap");

      if(TPC_outer_radius < tpcEndplateServices_R[i] + tpcEndplateServices_r[i])
	throw std::runtime_error("TPC endplate services: ring dimensions exceed TPC_outer_radius");

      Torus ringSolid(tpcEndplateServices_R[i],
		      0,
		      tpcEndplateServices_r[i],
		      0, 2*M_PI);

      Volume ringLogical("ringLogical",ringSolid,
			 theDetector.material("copper"));

      double z_position = TPC_Ecal_Hcal_barrel_halfZ + 
	                  tpcEndplateServices_r[i];
      
#ifdef VERBOSE

G4cout << "Placing ring number " << i << " at z = " << z_position << G4endl;

#endif
      //first placement at -Z
      Position ringLogical_pos1(0,0,-z_position);
      pVol = envelope.placeVolume(ringLogical,ringLogical_pos1);

      //second placement at +Z
      Position ringLogical_pos2(0,0,z_position);
      pVol = envelope.placeVolume(ringLogical,ringLogical_pos2);
   }
 
  // Print out some parameters
  std::cout <<"\n   - TPC_Ecal_Hcal_barrel_halfZ = "
	    <<TPC_Ecal_Hcal_barrel_halfZ
	    <<"\n   - TPC_inner_radius = "
	    <<TPC_inner_radius
	    <<"\n   - TPC_outer_radius = "
	    <<TPC_outer_radius
	    <<"\n   - TUBE_IPOuterBulge_end_radius = "
	    <<TUBE_IPOuterBulge_end_radius
	    <<"\n   - TUBE_IPOuterBulge_end_z = "
	    <<TUBE_IPOuterBulge_end_z
	    <<"\n   - tpcEndplateServicesRing1_R = "
	    <<tpcEndplateServices_R[0]
	    <<"\n   - tpcEndplateServicesRing1_ro = "
	    <<tpcEndplateServices_r[0]
	    <<"\n   - tpcEndplateServicesRing2_R = "
	    <<tpcEndplateServices_R[1]
	    <<"\n   - tpcEndplateServicesRing2_ro = "
	    <<tpcEndplateServices_r[1]
	    <<"\n   - tpcEndplateServicesRing3_R = "
	    <<tpcEndplateServices_R[2]
	    <<"\n   - tpcEndplateServicesRing3_ro = "
	    <<tpcEndplateServices_r[2]
	    <<"\n   - tpcEndplateServicesRing4_R = "
	    <<tpcEndplateServices_R[3]
	    <<"\n   - tpcEndplateServicesRing4_ro = "
	    <<tpcEndplateServices_r[3]
	    <<"\n   - tpcEndplateServicesRing5_R = "
	    <<tpcEndplateServices_R[4]
	    <<"\n   - tpcEndplateServicesRing5_ro = "
	    <<tpcEndplateServices_r[4]
	    <<"\n   - tpcEndplateServicesRing6_R = "
	    <<tpcEndplateServices_R[5]
	    <<"\n   - tpcEndplateServicesRing6_ro = "
	    <<tpcEndplateServices_r[5]
	    <<"\n   - tpcEndplateServicesRing7_R = "
	    <<tpcEndplateServices_R[6]
	    <<"\n   - tpcEndplateServicesRing7_ro = "
	    <<tpcEndplateServices_r[6]
	    <<"\n   - tpcEndplateServicesRing8_R = "
	    <<tpcEndplateServices_R[7]
	    <<"\n   - tpcEndplateServicesRing8_ro = "
	    <<tpcEndplateServices_r[7]
	    <<std::endl;

   return true;
}

bool
SServices00_v01::
BuildEcalBarrelServices(PlacedVolume &pVol,
			Assembly &envelope,
			Detector &theDetector)
{

  Box ContainerSolid(top_dim_x/2.,
		     RailHeight/2.,
		     TPC_Ecal_Hcal_barrel_halfZ); 

  Volume containerLogical("EcalBarrelServicesContainerLogical",
			  ContainerSolid,
			  theDetector.material("air"));

  if(!FillEcalBarrelServicesContainer(pVol,containerLogical,theDetector))
	return false;

  for (int stave_id = 1; stave_id < 9 ; stave_id++)
  {
	double phirot = (stave_id-1) * M_PI/4;
	RotationZYX staveRot(-phirot,0,0);

	Position  stavePosition((Ecal_outer_radius + RailHeight/2.)*sin(phirot) + module_thickness * sin(M_PI/4.)*cos(-phirot), 
				(Ecal_outer_radius + RailHeight/2.)*cos(phirot) + module_thickness * sin(M_PI/4.)*sin(-phirot), 
				0 );
	Transform3D tran3D(staveRot,stavePosition);
	pVol = envelope.placeVolume(containerLogical,tran3D);
  }

  return true;
}

bool
SServices00_v01::
FillEcalBarrelServicesContainer(PlacedVolume &pVol,
				Volume &pContainerLogical,
				Detector &theDetector)
{
  double RailDistanceToRight = theDetector.constant<double>("EcalBarrelServices_RailDistanceToRight");

  double RailSeparation = theDetector.constant<double>("EcalBarrelServices_RailSeparation");

  double RailWidth = theDetector.constant<double>("EcalBarrelServices_RailWidth");

  bool RailSeparationChanged = false;
  if(fabs(
         top_dim_x/2. - (RailDistanceToRight + RailWidth*1.5 + RailSeparation)
         ) >= RailWidth*0.5)
  {
	RailSeparation = top_dim_x/2. - (RailDistanceToRight + RailWidth*1.5);
  	if(RailSeparation <= 0.)
		throw std::runtime_error("EcalBarrel services:no room for rails, cables, etc...");

        RailSeparationChanged = true;
  }

  Box railSolid(RailWidth/2.,
		RailHeight/2.,
		TPC_Ecal_Hcal_barrel_halfZ); 

  Volume railLogical("RailLogical",railSolid,
		     theDetector.material("aluminium"));
 
  double railPosition = top_dim_x/2. - RailDistanceToRight - RailWidth/2.;

  for(int i=1; i<=3; i++)
  {

    Position rail_pos(railPosition,0,0);

    pVol = pContainerLogical.placeVolume(railLogical,rail_pos);

    railPosition -= (RailWidth + RailSeparation);
  }

//-Z thicknesses 
  double ZMinus_FirstInterrail_PE_Thickness = theDetector.constant<double>("EcalBarrelServices_ZMinus_FirstInterrail_PE_Thickness");

  double ZMinus_FirstInterrail_Cu_Thickness = theDetector.constant<double>("EcalBarrelServices_ZMinus_FirstInterrail_Cu_Thickness");

  double ZMinus_SecondInterrail_Cu_Thickness = theDetector.constant<double>("EcalBarrelServices_ZMinus_SecondInterrail_Cu_Thickness");

//+Z thicknesses 
  double ZPlus_FirstInterrail_PE_Thickness = theDetector.constant<double>("EcalBarrelServices_ZPlus_FirstInterrail_PE_Thickness");

  double ZPlus_FirstInterrail_Cu_Thickness = theDetector.constant<double>("EcalBarrelServices_ZPlus_FirstInterrail_Cu_Thickness");

  double ZPlus_SecondInterrail_Cu_Thickness = theDetector.constant<double>("EcalBarrelServices_ZPlus_SecondInterrail_Cu_Thickness");

  if(RailSeparationChanged)
  {
     double OldRailSeparation = theDetector.constant<double>("EcalBarrelServices_RailSeparation");

     double fraction = OldRailSeparation / RailSeparation;
     ZMinus_FirstInterrail_PE_Thickness *= fraction;
     ZMinus_FirstInterrail_Cu_Thickness *= fraction;
     ZMinus_SecondInterrail_Cu_Thickness *= fraction;

     ZPlus_FirstInterrail_PE_Thickness *= fraction;
     ZPlus_FirstInterrail_Cu_Thickness *= fraction;
     ZPlus_SecondInterrail_Cu_Thickness *= fraction;
  }

  double moduleLength = TPC_Ecal_Hcal_barrel_halfZ * 2 / 5.;

  //First place ingredients at -Z:
  for(int i=0; i<3; i++)
  {
     double x_half_dim = (RailSeparation  * (3 - i) / 3.) / 2.;
     Box PESolid(x_half_dim,
		 ZMinus_FirstInterrail_PE_Thickness / 2.,
		 moduleLength / 2.); 

     Volume PELogical("PELogical",PESolid,
		      theDetector.material("polyethylene"));

     Position posPE(top_dim_x/2.-RailDistanceToRight-RailWidth-RailSeparation+x_half_dim,
		    -RailHeight/2. + ZMinus_FirstInterrail_PE_Thickness/2.,
		    -TPC_Ecal_Hcal_barrel_halfZ + moduleLength*(i+1/2.));

     pVol = pContainerLogical.placeVolume(PELogical,posPE);

     Box Cu_1_Solid(x_half_dim,
		    ZMinus_FirstInterrail_Cu_Thickness / 2.,
		    moduleLength / 2.); 

     Volume Cu_1_Logical("Cu_1_Logical",Cu_1_Solid,
                        theDetector.material("copper"));

     Position posCu1(top_dim_x/2.-RailDistanceToRight-RailWidth-RailSeparation+x_half_dim,
		     -RailHeight/2. + ZMinus_FirstInterrail_PE_Thickness+ZMinus_FirstInterrail_Cu_Thickness/2.,
		     -TPC_Ecal_Hcal_barrel_halfZ + moduleLength*(i+1/2.));
    
     pVol = pContainerLogical.placeVolume(Cu_1_Logical,posCu1);

     Box Cu_2_Solid(x_half_dim,
		    ZMinus_SecondInterrail_Cu_Thickness / 2.,
		    moduleLength / 2.); 

     Volume Cu_2_Logical("Cu_2_Logical",Cu_2_Solid,
                        theDetector.material("copper"));

     Position posCu2(top_dim_x/2.-RailDistanceToRight-2*RailWidth-2*RailSeparation+x_half_dim,
		     -RailHeight/2. + ZMinus_SecondInterrail_Cu_Thickness/2.,
		     -TPC_Ecal_Hcal_barrel_halfZ + moduleLength*(i+1/2.));
    
     pVol = pContainerLogical.placeVolume(Cu_2_Logical,posCu2);
  }
 
  //Now place ingredients at +Z:
  for(int i=0; i<2; i++)
  {
     double x_half_dim = (RailSeparation  * (2 - i) / 2.) / 2.;
     Box PESolid(x_half_dim,
		 ZPlus_FirstInterrail_PE_Thickness / 2.,
		 moduleLength / 2.); 

     Volume PELogical("PELogical",PESolid,
		      theDetector.material("polyethylene"));

     Position posPE(top_dim_x/2.-RailDistanceToRight-RailWidth-RailSeparation+x_half_dim,
		    -RailHeight/2. + ZPlus_FirstInterrail_PE_Thickness/2.,
		    TPC_Ecal_Hcal_barrel_halfZ - moduleLength*(i+1/2.));
    
     pVol = pContainerLogical.placeVolume(PELogical,posPE);


     Box Cu_1_Solid(x_half_dim,
		    ZPlus_FirstInterrail_Cu_Thickness / 2.,
		    moduleLength / 2.); 

     Volume Cu_1_Logical("Cu_1_Logical",Cu_1_Solid,
                        theDetector.material("copper"));

     Position posCu1(top_dim_x/2.-RailDistanceToRight-RailWidth-RailSeparation+x_half_dim,
		     -RailHeight/2. + ZPlus_FirstInterrail_PE_Thickness+ZPlus_FirstInterrail_Cu_Thickness/2.,
		     TPC_Ecal_Hcal_barrel_halfZ - moduleLength*(i+1/2.));

     pVol = pContainerLogical.placeVolume(Cu_1_Logical,posCu1);


     Box Cu_2_Solid(x_half_dim,
		    ZPlus_SecondInterrail_Cu_Thickness / 2.,
		    moduleLength / 2.); 

     Volume Cu_2_Logical("Cu_2_Logical",Cu_2_Solid,
                        theDetector.material("copper"));

     Position posCu2(top_dim_x/2.-RailDistanceToRight-2*RailWidth-2*RailSeparation+x_half_dim,
		     -RailHeight/2. + ZPlus_SecondInterrail_Cu_Thickness/2.,
		     TPC_Ecal_Hcal_barrel_halfZ - moduleLength*(i+1/2.));

     pVol = pContainerLogical.placeVolume(Cu_2_Logical,posCu2);

  }

  // Print out some parameters
  std::cout <<"\n   - EcalBarrelServices_RailDistanceToRight = "
	    <<RailDistanceToRight
	    <<"\n   - EcalBarrelServices_RailHeight = "
	    <<RailHeight
	    <<"\n   - EcalBarrelServices_RailSeparation = "
	    <<RailSeparation
	    <<"\n   - EcalBarrelServices_RailWidth = "
	    <<RailWidth
	    <<"\n   - EcalBarrelServices_ZMinus_FirstInterrail_Cu_Thickness = "
	    <<ZMinus_FirstInterrail_Cu_Thickness
	    <<"\n   - EcalBarrelServices_ZMinus_FirstInterrail_PE_Thickness = "
	    <<ZMinus_FirstInterrail_PE_Thickness
	    <<"\n   - EcalBarrelServices_ZMinus_SecondInterrail_Cu_Thickness = "
	    <<ZMinus_SecondInterrail_Cu_Thickness
	    <<"\n   - EcalBarrelServices_ZPlus_FirstInterrail_Cu_Thickness = "
	    <<ZPlus_FirstInterrail_Cu_Thickness
	    <<"\n   - EcalBarrelServices_ZPlus_FirstInterrail_PE_Thickness = "
	    <<ZPlus_FirstInterrail_PE_Thickness
	    <<"\n   - EcalBarrelServices_ZPlus_SecondInterrail_Cu_Thickness = "
	    <<ZPlus_SecondInterrail_Cu_Thickness 
	    <<std::endl;


  return true;
}

bool
SServices00_v01::
BuildEcalBarrel_EndCapServices(PlacedVolume &pVol,
			       Assembly &envelope,
			       Detector &theDetector)
{
//-Z thicknesses 
  double ZMinus_PE_Thickness = theDetector.constant<double>("EcalBarrel_EndCapServices_ZMinus_PE_Thickness");

  double ZMinus_Cu_Thickness = theDetector.constant<double>("EcalBarrel_EndCapServices_ZMinus_Cu_Thickness");

//+Z thicknesses 
  double ZPlus_PE_Thickness = theDetector.constant<double>("EcalBarrel_EndCapServices_ZPlus_PE_Thickness");

  double ZPlus_Cu_Thickness = theDetector.constant<double>("EcalBarrel_EndCapServices_ZPlus_Cu_Thickness");

  double containerThickness = ZMinus_PE_Thickness + ZMinus_Cu_Thickness;
  double z_position = -TPC_Ecal_Hcal_barrel_halfZ -containerThickness/2.;

  double Cu_Thickness = ZMinus_Cu_Thickness;
  double PE_Thickness = ZMinus_PE_Thickness;

  double container_x_dim = top_dim_x/2. + module_thickness*sin(M_PI/4.)-
		InnerServicesWidth/2.;

  for(int i=0; i<=1; i++)
  {

    if(Ecal_cables_gap < containerThickness)
      throw std::runtime_error("EcalBarrel_EndCap services: Cu+PE thicknesses exceed Ecal_cables_gap");

    Box ContainerSolid(container_x_dim/2.,
		       module_thickness/2.,
		       containerThickness/2.); 

    Volume containerLogical("EcalBarrel_EndCapServicesContainerLogical",ContainerSolid,
			    theDetector.material("air"));


    Box PESolid(container_x_dim/2.,
		module_thickness/2.,
		PE_Thickness/2.); 

    Volume PELogical("EcalBarrel_EndCap_PELogical",PESolid,
		     theDetector.material("polyethylene"));

    Position  PosPE(0,0,containerThickness/2. - PE_Thickness/2.);
      
    pVol = containerLogical.placeVolume(PELogical,PosPE);
 

    Box Cu_Solid(container_x_dim/2.,
		 module_thickness/2.,
		 Cu_Thickness/2.); 

    Volume Cu_Logical("EcalBarrel_EndCap_Cu_Logical",Cu_Solid,
		      theDetector.material("copper"));

    Position  PosCu(0,0,-containerThickness/2. + Cu_Thickness/2.);
      
    pVol = containerLogical.placeVolume(Cu_Logical,PosCu);


    for (int stave_id = 1; stave_id < 9 ; stave_id++)
    {
	double phirot = (stave_id-1) * M_PI/4;

	RotationZYX rot;
	Position stavePosition;

        if(z_position > 0) 
	  {
	    RotationZYX rot_P(-phirot,0,0);
	    rot = rot_P;
	    Position stavePosition_P((Ecal_outer_radius - module_thickness/2.)*sin(phirot)+(InnerServicesWidth/2. + container_x_dim/2.)*cos(-phirot),
				     (Ecal_outer_radius - module_thickness/2.)*cos(phirot)+(InnerServicesWidth/2. + container_x_dim/2.)*sin(-phirot), 
				     z_position);

	    stavePosition = stavePosition_P;
	  }
	else 
	  {
	    RotationZYX rot_M(phirot,M_PI,0);
	    rot = rot_M;
	    Position stavePosition_M((Ecal_outer_radius - module_thickness/2.)*sin(phirot)+(-InnerServicesWidth/2. - container_x_dim/2.)*cos(-phirot),
				     (Ecal_outer_radius - module_thickness/2.)*cos(phirot)+(-InnerServicesWidth/2. - container_x_dim/2.)*sin(-phirot),
				     z_position);
	    
	    stavePosition = stavePosition_M;
	  }

	Transform3D tran3D(rot,stavePosition);
	pVol = envelope.placeVolume(containerLogical,tran3D);

    }

    containerThickness = ZPlus_PE_Thickness + ZPlus_Cu_Thickness;
    z_position = TPC_Ecal_Hcal_barrel_halfZ + containerThickness/2.;
    
    Cu_Thickness = ZPlus_Cu_Thickness;
    PE_Thickness = ZPlus_PE_Thickness;

  }

  // Print out some parameters
  std::cout <<"\n   - EcalBarrel_EndCapServices_ZMinus_Cu_Thickness = "
	    <<ZMinus_Cu_Thickness
	    <<"\n   - EcalBarrel_EndCapServices_ZMinus_PE_Thickness = "
	    <<ZMinus_PE_Thickness
	    <<"\n   - EcalBarrel_EndCapServices_ZPlus_Cu_Thickness = "
	    <<ZPlus_Cu_Thickness
	    <<"\n   - EcalBarrel_EndCapServices_ZPlus_PE_Thickness = "
	    <<ZPlus_PE_Thickness
	    <<"\n   - Ecal_Tpc_gap = "
	    <<theDetector.constant<double>("Ecal_Tpc_gap")
	    <<"\n   - Ecal_cables_gap = "
	    <<Ecal_cables_gap
	    <<"\n   - Ecal_outer_radius = "
	    <<Ecal_outer_radius
	    <<std::endl;


  return true;
}

bool
SServices00_v01::
BuildHcalBarrel_EndCapServices(PlacedVolume &pVol,
			       Assembly &envelope,
			       Detector &theDetector)
{
  double Hcal_stave_gaps = 
	theDetector.constant<double>("Hcal_stave_gaps");
  double Hcal_inner_radius = Ecal_outer_radius + 
	theDetector.constant<double>("Hcal_Ecal_gap");
  double Hcal_R_max = 
	theDetector.constant<double>("Hcal_R_max");

  Hcal_total_dim_y = Hcal_R_max * cos(M_PI/16) - Hcal_inner_radius;
  double Hcal_module_radius = Hcal_inner_radius + Hcal_total_dim_y;
  Hcal_y_dim2_for_x  = 
	(Hcal_module_radius - Hcal_module_radius*cos(M_PI/8));
  double Hcal_y_dim1_for_x  = Hcal_total_dim_y - Hcal_y_dim2_for_x;
  Hcal_bottom_dim_x  = 
	2.*Hcal_inner_radius*tan(M_PI/8.)- Hcal_stave_gaps;
  Hcal_midle_dim_x   = 
	Hcal_bottom_dim_x + 2* Hcal_y_dim1_for_x*tan(M_PI/8.);
  Hcal_top_dim_x     = 
	Hcal_midle_dim_x - 2 * Hcal_y_dim2_for_x/tan(M_PI/8.);

#ifdef VERBOSE
  double Hcal_outer_radius = Hcal_inner_radius + Hcal_total_dim_y;

  G4cout << "BuildHcalBarrel_EndCapServices information: "
         << "\n                       Hcal_outer_radius = "
         << Hcal_outer_radius
         << "\n                       module thickness = "
         << Hcal_total_dim_y
         << "\n                       Hcal_R_max = "
         << Hcal_R_max
         << "\n                       Hcal_bottom_dim_x = "
         << Hcal_bottom_dim_x
	 << G4endl;
#endif

  double BHX  = Hcal_bottom_dim_x /2. - SurfaceTolerance;
  double MHX  = Hcal_midle_dim_x / 2. - SurfaceTolerance;
  double THX  = Hcal_top_dim_x / 2. - SurfaceTolerance;
  double YX1H = Hcal_y_dim1_for_x / 2.;
  double YX2H = Hcal_y_dim2_for_x / 2.;
  double DHZ  = Ecal_cables_gap / 2. - SurfaceTolerance;

  // Attention: on bâtit le module dans la verticale
  // à cause du G4Trd et on le tourne avant de le positioner
  Trapezoid Bottom(BHX, MHX, DHZ, DHZ, YX1H);

  Trapezoid Top(MHX, THX, DHZ, DHZ, YX2H);

  Position pos(0, 0, YX1H + YX2H);
  UnionSolid ModuleSolid( Bottom, Top, pos);


  Volume ModuleLogicalZMinus("ServicesHcalModuleZMinus", ModuleSolid,
			     theDetector.material("air"));


  Volume  ModuleLogicalZPlus("ServicesHcalModuleZPlus", ModuleSolid,
			     theDetector.material("air"));

  //First place layer models of services coming from Ecal and TPC:
  if(!FillHcalServicesModuleWithInnerServices(pVol,
		      ModuleLogicalZMinus,ModuleLogicalZPlus,theDetector))
      return false;

  //Then place layer models of HCAL electronics interface:
  int BuildHcalElectronicsInterface = 
    theDetector.constant<int>("BuildHcalElectronicsInterface");

  if(!(BuildHcalElectronicsInterface == 0))
    if(!FillHcalServicesModuleWithHcalElectronicsInterface(pVol,
		        ModuleLogicalZMinus,ModuleLogicalZPlus,theDetector))
      return false;

  double Y = Hcal_inner_radius + YX1H;
  double stave_phi_offset = 0;

  for (int stave_id = 1;
       stave_id <= 8;
       stave_id++)
    {
          double module_z_offset = 
		- TPC_Ecal_Hcal_barrel_halfZ - Ecal_cables_gap/2.;

          double phirot = stave_phi_offset;

	  RotationZYX rot(0,phirot,-M_PI*0.5);
	  Position stave_pos1(Y*sin(phirot),
			      Y*cos(phirot),
			      module_z_offset);
	  Transform3D tran3D1(rot,stave_pos1);
	  pVol = envelope.placeVolume(ModuleLogicalZMinus,tran3D1);

	  
          module_z_offset = - module_z_offset;

	  Position stave_pos2(Y*sin(phirot),
			      Y*cos(phirot),
			      module_z_offset);
	  Transform3D tran3D2(rot,stave_pos2);
	  pVol = envelope.placeVolume(ModuleLogicalZPlus,tran3D2);
	  	  
          stave_phi_offset -=  M_PI/4;
    }

   return true;
}

bool
SServices00_v01::
FillHcalServicesModuleWithInnerServices(PlacedVolume &pVol,
					Volume &ModuleLogicalZMinus,
					Volume &ModuleLogicalZPlus,
					Detector &theDetector)
{
//-Z thicknesses 
  double ZMinus_StainlessSteel_Thickness = theDetector.constant<double>("HcalServicesModule_ZMinus_StainlessSteel_Thickness");

  double ZMinus_PE_Thickness = theDetector.constant<double>("HcalServicesModule_ZMinus_PE_Thickness");

  double ZMinus_Cu_Thickness = theDetector.constant<double>("HcalServicesModule_ZMinus_Cu_Thickness");

//+Z thicknesses 
  double ZPlus_StainlessSteel_Thickness = theDetector.constant<double>("HcalServicesModule_ZPlus_StainlessSteel_Thickness");

  double ZPlus_PE_Thickness = theDetector.constant<double>("HcalServicesModule_ZPlus_PE_Thickness");

  double ZPlus_Cu_Thickness = theDetector.constant<double>("HcalServicesModule_ZPlus_Cu_Thickness");

  double StainlessSteel_Thickness = ZMinus_StainlessSteel_Thickness;
  double Cu_Thickness = ZMinus_Cu_Thickness;
  double PE_Thickness = ZMinus_PE_Thickness;

  Volume motherLogical = ModuleLogicalZMinus;

  double yShift = Hcal_y_dim2_for_x/2.;

  for(int i=0; i<=1; i++)
  {
    Position layerPositionSteel(0,yShift,
			   Ecal_cables_gap/2. - StainlessSteel_Thickness/2.);

    if(!PlaceHcalInnerServicesLayer(pVol,motherLogical,
		theDetector.material("stainless_steel"),
		StainlessSteel_Thickness, layerPositionSteel))

	return false;

    Position layerPositionPoly(0,yShift,
			   Ecal_cables_gap/2. -StainlessSteel_Thickness - PE_Thickness/2.);

    if(!PlaceHcalInnerServicesLayer(pVol,motherLogical,
		theDetector.material("polyethylene"),
		PE_Thickness, layerPositionPoly))

	return false;

    Position layerPositionCopper(0,yShift,
		Ecal_cables_gap/2.-StainlessSteel_Thickness-PE_Thickness -
		Cu_Thickness/2.);

    if(!PlaceHcalInnerServicesLayer(pVol,motherLogical,
		theDetector.material("copper"),
		Cu_Thickness, layerPositionCopper))

	return false;

    StainlessSteel_Thickness = ZPlus_StainlessSteel_Thickness;
    Cu_Thickness = ZPlus_Cu_Thickness;
    PE_Thickness = ZPlus_PE_Thickness;

    yShift = -Hcal_y_dim2_for_x/2.;

    motherLogical = ModuleLogicalZPlus;
  }

  // Print out some parameters
  std::cout <<"\n   - HcalServicesModule_InnerServicesWidth = "
	    <<InnerServicesWidth
	    <<"\n   - HcalServicesModule_ZMinus_Cu_Thickness = "
	    <<ZMinus_Cu_Thickness
	    <<"\n   - HcalServicesModule_ZMinus_PE_Thickness = "
	    <<ZMinus_PE_Thickness
	    <<"\n   - HcalServicesModule_ZMinus_StainlessSteel_Thickness = "
	    <<ZMinus_StainlessSteel_Thickness
	    <<"\n   - HcalServicesModule_ZPlus_Cu_Thickness = "
	    <<ZPlus_Cu_Thickness
	    <<"\n   - HcalServicesModule_ZPlus_PE_Thickness = "
	    <<ZPlus_PE_Thickness
	    <<"\n   - HcalServicesModule_ZPlus_StainlessSteel_Thickness = "
	    <<ZPlus_StainlessSteel_Thickness
	    <<"\n   - HcalServices_outer_Cu_thickness = "
	    <<theDetector.constant<double>("HcalServices_outer_Cu_thickness")
	    <<"\n   - HcalServices_outer_FR4_thickness = "
	    <<theDetector.constant<double>("HcalServices_outer_FR4_thickness")
	    <<"\n   - Hcal_Ecal_gap = "
	    <<theDetector.constant<double>("Hcal_Ecal_gap")
	    <<"\n   - Hcal_R_max = "
	    <<theDetector.constant<double>("Hcal_R_max")
	    <<"\n   - Hcal_back_plate_thickness = "
	    <<theDetector.constant<double>("Hcal_back_plate_thickness")
	    <<"\n   - Hcal_nlayers = "
	    <<theDetector.constant<int>("Hcal_nlayers")
	    <<"\n   - Hcal_radiator_thickness = "
	    <<theDetector.constant<double>("Hcal_radiator_thickness")
	    <<"\n   - Hcal_stave_gaps = "
	    <<theDetector.constant<int>("Hcal_stave_gaps")
	    <<"\n   - Hcal_steel_cassette_thickness = "
	    <<theDetector.constant<double>("Hcal_steel_cassette_thickness")
	    <<std::endl;

   return true;
}

bool
SServices00_v01::
FillHcalServicesModuleWithHcalElectronicsInterface(PlacedVolume &pVol,
						   Volume &ModuleLogicalZMinus,
						   Volume &ModuleLogicalZPlus,
						   Detector &theDetector)
{
   std::cout<<"   - BuildHcalElectronicsInterface = true " <<std::endl;

   double Hcal_back_plate_thickness = theDetector.constant<double>("Hcal_back_plate_thickness");
  
   double Hcal_nlayers = theDetector.constant<int>("Hcal_nlayers");

   double Hcal_radiator_thickness = theDetector.constant<double>("Hcal_radiator_thickness");

   double Hcal_layer_thickenss = 
	(Hcal_total_dim_y - Hcal_back_plate_thickness) / Hcal_nlayers;

   double Hcal_chamber_thickness = 
	Hcal_layer_thickenss - Hcal_radiator_thickness;

   if(Hcal_chamber_thickness <= 0)
      throw std::runtime_error("Hcal Barrel-EndCap services: Hcal chamber thicknesses  <= 0");

#ifdef VERBOSE
G4cout << "Hcal Barrel-EndCap services: Hcal_chamber_thickness = " <<
	Hcal_chamber_thickness << G4endl;
#endif

   double Hcal_y_dim1_for_x  = Hcal_total_dim_y - Hcal_y_dim2_for_x;

   double layer_x_dim = 0;
   double layer_y_offset = 0;

   for(int layer_id=1; layer_id<=Hcal_nlayers; layer_id++) {

	layer_y_offset = (layer_id - 1)*Hcal_layer_thickenss + 
                        Hcal_radiator_thickness;

        //---- bottom barrel----
	if(layer_id * Hcal_layer_thickenss  < Hcal_y_dim1_for_x )
		layer_x_dim = Hcal_bottom_dim_x + 2 * 
				layer_y_offset * tan(M_PI/8.);
	else   //----- top barrel ---
		layer_x_dim = Hcal_midle_dim_x - 2*
			( layer_y_offset + Hcal_chamber_thickness -
				Hcal_y_dim1_for_x ) / tan(M_PI/8.);

	if(!FillHcalElectronicsInterfaceLayer(pVol,
                ModuleLogicalZMinus, ModuleLogicalZPlus,
		layer_y_offset, layer_x_dim,theDetector))
		
		return false;
 
  }

   return true;
}

bool
SServices00_v01::
FillHcalElectronicsInterfaceLayer(PlacedVolume &pVol,
				  Volume &ModuleLogicalZMinus,
				  Volume &ModuleLogicalZPlus,
				  double layer_y_offset, double layer_x_dim,
				  Detector &theDetector)
{
   double Hcal_y_dim1_for_x  = Hcal_total_dim_y - Hcal_y_dim2_for_x;
   double Hcal_steel_cassette_thickness = theDetector.constant<double>("Hcal_steel_cassette_thickness");

   double y_position = -Hcal_y_dim1_for_x/2. + layer_y_offset +
		Hcal_steel_cassette_thickness/2.;

   if(!PlaceHcalElectronicsInterfaceComponent(pVol,
					      ModuleLogicalZMinus,
					      ModuleLogicalZPlus,
					      theDetector.material("S235"),
					      Hcal_steel_cassette_thickness,y_position,layer_x_dim)
      )
		return false;

   double FR4_thickness = theDetector.constant<double>("HcalServices_outer_FR4_thickness");

   y_position += (Hcal_steel_cassette_thickness/2. +
			FR4_thickness/2.);

   if(!PlaceHcalElectronicsInterfaceComponent(pVol,
					      ModuleLogicalZMinus,
					      ModuleLogicalZPlus,
					      theDetector.material("PCB"),
					      FR4_thickness,y_position,layer_x_dim)
      )
		return false;

   double Cu_thickness = theDetector.constant<double>("HcalServices_outer_Cu_thickness");

   y_position += (FR4_thickness/2. + Cu_thickness/2.);

   if(!PlaceHcalElectronicsInterfaceComponent(pVol,
					      ModuleLogicalZMinus,
					      ModuleLogicalZPlus,
					      theDetector.material("copper"),
					      Cu_thickness,y_position,layer_x_dim)
      )
		return false;

   return true;
}

bool
SServices00_v01::
PlaceHcalElectronicsInterfaceComponent(PlacedVolume &pVol,
				       Volume &ModuleLogicalZMinus,
				       Volume &ModuleLogicalZPlus,
				       Material layerMaterial,
				       double layerThickness,
				       double y_position, double layer_x_dim)
{
   Box layerBox (layer_x_dim/2. - 4*(SurfaceTolerance),
		 layerThickness/2. - 4*(SurfaceTolerance),
		 Ecal_cables_gap/2. - 4*(SurfaceTolerance)); 

   Solid layerSolid(layerBox);

   Solid cutSolid = CutLayer(layerSolid,y_position);

   Volume layerLogical("ElectronicsInterfaceLayerLogical",cutSolid,
		       layerMaterial);

   
   RotationZYX rot(0,0,M_PI*0.5);
   Position pos(0,0, y_position);
   Transform3D tran3D(rot,pos);

   pVol = ModuleLogicalZMinus.placeVolume(layerLogical,tran3D);
   

   pVol = ModuleLogicalZPlus.placeVolume(layerLogical,tran3D);

   return true;
}

Solid  SServices00_v01::CutLayer(Solid &layerSolid, double y_position)
{
    double tolerance = 2 * dd4hep::mm;

    Box Solid_Side(InnerServicesWidth/2. + tolerance,
		   Hcal_total_dim_y + 8*(SurfaceTolerance),
		   Ecal_cables_gap + 8*(SurfaceTolerance)); 

    double xShift = (Hcal_bottom_dim_x + Hcal_total_dim_y * tan(M_PI/8.)) / 2.;
    double yShift = Hcal_y_dim2_for_x/2. - y_position;

    Position rightSideLayerPosition(xShift,yShift,0);

    RotationZYX rightRotationSide(M_PI/8.,0,0);
    Transform3D tran3Dright(rightRotationSide,rightSideLayerPosition);

    SubtractionSolid subtractionRight(layerSolid,Solid_Side,tran3Dright);

    Position leftSideLayerPosition(-xShift,yShift,0);
    RotationZYX leftRotationSide(-M_PI/8.,0,0);
    Transform3D tran3Dleft(leftRotationSide,leftSideLayerPosition);

    SubtractionSolid subtractionLeft(subtractionRight,Solid_Side,tran3Dleft);

    Position centerLayerPosition(0,yShift,0);

    SubtractionSolid theFinalCut(subtractionLeft,Solid_Side,centerLayerPosition);

    return theFinalCut;

}

bool
SServices00_v01::
PlaceHcalInnerServicesLayer(PlacedVolume &pVol,Volume &motherLogical,
			    Material layerMaterial,
			    double layerThickness,
			    Position &layerPosition)
{
    RotationZYX rot(0,0,M_PI*0.5);

    Box Center_Solid(InnerServicesWidth/2. - 2*(SurfaceTolerance),
		     Hcal_total_dim_y/2. - 2*(SurfaceTolerance),
		     layerThickness/2. - 2*(SurfaceTolerance)); 

    Volume Center_Logical("HcalBarrel_EndCap_Center_Logical",Center_Solid,
			  layerMaterial);

    Transform3D tran3D(rot,layerPosition);
    pVol = motherLogical.placeVolume(Center_Logical,tran3D);
        
    Box Box_Side(InnerServicesWidth/2. - 2*(SurfaceTolerance),
		 Hcal_total_dim_y - 4*(SurfaceTolerance),
		 layerThickness/2. - 2*(SurfaceTolerance)); 

    Solid Solid_Side(Box_Side);

    Solid motherSolid = motherLogical.solid();

    double xShift = (Hcal_bottom_dim_x + Hcal_total_dim_y * tan(M_PI/8.)) / 2.;

    Position sideLayerPosition = layerPosition;
    sideLayerPosition.SetX(xShift);
 
    RotationZYX  rotRight(-M_PI/8.,0,M_PI*0.5);
    Transform3D tran3DRight(rotRight,sideLayerPosition);

    IntersectionSolid intersectionSolidRight(motherSolid,Solid_Side,tran3DRight);

    Volume Logical_Right("HcalBarrel_EndCap_Logical_Right",intersectionSolidRight,
			layerMaterial);

    Position pos(0,0,0);
    pVol = motherLogical.placeVolume(Logical_Right,pos);
    
    sideLayerPosition = layerPosition;
    sideLayerPosition.SetX(-xShift);

    RotationZYX  rotLeft(M_PI/8.,0,M_PI*0.5);
    Transform3D tran3DLeft(rotLeft,sideLayerPosition);

    IntersectionSolid intersectionSolidLeft(motherSolid,Solid_Side,tran3DLeft);

    Volume Logical_Left("HcalBarrel_EndCap_Logical_Left",intersectionSolidLeft,
			layerMaterial);

    pVol = motherLogical.placeVolume(Logical_Left,pos);
    
   return true;
}


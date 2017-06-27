#include <utility>
#include <vector>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"


//=====================================
// class for  BuildTPCEndplateServices
//=====================================
//using vector of double pair to hold the value 
//of "tpcEndplateServicesRing_R" and "tpcEndplateServicesRing_ro".
//Which provide more flexible for handling the modification in the number
//of TPC cooling rings. This class can handle different number of rings.

class BuildTPCEndplateServices{
public:
  BuildTPCEndplateServices(){};
  ~BuildTPCEndplateServices(){};

  void setMaterial (dd4hep::Material copper4cooling) {cooling_Material = copper4cooling;};
  void sethalfZ (double TPC_Ecal_Hcal_barrel_halfZ){TPC_barrel_halfZ = TPC_Ecal_Hcal_barrel_halfZ ;};
  void settpcEndplateServicesRing_R_ro( double tpcEndplateServicesRing_R, double tpcEndplateServicesRing_ro) {
    tpcEndplateServicesRing_R_ro.push_back(std::make_pair(tpcEndplateServicesRing_R,tpcEndplateServicesRing_ro));
  };

  // to build TPC cooling rings into the service assembly 
  bool DoBuildTPCEndplateServices(dd4hep::PlacedVolume &pVol,dd4hep::Assembly &envelope){
    for( std::vector< std::pair<double,double> >::const_iterator it = tpcEndplateServicesRing_R_ro.begin() ;
	 it != tpcEndplateServicesRing_R_ro.end(); ++it){
      dd4hep::Torus solidTube(it->first, 0, it->second, 0, 2*M_PI);
      double z_position = TPC_barrel_halfZ + it->second;
      dd4hep::Volume EnvTube("service_tube",solidTube,cooling_Material);
      //EnvTube.setVisAttributes("RedVis");
      dd4hep::Position pos1(0,0,z_position);
      pVol = envelope.placeVolume(EnvTube,pos1);
      dd4hep::Position pos2(0,0,-z_position);
      pVol = envelope.placeVolume(EnvTube,pos2);
    } 
    return true;
  };

private:
  dd4hep::Material cooling_Material;
  double TPC_barrel_halfZ;
  std::vector< std::pair<double,double> > tpcEndplateServicesRing_R_ro;
};



//=====================================
// class for  BuildEcalBarrelServices
//=====================================
class BuildEcalBarrelServices {

public:
  BuildEcalBarrelServices(){};
  ~BuildEcalBarrelServices(){};

  void setMaterialAir (dd4hep::Material air4container) {air = air4container;};
  void setMaterialAluminium (dd4hep::Material al) {aluminium = al;};
  void setMaterialPolyethylene (dd4hep::Material PE) {polyethylene = PE;};
  void setMaterialCopper (dd4hep::Material Cu) {copper = Cu;};

  void sethalfZ (double TPC_Ecal_Hcal_barrel_halfZ){Ecal_barrel_halfZ = TPC_Ecal_Hcal_barrel_halfZ ;};
  void setTopDimX (double dim_x){top_dim_x = dim_x ;};
  void setRailHeight (double height){RailHeight = height ;};
  void setOutRadius (double outer_radius){Ecal_outer_radius = outer_radius ;};
  void setModuleThickness (double thickness){ module_thickness = thickness ;};

  void setRailSeparation (double separation) {RailSeparation = separation;};
  void setRailWidth (double width) {RailWidth = width;};

  void setZMinus_FirstInterrail_PE_Thickness  (double PE_Thickness_ZMF) {ZMinus_FirstInterrail_PE_Thickness = PE_Thickness_ZMF;};
  void setZMinus_FirstInterrail_Cu_Thicknesst  (double Cu_Thickness_ZMF) {ZMinus_FirstInterrail_Cu_Thickness = Cu_Thickness_ZMF;};

  void setZPlus_FirstInterrail_PE_Thickness  (double PE_Thickness_ZPF) {ZPlus_FirstInterrail_PE_Thickness = PE_Thickness_ZPF;};
  void setZPlus_FirstInterrail_Cu_Thickness  (double Cu_Thickness_ZPF) {ZPlus_FirstInterrail_Cu_Thickness = Cu_Thickness_ZPF;};

  void setRailSeparationChanged(void) {RailSeparationChanged = false;};

  // to fill the detail service layers into the container
  bool FillEcalBarrelServicesContainer(dd4hep::PlacedVolume &pVol, dd4hep::Volume &pContainerLogical){

    const int NRAILS=2;
    if ( RailSeparation*(NRAILS-1) > top_dim_x - RailWidth ) {
      std::cout << "WARNING, requested rail separation too large! " << std::endl;
      std::cout << "Must be smaller than " << (top_dim_x - RailWidth)/(NRAILS-1) << " cm" << std::endl;
      assert(0);
    }

    float railPositions[NRAILS]; // positions wrt module centre: equally spaced, centered on module
    for (int i=0; i<NRAILS; i++) {
      railPositions[i] = i*RailSeparation - RailSeparation*(NRAILS-1)/2.;
    }
   
    dd4hep::Box railSolid(RailWidth/2., RailHeight/2., Ecal_barrel_halfZ); 

    dd4hep::Volume railLogical("railLogical", railSolid, aluminium);

    for(int i=0; i<NRAILS; i++)
      {
	dd4hep::Position pos(railPositions[i],0,0);
	pVol = pContainerLogical.placeVolume(railLogical,pos);
      }
    
  double moduleLength = Ecal_barrel_halfZ * 2. / 5.;

  double internalZoneWidth = RailSeparation - RailWidth - 0.1; // space between rails (in x-dir) (minus a small gap)

  // now fit fit service materials into the zone between rails

  const int NMODULES_IN_Z=5;
  for(int i=0; i<NMODULES_IN_Z; i++){ // these are the 5 modules in z
    
    float pe_thick = ZMinus_FirstInterrail_PE_Thickness;
    float cu_thick = ZMinus_FirstInterrail_Cu_Thickness;
    if ( i>2 ) {
      pe_thick = ZPlus_FirstInterrail_PE_Thickness;
      cu_thick = ZPlus_FirstInterrail_Cu_Thickness;
    }

    // polyethylene
    dd4hep::Box PESolid( internalZoneWidth/2. , pe_thick / 2., moduleLength / 2.); 
    
    std::string PELogical_name  = "PELogical"+dd4hep::_toString(i,"_%d");
    dd4hep::Volume PELogical(PELogical_name,PESolid,polyethylene);
    //PELogical.setVisAttributes("GreenVis");

    dd4hep::Position posPE( 0 , 
                            -RailHeight/2. + pe_thick/2.,
                            -Ecal_barrel_halfZ + moduleLength*(i+1/2.));
    
    pVol = pContainerLogical.placeVolume(PELogical,posPE);
    
    // Cu_1
    dd4hep::Box Cu_1_Solid( internalZoneWidth/2., cu_thick / 2., moduleLength / 2.); 
    
    std::string Cu_1_Logical_name  = "Cu_1_Logical"+dd4hep::_toString(i,"_%d");
    dd4hep::Volume Cu_1_Logical(Cu_1_Logical_name, Cu_1_Solid,copper);
    //Cu_1_Logical.setVisAttributes("BlueVis");

    dd4hep::Position posCu1( 0 , 
                             -RailHeight/2. + pe_thick + cu_thick/2.,
                             -Ecal_barrel_halfZ + moduleLength*(i+1/2.));
    
    pVol = pContainerLogical.placeVolume(Cu_1_Logical,posCu1);

  }

    return true;
  };

  // to build Ecal Barrel service into the service assembly 
  bool DoBuildEcalBarrelServices(dd4hep::PlacedVolume &pVol,dd4hep::Assembly &envelope){

    dd4hep::Box ContainerSolid(top_dim_x/2., RailHeight/2., Ecal_barrel_halfZ); 
    dd4hep::Volume containerLogical("EcalBarrelServicesContainerLogical",ContainerSolid,air);
    //containerLogical.setVisAttributes("SeeThrough");
    //containerLogical.setVisAttributes("MagentaVis");

    if(!FillEcalBarrelServicesContainer(pVol,containerLogical))
      return false;

    for (int stave_id = 1; stave_id < 9 ; stave_id++)
      {
	double phirot = (stave_id-1) * M_PI/4.;
	dd4hep::RotationZYX rot(-phirot,0,0);
	//Position  stavePosition(module_thickness * sin(M_PI/4.), Ecal_outer_radius + RailHeight/2., 0);
	dd4hep::Position stavePosition((Ecal_outer_radius + RailHeight/2.)*sin(phirot) - module_thickness * sin(M_PI/4.)*cos(-phirot), 
                                       (Ecal_outer_radius + RailHeight/2.)*cos(phirot) - module_thickness * sin(M_PI/4.)*sin(-phirot), 
				0 );
	dd4hep::Transform3D tran3D(rot,stavePosition);
	pVol = envelope.placeVolume(containerLogical,tran3D);
      }
    
    return true;
  };

private:
  double Ecal_barrel_halfZ;
  double top_dim_x;
  double RailHeight;
  double Ecal_outer_radius;
  double module_thickness;

  //  double RailDistanceToRight;
  double RailSeparation;
  double RailWidth;
  bool   RailSeparationChanged;

//-Z thicknesses
  double ZMinus_FirstInterrail_PE_Thickness;
  double ZMinus_FirstInterrail_Cu_Thickness;

//+Z thicknesses 
  double ZPlus_FirstInterrail_PE_Thickness;
  double ZPlus_FirstInterrail_Cu_Thickness;

  dd4hep::Material air;
  dd4hep::Material aluminium;
  dd4hep::Material polyethylene;
  dd4hep::Material copper;

};






//=====================================
// class for BuildEcalBarrel_EndCapServices
//=====================================
class BuildEcalBarrel_EndCapServices {

public:
  BuildEcalBarrel_EndCapServices(){};
  ~BuildEcalBarrel_EndCapServices(){};

  void setMaterialAir (dd4hep::Material air4container) {air = air4container;};
  void setMaterialPolyethylene (dd4hep::Material PE) {polyethylene = PE;};
  void setMaterialCopper (dd4hep::Material Cu) {copper = Cu;};

  void sethalfZ (double TPC_Ecal_Hcal_barrel_halfZ){Ecal_barrel_halfZ = TPC_Ecal_Hcal_barrel_halfZ ;};
  void setTopDimX (double dim_x){top_dim_x = dim_x ;};
  void setModuleThickness (double thickness){ module_thickness = thickness ;};
  void setInnerServicesWidth (double width){InnerServicesWidth = width;};
  void setEcal_cables_gap(double gap){Ecal_cables_gap = gap;};
  void setOutRadius (double outer_radius){Ecal_outer_radius = outer_radius ;};

  void setZMinus_PE_Thickness (double PE_Thickness_ZM) {ZMinus_PE_Thickness = PE_Thickness_ZM;};
  void setZMinus_Cu_Thickness (double Cu_Thickness_ZM) {ZMinus_Cu_Thickness = Cu_Thickness_ZM;};

  void setZPlus_PE_Thickness (double PE_Thickness_ZP) {ZPlus_PE_Thickness = PE_Thickness_ZP;};
  void setZPlus_Cu_Thickness (double Cu_Thickness_ZP) {ZPlus_Cu_Thickness = Cu_Thickness_ZP;};


  // to build Ecal Barrel service into the service assembly 
  bool DoBuildEcalBarrel_EndCapServices(dd4hep::PlacedVolume &pVol,dd4hep::Assembly &envelope){

    double containerThickness = ZMinus_PE_Thickness + ZMinus_Cu_Thickness;
    double z_position = -Ecal_barrel_halfZ -containerThickness/2.;
    double container_x_dim = top_dim_x/2. + module_thickness*sin(M_PI/4.) - InnerServicesWidth/2.;
    double Cu_Thickness = ZMinus_Cu_Thickness;
    double PE_Thickness = ZMinus_PE_Thickness;

    for(int i=0; i<=1; i++){

      if(Ecal_cables_gap < containerThickness) return false;
	//Control::Abort(
	//	       "EcalBarrel_EndCap services: Cu+PE thicknesses exceed Ecal_cables_gap",
	//	       MOKKA_ERROR_BAD_DATABASE_PARAMETERS);

      dd4hep::Box ContainerSolid(container_x_dim/2., module_thickness/2., containerThickness/2.); 

      std::string containerLogical_name  = "containerLogical"+dd4hep::_toString(i,"_%d");
      dd4hep::Volume containerLogical(containerLogical_name,ContainerSolid,air);


      dd4hep::Box PESolid(container_x_dim/2., module_thickness/2., PE_Thickness/2.);
 
      std::string PELogical_name  = "PELogical"+dd4hep::_toString(i,"_%d");
      dd4hep::Volume PELogical(PELogical_name,PESolid,polyethylene);

      dd4hep::Position  PosPE(0,0,containerThickness/2. - PE_Thickness/2.);
      
      pVol = containerLogical.placeVolume(PELogical,PosPE);
      
      dd4hep::Box Cu_Solid(container_x_dim/2., module_thickness/2., Cu_Thickness/2.); 
      
      std::string Cu_Logical_name  = "Cu_Logical"+dd4hep::_toString(i,"_%d");
      dd4hep::Volume Cu_Logical(Cu_Logical_name,Cu_Solid,copper);
      
      dd4hep::Position  PosCu(0,0,-containerThickness/2. + Cu_Thickness/2.);
      
      pVol = containerLogical.placeVolume(Cu_Logical,PosCu);

      
      for (int stave_id = 1; stave_id < 9 ; stave_id++){
	
	double phirot = (stave_id-1) * M_PI/4;
	
	dd4hep::RotationZYX rot;
	dd4hep::Position stavePosition;
	
        if(z_position > 0) {
	  dd4hep::RotationZYX rot_P(-phirot,0,0);
	  rot = rot_P;
	  dd4hep::Position stavePosition_P((Ecal_outer_radius - module_thickness/2.)*sin(phirot)+(InnerServicesWidth/2. + container_x_dim/2.)*cos(-phirot),
                                           (Ecal_outer_radius - module_thickness/2.)*cos(phirot)+(InnerServicesWidth/2. + container_x_dim/2.)*sin(-phirot), 
                                           z_position);

	  stavePosition = stavePosition_P;
	}
	else
	  {
	    dd4hep::RotationZYX rot_M(phirot,M_PI,0);
	    rot = rot_M;
	    dd4hep::Position stavePosition_M((Ecal_outer_radius - module_thickness/2.)*sin(phirot)+(-InnerServicesWidth/2. - container_x_dim/2.)*cos(-phirot),
                                             (Ecal_outer_radius - module_thickness/2.)*cos(phirot)+(-InnerServicesWidth/2. - container_x_dim/2.)*sin(-phirot),
                                             z_position);
	    
	    stavePosition = stavePosition_M;
	  }
	
	dd4hep::Transform3D tran3D(rot,stavePosition);
	pVol = envelope.placeVolume(containerLogical,tran3D);
	
      }
      
      containerThickness = ZPlus_PE_Thickness + ZPlus_Cu_Thickness;
      z_position = Ecal_barrel_halfZ + containerThickness/2.;
      
      Cu_Thickness = ZPlus_Cu_Thickness;
      PE_Thickness = ZPlus_PE_Thickness;
      
    }


    return true;
  };


private:
  double Ecal_barrel_halfZ;
  double top_dim_x;
  double module_thickness;
  double InnerServicesWidth;
  double Ecal_cables_gap;
  double Ecal_outer_radius;

//-Z thicknesses
  double ZMinus_PE_Thickness;
  double ZMinus_Cu_Thickness;

//+Z thicknesses 
  double ZPlus_PE_Thickness;
  double ZPlus_Cu_Thickness;

  dd4hep::Material air;
  dd4hep::Material polyethylene;
  dd4hep::Material copper;


};

//=====================================
// class for BuildSitCables
//=====================================
class BuildSitCables {

public:
  BuildSitCables(){};
  ~BuildSitCables(){};

  void setMaterialAluminium (dd4hep::Material al) {aluminium = al;};
  void setSit_cables_cylinder_thickness(double cables_cylinder_thickness){Sit_cables_cylinder_thickness = cables_cylinder_thickness ;};
  void sethalfZ (double TPC_barrel_halfZ){TPC_Ecal_Hcal_barrel_halfZ = TPC_barrel_halfZ ;};
  void setTPC_Inner_Radius (double TPCInnerRadius){TPC_inner_radius = TPCInnerRadius ;};
  void setftd4to7_tpc_radial_gap (double gap){ftd4to7_tpc_radial_gap = gap ;};
  void setz2_position_ReltoTPCLength (double z2_position){z2_position_ReltoTPCLength = z2_position ;};
  void setz3_position_ReltoTPCLength (double z3_position){z3_position_ReltoTPCLength = z3_position ;};
  void setpetal_cp_support_thickness (double support_thickness){petal_cp_support_thickness = support_thickness ;};
  void setdisk_si_thickness (double si_thickness){disk_si_thickness = si_thickness ;};
  void setpetal_support_zoffset (double support_zoffset){petal_support_zoffset = support_zoffset ;};
  void setftd2_sit1_radial_diff (double ftd2_sit1_radial){ftd2_sit1_radial_diff = ftd2_sit1_radial ;};
  void setftd3_sit2_radial_diff (double ftd3_sit2_radial){ftd3_sit2_radial_diff = ftd3_sit2_radial ;};
  void setSIT1_Radius (double SIT1Radius){SIT1_Radius =  SIT1Radius;};
  void setSIT2_Radius (double SIT2Radius){SIT2_Radius =  SIT2Radius;};
  void setFTD2_cone_thickness (double FTD2_thickness){FTD2_cone_thickness =  FTD2_thickness;};
  void setFTD3_cone_thickness (double FTD3_thickness){FTD3_cone_thickness = FTD3_thickness ;};
  void setTUBE_IPOuterBulge_end_radius (double end_radius){TUBE_IPOuterBulge_end_radius = end_radius ;};
  void setTUBE_IPOuterBulge_end_z (double end_z){TUBE_IPOuterBulge_end_z = end_z;};
  void setSServices_FTD7_cables_thickness (double cable_thickness){SServices_FTD7_cables_thickness = cable_thickness;};

  // to build Ecal Barrel service into the service assembly
  bool DoBuildSitCables(dd4hep::PlacedVolume &pVol, dd4hep::Assembly &envelope){

    //First place the Cylinder

    if(ftd4to7_tpc_radial_gap < Sit_cables_cylinder_thickness)
      throw std::runtime_error("SServices00: the ftd-tpc radial gap is less than Sit cable thickness");

    double SitTube_inner_radius = TPC_inner_radius -
      ftd4to7_tpc_radial_gap/2.;

    double z_start_3 = TPC_Ecal_Hcal_barrel_halfZ * z3_position_ReltoTPCLength;

    double petalairthickness_half = 0.5 * ( petal_cp_support_thickness + 2.0*disk_si_thickness);

    double max_half_thickness_disk_3 = petal_support_zoffset + petalairthickness_half ;

    double z_half_len = (TPC_Ecal_Hcal_barrel_halfZ -
                        z_start_3) / 2.;

    dd4hep::Tube SitTubeSolid (SitTube_inner_radius,
			       SitTube_inner_radius + Sit_cables_cylinder_thickness,
			       z_half_len,
			       0., 2 * M_PI);

    dd4hep::Volume SitTubeLog("SitTubeLog",SitTubeSolid, aluminium);

    dd4hep::Position Cylinder_pos1(0, 0, z_start_3 + z_half_len);
    pVol = envelope.placeVolume(SitTubeLog,Cylinder_pos1);

    dd4hep::Position Cylinder_pos2(0, 0, -z_start_3 - z_half_len);
    pVol = envelope.placeVolume(SitTubeLog,Cylinder_pos2);


    //Then place the cone

    double z_start_2 = TPC_Ecal_Hcal_barrel_halfZ * z2_position_ReltoTPCLength;

    petalairthickness_half = 0.5 * ( petal_cp_support_thickness + 2.0*disk_si_thickness);

    double max_half_thickness_disk_2 = petal_support_zoffset + petalairthickness_half ;

    double FTD2_outer_radius = SIT1_Radius + ftd2_sit1_radial_diff;

    double FTD3_outer_radius = SIT2_Radius + ftd3_sit2_radial_diff;

    double zPlane[2];
    zPlane[0] = z_start_2 + max_half_thickness_disk_2 + SurfaceTolerance;

    zPlane[1] = z_start_3 - max_half_thickness_disk_3 - SurfaceTolerance;

    double rInner[2];
    rInner[0] = FTD2_outer_radius + SurfaceTolerance;
    rInner[1] = FTD3_outer_radius + SurfaceTolerance;

    double rOuter[2];
    rOuter[0] = FTD2_outer_radius + FTD2_cone_thickness;
    rOuter[1] = FTD3_outer_radius + FTD3_cone_thickness;

    vector<double> rmin;
    vector<double> rmax;
    vector<double> z;
    for(int i=0;i<=1;i++) {
      rmin.push_back(rInner[i]);
      rmax.push_back(rOuter[i]);
      z.push_back(zPlane[i]);
    }

    dd4hep::Polycone SitConeSolid (0., 2.0 * M_PI, rmin, rmax, z);

    dd4hep::Volume SitConeLog("SitConeLog", SitConeSolid, aluminium);

    dd4hep::Position cone_pos(0, 0, 0);
    pVol = envelope.placeVolume(SitConeLog,cone_pos);

    dd4hep::RotationZYX rot(0,M_PI,0);
    dd4hep::Transform3D tran3D(rot,cone_pos);
    pVol = envelope.placeVolume(SitConeLog,tran3D);


    //Then place the disk
    double Sit_cables_disk_thickness =((1+TUBE_IPOuterBulge_end_radius/TPC_inner_radius)*0.5*
                                      SServices_FTD7_cables_thickness );

    dd4hep::Tube SitDiskSolid(TUBE_IPOuterBulge_end_radius + SurfaceTolerance,
			      TPC_inner_radius,
			      Sit_cables_disk_thickness/2. - SurfaceTolerance,
			      0., 2 * M_PI);

    dd4hep::Volume SitDiskLog("SitDiskLog", SitDiskSolid, aluminium);

    dd4hep::Position disk_pos1(0, 0, TUBE_IPOuterBulge_end_z + Sit_cables_disk_thickness/2.);
    pVol = envelope.placeVolume(SitDiskLog,disk_pos1);

    dd4hep::Position disk_pos2(0, 0, -TUBE_IPOuterBulge_end_z - Sit_cables_disk_thickness/2.);
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
	      <<"\n   - TUBE_IPOuterBulge_end_z = "
	      <<TUBE_IPOuterBulge_end_z
	      <<"\n   - SServices_FTD7_cables_thickness = "
	      <<SServices_FTD7_cables_thickness
	      <<std::endl;


    return true;
  };

private:
  double Sit_cables_cylinder_thickness;
  double TPC_inner_radius;
  double TPC_Ecal_Hcal_barrel_halfZ;
  double SIT1_Radius;
  double SIT2_Radius;
  double FTD2_cone_thickness;
  double FTD3_cone_thickness;
  double TUBE_IPOuterBulge_end_radius;
  double SServices_FTD7_cables_thickness;
  double TUBE_IPOuterBulge_end_z;
  double z2_position_ReltoTPCLength;
  double z3_position_ReltoTPCLength;
  double petal_cp_support_thickness;
  double disk_si_thickness;
  double petal_support_zoffset;
  double ftd2_sit1_radial_diff;
  double ftd3_sit2_radial_diff;
  double ftd4to7_tpc_radial_gap;

  const double SurfaceTolerance = 0.0001*dd4hep::mm;

  dd4hep::Material aluminium;

};

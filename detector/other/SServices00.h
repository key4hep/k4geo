#include <utility>
#include <vector>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetElement.h"

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

  /// to build TPC cooling rings into the service assembly 
  bool DoBuildTPCEndplateServices(dd4hep::PlacedVolume &pVol,dd4hep::Assembly &envelope) ;

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
  void setenv_safety(double safety) {
    env_safety = safety ;
    RailHeight -= env_safety*4.0 ;
  };

  /// to fill the detail service layers into the container
  bool FillEcalBarrelServicesContainer(dd4hep::PlacedVolume &pVol, dd4hep::Volume &pContainerLogical) ;


  /// to build Ecal Barrel service into the service assembly 
  bool DoBuildEcalBarrelServices(dd4hep::PlacedVolume &pVol,dd4hep::Assembly &envelope) ;

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

  double env_safety;
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

  void setenv_safety(double safety) {
    env_safety = safety ;
    Ecal_cables_gap -= env_safety*4.0 ;
  };

  /// to build Ecal Barrel service into the service assembly 
  bool DoBuildEcalBarrel_EndCapServices(dd4hep::PlacedVolume &pVol,dd4hep::Assembly &envelope) ;

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

  double env_safety;
};

//=====================================
// class for BuildSitCables
//=====================================
class BuildSitCables {

public:
  BuildSitCables(dd4hep::DetElement det) : detElt( det) { };
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
  void setTUBE_IPOuterBulge_end_envradius (double end_envradius){TUBE_IPOuterBulge_end_envradius = end_envradius ;};
  void setTUBE_IPOuterBulge_end_z (double end_z){TUBE_IPOuterBulge_end_z = end_z;};
  void setSServices_FTD7_cables_thickness (double cable_thickness){SServices_FTD7_cables_thickness = cable_thickness;};

  /// to build Ecal Barrel service into the service assembly
  bool DoBuildSitCables(dd4hep::PlacedVolume &pVol, dd4hep::Assembly &envelope) ;


private:
  double Sit_cables_cylinder_thickness;
  double TPC_inner_radius;
  double TPC_Ecal_Hcal_barrel_halfZ;
  double SIT1_Radius;
  double SIT2_Radius;
  double FTD2_cone_thickness;
  double FTD3_cone_thickness;
  double TUBE_IPOuterBulge_end_radius;
  double TUBE_IPOuterBulge_end_envradius;
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
  dd4hep::DetElement detElt ;
};



//=====================================
// class for BuildHcalBarrel_EndCapServices
//=====================================
class BuildHcalBarrel_EndCapServices {

public:
  BuildHcalBarrel_EndCapServices(){};
  ~BuildHcalBarrel_EndCapServices(){};

public:

  void setMaterialAir (dd4hep::Material Air) {air = Air;};
  void setMaterialPolyethylene (dd4hep::Material PE) {polyethylene = PE;};
  void setMaterialCopper (dd4hep::Material Cu) {copper = Cu;};
  void setMaterialstainless_steel (dd4hep::Material Fe) {stainless_steel = Fe;};
  void setMaterialS235 (dd4hep::Material steel) {S235 = steel;};
  void setMaterialPCB (dd4hep::Material pcb) {PCB = pcb;};

  void setHcal_stave_gaps (double stave_gaps){Hcal_stave_gaps = stave_gaps ;};
  void setHcal_inner_radius (double inner_radius){Hcal_inner_radius = inner_radius ;};
  void setHcal_R_max (double R_max){Hcal_R_max = R_max ;};
  void setBuildHcalElectronicsInterface (int build){BuildHcalElectronicsInterface = build ;};
  void sethalfZ (double Hcal_barrel_halfZ){TPC_Ecal_Hcal_barrel_halfZ = Hcal_barrel_halfZ ;};
  void setEcal_cables_gap(double gap){Ecal_cables_gap = gap;};

//-Z thicknesses
  void setZMinus_StainlessSteel_Thickness (double ZMinus_StainlessSteel){ZMinus_StainlessSteel_Thickness = ZMinus_StainlessSteel ;};
  void setZMinus_PE_Thickness (double ZMinus_PE){ZMinus_PE_Thickness = ZMinus_PE ;};
  void setZMinus_Cu_Thickness (double ZMinus_Cu){ZMinus_Cu_Thickness = ZMinus_Cu ;};
//+Z thicknesses
  void setZPlus_StainlessSteel_Thickness (double ZPlus_StainlessSteel){ZPlus_StainlessSteel_Thickness = ZPlus_StainlessSteel ;};
  void setZPlus_PE_Thickness (double ZPlus_PE){ZPlus_PE_Thickness = ZPlus_PE ;};
  void setZPlus_Cu_Thickness (double ZPlus_Cu){ZPlus_Cu_Thickness = ZPlus_Cu ;};

  void setInnerServicesWidth (double ServicesWidth){InnerServicesWidth = ServicesWidth ;};
  void setHcal_back_plate_thickness (double back_plate_thickness){Hcal_back_plate_thickness = back_plate_thickness ;};
  void setHcal_nlayers (int nlayer){Hcal_nlayers = nlayer;};
  void setHcal_radiator_thickness (double radiator_thickness){Hcal_radiator_thickness = radiator_thickness ;};

  void setHcal_steel_cassette_thickness (double cassette_thickness){Hcal_steel_cassette_thickness = cassette_thickness;};
  void setHcalServices_outer_FR4_thickness (double PCB_thickness){HcalServices_outer_FR4_thickness = PCB_thickness;};
  void setHcalServices_outer_Cu_thickness (double copper_thickness){HcalServices_outer_Cu_thickness = copper_thickness;};

  void setenv_safety(double safety) {
    env_safety = safety ;
    Ecal_cables_gap -= env_safety*4.0 ;
  };

  bool FillHcalServicesModuleWithInnerServices(dd4hep::PlacedVolume &pVol,
					       dd4hep::Volume &ModuleLogicalZMinus,
					       dd4hep::Volume &ModuleLogicalZPlus);




  bool PlaceHcalInnerServicesLayer(dd4hep::PlacedVolume &pVol,
				   dd4hep::Volume &motherLogical,
				   dd4hep::Material layerMaterial,
				   double layerThickness,
				   dd4hep::Position &layerPosition);



  bool FillHcalServicesModuleWithHcalElectronicsInterface(dd4hep::PlacedVolume &pVol,
							  dd4hep::Volume &ModuleLogicalZMinus,
							  dd4hep::Volume &ModuleLogicalZPlus);




  bool FillHcalElectronicsInterfaceLayer(dd4hep::PlacedVolume &pVol,
					 dd4hep::Volume &ModuleLogicalZMinus,
					 dd4hep::Volume &ModuleLogicalZPlus,
					 double layer_y_offset, double layer_x_dim);




  bool PlaceHcalElectronicsInterfaceComponent(dd4hep::PlacedVolume &pVol,
					      dd4hep::Volume &ModuleLogicalZMinus,
					      dd4hep::Volume &ModuleLogicalZPlus,
					      dd4hep::Material layerMaterial,
					      double layerThickness,
					      double y_position, double layer_x_dim){

    dd4hep::Box layerBox (layer_x_dim/2. - 4*(SurfaceTolerance),
			  layerThickness/2. - 4*(SurfaceTolerance),
			  Ecal_cables_gap/2. - 4*(SurfaceTolerance));

    dd4hep::Solid layerSolid(layerBox);

    dd4hep::Solid cutSolid = CutLayer(layerSolid,y_position);

    dd4hep::Volume layerLogical("ElectronicsInterfaceLayerLogical",cutSolid,
				layerMaterial);


    dd4hep::RotationZYX rot(0,0,M_PI*0.5);
    dd4hep::Position pos(0,0, y_position);
    dd4hep::Transform3D tran3D(rot,pos);

    pVol = ModuleLogicalZMinus.placeVolume(layerLogical,tran3D);

    pVol = ModuleLogicalZPlus.placeVolume(layerLogical,tran3D);

    return true;

  };


  dd4hep::Solid CutLayer(dd4hep::Solid &layerSolid, double y_position) ;



  bool DoBuildHcalBarrel_EndCapServices(dd4hep::PlacedVolume &pVol,
					dd4hep::Assembly &envelope) ;



private:
  double Hcal_stave_gaps;
  double Hcal_inner_radius;
  double Hcal_R_max;
  double Ecal_cables_gap;
  double TPC_Ecal_Hcal_barrel_halfZ;
//-Z thicknesses
  double ZMinus_StainlessSteel_Thickness;
  double ZMinus_PE_Thickness;
  double ZMinus_Cu_Thickness;
//+Z thicknesses
  double ZPlus_StainlessSteel_Thickness;
  double ZPlus_PE_Thickness;
  double ZPlus_Cu_Thickness;

  double InnerServicesWidth;

  double Hcal_total_dim_y;
  double Hcal_y_dim2_for_x;
  double Hcal_bottom_dim_x;
  double Hcal_midle_dim_x;
  double Hcal_top_dim_x;

  double Hcal_back_plate_thickness;
  double Hcal_radiator_thickness;
  int Hcal_nlayers;

  double Hcal_steel_cassette_thickness;
  double HcalServices_outer_FR4_thickness;
  double HcalServices_outer_Cu_thickness;

  int BuildHcalElectronicsInterface;

  const double SurfaceTolerance = 0.0001*dd4hep::mm;
  double env_safety;

  dd4hep::Material air;
  dd4hep::Material polyethylene;
  dd4hep::Material copper;
  dd4hep::Material stainless_steel;
  dd4hep::Material S235;
  dd4hep::Material PCB;
};


//=====================================
// class for BuildVXDCables
//=====================================
class BuildVXDCables {

public:
  BuildVXDCables(dd4hep::DetElement det) : detElt( det) { };
  ~BuildVXDCables(){};

  void setMaterialCopper (dd4hep::Material Cu) {copper = Cu;};
  void setVXD_cable_cross_section_area(double cable_area){VXD_cable_cross_section_area = cable_area ;};
  void setVXD_cable_z_start(double z_start){VXD_cable_z_start = z_start ;};
  void setVXD_cable_z_end(double z_end){VXD_cable_z_end = z_end ;};
  void setVXD_cable_inner1_radius(double r1_inner){VXD_cable_inner1_radius = r1_inner ;};
  void setVXD_cable_inner2_radius(double r2_inner){VXD_cable_inner2_radius = r2_inner ;};

  /// to build VXD calbe cone into the service assembly
  bool DoBuildVXDCables(dd4hep::PlacedVolume &pVol, dd4hep::Assembly &envelope) ;

private:
  double VXD_cable_cross_section_area;
  double VXD_cable_inner1_radius;
  double VXD_cable_inner2_radius;
  double VXD_cable_cone_middle_thickness;
  double VXD_cable_outer1_radius;
  double VXD_cable_outer2_radius;
  double VXD_cable_z_start;
  double VXD_cable_z_end;

  const double SurfaceTolerance = 0.0001*dd4hep::mm;

  dd4hep::Material copper;
  dd4hep::DetElement detElt;
};

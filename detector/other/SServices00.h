#include <utility>
#include <vector>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
 
using namespace std;
using namespace DD4hep;
using namespace dd4hep;
using namespace DD4hep::Geometry;

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

  void setMaterial (Material copper4cooling) {cooling_Material = copper4cooling;};
  void sethalfZ (double TPC_Ecal_Hcal_barrel_halfZ){TPC_barrel_halfZ = TPC_Ecal_Hcal_barrel_halfZ ;};
  void settpcEndplateServicesRing_R_ro( double tpcEndplateServicesRing_R, double tpcEndplateServicesRing_ro) {
    tpcEndplateServicesRing_R_ro.push_back(make_pair(tpcEndplateServicesRing_R,tpcEndplateServicesRing_ro));
  };

  // to build TPC cooling rings into the service assembly 
  bool DoBuildTPCEndplateServices(PlacedVolume *pVol,Assembly *envelope){
    for( vector< pair<double,double> >::const_iterator it = tpcEndplateServicesRing_R_ro.begin() ; 
	 it != tpcEndplateServicesRing_R_ro.end(); ++it){
      Torus solidTube(it->first, 0, it->second, 0, 2*M_PI);
      double z_position = TPC_barrel_halfZ + it->second;
      Volume EnvTube("service_tube",solidTube,cooling_Material);
      Position pos1(0,0,z_position);
      *pVol = envelope->placeVolume(EnvTube,pos1);
      Position pos2(0,0,-z_position);
      *pVol = envelope->placeVolume(EnvTube,pos2);
    } 
    return true;
  };

private:
  Material cooling_Material;
  double TPC_barrel_halfZ;
  vector< pair<double,double> > tpcEndplateServicesRing_R_ro;  
};



//=====================================
// class for  BuildEcalBarrelServices
//=====================================
class BuildEcalBarrelServices {

public:
  BuildEcalBarrelServices(){};
  ~BuildEcalBarrelServices(){};

  void setMaterial (Material air4container) {container_Material = air4container;};
  void sethalfZ (double TPC_Ecal_Hcal_barrel_halfZ){Ecal_barrel_halfZ = TPC_Ecal_Hcal_barrel_halfZ ;};
  void setTopDimX (double dim_x){top_dim_x = dim_x ;};
  void setRailHeight (double height){RailHeight = height ;};
  void setOutRadius (double outer_radius){Ecal_outer_radius = outer_radius ;};
  void setModuleThickness (double thickness){ module_thickness = thickness ;};

  // to build Ecal Barrel service into the service assembly 
  bool DoBuildEcalBarrelServices(PlacedVolume *pVol,Assembly *envelope){

    Box ContainerSolid(top_dim_x/2., RailHeight/2., Ecal_barrel_halfZ); 
    Volume containerLogical("EcalBarrelServicesContainerLogical",ContainerSolid,container_Material);

    for (int stave_id = 1; stave_id < 9 ; stave_id++)
      {
	double phirot = (stave_id-1) * M_PI/4.;
	RotationZYX rot(-phirot,0,0);
	//Position  stavePosition(module_thickness * sin(M_PI/4.), Ecal_outer_radius + RailHeight/2., 0);
	Position  stavePosition((Ecal_outer_radius + RailHeight/2.)*sin(phirot) + module_thickness * sin(M_PI/4.)*cos(-phirot), 
				(Ecal_outer_radius + RailHeight/2.)*cos(phirot)+module_thickness * sin(M_PI/4.)*sin(-phirot), 
				0 );
	Transform3D tran3D(rot,stavePosition);
	*pVol = envelope->placeVolume(containerLogical,tran3D);
      }
    
    
    return true;
  };

private:
  Material container_Material;
  double Ecal_barrel_halfZ;
  double top_dim_x;
  double RailHeight;
  double Ecal_outer_radius;
  double module_thickness;
};

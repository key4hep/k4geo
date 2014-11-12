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


//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for YokeEndcaps
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================
// *********************************************************
// *                         Mokka                         *
// *    -- A Detailed Geant 4 Simulation for the ILC --    *
// *                                                       *
// *  polywww.in2p3.fr/geant4/tesla/www/mokka/mokka.html   *
// *********************************************************
//
// $Id: YokePL04.cc,v 1.3 2008/10/31 15:41:52 frank Exp $
// $Name:  $
//
// History:  
// - first implementation P. Mora de Freitas (May 2001)
// - selectable symmetry, self-scaling, removed pole tips 
// - Adrian Vogel, 2006-03-17
// - muon system plus
//   instrumented pole tip back for TESLA models   
// - Predrag Krstonosic , 2006-08-30
// - added barrelEndcapGap, gear parameters, made barrel 
//   and endcap same thickness, made plug insensitive,  
// - F.Gaede, DESY 2008-10-04
//

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "TGeoTrd2.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

#define VERBOSE 1

static Ref_t create_detector(LCDD& lcdd, xml_h element, SensitiveDetector sens)  {
  static double tolerance = 0e0;

  xml_det_t     x_det     = element;
  string        det_name  = x_det.nameStr();
  Layering      layering (element);
 
  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();

  Material      air       = lcdd.air();
  Material      vacuum    = lcdd.vacuum();

  Material      yokeMaterial  = lcdd.material(x_det.materialStr());;

  int           det_id    = x_det.id();
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  Assembly envelope_assembly( det_name + "assembly"  ) ;  
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope_assembly);

  env_phv.addPhysVolID("system",det_id);
  sdet.setPlacement(env_phv);

  sens.setType("yoke");

 
//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build Yoke05Endcaps
//
//====================================================================
  double Yoke_barrel_inner_radius           = lcdd.constant<double>("Yoke_barrel_inner_radius");
  double Yoke_endcap_inner_radius           = lcdd.constant<double>("Yoke_endcap_inner_radius");
  //double Yoke_thickness                     = lcdd.constant<double>("Yoke_thickness");
  //double Yoke_Barrel_Half_Z                 = lcdd.constant<double>("Yoke_Barrel_Half_Z");  
  double Yoke_Z_start_endcaps               = lcdd.constant<double>("Yoke_Z_start_endcaps");


  //Database *db = new Database(env.GetDBName());
  //db->exec("SELECT * FROM `yoke`;");
  //db->getTuple();
  ////... Geometry parameters from the environment and from the database
  //symmetry = db->fetchInt("symmetry");
  //const G4double rInnerBarrel  = 
  //  env.GetParameterAsDouble("Yoke_barrel_inner_radius");
  //const G4double rInnerEndcap  = 
  //  env.GetParameterAsDouble("Yoke_endcap_inner_radius");
  //const G4double zStartEndcap  = 
  //  env.GetParameterAsDouble("Yoke_Z_start_endcaps");
  //
  //db->exec("SELECT * FROM `muon`;");
  //db->getTuple();
  //iron_thickness               = db->fetchDouble("iron_thickness");
  //G4double gap_thickness       = db->fetchDouble("layer_thickness");
  //number_of_layers             = db->fetchInt("number_of_layers");
  //G4double yokeBarrelEndcapGap = db->fetchInt("barrel_endcap_gap");
  //G4double cell_dim_x          = db->fetchDouble("cell_size");
  //G4double cell_dim_z          = db->fetchDouble("cell_size"); 
  //G4double chamber_thickness   = 10*mm;  

  double yokeBarrelEndcapGap     = 2.5;// ?? lcdd.constant<double>("barrel_endcap_gap"); //25.0*mm


//====================================================================
//
// general calculated parameters
//
//====================================================================

  //port from Mokka Yoke05, the following parameters used by Yoke05
  int    symmetry            = nsides;
  double rInnerBarrel        = Yoke_barrel_inner_radius;
  double zStartEndcap        = Yoke_Z_start_endcaps; // has been updated to 4072.0*mm by driver SCoil02 

  //TODO: put all magic numbers into ILD_o1_v05.xml file.
  double gap_thickness = 4.0;
  double iron_thickness = 10.0; //10.0 cm
  int number_of_layers = 10;

  //... Barrel parameters: 
  //... tolerance 1 mm
  double yokeBarrelThickness    = gap_thickness 
    + number_of_layers*(iron_thickness + gap_thickness) 
    + 3*(5.6*iron_thickness + gap_thickness) 
    + 0.1; // the tolerance 1 mm

  double yokeEndcapThickness    =   number_of_layers*(iron_thickness  + gap_thickness)
    + 2*(5.6*iron_thickness + gap_thickness);

  double rInnerEndcap           =    Yoke_endcap_inner_radius;
  double rOuterEndcap           =    rInnerBarrel + yokeBarrelThickness;
  double z_halfBarrel           =    zStartEndcap - yokeBarrelEndcapGap;    


  // In this release the number of modules is fixed to 3
  double Yoke_Endcap_module_dim_z =  yokeEndcapThickness+0.1; //1.0*mm;

  cout<<" Build the yoke within this dimension "<<endl;
  cout << "  ...Yoke  db: symmetry             " << symmetry <<endl;
  cout << "  ...Yoke  db: rInnerEndcap         " << rInnerEndcap <<endl;
  cout << "  ...Yoke  db: rOuterEndcap         " << rOuterEndcap <<endl;
  cout << "  ...Yoke  db: zStartEndcap         " << zStartEndcap <<endl;

  cout << "  ...Muon  db: iron_thickness       " << iron_thickness <<endl;
  cout << "  ...Muon  db: gap_thickness        " << gap_thickness <<endl;
  cout << "  ...Muon  db: number_of_layers     " << number_of_layers <<endl;

  cout << "  ...Muon par: yokeEndcapThickness  " << yokeEndcapThickness <<endl;
  cout << "  ...Muon par: Barrel_half_z        " << z_halfBarrel <<endl;




  PolyhedraRegular YokeEndcapSolid( symmetry, M_PI/symmetry, rInnerEndcap, rOuterEndcap,  Yoke_Endcap_module_dim_z);

  Volume mod_vol(det_name+"_module", YokeEndcapSolid, yokeMaterial);

  mod_vol.setVisAttributes(lcdd.visAttributes(x_det.visStr()));
     

//====================================================================
// Place Yoke05 Endcaps module into the world volume
//====================================================================

  double zEndcap          =   zStartEndcap + yokeEndcapThickness/2;
  
  for(int module_num=0;module_num<2;module_num++) {

    int module_id = ( module_num == 0 ) ? 0:6;
    double this_module_z_offset = ( module_id == 0 ) ? - zEndcap : zEndcap; 
    double this_module_rotY = ( module_id == 0 ) ? M_PI:0; 
  
    Position xyzVec(0,0,this_module_z_offset);
    RotationZYX rot(0,this_module_rotY,0);
    Rotation3D rot3D(rot);
    Transform3D tran3D(rot3D,xyzVec);

    PlacedVolume pv = envelope_assembly.placeVolume(mod_vol,tran3D);
    pv.addPhysVolID("module",module_id); // z: -/+ 0/6

    string m_name = _toString(module_id,"module%d");
    DetElement sd (m_name,det_id);
    sd.setPlacement(pv);
    sdet.add(sd);

  }

  
  return sdet;
}

DECLARE_DETELEMENT(Yoke05_Endcaps,create_detector)

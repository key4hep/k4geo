//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for YokeBarrel
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
 
  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();

  Material      air       = lcdd.air();
  Material      vacuum    = lcdd.vacuum();

  Material      yokeMaterial  = lcdd.material(x_det.materialStr());;

  int           det_id    = x_det.id();
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  Assembly envelope_assembly( det_name + "assembly"  ) ;  
  PlacedVolume pv;

  sens.setType("yoke");

  DetElement    module_det("module0",det_id);


//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build Yoke05Barrel
//
//====================================================================
  double Yoke_barrel_inner_radius           = lcdd.constant<double>("Yoke_barrel_inner_radius");
  double Yoke_thickness                     = lcdd.constant<double>("Yoke_thickness");
  double Yoke_Barrel_Half_Z                 = lcdd.constant<double>("Yoke_Barrel_Half_Z");  
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
  double Yoke_outer_radius                  = Yoke_barrel_inner_radius + Yoke_thickness;


  //port from Mokka Yoke05, the following parameters used by Yoke05
  int    symmetry            = nsides;
  double rInnerBarrel        = Yoke_barrel_inner_radius;
  double zStartEndcap        = Yoke_Z_start_endcaps; // has been updated to 4072.0*mm by driver SCoil02 
  //double yokeBarrelEndcapGap = barrel_endcap_gap;


  //... Barrel parameters: 
  //... tolerance 1 mm
  //double yokeBarrelThickness    = gap_thickness 
  //  + number_of_layers*(iron_thickness + gap_thickness) 
  //  + 3*(5.6*iron_thickness + gap_thickness) 
  //  + 1*mm;
  double yokeBarrelThickness    =    Yoke_thickness;

  double rOuterBarrel           =    rInnerBarrel + yokeBarrelThickness;    
  double z_halfBarrel           =    zStartEndcap - yokeBarrelEndcapGap;    


  // In this release the number of modules is fixed to 3
  double Yoke_Barrel_module_dim_z = 2.0*(zStartEndcap-yokeBarrelEndcapGap)/3.0 ;



  cout<<" Build the yoke within this dimension "<<endl;
  cout << "  ...Yoke  db: symmetry             " << symmetry <<endl;
  cout << "  ...Yoke  db: rInnerBarrel         " << rInnerBarrel <<endl;
  cout << "  ...Yoke  db: zStartEndcap         " << zStartEndcap <<endl;
 




  PolyhedraRegular YokeBarrelSolid( symmetry, M_PI/symmetry, rInnerBarrel, rOuterBarrel,  Yoke_Barrel_module_dim_z);

  Volume mod_vol(det_name+"_module", YokeBarrelSolid, yokeMaterial);

  mod_vol.setVisAttributes(lcdd.visAttributes("YellowVis"));
     
    
//====================================================================
// Place Yoke05 Barrel stave module into the assembly envelope
//====================================================================

  for (int module_id = 1; module_id < 4; module_id++)
    {
      double module_z_offset =  (module_id-2) * Yoke_Barrel_module_dim_z;
      
      Position pos(0,0,module_z_offset);
      
      PlacedVolume pv = envelope_assembly.placeVolume(mod_vol,pos);
      pv.addPhysVolID("module",module_id);
      DetElement sd = (module_id==0) ? module_det : module_det.clone(_toString(module_id,"module%d"));
      sd.setPlacement(pv);
      sdet.add(sd);
      
    }



//====================================================================
// Place Coil into the world volume
//====================================================================
  
  pv = motherVol.placeVolume(envelope_assembly);
  pv.addPhysVolID("system", sdet.id());
  sdet.setPlacement(pv);
  
  return sdet;
}

DECLARE_DETELEMENT(Yoke05_Barrel,create_detector)

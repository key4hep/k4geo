//====================================================================
//  DDSim - LC detector models in DD4hep 
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

  sens.setType("calorimeter");

 
//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build Yoke05Endcaps
//
//====================================================================
  double Yoke_barrel_inner_radius           = lcdd.constant<double>("Yoke_barrel_inner_radius");
  double Yoke_endcap_inner_radius           = lcdd.constant<double>("Yoke_endcap_inner_radius");
  double Yoke_Z_start_endcaps               = lcdd.constant<double>("Yoke_Z_start_endcaps");
  double HCAL_R_max                         = lcdd.constant<double>("Hcal_outer_radius");

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

  //Port from Mokka: 
  // Endcap Thickness has no tolerance
  // But the placement should shift_middle by -0.05 (0.5*mm) later
  double Yoke_Endcap_module_dim_z =  yokeEndcapThickness;

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
// Build chamber volume
//====================================================================
  //double gap_thickness       = db->fetchDouble("layer_thickness");

  //-------------------- start loop over Yoke layers ----------------------
    // Loop over the sets of layer elements in the detector.
    int l_num = 1;
    for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
      xml_comp_t x_layer = li;
      int repeat = x_layer.repeat();

      // Loop over number of repeats for this layer.
      for (int i=0; i<repeat; i++)    {
	//if(i>11) continue;
	string l_name = _toString(l_num,"layer%d");
	double l_thickness = layering.layer(i)->thickness();  // Layer's thickness.
	
	PolyhedraRegular ChamberSolid( symmetry, M_PI/symmetry, rInnerEndcap + tolerance, rOuterEndcap - tolerance,  l_thickness);
	Volume     ChamberLog(det_name+"_"+l_name,ChamberSolid,air);
	DetElement layer(l_name, det_id);

	ChamberLog.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));

	// Loop over the sublayers or slices for this layer.
	int s_num = 1;
	double s_pos_z = -(l_thickness / 2);


	
	//--------------------------------------------------------------------------------
	// Build Layer, Sensitive Scintilator in the middle, and Air tolorance at two sides 
	//--------------------------------------------------------------------------------
	for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	  xml_comp_t x_slice = si;
	  string     s_name  =  _toString(s_num,"slice%d");
	  double     s_thickness = x_slice.thickness();

	  PolyhedraRegular sliceSolid( symmetry, M_PI/symmetry, rInnerEndcap + tolerance, rOuterEndcap - tolerance,  s_thickness);
	  Volume     s_vol(det_name+"_"+l_name+"_"+s_name,sliceSolid,lcdd.material(x_slice.materialStr()));
          DetElement slice(layer,s_name,det_id);

	  if ( x_slice.isSensitive() ) {
	    s_vol.setSensitiveDetector(sens);
	  }
	  // Set region, limitset, and vis.
	  s_vol.setAttributes(lcdd,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());

	  s_pos_z += s_thickness/2.;

	  Position   s_pos(0,0,s_pos_z);      // Position of the layer.
	  PlacedVolume  s_phv = ChamberLog.placeVolume(s_vol,s_pos);
	  slice.setPlacement(s_phv);

	  if ( x_slice.isSensitive() ) {
	    s_phv.addPhysVolID("slice",s_num);
	  }

	  slice.setPlacement(s_phv);
	  // Increment x position for next slice.
	  s_pos_z += s_thickness/2.;

	  ++s_num;

	}
	
	++l_num;


	double shift_middle    = - yokeEndcapThickness/2 - 0.05 //0.5*mm 
	  + iron_thickness*(i+1) 
	  + (i+0.5)*gap_thickness; 
	
	if( i>= 10)
	  {
	    shift_middle    = - yokeEndcapThickness/2 - 0.05 //0.5*mm 
	      + iron_thickness*(i+1+(i-9)*4.6) + (i+0.5)*gap_thickness; 
	  }	

	Position xyzVec(0,0,shift_middle);
	
	PlacedVolume layer_phv =  mod_vol.placeVolume(ChamberLog,xyzVec);
	layer_phv.addPhysVolID("layer", l_num);
	//string stave_name  = "stave1";
	string stave_layer_name = "stave1"+_toString(l_num,"layer%d");
	DetElement stave(stave_layer_name,det_id);;
	stave.setPlacement(layer_phv);
	sdet.add(stave);

	
      }

    }  

//====================================================================
// Check Yoke05 plug module
//====================================================================
    bool   build_plug = false;
    double HCAL_z         = 393.7;
    double HCAL_plug_gap  = 4.5;
    double plug_thickness = zStartEndcap-HCAL_z-HCAL_plug_gap;
    
    double rInnerPlug           =    Yoke_endcap_inner_radius;
    double rOuterPlug           =    HCAL_R_max;
    double Yoke_Plug_module_dim_z =  plug_thickness;

    // Is there a space to build Yoke plug
    if( Yoke_Plug_module_dim_z > 0 ) 
      {
	build_plug = true;

	cout << "  ...Plug par: build_plug is true, there is space to build yoke plug" <<endl;
	cout << "  ...Plug par: HCAL_half_z          " << HCAL_z <<endl;
	cout << "  ...Plug par: HCAL_Plug_Gap        " << HCAL_plug_gap <<endl;
	cout << "  ...Plug par: Plug Thickness       " << plug_thickness <<endl;
	cout << "  ...Plug par: Plug Radius          " << HCAL_R_max <<endl;

      }

//====================================================================
// Place Yoke05 Endcaps module into the world volume
//====================================================================

  double zEndcap          =   zStartEndcap + yokeEndcapThickness/2.0 + 0.1; // Need 0.1 (1.0*mm) according to the Mokka Yoke05 driver.
  double zPlug            =   zStartEndcap - plug_thickness/2.0 -0.05; //  Need 0.05 (0.5*mm) according to the Mokka Yoke05 driver.
  
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

    //====================================================================
    // If build_plug is true, Place the plug module into the world volume
    //====================================================================
    if(build_plug == true){
      PolyhedraRegular YokePlugSolid( symmetry, M_PI/symmetry, rInnerPlug, rOuterPlug,  Yoke_Plug_module_dim_z);
      Volume plug_vol(det_name+"_plug", YokePlugSolid, yokeMaterial);
      plug_vol.setVisAttributes(lcdd.visAttributes(x_det.visStr()));

      double this_plug_z_offset = ( module_id == 0 ) ? - zPlug : zPlug; 
      Position   plug_pos(0,0,this_plug_z_offset);
      PlacedVolume  plug_phv = envelope_assembly.placeVolume(plug_vol,plug_pos);
      string plug_name = _toString(module_id,"plug%d");
      DetElement plug (plug_name,det_id);
      plug.setPlacement(plug_phv);
    }

  }

  
  return sdet;
}

DECLARE_DETELEMENT(Yoke05_Endcaps,create_detector)

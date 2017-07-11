//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for YokeEndcaps
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id$
//====================================================================
// *********************************************************
// *                         Mokka                         *
// *    -- A Detailed Geant 4 Simulation for the ILC --    *
// *                                                       *
// *  polywww.in2p3.fr/geant4/tesla/www/mokka/mokka.html   *
// *********************************************************
//
// $Id$
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
#include "DD4hep/DetType.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::DetType;
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

#define VERBOSE 1

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {
  static double tolerance = 0e0;

  xml_det_t     x_det     = element;
  string        det_name  = x_det.nameStr();
  Layering      layering (element);
 
  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();

  Material      air       = theDetector.air();
  //unused: Material      vacuum    = theDetector.vacuum();

  Material      yokeMaterial  = theDetector.material(x_det.materialStr());;

  int           det_id    = x_det.id();
  DetElement    sdet      (det_name,det_id);

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  sens.setType("calorimeter");

 
//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build Yoke05Endcaps
//
//====================================================================
  double Yoke_barrel_inner_radius           = theDetector.constant<double>("Yoke_barrel_inner_radius");
  double Yoke_endcap_inner_radius           = theDetector.constant<double>("Yoke_endcap_inner_radius");
  double Yoke_Z_start_endcaps               = theDetector.constant<double>("Yoke_Z_start_endcaps");
  double HCAL_R_max                         = theDetector.constant<double>("Hcal_outer_radius");
  //double Yoke_cells_size                    = theDetector.constant<double>("Yoke_cells_size");

  double yokeBarrelEndcapGap     = 2.5;// ?? theDetector.constant<double>("barrel_endcap_gap"); //25.0*mm


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
    + 2*(5.6*iron_thickness + gap_thickness) + gap_thickness;

  double rInnerEndcap           =    Yoke_endcap_inner_radius;
  double rOuterEndcap           =    rInnerBarrel + yokeBarrelThickness;
  double z_halfBarrel           =    zStartEndcap - yokeBarrelEndcapGap;    

  //Port from Mokka: 
  // Endcap Thickness has no tolerance
  // But the placement should shift_middle by -0.05 (0.5*mm) later
  double Yoke_Endcap_module_dim_z =  yokeEndcapThickness;

  //double Yoke_cell_dim_x        = rOuterEndcap*2.0 / floor (rOuterEndcap*2.0/Yoke_cells_size);
  //double Yoke_cell_dim_y        = Yoke_cell_dim_x;

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

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];
  
  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = symmetry  ;
  caloData->outer_symmetry = symmetry  ;
  caloData->phi0 = 0 ; // hardcoded

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rInnerEndcap ;
  caloData->extent[1] = rOuterEndcap ;
  caloData->extent[2] = zStartEndcap ;
  caloData->extent[3] = zStartEndcap + Yoke_Endcap_module_dim_z ;



  PolyhedraRegular YokeEndcapSolid( symmetry, M_PI/symmetry, rInnerEndcap, rOuterEndcap,  Yoke_Endcap_module_dim_z);

  Volume mod_vol(det_name+"_module", YokeEndcapSolid, yokeMaterial);

  mod_vol.setVisAttributes(theDetector.visAttributes(x_det.visStr()));
     

//====================================================================
// Build chamber volume
//====================================================================
  //double gap_thickness       = db->fetchDouble("layer_thickness");

  //-------------------- start loop over Yoke layers ----------------------
  // Loop over the sets of layer elements in the detector.

  double nRadiationLengths=0.;
  double nInteractionLengths=0.;
  double thickness_sum=0;


    int l_num = 1;
    for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
      xml_comp_t x_layer = li;
      int repeat = x_layer.repeat();

      // Loop over number of repeats for this layer.
      for (int i=0; i<repeat; i++)    {
	//if(i>11) continue;
	string l_name = _toString(l_num,"layer%d");
	double l_thickness = layering.layer(i)->thickness();  // Layer's thickness.

	LayeredCalorimeterData::Layer caloLayer ;
	caloLayer.cellSize0 = cell_sizeX;
	caloLayer.cellSize1 = cell_sizeY;
	
	PolyhedraRegular ChamberSolid( symmetry, M_PI/symmetry, rInnerEndcap + tolerance, rOuterEndcap - tolerance,  l_thickness);
	Volume     ChamberLog(det_name+"_"+l_name,ChamberSolid,air);
	DetElement layer(l_name, det_id);

	ChamberLog.setVisAttributes(theDetector.visAttributes(x_layer.visStr()));

	// Loop over the sublayers or slices for this layer.
	int s_num = 1;
	double s_pos_z = -(l_thickness / 2);


	
	//--------------------------------------------------------------------------------
	// Build Layer, Sensitive Scintilator in the middle, and Air tolorance at two sides 
	//--------------------------------------------------------------------------------
	double radiator_thickness = 0.05 + 0.5*gap_thickness + iron_thickness - l_thickness/2.0 ;
	if ( i>0 )   radiator_thickness = gap_thickness + iron_thickness - l_thickness ;
	if ( i>=10 ) radiator_thickness = gap_thickness + 5.6*iron_thickness - l_thickness ;

	nRadiationLengths   = radiator_thickness/(yokeMaterial.radLength());
	nInteractionLengths = radiator_thickness/(yokeMaterial.intLength());
	thickness_sum       = radiator_thickness;
	
	for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	  xml_comp_t x_slice = si;
	  string     s_name  =  _toString(s_num,"slice%d");
	  double     s_thickness = x_slice.thickness();
	  Material slice_material  = theDetector.material(x_slice.materialStr());
	  
	  PolyhedraRegular sliceSolid( symmetry, M_PI/symmetry, rInnerEndcap + tolerance, rOuterEndcap - tolerance,  s_thickness);
	  Volume     s_vol(det_name+"_"+l_name+"_"+s_name,sliceSolid,slice_material);
          DetElement slice(layer,s_name,det_id);

	  nRadiationLengths   += s_thickness/(2.*slice_material.radLength());
	  nInteractionLengths += s_thickness/(2.*slice_material.intLength());
	  thickness_sum       += s_thickness/2;

	  if ( x_slice.isSensitive() ) {
	    s_vol.setSensitiveDetector(sens);

#if DD4HEP_VERSION_GE( 0, 15 )
	  //Store "inner" quantities
	  caloLayer.inner_nRadiationLengths   = nRadiationLengths;
	  caloLayer.inner_nInteractionLengths = nInteractionLengths;
	  caloLayer.inner_thickness           = thickness_sum;
	  //Store scintillator thickness
	  caloLayer.sensitive_thickness       = s_thickness;
#endif
	  //Reset counters to measure "outside" quantitites
	  nRadiationLengths=0.;
	  nInteractionLengths=0.;
	  thickness_sum = 0.;
	  }

	  nRadiationLengths   += s_thickness/(2.*slice_material.radLength());
	  nInteractionLengths += s_thickness/(2.*slice_material.intLength());
	  thickness_sum       += s_thickness/2;
	  
	  // Set region, limitset, and vis.
	  s_vol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());

	  s_pos_z += s_thickness/2.;

	  Position   s_pos(0,0,s_pos_z);      // Position of the layer.
	  PlacedVolume  s_phv = ChamberLog.placeVolume(s_vol,s_pos);
	  slice.setPlacement(s_phv);

	  // Increment x position for next slice.
	  s_pos_z += s_thickness/2.;

	  ++s_num;

	}

#if DD4HEP_VERSION_GE( 0, 15 )
	//Store "outer" quantities
	caloLayer.outer_nRadiationLengths   = nRadiationLengths;
	caloLayer.outer_nInteractionLengths = nInteractionLengths;
	caloLayer.outer_thickness           = thickness_sum;
#endif
	
	++l_num;


	double shift_middle    = - yokeEndcapThickness/2 + 0.05 //0.5*mm 
	  + iron_thickness*(i+1) 
	  + (i+0.5)*gap_thickness; 
	
	if( i>= 10)
	  {
	    shift_middle    = - yokeEndcapThickness/2 + 0.05 //0.5*mm 
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

      //-----------------------------------------------------------------------------------------
	

	caloLayer.distance = zStartEndcap + yokeEndcapThickness/2.0 + shift_middle
	  - caloLayer.inner_thickness ;
	caloLayer.absorberThickness = radiator_thickness ;
	
	caloData->layers.push_back( caloLayer ) ;

      //-----------------------------------------------------------------------------------------
	
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

    PlacedVolume pv = envelope.placeVolume(mod_vol,tran3D);
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
      plug_vol.setVisAttributes(theDetector.visAttributes(x_det.visStr()));

      double this_plug_z_offset = ( module_id == 0 ) ? - zPlug : zPlug; 
      Position   plug_pos(0,0,this_plug_z_offset);
      PlacedVolume  plug_phv = envelope.placeVolume(plug_vol,plug_pos);
      string plug_name = _toString(module_id,"plug%d");
      DetElement plug (plug_name,det_id);
      plug.setPlacement(plug_phv);
    }

  }

  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
  
  return sdet;
}

DECLARE_DETELEMENT(Yoke05_Endcaps,create_detector)

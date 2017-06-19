//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for YokeBarrel
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
#include "TGeoTrd2.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
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

  xml_comp_t    x_staves  = x_det.staves();
  Material      yokeMaterial  = theDetector.material(x_staves.materialStr());

  //unused: Material      env_mat     = theDetector.material(x_dim.materialStr());

  xml_comp_t    env_pos     = x_det.position();
  xml_comp_t    env_rot     = x_det.rotation();

  Position      pos(env_pos.x(),env_pos.y(),env_pos.z());
  RotationZYX   rotZYX(env_rot.z(),env_rot.y(),env_rot.x());

  Transform3D   tr(rotZYX,pos);

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
// Use them to build Yoke05Barrel
//
//====================================================================
  double Yoke_barrel_inner_radius           = theDetector.constant<double>("Yoke_barrel_inner_radius");
  //double Yoke_thickness                     = theDetector.constant<double>("Yoke_thickness");
  //double Yoke_Barrel_Half_Z                 = theDetector.constant<double>("Yoke_Barrel_Half_Z");  
  double Yoke_Z_start_endcaps               = theDetector.constant<double>("Yoke_Z_start_endcaps");
  //double Yoke_cells_size                    = theDetector.constant<double>("Yoke_cells_size");

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

  double rOuterBarrel           =    rInnerBarrel + yokeBarrelThickness;    
  double z_halfBarrel           =    zStartEndcap - yokeBarrelEndcapGap;    


  // In this release the number of modules is fixed to 3
  double Yoke_Barrel_module_dim_z = 2.0*(zStartEndcap-yokeBarrelEndcapGap)/3.0 ;

  //double Yoke_cell_dim_x        = Yoke_cells_size;
  //double Yoke_cell_dim_z        = Yoke_Barrel_module_dim_z / floor (Yoke_Barrel_module_dim_z/Yoke_cell_dim_x);

  cout<<" Build the yoke within this dimension "<<endl;
  cout << "  ...Yoke  db: symmetry             " << symmetry <<endl;
  cout << "  ...Yoke  db: rInnerBarrel         " << rInnerBarrel <<endl;
  cout << "  ...Yoke  db: zStartEndcap         " << zStartEndcap <<endl;

  cout << "  ...Muon  db: iron_thickness       " << iron_thickness <<endl;
  cout << "  ...Muon  db: gap_thickness        " << gap_thickness <<endl;
  cout << "  ...Muon  db: number_of_layers     " << number_of_layers <<endl;

  cout << "  ...Muon par: yokeBarrelThickness  " << yokeBarrelThickness <<endl;
  cout << "  ...Muon par: Barrel_half_z        " << z_halfBarrel <<endl;

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];
  
  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::BarrelLayout ;
  caloData->inner_symmetry = symmetry  ;
  caloData->outer_symmetry = symmetry  ;
  caloData->phi0 = 0 ; // also hardcoded below

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rInnerBarrel ;
  caloData->extent[1] = rOuterBarrel ;
  caloData->extent[2] = 0. ;
  caloData->extent[3] = z_halfBarrel ;


// ========= Create Yoke Barrel module   ====================================
  PolyhedraRegular YokeBarrelSolid( symmetry, M_PI/2.0-M_PI/symmetry, rInnerBarrel, rOuterBarrel,  Yoke_Barrel_module_dim_z);

  Volume mod_vol(det_name+"_module", YokeBarrelSolid, yokeMaterial);

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
	//double l_thickness = layering.layer(l_num-1)->thickness();  // Layer's thickness.
	double l_thickness = layering.layer(i)->thickness();  // Layer's thickness.
	
	//double gap_thickness = l_thickness;
	//double iron_thickness = 10.0; //10.0 cm

	double radius_low = rInnerBarrel+ 0.05 + i*gap_thickness + i*iron_thickness; 
	//rInnerBarrel+ 0.5*mm + i*gap_thickness + i*iron_thickness; 
	//double radius_mid      = radius_low+0.5*gap_thickness;  
	//double radius_sensitive = radius_mid;

	if( i>=10 ) radius_low =  rInnerBarrel + 0.05 + i*gap_thickness  + (i+(i-10)*4.6)*iron_thickness;
	//{ radius_low =  
	//    rInnerBarrel + 0.5*mm + i*gap_thickness 
	//    + (i+(i-10)*4.6)*iron_thickness;
	//radius_mid       = radius_low+0.5*gap_thickness;  
	//radius_sensitive = radius_mid;
	//}
	
	//... safety margines of 0.1 mm for x,y of chambers
	//double dx = radius_low*tan(Angle2)-0.1*mm;
	//double dy = (zStartEndcap-yokeBarrelEndcapGap)/3.0-0.1*mm; 

	double Angle2 = M_PI/symmetry;
	double dx = radius_low*tan(Angle2)-0.01;
	double dy = (zStartEndcap-yokeBarrelEndcapGap)/3.0-0.01; 
	//Box ChamberSolid(dx,gap_thickness/2.,dy);
	//Volume ChamberLog("muonSci",ChamberSolid,air);

	LayeredCalorimeterData::Layer caloLayer ;
	caloLayer.cellSize0 = cell_sizeX;
	caloLayer.cellSize1 = cell_sizeY;
      
	Box        ChamberSolid(dx,l_thickness/2.0, dy);
	Volume     ChamberLog(det_name+"_"+l_name,ChamberSolid,air);
	DetElement layer(l_name, det_id);

	ChamberLog.setVisAttributes(theDetector.visAttributes(x_layer.visStr()));

	// Loop over the sublayers or slices for this layer.
	int s_num = 1;
	double s_pos_y = -(l_thickness / 2);


	
	//--------------------------------------------------------------------------------
	// Build Layer, Sensitive Scintilator in the middle, and Air tolorance at two sides 
	//--------------------------------------------------------------------------------
	double radiator_thickness = 0.05; // Yoke05 Barrel: No radiator before first sensitive layer.
      	if ( i>0 )   radiator_thickness = gap_thickness + iron_thickness - l_thickness;
	if ( i>=10 ) radiator_thickness = gap_thickness + 5.6*iron_thickness - l_thickness;

	nRadiationLengths   = radiator_thickness/(yokeMaterial.radLength());
	nInteractionLengths = radiator_thickness/(yokeMaterial.intLength());
	thickness_sum       = radiator_thickness;


	for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	  xml_comp_t x_slice = si;
	  string     s_name  =  _toString(s_num,"slice%d");
	  double     s_thickness = x_slice.thickness();
	  Material slice_material  = theDetector.material(x_slice.materialStr());

	  double slab_dim_x = dx-tolerance;
	  double slab_dim_y = s_thickness/2.;
	  double slab_dim_z = dy-tolerance;

	  Box        s_box(slab_dim_x,slab_dim_y,slab_dim_z);
	  Volume     s_vol(det_name+"_"+l_name+"_"+s_name,s_box,slice_material);
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

	  s_pos_y += s_thickness/2.;

	  Position   s_pos(0,s_pos_y,0);      // Position of the layer.
	  PlacedVolume  s_phv = ChamberLog.placeVolume(s_vol,s_pos);
	  slice.setPlacement(s_phv);

	  // Increment x position for next slice.
	  s_pos_y += s_thickness/2.;

	  ++s_num;

	}

#if DD4HEP_VERSION_GE( 0, 15 )
	//Store "outer" quantities
	caloLayer.outer_nRadiationLengths   = nRadiationLengths;
	caloLayer.outer_nInteractionLengths = nInteractionLengths;
	caloLayer.outer_thickness           = thickness_sum;
#endif
      
	++l_num;


	
	double phirot = 0;

	for(int j=0;j<symmetry;j++)
	  {
	    double Y = radius_low + l_thickness/2.0;
	    Position xyzVec(-Y*sin(phirot), Y*cos(phirot), 0);

	    RotationZYX rot(phirot,0,0);
	    Rotation3D rot3D(rot);

	    Transform3D tran3D(rot3D,xyzVec); 
	    PlacedVolume layer_phv =  mod_vol.placeVolume(ChamberLog,tran3D);
	    layer_phv.addPhysVolID("layer", l_num).addPhysVolID("stave",j+1);
	    string     stave_name  =  _toString(j+1,"stave%d");
	    string stave_layer_name = stave_name+_toString(l_num,"layer%d");
	    DetElement stave(stave_layer_name,det_id);;
	    stave.setPlacement(layer_phv);
	    sdet.add(stave);
	    phirot -= M_PI/symmetry*2.0;

	  }

	//-----------------------------------------------------------------------------------------


	caloLayer.distance = radius_low - radiator_thickness ;
	caloLayer.absorberThickness = radiator_thickness ;
	
	caloData->layers.push_back( caloLayer ) ;

	//-----------------------------------------------------------------------------------------
	
      }

    }  


//====================================================================
// Place Yoke05 Barrel stave module into the world volume
//====================================================================

  for (int module_id = 1; module_id < 4; module_id++)
    {
      double module_z_offset =  (module_id-2) * Yoke_Barrel_module_dim_z;
      
      Position mpos(0,0,module_z_offset);
      
      PlacedVolume m_phv = envelope.placeVolume(mod_vol,mpos);
      m_phv.addPhysVolID("module",module_id).addPhysVolID("system", det_id);
      m_phv.addPhysVolID("tower", 1);// Not used
      string m_name = _toString(module_id,"module%d");
      DetElement sd (m_name,det_id);
      sd.setPlacement(m_phv);
      sdet.add(sd);
      
    }

  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
  
  return sdet;
}

DECLARE_DETELEMENT(Yoke05_Barrel,create_detector)

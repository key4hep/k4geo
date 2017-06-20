//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for SiWEcalBarrel
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id: SEcal04_Barrel_v01.cpp 876 2016-01-26 15:20:44Z shaojun $
//====================================================================

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
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Trapezoid;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

//#define VERBOSE 1

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

/** SEcal04.cc
 *
 *  @author: Shaojun Lu, DESY
 *  @version $Id: SEcal04_Barrel.cpp 876 2016-01-26 15:20:44Z shaojun $
 *              Ported from Mokka SEcal04 Barrel part. Read the constants from XML
 *              instead of the DB. Then build the Barrel in the same way with DD4hep
 * 		construct.
 *
 * @history: F.Gaede, CERN/DESY, Nov. 10, 2014
 *              added information for reconstruction: LayeringExtension and surfaces (experimental)
 *              removed DetElement for slices (not needed) increased multiplicity for layer DetElement
 *              along tower index  
 *   F.Gaede: 03/2015: 
 *            create the envelope volume with create_placed_envelope() using the xml 
 *   S.Lu:    02/2016:
 *            simplify the geometry by removing the Si wafer structure.
 *            Todo: thinking about a virtual cell and gap in Si layer.
 */

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {
  static double tolerance = 0e0;

  xml_det_t     x_det     = element;
  string        det_name  = x_det.nameStr();
  Layering      layering (element);

  Material      air       = theDetector.air();
  //unused: Material      vacuum    = theDetector.vacuum();

  int           det_id    = x_det.id();
  xml_comp_t    x_staves  = x_det.staves();
  DetElement    sdet      (det_name,det_id);

  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();
  double        dphi      = (2*M_PI/nsides);
  double        hphi      = dphi/2;


 // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;

  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------


  sens.setType("calorimeter");

  Material stave_material  = theDetector.material(x_staves.materialStr());

  DetElement    stave_det("module0stave0",det_id);

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];

//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build HcalBarrel
//
//====================================================================

  int N_FIBERS_W_STRUCTURE = 2; 
  int N_FIBERS_ALVOULUS = 3;

  //  read parametere from compact.xml file
  double Ecal_Alveolus_Air_Gap              = theDetector.constant<double>("Ecal_Alveolus_Air_Gap");
  double Ecal_Slab_shielding                = theDetector.constant<double>("Ecal_Slab_shielding");
  double Ecal_Slab_copper_thickness         = theDetector.constant<double>("Ecal_Slab_copper_thickness");
  double Ecal_Slab_PCB_thickness            = theDetector.constant<double>("Ecal_Slab_PCB_thickness");
  double Ecal_Slab_glue_gap                 = theDetector.constant<double>("Ecal_Slab_glue_gap");
  double Ecal_Slab_ground_thickness         = theDetector.constant<double>("Ecal_Slab_ground_thickness");
  double Ecal_fiber_thickness               = theDetector.constant<double>("Ecal_fiber_thickness");
  double Ecal_Si_thickness                  = theDetector.constant<double>("Ecal_Si_thickness");
  
  double Ecal_inner_radius                  = theDetector.constant<double>("TPC_outer_radius") +theDetector.constant<double>("Ecal_Tpc_gap");
  double Ecal_radiator_thickness1           = theDetector.constant<double>("Ecal_radiator_layers_set1_thickness");
  double Ecal_radiator_thickness2           = theDetector.constant<double>("Ecal_radiator_layers_set2_thickness");
  double Ecal_radiator_thickness3           = theDetector.constant<double>("Ecal_radiator_layers_set3_thickness");
  double Ecal_Barrel_halfZ                  = theDetector.constant<double>("Ecal_Barrel_halfZ");
  
  double Ecal_support_thickness             = theDetector.constant<double>("Ecal_support_thickness");
  double Ecal_front_face_thickness          = theDetector.constant<double>("Ecal_front_face_thickness");
  double Ecal_lateral_face_thickness        = theDetector.constant<double>("Ecal_lateral_face_thickness");
  double Ecal_Slab_H_fiber_thickness        = theDetector.constant<double>("Ecal_Slab_H_fiber_thickness");

  double Ecal_Slab_Sc_PCB_thickness         = theDetector.constant<double>("Ecal_Slab_Sc_PCB_thickness");
  double Ecal_Sc_thickness                  = theDetector.constant<double>("Ecal_Sc_thickness");
  double Ecal_Sc_reflector_thickness        = theDetector.constant<double>("Ecal_Sc_reflector_thickness");

  int    Ecal_nlayers1                      = theDetector.constant<int>("Ecal_nlayers1");
  int    Ecal_nlayers2                      = theDetector.constant<int>("Ecal_nlayers2");
  int    Ecal_nlayers3                      = theDetector.constant<int>("Ecal_nlayers3");
  int    Ecal_barrel_number_of_towers       = theDetector.constant<int>("Ecal_barrel_number_of_towers");
  
  //double Ecal_guard_ring_size               = theDetector.constant<double>("Ecal_guard_ring_size");
  
//====================================================================
//
// general calculated parameters
//
//====================================================================
  
  double Ecal_total_SiSlab_thickness = 
    Ecal_Slab_shielding + 
    Ecal_Slab_copper_thickness + 
    Ecal_Slab_PCB_thickness +
    Ecal_Slab_glue_gap + 
    Ecal_Si_thickness + 
    Ecal_Slab_ground_thickness +
    Ecal_Alveolus_Air_Gap / 2;
#ifdef VERBOSE
  std::cout << " Ecal_total_SiSlab_thickness = " << Ecal_total_SiSlab_thickness  << std::endl;
#endif
  
  

  double Ecal_total_ScSlab_thickness = 
    Ecal_Slab_shielding + 
    Ecal_Slab_copper_thickness + 
    Ecal_Slab_Sc_PCB_thickness +
    Ecal_Sc_thickness + 
    Ecal_Sc_reflector_thickness * 2 +
    Ecal_Alveolus_Air_Gap / 2;
#ifdef VERBOSE
  std::cout << " Ecal_total_ScSlab_thickness = " << Ecal_total_ScSlab_thickness  << std::endl;
#endif
  

    int Number_of_Si_Layers_in_Barrel = 0;
    int Number_of_Sc_Layers_in_Barrel = 0;


#ifdef VERBOSE
  std::cout << " Ecal total number of Silicon layers = " << Number_of_Si_Layers_in_Barrel  << std::endl;
  std::cout << " Ecal total number of Scintillator layers = " << Number_of_Sc_Layers_in_Barrel  << std::endl;
#endif
  
  // In this release the number of modules is fixed to 5
  double Ecal_Barrel_module_dim_z = 2 * Ecal_Barrel_halfZ / 5. ;
#ifdef VERBOSE
  std::cout << "Ecal_Barrel_module_dim_z  = " << Ecal_Barrel_module_dim_z  << std::endl;
#endif

  // The alveolus size takes in account the module Z size
  // but also 4 fiber layers for the alveoulus wall, the all
  // divided by the number of towers
  double alveolus_dim_z = 
    (Ecal_Barrel_module_dim_z - 2. * Ecal_lateral_face_thickness) /
    Ecal_barrel_number_of_towers - 
    2 * N_FIBERS_ALVOULUS  * Ecal_fiber_thickness  - 
    2 * Ecal_Slab_H_fiber_thickness -
    2 * Ecal_Slab_shielding;


#ifdef VERBOSE
  std::cout << "alveolus_dim_z = " <<  alveolus_dim_z << std::endl;
#endif

  int n_total_layers = Ecal_nlayers1 + Ecal_nlayers2 + Ecal_nlayers3;
  Number_of_Si_Layers_in_Barrel = n_total_layers+1;

  double module_thickness = 
    Ecal_nlayers1 * Ecal_radiator_thickness1 +
    Ecal_nlayers2 * Ecal_radiator_thickness2 +
    Ecal_nlayers3 * Ecal_radiator_thickness3 +
    
    int(n_total_layers/2) * // fiber around W struct layers
    (N_FIBERS_W_STRUCTURE * 2 *  Ecal_fiber_thickness) +
    
    Number_of_Si_Layers_in_Barrel * // Silicon slabs plus fiber around and inside
    (Ecal_total_SiSlab_thickness +
     (N_FIBERS_ALVOULUS + 1 ) * Ecal_fiber_thickness) +
    
    Number_of_Sc_Layers_in_Barrel * // Scintillator slabs plus fiber around and inside
    (Ecal_total_ScSlab_thickness +
     (N_FIBERS_ALVOULUS + 1 ) * Ecal_fiber_thickness) +
    
    Ecal_support_thickness + Ecal_front_face_thickness;
  
#ifdef VERBOSE
  std::cout << "For information : module_thickness = " << module_thickness  << std::endl;
#endif
  
  // module barrel key parameters
  double  bottom_dim_x = 2. * tan(M_PI/8.) * Ecal_inner_radius +
    module_thickness/sin(M_PI/4.);
  
  double top_dim_x = bottom_dim_x - 2 * module_thickness;

  //------------------------------------------------------------------------------------

  LayeredCalorimeterData::Layer caloLayer ;
  caloLayer.cellSize0 = cell_sizeX;
  caloLayer.cellSize1 = cell_sizeY;

  //== For Wafer ===  
  //double cell_dim_x = caloLayer.cellSize0;
  //double total_Si_dim_z = alveolus_dim_z;

  //double util_SI_wafer_dim_z = 
  //  total_Si_dim_z/2 -  2 * Ecal_guard_ring_size;

  //double cell_dim_z =  util_SI_wafer_dim_z/ 
  //  floor(util_SI_wafer_dim_z/
  //	  cell_dim_x);

  //int N_cells_in_Z = int(util_SI_wafer_dim_z/cell_dim_z);
  //int N_cells_in_X = N_cells_in_Z;
  
  //cell_dim_x = cell_dim_z;


  
#ifdef VERBOSE
  std::cout << " bottom_dim_x = " << bottom_dim_x  << std::endl;
  std::cout << " top_dim_x = " << top_dim_x << std::endl;
  std::cout << " Ecal total number of Silicon layers = " << Number_of_Si_Layers_in_Barrel  << std::endl;
  std::cout << " Ecal total number of Scintillator layers = " << Number_of_Sc_Layers_in_Barrel  << std::endl;
#endif
  


// ========= Create Ecal Barrel stave   ====================================
//  It will be the volume for palcing the Ecal Barrel alveolus(i.e. Layers).
//  And the structure W plate.
//  Itself will be placed into the world volume.
// ==========================================================================

  // The TOP_X and BOTTOM_X is different in Mokka and DD4hep
  Trapezoid trd(top_dim_x / 2,
		bottom_dim_x / 2, 
		Ecal_Barrel_module_dim_z / 2,
		Ecal_Barrel_module_dim_z / 2,
		module_thickness/2);

  Volume mod_vol(det_name+"_module",trd,air);



  // We count the layers starting from IP and from 1,
  // so odd layers should be inside slabs and
  // even ones on the structure.
  // The structure W layers are here big plans, as the 
  // gap between each W plate is too small to create problems 
  // The even W layers are part of H structure placed inside
  // the alveolus.

  // ############################
  //  Dimension of radiator wLog
  //  slice provide the thickness
  // ############################


  double y_floor = 
    Ecal_front_face_thickness +
    N_FIBERS_ALVOULUS * Ecal_fiber_thickness;


  // ############################
  //  Dimension of alveolus
  //  slice provide the thickness
  // ############################
  
    // =====  build Si Slab and put into the Layer volume =====
    // =====  place the layer into the module 5 time for one full layer into the trd module ====
    // =====  build and place barrel structure into trd module ====
    // Parameters for computing the layer X dimension:
    double stave_z  =(Ecal_Barrel_module_dim_z - 2. * Ecal_lateral_face_thickness) / Ecal_barrel_number_of_towers/2.;
    double l_dim_x  = bottom_dim_x/2.;                            // Starting X dimension for the layer.
    double l_pos_z  = module_thickness/2;

    l_dim_x -= y_floor;
    l_pos_z -= y_floor;


    // ------------- create extension objects for reconstruction -----------------
    
    //========== fill data for reconstruction ============================
    LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
    caloData->layoutType = LayeredCalorimeterData::BarrelLayout ;
    caloData->inner_symmetry = nsides  ;
    //added by Thorben Quast
    caloData->outer_symmetry = nsides  ;
    caloData->phi0 = 0 ; // hardcoded 
    
    /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
    caloData->extent[0] = Ecal_inner_radius ;
    //line fixed by Thorben Quast since actual conversion is made during the drawing
    caloData->extent[1] = ( Ecal_inner_radius + module_thickness );
    //caloData->extent[1] = ( Ecal_inner_radius + module_thickness ) / cos( M_PI/8. ) ;
    caloData->extent[2] = 0. ;
    caloData->extent[3] = Ecal_Barrel_halfZ ;

    // // base vectors for surfaces:
    // dd4hep::rec::Vector3D u(1,0,0) ;
    // dd4hep::rec::Vector3D v(0,1,0) ;
    // dd4hep::rec::Vector3D n(0,0,1) ;


    //-------------------- start loop over ECAL layers ----------------------
    // Loop over the sets of layer elements in the detector.

    double nRadiationLengths   = 0. ;
    double nInteractionLengths = 0. ;
    double thickness_sum       = 0. ;

    nRadiationLengths   = Ecal_radiator_thickness1/(stave_material.radLength())
      + y_floor/air.radLength();
    nInteractionLengths = Ecal_radiator_thickness1/(stave_material.intLength())
      + y_floor/air.intLength();
    thickness_sum       = Ecal_radiator_thickness1 + y_floor;

    int l_num = 1;
    bool isFirstSens = true;
    int myLayerNum = 0 ;

    for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
      xml_comp_t x_layer = li;
      int repeat = x_layer.repeat();
      // Loop over number of repeats for this layer.
      for (int j=0; j<repeat; j++)    {
	string l_name = _toString(l_num,"layer%d");
	double l_thickness = layering.layer(l_num-1)->thickness();  // Layer's thickness.
	double xcut = (l_thickness);                     // X dimension for this layer.
	l_dim_x -= xcut;


	Box        l_box(l_dim_x-tolerance,stave_z-tolerance,l_thickness/2.0-tolerance);
	Volume     l_vol(det_name+"_"+l_name,l_box,air);

	l_vol.setVisAttributes(theDetector.visAttributes(x_layer.visStr()));

	//fg: need vector of DetElements for towers ! 
	//    DetElement layer(stave_det, l_name, det_id);
	std::vector< DetElement > layers( Ecal_barrel_number_of_towers )  ;
	
	// place layer 5 times in module. at same layer position (towers !)
	double l_pos_y = Ecal_Barrel_module_dim_z / 2. 
	  - ( Ecal_lateral_face_thickness +
	      Ecal_fiber_thickness * N_FIBERS_ALVOULUS +
	      Ecal_Slab_shielding + 
	      Ecal_Slab_H_fiber_thickness +
	      alveolus_dim_z /2.);						  
       	for (int i=0; i<Ecal_barrel_number_of_towers; i++){ // need four clone

	  layers[i] = DetElement( stave_det, l_name+_toString(i,"tower%02d") , det_id ) ;

	  Position   l_pos(0,l_pos_y,l_pos_z-l_thickness/2.);      // Position of the layer.
	  PlacedVolume layer_phv = mod_vol.placeVolume(l_vol,l_pos);
	  // layer_phv.addPhysVolID("layer", l_num);
	  layer_phv.addPhysVolID("tower", i);

	  layers[i].setPlacement(layer_phv);
	  l_pos_y -= (alveolus_dim_z + 
		      2. * Ecal_fiber_thickness * N_FIBERS_ALVOULUS +
		      2. * Ecal_Slab_H_fiber_thickness +
		      2. * Ecal_Slab_shielding);
	}




	// Loop over the sublayers or slices for this layer.
	int s_num = 1;
	double s_pos_z = l_thickness / 2.;


	//--------------------------------------------------------------------------------
	// BuildBarrelAlveolus: BuildSiliconSlab:
	//--------------------------------------------------------------------------------
	double radiator_dim_y = Ecal_radiator_thickness1; //to be updated with slice radiator thickness 

	for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	  xml_comp_t x_slice = si;
	  string     s_name  =  _toString(s_num,"slice%d");
	  double     s_thick = x_slice.thickness();
	  Material slice_material  = theDetector.material(x_slice.materialStr());
#ifdef VERBOSE
	  std::cout<<"Ecal_barrel_number_of_towers: "<< Ecal_barrel_number_of_towers <<std::endl;
#endif
	  double slab_dim_x = l_dim_x-tolerance;
	  double slab_dim_y = s_thick/2.;
	  double slab_dim_z = stave_z-tolerance;

	  Box        s_box(slab_dim_x,slab_dim_z,slab_dim_y);
	  Volume     s_vol(det_name+"_"+l_name+"_"+s_name,s_box,slice_material);
	  //fg: not needed          DetElement slice(layer,s_name,det_id);

	  s_vol.setVisAttributes(theDetector.visAttributes(x_slice.visStr()));
#ifdef VERBOSE
	  std::cout<<"x_slice.materialStr(): "<< x_slice.materialStr() <<std::endl;
#endif
	  if (x_slice.materialStr().compare(x_staves.materialStr()) == 0){
	    radiator_dim_y = s_thick;
	  // W StructureLayer has the same thickness as W radiator layer in the Alveolus layer
	    
#if DD4HEP_VERSION_GE( 0, 15 )
	    caloLayer.outer_nRadiationLengths   = nRadiationLengths;
	    caloLayer.outer_nInteractionLengths = nInteractionLengths;
	    caloLayer.outer_thickness           = thickness_sum;

	    if (!isFirstSens){ caloData->layers.push_back( caloLayer ) ;
#ifdef VERBOSE
	    std::cout<<" caloLayer.distance: "<< caloLayer.distance <<std::endl;

	    std::cout<<" caloLayer.inner_nRadiationLengths: "<< caloLayer.inner_nRadiationLengths <<std::endl;
	    std::cout<<" caloLayer.inner_nInteractionLengths: "<< caloLayer.inner_nInteractionLengths <<std::endl;
	    std::cout<<" caloLayer.inner_thickness: "<< caloLayer.inner_thickness <<std::endl;
	    std::cout<<" caloLayer.sensitive_thickness: "<< caloLayer.sensitive_thickness <<std::endl;

	    std::cout<<" caloLayer.outer_nRadiationLengths: "<< caloLayer.outer_nRadiationLengths <<std::endl;
	    std::cout<<" caloLayer.outer_nInteractionLengths: "<< caloLayer.outer_nInteractionLengths <<std::endl;
	    std::cout<<" caloLayer.outer_thickness: "<< caloLayer.outer_thickness <<std::endl;

	    std::cout<<" EcalBarrel[1]==>caloLayer.inner_thickness + caloLayer.outer_thickness: "
		     << caloLayer.inner_thickness + caloLayer.outer_thickness <<std::endl;
#endif
	    }
#endif
	    // Init for inner
	    nRadiationLengths   = 0. ;
	    nInteractionLengths = 0. ;
	    thickness_sum       = 0. ;	    
	    isFirstSens         = false;

	  }
	  nRadiationLengths   += s_thick/(2.*slice_material.radLength());
	  nInteractionLengths += s_thick/(2.*slice_material.intLength());
	  thickness_sum       += s_thick/2.;

          if ( x_slice.isSensitive() ) {
	    s_vol.setSensitiveDetector(sens);

#if DD4HEP_VERSION_GE( 0, 15 )
	    //Store "inner" quantities
	    caloLayer.inner_nRadiationLengths   = nRadiationLengths ;
	    caloLayer.inner_nInteractionLengths = nInteractionLengths ;
	    caloLayer.inner_thickness           = thickness_sum ;
	    //Store sensitive slice thickness
	    caloLayer.sensitive_thickness       = s_thick ;
#ifdef VERBOSE	    
	    std::cout<<" l_num: "<<l_num <<std::endl;
	    std::cout<<" s_num: "<<s_num <<std::endl;
	    std::cout<<" Ecal_inner_radius: "<< Ecal_inner_radius <<std::endl;
	    std::cout<<" module_thickness: "<< module_thickness <<std::endl;
	    std::cout<<" l_pos_z: "<< l_pos_z <<std::endl;
	    std::cout<<" l_thickness: "<< l_thickness <<std::endl;
	    std::cout<<" s_pos_z: "<< s_pos_z <<std::endl;
	    std::cout<<" s_thick: "<< s_thick <<std::endl;
	    std::cout<<" radiator_dim_y: "<< radiator_dim_y <<std::endl;
#endif	    
	    //-----------------------------------------------------------------------------------------
	    caloLayer.distance  = Ecal_inner_radius + module_thickness/2.0 - l_pos_z + l_thickness/2. + (s_pos_z+s_thick/2.) 
	      - caloLayer.inner_thickness;
	    caloLayer.absorberThickness = radiator_dim_y ;

	    //-----------------------------------------------------------------------------------------
#endif
	    // Init for outer
	    nRadiationLengths   = 0. ;
	    nInteractionLengths = 0. ;
	    thickness_sum       = 0. ;

	  }

	  nRadiationLengths   += s_thick/(2.*slice_material.radLength());
	  nInteractionLengths += s_thick/(2.*slice_material.intLength());
	  thickness_sum       += s_thick/2;


          // Slice placement.
          PlacedVolume slice_phv = l_vol.placeVolume(s_vol,Position(0,0,s_pos_z-s_thick/2));

          if ( x_slice.isSensitive() ) {
	    slice_phv.addPhysVolID("layer", myLayerNum++ );
	  }

          // Increment Z position of slice.
          s_pos_z -= s_thick;
                                        
          // Increment slice number.
          ++s_num;

        }        

#if DD4HEP_VERSION_GE( 0, 15 )
	caloLayer.outer_nRadiationLengths   = nRadiationLengths
	  + (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE))/air.radLength();
	caloLayer.outer_nInteractionLengths = nInteractionLengths
	  + (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE))/air.intLength();
	caloLayer.outer_thickness           = thickness_sum 
	  + (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));

	//Only fill the layers information into DDRec after second layer as Mokka Gear.
	if (!isFirstSens) caloData->layers.push_back( caloLayer ) ;
#ifdef VERBOSE
	std::cout<<" caloLayer.distance: "<< caloLayer.distance <<std::endl;
	
	std::cout<<" caloLayer.inner_nRadiationLengths: "<< caloLayer.inner_nRadiationLengths <<std::endl;
	std::cout<<" caloLayer.inner_nInteractionLengths: "<< caloLayer.inner_nInteractionLengths <<std::endl;
	std::cout<<" caloLayer.inner_thickness: "<< caloLayer.inner_thickness <<std::endl;
	std::cout<<" caloLayer.sensitive_thickness: "<< caloLayer.sensitive_thickness <<std::endl;

	std::cout<<" caloLayer.outer_nRadiationLengths: "<< caloLayer.outer_nRadiationLengths <<std::endl;
	std::cout<<" caloLayer.outer_nInteractionLengths: "<< caloLayer.outer_nInteractionLengths <<std::endl;
	std::cout<<" caloLayer.outer_thickness: "<< caloLayer.outer_thickness <<std::endl;

	std::cout<<" EcalBarrel[2]==>caloLayer.inner_thickness + caloLayer.outer_thickness: "
		 << caloLayer.inner_thickness + caloLayer.outer_thickness <<std::endl;
#endif
#endif
	// Init for next double layer
	nRadiationLengths   = radiator_dim_y/(stave_material.radLength())
	  + (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE))/air.radLength();
	nInteractionLengths = radiator_dim_y/(stave_material.intLength())
	  + (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE))/air.intLength();
	thickness_sum       = radiator_dim_y
	  + (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));
	
	if(radiator_dim_y <= 0) {
	  stringstream err;
	  err << " \n ERROR: The subdetector " << x_det.nameStr() << " geometry parameter -- radiator_dim_y = " << radiator_dim_y ;
	  err << " \n Please check the radiator material name in the subdetector xml file";
	  throw runtime_error(err.str());
	}

	// #########################
	// BuildBarrelStructureLayer
	// #########################


	l_dim_x -=  (radiator_dim_y +  Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));
	double radiator_dim_x = l_dim_x*2.;

#ifdef VERBOSE
	std::cout << "radiator_dim_x = " << radiator_dim_x << std::endl;
#endif  

	double radiator_dim_z =
	  Ecal_Barrel_module_dim_z -
	  2 * Ecal_lateral_face_thickness -
	  2 * N_FIBERS_W_STRUCTURE * Ecal_fiber_thickness;
	
	string bs_name="bs";
	Box        barrelStructureLayer_box(radiator_dim_x/2.,radiator_dim_z/2.,radiator_dim_y/2.);
	Volume     barrelStructureLayer_vol(det_name+"_"+l_name+"_"+bs_name,barrelStructureLayer_box,stave_material);

	barrelStructureLayer_vol.setVisAttributes(theDetector.visAttributes(x_layer.visStr()));	


        // Increment to next layer Z position.
        l_pos_z -= l_thickness;          

	// Without last W StructureLayer, the last part is Si SD even layer.
	// the last number of  Ecal_nlayers1, Ecal_nlayers2 and  Ecal_nlayers3 is odd.
	int even_layer = l_num*2;
	if(even_layer > Ecal_nlayers1 + Ecal_nlayers2 + Ecal_nlayers3) continue;
	//if ( Number_of_Si_Layers_in_Barrel > n_total_layers ) continue;

	double bsl_pos_z = l_pos_z - (radiator_dim_y/2. + Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));
	l_pos_z -= (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));

	Position   bsl_pos(0,0,bsl_pos_z);      // Position of the layer.
	// PlacedVolume  barrelStructureLayer_phv =
	  mod_vol.placeVolume(barrelStructureLayer_vol,bsl_pos);

	l_dim_x -=  (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));	
	l_pos_z -= (radiator_dim_y + Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));
	
        ++l_num;

      }
    }
  

    // Set stave visualization.
    if (x_staves)   {
      mod_vol.setVisAttributes(theDetector.visAttributes(x_staves.visStr()));
    }
    
    
    
//====================================================================
// Place ECAL Barrel stave module into the envelope volume
//====================================================================
    double X,Y;
    X = module_thickness * sin(M_PI/4.);
    Y = Ecal_inner_radius + module_thickness / 2.;
    
    for (int stave_id = 1; stave_id <= nsides ; stave_id++)
      for (int module_id = 1; module_id < 6; module_id++)
	{
	  double phirot =  (stave_id-1) * dphi - hphi;
	  double module_z_offset =  (2 * module_id-6) * Ecal_Barrel_module_dim_z/2.;
	  
	  // And the rotation in Mokka is right hand rule, and the rotation in DD4hep is clockwise rule
	  // So there is a negitive sign when port Mokka into DD4hep 
	  Transform3D tr(RotationZYX(0,phirot,M_PI*0.5),Translation3D(X*cos(phirot)-Y*sin(phirot),
								      X*sin(phirot)+Y*cos(phirot),
								      module_z_offset));
	  PlacedVolume pv = envelope.placeVolume(mod_vol,tr);
	  pv.addPhysVolID("module",module_id);
	  pv.addPhysVolID("stave",stave_id);
	  DetElement sd = (module_id==0&&stave_id==0) ? stave_det : stave_det.clone(_toString(module_id,"module%d")+_toString(stave_id,"stave%d"));
	  sd.setPlacement(pv);
	  sdet.add(sd);
	  
	}
    
    // Set envelope volume attributes.
    envelope.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
    
    sdet.addExtension< LayeredCalorimeterData >( caloData ) ; 


    return sdet;
}

DECLARE_DETELEMENT(SEcal04_Barrel_v01,create_detector)

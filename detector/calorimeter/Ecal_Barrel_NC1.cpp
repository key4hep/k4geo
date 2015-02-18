//====================================================================
//  DDSim - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for SiWEcalBarrel
//  Ported from Mokka, adapted for DDRec
//--------------------------------------------------------------------
//  Mokka port - S.Lu, DESY 
//  DD4Rec adaptation - D.Protopopescu
//  $Id$
//====================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "TGeoTrd2.h"

#include "DDRec/Extensions/LayeringExtensionImpl.h"
#include "DDRec/Extensions/SubdetectorExtensionImpl.h"
#include "DDRec/Surface.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

#define VERBOSE 1

static Ref_t create_ecal_barrel(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
  static double tolerance = 0e0;
  Material      air       = lcdd.air();
  //unused:  Material      vacuum    = lcdd.vacuum();

  xml_det_t     x_det     = e;
  Layering      layering (e);
  string        det_name  = x_det.nameStr();
  int           det_id    = x_det.id();
  xml_dim_t     dim       = x_det.dimensions();
  int           numsides  = dim.numsides();
  double        rmin      = dim.rmin();

  xml_comp_t    x_staves  = x_det.staves();
  xml_dim_t     s_dim     = x_staves.dimensions();
  int           ntowers   = s_dim.nmodules();
  double        tower_z   = 2.*dim.z()/(ntowers);  

  xml_comp_t    towers    = x_staves.trd();
  xml_dim_t     t_dim     = towers.dimensions();
  int           nlayers   = t_dim.nmodules();

  DetElement    sdet (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);

  Assembly envelope_assembly( det_name + "assembly"  ) ;  
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope_assembly);

  env_phv.addPhysVolID("system",det_id);
  sdet.setPlacement(env_phv);

  sens.setType("calorimeter");

  DetElement stave_det("module0stave0",det_id);

  // cannot get rid of those with the current layers/slices structure !!
  int N_FIBERS_W_STRUCTURE            = lcdd.constant<int>("Ecal_barrel_number_of_fibers_W_structure");//2; 
  int N_FIBERS_ALVOULUS               = lcdd.constant<int>("Ecal_barrel_number_of_fibers_alveolus");//3;
  double Ecal_fiber_thickness         = lcdd.constant<double>("Ecal_fiber_thickness");
  double Ecal_front_face_thickness    = lcdd.constant<double>("Ecal_front_face_thickness");
  double Ecal_lateral_face_thickness  = lcdd.constant<double>("Ecal_lateral_face_thickness");
  //unused:  double Ecal_Slab_H_fiber_thickness  = lcdd.constant<double>("Ecal_Slab_H_fiber_thickness");

#ifdef VERBOSE
  std::cout << " barrel_halfz  = " << dim.z()  << std::endl;
  std::cout << " tower_dim_z  = " << tower_z  << std::endl;
#endif

  // module barrel key parameters
  double module_thickness = s_dim.thickness();
  double si_slab_thickness = t_dim.thickness(); 
  double bottom_dim_x = 2.*tan(M_PI/numsides)*rmin + module_thickness/sin(2.*M_PI/numsides);
  double top_dim_x = bottom_dim_x - 2.*module_thickness;

#ifdef VERBOSE
  std::cout << " si_slab_thickness = " << si_slab_thickness  << std::endl;
  std::cout << " module_thickness = " << module_thickness  << std::endl;
  std::cout << " bottom_dim_x = " << bottom_dim_x  << std::endl;
  std::cout << " top_dim_x = " << top_dim_x << std::endl;
#endif
  

// ========= Create Ecal Barrel stave   ====================================
//  It will be the volume for palcing the Ecal Barrel alveolus(i.e. Layers).
//  And the structure W plate.
//  Itself will be placed into the world volume.
// ==========================================================================

  // top_x and bottom_x differ between Mokka and DD4hep
  Trapezoid trd(top_dim_x / 2,
		bottom_dim_x / 2, 
		tower_z / 2,
		tower_z / 2,
		module_thickness/2);

  Volume mod_vol(det_name+"_module",trd,air);

  // We count the layers starting from IP and from 1,
  // so odd layers should be inside slabs and even ones on the structure.
  // The structure W layers are here big plans, as the 
  // gap between each W plate is too small to create problems 
  // The even W layers are part of H structure placed inside the alveolus.

  //  Dimension of radiator wLog slice provides the thickness
  double y_floor =  Ecal_front_face_thickness + N_FIBERS_ALVOULUS * Ecal_fiber_thickness;

  // ############################
  //  Dimension of alveolus slice provides the thickness
  // ############################
  
    // =====  build Si Slab and put into the Layer volume =====
    // =====  place the layer into the module 5 time for one full layer into the trd module ====
    // =====  build and place barrel structure into trd module ====

    // Parameters for computing the layer X dimension:
    double stave_z  = (tower_z - 2.*Ecal_lateral_face_thickness) / (2.*ntowers);
    double l_dim_x  = bottom_dim_x/2.;     // Starting X dimension for the layer
    double l_pos_z  = module_thickness/2.;

    l_dim_x -= y_floor;
    l_pos_z -= y_floor;

    // ------------------ create extension objects for reconstruction ------------------
    //unused:   DDRec::LayeringExtensionImpl* layeringExtension = new DDRec::LayeringExtensionImpl ;
    DDRec::SubdetectorExtensionImpl* subDetExtension = new DDRec::SubdetectorExtensionImpl( sdet )  ;
    Position layerNormal(0,0,1); //fg: defines direction of thickness in Box for layer slices
    
    subDetExtension->setIsBarrel(true);
    subDetExtension->setNSides( numsides );
    subDetExtension->setRMin( rmin );
    subDetExtension->setRMax( (rmin + module_thickness )/cos(M_PI/numsides) );
    subDetExtension->setZMin( 0. );
    subDetExtension->setZMax( dim.z() ); //Ecal_Barrel_halfZ
    
    // base vectors for surfaces:
    DDSurfaces::Vector3D u(1,0,0) ;
    DDSurfaces::Vector3D v(0,1,0) ;
    DDSurfaces::Vector3D n(0,0,1) ;

    //------------------------------------------------------------------------------------

    //-------------------- start loop over ECAL layers ----------------------
    // Loop over the sets of layer elements in the detector.
    int l_num = 1;
    for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
      xml_comp_t x_layer = li;
      int repeat = x_layer.repeat();
      // Loop over number of repeats for this layer.
      for (int j=0; j< repeat; j++)    {
	string l_name = _toString(l_num,"layer%d");
	double l_thickness = layering.layer(l_num-1)->thickness();  // layer thickness
	double xcut = (l_thickness);                     // X dimension for this layer
	l_dim_x -= xcut;

#ifdef VERBOSE
	std::cout << " layer "<< l_num << std::endl;
#endif

	Box        l_box(l_dim_x-tolerance,stave_z-tolerance,l_thickness/2.0-tolerance);
	Volume     l_vol(det_name+"_"+l_name,l_box,air);

	l_vol.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));

	// fg: need vector of DetElements for towers ! 
	// DetElement layer(stave_det, l_name, det_id);
	std::vector< DetElement > layers( ntowers );
	
	// place layer 5 times in module. at same layer position (towers !)
	double l_pos_y = tower_z/2.;
       	for (int i=0; i < ntowers; i++){ // need four clones

	  layers[i] = DetElement( stave_det, l_name+_toString(i,"tower%02d") , det_id ) ;
	  l_pos_y -= stave_z;
	  Position   l_pos(0, l_pos_y, l_pos_z - l_thickness/2.);  // Position of the layer
	  PlacedVolume layer_phv = mod_vol.placeVolume(l_vol,l_pos);
	  layer_phv.addPhysVolID("layer", l_num);
	  layer_phv.addPhysVolID("tower", i);

	  layers[i].setPlacement(layer_phv);
	  l_pos_y -= stave_z;
	}

	//--------------------------------------------------------------------------------
        // BuildBarrelAlveolus: BuildSiliconSlab:
	//--------------------------------------------------------------------------------
	int    s_num = 1;
	double s_pos_z = -( l_thickness / 2.);
	double radiator_dim_y = -1.0; // to be updated with slice radiator thickness 

	// Loop over slices for this layer
	for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	  xml_comp_t x_slice = si;
	  string     s_name  = _toString(s_num,"slice%d");
	  double     s_thick = x_slice.thickness();

	  double slab_dim_x = l_dim_x-tolerance;
	  double slab_dim_y = s_thick/2.;
	  double slab_dim_z = stave_z-tolerance;

	  Box        s_box(slab_dim_x, slab_dim_z, slab_dim_y);
	  Volume     s_vol(det_name+"_"+l_name+"_"+s_name,s_box,lcdd.material(x_slice.materialStr()));

	  s_vol.setVisAttributes(lcdd.visAttributes(x_slice.visStr()));

#ifdef VERBOSE
	  std::cout << "  x_slice.materialStr(): "<< x_slice.materialStr() << std::endl;
#endif
	  if (x_slice.materialStr().compare(x_staves.materialStr()) == 0){
            // W StructureLayer has the same thickness as W radiator layer in the Alveolus layer
	    radiator_dim_y = s_thick;
	  }
          if ( x_slice.isSensitive() ) {
	    s_vol.setSensitiveDetector(sens);
	  }

          PlacedVolume slice_phv = l_vol.placeVolume(s_vol,Position(0,0,s_pos_z+s_thick/2));

          if ( x_slice.isSensitive() ) {
	    slice_phv.addPhysVolID("slice",s_num);

	    // add a measurement surface to the layer for every sensitive slice:
	    DDRec::VolPlane surf( s_vol , DDSurfaces::SurfaceType(DDSurfaces::SurfaceType::Sensitive) , 
                                  slab_dim_y , slab_dim_y , u,v,n ) ; //,o ) ;

	    // add them to the layers of all towers
	    for (int i=0; i< ntowers; i++){
	      DDRec::volSurfaceList( layers[i] )->push_back( surf ) ;
	    }
	  }

          // Increment Z position of slice.
          s_pos_z += s_thick;
                                        
          // Increment slice number.
          ++s_num;
        }        

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
	std::cout << " radiator_dim_x = " << radiator_dim_x << std::endl;
#endif  

	double radiator_dim_z = tower_z -
	  2 * Ecal_lateral_face_thickness -
	  2 * N_FIBERS_W_STRUCTURE * Ecal_fiber_thickness;
	
	string bs_name = "bs";
	
	Box     barrelStructureLayer_box(radiator_dim_x/2., radiator_dim_z/2., radiator_dim_y/2.);
	Volume  barrelStructureLayer_vol(det_name+"_"+l_name+"_"+bs_name,barrelStructureLayer_box,lcdd.material(x_staves.materialStr()));

	barrelStructureLayer_vol.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));	

        // Increment to next layer Z position.
        l_pos_z -= l_thickness;          

	// Without last W StructureLayer, the last part is Si SD even layer.
	// the last number of Ecal_nlayers1, Ecal_nlayers2 and  Ecal_nlayers3 is odd.
	int even_layer = l_num*2;
	if(even_layer > nlayers ) continue;

	double bsl_pos_z = l_pos_z - (radiator_dim_y/2. + Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));
	l_pos_z -= (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));

	Position   bsl_pos(0,0,bsl_pos_z);      // Position of the layer.
	//unused: PlacedVolume  barrelStructureLayer_phv = 
	mod_vol.placeVolume(barrelStructureLayer_vol,bsl_pos);

	l_dim_x -=  (Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));	
	l_pos_z -= (radiator_dim_y + Ecal_fiber_thickness * (N_FIBERS_ALVOULUS + N_FIBERS_W_STRUCTURE));
	
        ++l_num;
      }
    }
  
    // Set stave visualization.
    if (x_staves)   {
      mod_vol.setVisAttributes(lcdd.visAttributes(x_staves.visStr()));
    }

        
//====================================================================
// Place ECAL Barrel stave module into the assembly envelope
//====================================================================
    double X,Y;
    X = module_thickness * sin(2.*M_PI/numsides);
    Y = rmin + module_thickness / 2.;
    
    for (int stave_id = 1; stave_id <= numsides ; stave_id++)
      for (int module_id = 1; module_id <= ntowers; module_id++)
	{
	  double phirot =  (stave_id - 1)*(2.*M_PI/numsides);
	  double module_z_offset =  (2*module_id - ntowers - 1)*tower_z/2.;
	  
	  // The rotation in Mokka is right hand rule, and the rotation in DD4hep is clockwise
	  // So there is a negative sign when port Mokka into DD4hep 
	  Transform3D tr(RotationZYX(0, phirot, M_PI/2.),Translation3D(X*cos(phirot)-Y*sin(phirot),
								      X*sin(phirot)+Y*cos(phirot),
								      module_z_offset));
	  PlacedVolume pv = envelope_assembly.placeVolume(mod_vol,tr);
	  pv.addPhysVolID("system", det_id);
	  pv.addPhysVolID("module", module_id);
	  pv.addPhysVolID("stave", stave_id);
	  DetElement sd = (module_id==0&&stave_id==0) ? stave_det : stave_det.clone(
                           _toString(module_id,"module%d") + _toString(stave_id,"stave%d"));
	  sd.setPlacement(pv);
	  sdet.add(sd);
	  
	}
    
    // Set envelope_assembly volume attributes.
    envelope_assembly.setAttributes(lcdd, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
    return sdet;
}

DECLARE_DETELEMENT(Ecal_Barrel_NC1, create_ecal_barrel)

//====================================================================
//  DDSim - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for HcalEndcapRing
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id$
//====================================================================

 /* History:  
  
   initial version: 
   F.Gaede: identical to  Hcal03 driver except that an additional gap
            for the fibres is introduced between scintillator and steel
   PMoradeFreitas: Super driver without tmp database and able
                   to build Hcal barrel with just two modules in stave.
   Angela Lucaci: similar to SHcal03, just that the drift chambers option is 
                  not considered anymore, only the scintillator one. In addition,
                  a gap in the middle of the stave is build, and fractional cells
                  at the edges are introduced (see SDHcalBarrel.cc)
   Ralf Diener: Corrected use of GEAR interface
   Andre Sailer: Added Tungsten
   Andre Sailer: Added possibility for different endcap/barrel material
                 Which needs three additional parameters: 
		 endcap_radiator_thickness, endcap_radiator_material, endcap_layers
                 Also added parameters to the database for this driver
   Shaojun Lu:  Barrel driver has been changed for the new engineering design shape.
                Endcaps driver has been changed for the new engineering design shape.
		Updated for flexibility.
   Shaojun Lu:  Ported from Mokka SHcalSc04 EndcapRing part. Read the constants from XML
                instead of the DB. Then build the EndcapRing in the same way with DD4hep
		construct.
 */

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "DD4hep/Shapes.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

//#define VERBOSE 1

static Ref_t create_detector(LCDD& lcdd, xml_h element, SensitiveDetector sens)  {
  //unused:  static double tolerance = 0e0;

  xml_det_t   x_det     = element;
  string      det_name    = x_det.nameStr();
  Layering    layering(x_det);

  Material    air         = lcdd.air();
  Material    Steel235    = lcdd.material(x_det.materialStr());

  int           det_id    = x_det.id();
  xml_comp_t    x_staves  = x_det.staves();
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = lcdd.pickMotherVolume(sdet);


  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = XML::createPlacedEnvelope( lcdd,  element , sdet ) ;
  
  if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  sens.setType("calorimeter");

  DetElement    module_det("module0",det_id);

  // Some verbose output
  cout << " \n\n\n CREATE DETECTOR: SHcalSC04_EndcapRing" << endl;




//====================================================================
//
// Read all the constant from ILD_o1_v05.xml
// Use them to build HcalEndcapRing
//
//====================================================================
  // The way to read constant from XML/LCDD file.
  double      Hcal_radiator_thickness          = lcdd.constant<double>("Hcal_radiator_thickness");
  double      Hcal_chamber_thickness           = lcdd.constant<double>("Hcal_chamber_thickness");
  double      Hcal_inner_radius                = lcdd.constant<double>("Hcal_inner_radius");
  double      Hcal_back_plate_thickness        = lcdd.constant<double>("Hcal_back_plate_thickness");
  double      Hcal_lateral_plate_thickness     = lcdd.constant<double>("Hcal_lateral_structure_thickness");
  double      Hcal_stave_gaps                  = lcdd.constant<double>("Hcal_stave_gaps");
  double      Hcal_modules_gap                 = lcdd.constant<double>("Hcal_modules_gap"); 
  double      TPC_Ecal_Hcal_barrel_halfZ       = lcdd.constant<double>("TPC_Ecal_Hcal_barrel_halfZ"); 
  //unused:  double      Hcal_middle_stave_gaps           = lcdd.constant<double>("Hcal_middle_stave_gaps");
  //unused:  double      Hcal_layer_air_gap               = lcdd.constant<double>("Hcal_layer_air_gap");

  double      Hcal_module_radius               = lcdd.constant<double>("Hcal_outer_radius");

  int         Hcal_nlayers                     = lcdd.constant<int>("Hcal_nlayers");

  //double      Hcal_inner_radius                = lcdd.constant<double>("HcalBarrel_rmin");
  //double      Hcal_outer_radius                = lcdd.constant<double>("Hcal_outer_radius");
  //double      Hcal_Barrel_rotation             = lcdd.constant<double>("Hcal_Barrel_rotation");

  double      Ecal_endcap_zmin                 = lcdd.constant<double>("Ecal_endcap_zmin");
  double      Hcal_endcap_radiator_thickness   = lcdd.constant<double>("Hcal_endcap_radiator_thickness");
  double      Hcal_radial_ring_inner_gap       = lcdd.constant<double>("Hcal_radial_ring_inner_gap");
  double      Hcal_endcap_cables_gap           = lcdd.constant<double>("Hcal_endcap_cables_gap");
  double      Hcal_endcap_ecal_gap             = lcdd.constant<double>("Hcal_endcap_ecal_gap");

  double      Hcal_Cu_thickness                = lcdd.constant<double>("Hcal_Cu_thickness");
  double      Hcal_PCB_thickness               = lcdd.constant<double>("Hcal_PCB_thickness");

  int         Hcal_endcap_nlayers              = lcdd.constant<int>("Hcal_endcap_nlayers");

  //unused:  int         Hcal_ring                        = lcdd.constant<int>("Hcal_ring");

  //TODO: thinking about how to pass the updated value at runtime from other inner drivers? 
  double      Ecal_endcap_zmax                 = 263.5; //cm //= lcdd.constant<double>("Ecal_endcap_zmax");
  double      Ecal_endcap_outer_radius         = 208.88;//cm //= lcdd.constant<double>("Ecal_endcap_outer_radius");
  double      Hcal_cells_size                  = lcdd.constant<double>("Hcal_cells_size");

//====================================================================
//
// general calculated parameters
//
//====================================================================

  double Hcal_total_dim_y   = Hcal_nlayers * (Hcal_radiator_thickness + Hcal_chamber_thickness) 
                            + Hcal_back_plate_thickness;

  
  // Moved the calculation into ILD_o1_v05.xml
  // double Hcal_module_radius = Hcal_inner_radius /  cos(M_PI/8.)
  //   + ( Hcal_radiator_thickness + Hcal_chamber_thickness ) * Hcal_nlayers
  //   + Hcal_back_plate_thickness;
  // std::cout << " *********** hcal outer radius : " << Hcal_module_radius  << std::endl ;


  double Hcal_y_dim1_for_x  = Hcal_module_radius*cos(M_PI/8.) - Hcal_inner_radius;
  //unused:  double Hcal_bottom_dim_x  = 2.*Hcal_inner_radius*tan(M_PI/8.)- Hcal_stave_gaps;
  double Hcal_normal_dim_z  = (2 * TPC_Ecal_Hcal_barrel_halfZ - Hcal_modules_gap)/2.;

 //only the middle has the steel plate.
  //unused:  double Hcal_regular_chamber_dim_z = Hcal_normal_dim_z - Hcal_lateral_plate_thickness;


  //double totalThickness_Hcal_Barrel       = (Hcal_radiator_thickness + Hcal_chamber_thickness) 
  //                                                 * HcalBarrel_layers + Hcal_back_plate_thickness;
  //double Hcal_module_radius               = Hcal_inner_radius + totalThickness_Hcal_Barrel;

  //double Hcal_cell_dim_x            = Hcal_cells_size;
  //double Hcal_cell_dim_z            = Hcal_regular_chamber_dim_z / floor (Hcal_regular_chamber_dim_z/Hcal_cell_dim_x);
  //double Hcal_digi_cell_dim_x       = Hcal_digi_cells_size;
  //double Hcal_digi_cell_dim_z       = Hcal_regular_chamber_dim_z / floor (Hcal_regular_chamber_dim_z/Hcal_digi_cell_dim_x);

 

  // just two modules per stave
  //double Hcal_normal_dim_z = (2 * TPC_Ecal_Hcal_barrel_halfZ - Hcal_modules_gap)/2.;
  double Hcal_start_z =  Hcal_normal_dim_z + Hcal_modules_gap / 2. + Hcal_endcap_cables_gap;

  cout<<"    HCAL endcap rings: Hcal_start_z:  "<< Hcal_start_z <<endl<<endl;

  // Hcal_start_z is the Hcal Endcap boundary coming from the IP
  // Test Hcal_start_z against Ecal_endcap_zmax + Hcal_endcap_ecal_gap
  // to avoid overlap problems with Ecal if scaled.
  //
  double Hcal_endcap_rmax   = Hcal_inner_radius + Hcal_y_dim1_for_x;

  if( Hcal_start_z < Ecal_endcap_zmax + Hcal_endcap_ecal_gap )
    Hcal_start_z = Ecal_endcap_zmax + Hcal_endcap_ecal_gap;

  cout<<"    HCAL endcap rings: Hcal_start_z:  "<< Hcal_start_z <<endl<<endl;






// ========= Create Hcal end cap ring   ====================================
//  It will be the volume for palcing the Hcal endcaps Ring Layers.
//  And the absorber plate.
//  Itself will be placed into the world volume.
// ==========================================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~              EndcapRings                          ~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // old parameters from database
  double pRMax, pDz, pRMin;
  pRMax = Hcal_endcap_rmax;

  // The rings start from inner Ecal endcap boundary
  // and finish at inner Hcal endcap one.
  double start_z, stop_z;
  start_z = Ecal_endcap_zmin;
  double SpaceForLayers = Hcal_start_z - Hcal_endcap_ecal_gap
    - Ecal_endcap_zmin - Hcal_back_plate_thickness;

  int MaxNumberOfLayers = (int) (SpaceForLayers /
				     (Hcal_chamber_thickness + Hcal_endcap_radiator_thickness));


  cout<<"    HCAL endcap rings: Hcal_start_z:  "<< Hcal_start_z <<endl;
  cout<<"    HCAL endcap rings: Hcal_endcap_rmax: "<< Hcal_endcap_rmax <<endl;
  cout<<"    HCAL endcap rings: SpaceForLayers:  "<< SpaceForLayers <<endl;
  cout<<"    HCAL endcap rings will have "<< MaxNumberOfLayers << " layers."<<endl;

  stop_z = start_z + MaxNumberOfLayers * (Hcal_chamber_thickness + Hcal_endcap_radiator_thickness)
    + Hcal_back_plate_thickness;

  pDz = (stop_z - start_z) / 2.;

  pRMin = Ecal_endcap_outer_radius
    + Hcal_radial_ring_inner_gap;

  //double zPlane[2];
  //zPlane[0]=-pDz;
  //zPlane[1]=-zPlane[0];

  double zlen = pDz*2.;

  //double rInner[2],rOuter[2];
  //rInner[0]=rInner[1]=pRMin;
  //rOuter[0]=rOuter[1]=pRMax;

  double rmin = pRMin;
  double rmax = pRMax;
  int numSide = 8;

  PolyhedraRegular HcalEndCapRingSolid( numSide, M_PI/8., rmin, rmax,  zlen);

  Volume  HcalEndCapRingLogical("HcalEndCapRingLogical",HcalEndCapRingSolid, Steel235);


  //========== fill data for reconstruction ============================
  DDRec::LayeredCalorimeterData* caloData = new DDRec::LayeredCalorimeterData ;
  caloData->layoutType = DDRec::LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = numSide  ;
  caloData->outer_symmetry = numSide  ;
  caloData->phi0 = 0 ; // hardcoded 

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rmin ;
  caloData->extent[1] = rmax ;
  caloData->extent[2] = start_z ;
  caloData->extent[3] = stop_z  ;


  //==============================================================
  // build the layer and place into the Hcal EndcapRing
  //==============================================================



  double lpRMax, lpDz, lpRMin;

  //lpRMax = rmax * cos(M_PI/numSide) - Hcal_lateral_plate_thickness;
  lpRMax = rmax - Hcal_lateral_plate_thickness;

  lpDz = Hcal_chamber_thickness / 2.;

  //lpRMin = rmin * cos(M_PI/numSide) + Hcal_lateral_plate_thickness;
  lpRMin = rmin + Hcal_lateral_plate_thickness;

  // G4Polyhedra Envelope parameters
  double phiStart = M_PI/8.;
  //double phiTotal = 2.*pi;
  //int numSide     = motherPolyhedraParameters->numSide;
  //int numZPlanes  = 2;

  //double zPlane[2];
  //zPlane[0] = - pDz;
  //zPlane[1] = - zPlane[0];

  //double rInner[2],rOuter[2];
  //rInner[0] = rInner[1] = pRMin;
  //rOuter[0] = rOuter[1] = pRMax;

  /*
#ifdef SHCALSC04_DEBUG
  if(rings==true){
    G4cout<<"    EndcapRingsSensitive: Rinner="<<pRMin<<"mm Router="<<pRMax<<G4endl<<G4endl;
  }
#endif
  */

  double lzlen = lpDz*2.;

  PolyhedraRegular HcalEndCapRingChamberSolid( numSide, phiStart, lpRMin, lpRMax,  lzlen);

  Box IntersectionStaveBox( lpRMax, lpRMax, Hcal_total_dim_y);

  // set up the translation and rotation for the intersection process 
  // this happens in the mother volume coordinate system, so a coordinate transformation is needed
  Position IntersectXYZtrans((pRMax + (Hcal_stave_gaps/2.)), 
			     (pRMax + (Hcal_stave_gaps/2.)),
			     (Hcal_total_dim_y/2.));

  RotationZYX rot(0.,0.,0.);
  Transform3D tran3D(rot,IntersectXYZtrans);
  // intersect the octagonal layer with a square to get only one quadrant
  IntersectionSolid  HcalEndCapRingStaveSolid( HcalEndCapRingChamberSolid, IntersectionStaveBox, tran3D); 

  int EC_Number_of_towers = 0; 



  // chamber placements
  int number_of_chambers = Hcal_endcap_nlayers;
  int possible_number_of_chambers = (int) floor( zlen / (Hcal_chamber_thickness + Hcal_endcap_radiator_thickness));
  if(possible_number_of_chambers < number_of_chambers)
    number_of_chambers = possible_number_of_chambers;


  //place the four staves in their right positions
  for (int stave_id = 1;
       stave_id <= 4;
       stave_id++)
    {
      
      for (int layer_id = 1;
	   layer_id <= number_of_chambers;
	   layer_id++)
	{
	  
	  double Zoff = - ( zlen/2.)
	    + (layer_id-1) *(Hcal_chamber_thickness + Hcal_endcap_radiator_thickness)
	    + Hcal_endcap_radiator_thickness 
	    + (Hcal_chamber_thickness - Hcal_PCB_thickness - Hcal_Cu_thickness)/2.;
	  
	  double layer_thickness = lzlen;
	  
	  //====================================================================
	  // Create Hcal EndcapRing Chamber without radiator
	  // Place into the Hcal EndcapRing logical, after each radiator 
	  //====================================================================
	  xml_coll_t c(x_det,_U(layer));
	  xml_comp_t   x_layer = c;
	  string layer_name      = det_name+_toString(layer_id,"_layer%d");
	  
	  Volume HcalEndCapRingStaveLogical("HcalEndCapRingStaveLogical",HcalEndCapRingStaveSolid, air);

	  // Create the slices (sublayers) within the Hcal Barrel Chamber.
	  double slice_pos_z = -(layer_thickness/2.);
	  int slice_number = 0;
	  for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
	    xml_comp_t x_slice = k;
	    string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
	    double   slice_thickness = x_slice.thickness();
	    Material slice_material  = lcdd.material(x_slice.materialStr());
	    DetElement slice(layer_name,_toString(slice_number,"slice%d"),x_det.id());
	    
	    slice_pos_z += slice_thickness/2.;
	    
	    // Slice volume
	    PolyhedraRegular slicePolyhedraRegularSolid( numSide, phiStart, lpRMin, lpRMax,  slice_thickness);
	    
	    Box sliceIntersectionStaveBox( lpRMax, lpRMax, Hcal_total_dim_y);
	    
	    // set up the translation and rotation for the intersection process 
	    // this happens in the mother volume coordinate system, so a coordinate transformation is needed
	    Position IntersectXYZtrans((pRMax + (Hcal_stave_gaps/2.)), 
				       (pRMax + (Hcal_stave_gaps/2.)),
				       (Hcal_total_dim_y/2.));
	    
	    RotationZYX rot(0.,0.,0.);
	    Transform3D tran3D(rot,IntersectXYZtrans);
	    // intersect the octagonal layer with a square to get only one quadrant
	    IntersectionSolid  slice_Solid( slicePolyhedraRegularSolid, sliceIntersectionStaveBox, tran3D); 
	    
	    //Volume HcalEndCapRingStaveLogical("HcalEndCapRingStaveLogical",HcalEndCapRingStaveSolid, air);
	    
	    Volume slice_vol(slice_name,slice_Solid,slice_material);
	    
	    if ( x_slice.isSensitive() ) {
	      slice_vol.setSensitiveDetector(sens);
	    }
	    // Set region, limitset, and vis.
	    slice_vol.setAttributes(lcdd,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	    // slice PlacedVolume
	    PlacedVolume slice_phv = HcalEndCapRingStaveLogical.placeVolume(slice_vol,Position(0.,0.,slice_pos_z));
	    if ( x_slice.isSensitive() ) {
	      slice_phv.addPhysVolID("layer",layer_id).addPhysVolID("slice",slice_number);
	      cout<<"  layer_id:  "<< layer_id<<"   slice_id:  "<< slice_number<<endl;
	    }
	    
	    slice.setPlacement(slice_phv);
	    // Increment x position for next slice.
	    slice_pos_z += slice_thickness/2.;
	    // Increment slice number.
	    ++slice_number;             
	  }
	  
	  
	  
	  string l_name = _toString(layer_id,"layer%d");
	  string stave_name = _toString(stave_id,"stave%d");
	  DetElement layer(module_det, l_name+stave_name, det_id);
	  
	  double angle_module = M_PI/2. * ( stave_id );
	  
	  Position l_pos(0., 0., Zoff);
	  RotationZ rotz(angle_module);
	  Position l_new = rotz*l_pos;  
	  
	  RotationZYX rot(angle_module,0,0);
	  Transform3D tran3D(rot,l_new);
	      
	  PlacedVolume layer_phv = HcalEndCapRingLogical.placeVolume(HcalEndCapRingStaveLogical,tran3D);
	  layer_phv.addPhysVolID("layer", layer_id);
	  layer_phv.addPhysVolID("tower", EC_Number_of_towers);
	  layer_phv.addPhysVolID("stave", stave_id);
	  layer.setPlacement(layer_phv);


	  //-----------------------------------------------------------------------------------------
	  if ( caloData->layers.size() <= (unsigned int)number_of_chambers ) {
	    DDRec::LayeredCalorimeterData::Layer caloLayer ;

	    caloLayer.distance = Ecal_endcap_zmin + pDz + Zoff;
	    caloLayer.thickness = Hcal_chamber_thickness + Hcal_endcap_radiator_thickness ;
	    caloLayer.absorberThickness = Hcal_endcap_radiator_thickness ;
	    caloLayer.cellSize0 = Hcal_cells_size ;
	    caloLayer.cellSize1 = Hcal_cells_size ;
	    
	    caloData->layers.push_back( caloLayer ) ;
	  }
	  //-----------------------------------------------------------------------------------------

	}


      /*
     theSD->AddLayer(layer_id,
		      0,
		      0,
		      Zoff);

#ifdef MOKKA_GEAR
      // count the layers
      if(rings==true){
	helpEndcapRing.count += 1;

	// position of layer as offset + half thickness
	helpEndcapRing.layerPos.push_back(Zoff + std::abs(zPlane[0]) + (Hcal_PCB_thickness + Hcal_Cu_thickness)/2 ) ;

	helpEndcapRing.sensThickness.push_back( Hcal_scintillator_thickness );
	helpEndcapRing.steelCassetteThickness.push_back( Hcal_steel_cassette_thickness );
	helpEndcapRing.PCBThickness.push_back(Hcal_PCB_thickness);
	helpEndcapRing.CuThickness.push_back(Hcal_Cu_thickness);
      }


#endif
      */

    }  







  // Set stave visualization.
  if (x_staves)   {
    HcalEndCapRingLogical.setVisAttributes(lcdd.visAttributes(x_staves.visStr()));
   }
  

  //====================================================================
  // Place Hcal Endcap Ring module into the assembly envelope volume
  //====================================================================
  
  double endcap_z_offset = Ecal_endcap_zmin + pDz;
  
  for(int module_num=0;module_num<2;module_num++) {

    int module_id = ( module_num == 0 ) ? 0:6;
    double this_module_z_offset = ( module_id == 0 ) ? -endcap_z_offset : endcap_z_offset; 
    double this_module_rotY = ( module_id == 0 ) ? 0:M_PI; 
  
    Position xyzVec(0,0,this_module_z_offset);
    RotationZYX rot(0,this_module_rotY,0);
    Rotation3D rot3D(rot);
    Transform3D tran3D(rot3D,xyzVec);

    PlacedVolume pv = envelope.placeVolume(HcalEndCapRingLogical,tran3D);
    pv.addPhysVolID("module",module_id); // z: +/-

    DetElement sd = (module_num==0) ? module_det : module_det.clone(_toString(module_num,"module%d"));
    sd.setPlacement(pv);

  }
  
  sdet.addExtension< DDRec::LayeredCalorimeterData >( caloData ) ;

  return sdet;
  
}



DECLARE_DETELEMENT(SHcalSc04_EndcapRing, create_detector)

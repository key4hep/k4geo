//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  DD4hep Geometry driver for HcalEndcapRing
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id: SHcalSc04_EndcapRing.cpp 1034 2016-07-26 15:45:32Z shaojun $
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
   Shaojun Lu:  Adapt to the DD4hep envelope.

   Tibor Kurca: modify for SDHcal 
                ... doesn't correspond yet exactly to the technical design: 8 staves, 
                inner/outer symmetries 8/16 
                - check if requested number of chambers fits into the dedicated space
                - forsee different inner/outer symmetries corresponding to EcalEndcap/HcalEndcap
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "XML/Layering.h"
#include "DD4hep/Shapes.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::IntersectionSolid;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationZ;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

//#define VERBOSE 1

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {
  //unused:  static double tolerance = 0e0;

  xml_det_t   x_det     = element;
  string      det_name    = x_det.nameStr();
  Layering    layering(x_det);

  Material    air         = theDetector.air();
  Material    stavesMaterial    = theDetector.material(x_det.materialStr());

  int           det_id    = x_det.id();
  xml_comp_t    x_staves  = x_det.staves();
  DetElement    sdet      (det_name,det_id);


  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  sens.setType("calorimeter");

  DetElement    module_det("module0",det_id);

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];

  // Some verbose output
  cout << " \n\n\n CREATE DETECTOR: Hcal_EndcapRing_SD_v01" << endl;



//====================================================================
//
// Read all the constant from ILD_....xml
// Use them to build HcalEndcapRing
//
//====================================================================
  // The way to read constants from XML/Detector file.
//  double      Hcal_cells_size                  = theDetector.constant<double>("HcalSD_cells_size");
  double      Hcal_radiator_thickness          = theDetector.constant<double>("HcalSD_radiator_thickness");
  double      Hcal_back_plate_thickness        = theDetector.constant<double>("HcalSD_back_plate_thickness");
  double      Hcal_lateral_plate_thickness     = theDetector.constant<double>("HcalSD_lateral_structure_thickness");
  double      Hcal_stave_gaps                  = theDetector.constant<double>("HcalSD_stave_gaps");
  int         HcalEndcapRing_nlayers           = theDetector.constant<int>("HcalEndcapRingSD_nlayers");
  int         HcalEndcapRing_inner_symmetry    = theDetector.constant<int>("HcalEndcapRingSD_inner_symmetry");
//  int         HcalEndcapRing_outer_symmetry    = theDetector.constant<int>("HcalEndcapRingSD_outer_symmetry");
  double      HcalEndcapRing_inner_radius      = theDetector.constant<double>("HcalEndcapRing_inner_radius");
  double      HcalEndcapRing_outer_radius      = theDetector.constant<double>("HcalEndcapRing_outer_radius");
  double      HcalEndcapRing_min_z             = theDetector.constant<double>("HcalEndcapRing_min_z");
  double      HcalEndcapRing_max_z             = theDetector.constant<double>("HcalEndcapRing_max_z");
//====================================================================
//
// general calculated parameters
//
//====================================================================

  double layer_thickness = 0.0;
  xml_coll_t c(x_det,_U(layer));
  xml_comp_t   x_layer = c;

  for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
     xml_comp_t x_slice = k;
     layer_thickness += x_slice.thickness();
    }

  double Hcal_total_dim_z   = HcalEndcapRing_nlayers * layer_thickness + Hcal_back_plate_thickness;
//  double Hcal_total_dim_z   = HcalEndcapRing_nlayers * layer_thickness;


// ========= Create Hcal end cap ring   ====================================
//  It will be the volume for placing the Hcal endcaps Ring Layers.
//  Itself will be placed into the world volume.
// ==========================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~              EndcapRings                          ~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  double pRMax, pDz, pRMin;
  pRMax = HcalEndcapRing_outer_radius;

  double SpaceForLayers = HcalEndcapRing_max_z -HcalEndcapRing_min_z
                          - Hcal_back_plate_thickness;

  int MaxNumberOfLayers = (int) (SpaceForLayers / layer_thickness);
  int number_of_chambers = HcalEndcapRing_nlayers;
  if(MaxNumberOfLayers < number_of_chambers)
    number_of_chambers = MaxNumberOfLayers;

  cout<<"  cell_sizeX= "<< cell_sizeX << "  cell_sizeY= "<< cell_sizeY << endl;
  cout<<"  layer_thickness (from slices) = "<<layer_thickness<<endl;
  cout<<"  Hcal_radiator_thickness: "<< Hcal_radiator_thickness <<endl;
  cout<<"  Hcal_back_plate_thickness: "<< Hcal_back_plate_thickness <<endl;
  cout<<"  SpaceForLayers: Available "<< SpaceForLayers <<" Requested: "<<Hcal_total_dim_z<<endl;
  cout<<"  Possible / Requested #of layers= "<< MaxNumberOfLayers << "/"<<HcalEndcapRing_nlayers<<endl;
  cout<<"  Remaning Space (free):  "<<  SpaceForLayers - MaxNumberOfLayers * layer_thickness <<endl;
  cout<<"  HcalEndcapRing_min_z:  "<< HcalEndcapRing_min_z <<endl;
  cout<<"  HcalEndcapRing_max_z:  "<< HcalEndcapRing_max_z <<endl;
  cout<<"  HcalEndcapRing_inner_radius: "<< HcalEndcapRing_inner_radius <<endl;
  cout<<"  HcalEndcapRing_outer_radius: "<< HcalEndcapRing_outer_radius <<endl;

  double Dzfree = SpaceForLayers - MaxNumberOfLayers*layer_thickness;
 
// remaining free space "Dzfree" is put at the start of the detector, before the 1st layer
// --> prefer the space filled by air than by steel
  HcalEndcapRing_min_z = HcalEndcapRing_min_z + Dzfree ;

  pDz = (MaxNumberOfLayers*layer_thickness + Hcal_back_plate_thickness)/2.;

  pRMin =HcalEndcapRing_inner_radius;

  double zlen = pDz*2.;

  double rmin = pRMin;
  double rmax = pRMax;
  int numSide = HcalEndcapRing_inner_symmetry;

  PolyhedraRegular HcalEndCapRingSolid( numSide, M_PI/numSide, rmin, rmax,  zlen);

  Volume  HcalEndCapRingLogical("HcalEndCapRingLogical",HcalEndCapRingSolid, stavesMaterial);


  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = numSide  ;
  caloData->outer_symmetry = numSide  ;
  caloData->phi0 = 0 ; // hardcoded 

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rmin ;
  caloData->extent[1] = rmax ;
  caloData->extent[2] = HcalEndcapRing_min_z ;
  caloData->extent[3] = HcalEndcapRing_max_z ;


  //==============================================================
  // build the layer and place into the Hcal EndcapRing
  //==============================================================



  double lpRMax, lpDz, lpRMin;

  lpRMax = rmax - Hcal_lateral_plate_thickness;

  lpDz = layer_thickness / 2.;

  lpRMin = rmin + Hcal_lateral_plate_thickness;

  // G4Polyhedra Envelope parameters
  double phiStart = M_PI/numSide;

  double lzlen = lpDz*2.;

  PolyhedraRegular HcalEndCapRingChamberSolid( numSide, phiStart, lpRMin, lpRMax,  lzlen);

//  Box IntersectionStaveBox( lpRMax, lpRMax, Hcal_total_dim_z);
  Box IntersectionStaveBox( lpRMax, lpRMax, zlen);

  // set up the translation and rotation for the intersection process 
  // this happens in the mother volume coordinate system, so a coordinate transformation is needed
  Position IntersectXYZtrans((pRMax + (Hcal_stave_gaps/2.)), 
			     (pRMax + (Hcal_stave_gaps/2.)),
			     (zlen/2.));

  RotationZYX rot(0.,0.,0.);
  Transform3D tran3D(rot,IntersectXYZtrans);
  // intersect the octagonal layer with a square to get only one quadrant
  IntersectionSolid  HcalEndCapRingStaveSolid( HcalEndCapRingChamberSolid, IntersectionStaveBox, tran3D); 

  int EC_Number_of_towers = 0; 



  // chamber placements


  //place the four staves in their right positions
  for (int stave_id = 0;
       stave_id < 4;
       stave_id++)
    {
      
      for (int layer_id = 1;
	   layer_id <= number_of_chambers;
	   layer_id++)
	{
// Insert layers from zmin -> zmax  BUT! be sure to end the detector with  "last_layer + back_plate"  ... and no free space
	  double Zoff = -zlen/2.0 + (layer_id - 0.5) *layer_thickness ;
          if(stave_id==0)
            cout<<"Ring:layerID, Zoff, staveID= "<<layer_id<<" "<<Zoff<<" "<<stave_id<<endl;
          
	  //====================================================================
	  // Create Hcal EndcapRing Chamber without radiator
	  // Place into the Hcal EndcapRing logical, after each radiator 
	  //====================================================================
	  string layer_name      = det_name+_toString(layer_id,"_layer%d");
	  
	  Volume HcalEndCapRingStaveLogical("HcalEndCapRingStaveLogical",HcalEndCapRingStaveSolid, air);

	  LayeredCalorimeterData::Layer caloLayer ;
	  caloLayer.cellSize0 = cell_sizeX;
	  caloLayer.cellSize1 = cell_sizeY;

	  // Create the slices (sublayers) within the Hcal Barrel Chamber.
	  double slice_pos_z = -layer_thickness/2.;
	  int slice_number = 0;

	  double nRadiationLengths=0.;
	  double nInteractionLengths=0.;
	  double thickness_sum=0;

	  nRadiationLengths   = Hcal_radiator_thickness/(stavesMaterial.radLength());
	  nInteractionLengths = Hcal_radiator_thickness/(stavesMaterial.intLength());

	  for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
	    xml_comp_t x_slice = k;
	    string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
	    double   slice_thickness = x_slice.thickness();
	    Material slice_material  = theDetector.material(x_slice.materialStr());
	    DetElement slice(layer_name,_toString(slice_number,"slice%d"),x_det.id());
	    
	    slice_pos_z += slice_thickness/2.;
            if(stave_id==0 && layer_id==1)
              cout<<"Layer_slice: "<<slice_name<<" slice_thickness: "<<slice_thickness<<" slice_position_z: "<<slice_pos_z<<endl;
	    
	    // Slice volume
	    PolyhedraRegular slicePolyhedraRegularSolid( numSide, phiStart, lpRMin, lpRMax,  slice_thickness);
	    
//	    Box sliceIntersectionStaveBox( lpRMax, lpRMax, Hcal_total_dim_z);
	    Box sliceIntersectionStaveBox( lpRMax, lpRMax, zlen);
	    
	    // set up the translation and rotation for the intersection process 
	    // this happens in the mother volume coordinate system, so a coordinate transformation is needed
	    Position sIntersectXYZtrans((pRMax + (Hcal_stave_gaps/2.)), 
					(pRMax + (Hcal_stave_gaps/2.)),
					(zlen/2.));
//					(Hcal_total_dim_z/2.));
	    
	    RotationZYX srot(0.,0.,0.);
	    Transform3D stran3D(srot,sIntersectXYZtrans);
	    // intersect the octagonal layer with a square to get only one quadrant
	    IntersectionSolid  slice_Solid( slicePolyhedraRegularSolid, sliceIntersectionStaveBox, stran3D); 
	    
	    Volume slice_vol(slice_name,slice_Solid,slice_material);
	    
	    nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
	    nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
	    thickness_sum       += slice_thickness/2;

	    if ( x_slice.isSensitive() ) {
	      slice_vol.setSensitiveDetector(sens);

#if DD4HEP_VERSION_GE( 0, 15 )
	      //Store "inner" quantities
	      caloLayer.inner_nRadiationLengths = nRadiationLengths;
	      caloLayer.inner_nInteractionLengths = nInteractionLengths;
	      caloLayer.inner_thickness = thickness_sum;
	      //Store gas layer thickness
	      caloLayer.sensitive_thickness = slice_thickness;
#endif
	      //Reset counters to measure "outside" quantitites
	      nRadiationLengths=0.;
	      nInteractionLengths=0.;
	      thickness_sum = 0.;
	    }

	    nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
	    nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
	    thickness_sum += slice_thickness/2;

	    // Set region, limitset, and vis.
	    slice_vol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	    // slice PlacedVolume
	    PlacedVolume slice_phv = HcalEndCapRingStaveLogical.placeVolume(slice_vol,Position(0.,0.,slice_pos_z));
	    if ( x_slice.isSensitive() ) {
	      slice_phv.addPhysVolID("layer",layer_id);
	    }
	    
	    slice.setPlacement(slice_phv);
	    // Increment x position for next slice.
	    slice_pos_z += slice_thickness/2.;
	    // Increment slice number.
	    ++slice_number;             
	  }
	  
#if DD4HEP_VERSION_GE( 0, 15 )
	  //Store "outer" quantities
	  caloLayer.outer_nRadiationLengths = nRadiationLengths;
	  caloLayer.outer_nInteractionLengths = nInteractionLengths;
	  caloLayer.outer_thickness = thickness_sum;
#endif  
	  
	  string l_name = _toString(layer_id,"layer%d");
	  string stave_name = _toString(stave_id,"stave%d");
	  DetElement layer(module_det, l_name+stave_name, det_id);
	  
	  double angle_module = M_PI/2. *  stave_id ;
	  
	  Position l_pos(0., 0., Zoff);
	  RotationZ lrotz(angle_module);
	  Position l_new = lrotz*l_pos;  
	  
	  RotationZYX lrot(angle_module,0,0);
	  Transform3D ltran3D(lrot,l_new);
	      
	  PlacedVolume layer_phv = HcalEndCapRingLogical.placeVolume(HcalEndCapRingStaveLogical,ltran3D);
	  layer_phv.addPhysVolID("tower", EC_Number_of_towers);
	  layer_phv.addPhysVolID("stave", stave_id);
	  layer.setPlacement(layer_phv);


	  //-----------------------------------------------------------------------------------------
	  if ( caloData->layers.size() < (unsigned int)number_of_chambers ) {

	    // distance to the surface of the radiator.
	    caloLayer.distance = HcalEndcapRing_min_z + (layer_id - 1)* layer_thickness  ;
//	    caloLayer.thickness = caloLayer.inner_thickness +caloLayer.outer_thickness ;
	    caloLayer.absorberThickness = Hcal_radiator_thickness ;
	    
	    caloData->layers.push_back( caloLayer ) ;
	  }
	  //-----------------------------------------------------------------------------------------

	}

    }  

  // Set stave visualization.
  if (x_staves)   {
    HcalEndCapRingLogical.setVisAttributes(theDetector.visAttributes(x_staves.visStr()));
   }

  //====================================================================
  // Place Hcal Endcap Ring module into the assembly envelope volume
  //====================================================================
  
  double endcap_z_offset = HcalEndcapRing_min_z + pDz;
  
  for(int module_num=0;module_num<2;module_num++) {

    int module_id = ( module_num == 0 ) ? 0:6;
    double this_module_z_offset = ( module_id == 0 ) ? -endcap_z_offset : endcap_z_offset; 
    double this_module_rotY = ( module_id == 0 ) ? M_PI:0.0; 
    Position mxyzVec(0,0,this_module_z_offset);
    RotationZYX mrot(0,this_module_rotY,0);
    Rotation3D mrot3D(mrot);
    Transform3D mtran3D(mrot3D,mxyzVec);

    PlacedVolume pv = envelope.placeVolume(HcalEndCapRingLogical,mtran3D);
    pv.addPhysVolID("module",module_id); // z: +/-

    DetElement sd = (module_num==0) ? module_det : module_det.clone(_toString(module_num,"module%d"));
    sd.setPlacement(pv);

  }
  
  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;

  return sdet;
  
}



DECLARE_DETELEMENT(Hcal_EndcapRing_SD_v01, create_detector)

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
 
   Tibor Kurca: Modify for SemiDigital Hcal_Endcaps &
                - calculate chamber thickness from slice-thicknesses values in xml 
                - currently all volume filled with layers, not corresponding to an engineering design 
                - inner, outer symmetries should be also different ????!!!
                - number of staves should be defined also by parameter in xml file????
                
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
using dd4hep::SubtractionSolid;
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
  string      det_type    = x_det.typeStr();
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
  DetElement    stave_det("module0stave0",det_id);

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];


//====================================================================
//
// Read all the constants from hcal_defs.xml
// Use them to build HcalEndcaps
//
//====================================================================
  // The way to read constant from XML/THEDETECTOR file.
  double      Hcal_radiator_thickness          = theDetector.constant<double>("HcalSD_radiator_thickness");
  double      HcalEndcap_nlayers               = theDetector.constant<double>("HcalEndcapSD_nlayers");

  double      HcalEndcap_inner_radius          = theDetector.constant<double>("HcalEndcap_inner_radius");
  double      HcalEndcap_outer_radius          = theDetector.constant<double>("HcalEndcap_outer_radius");
  double      HcalEndcap_min_z                 = theDetector.constant<double>("HcalEndcap_min_z");
  double      HcalEndcap_max_z                 = theDetector.constant<double>("HcalEndcap_max_z");

  double      Hcal_stave_gaps                  = theDetector.constant<double>("HcalSD_stave_gaps");
  double      Hcal_lateral_plate_thickness     = theDetector.constant<double>("HcalSD_endcap_lateral_structure_thickness");
  double      Hcal_back_plate_thickness        = theDetector.constant<double>("HcalSD_back_plate_thickness");
  int         HcalEndcap_symmetry              = theDetector.constant<int>("HcalEndcap_symmetry");

// Some verbose output
  cout << " \n\n\n CREATE DETECTOR: Hcal_Endcaps_SD_v01" << endl;
  cout << " \n\n\n Detector Type " << det_type<< endl;
  cout<<"  cell_sizeX, cell_sizeY:  "<<cell_sizeX <<" "<<cell_sizeY <<endl;
  cout<<"  HcalEndcap_inner_radius: "<< HcalEndcap_inner_radius <<endl;
  cout<<"  HcalEndcap_outer_radius: "<< HcalEndcap_outer_radius <<endl;
  cout<<"  HcalEndcap_min_z:  "<< HcalEndcap_min_z <<endl;
  cout<<"  HcalEndcap_max_z:  "<< HcalEndcap_max_z <<endl;
  cout<<"  Hcal_stave_gaps:  "<< Hcal_stave_gaps <<endl;

// First get the layer (chamber) thickness caluculated from the Hcal_Endcaps_SD_v0x.xml file
// not from the fixed, defined value "Hcal_chamber_thickness" 
  xml_coll_t c(x_det,_U(layer));
  xml_comp_t   x_layer = c;

  double layer_thickness = 0.0;

  for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
    xml_comp_t x_slice = k;
    layer_thickness += x_slice.thickness();
  }
  cout<<" layer_thickness (from slices) = "<<layer_thickness<<endl; 
//====================================================================
//
// general calculated parameters
//
//====================================================================


// ========= Create Hcal end cap ring   ====================================
//  It will be the volume for placing the Hcal endcaps Ring Layers.
//  And the absorber plate.
//  Itself will be placed into the world volume.
// ==========================================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~              Endcaps                          ~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  double pRMax, pDz, pRMin;
  pRMax = HcalEndcap_outer_radius;
  pRMin = HcalEndcap_inner_radius;

  double SpaceForLayers = HcalEndcap_max_z -HcalEndcap_min_z
                        - Hcal_back_plate_thickness;

  pDz = (HcalEndcap_max_z -HcalEndcap_min_z) / 2.;


  double zlen = pDz*2.;

  double rmin = pRMin;
  double rmax = pRMax;
  int numSide = HcalEndcap_symmetry;
  //double phiStart = M_PI/numSide;
  double phiStart = 0.0;

  PolyhedraRegular HcalEndCapSolidP( numSide, phiStart, rmin, rmax,  zlen);
  Box CenterHcalBox(rmin, rmin, zlen);

  SubtractionSolid HcalEndCapSolid(HcalEndCapSolidP, CenterHcalBox);

  Volume  HcalEndCapLogical("HcalEndCapLogical",HcalEndCapSolid, stavesMaterial);


  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = 4  ;
  caloData->outer_symmetry = 0  ;
//  caloData->outer_symmetry = numSide  ;
  caloData->phi0 = 0 ; // hardcoded 

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rmin ;
  caloData->extent[1] = rmax ;
  caloData->extent[2] = HcalEndcap_min_z ;
  caloData->extent[3] = HcalEndcap_max_z ;


  //==============================================================
  // build the layer and place into the Hcal Endcap
  //==============================================================



  double lpRMax, lpDz, lpRMin;

  lpRMax = rmax - Hcal_lateral_plate_thickness;

  lpDz = layer_thickness / 2.;

  lpRMin = rmin + Hcal_lateral_plate_thickness;

  // G4Polyhedra Envelope parameters

  double lzlen = lpDz*2.;

  PolyhedraRegular HcalEndCapChamberSolidC( numSide, phiStart, lpRMin, lpRMax,  lzlen);
  
  SubtractionSolid HcalEndCapChamberSolid(HcalEndCapChamberSolidC, CenterHcalBox);
  Box IntersectionStaveBox( lpRMax, lpRMax, zlen);

  // set up the translation and rotation for the intersection process 
  // this happens in the mother volume coordinate system, so a coordinate transformation is needed
  Position IntersectXYZtrans((pRMax + (Hcal_stave_gaps/2.)), 
			     (pRMax + (Hcal_stave_gaps/2.)),
			     (zlen/2.));

  RotationZYX rot(0.,0.,0.);
  Transform3D tran3D(rot,IntersectXYZtrans);
  // intersect the octagonal layer with a square to get only one quadrant
  IntersectionSolid  HcalEndCapStaveSolid( HcalEndCapChamberSolid, IntersectionStaveBox, tran3D); 

  int EC_Number_of_towers = 0; 



  // chamber placements
  int number_of_chambers = HcalEndcap_nlayers;
  int possible_number_of_chambers = (int) floor( SpaceForLayers/ layer_thickness );
  if(possible_number_of_chambers < number_of_chambers)
    number_of_chambers = possible_number_of_chambers;
   
  cout<<"  Possible/ Requested #of chambers= "<<possible_number_of_chambers<<" / "<<HcalEndcap_nlayers<<endl; 
  //place the four staves in their right positions
  for (int stave_id = 0;
       stave_id < 4;
       stave_id++)
    {
      
      for (int layer_id = 1;
	   layer_id <= number_of_chambers;
	   layer_id++)
	{
	  double Zoff = -zlen/2.0 + (layer_id - 0.5)* layer_thickness ; 
	  
	  //====================================================================
	  // Create Hcal Endcap Chamber without radiator ...!! TK  with radiator defined as slice in xml 
	  // Place into the Hcal Endcap logical
	  //====================================================================
	  
          string layer_name      = det_name+_toString(layer_id,"_layer%d");
	  Volume HcalEndCapStaveLogical("HcalEndCapStaveLogical",HcalEndCapStaveSolid, air);

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
            if(stave_id==1 && layer_id==1)
	      cout<<"  Layer_slice:  "<<  slice_name<<" slice_thickness:  "<< slice_thickness<< endl;
	    DetElement slice(layer_name,_toString(slice_number,"slice%d"),x_det.id());
	    
	    slice_pos_z += slice_thickness/2.;
	    
	    // Slice volume
	    PolyhedraRegular slicePolyhedraRegularSolidS( numSide, phiStart, lpRMin, lpRMax,  slice_thickness);
            SubtractionSolid slicePolyhedraRegularSolid(slicePolyhedraRegularSolidS, CenterHcalBox);
	    
	    Box sliceIntersectionStaveBox( lpRMax, lpRMax, zlen);
	    
	    // set up the translation and rotation for the intersection process 
	    // this happens in the mother volume coordinate system, so a coordinate transformation is needed
	    Position sIntersectXYZtrans((pRMax + (Hcal_stave_gaps/2.)), 
					(pRMax + (Hcal_stave_gaps/2.)),
					(zlen/2.));
	    
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
	    PlacedVolume slice_phv = HcalEndCapStaveLogical.placeVolume(slice_vol,Position(0.,0.,slice_pos_z));
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
	  
// git +	  double angle_module = -M_PI/2. *  stave_id;
	  double angle_module = M_PI/2. *  stave_id;
	  
	  Position l_pos(0., 0., Zoff);
	  RotationZ lrotz(angle_module);
	  Position l_new = lrotz*l_pos;  
	  
	  RotationZYX lrot(angle_module,0,0);
	  Transform3D ltran3D(lrot,l_new);
	      
	  PlacedVolume layer_phv = HcalEndCapLogical.placeVolume(HcalEndCapStaveLogical,ltran3D);
	  layer_phv.addPhysVolID("tower", EC_Number_of_towers);
	  layer_phv.addPhysVolID("stave", stave_id);
	  layer.setPlacement(layer_phv);


	  //-----------------------------------------------------------------------------------------
	  if ( caloData->layers.size() < (unsigned int)number_of_chambers ) {

	    // distance to the surface of the radiator.
	    caloLayer.distance = HcalEndcap_min_z + (layer_id - 1)* layer_thickness  ;
//	    caloLayer.thickness = caloLayer.inner_thickness +caloLayer.outer_thickness ;
	    caloLayer.absorberThickness = Hcal_radiator_thickness ;
	    
	    caloData->layers.push_back( caloLayer ) ;
	  }
	  //-----------------------------------------------------------------------------------------

	}

    }  

  // Set stave visualization.
  if (x_staves)   {
    HcalEndCapLogical.setVisAttributes(theDetector.visAttributes(x_staves.visStr()));
   }

  //====================================================================
  // Place Hcal Endcap module into the assembly envelope volume
  //====================================================================
  
  double endcap_z_offset = HcalEndcap_min_z + pDz;
  
  for(int module_num=0;module_num<2;module_num++) {

    int module_id = ( module_num == 0 ) ? 0:6;
    double this_module_z_offset = ( module_id == 0 ) ? -endcap_z_offset : endcap_z_offset; 
    double this_module_rotY = ( module_id == 0 ) ? M_PI:0.0; 
    cout <<"  Hcal_Endcaps: module_id, module_num, z_offset, roty=  "<< module_id <<" "<<module_num<<" "<<this_module_z_offset<< " "<<this_module_rotY <<endl;
  
    Position mxyzVec(0,0,this_module_z_offset);
    RotationZYX mrot(0,this_module_rotY,0);
    Rotation3D mrot3D(mrot);
    Transform3D mtran3D(mrot3D,mxyzVec);

    PlacedVolume pv = envelope.placeVolume(HcalEndCapLogical,mtran3D);
    pv.addPhysVolID("module",module_id); // z: +/-

    DetElement sd = (module_num==0) ? module_det : module_det.clone(_toString(module_num,"module%d"));
    sd.setPlacement(pv);

  }
  
  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;

  return sdet;
  
}



DECLARE_DETELEMENT(Hcal_Endcaps_SD_v01, create_detector)

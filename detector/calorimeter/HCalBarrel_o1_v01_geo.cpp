//====================================================================
//  AIDA Detector description implementation for LCD
//--------------------------------------------------------------------
//
//  Author     : N. Nikiforou
//               Adapted from PolyhedraBarrelCalorimeter by M. Frank
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Layer;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
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

static void placeStaves(DetElement& parent, DetElement& stave, double rmin, int numsides, double total_thickness,
			Volume envelopeVolume, double innerAngle, Volume sectVolume) {
  double innerRotation = innerAngle;
  double offsetRotation = -innerRotation / 2;
  double sectCenterRadius = rmin + total_thickness / 2;
  double rotX = M_PI / 2;
  double rotY = -offsetRotation;
  double posX = -sectCenterRadius * std::sin(rotY);
  double posY = sectCenterRadius * std::cos(rotY);

  for (int module = 1; module <= numsides; ++module) {
    DetElement det = module > 1 ? stave.clone(_toString(module,"stave%d")) : stave;
    Transform3D trafo(RotationZYX(0, rotY, rotX), Translation3D(-posX, -posY, 0));
    PlacedVolume pv = envelopeVolume.placeVolume(sectVolume,trafo);
    // Not a valid volID: pv.addPhysVolID("stave", 0);
    pv.addPhysVolID("module", module);
    det.setPlacement(pv);
    parent.add(det);
    rotY -= innerRotation;
    posX = -sectCenterRadius * std::sin(rotY);
    posY = sectCenterRadius * std::cos(rotY);
  }
}

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  Layering layering(x_det);
  xml_comp_t staves = x_det.staves();
  xml_dim_t dim = x_det.dimensions();
  string det_name = x_det.nameStr();
  string det_type = x_det.typeStr();
  Material air = theDetector.air();
  double totalThickness = layering.totalThickness();
  int totalRepeat = 0;
  int totalSlices = 0;
  double gap = xml_dim_t(x_det).gap();
  int nsides = dim.numsides();

  double detZ = dim.z();
  double rmin = dim.rmin();
  DetElement sdet(det_name, x_det.id());
  DetElement stave("stave1", x_det.id());
  //  Volume motherVol = theDetector.pickMotherVolume(sdet);
  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];
  
  //Create caloData object to extend driver with data required for reconstruction
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::BarrelLayout ;
  caloData->inner_symmetry = nsides;
  caloData->outer_symmetry = nsides; 
  
  /** NOTE: phi0=0 means lower face flat parallel to experimental floor
   *  This is achieved by rotating the modules with respect to the envelope
   *  which is assumed to be a Polyhedron and has its axes rotated with respect
   *  to the world by 180/nsides. In any other case (e.g. if you want to have
   *  a tip of the calorimeter touching the ground) this value needs to be computed
   */
  
  caloData->inner_phi0 = 0.; 
  caloData->outer_phi0 = 0.; 
  caloData->gap0 = 0.; //FIXME
  caloData->gap1 = 0.; //FIXME
  caloData->gap2 = 0.; //FIXME  
  
  
  for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
    xml_comp_t x_layer = c;
    int repeat = x_layer.repeat();
    totalRepeat += repeat;
    totalSlices += x_layer.numChildren(_U(slice));
  }

 // --- create an envelope volume and position it into the world ---------------------

  Volume envelopeVol = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

//  PolyhedraRegular polyhedra(nsides, rmin, rmin + totalThickness, detZ);
//  Volume envelopeVol(det_name, polyhedra, air);

// std::cout<<"!!!!!!!!!!"<<std::setprecision(16)<<rmin + totalThickness<<std::endl;

  // Add the subdetector envelope to the structure.
  double innerAngle = 2 * M_PI / nsides;
  double halfInnerAngle = innerAngle / 2;
  double tan_inner = std::tan(halfInnerAngle) * 2;
  double innerFaceLen = rmin * tan_inner;
  double outerFaceLen = (rmin + totalThickness) * tan_inner;
  double staveThickness = totalThickness;
  
  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rmin ;
  caloData->extent[1] = rmin + totalThickness ;
  caloData->extent[2] = 0. ;
  caloData->extent[3] = detZ/2.0 ;
  
  Trapezoid staveTrdOuter(innerFaceLen / 2, outerFaceLen / 2, detZ / 2, detZ / 2, staveThickness / 2);
  Volume staveOuterVol("stave_outer", staveTrdOuter, air);

  Trapezoid staveTrdInner(innerFaceLen / 2 - gap, outerFaceLen / 2 - gap, detZ / 2, detZ / 2, staveThickness / 2);
  Volume staveInnerVol("stave_inner", staveTrdInner, air);

  double layerOuterAngle = (M_PI - innerAngle) / 2;
  double layerInnerAngle = (M_PI / 2 - layerOuterAngle);
  double layer_pos_z = -(staveThickness / 2);
  double layer_dim_x = innerFaceLen / 2 - gap * 2;
  int layer_num = 1;

  //#### LayeringExtensionImpl* layeringExtension = new LayeringExtensionImpl();
  //#### Position layerNormal(0,0,1);
      std::cout<<"XXX "<<std::setprecision(8);

  for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
    xml_comp_t x_layer = c;
    int repeat = x_layer.repeat();            // Get number of times to repeat this layer.
    const Layer* lay = layering.layer(layer_num - 1); // Get the layer from the layering engine.
    // Loop over repeats for this layer.
    for (int j = 0; j < repeat; j++) {
      string layer_name = _toString(layer_num, "layer%d");
      double layer_thickness = lay->thickness();
      DetElement layer(stave, layer_name, layer_num);
      //### layeringExtension->setLayer(layer_num, layer, layerNormal);
      
      LayeredCalorimeterData::Layer caloLayer ;
      
      
      // Layer position in Z within the stave.
      layer_pos_z += layer_thickness / 2;
      // Layer box & volume
      Volume layer_vol(layer_name, Box(layer_dim_x, detZ / 2, layer_thickness / 2), air);

      // Create the slices (sublayers) within the layer.
      double slice_pos_z = -(layer_thickness / 2);
      int slice_number = 0;
      
      double totalAbsorberThickness=0.;
      double sens_pos= 0.;

      double th_i(0.), th_o(-1.) ;
      for (xml_coll_t k(x_layer, _U(slice)); k; ++k) {
        xml_comp_t x_slice = k;
        string slice_name = _toString(slice_number, "slice%d");
        double slice_thickness = x_slice.thickness();
        Material slice_material = theDetector.material(x_slice.materialStr());
        DetElement slice(layer, slice_name, slice_number);
        
        slice_pos_z += slice_thickness / 2;
        // Slice volume & box
        Volume slice_vol(slice_name, Box(layer_dim_x, detZ / 2, slice_thickness / 2), slice_material);
        
        if( x_slice.isRadiator() ==true)
          totalAbsorberThickness+= slice_thickness;
        
          
        
        if (x_slice.isSensitive()) {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
          sens_pos = slice_pos_z;
	  th_i += slice_thickness / 2. ;
	  th_o  = slice_thickness / 2. ;
	} else {
	  if( th_o < 0. ){
	    th_i += slice_thickness;
	  } else {
	    th_o += slice_thickness;
	  }
        } 
        // Set region, limitset, and vis.
        slice_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
        // slice PlacedVolume
        PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol, Position(0, 0, slice_pos_z));
        
        slice.setPlacement(slice_phv);
        // Increment Z position for next slice.
        slice_pos_z += slice_thickness / 2;
        // Increment slice number.
        ++slice_number;
      }
      // Set region, limitset, and vis.
      layer_vol.setAttributes(theDetector, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

      // Layer physical volume.
      PlacedVolume layer_phv = staveInnerVol.placeVolume(layer_vol, Position(0, 0, layer_pos_z));
      layer_phv.addPhysVolID("layer", layer_num);
      layer.setPlacement(layer_phv);

      caloLayer.distance = rmin+layer_pos_z + staveThickness / 2 + sens_pos; //NEED TO START FROM ORIGIN
      
      std::cout<<" "<< caloLayer.distance/dd4hep::mm;

      caloLayer.inner_thickness = th_i ;
      caloLayer.outer_thickness = th_o ;
      caloLayer.absorberThickness = totalAbsorberThickness;
      caloLayer.cellSize0 = cell_sizeX;
      caloLayer.cellSize1 = cell_sizeY;
      
      caloData->layers.push_back( caloLayer ) ;
      
      // Increment the layer X dimension.
      layer_dim_x += layer_thickness * std::tan(layerInnerAngle);    // * 2;
      // Increment the layer Z position.
      layer_pos_z += layer_thickness / 2;
      // Increment the layer number.
      ++layer_num;
    }
  }
  
  std::cout<<std::endl;

  // Add stave inner physical volume to outer stave volume.
  staveOuterVol.placeVolume(staveInnerVol);
  if ( staves )  {
    // Set the vis attributes of the outer stave section.
    stave.setVisAttributes(theDetector, staves.visStr(), staveInnerVol);
    stave.setVisAttributes(theDetector, staves.visStr(), staveOuterVol);
  }
  // Place the staves.
  placeStaves(sdet, stave, rmin, nsides, totalThickness, envelopeVol, innerAngle, staveOuterVol);
  // Set envelope volume attributes.
  envelopeVol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  

  //FOR NOW, USE A MORE "SIMPLE" VERSION OF EXTENSIONS, INCLUDING NECESSARY GEAR PARAMETERS
  //Copied from Frank's SHcalSc04 Implementation
  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
  
  return sdet;
}

DECLARE_DETELEMENT(HCalBarrel_o1_v01, create_detector)

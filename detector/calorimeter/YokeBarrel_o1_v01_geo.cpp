// $Id: PolyhedraBarrelCalorimeter2_geo.cpp 1390 2014-11-14 16:32:48Z Christian.Grefe@cern.ch $
//====================================================================
//  AIDA Detector description implementation for LCD
//--------------------------------------------------------------------
//
//  Author     : M.Frank
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

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
    pv.addPhysVolID("stave", 0);
    pv.addPhysVolID("module", module);
    det.setPlacement(pv);
    parent.add(det);
    rotY -= innerRotation;
    posX = -sectCenterRadius * std::sin(rotY);
    posY = sectCenterRadius * std::cos(rotY);
  }
}

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  Layering layering(x_det);
  xml_comp_t staves = x_det.staves();
  xml_dim_t dim = x_det.dimensions();
  string det_name = x_det.nameStr();
  string det_type = x_det.typeStr();
  Material air = lcdd.air();
  double totalThickness = layering.totalThickness();
  int totalRepeat = 0;
  int totalSlices = 0;
  double gap = xml_dim_t(x_det).gap();
  int numSides = dim.numsides();
  double detZ = dim.z();
  double rmin = dim.rmin();
  DetElement sdet(det_name, x_det.id());
  DetElement stave("stave1", x_det.id());
  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];  
//  Volume motherVol = lcdd.pickMotherVolume(sdet);

  //Create caloData object to extend driver with data required for reconstruction
  DDRec::LayeredCalorimeterData* caloData = new DDRec::LayeredCalorimeterData ;
  caloData->layoutType = DDRec::LayeredCalorimeterData::BarrelLayout ;
  caloData->inner_symmetry = numSides;
  caloData->outer_symmetry = numSides; 
  caloData->phi0 = 0;
  
  for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
    xml_comp_t x_layer = c;
    int repeat = x_layer.repeat();
    totalRepeat += repeat;
    totalSlices += x_layer.numChildren(_U(slice));
  }

 // --- create an envelope volume and position it into the world ---------------------

  Volume envelopeVol = XML::createPlacedEnvelope( lcdd,  e , sdet ) ;

  if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

//  PolyhedraRegular polyhedra(numSides, rmin, rmin + totalThickness, detZ);
//  Volume envelopeVol(det_name, polyhedra, air);

std::cout<<"!!!!!!!!!!"<<std::setprecision(16)<<rmin + totalThickness<<std::endl;

  // Add the subdetector envelope to the structure.
  double innerAngle = 2 * M_PI / numSides;
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

      DDRec::LayeredCalorimeterData::Layer caloLayer ;
      
      // Layer position in Z within the stave.
      layer_pos_z += layer_thickness / 2;
      // Layer box & volume
      Volume layer_vol(layer_name, Box(layer_dim_x, detZ / 2, layer_thickness / 2), air);

      // Create the slices (sublayers) within the layer.
      double slice_pos_z = -(layer_thickness / 2);
      int slice_number = 1;
      double totalAbsorberThickness=0.;
      int sensor_number=0;
      for (xml_coll_t k(x_layer, _U(slice)); k; ++k) {
        xml_comp_t x_slice = k;
        string slice_name = _toString(slice_number, "slice%d");
        double slice_thickness = x_slice.thickness();
        Material slice_material = lcdd.material(x_slice.materialStr());
        DetElement slice(layer, slice_name, slice_number);

        slice_pos_z += slice_thickness / 2;
        // Slice volume & box
        Volume slice_vol(slice_name, Box(layer_dim_x, detZ / 2, slice_thickness / 2), slice_material);
        char val = x_slice.hasAttr(_U(radiator)) ? x_slice.attr < string > (_U(radiator))[0] : 'f';
        val = std::toupper(val);
        bool isAbsorber =  (val == 'T' || val == 'Y');
        
        if( isAbsorber ==true)
          totalAbsorberThickness+= slice_thickness;
        

        // Set region, limitset, and vis.
        slice_vol.setAttributes(lcdd, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
        // slice PlacedVolume
        PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol, Position(0, 0, slice_pos_z));
        slice.setPlacement(slice_phv);
        
        if (x_slice.isSensitive()) {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
          slice_phv.addPhysVolID("submodule", sensor_number); 
          
          sensor_number++;
        }
        
        // Increment Z position for next slice.
        slice_pos_z += slice_thickness / 2;
        // Increment slice number.
        ++slice_number;
      }
      // Set region, limitset, and vis.
      layer_vol.setAttributes(lcdd, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

      // Layer physical volume.
      PlacedVolume layer_phv = staveInnerVol.placeVolume(layer_vol, Position(0, 0, layer_pos_z));
      layer_phv.addPhysVolID("layer", layer_num);
      layer_phv.addPhysVolID("side", 0); ///FIXME! Should try to set in parent
      
      layer.setPlacement(layer_phv);

      
      caloLayer.distance = layer_pos_z;
      caloLayer.thickness = layer_thickness;
      caloLayer.absorberThickness = totalAbsorberThickness;
      caloLayer.cellSize0 = cell_sizeX;
      caloLayer.cellSize1 = cell_sizeY;
      
      caloData->layers.push_back( caloLayer ) ; //PROBLEM WHEN HAVING TWO ACTIVE ELEMENTS PER LAYER
      
      // Increment the layer X dimension.
      layer_dim_x += layer_thickness * std::tan(layerInnerAngle);    // * 2;
      // Increment the layer Z position.
      layer_pos_z += layer_thickness / 2;
      // Increment the layer number.
      ++layer_num;
    }
  }

  // Add stave inner physical volume to outer stave volume.
  staveOuterVol.placeVolume(staveInnerVol);
  if ( staves )  {
    // Set the vis attributes of the outer stave section.
    stave.setVisAttributes(lcdd, staves.visStr(), staveInnerVol);
    stave.setVisAttributes(lcdd, staves.visStr(), staveOuterVol);
  }
  // Place the staves.
  placeStaves(sdet, stave, rmin, numSides, totalThickness, envelopeVol, innerAngle, staveOuterVol);
  // Set envelope volume attributes.
  envelopeVol.setAttributes(lcdd, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());



  sdet.addExtension< DDRec::LayeredCalorimeterData >( caloData ) ;
  
  
  return sdet;
}

DECLARE_DETELEMENT(YokeBarrel_o1_v01, create_detector)

//====================================================================
//  DD4hep Geometry driver for Sampling Calo BOX prototype
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:  $
//====================================================================
#include "DD4hep/Printout.h"
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "LcgeoExceptions.h"

#include <iostream>
#include <vector>

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Layer;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::SensitiveDetector;
using dd4hep::Volume;
using dd4hep::_toString;

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif


static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {

  xml_det_t   x_det       = element;
  string      det_name    = x_det.nameStr();
  DetElement  sdet( det_name,x_det.id() );

  Layering    layering(x_det);
  xml_dim_t dim = x_det.dimensions();

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;
  
  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  Material      air               = theDetector.air();

  sens.setType("calorimeter");

//====================================================================
//
// Read all the dimensions from compact.xml, user can update the value.
// Use them to build a calo box prototye.
//
//====================================================================

  double      Calo_dim_x          = dim.x();
  double      Calo_dim_y          = dim.y();
  double      Calo_dim_z          = dim.z();

  printout( dd4hep::DEBUG,  "building SamplingCaloBoxPrototype_v01",
	    "Calo_dim_x : %e    Calo_dim_y: %e    Calo_dim_z: %e ",
  	    Calo_dim_x, Calo_dim_y, Calo_dim_z) ;


//====================================================================
//
// general calculated parameters
//
//====================================================================
 
  //calorimeter dimensions
  double cal_hx = Calo_dim_x/2.0;
  double cal_hy = Calo_dim_y/2.0;
  double cal_hz = Calo_dim_z/2.0;


//====================================================================
//
// build sampling layers in the CaloBox
//
//====================================================================

    int layer_num = 0;
    int layerType   = 0;

    double layer_pos_z = - cal_hz;

    for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
        xml_comp_t x_layer = c;
        int repeat = x_layer.repeat();
        const Layer* lay = layering.layer(layer_num);
        double layer_thickness = lay->thickness();
        string layer_type_name   = _toString(layerType,"layerType%d");

        // Loop over repeats for this layer.
        for (int j = 0; j < repeat; j++) {
            string layer_name = _toString(layer_num, "layer%d");
            DetElement layer(layer_name, layer_num);

            // Layer box & volume
            Volume layer_vol(layer_type_name, Box(cal_hx, cal_hy, layer_thickness / 2), air);
            
            // Create the slices (sublayers) within the layer.
            double slice_pos_z = -(layer_thickness / 2);
            int slice_number = 0;
            
            for (xml_coll_t k(x_layer, _U(slice)); k; ++k) {
                xml_comp_t x_slice = k;
                string slice_name = _toString(slice_number, "slice%d");
                double slice_thickness = x_slice.thickness();
                Material slice_material = theDetector.material(x_slice.materialStr());
                DetElement slice(layer, slice_name, slice_number);

                slice_pos_z += slice_thickness / 2;
                // Slice volume & box
                Volume slice_vol(slice_name, Box(cal_hx, cal_hy, slice_thickness / 2), slice_material);
                
                
                if (x_slice.isSensitive()) {
                    sens.setType("calorimeter");
                    slice_vol.setSensitiveDetector(sens);                   
                } 

                
                // Set region, limitset, and vis.
                slice_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
                // slice PlacedVolume
                PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol, Position(0, 0, slice_pos_z));
                slice_phv.addPhysVolID("slice", slice_number);
                slice.setPlacement(slice_phv);

                // Increment Z position for next slice.
                slice_pos_z += slice_thickness / 2;
                // Increment slice number.
                ++slice_number;
            }
            
            
            // Set region, limitset, and vis.
            layer_vol.setAttributes(theDetector, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
           
            // Layer position in Z within the stave.
            layer_pos_z += layer_thickness / 2;
            // Layer physical volume.
            PlacedVolume layer_phv = envelope.placeVolume(layer_vol, Position(0, 0, layer_pos_z));
            layer_phv.addPhysVolID("layer", layer_num);
            layer.setPlacement(layer_phv);
            
            // Increment the layer Z position.
            layer_pos_z += layer_thickness / 2;
            // Increment the layer number.
            ++layer_num;
        }
        
        ++layerType;
    }

 
    return sdet;

}

DECLARE_DETELEMENT(CaloPrototype_v01, create_detector)

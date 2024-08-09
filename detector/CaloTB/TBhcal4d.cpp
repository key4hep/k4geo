//====================================================================
//  DD4hep Geometry driver for TB AHCAL 2015 prototype
//--------------------------------------------------------------------
//  E.Brianne, DESY
//  $Id:  $
//====================================================================
#include "DD4hep/Printout.h"
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"

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
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Segmentation;
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

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;
  
  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  Material air = theDetector.air();

  sens.setType("calorimeter");

  //====================================================================
  //
  // Read all the constant from compact.xml, user can update the value.
  // Use them to build a calo box prototye.
  //
  //====================================================================

  //unused: double HBU_SSF_half_x = theDetector.constant<double>("HBU_SSF_dim_x")/2.0;
  //unused: double HBU_SSF_half_y = theDetector.constant<double>("HBU_SSF_dim_y")/2.0;

  int HCAL_SSF_ncell_x = theDetector.constant<int>("HCAL_SSF_ncell_x");
  int HCAL_SSF_ncell_y = theDetector.constant<int>("HCAL_SSF_ncell_y");

  //unused: double HBU_BL_half_x = theDetector.constant<double>("HBU_BL_dim_x")/2.0;
  //unused: double HBU_BL_half_y = theDetector.constant<double>("HBU_BL_dim_y")/2.0;

  int HCAL_BL_ncell_x = theDetector.constant<int>("HCAL_BL_ncell_x");
  int HCAL_BL_ncell_y = theDetector.constant<int>("HCAL_BL_ncell_y");

  int HCAL_SSF_nlayers = theDetector.constant<int>("HCAL_SSF_nlayers");
  int HCAL_BL_nlayers = theDetector.constant<int>("HCAL_BL_nlayers");
  double Hcal_layer_thickness = theDetector.constant<double>("Hcal_layer_thickness");

  Readout  readout = sens.readout();
  Segmentation seg = readout.segmentation();

  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
  double cell_sizeX = cellSizeVector[0];
  double cell_sizeY = cellSizeVector[1];

  //====================================================================
  //
  // general calculated parameters
  //
  //====================================================================
 
  //calorimeter dimensions
  double cal_SSF_hx = (double) (HCAL_SSF_ncell_x * cell_sizeX)/2.;
  double cal_SSF_hy = (double) (HCAL_SSF_ncell_y * cell_sizeY)/2.;
  double cal_SSF_hz = (double) (Hcal_layer_thickness * HCAL_SSF_nlayers)/2;

  double cal_BL_hx = (double) (HCAL_BL_ncell_x * cell_sizeX)/2.;
  double cal_BL_hy = (double) (HCAL_BL_ncell_y * cell_sizeY)/2.;
  double cal_BL_hz = (double) (Hcal_layer_thickness * HCAL_BL_nlayers)/2;

  //====================================================================
  //
  // Chambers in the CaloBox
  //
  //====================================================================

  int layer_num = 0;
  int layerType = 0;

  double layer_pos_z = - cal_SSF_hz - cal_BL_hz;

  for (xml_coll_t c(x_det, _U(layer)); c; ++c) 
    {
      xml_comp_t x_layer = c;
      int repeat = x_layer.repeat();                // Get number of times to repeat this layer.
      const Layer* lay = layering.layer(layer_num); // Get the layer from the layering engine.
      double layer_thickness = lay->thickness();
      string layer_type_name   = _toString(layerType,"layerType%d");

      if(layer_num < HCAL_SSF_nlayers)
	{
	  // Loop over repeats for this layer.
	  for (int j = 0; j < repeat; j++) 
	    {
	      string layer_name = _toString(layer_num, "layer%d");
	      DetElement layer_SSF(layer_name, layer_num);

	      // Layer box & volume
	      Volume layer_SSF_vol(layer_type_name, Box(cal_BL_hx, cal_BL_hy, layer_thickness / 2), air);
            
	      // Create the slices (sublayers) within the layer.
	      double slice_SSF_pos_z = -(layer_thickness / 2);
	      int slice_number = 0;
            
	      for (xml_coll_t k(x_layer, _U(slice)); k; ++k) 
		{
		  xml_comp_t x_slice = k;
		  string slice_name = _toString(slice_number, "slice%d");
		  double slice_thickness = x_slice.thickness();
		  Material slice_material = theDetector.material(x_slice.materialStr());
                
		  slice_SSF_pos_z += slice_thickness / 2;

		  //Case of absorber make it bigger than the actual layer (*2 for HBU SSF)
		  if(slice_number == 0)
		    {
		      // Slice volume & box
		      Volume slice_SSF_vol(slice_name, Box(cal_BL_hx, cal_BL_hy, slice_thickness / 2), slice_material);

		      // Set region, limitset, and vis.
		      slice_SSF_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
		      // slice PlacedVolume
		      layer_SSF_vol.placeVolume(slice_SSF_vol, Position(0, 0, slice_SSF_pos_z));
		    }
		  else
		    {
		      // Slice volume & box
		      Volume slice_SSF_vol(slice_name, Box(cal_SSF_hx, cal_SSF_hy, slice_thickness / 2), slice_material);

		      if (x_slice.isSensitive()) 
			{
			  sens.setType("calorimeter");
			  slice_SSF_vol.setSensitiveDetector(sens);                   
			} 

		      // Set region, limitset, and vis.
		      slice_SSF_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
		      // slice PlacedVolume
		      layer_SSF_vol.placeVolume(slice_SSF_vol, Position(0, 0, slice_SSF_pos_z));
		    }

		  // Increment Z position for next slice.
		  slice_SSF_pos_z += slice_thickness / 2;
		  // Increment slice number.
		  ++slice_number;
		}
            
	      // Set region, limitset, and vis.
	      layer_SSF_vol.setAttributes(theDetector, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
           
	      // Layer position in Z within the stave.
	      layer_pos_z += layer_thickness / 2;
	      // Layer physical volume.
	      PlacedVolume layer_SSF_phv = envelope.placeVolume(layer_SSF_vol, Position(0, 0, layer_pos_z));
	      //layer_phv.addPhysVolID("layer", layer_num);
	      layer_SSF_phv.addPhysVolID("K", layer_num);
	      layer_SSF.setPlacement(layer_SSF_phv);
            
	      // Increment the layer Z position.
	      layer_pos_z += layer_thickness / 2;
	      // Increment the layer number.
	      ++layer_num;
	    }
	}
      else
	{
	  // Loop over repeats for this layer.
	  for (int j = 0; j < repeat; j++) 
	    {
	      string layer_name = _toString(layer_num, "layer%d");
	      DetElement layer_BL(layer_name, layer_num);

	      // Layer box & volume
	      Volume layer_BL_vol(layer_type_name, Box(cal_BL_hx, cal_BL_hy, layer_thickness / 2), air);
            
	      // Create the slices (sublayers) within the layer.
	      double slice_BL_pos_z = -(layer_thickness / 2);
	      int slice_number = 0;
            
	      for (xml_coll_t k(x_layer, _U(slice)); k; ++k) 
		{
		  xml_comp_t x_slice = k;
		  string slice_name = _toString(slice_number, "slice%d");
		  double slice_thickness = x_slice.thickness();
		  Material slice_material = theDetector.material(x_slice.materialStr());
                
		  slice_BL_pos_z += slice_thickness / 2;

		  // Slice volume & box
		  Volume slice_BL_vol(slice_name, Box(cal_BL_hx, cal_BL_hy, slice_thickness / 2), slice_material);

		  if (x_slice.isSensitive()) 
		    {
		      sens.setType("calorimeter");
		      slice_BL_vol.setSensitiveDetector(sens);                   
		    } 

		  // Set region, limitset, and vis.
		  slice_BL_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
		  // slice PlacedVolume
		  layer_BL_vol.placeVolume(slice_BL_vol, Position(0, 0, slice_BL_pos_z));

		  // Increment Z position for next slice.
		  slice_BL_pos_z += slice_thickness / 2;
		  // Increment slice number.
		  ++slice_number;
		}
            
	      // Set region, limitset, and vis.
	      layer_BL_vol.setAttributes(theDetector, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
           
	      // Layer position in Z within the stave.
	      layer_pos_z += layer_thickness / 2;
	      // Layer physical volume.
	      PlacedVolume layer_BL_phv = envelope.placeVolume(layer_BL_vol, Position(0, 0, layer_pos_z));
	      //layer_phv.addPhysVolID("layer", layer_num);
	      layer_BL_phv.addPhysVolID("K", layer_num);
	      layer_BL.setPlacement(layer_BL_phv);
            
	      // Increment the layer Z position.
	      layer_pos_z += layer_thickness / 2;
	      // Increment the layer number.
	      ++layer_num;
	    }
	}
        
      ++layerType;
    }
 
  return sdet;

}

DECLARE_DETELEMENT(TBhcal4d, create_detector)

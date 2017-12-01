//====================================================================
//  $Id: SteppedMuonBarrel_o2_v02_geo.cpp 1390 2017-11-21 15:32:48Z protopop@cern.ch $
//--------------------------------------------------------------------
//
//  Author : D.Protopopescu
//           Stepped yoke geometry implemented Nov 2017 
//  Author : N. Nikiforou
//           Adapted from PolyhedraBarrelCalorimeter by M. Frank
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

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens) {
    xml_det_t x_det = e;
    Layering layering(x_det);
    //xml_comp_t staves = x_det.staves();
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
    // For stepped geometry we need to redefine rmin in terms of rmax 
    double rmax = dim.rmax();//*std::cos(M_PI/nsides)
    double rmin = rmax - totalThickness;

    std::cout << "*** In XML rmin=" << dim.rmin() << " but using rmin=" << rmin << std::endl;

    DetElement sdet(det_name, x_det.id());
    //DetElement stave("stave0", x_det.id());
    //Volume motherVol = theDetector.pickMotherVolume(sdet);
    Readout readout = sens.readout();
    Segmentation seg = readout.segmentation();
    
    std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
    double cell_sizeX = cellSizeVector[0];
    double cell_sizeY = cellSizeVector[1];
    
    //Create caloData object to extend driver with data required for reconstruction
    LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
    caloData->layoutType = LayeredCalorimeterData::BarrelLayout ;
    caloData->inner_symmetry = nsides;
    caloData->outer_symmetry = nsides; 
    
    /*  NOTE: phi0=0 means lower face flat parallel to experimental floor
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
    
    // CHECK THIS!
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelopeVol = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;
    
    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    
    // Extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
    caloData->extent[0] = rmin ;
    caloData->extent[1] = rmin + totalThickness;
    caloData->extent[2] = 0. ;

    double innerAngle = 2 * M_PI / nsides;
    double tan_inner = std::tan(innerAngle/2.) * 2;
    double layer_pos_z = -(totalThickness / 2);
    double layer_dim_x = rmin * tan_inner/2. - gap * 2;
  
    int layer_num = 0;
    int layerType   = 0;

    for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
        xml_comp_t x_layer = c;
        int repeat = x_layer.repeat();            // Get number of times to repeat this layer.
        const Layer* lay = layering.layer(layer_num); // Get the layer from the layering engine.
        double layer_thickness = lay->thickness();
        string layer_type_name = _toString(layerType, "layerType%d");

	// This needs checked!
	caloData->extent[3] = detZ/2.0 + floor(repeat/2)*totalThickness/repeat;
	
        LayeredCalorimeterData::Layer caloLayer ;
        caloLayer.cellSize0 = cell_sizeX;
        caloLayer.cellSize1 = cell_sizeY;

        // Loop over repeats for this layer.
        for (int j = 0; j < repeat; j++) {
	  
	    // Stepped geometry is implemented like this
	    double layer_dz = detZ/2 + floor((j+1)/2)*layer_thickness;
	    //caloData.extent[3] = layer_dz;

            string layer_name = _toString(layer_num, "layer%d");
            DetElement layer(sdet, layer_name, layer_num);

            // Layer box & volume
            Volume layer_vol(layer_type_name, Box(layer_dim_x, layer_dz, layer_thickness / 2), air);
            
            // Create the slices (sublayers) within the layer.
            double slice_pos_z = -(layer_thickness / 2);
            int slice_number = 0;
            
            double nRadiationLengths=0.;
            double nInteractionLengths=0.;
            double thickness_sum=0;
            
	    int sensor = 0;
            for (xml_coll_t k(x_layer, _U(slice)); k; ++k) {
                xml_comp_t x_slice = k;
                string slice_name = _toString(slice_number, "slice%d");
                double slice_thickness = x_slice.thickness();
                Material slice_material = theDetector.material(x_slice.materialStr());
                
                slice_pos_z += slice_thickness / 2;

                // Slice volume & box
                Volume slice_vol(slice_name, Box(layer_dim_x, layer_dz, slice_thickness / 2), slice_material);
		// Set region, limitset, and vis.
                slice_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
		// slice PlacedVolume
                PlacedVolume slice_pv = layer_vol.placeVolume(slice_vol, Position(0, 0, slice_pos_z));
                
                nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
                nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
                thickness_sum += slice_thickness/2;
                
                if (x_slice.isSensitive()) {
                    sens.setType("calorimeter");
                    slice_vol.setSensitiveDetector(sens);
		    slice_pv.addPhysVolID("submodule", sensor);

#if DD4HEP_VERSION_GE( 0, 15 )
                    //Store "inner" quantities
                    caloLayer.inner_nRadiationLengths = nRadiationLengths;
                    caloLayer.inner_nInteractionLengths = nInteractionLengths;
                    caloLayer.inner_thickness = thickness_sum;

                    //Store scintillator thickness
                    caloLayer.sensitive_thickness = slice_thickness;
#endif
                    //Reset counters to measure "outside" quantitites
                    nRadiationLengths=0.;
                    nInteractionLengths=0.;
                    thickness_sum = 0.;
		    sensor++;
                } 
                
                nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
                nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
                thickness_sum += slice_thickness/2;
                                
                // Increment Z position for next slice.
                slice_pos_z += slice_thickness / 2;
                // Increment slice number.
                ++slice_number;
            }
            
#if DD4HEP_VERSION_GE( 0, 15 )
            //Store "outer" quantities
            caloLayer.outer_nRadiationLengths = nRadiationLengths;
            caloLayer.outer_nInteractionLengths = nInteractionLengths;
            caloLayer.outer_thickness = thickness_sum;
#endif        
            
            // Set region, limitset, and vis.
            layer_vol.setAttributes(theDetector, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

            // Layer position in Z within the stave.
            layer_pos_z += layer_thickness / 2;

            //The rest of the data is constant; only the distance needs to be updated  
            //Store the position up to the inner face of the layer
            caloLayer.distance = rmin + layer_pos_z + totalThickness / 2 - layer_thickness / 2; 
            std::cout<<"MuonBarrel layer: "<<layer_num<<" Rmin: "<<rmin<<" layer_pos_z: "<<layer_pos_z<<" Dist: "<<caloLayer.distance
#if DD4HEP_VERSION_GE( 0, 15 )
		 <<" inner_thickness: "<<caloLayer.inner_thickness<<" outer_thickness: "<<caloLayer.outer_thickness
#endif
		 <<std::endl;

            //Push back a copy to the caloData structure
            caloData->layers.push_back( caloLayer ) ;

	    // Placement directly within the envelope
            double innerRotation = innerAngle;
	    double offsetRotation = -innerRotation / 2;
            double sectCenterRadius = rmin + totalThickness/2 + layer_pos_z;
            double rotX = M_PI / 2;
            double rotY = -offsetRotation;
            double posX = -sectCenterRadius * std::sin(rotY) + 0;
            double posY = sectCenterRadius * std::cos(rotY) + 0;

	    for (int module = 0; module < nsides; ++module) {
	      DetElement mod = module > 0 ? layer.clone(_toString(module, "layer%d")) : layer;
	      Transform3D trafo(RotationZYX(0, rotY, rotX), Translation3D(-posX, -posY, 0));
	      PlacedVolume mod_phv = envelopeVol.placeVolume(layer_vol, trafo);
	      mod_phv.addPhysVolID("layer", layer_num);
	      mod_phv.addPhysVolID("module", module);
	      mod.setPlacement(mod_phv);

	      rotY -= innerRotation;
	      posX = -sectCenterRadius * std::sin(rotY);
	      posY = sectCenterRadius * std::cos(rotY);
	    }
            
            // Increment the layer X dimension.
            layer_dim_x += layer_thickness * std::tan(innerAngle/2.);
            // Increment the layer Z position.
            layer_pos_z += layer_thickness / 2;
            // Increment the layer number.
            ++layer_num;
        }
        
        ++layerType;
    }
    
    // Set envelope volume attributes.
    envelopeVol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
        
    //FOR NOW, USE A MORE "SIMPLE" VERSION OF EXTENSIONS, INCLUDING NECESSARY GEAR PARAMETERS
    //Copied from Frank's SHcalSc04 Implementation
    sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
    
    return sdet;
}

DECLARE_DETELEMENT(SteppedMuonBarrel_o2_v02, create_detector)

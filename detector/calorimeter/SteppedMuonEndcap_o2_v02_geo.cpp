//====================================================================
//  $Id: SteppedMuonEndcap_o2_v02_geo.cpp 1390 2017-11-30 15:32:48Z protopop@cern.ch $
//--------------------------------------------------------------------
//
//  Author : D.Protopopescu
//           Stepped yoke geometry implemented Nov 2017
//  Author : N.Nikiforou
//           Adapted from PolyhedraBarrel Calorimeter by M.Frank
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
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
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens)  {
    xml_det_t   x_det     = e;
    int         det_id    = x_det.id();
    string      det_name  = x_det.nameStr();
    DetElement  sdet      (det_name,det_id);
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelopeVol = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;
    
    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    
    xml_dim_t   dim       = x_det.dimensions();
    Material    air       = theDetector.air();
    int         nsides_inner = dim.nsides_inner();
    int         nsides_outer = dim.nsides_outer();
    double      rmin      = dim.rmin();
    double      rmax      = dim.rmax();//*std::cos(M_PI/nsides_outer); 
    double      zmin      = dim.zmin();

    cout << "Using rmax=" << rmax << " and rmin=" << rmin << " and diff=" << rmax - rmin << endl;
        
    Layering    layering(x_det);
    double      totalThickness = layering.totalThickness();
    Readout     readout = sens.readout();
    Segmentation seg = readout.segmentation();
    
    std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
    double cell_sizeX = cellSizeVector[0];
    double cell_sizeY = cellSizeVector[1];

    int layer_num = 0;
    int layerType = 0;
    double layerZ = zmin;
    
    //Create caloData object to extend driver with data required for reconstruction
    LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
    caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
    caloData->inner_symmetry = nsides_inner;
    caloData->outer_symmetry = nsides_outer; 
    
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
    
    /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
    caloData->extent[0] = rmin ;
    caloData->extent[1] = rmax ; //CHECK THIS!
    caloData->extent[2] = zmin ;
    caloData->extent[3] = zmin + totalThickness;
    
    //endcapVol.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

    for(xml_coll_t c(x_det,_U(layer)); c; ++c)  {
      xml_comp_t  x_layer = c;
      double      layer_thick = layering.layer(layer_num)->thickness();
      string      layer_type_name = _toString(layerType,"layerType%d");
      int         layer_repeat = x_layer.repeat();
      double      layer_rcutout = x_layer.hasAttr(_U(gap)) ? x_layer.gap() : 0;
      
      std::cout << "Number of layers in group " << layerType << " : " << layer_repeat <<std::endl; 
      std::cout << "MuonEndcap total thickness: " << layer_thick*layer_repeat << std::endl;

      for(int j=0; j<layer_repeat; ++j) {
	
	// Stepped geometry is implemented like this
	double layer_rmax = rmax + (j<5 ? j*2*layer_thick - 5*2*layer_thick : 0.0);
	std::cout << "MuonEndcap layer " << j << " with layer_rmax = " << layer_rmax << std::endl;
	
	PolyhedraRegular polyVolume(nsides_outer, rmin+layer_rcutout, layer_rmax, layer_thick);
	Volume layer_vol(layer_type_name, polyVolume, air);

        int slice_num = 0;
        double sliceZ = -layer_thick/2;
        
        //Create a caloLayer struct for this layer type to store copies of in the parent struct
        LayeredCalorimeterData::Layer caloLayer ;
        caloLayer.cellSize0 = cell_sizeX;
        caloLayer.cellSize1 = cell_sizeY; 
        
        double nRadiationLengths=0.;
        double nInteractionLengths=0.;
        double thickness_sum=0;
        
	int sensor = 0;
        for(xml_coll_t s(x_layer,_U(slice)); s; ++s)  {
            xml_comp_t x_slice = s;
            string     slice_name  = _toString(slice_num,"slice%d");
            double     slice_thickness = x_slice.thickness();
            Material   slice_material   = theDetector.material(x_slice.materialStr());
            Volume     slice_vol(slice_name,PolyhedraRegular(nsides_outer,rmin+layer_rcutout,layer_rmax,slice_thickness),slice_material);
            
            slice_vol.setVisAttributes(theDetector.visAttributes(x_slice.visStr()));
            sliceZ += slice_thickness/2;
            PlacedVolume  slice_pv = layer_vol.placeVolume(slice_vol,Position(0,0,sliceZ));
            
            nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
            nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
            thickness_sum += slice_thickness/2;
            
            if ( x_slice.isSensitive() )  {
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

            sliceZ += slice_thickness/2;
            slice_num++;
        }
        
#if DD4HEP_VERSION_GE( 0, 15 )
        //Store "outer" quantities
        caloLayer.outer_nRadiationLengths = nRadiationLengths;
        caloLayer.outer_nInteractionLengths = nInteractionLengths;
        caloLayer.outer_thickness = thickness_sum;
#endif        
	// Set region, limitset, and vis.
	layer_vol.setAttributes(theDetector, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
        
        if (layer_repeat <= 0) 
	  throw std::runtime_error(x_det.nameStr()+"> Invalid repeat value");
            
	//Store the position up to the inner face of the layer
	caloLayer.distance = layerZ;

	//Push back a copy to the caloData structure
	caloData->layers.push_back( caloLayer );

	string phys_lay = _toString(layer_num, "layer%d");
	double z_pos = layerZ + layer_thick/2;

	DetElement    layer_sideA(sdet, phys_lay, layer_num);
	PlacedVolume  pv = envelopeVol.placeVolume(layer_vol, Position(0, 0, z_pos));
	pv.addPhysVolID("layer", layer_num);
	pv.addPhysVolID("side", 1);
	layer_sideA.setPlacement(pv);

	// Reflect it
	DetElement layer_sideB = layer_sideA.clone("layer_sideB");
	pv = envelopeVol.placeVolume(layer_vol, Transform3D(RotationZYX(0,M_PI,0),
							 Position(0,0,-z_pos)));
	pv.addPhysVolID("side", -1);
	pv.addPhysVolID("layer", layer_num);
	layer_sideB.setPlacement(pv);

	layerZ += layer_thick;
	++layer_num;
      }      
      ++layerType;    
    }

    // Set envelope volume attributes.
    envelopeVol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

    sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
    
    return sdet;    
}

DECLARE_DETELEMENT(SteppedMuonEndcap_o2_v02,create_detector)


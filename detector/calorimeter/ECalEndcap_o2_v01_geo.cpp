
//====================================================================
//  Modified ECal Endcap Driver for the CLIC detector
//--------------------------------------------------------------------
//
//  Author     : N.Nikiforou
//  Adapted from PolyhedraBarrel Calorimeter by M. Frank
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
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens) {

  xml_det_t   x_det     = e;
  int         det_id    = x_det.id();
  string      det_name  = x_det.nameStr();
  DetElement    sdet      (det_name,det_id);

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
  
  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
  
  //-----------------------------------------------------------------------------------
  

    
    
//   std::cout<<"Building ECal EndCap inside envelope."<<std::endl;
  xml_dim_t   dim       = x_det.dimensions();
  Material    air       = theDetector.air();
  int         nsides_inner = dim.nsides_inner();
  int         nsides_outer = dim.nsides_outer();
  double      rmin      = dim.rmin();
  double      rmax      = dim.rmax(); /// FIXME: IS THIS RIGHT?
  double      zmin      = dim.zmin();
  
  
  Layering    layering(x_det);
  double      totalThickness = layering.totalThickness();
  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];
  
  PolyhedraRegular polyVolume(nsides_outer,rmin,rmax,totalThickness);

  
  Volume      endcapVol("endcap",polyVolume,air);
  
  DetElement  endcapA(sdet,"endcap",det_id);
  Ref_t(endcapA)->SetName((det_name+"_A").c_str());
  
  int l_num = 1;
  int layerType   = 0;
  double layerZ   = -totalThickness/2;
  
  //Create caloData object to extend driver with data required for reconstruction
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = nsides_inner;
  caloData->outer_symmetry = nsides_outer; 
  
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
  
  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rmin ;
  caloData->extent[1] = rmax ; ///FIXME: CHECK WHAT IS NEEDED (EXSCRIBED?)
  caloData->extent[2] = zmin ;
  caloData->extent[3] = zmin + totalThickness;
  
  endcapVol.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  
  for(xml_coll_t c(x_det,_U(layer)); c; ++c)  {
    xml_comp_t       x_layer  = c;
    double           l_thick  = layering.layer(l_num-1)->thickness();
    string           l_name   = _toString(layerType,"layer%d");
    int              l_repeat = x_layer.repeat();
    
    std::cout<<"Number of layers in group "<<layerType<<" : "<<l_repeat<<std::endl; 
    
    Volume           l_vol(l_name,PolyhedraRegular(nsides_outer,rmin,rmax,l_thick),air);
    
    vector<PlacedVolume> sensitives;
    
    int s_num = 1;
    double sliceZ = -l_thick/2;
    double totalAbsorberThickness=0.;

    double th_i(0.), th_o(-1.) ;
    for(xml_coll_t s(x_layer,_U(slice)); s; ++s)  {
      xml_comp_t x_slice = s;
      string     s_name  = _toString(s_num,"slice%d");
      double     s_thick = x_slice.thickness();
      Material   s_mat   = theDetector.material(x_slice.materialStr());
      Volume     s_vol(s_name,PolyhedraRegular(nsides_outer,rmin,rmax,s_thick),s_mat);
      
      s_vol.setVisAttributes(theDetector.visAttributes(x_slice.visStr()));
      sliceZ += s_thick/2;
      PlacedVolume s_phv = l_vol.placeVolume(s_vol,Position(0,0,sliceZ));
      if ( x_slice.isSensitive() )  {
        sens.setType("calorimeter");
        s_vol.setSensitiveDetector(sens);
        sensitives.push_back(s_phv);
	th_i += s_thick / 2. ;
	th_o  = s_thick / 2. ;
      } else {
	if( th_o < 0. ){
	  th_i += s_thick;
	} else {
	  th_o += s_thick;
	}
      }
      
      if( x_slice.isRadiator() == true)
        totalAbsorberThickness+= s_thick;
      
      sliceZ += s_thick/2;
      s_num++;
    }
    l_vol.setVisAttributes(theDetector.visAttributes(x_layer.visStr()));
    if ( l_repeat <= 0 ) throw std::runtime_error(x_det.nameStr()+"> Invalid repeat value");
    for(int j=0; j<l_repeat; ++j) {
      string phys_lay = _toString(l_num,"layer%d");
      layerZ += l_thick/2;
      DetElement    layer_elt(endcapA, phys_lay, l_num);
      PlacedVolume  pv = endcapVol.placeVolume(l_vol,Position(0,0,layerZ));
      pv.addPhysVolID("layer", l_num);
      layer_elt.setPlacement(pv);
      for(size_t ic=0; ic<sensitives.size(); ++ic)  {
        PlacedVolume sens_pv = sensitives[ic];
        DetElement comp_elt(layer_elt,sens_pv.volume().name(),l_num);
        comp_elt.setPlacement(sens_pv);
      }
      
      ///FIXME: IS ORIENTATION RIGHT? WHICH SIDE DO WE NEED TO ADD TO STRUCTURE?
      LayeredCalorimeterData::Layer caloLayer ;
      caloLayer.distance = zmin + totalThickness/2 + layerZ;
      caloLayer.inner_thickness = th_i ;
      caloLayer.outer_thickness = th_o ;
      caloLayer.absorberThickness = totalAbsorberThickness;
      caloLayer.cellSize0 = cell_sizeX; 
      caloLayer.cellSize1 = cell_sizeY; 
      
      caloData->layers.push_back( caloLayer ) ;
//       std::cout<<"Layer "<<j<<" distance= " <<caloLayer.distance << " layerZ= " << layerZ<<std::endl;
      
      layerZ += l_thick/2;
      ++l_num;
    }
    ++layerType;
  }
  
  double z_pos = zmin+totalThickness/2;
  PlacedVolume pv;
  // Reflect it.

  DetElement  endcapB = endcapA.clone(det_name+"_B",x_det.id());
  
  //Removed rotations to align with envelope
  pv = envelope.placeVolume(endcapVol,Transform3D(RotationZYX(0,0,0),
                                                  Position(0,0,z_pos)));
  pv.addPhysVolID("side", 1);
  endcapA.setPlacement(pv);
  
  //Removed rotations
  pv = envelope.placeVolume(endcapVol,Transform3D(RotationZYX(0,M_PI,0),
                                                  Position(0,0,-z_pos)));
  pv.addPhysVolID("side", 2);
  endcapB.setPlacement(pv);
  

  sdet.add(endcapB);
  
  sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
  
  return sdet;
  
}

DECLARE_DETELEMENT(ECalEndcap_o2_v01,create_detector)


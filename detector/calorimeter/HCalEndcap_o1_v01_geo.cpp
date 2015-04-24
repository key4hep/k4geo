
//====================================================================
//  Modified HCal Endcap Driver for the CLIC detector
//--------------------------------------------------------------------
//
//  Author     : N.Nikiforou
//  Adapted from PolyhedraBarrel Calorimeter by M. Frank
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
  xml_det_t   x_det     = e;
  int         det_id    = x_det.id();
  string      det_name  = x_det.nameStr();
  DetElement    sdet      (det_name,det_id);

    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = XML::createPlacedEnvelope( lcdd,  e , sdet ) ;
    
    if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    
    
//    return sdet; //just temporary
    
  std::cout<<"Building HCal EndCap inside envelope."<<std::endl;
  xml_dim_t   dim       = x_det.dimensions();
  Material    air       = lcdd.air();
  int         numsides  = dim.numsides();
  double      rmin      = dim.rmin();
  double      rmax      = dim.rmax()*std::cos(M_PI/numsides);
  double      zmin      = dim.zmin();
  float       rcutout = dim.hasAttr(_U(rmin2)) ? dim.rmin2() : 0;
  float       zcutout = dim.hasAttr(_U(z2)) ? dim.z2() : 0;
  
  
  Layering    layering(x_det);
  double      totalThickness = layering.totalThickness();
  PolyhedraRegular polyVolume(numsides,rmin,rmax,totalThickness);
  PolyhedraRegular cutoutPolyVolume(numsides,0,rmin+rcutout,zcutout);
  Position cutoutPos(0,0,(zcutout-totalThickness)/2.0);
  std::cout<<"Cutout z width will be  "<<zcutout<<std::endl; 
  Volume      endcapVol("endcap",SubtractionSolid(polyVolume,cutoutPolyVolume,cutoutPos),air);
  DetElement  endcapA(sdet,"endcap",det_id);
  Ref_t(endcapA)->SetName((det_name+"_A").c_str());
  
  int l_num = 1;
  int layerType   = 0;
  double layerZ   = -totalThickness/2;
  
  endcapVol.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  
  for(xml_coll_t c(x_det,_U(layer)); c; ++c)  {
    xml_comp_t       x_layer  = c;
    double           l_thick  = layering.layer(l_num-1)->thickness();
    string           l_name   = _toString(layerType,"layer%d");
    int              l_repeat = x_layer.repeat();
    float            l_rcutout = x_layer.hasAttr(_U(gap)) ? x_layer.gap() : 0;
    
    std::cout<<"Number of layers in group "<<layerType<<" : "<<l_repeat<<std::endl; 
    
    Volume           l_vol(l_name,PolyhedraRegular(numsides,rmin+l_rcutout,rmax,l_thick),air);
    vector<PlacedVolume> sensitives;
    
    int s_num = 1;
    double sliceZ = -l_thick/2;
    for(xml_coll_t s(x_layer,_U(slice)); s; ++s)  {
      xml_comp_t x_slice = s;
      string     s_name  = _toString(s_num,"slice%d");
      double     s_thick = x_slice.thickness();
      Material   s_mat   = lcdd.material(x_slice.materialStr());
      Volume     s_vol(s_name,PolyhedraRegular(numsides,rmin+l_rcutout,rmax,s_thick),s_mat);
      
      s_vol.setVisAttributes(lcdd.visAttributes(x_slice.visStr()));
      sliceZ += s_thick/2;
      PlacedVolume s_phv = l_vol.placeVolume(s_vol,Position(0,0,sliceZ));
      s_phv.addPhysVolID("slice",s_num);
      if ( x_slice.isSensitive() )  {
        sens.setType("calorimeter");
        s_vol.setSensitiveDetector(sens);
        sensitives.push_back(s_phv);
      }
      sliceZ += s_thick/2;
      s_num++;
    }
    l_vol.setVisAttributes(lcdd.visAttributes(x_layer.visStr()));
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
      layerZ += l_thick/2;
      ++l_num;
    }
    ++layerType;
  }
  
  double z_pos = zmin+totalThickness/2;
  PlacedVolume pv;
  // Reflect it.

  DetElement  endcapB = endcapA.clone(det_name+"_B",x_det.id());
  
  pv = envelope.placeVolume(endcapVol,Transform3D(RotationZYX(M_PI/numsides,0,0),
                                                  Position(0,0,z_pos)));
  pv.addPhysVolID("barrel", 1);
  endcapA.setPlacement(pv);
  
  pv = envelope.placeVolume(endcapVol,Transform3D(RotationZYX(M_PI/numsides,M_PI,0),
                                                  Position(0,0,-z_pos)));
  pv.addPhysVolID("barrel", 2);
  endcapB.setPlacement(pv);
  
//   pv.addPhysVolID("system", det_id);
//   both_endcaps.setPlacement(pv);
//   sdet.add(sdetA);
  sdet.add(endcapB);
  return sdet;
  
}

DECLARE_DETELEMENT(HCalEndcap_o1_v01,create_detector)


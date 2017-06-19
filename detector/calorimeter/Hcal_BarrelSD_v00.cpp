//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  Geometry driver for semi-digitia HcalBarrel
//--------------------------------------------------------------------
//  G.Grenier 
//  $Id$
//====================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "LcgeoExceptions.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::PlacedVolume;
using dd4hep::Ref_t;
using dd4hep::SensitiveDetector;

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens) {

  // static double tolerance = 0e0;

  xml_det_t   x_det       = element;
  string      det_name    = x_det.nameStr();

  // xml_comp_t   x_dim      = x_det.dimensions();
  // int          nsides     = x_dim.numsides();

  // Hcal Barrel module shapers' parameters
  // double  Hcal_inner_radius   = x_dim.rmin();
  // double  Hcal_outer_radius   = x_dim.rmax();
  // double  detZ                = x_dim.zhalf();

  DetElement  sdet( det_name,x_det.id() );


  // --- create an envelope volume and position it into the world ---------------------
  
  // Volume envelope =
  dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  PlacedVolume pv;

  sens.setType("calorimeter");


  // Some verbose output
  cout << " \n\n\n CREATE DETECTOR: Hcal_BarrelSD_v00" << endl;


  //-----------------------------------------------------------------------------------
  //  
  //  ...  construct the calorimeter in its envelope
  // 
  //-----------------------------------------------------------------------------------



  //  sdet.addExtension< DDRec::LayeredCalorimeterData >( caloData ) ;
  return sdet;
}
DECLARE_DETELEMENT(Hcal_BarrelSD_v00, create_detector)

//=========================================================================
//  Barrel ECal driver implementation for the CLIC NDM
//-------------------------------------------------------------------------
//  Baseon on
//   * S. Lu's driver for SiWEcalBarrel, ported from Mokka
//   * M. Frank's generic driver EcalBarrel_geo.cpp
//-------------------------------------------------------------------------
//  D.Protopopescu, Glasgow
//  M.Frank, CERN
//  S.Lu, DESY
//  $Id: ECalBarrel_o1_v01_geo.cpp 328 2015-03-18 14:40:50Z /C=UK/O=eScience/OU=Glasgow/L=Compserv/CN=dan protopopescu $
//=========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Layering.h"
#include "TGeoTrd2.h"
#include "XML/Utilities.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Ref_t;
using dd4hep::SensitiveDetector;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector ) {

    // static double tolerance = 0e0;
    
    xml_det_t     x_det     = e;
    int           det_id    = x_det.id();
    string        det_name  = x_det.nameStr();
    DetElement    sdet      (det_name,det_id);
    
    // --- create an envelope volume and position it into the world ---------------------
    
    // Volume envelope =
    dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;
    
    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    
    
    return sdet; //just temporary
    
}

DECLARE_DETELEMENT(ECalPlug_o1_v01,create_detector)

//====================================================================
//  LCGeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  Shaojun Lu, DESY
//  $Id$
//====================================================================

#include "SServices00_v01.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"
//#include "XMLHandlerDB.h"
#include <cmath>
#include <map>
#include <string>

using namespace std;

using dd4hep::Assembly;
using dd4hep::ConeSegment;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::RotateY;
using dd4hep::RotationY;
using dd4hep::Solid;
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Tube;
using dd4hep::Volume;
using dd4hep::Ref_t;
using dd4hep::Detector;
using dd4hep::DetElement;

static Ref_t create_element(Detector& lcdd, xml_h element, Ref_t)  {

  //unused:  static double tolerance = 0e0;

  xml_det_t x_det = element;
  string det_name = x_det.nameStr();
  
  int det_id = x_det.id();
  DetElement sdet (det_name,det_id);
  Volume motherVol = lcdd.pickMotherVolume(sdet);
  
  Assembly envelope_assembly( det_name + "assembly" ) ;
  PlacedVolume pv;

  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;

//====================================================================
// build all services
//====================================================================
  std::cout << "\nBuilding SServices00_v01"<< std::endl;

  SServices00_v01 services;
  services.BuildServices(pv,envelope_assembly,lcdd);

//====================================================================
// Place services into the world volume
//====================================================================

  pv = motherVol.placeVolume(envelope_assembly);
  pv.addPhysVolID("system", det_id);
  sdet.setVisAttributes( lcdd, x_det.visStr(),  envelope_assembly);
  sdet.setPlacement(pv);

  return sdet;
}

DECLARE_DETELEMENT(ILDSServices,create_element)



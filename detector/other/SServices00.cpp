//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
// DD4hep Geometry driver for SService00
// Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:$
//====================================================================
//
//*******************************************************
//*                                                     *
//*                      Mokka                          * 
//*   - the detailed geant4 simulation for ILC   -      *
//*                                                     *
//* For more information about Mokka, visit the         *
//*                                                     *
//*  Mokka.in2p3.fr  Mokka home page.                   *
//*                                                     *
//*******************************************************
//
// $Id: SServices00.cc,v 1.1 2009/06/04 11:40:42 musat Exp $
// $Name: mokka-07-00 $
//
// 
//
// SServices00.cc
//

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "XMLHandlerDB.h"
 
using namespace std;
using namespace DD4hep;
using namespace dd4hep;
using namespace DD4hep::Geometry;

bool BuildTPCEndplateServices(PlacedVolume *pVol,Assembly *envelope)
{

  return true;
}

bool BuildEcalBarrelServices(PlacedVolume *pVol,Assembly *envelope)
{

  return true;
}

bool BuildEcalBarrel_EndCapServices(PlacedVolume *pVol,Assembly *envelope)
{

  return true;
}

bool BuildHcalBarrel_EndCapServices(PlacedVolume *pVol,Assembly *envelope)
{

  return true;
}

bool BuildSitCables(PlacedVolume *pVol,Assembly *envelope)
{

  return true;
}


static Ref_t create_element(LCDD& lcdd, xml_h element, SensitiveDetector sens)  {

  static double tolerance = 0e0;

  xml_det_t x_det = element;
  string det_name = x_det.nameStr();
  
  int det_id = x_det.id();
  DetElement sdet (det_name,det_id);
  Volume motherVol = lcdd.pickMotherVolume(sdet);
  
  Assembly envelope_assembly( det_name + "assembly" ) ;
  PlacedVolume pv;


//====================================================================
// build all services
//====================================================================
  cout << "\nBuilding SServices00"<< endl;

  //BuildTPCEndplateServices(&pv,&envelope_assembly);

  //BuildEcalBarrelServices(&pv,&envelope_assembly);

  //BuildEcalBarrel_EndCapServices(&pv,&envelope_assembly);

  //BuildHcalBarrel_EndCapServices(&pv,&envelope_assembly);

  //BuildSitCables(&pv,&envelope_assembly);


//====================================================================
// Place services into the world volume
//====================================================================

  pv = motherVol.placeVolume(envelope_assembly);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  return sdet;
}

DECLARE_DETELEMENT(SServices00,create_element)

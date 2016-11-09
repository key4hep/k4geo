// $Id: $
//==========================================================================
//  AIDA Detector description implementation for LCD
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
// Modified by: D.Protopopescu 
// Includes extensions as implemented by N. Nikiforou in TrackerEndcap_o2_v04_geo.cpp
//
//==========================================================================
//
// Specialized generic detector constructor
// 
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include <map>
#include "DDRec/DetectorData.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
  typedef vector<PlacedVolume> Placements;
  xml_det_t   x_det     = e;
  Material    vacuum    = lcdd.vacuum();
  int         det_id    = x_det.id();
  string      det_name  = x_det.nameStr();
  bool        reflect   = x_det.reflect(false);
  DetElement  sdet        (det_name,det_id);
  int         m_id=0, c_id=0, n_sensor=0;
  map<string,Volume> modules;
  map<string, Placements>  sensitives;
  PlacedVolume pv;

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = XML::createPlacedEnvelope( lcdd,  e , sdet ) ;
  XML::setDetectorTypeFlag( e, sdet ) ;
  
  if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;
  
  //-----------------------------------------------------------------------------------
  
  DDRec::ZDiskPetalsData*  zDiskPetalsData = new DDRec::ZDiskPetalsData ;
  
  envelope.setVisAttributes(lcdd.invisible());
  sens.setType("tracker");

  for(xml_coll_t mi(x_det,_U(module)); mi; ++mi, ++m_id)  {
    xml_comp_t x_mod   = mi;
    string     m_nam   = x_mod.nameStr();
    xml_comp_t trd     = x_mod.trd();
    double     posY;
    double     x1      = trd.x1();
    double     x2      = trd.x2();
    double     z       = trd.z();
    double     y1, y2, total_thickness=0.;
    xml_coll_t ci(x_mod,_U(module_component));
    for(ci.reset(), total_thickness=0.0; ci; ++ci)
      total_thickness += xml_comp_t(ci).thickness();
      
    y1 = y2 = total_thickness / 2;
    Volume  m_volume(m_nam, Trapezoid(x1, x2, y1, y2, z), vacuum);      
    m_volume.setVisAttributes(lcdd.visAttributes(x_mod.visStr()));

    // Loop over slices 
    // The first slice (top in the xml) is placed at the "bottom" of the module
    for(ci.reset(), n_sensor=1, c_id=0, posY=-y1; ci; ++ci, ++c_id)  {
      xml_comp_t c       = ci;
      double     c_thick = c.thickness();
      Material   c_mat   = lcdd.material(c.materialStr());
      string     c_name  = _toString(c_id,"component%d");
      Volume     c_vol(c_name, Trapezoid(x1,x2,c_thick/2e0,c_thick/2e0,z), c_mat);

      c_vol.setVisAttributes(lcdd.visAttributes(c.visStr()));
      pv = m_volume.placeVolume(c_vol,Position(0,posY+c_thick/2,0));
      if ( c.isSensitive() ) {
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        ++n_sensor;
      }
      posY += c_thick;
    }
    modules[m_nam] = m_volume;
  }
  
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    xml_comp_t  x_layer(li);
    int l_id    = x_layer.id();
    int mod_num = 0;
    int ring_num = 0;

    double sumZ(0.), innerR(1e100), outerR(0.);

    for(xml_coll_t ri(x_layer,_U(ring)); ri; ++ri)  {
      xml_comp_t x_ring = ri;
      double r        = x_ring.r();
      double phi0     = x_ring.phi0(0);
      double zstart   = x_ring.zstart();
      double dz       = x_ring.dz(0);
      int    nmodules = x_ring.nmodules();
      string m_nam    = x_ring.moduleStr();
      Volume m_vol    = modules[m_nam];
      double iphi     = 2*M_PI/nmodules;
      double phi      = phi0;
      Placements& sensVols = sensitives[m_nam];

      DD4hep::Geometry::Box mod_shape(m_vol.solid());
      if(r - mod_shape->GetDZ()<innerR)
	innerR = r - mod_shape->GetDZ();
      if(r + mod_shape->GetDZ()>outerR)
	outerR = r + mod_shape->GetDZ();
      sumZ += zstart;
      r = r + mod_shape->GetDY();

      for(int k=0; k<nmodules; ++k) {
        string m_base = _toString(l_id,"layer%d") + _toString(mod_num,"_module%d") + _toString(k,"_sensor%d");
        double x = -r*std::cos(phi);
        double y = -r*std::sin(phi);
        DetElement module(sdet,m_base+"_pos",det_id);
        pv = envelope.placeVolume(m_vol,Transform3D(RotationZYX(0,-M_PI/2-phi,-M_PI/2),Position(x,y,zstart+dz)));
        pv.addPhysVolID("side",1).addPhysVolID("layer", l_id).addPhysVolID("module",mod_num).addPhysVolID("sensor",k);
        module.setPlacement(pv);
        for(size_t ic=0; ic<sensVols.size(); ++ic)  {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement comp_elt(module,sens_pv.volume().name(),mod_num);
          comp_elt.setPlacement(sens_pv);
        }

        if ( reflect ) {
          pv = envelope.placeVolume(m_vol,Transform3D(RotationZYX(M_PI,-M_PI/2-phi,-M_PI/2),Position(x,y,-zstart-dz)));
          pv.addPhysVolID("side",-1).addPhysVolID("layer",l_id).addPhysVolID("module",mod_num).addPhysVolID("sensor",k);
          DetElement r_module(sdet,m_base+"_neg",det_id);
          r_module.setPlacement(pv);
          for(size_t ic=0; ic<sensVols.size(); ++ic)  {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement comp_elt(r_module,sens_pv.volume().name(),mod_num);
            comp_elt.setPlacement(sens_pv);
          }
        }
        dz   = -dz;
        phi += iphi;
      }
      ++mod_num;
      ++ring_num;
    }

    DDRec::ZDiskPetalsData::LayerLayout thisLayer ;

    // Only filling what is needed for CED/DDMarlinPandora
    thisLayer.zPosition = sumZ/ring_num; // average z
    thisLayer.distanceSensitive = innerR;
    thisLayer.lengthSensitive = outerR - innerR;
    thisLayer.petalNumber = ring_num; // number of rings in petalNumber, needed for tracking
    zDiskPetalsData->layers.push_back(thisLayer);

  }

  sdet.addExtension< DDRec::ZDiskPetalsData >( zDiskPetalsData ) ;

  return sdet;
}

DECLARE_DETELEMENT(SiTrackerEndcap_o2_v01ext,create_detector)

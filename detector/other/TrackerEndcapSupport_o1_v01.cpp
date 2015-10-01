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
//
//==========================================================================
//
// Specialized generic detector constructor
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t     x_det     = e;
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  DetElement    sdet(det_name, det_id);

  Volume envelope = XML::createPlacedEnvelope(lcdd,  e , sdet) ;

  if (lcdd.buildType() == BUILD_ENVELOPE) return sdet ;

  Material   air       = lcdd.air();
  bool       reflect   = x_det.reflect();

  PlacedVolume pv;
  int l_num = 0;

  for (xml_coll_t i(x_det, _U(layer)); i; ++i, ++l_num)  {
    xml_comp_t x_layer = i;
    string l_nam = det_name + _toString(l_num, "_layer%d");
    double  zmin = x_layer.inner_z();
    double  rmin = x_layer.inner_r();
    double  rmax = x_layer.outer_r();
    double  z    = zmin, layerWidth = 0.;
    int     s_num = 0;

    for (xml_coll_t j(x_layer, _U(slice)); j; ++j)  {
      double thickness = xml_comp_t(j).thickness();
      layerWidth += thickness;
    }
    Tube    l_tub(rmin, rmax, layerWidth, 2 * M_PI);
    Volume  l_vol(l_nam, l_tub, air);
    l_vol.setVisAttributes(lcdd, x_layer.visStr());
    for (xml_coll_t j(x_layer, _U(slice)); j; ++j, ++s_num)  {
      xml_comp_t x_slice = j;
      double thick = x_slice.thickness();
      Material mat = lcdd.material(x_slice.materialStr());
      string s_nam = l_nam + _toString(s_num, "_slice%d");
      Volume s_vol(s_nam, Tube(rmin, rmax, thick), mat);

      s_vol.setAttributes(lcdd, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
      pv = l_vol.placeVolume(s_vol, Position(0, 0, z - zmin - layerWidth / 2 + thick / 2));
      pv.addPhysVolID("sensor", s_num);
    }

    DetElement layer(sdet, l_nam + "_pos", l_num);
    pv = envelope.placeVolume(l_vol, Position(0, 0, zmin + layerWidth / 2.));
    pv.addPhysVolID("layer", l_num);
    pv.addPhysVolID("barrel", 1);
    layer.setPlacement(pv);
    if (reflect)  {
      pv = envelope.placeVolume(l_vol, Transform3D(RotationY(M_PI), Position(0, 0, -zmin - layerWidth / 2)));
      pv.addPhysVolID("layer", l_num);
      pv.addPhysVolID("barrel", 2);
      DetElement layerR = layer.clone(l_nam + "_neg");
      sdet.add(layerR.setPlacement(pv));
    }
  }


  return sdet;
}

DECLARE_DETELEMENT(TrackerEndcapSupport_o1_v01, create_detector)

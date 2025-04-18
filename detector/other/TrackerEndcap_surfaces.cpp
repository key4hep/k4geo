// $Id$
//==========================================================================
//  Surface installer plugin for TrackerEndcap geometry driver
//--------------------------------------------------------------------------
//
// Author     : N. Nikiforou, forked from DD4hep/SiTrackerEndcap_surfaces.cpp
//              by M. Frank
//==========================================================================
//
// Specialized generic detector constructor
//
//==========================================================================
// Framework include files
#define DD4HEP_USE_SURFACEINSTALL_HELPER TrackerEndcapSurfacePlugin
#include "DD4hep/SurfaceInstaller.h"

/// Install measurement surfaces
template <typename UserData>
void Installer<UserData>::install(dd4hep::DetElement component, dd4hep::PlacedVolume pv) {
  dd4hep::Volume comp_vol = pv.volume();
  if (comp_vol.isSensitive()) {
    dd4hep::Volume mod_vol = parentVolume(component);
    dd4hep::Trapezoid comp_shape(comp_vol.solid()), mod_shape(mod_vol.solid());

    if (!comp_shape.isValid() || !mod_shape.isValid()) {
      invalidInstaller("Components and/or modules are not Trapezoid -- invalid shapes");
    } else if (!handleUsingCache(component, comp_vol)) {
      const double* trans = placementTranslation(component);
      double half_mod_thickness = mod_shape->GetDy1();
      // double half_comp_thickness = comp_shape->GetDy1();
      double si_position = trans[1];
      double outer_thickness = half_mod_thickness - si_position;
      double inner_thickness = half_mod_thickness + si_position;
      // After the rotation of the module and placement to the endcap perpendicular to z
      // the orientation should be n parallel to z, u along r and v along rPhi
      Vector3D u(0., 0., 1.), v(1., 0., 0.), n(0., 1., 0.), o(0., 0., 0.);

      VolPlane surf(comp_vol, Type(Type::Sensitive, Type::Measurement1D), inner_thickness, outer_thickness, u, v, n, o);
      addSurface(component, surf);
    }
  }
}

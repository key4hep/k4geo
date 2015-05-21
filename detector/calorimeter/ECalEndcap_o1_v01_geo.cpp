//=========================================================================
//  EndCap ECal driver implementation for the CLIC NDM
//-------------------------------------------------------------------------
//  Based on
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
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
    
    xml_det_t     x_det     = e;
    int           det_id    = x_det.id();
    string        det_name  = x_det.nameStr();
    DetElement    sdet      (det_name,det_id);
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = XML::createPlacedEnvelope(lcdd, e, sdet);
    
    if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
	return sdet; //Temporary fox for the duplicated volumes -> Please fix it Dan.    
    DD4hep::XML::Component xmlParameter = x_det.child(_Unicode(parameter));
    const double tolerance = xmlParameter.attr<double>(_Unicode(ecal_endcap_tolerance));
    const double faceThickness = xmlParameter.attr<double>(_Unicode(FaceThickness));
    Material WDens24 = lcdd.material(xmlParameter.attr<string>(_Unicode(RadiatorMaterial)));
    
    Layering      layering (e);
    Material      air       = lcdd.air();
    xml_comp_t    x_sectors = x_det.staves();
    xml_comp_t    x_dim     = x_det.dimensions();
    int           nsides    = x_dim.numsides();
    double        dphi      = (2.*M_PI/nsides);
    double        hphi      = dphi/2.;
    double        inner_r   = x_dim.rmin()*std::sqrt(2);   // inscribing circle radius
    double        outer_r   = x_dim.rmax()/std::cos(hphi); // circumscribing circle radius
    double        z_min     = x_dim.zmin();
    double        z_max     = x_dim.zmax();
    double        mod_z     = layering.totalThickness();

    // The endcap has a square hole for the ECal plug, but the 12/4 symmetry is awkward
    // Will build 12 identical sectors, and each made of trapezoidal layers
    
    std::cout << "Sector z-thickness: " << mod_z << std::endl;
    std::cout << "Sector r_out=" << outer_r << " and r_in=" << inner_r << std::endl;
    std::cout << "Sector z_min=" << z_min << " and z_max=" << z_max << std::endl;
   
    if((z_max - z_min) < mod_z) { // throw exception
      throw std::runtime_error("ERROR: Layers don't fit within the envelope volume!");
    }

    DetElement sector_det("sector0", det_id);

    // Compute sector dimensions
    double trd_x1 = mod_z/2.;
    double trd_x2 = mod_z/2.;
    double trd_y1 = inner_r*std::sin(hphi) - tolerance;
    double trd_y2 = outer_r*std::sin(hphi) - tolerance;
    double trd_z  = (outer_r - inner_r)*std::cos(hphi)/2 - tolerance;

    // Create the trapezoid for the sector
    Trapezoid sec(trd_x1, 
                  trd_x2, 
                  trd_y1, 
                  trd_y2, 
                  trd_z); 
    Volume sector_vol("sector", sec, WDens24); // solid tungsten volume

    sens.setType("calorimeter");

    {// Build the sector, layer by layer and slice by slice
      
      // Parameters for computing the layer X dimension:
      double l_y1 = trd_y1 - faceThickness/std::sin(hphi);  
      double l_y2 = trd_y2 - faceThickness/std::sin(hphi);
      double l_zh = trd_z;
      double l_pos_z  = -trd_x1 + tolerance;
      
      // Loop over the sets of layer elements in the detector
      int l_num = 1;
      for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
	
	xml_comp_t x_layer = li;
	int repeat = x_layer.repeat();
	// Loop over number of repeats for this layer.
	for (int j=0; j<repeat; j++)  {
	  
	  string l_name = _toString(l_num, "layer%d");
	  double l_thickness = layering.layer(j)->thickness();  // this layer's thickness
	  
	  Position   l_pos(l_pos_z + l_thickness/2., 0. , 0.);   // position of the layer
	  Trapezoid  l_box(l_thickness/2., l_thickness/2., l_y1 - tolerance, l_y2 - tolerance, l_zh - tolerance);
	  Volume     l_vol(l_name, l_box, air);
	  DetElement layer(sector_det, l_name, det_id);
	  
	  // Loop over the slices for this layer
	  int s_num = 1;
	  double s_pos_z = -(l_thickness/2);
	  for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	    
	    xml_comp_t x_slice = si;
	    string     s_name  = _toString(s_num, "slice%d");
	    double     s_thick = x_slice.thickness();
	    //std::cout << s_name << " thickness: " << s_thick << std::endl;
            Trapezoid  s_box(s_thick/2., s_thick/2., l_y1 - tolerance, l_y2 - tolerance, l_zh - tolerance);
	    Volume     s_vol(s_name, s_box, lcdd.material(x_slice.materialStr()));
	    DetElement slice(layer, s_name, det_id);
	    
	    if ( x_slice.isSensitive() ) {
	      s_vol.setSensitiveDetector(sens);
	    }
	    slice.setAttributes(lcdd, s_vol, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	    
	    // Slice placement
	    PlacedVolume slice_phv = l_vol.placeVolume(s_vol, Position(s_pos_z + s_thick/2, 0, 0));
	    slice_phv.addPhysVolID("slice", s_num);
	    slice.setPlacement(slice_phv);
	    // Increment Z position of slice
	    s_pos_z += s_thick;
	    
	    // Increment slice number
	    ++s_num;
	  }
	  
	  // Set region, limitset, and visibility of layer
	  layer.setAttributes(lcdd, l_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
	  
	  PlacedVolume layer_phv = sector_vol.placeVolume(l_vol, l_pos);
	  layer_phv.addPhysVolID("layer", l_num);
	  layer.setPlacement(layer_phv);
	  // Increment to next layer Z position
	  l_pos_z += l_thickness;
	  ++l_num;
	}
      }
    }

    // Set visualization
    if ( x_sectors )   {
        sector_vol.setVisAttributes(lcdd.visAttributes(x_sectors.visStr()));
    }

    // Phi start for sectors
    double phi = hphi;
    double r_off = (inner_r + outer_r)*std::cos(hphi)/2;    
    double z_off = (z_max + z_min)/2;

    for (int j = -1; j < 2; j += 2){ 
      for (int i = 0; i < nsides; i++, phi -= dphi) { 
	// Compute the sector position in envelope coordinates
	double m_pos_x = r_off * std::sin(phi);
	double m_pos_y = r_off * std::cos(phi);
	double m_pos_z = j*z_off;
        Transform3D tr(RotationZYX(M_PI/2, (j+1)*M_PI/2 - j*phi, j*M_PI/2.), Translation3D(m_pos_x, m_pos_y, m_pos_z));
	PlacedVolume pv = envelope.placeVolume(sector_vol, tr);
	int mid = (nsides+1)*(2+j)+i;
	pv.addPhysVolID("barrel", j); 
	pv.addPhysVolID("sector", i);
	DetElement sd = mid==0 ? sector_det : sector_det.clone(_toString(mid,"sector%d"));
	sd.setPlacement(pv);
	sdet.add(sd);
      }
    }

    // Set envelope volume attributes
    envelope.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
    return sdet;
}

DECLARE_DETELEMENT(ECalEndcap_o1_v01,create_detector)

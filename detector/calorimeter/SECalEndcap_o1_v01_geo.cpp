//=========================================================================
//  EndCap SiWECal driver implementation for (the CLIC NDM) and ILD
//-------------------------------------------------------------------------
//  Based on
//   * D.Protopopescu's ECalEndcap_o1_v01_geo.cpp
//   * https://www.evernote.com/l/AJ0hig1OWixLH7McEetl6XiioRlf4ycRPzw
//-------------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:  $
//=========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/DetType.h"
#include "XML/Layering.h"
#include "TGeoTrd2.h"
#include "TMath.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::DetType;
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Trap;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;


static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens)  {
    
    xml_det_t     x_det     = e;
    int           det_id    = x_det.id();
    string        det_name  = x_det.nameStr();
    DetElement    sdet      (det_name,det_id);
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, e, sdet);
    
    sdet.setTypeFlag( DetType::CALORIMETER |  DetType::ENDCAP  | DetType::ELECTROMAGNETIC ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    dd4hep::xml::Component xmlParameter = x_det.child(_Unicode(parameter));
    const double tolerance = xmlParameter.attr<double>(_Unicode(ecal_endcap_tolerance));
    const int insides = xmlParameter.attr<double>(_Unicode(InSides));
    Material Composite = theDetector.material(xmlParameter.attr<string>(_Unicode(CarbonFiber)));

    Layering      layering (e);
    Material      air       = theDetector.air();
    xml_comp_t    x_panels  = x_det.staves();
    xml_comp_t    x_dim     = x_det.dimensions();
    int           outsides  = x_dim.numsides();
    double        dphi      = (2.*M_PI/outsides);
    double        hphi      = dphi/2.;
    double        inner_r   = x_dim.rmin();  // inscribed circle radius
    double        outer_r   = x_dim.rmax();  // inscribed circle radius
    double        z_min     = x_dim.zmin();
    double        z_max     = x_dim.zmax();
    double        mod_z     = layering.totalThickness() + 2*tolerance;
    double eCal_cell_size = theDetector.constant<double>("ECal_cell_size"); //Should try to obtain from segmentation
    
    
    
    //Create caloData object to extend driver with data required for reconstruction
    LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
    caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
    caloData->inner_symmetry = insides;
    caloData->outer_symmetry = outsides; 
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
    caloData->extent[0] = inner_r;
    caloData->extent[1] = outer_r ; ///FIXME: CHECK WHAT IS NEEDED (EXSCRIBED?)
    caloData->extent[2] = z_min ;
    caloData->extent[3] = z_min + mod_z;
    
    
    
    if(insides != 4){ // throw exception
      throw std::runtime_error("ERROR: This driver is for 4 inner symmetry!");
    }
    if(outsides%4 != 0){ // throw exception
      throw std::runtime_error("ERROR: This driver is for 4*N (N=2,3) outer symmetry!");
    }

    sens.setType("calorimeter");
    
    std::cout << "Global tolerance: " << tolerance << " cm" << std::endl;
    std::cout << "Panels z-thickness: " << mod_z << " cm" << std::endl;
    std::cout << "EndCap z_min = " << z_min << " cm and z_max = " << z_max << " cm" << std::endl;

    if((z_max - z_min) < mod_z) { // throw exception
      throw std::runtime_error("ERROR: Layers don't fit within the envelope volume!");
    } else std::cout << "Got z_max - z_min = " << z_max - z_min << " cm" << std::endl;

    // The actual panel shapes depend on whether the side of the inner cutout (l_i)
    // is bigger or smaller than the side of the outer polygone l_o
    double l_i = 2.*inner_r; 
    double l_o = 2.*outer_r*std::tan(hphi);
    std::cout << "Got l_i = " << l_i << " cm, and l_o = " << l_o << " cm" << std::endl; 

    if(outsides == 12 && l_i < l_o) { // throw exception
      throw std::runtime_error("ERROR: This driver is for the l_i > l_o case when outer symmetry is 12!");
    } else {
      std::cout << "EndCap r_out = " << outer_r << " cm (inscribed circle radius)" << std::endl;
      std::cout << "EndCap r_in = " << inner_r << " cm (inscribed circle radius)" << std::endl;
    }
    
    // This driver is for a 12/4 symmetry, with 4 sectors containing
    // 3 right trapezoidal modules each (l_i > l_o)
    // Dimensions based on the sketch from http://cern.ch/go/R89p
    double ta = abs(l_i - l_o)/2;        
    double tc = ((outer_r - inner_r)>(l_o*std::sin(dphi)))? (l_o*std::sin(dphi)):(outer_r - inner_r);
    double tb = (outer_r/std::cos(hphi)*std::cos(hphi+2*dphi) > inner_r)?
      (outer_r - outer_r/std::cos(hphi)*std::cos(hphi+2*dphi) - tc):(outer_r - inner_r - tc);
    double td = l_o  + tc/std::tan(dphi);
    double te = td + tb/std::tan(2*dphi);
    double tf = tc + tb;
    double tg = tf - ta/std::tan(2*dphi);

    std::cout << "Set a=" << ta << " b=" << tb << " c=" << tc << " d=" << 
      td << " e=" << te << " f=" << tf << " g=" << tg << " (cm)" << std::endl;

    // Loop over panel modules
    for(int m=1; m<3; m++){

      DetElement panel_det("panel100", det_id);

      // Function parameter names are preserved here for readability
      // Trap(z, theta, phi, y1, x1, x2, alpha1, y2, x3, x4, alpha2)
      // but for this y2 = y1, x3 = x1, x4 = x2, and alpha2 = alpha1
      double z = mod_z/2;
      double theta = 0.0;
      double phi = M_PI/2.;
      double y1, x1, x2, x_off, y_off;
      double tol = 0.9999; // set at 0.96 for visulalising the panels, else ~1
      switch(m){
      case 1:
	y1 = tol*tc/2; // all half-sizes
	x1 = (l_i > l_o)? (tol*td/2):(tol*td/2 - ta/2);
	x2 = (l_i > l_o)? (tol*l_o/2.):(tol*(l_o/2. - ta/2));
	// Offsets in definition coordinates
        x_off = (l_i > l_o)? ((td/2 - l_o/2)/2):((td/2 - l_o/2)/2 + ta/2);
        y_off = (outer_r/std::cos(hphi)*std::cos(hphi+2*dphi) > inner_r)? 
	  (l_i/2 + tb + tc/2 + outer_r/std::cos(hphi)*std::cos(hphi+2*dphi) - inner_r):(l_i/2 + tb + tc/2);
	break;
      case 2:
	y1 = tol*tb/2;
	x1 = (l_i > l_o)? (tol*te/2):(tol*te/2 - ta/2);
	x2 = (l_i > l_o)? (tol*td/2):(tol*td/2 - ta/2);
	x_off = (l_i > l_o)? ((te + td)/4 - l_o/2):((te + td)/4 - l_o/2 + ta/2);
        y_off = (outer_r/std::cos(hphi)*std::cos(hphi+2*dphi) > inner_r)? 
	  (l_i/2 + tb/2 + outer_r/std::cos(hphi)*std::cos(hphi+2*dphi) - inner_r):(l_i/2 + tb/2);
	break;
      case 3:
      default:
	std::cout << "With the current dimensions module 3 is too small to build" << std::endl;
	exit(-3);
      }   
      double alpha1 = -TMath::ATan2(x1 - x2, 2*y1); // in radians
      Trap trm(z, theta, phi, y1, x1, x2, alpha1, y1, x1, x2, alpha1);
      Volume panel_vol(_toString(100*m, "panel%d"), trm, Composite);

      // Build the panels, layer by layer and slice by slice
      double l_pos_z  = -z + tolerance;
	
      // Loop over the sets of layer elements in the detector
      double check_thick = 0.0;
      int l_num = 1;
      int l_set = 0;
      for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
	
	xml_comp_t x_layer = li;
	int repeat = x_layer.repeat();

	// Loop over number of repeats for this layer.
	for (int k=0; k<repeat; k++)  {
    LayeredCalorimeterData::Layer caloLayer ;
    
	  string l_name = _toString(100*m + l_num, "layer%d");
	  double l_thickness = layering.layer(l_num-1)->thickness();  // this layer's thickness	  
	  Position   l_pos(0, 0, l_pos_z + l_thickness/2.);           // position of the layer
	  Trap       l_box(l_thickness/2., theta, phi, y1 - tolerance, x1 - tolerance, x2 - tolerance, alpha1, 
			                               y1 - tolerance, x1 - tolerance, x2 - tolerance, alpha1);
	  Volume     l_vol(l_name, l_box, air);
	  DetElement layer(panel_det, l_name, det_id);
	  
	  // Loop over the slices for this layer
	  int s_num = 1;
	  double s_pos_z = -(l_thickness/2);
	  double totalAbsorberThickness=0.;
    
	  double th_i(0.), th_o(-1.) ;
	  for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	    
	    xml_comp_t x_slice = si;
	    string     s_name  = _toString(s_num, "slice%d");
	    double     s_thick = x_slice.thickness();     // slice thickness
	    Position   s_pos(0, 0, s_pos_z + s_thick/2);  // slice position
	    double     xtol = tolerance + 0.00001;
	    Trap       s_box(s_thick/2., theta, phi, y1 - xtol, x1 - xtol, x2 - xtol, alpha1, 
                                                     y1 - xtol, x1 - xtol, x2 - xtol, alpha1);
	    Volume     s_vol(s_name, s_box, theDetector.material(x_slice.materialStr()));
	    DetElement slice(layer, s_name, det_id);
	    
	    if ( x_slice.isSensitive() ) {
	      s_vol.setSensitiveDetector(sens);
	      th_i += s_thick / 2. ;
	      th_o  = s_thick / 2. ;
	    } else {
	      if( th_o < 0. ){
		th_i += s_thick;
	      } else {
		th_o += s_thick;
	      }
	    }
	    slice.setAttributes(theDetector, s_vol, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	    
      char val = x_slice.hasAttr(_U(radiator)) ? x_slice.attr < string > (_U(radiator))[0] : 'f';
      val = std::toupper(val);
      bool isAbsorber =  (val == 'T' || val == 'Y');
      
      if( isAbsorber ==true)
        totalAbsorberThickness+= s_thick;
      
      
	    // Slice placement
	    PlacedVolume slice_phv = l_vol.placeVolume(s_vol, Position(0, 0, s_pos_z + s_thick/2));
	    slice_phv.addPhysVolID("slice", s_num);
	    slice.setPlacement(slice_phv);

	    // Check
	    check_thick += s_thick;
	    /* std::cout << "Panel " << m << ", layer " <<  l_set 
                      << ", repeat " << k << " (layer" << 100*m + l_num
                      << "), slice " << s_num 
		      << " (dz=" << s_thick << "cm, z=" << s_pos_z + s_thick/2 << "cm)" << std::endl;
	    */

	    // Increment Z position of slice
	    s_pos_z += s_thick;	    
	    // Increment slice number
	    ++s_num;
	  }
	  
	  // Set region, limitset, and visibility of layer
	  layer.setAttributes(theDetector, l_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

	  // std::cout << "Layer " << 100*m + l_num << std::endl;
	  PlacedVolume layer_phv = panel_vol.placeVolume(l_vol, l_pos);
	  layer_phv.addPhysVolID("layer", 100*m + l_num);
	  layer.setPlacement(layer_phv);
    
	  caloLayer.distance = l_pos_z + l_thickness/2.;
	  caloLayer.inner_thickness = th_i ;
	  caloLayer.outer_thickness = th_o ;
	  caloLayer.absorberThickness = totalAbsorberThickness;
	  caloLayer.cellSize0 = eCal_cell_size;
	  caloLayer.cellSize1 = eCal_cell_size;

	  // Increment to next layer Z position
	  l_pos_z += l_thickness;
	  // Increment layer number
	  ++l_num;
	}
	l_set++;
      }

      std::cout << "Total slice thickness check " << check_thick << " cm" << std::endl;
     
      // Set visualization
      if ( x_panels ) {
      	panel_vol.setVisAttributes(theDetector.visAttributes(x_panels.visStr()));
      }

      // Placement
      double zvis = 1.0; // for visualisation set to 1.1, else 1
      double z_off = (z_max + zvis*z_min)/2.;
      double r_off = std::sqrt(x_off*x_off + y_off*y_off);
      double gamma = std::atan(x_off/y_off);

      // Loop over endcaps
      for (int j = -1; j < 2; j += 2){

	// Loop over sectors
	for (int i = 0; i < insides; i++) { 

	  // Compute the panel position in envelope coordinates
	  double m_pos_x = j*r_off*std::sin(gamma);
	  double m_pos_y = r_off*std::cos(gamma);
	  double m_pos_z = j*z_off;
	
	  Transform3D tr(RotationZYX(0, (j-1)*M_PI/2, 0.), Translation3D(m_pos_x, m_pos_y, m_pos_z));
	  PlacedVolume pv = envelope.placeVolume(panel_vol, RotationZYX(i*M_PI/2. + hphi, 0., 0.)*tr);
	  
	  int mid = 100*(3+j)/2 + 10*i + m; // unique composed panel ID
	  // std::cout << "Unique module ID = " << mid << std::endl;
	  pv.addPhysVolID("barrel", j); 
	  pv.addPhysVolID("sector", i);
	  pv.addPhysVolID("panel", m);
	  
	  DetElement sd = mid==100 ? panel_det : panel_det.clone(_toString(mid,"panel%d"));
	  sd.setPlacement(pv);
	  sdet.add(sd);
	  
	}
      }
    }

    // Set envelope volume attributes
    envelope.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
    
    sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
    
    return sdet;
}

DECLARE_DETELEMENT(SECalEndcap_o1_v01,create_detector)

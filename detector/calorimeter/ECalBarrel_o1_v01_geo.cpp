//=========================================================================
//  Barrel ECal driver implementation for the CLIC NDM
//-------------------------------------------------------------------------
//  Based on
//   * S. Lu's driver for SiWEcalBarrel, ported from Mokka
//   * M. Frank's generic driver EcalBarrel_geo.cpp
//-------------------------------------------------------------------------
//  D.Protopopescu, Glasgow
//  M.Frank, CERN
//  S.Lu, DESY
//  $Id$
//=========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Layering.h"
#include "TGeoTrd2.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Trapezoid;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens)  {
    
    xml_det_t     x_det     = e;
    int           det_id    = x_det.id();
    string        det_name  = x_det.nameStr();
    DetElement    sdet      (det_name, det_id);
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, e, sdet);
    
    if(theDetector.buildType() == BUILD_ENVELOPE) return sdet;
    
    //-----------------------------------------------------------------------------------

    dd4hep::xml::Component xmlParameter = x_det.child(_Unicode(parameter));
    const double tolerance = xmlParameter.attr<double>(_Unicode(ecal_barrel_tolerance));
    const int n_towers = (int) xmlParameter.attr<double>(_Unicode(num_towers));
    const int n_rails = (int) xmlParameter.attr<double>(_Unicode(num_rails));
    const int n_stacks = (int) xmlParameter.attr<double>(_Unicode(stacks_per_tower));
    const double towersAirGap = xmlParameter.attr<double>(_Unicode(TowersAirGap));
    const double faceThickness = xmlParameter.attr<double>(_Unicode(TowersFaceThickness));
    const double supportRailCS = xmlParameter.attr<double>(_Unicode(supportRailCrossSection));
    Material Composite = theDetector.material(xmlParameter.attr<string>(_Unicode(AlveolusMaterial)));
    Material Steel = theDetector.material(xmlParameter.attr<string>(_Unicode(supportRailMaterial)));

    Readout readout = sens.readout();
    Segmentation seg = readout.segmentation();
    
    std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
    double cell_sizeX      = cellSizeVector[0];
    double cell_sizeY      = cellSizeVector[1];
    
    Layering      layering (e);
    Material      air       = theDetector.air();
    xml_comp_t    x_staves  = x_det.staves();
    xml_comp_t    x_dim     = x_det.dimensions();
    int           nsides    = x_dim.numsides();
    double        inner_r   = x_dim.rmin();                     // inscribed cylinder
    double        dphi      = (2.*M_PI/nsides);
    double        hphi      = dphi/2.;
    double        mod_z     = layering.totalThickness();
    double        r_max     = x_dim.rmax()/std::cos(hphi);      // circumscribed cylinder, see http://cern.ch/go/r9mZ
    double        ry        = sqrt(supportRailCS)/2;            // support rail thickness
    double        outer_r   = (inner_r + mod_z)/std::cos(hphi); // r_outer of actual module
    
    // Create caloData object to extend driver with data required for reconstruction
    // Added by Nikiforos, 28/05/15. Should not interfere with driver but it does
    LayeredCalorimeterData* caloData = new LayeredCalorimeterData;
    caloData->layoutType = LayeredCalorimeterData::BarrelLayout;
    caloData->inner_symmetry = nsides;
    caloData->outer_symmetry = nsides; 
    
    /** NOTE: phi0=0 means lower face flat parallel to experimental floor
     *  This is achieved by rotating the modules with respect to the envelope
     *  which is assumed to be a Polyhedron and has its axes rotated with respect
     *  to the world by 180/nsides. In any other case (e.g. if you want to have
     *  a tip of the calorimeter touching the ground) this value needs to be computed
     */
    
    caloData->inner_phi0 = 0.; 
    caloData->outer_phi0 = 0.; 
    // One should ensure that these sensitivity gaps are correctly used
    caloData->gap0 = towersAirGap;  // the 4 gaps between the 5 towers, along z 
    caloData->gap1 = faceThickness; // gaps between stacks in a module, along z
    caloData->gap2 = tolerance*std::cos(hphi); // gaps where the staves overlap

    // extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
    caloData->extent[0] = inner_r;
    caloData->extent[1] = outer_r; // or r_max ?
    caloData->extent[2] = 0.;      // NN: for barrel detectors this is 0
    caloData->extent[3] = x_dim.z()/2;

    // This detector is now composed of 12 staves,
    // each stave made of 5 towers which contain each 3 stacks of 26 sensitive layers.
    // See schematics in http://cern.ch/go/MgF6

    // Compute stave dimensions
    double dx = mod_z*(std::tan(hphi) + 1./std::tan(dphi));       // lateral shift
    double trd_x1 = inner_r*std::tan(hphi) + dx/2. - tolerance;   // inner, long side
    double trd_x2 = outer_r*std::sin(hphi) - dx/2. - tolerance;   // outer, short side
    double trd_y1s = x_dim.z()/2. - tolerance;                    // and y2 = y1
    double trd_z  = mod_z/2.;

    // Create the trapezoid for the stave. For coordinate system, see
    // https://root.cern.ch/root/html/TGeoTrd2.html
    // This stave is slightly bigger to contain the support rails
    Trapezoid stv(trd_x1, // Inner side, i.e. "long" X side
                  trd_x2 - 2*ry/std::tan(dphi), // Outer side, i.e. "short" X side
                  trd_y1s,       // Depth at bottom face 
                  trd_y1s,       // Depth at top face y2=y1
                  trd_z + ry);   // Height of the stave, when rotated 
    Volume stave_vol("stave", stv, air);

    std::cout << "Stave height/thickness: " << 2*trd_z << std::endl;
    std::cout << "Got r_inner = " << inner_r << std::endl;
    std::cout << "Got r_outer = " << outer_r << " and r_max = " << r_max << std::endl;
    if(outer_r > r_max) { // throw exception 
      throw std::runtime_error("ERROR: Layers don't fit within the envelope volume!");
    }

    // Create the trapezoid for the towers (5 of them per stave) 
    double trd_y1t = (x_dim.z() - (n_towers-1)*towersAirGap)/n_towers/2 - tolerance;
    Trapezoid twr(trd_x1, // Inner side, i.e. the "long" X side
                  trd_x2, // Outer side, i.e. the "short" X side
                  trd_y1t, // Depth at bottom face
                  trd_y1t, // Depth at top face y2=y1
                  trd_z);  // Height of the tower, when rotated
    
    // Create support rails
    Box rail(ry, trd_y1s, ry); // y=trd_y1t for 3x5

    // Create trapezoids for each of the stacks (3 of them per tower)
    double trd_y1 = (trd_y1t - (n_stacks - 1)*faceThickness)/n_stacks; // add end faces + tolerance
    Trapezoid stk(trd_x1, // Inner side, i.e. the "long" X side
                  trd_x2, // Outer side, i.e. the "short" X side
                  trd_y1, // Depth at bottom face
                  trd_y1, // Depth at top face y2=y1
                  trd_z); // Height of the stack, when rotated

    sens.setType("calorimeter");

    // Start position and angle for the module
    double phi = M_PI/2.0 - M_PI/nsides;                // following the envelope rotation, ECalBarrel_o1_v01_01.xml
    double mod_x_off = dx/2 + tolerance;                // Module X offset
    //double mod_y_off = (inner_r + r_max*cos(hphi))/2; // Module Y offset - mid-envelope placement
    double mod_y_off = inner_r + ry + trd_z + tolerance;// Module Y offset - inner r placement
    std::cout << "Module y offset: " << mod_y_off << std::endl;

    DetElement stave_det("stave0", det_id);    

    for (int t = 0; t < n_towers; t++) {
      double t_pos_y = (t - (n_towers-1)/2)*(2*trd_y1t + towersAirGap);
      string t_name = _toString(t+1, "module%d");
      Volume tower_vol(t_name, twr, Composite); // solid C composite
      DetElement tower_det(stave_det, t_name, det_id);

      for (int m = 0; m < n_stacks; m++) { 
	double m_pos_y = (m - (n_stacks-1)/2)*(2*trd_y1 + faceThickness);
	string k_name = _toString(m+1, "submodule%d");
	Volume stack_vol(k_name, stk, Composite); // solid C composite
        DetElement stack_det(tower_det, k_name, det_id);

        { // Parameters for computing the layer X dimension:
	  double l_dim_y  = trd_y1 - 2.*faceThickness;
	  double l_dim_x  = trd_x1; // Starting X dimension for the layer
	  double tan_beta = std::tan(dphi);
	  double l_pos_z  = -trd_z;
	  
	  // Loop over the sets of layer elements in the module
	  int l_num = 1;
	  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
	    
	    xml_comp_t x_layer = li;
	    int repeat = x_layer.repeat();
	    
	    // Loop over number of repeats for this layer.
	    for (int j=0; j<repeat; j++) {
	      
	      LayeredCalorimeterData::Layer caloLayer ;
	      
	      string l_name = _toString(l_num, "layer%d");
	      double l_thickness = layering.layer(l_num-1)->thickness(); // layer thickness          
	      l_dim_x -= l_thickness/tan_beta;                        // decreasing width 
	      
	      Position   l_pos(0., 0., l_pos_z + l_thickness/2.);     // layer position 
	      Box        l_box(l_dim_x - tolerance, l_dim_y - tolerance, l_thickness/2.);
	      Volume     l_vol(l_name, l_box, air);
	      DetElement layer(stack_det, l_name, det_id);
	      
	      // For caloData
	      double totalAbsorberThickness = 0.;
	      
	      // Loop over the sublayers or slices for this layer
              int s_num = 1;
              double s_pos_z = -(l_thickness/2);

	      double th_i(0.), th_o(-1.) ;

	      for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
		
		xml_comp_t x_slice = si;
		string     s_name  = _toString(s_num, "slice%d");
		double     s_thick = x_slice.thickness();
		Box        s_box(l_dim_x - tolerance, l_dim_y - tolerance, s_thick/2.);
		Material   slice_material  = theDetector.material(x_slice.materialStr());
		Volume     s_vol(s_name, s_box, slice_material);
		DetElement slice(layer, s_name, det_id);

		if ( x_slice.isSensitive() ) {
		  s_vol.setSensitiveDetector(sens);
		  caloLayer.distance = s_pos_z + s_thick/.2 + l_pos_z + l_thickness/2. + mod_y_off - ry;

		  th_i += s_thick / 2. ;
		  th_o  = s_thick / 2. ;

		} else {

		  if( th_o < 0. ){
		    th_i += s_thick;
		  } else {
		    th_o += s_thick;
		  }

		} 
		

		if( x_slice.isRadiator() == true) {
		  totalAbsorberThickness += s_thick;
		}

		slice.setAttributes(theDetector, s_vol, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
		
		// Slice placement
		PlacedVolume slice_pv = l_vol.placeVolume(s_vol, Position(0., 0., s_pos_z + s_thick/2.));
		slice.setPlacement(slice_pv);	  
		// Increment Z position of slice
		s_pos_z += s_thick;
		
		// Increment slice number
		++s_num;
	      }
	      
	      // Set region, limitset, and visibility of layer
	      layer.setAttributes(theDetector, l_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
	      
	      PlacedVolume layer_pv = stack_vol.placeVolume(l_vol, l_pos);
	      layer_pv.addPhysVolID("layer", l_num);
	      layer.setPlacement(layer_pv);
	      
	      caloLayer.inner_thickness = th_i ;
	      caloLayer.outer_thickness = th_o ;
	      caloLayer.absorberThickness = totalAbsorberThickness;
	      caloLayer.cellSize0 = cell_sizeX;
	      caloLayer.cellSize1 = cell_sizeY;

	      caloData->layers.push_back( caloLayer ) ;
	      
	      // Increment to next layer Z position
	      l_pos_z += l_thickness;
	      ++l_num;
	    }
	  }
	}
	//stack_det.setAttributes(theDetector, stack_vol, x_staves.regionStr(), x_staves.limitsStr(), x_staves.visStr());
	stack_vol.setVisAttributes(theDetector.visAttributes(x_staves.visStr())); // Fix this!
	
	PlacedVolume stack_pv = tower_vol.placeVolume(stack_vol, Position(0., m_pos_y, 0.));
	stack_pv.addPhysVolID("submodule", m); 
	stack_det.setPlacement(stack_pv);
      }
      //tower_det.setAttributes(theDetector, tower_vol, x_staves.regionStr(), x_staves.limitsStr(), x_staves.visStr());
      tower_vol.setVisAttributes(theDetector.visAttributes(x_staves.visStr())); // Fix this!
      
      PlacedVolume tower_pv = stave_vol.placeVolume(tower_vol, Position(0., t_pos_y, -ry));
      tower_pv.addPhysVolID("module", t);  
      tower_det.setPlacement(tower_pv);
    }
      
    // Create support rails (move in the loop for 3x5, i.e. 3 rails per tower)
    double rail_spacing = trd_x2/2.;
    for (int r = 0; r < n_rails; r++){
      string r_name = _toString(r, "support_rail%d");
      Volume rail_vol(r_name, rail, Steel);
      DetElement rail_support(stave_det, r_name, det_id);
      double r_pos_y = (r - (n_rails - 1)/2)*rail_spacing;
      //PlacedVolume rail_pv = stave_vol.placeVolume(rail_vol, Position(r_pos_y, t_pos_y, trd_z));//3x5
      PlacedVolume rail_pv = stave_vol.placeVolume(rail_vol, Position(r_pos_y, 0., trd_z));//3
      //rail_vol.setVisAttributes(theDetector.visAttributes(x_staves.visStr())); // fix!
      rail_support.setPlacement(rail_pv);
    }
    
    // Set visualization
    if ( x_staves ) {
        stave_vol.setVisAttributes(theDetector.visAttributes(x_staves.visStr()));
    }
    
    // Place the staves
    for (int i = 0; i < nsides; i++, phi -= dphi) { // i is the stave number
      double s_pos_x = mod_x_off * std::cos(phi) - mod_y_off * std::sin(phi);
      double s_pos_y = mod_x_off * std::sin(phi) + mod_y_off * std::cos(phi);
      double s_pos_z = 0.;
      Transform3D tr(RotationZYX(0., phi, M_PI/2.), Translation3D(-s_pos_x, -s_pos_y, -s_pos_z));
      PlacedVolume pv = envelope.placeVolume(stave_vol, tr);
      pv.addPhysVolID("stave", i);
      pv.addPhysVolID("side", 0); //should set once in parent volume placement
      //pv.addPhysVolID("system",det_id); // not needed (?)
      DetElement sd = i==0 ? stave_det : stave_det.clone(_toString(i, "stave%d"));
      sd.setPlacement(pv);
      sdet.add(sd);
    }    

    // Set envelope volume attributes
    envelope.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
    
    //FOR NOW, USE A MORE "SIMPLE" VERSION OF EXTENSIONS, INCLUDING NECESSARY GEAR PARAMETERS
    //Copied from Frank's SHcalSc04 Implementation
    sdet.addExtension< LayeredCalorimeterData >( caloData ) ;
    
    return sdet;
}

DECLARE_DETELEMENT(ECalBarrel_o1_v01,create_detector)

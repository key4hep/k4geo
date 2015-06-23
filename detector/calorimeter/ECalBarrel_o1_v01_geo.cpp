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
using namespace DD4hep;
using namespace DD4hep::Geometry;

static Ref_t create_detector(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
    
    xml_det_t     x_det     = e;
    int           det_id    = x_det.id();
    string        det_name  = x_det.nameStr();
    DetElement    sdet      (det_name, det_id);
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = XML::createPlacedEnvelope(lcdd, e, sdet);
    
    if(lcdd.buildType() == BUILD_ENVELOPE) return sdet;
    
    //-----------------------------------------------------------------------------------

    DD4hep::XML::Component xmlParameter = x_det.child(_Unicode(parameter));
    const double tolerance = xmlParameter.attr<double>(_Unicode(ecal_barrel_tolerance));
    const int n_towers = (int) xmlParameter.attr<double>(_Unicode(num_towers));
    const double towersAirGap = xmlParameter.attr<double>(_Unicode(TowersAirGap));
    const double faceThickness = xmlParameter.attr<double>(_Unicode(TowersFaceThickness));
    Material Composite = lcdd.material(xmlParameter.attr<string>(_Unicode(AlveolusMaterial)));
    Readout readout = sens.readout();
    Segmentation seg = readout.segmentation();
    
    std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
    double cell_sizeX      = cellSizeVector[0];
    double cell_sizeY      = cellSizeVector[1];
    
    Layering      layering (e);
    Material      air       = lcdd.air();
    xml_comp_t    x_staves  = x_det.staves();
    xml_comp_t    x_dim     = x_det.dimensions();
    int           nsides    = x_dim.numsides();
    double        inner_r   = x_dim.rmin();                     // inscribing circle radius
    double        dphi      = (2.*M_PI/nsides);
    double        hphi      = dphi/2.;
    double        mod_z     = layering.totalThickness();
    double        outer_r   = (inner_r + mod_z)/std::cos(hphi); // circumscribing circle radius
    double        r_max     = x_dim.rmax()/std::cos(hphi);      // see http://cern.ch/go/r9mZ
    
    // Create caloData object to extend driver with data required for reconstruction
    // Added by Nikiforos, 28/05/15. Should not interfere with driver but it does
    DDRec::LayeredCalorimeterData* caloData = new DDRec::LayeredCalorimeterData ;
    caloData->layoutType = DDRec::LayeredCalorimeterData::BarrelLayout ;
    caloData->inner_symmetry = nsides;
    caloData->outer_symmetry = nsides; 
    caloData->phi0 = 0.; //FIXME
    /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
    caloData->extent[0] = inner_r ;
    caloData->extent[1] = x_dim.rmax(); // check !!
    caloData->extent[2] = 0. ;
    caloData->extent[3] = x_dim.z()/2.0 ;
    
    // The barrel is composed of 12 x 5 towers, each with 17 + 8 layers made of 7 slices
    // One could group 5 towers per stave and split each tower into 3 stacks

    std::cout << "Tower height/thickness: " << mod_z << std::endl;
    std::cout << "Got r_out = " << outer_r << " and r_max = " << r_max << std::endl;
    if(outer_r > r_max) { // throw exception 
      throw std::runtime_error("ERROR: Layers don't fit within the envelope volume!");
    }

    DetElement tower_det("tower0", det_id);
    
    // Compute stave dimensions
    double dx = mod_z*(std::tan(hphi) + 1./std::tan(dphi));               // lateral shift
    double trd_x1 = inner_r*std::tan(hphi) + dx/2. - tolerance;           // inner, long side
    double trd_x2 = outer_r*std::sin(hphi) - dx/2. - tolerance;           // outer, short side
    double trd_y1 = x_dim.z()/2. - tolerance;                             // and y2 = y1
    double trd_z  = mod_z/2.;
    
    /* Create the trapezoid for the stave. For coordinate system, see
    // https://root.cern.ch/root/html/TGeoTrd2.html    
    Trapezoid stv(trd_x1, // Inner side, i.e. the "long" X side
                  trd_x2, // Outer side, i.e. the "short" X side
                  trd_y1, // Depth at bottom face 
                  trd_y1, // Depth at top face y2=y1
                  trd_z); // Height of the stave, when rotated 
    Volume stave_vol("stave", stv, air);     // stave volume */

    // Create the trapezoid for the towers
    trd_y1 = (x_dim.z() - (n_towers-1)*towersAirGap)/n_towers/2. - tolerance;
    Trapezoid twr(trd_x1, // Inner side, i.e. the "long" X side
                  trd_x2, // Outer side, i.e. the "short" X side
                  trd_y1, // Depth at bottom face
                  trd_y1, // Depth at top face y2=y1
                  trd_z); // Height of the tower, when rotated
    Volume tower_vol("tower", twr, Composite); // solid tungsten volume
        
    sens.setType("calorimeter");

    {// Build the tower, layer by layer and slice by slice

        // Parameters for computing the layer X dimension:
        double l_dim_y  = trd_y1 - 2.*faceThickness;
        double l_dim_x  = trd_x1; // Starting X dimension for the layer
        double tan_beta = std::tan(dphi);
        double l_pos_z  = -trd_z;
        
        // Loop over the sets of layer elements in the detector
        int l_num = 1;
        for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
            
            xml_comp_t x_layer = li;
            int repeat = x_layer.repeat();

            // Loop over number of repeats for this layer.
            for (int j=0; j<repeat; j++)  {

                DDRec::LayeredCalorimeterData::Layer caloLayer ;
              
                string l_name = _toString(l_num, "layer%d");
                double l_thickness = layering.layer(l_num-1)->thickness();  // this layer's thickness          
		//std::cout << l_name << " thickness: " << l_thickness << std::endl;
                l_dim_x -= l_thickness/tan_beta;                      // decreasing width 

                Position   l_pos(0., 0., l_pos_z + l_thickness/2.);   // position of the layer
                Box        l_box(l_dim_x - tolerance, l_dim_y - tolerance, l_thickness/2.);
                Volume     l_vol(l_name, l_box, air);
                DetElement layer(tower_det, l_name, det_id);
                
                // Loop over the sublayers or slices for this layer
                int s_num = 1;
                double s_pos_z = -(l_thickness/2);
                double totalAbsorberThickness=0.;
                
                for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {

                    xml_comp_t x_slice = si;
                    string     s_name  = _toString(s_num,"slice%d");
                    double     s_thick = x_slice.thickness();
		    //std::cout << s_name << " thickness: " << s_thick <<std::endl;
                    Box        s_box(l_dim_x - tolerance, l_dim_y - tolerance, s_thick/2.);
                    Volume     s_vol(s_name, s_box, lcdd.material(x_slice.materialStr()));
                    DetElement slice(layer, s_name, det_id);
                    
                    if ( x_slice.isSensitive() ) {
                        s_vol.setSensitiveDetector(sens);
                    }
                    
                    char val = x_slice.hasAttr(_U(radiator)) ? x_slice.attr < string > (_U(radiator))[0] : 'f';
                    val = std::toupper(val);
                    bool isAbsorber =  (val == 'T' || val == 'Y');
                    
                    if( isAbsorber ==true)
                      totalAbsorberThickness+= s_thick;
                    
                    slice.setAttributes(lcdd, s_vol, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

                    // Slice placement
                    PlacedVolume slice_phv = l_vol.placeVolume(s_vol, Position(0., 0., s_pos_z + s_thick/2.));
                    slice.setPlacement(slice_phv);
                    // Increment Z position of slice
                    s_pos_z += s_thick;
                    
                    // Increment slice number
                    ++s_num;
                }
                
                // Set region, limitset, and visibility of layer
                layer.setAttributes(lcdd, l_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
                
                PlacedVolume layer_phv = tower_vol.placeVolume(l_vol, l_pos);
                layer_phv.addPhysVolID("layer", l_num);
                layer.setPlacement(layer_phv);
                
                caloLayer.distance = l_pos_z + l_thickness/2.;
                caloLayer.thickness = l_thickness;
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

    // Set visualization
    if ( x_staves )   {
        tower_vol.setVisAttributes(lcdd.visAttributes(x_staves.visStr()));
    }
    
    // Phi start for staves/towers
    double phi = M_PI/nsides;
    double mod_x_off = dx/2 + tolerance;               // Stave/tower X offset
    double mod_y_off = (inner_r + outer_r)/2;          // Stave/tower Y offset
    
    // Create the towers
    for (int t = 0; t < n_towers; t++)  {
        // Create nsides staves
        for (int i = 0; i < nsides; i++, phi -= dphi) { // i is the stave number
            // Compute the tower position in envelope coordinates
            double m_pos_x = mod_x_off * std::cos(phi) - mod_y_off * std::sin(phi);
            double m_pos_y = mod_x_off * std::sin(phi) + mod_y_off * std::cos(phi);
            double m_pos_z = -x_dim.z()/2 + (2*t+1)*(trd_y1 + tolerance) + t*towersAirGap;
            Transform3D tr(RotationZYX(0., phi, M_PI/2.), Translation3D(-m_pos_x, -m_pos_y, -m_pos_z));
            PlacedVolume pv = envelope.placeVolume(tower_vol, tr);
            int mid = (nsides+1)*t+i; 
            //pv.addPhysVolID("system",det_id); // not needed
            pv.addPhysVolID("side", 0); //should set once in parent volume placement
            pv.addPhysVolID("module", t); ///FIXME! DRIVER DEVELOPER PLEASE CHECK IF ASSIGNMENTS MAKE SENSE
            pv.addPhysVolID("stave", i);
            DetElement sd = mid==0 ? tower_det : tower_det.clone(_toString(mid,"tower%d"));
            sd.setPlacement(pv);
            sdet.add(sd);
	}
    }
    // Set envelope volume attributes
    envelope.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
    
    //FOR NOW, USE A MORE "SIMPLE" VERSION OF EXTENSIONS, INCLUDING NECESSARY GEAR PARAMETERS
    //Copied from Frank's SHcalSc04 Implementation
    sdet.addExtension< DDRec::LayeredCalorimeterData >( caloData ) ;
    
    return sdet;
}

DECLARE_DETELEMENT(ECalBarrel_o1_v01,create_detector)

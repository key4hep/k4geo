// Silicon Tracker Barrel implementation for the CLIC detector
//====================================================================
//--------------------------------------------------------------------
//
//  Author     : N. Nikiforou (forked from SiTrackerBarrel_geo.cpp
//
// Comment: Here each slice of the module has the same transversal
//          dimensions (x,y). Suitable for Surfaces.
//
// Comment: You have to use <include ref=... to define a module stack
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Utilities.h"
#include "XML/DocumentHandler.h"
#include "DDRec/DetectorData.h"

#include <UTIL/BitField64.h>
#include <UTIL/BitSet32.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

using dd4hep::Assembly;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::ERROR;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationX;
using dd4hep::RotationZ;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Tube;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::rec::NeighbourSurfacesData;
using dd4hep::rec::ZPlanarData;

void populateNeighbourData(NeighbourSurfacesData* neighbourSurfacesData, UTIL::BitField64& encoder, int module_idx, int sensor_idx, int nphi, int nz) {
    const dd4hep::CellID cellID = encoder.lowWord(); // 32 bits

    //compute neighbours 
    int n_neighbours_module = 1; // 1 gives the adjacent modules (i do not think we would like to change this)
    int n_neighbours_sensor = 1;

    int newmodule=0, newsensor=0;

    for(int imodule=-n_neighbours_module; imodule<=n_neighbours_module; imodule++){ // neighbouring modules
      for(int isensor=-n_neighbours_sensor; isensor<=n_neighbours_sensor; isensor++){ // neighbouring sensors
        
        if (imodule==0 && isensor==0) continue; // cellID we started with
        newmodule = module_idx + imodule;
        newsensor = sensor_idx + isensor;

        //compute special case at the boundary  
        //general computation to allow (if necessary) more then adjacent neighbours (ie: +-2)
        
        if (newmodule < 0) newmodule = nphi + newmodule;
        if (newmodule >= nphi) newmodule = newmodule - nphi;

        if (newsensor < 0 || newsensor >= nz) continue; //out of the stave

        //encoding
        encoder[lcio::LCTrackerCellID::module()] = newmodule;
        encoder[lcio::LCTrackerCellID::sensor()] = newsensor;
        
        neighbourSurfacesData->sameLayer[cellID].push_back(encoder.lowWord());

      }
    }

}

ZPlanarData::LayerLayout buildLayerLayout(double radius, int staves, double phi0, double z0, Volume module, PlacedVolume sensor) {
    ZPlanarData::LayerLayout thisLayer ;
    /// GET GEAR INFORMATION
    /// NOTE WORKS ONLY FOR ONE WAFER
    /// FIXME: does not take phi_dr and z_dr into account.
    Box mod_shape(module.solid()), comp_shape(sensor.volume().solid());

    const double* trans = sensor->GetMatrix()->GetTranslation();
    double half_module_thickness = mod_shape->GetDZ();
    double half_silicon_thickness = comp_shape->GetDZ();

    double sensitive_z_position  = trans[2];

    double inner_thickness = half_module_thickness - sensitive_z_position;

    thisLayer.distanceSupport  = radius;

    thisLayer.offsetSupport    =  0; 
    thisLayer.thicknessSupport = inner_thickness- half_silicon_thickness;
    thisLayer.zHalfSupport    = z0 + mod_shape->GetDY();
    thisLayer.widthSupport     = 2*mod_shape->GetDX(); 

    thisLayer.distanceSensitive = radius + sensitive_z_position; 
    thisLayer.offsetSensitive  = 0. ;
    thisLayer.thicknessSensitive = 2*half_silicon_thickness;//Assembled along Z
    //Changed by Thorben Quast (same applies to zHalfSupport)
    //z0 = center of most right sensor, comp_shape-GetDY() = half length of one sensitive are of the module
    thisLayer.zHalfSensitive    = z0 + comp_shape->GetDY();
    thisLayer.widthSensitive = 2*comp_shape->GetDX();
    thisLayer.ladderNumber = staves;
    thisLayer.phi0 =  phi0;
    
    return thisLayer;
}

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens)  {
    typedef std::vector<PlacedVolume> Placements;
    xml_det_t   x_det     = e;
    Material    air       = theDetector.air();
    int         det_id    = x_det.id();
    std::string      det_name  = x_det.nameStr();
    DetElement  sdet       (det_name,det_id);
    // Assembly    assembly   (det_name);
    std::map<std::string, Volume>    volumes;
    std::map<std::string, Placements>  sensitives;
    PlacedVolume pv;
    

    // for encoding
    std::string cellIDEncoding = sens.readout().idSpec().fieldDescription();
    UTIL::BitField64 encoder( cellIDEncoding );
    encoder.reset();
    encoder[lcio::LCTrackerCellID::subdet()] = det_id;
    encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::barrel;


    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;
    
    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    ZPlanarData*  zPlanarData = new ZPlanarData() ;
    NeighbourSurfacesData*  neighbourSurfacesData = new NeighbourSurfacesData() ;
    
    sens.setType("tracker");
    
    // TODO: refactor module creation code into a function...
    //NOTE modules are what is defined in compact. Later we call a "module" as a "sensor".
    for(xml_coll_t mi(x_det,_U(module)); mi; ++mi)  {
        xml_comp_t x_mod  = mi;
        xml_comp_t m_env  = x_mod.child(_U(module_envelope));
        std::string     m_nam  = x_mod.nameStr();
        
        if ( volumes.find(m_nam) != volumes.end() )   {
            printout(ERROR,"TrackerBarrel","Logics error in building modules.");
            throw std::runtime_error("Logic error in building modules.");
        }

        double module_thickness = 0;
        for(xml_coll_t incl(x_mod,_U(include)); incl; ++incl) {
            dd4hep::xml::DocumentHolder doc(dd4hep::xml::DocumentHandler().load(incl, incl.attr_value(_U(ref))));
            xml_h includes = doc.root();
            xml_det_t incl_stack = includes;
            for (xml_coll_t ci(incl_stack, _U(module_component)); ci; ++ci) {
                xml_comp_t x_comp = ci;
                module_thickness = module_thickness + x_comp.thickness();
            }
        }

        Volume     m_vol(m_nam,Box(m_env.width()/2.,m_env.length()/2.,module_thickness/2.),air);
        volumes[m_nam] = m_vol;
        m_vol.setVisAttributes(theDetector.visAttributes(x_mod.visStr()));
        
        
        int        ncomponents = 0; 
        
        //First component on top of the list is the innermost one. 
        double position_z= -module_thickness/2.;
        for(xml_coll_t incl(x_mod,_U(include)); incl; ++incl) {
            dd4hep::xml::DocumentHolder doc(dd4hep::xml::DocumentHandler().load(incl, incl.attr_value(_U(ref))));
            xml_h includes = doc.root();
            xml_det_t incl_stack = includes;
            for (xml_coll_t ci(incl_stack, _U(module_component)); ci; ++ci, ++ncomponents) {
                xml_comp_t x_comp = ci;
                std::string c_nam = _toString(ncomponents, "component%d");
                Box c_box(m_env.width() / 2.0, m_env.length() / 2.0, x_comp.thickness() / 2.0);
                Volume c_vol(c_nam, c_box, theDetector.material(x_comp.materialStr()));


                pv = m_vol.placeVolume(c_vol, Position(0, 0, position_z + x_comp.thickness() / 2.0));

                c_vol.setRegion(theDetector, x_comp.regionStr());
                c_vol.setLimitSet(theDetector, x_comp.limitsStr());
                c_vol.setVisAttributes(theDetector, x_comp.visStr());
                if (x_comp.isSensitive()) {
                    //         pv.addPhysVolID("wafer",wafer_number++);
                    c_vol.setSensitiveDetector(sens);
                    sensitives[m_nam].push_back(pv);
                }

                position_z += x_comp.thickness();
            }
        }
    }
    for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
        xml_comp_t x_layer  = li;
        xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
        xml_comp_t z_layout = x_layer.child(_U(z_layout));      // Get the <z_layout> element.
        int        lay_id   = x_layer.id();
        // int        type     = x_layer.type();
        std::string     m_nam    = x_layer.moduleStr();
        std::string     lay_nam  = _toString(x_layer.id(),"layer%d");

        Assembly   lay_vol   (lay_nam);         // Create the layer envelope volume.
        double     phi0     = x_layout.phi0();              // Starting phi of first sensor.
        double     phi_tilt = x_layout.phi_tilt();          // Phi tilt of a sensor.
        double     rc       = x_layout.rc();                // Radius of the sensor center.
        int        nphi     = x_layout.nphi();              // Number of sensors in phi.
        double     rphi_dr  = x_layout.dr();                // The delta radius of every other sensor.
        
        double     phi_incr = (M_PI * 2) / nphi;            // Phi increment for one sensor.
        double     z0       = z_layout.z0();                // Z position of first sensor in phi.
        double     nz       = z_layout.nz();                // Number of sensors to place in z.
        double     z_dr     = z_layout.dr();                // Radial displacement parameter, of every other sensor.
        Volume sensorVol = volumes[m_nam];
        DetElement lay_elt(sdet,_toString(x_layer.id(),"layer%d"),lay_id);
        Placements& waferVols = sensitives[m_nam];

        ZPlanarData::LayerLayout thisLayer = buildLayerLayout(rc, nphi, phi0, z0, sensorVol, waferVols[0]);
        zPlanarData->layers.push_back(thisLayer);

        Assembly stave("stave");
        // Z increment for sensor placement along Z axis.
        // Adjust for z0 at center of sensor rather than
        // the end of cylindrical envelope.
        DetElement staveElementTemplate("staveElementTemplate", 0);
        double z_incr = nz > 1 ? (2.0 * z0) / (nz - 1) : 0.0;
        for (int i = 0; i < nz; i++)
        {
            // sensor/module still lives in a world where:
            // x: width, y: length, z: thickness
            // we construct the stave in the same way,
            // i.e. the modules are placed along y from z0 to -z0
            // and dr is applied along z
            double y = z0 - i * z_incr;
            double z = i % 2 == 0 ? 0 : z_dr;
            Position pos(0., y, z);
            auto sensorPv = stave.placeVolume(sensorVol, pos);
            sensorPv.addPhysVolID("sensor", i);
            std::string sensorName = _toString(i, "sensor%d");
            DetElement sensorElement(sensorName, i);
            staveElementTemplate.add(sensorElement);
            sensorElement.setPlacement(sensorPv);


            for (const PlacedVolume& wafer_pv : waferVols) {
                DetElement comp_elt(sensorElement, wafer_pv.volume().name(), i);
                comp_elt.setPlacement(wafer_pv);
            }
        }

        // Loop over the number of staves in phi.
        for (int i = 0; i < nphi; i++)
        {
            double dr = i % 2 == 0 ? 0 : rphi_dr;
            double r = rc + dr;
            double phi = phi0 + i * phi_incr;
            double x = r * std::cos(phi); // Basic x stave position.
            double y = r * std::sin(phi); // Basic y stave position.

            // module: stave
            std::string moduleName = _toString(i, "module%d");
            // stave has x: width, y: length, z: thickness (pointing outwards)
            // rotate -pi/2 around x to align length with z
            // rotate -p/2 * phi around z to point outwards again
            Rotation3D rotation = RotationZ(-M_PI/2 + phi + phi_tilt) * RotationX(-M_PI/2);
            Transform3D transform(rotation, Position(x, y, 0.));
            auto stavePv = lay_vol.placeVolume(stave, transform);
            stavePv.addPhysVolID("module", i);
            DetElement staveElement = staveElementTemplate.clone(moduleName, i);
            staveElement.setPlacement(stavePv);
            lay_elt.add(staveElement);

            for (int j = 0; j < nz; j++)
            {

                ///////////////////

                //get cellID and fill map< cellID of surface, vector of cellID of neighbouring surfaces >

                //encoding

                encoder[lcio::LCTrackerCellID::layer()] = lay_id;
                encoder[lcio::LCTrackerCellID::module()] = i;
                encoder[lcio::LCTrackerCellID::sensor()] = j;

                populateNeighbourData(neighbourSurfacesData, encoder, i, j, nphi, nz);


                ///////////////////

            }

        }
        // Create the PhysicalVolume for the layer.
        pv = envelope.placeVolume(lay_vol); // Place layer in mother
        pv.addPhysVolID("layer", lay_id);       // Set the layer ID.
        lay_elt.setAttributes(theDetector,lay_vol,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());
        lay_elt.setPlacement(pv);

    }
    sdet.setAttributes(theDetector,envelope,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
    sdet.addExtension< ZPlanarData >( zPlanarData ) ;
    sdet.addExtension< NeighbourSurfacesData >( neighbourSurfacesData ) ;
    
    //envelope.setVisAttributes(theDetector.invisible());
    /*pv = theDetector.pickMotherVolume(sdet).placeVolume(assembly);
     pv.addPhysVolID("system", det_id);      // Set the subdetector system ID.
     pv.addPhysVolID("barrel", 0);           // Flag this as a barrel subdetector.
     sdet.setPlacement(pv);*/
    return sdet;
}

DECLARE_DETELEMENT(TrackerBarrel_o1_v06, create_detector)

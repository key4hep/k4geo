//====================================================================
//  Vertex Detector implementation for the FCC-ee IDEA detector
//--------------------------------------------------------------------
//
//  Based on VertexEndcap_o2_v06_geo.cpp from M. Petric, which was 
//  originaly forked form SiTrackerEndcap2 by M. Frank
//
//  Author     : A. Ilg
//
//  This code allows to build a tracker/vertex endcap made out of staves
//  as it is used in the IDEA vertex detector design by F. Palla and 
//  F. Bosi as of early-2023.
//  The staves are arranged in petals, and can feature any number of modules.
//  The modules can be built by smaller rectangular structures to represent
//  both sensitive and insensitive (periphery) parts so that e.g Quad
//  modules can be bult.
//====================================================================

#include "DD4hep/DetFactoryHelper.h"
#include <map>
#include "XML/Utilities.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"

using namespace std;

using dd4hep::Assembly;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::ERROR;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Trapezoid;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::rec::ZDiskPetalsData;
using dd4hep::Box;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens)  {
    typedef vector<PlacedVolume> Placements;
    xml_det_t   x_det     = e;
    Material    vacuum    = theDetector.vacuum();
    string      det_name  = x_det.nameStr();
    bool        reflect   = x_det.reflect(false);
    DetElement  sdet        (det_name,x_det.id());
    int         m_id=0;
    map<string,Volume> modules;
    map<string, Placements>  sensitives;
    PlacedVolume pv;
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    
    
    envelope.setVisAttributes(theDetector.invisible()); sens.setType("tracker");

    // -------- reconstruction parameters  ----------------
    ZDiskPetalsData*  zDiskPetalsData = new ZDiskPetalsData ;
    std::map< std::string, double > moduleSensThickness;

    // --- Module information struct ---
    struct module_information{
        string name;

        double readout_z_offset;
        vector<double> readout_thicknesses;
        vector<double> readout_widths;
        vector<double> readout_offsets; 
        vector<double> readout_z_offsets; 
        vector<Material> readout_materials;
        vector<string> readout_viss;

        double sensor_z_offset;
        double sensor_thickness;
        Material sensor_material;
        vector<bool> sensor_sensitives;
        vector<double> sensor_xmin;
        vector<double> sensor_xmax;
        vector<double> sensor_ymin;
        vector<double> sensor_ymax;
        double sensor_width;
        double sensor_length;
        vector<string> sensor_viss;
        vector<Volume> sensor_volumes;
        Volume m_volume;

        double support_z_offset;
        vector<double> support_thicknesses;
        vector<double> support_offsets;
        vector<double> support_z_offsets;
        vector<double> support_widths;
        vector<Material> support_materials;
        vector<string> support_viss;
    };
    list<module_information> module_information_list;

    // --- Collect module(s) information
    for(xml_coll_t mi(x_det,_U(module)); mi; ++mi, ++m_id)  {
        xml_comp_t x_mod   = mi;

        module_information m;
        m.name = x_mod.nameStr();

        // Readout
        xml_coll_t c_readout(x_mod,_U(readout));
        m.readout_z_offset = xml_comp_t(c_readout).z_offset();
        xml_coll_t c_component(c_readout,_U(component));
        for(c_component.reset(); c_component; ++c_component){
            xml_comp_t component = c_component;
            m.readout_thicknesses.push_back(component.thickness());
            m.readout_widths.push_back(component.width());
            m.readout_offsets.push_back(component.offset());
            m.readout_z_offsets.push_back(component.z_offset());
            m.readout_materials.push_back(theDetector.material(component.materialStr()));
            m.readout_viss.push_back(component.visStr());
        }

        // Sensor
        xml_coll_t c_sensor(x_mod,_U(sensor));
        m.sensor_z_offset = xml_comp_t(c_sensor).z_offset();
        m.sensor_thickness = xml_comp_t(c_sensor).thickness();
        m.sensor_material = theDetector.material(xml_comp_t(c_sensor).materialStr());
        c_component = xml_coll_t(c_sensor,_U(component));

        int iComponent=0;
        for(c_component.reset(); c_component; ++c_component){
            xml_comp_t component = c_component;
            m.sensor_sensitives.push_back(component.isSensitive());
            m.sensor_xmin.push_back(component.xmin());
            m.sensor_xmax.push_back(component.xmax());
            m.sensor_ymin.push_back(component.ymin());
            m.sensor_ymax.push_back(component.ymax());
            m.sensor_viss.push_back(component.visStr());

            // Already create volumes for all sensor components as this is independent of number of sensors per layer
            Box ele_box = Box( abs(component.xmax()-component.xmin())/2., abs(component.ymax()-component.ymin())/2., m.sensor_thickness/2.);
            m.sensor_volumes.push_back( Volume(m.name + _toString(iComponent, "_sensor_%d"), ele_box, m.sensor_material) );
            iComponent++;
        }
        m.sensor_width  = *max_element(m.sensor_xmax.begin(), m.sensor_xmax.end()) - *min_element(m.sensor_xmin.begin(), m.sensor_xmin.end());
        m.sensor_length = *max_element(m.sensor_ymax.begin(), m.sensor_ymax.end()) - *min_element(m.sensor_ymin.begin(), m.sensor_ymin.end());
        cout << "Module: " << m.name << ", sensor width: " << to_string(m.sensor_width)  << ", sensor length: " << to_string(m.sensor_length) << endl;

        Volume  m_volume(m.name, Box( m.sensor_width/2.0, m.sensor_length/2.0, m.sensor_thickness), vacuum);
        m_volume.setVisAttributes(theDetector.visAttributes(m.sensor_viss[0]));
        for(int i=0; i<m.sensor_volumes.size(); i++){
            double x_pos = m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.;
            double y_pos = m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2.;
            double z_pos = 0;
            pv = m_volume.placeVolume(m.sensor_volumes[i],Position(x_pos, y_pos, z_pos));
            if(m.sensor_sensitives[i]) {
                m.sensor_volumes[i].setSensitiveDetector(sens);
                sensitives[m.name].push_back(pv);
                moduleSensThickness[m.name] = m.sensor_thickness;
            }
        }
        m.m_volume = m_volume;


        // Support
        xml_coll_t c_support(x_mod,_U(support));
        m.support_z_offset = xml_comp_t(c_support).z_offset();
        c_component = xml_coll_t(c_support,_U(component));
        for(c_component.reset(); c_component; ++c_component){
            xml_comp_t component = c_component;
            m.support_thicknesses.push_back(component.thickness());
            m.support_offsets.push_back(component.offset());
            m.support_z_offsets.push_back(component.z_offset());
            m.support_widths.push_back(component.width());
            m.support_materials.push_back(theDetector.material(component.materialStr()));
            m.support_viss.push_back(component.visStr());
        }
        module_information_list.push_back(m);
    }
   
    vector<int> sides = {1};
    if(reflect){sides.push_back(-1);}
    for(auto & side : sides){
        string side_name = det_name + _toString(side,"_side%d");
        Assembly side_assembly(side_name);
        pv = envelope.placeVolume(side_assembly);
        pv.addPhysVolID("side", side);


        for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
            xml_comp_t  x_layer(li);
            int layer_id            = x_layer.id();
            double rmin         = x_layer.rmin();
            double dr           = x_layer.dr();
            double z            = x_layer.z();
            double layer_dz     = x_layer.dz();
            int nPetals         = x_layer.nPetals();
            double phi0_layer   = x_layer.phi0();

            int mod_num = 0;
                       
            // -------- reconstruction parameters  ----------------
            //NOTE: Mostly Dummy information for event display/DDMarlinPandora
            //FIXME: have to see how to properly accommodate spirals and the fact that there's only
            //one ring
            ZDiskPetalsData::LayerLayout thisLayer ;

            int numberOfRings=0; //check that only one ring is used
            
            string layer_name = side_name + _toString(layer_id,"_layer%d");
            Assembly layer_assembly(layer_name);
            // DetElement layerDE( sdet , layer_name, x_det.id() );
            pv = side_assembly.placeVolume(layer_assembly);
            pv.addPhysVolID("layer", layer_id );  
            // layerDE.setPlacement( pv ) ;

            for(int iPetal=0; iPetal<nPetals; iPetal++){
                double z_alternate_petal = (iPetal%2 == 0) ? 0.0 : layer_dz;

                string petal_name = layer_name + _toString(iPetal,"_petal%d");
                Assembly petal_assembly(petal_name);
                pv = layer_assembly.placeVolume(petal_assembly);

                int iStave = 0;
                int nStaves = 0;
                for(xml_coll_t ri(x_layer,_U(stave)); ri; ++ri)  
                    nStaves+=1;
        
                for(xml_coll_t ri(x_layer,_U(stave)); ri; ++ri)  {
                    if(numberOfRings>0){
                        printout(ERROR,"VertexEndcap","Driver (and ZDiskPetalsData structure) does not support more than one ring per layer!");
                        throw runtime_error("More than one ring per layer not supported by driver.");
                    }
                        
                    xml_comp_t x_stave = ri;

                    int    nmodules     = x_stave.nmodules();
                    double r_stave      = x_stave.r();
                    double z_offset     = x_stave.z_offset();
                    double stave_dz     = x_stave.dz(0);
                    double step         = x_stave.step();   // Spacing of modules
                    string moduleStr    = x_stave.moduleStr();
                    double phi0_stave   = x_stave.phi0();
    //                string m_nam    = x_ring.moduleStr();
    //                Volume m_vol    = modules[m_nam];

                    double phi     = 2*M_PI/nPetals*iPetal + phi0_layer + phi0_stave;

                    // Use the correct module
                    auto m = *find_if(module_information_list.cbegin(), module_information_list.cend(), [&moduleStr] (const module_information& module) {
                        return module.name == moduleStr;
                    });
            
                    Placements& sensVols = sensitives[m.name];

                    string stave_name = petal_name + _toString(iStave,"_stave%d");
                    Assembly stave_assembly(stave_name);
                    pv = petal_assembly.placeVolume(stave_assembly);
                    
                    // Place all components
                    RotationZYX rot( phi , 0, 0  );

                    double stave_length = nmodules*m.sensor_length + (nmodules-1)*step;
                    double r = rmin + m.sensor_width/2.0 + r_stave + ((iPetal%2 == 0) ? 0.0 : dr);

                    
                    // Place support
                    for(int i=0; i<m.support_thicknesses.size(); i++){
                        double x_pos = (r + m.support_offsets[i])*cos(phi);
                        double y_pos = r*sin(phi);
                        double z_pos = z + z_alternate_petal + z_offset + m.support_z_offset + m.support_z_offsets[i]; 
                        if(side == -1){z_pos = -z_pos;}
                        Position pos(x_pos, y_pos, z_pos);

                        Box ele_box = Box( m.support_widths[i]/2., stave_length/2., m.support_thicknesses[i]/2.);
                        Volume ele_vol = Volume( _toString(i, "support_%d"), ele_box, m.support_materials[i]);                    
                        ele_vol.setVisAttributes(theDetector.visAttributes(m.support_viss[i]));

                        pv = stave_assembly.placeVolume(ele_vol, Transform3D(rot, pos) );
                    }

                    // Place readout
                    for(int i=0; i<m.readout_thicknesses.size(); i++){
                        double x_pos = (r + m.readout_offsets[i])*cos(phi);
                        double y_pos = r*sin(phi);
                        double z_pos = z + z_alternate_petal + z_offset + m.readout_z_offset + m.readout_z_offsets[i];
                        if(side == -1){z_pos = -z_pos;}
                        Position pos(x_pos, y_pos, z_pos);

                        Box ele_box = Box( m.readout_widths[i]/2., stave_length/2., m.readout_thicknesses[i]/2.);
                        Volume ele_vol = Volume( _toString(i, "readout_%d"), ele_box, m.readout_materials[i]);                    
                        ele_vol.setVisAttributes(theDetector.visAttributes(m.readout_viss[i]));
                        pv = stave_assembly.placeVolume(ele_vol, Transform3D(rot, pos) );
                    }

                    for(int iModule=0; iModule<nmodules; iModule++){
                        double z_alternate_module = (iModule%2 == 0) ? 0.0 : stave_dz;
                        double x_pos = r*cos(phi) - (-(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + iModule*m.sensor_length + iModule*step)*sin(phi);
                        double y_pos = r*sin(phi) + (-(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + iModule*m.sensor_length + iModule*step)*cos(phi);
                        double z_pos = z + z_alternate_petal + z_offset + m.sensor_z_offset + z_alternate_module; 
                        if(side == -1){z_pos = -z_pos;}
                        Position pos(x_pos, y_pos, z_pos);

                        int iSensor=0;
                        string module_name = stave_name + _toString(mod_num,"_module%d");
                        DetElement module(sdet,module_name,x_det.id());
                        pv = envelope.placeVolume(m.m_volume,Transform3D(rot, pos));
                        pv.addPhysVolID("side",side).addPhysVolID("layer", layer_id).addPhysVolID("module",mod_num).addPhysVolID("sensor",iSensor);
                        module.setPlacement(pv);

                        for(size_t ic=0; ic<sensVols.size(); ++ic)  {
                            PlacedVolume sens_pv = sensVols[ic];
                            DetElement comp_elt(module,sens_pv.volume().name(),mod_num);
                            comp_elt.setPlacement(sens_pv);
                        }                        
                        // Assembly module_assembly(module_name);
                        // pv = layer_assembly.placeVolume(module_assembly);

                    //    // Place sensors: non-sensitive parts
                    //     for(int i=0; i<m.sensor_sensitives.size(); i++){
                    //         if(m.sensor_sensitives[i]){
                    //             continue;
                    //         }
                    //         Box ele_box = Box( abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2., abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2., m.sensor_thickness/2.);
                    //         Volume sensor_volume (m.name + _toString(i, "_sensorPassive_%d"), ele_box, m.sensor_material);

                    //         if(side == -1){z_pos = -z_pos;}
                    //         Position pos(x_pos, y_pos, z_pos);
                    //         sensor_volume.setVisAttributes(theDetector.visAttributes(m.sensor_viss[i]));
                    //         pv = module_assembly.placeVolume(sensor_volume, Transform3D(rot, pos) );

                    //         pv.addPhysVolID("side",side).addPhysVolID("layer", layer_id).addPhysVolID("module",mod_num).addPhysVolID("sensor",iSensor);
                    //         module.setPlacement(pv);

                    //     }
                    
                    //     // Place sensors: sensitive parts
                    //     for(int i=0; i<m.sensor_sensitives.size(); i++){
                    //         if(m.sensor_sensitives[i]==false){
                    //             continue;
                    //         }

                    //         string sensor_name = module_name + _toString(iSensor,"_sensor%d");

                    //         Box ele_box = Box( abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2., abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2., m.sensor_thickness/2.);
                    //         Volume sensor_volume (sensor_name, ele_box, m.sensor_material);

                    //         double z_alternate_module = (iModule%2 == 0) ? 0.0 : stave_dz;
                    //         double x_pos = (r + m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.)*cos(phi) - (-(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2. + iModule*m.sensor_length + iModule*step)*sin(phi);
                    //         double y_pos = (r + m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.)*sin(phi) + (-(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2. + iModule*m.sensor_length + iModule*step)*cos(phi);
                    //         double z_pos = z + z_alternate_petal + z_offset + m.sensor_z_offset + z_alternate_module; 
                    //         if(side == -1){z_pos = -z_pos;}
                    //         Position pos(x_pos, y_pos, z_pos);
                    //         sensor_volume.setVisAttributes(theDetector.visAttributes(m.sensor_viss[i]));
                    //         pv = module_assembly.placeVolume(sensor_volume, Transform3D(rot, pos) );

                    //         sensor_volume.setSensitiveDetector(sens);
                    //         moduleSensThickness[m.name] = m.sensor_thickness; //Assuming one sensitive slice per module
                    //         modules[m.name] = sensor_volume;


                    //         DetElement comp_elt(module, sensor_name, x_det.id());
                    //         comp_elt.setPlacement(pv);
                        iSensor++;
                        mod_num++;
                    }
                    iStave++;
                }

            }
        }
        
    }
    //attach data to detector
    sdet.addExtension< ZDiskPetalsData >( zDiskPetalsData ) ;
    
    std::cout<<"XXX Vertex endcap layers: "<<zDiskPetalsData->layers.size()<<std::endl;
    sdet.setAttributes(theDetector,envelope,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

    return sdet;
}

DECLARE_DETELEMENT(VertexEndcap_o1_v07,create_detector)

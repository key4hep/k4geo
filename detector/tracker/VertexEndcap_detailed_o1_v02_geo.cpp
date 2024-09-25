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
//  F. Bosi as of mid-2023.
//  The staves are arranged in petals, and can feature any number of modules.
//  The modules can be built by smaller rectangular structures to represent
//  both sensitive and insensitive (periphery) parts so that e.g Quad
//  modules can be bult. When using this detector constructor, in addition,
//  the DD4hep_GenericSurfaceInstallerPlugin plugin needs to be instantiated
//  in the xml compact file to define the sensitive surfaces.
//====================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "DD4hep/Printout.h"

using namespace std;

using dd4hep::Assembly;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::getAttrOrDefault;
using dd4hep::Box;
using dd4hep::Tube;
using dd4hep::DEBUG;
using dd4hep::INFO;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens)  {
    xml_det_t   x_det     = e;
    string      det_name  = x_det.nameStr();
    bool        reflect   = x_det.reflect(false);
    DetElement  sdet        (det_name,x_det.id());
    int         m_id=0;
    PlacedVolume pv;
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

    envelope.setVisAttributes(theDetector, x_det.visStr());
    sens.setType("tracker");

    // Struct to support multiple readout or support layers (both using this struct)
    struct componentsStruct{
        string name;
        double z_offset;
        double offset;
        vector<double> thicknesses;
        vector<double> widths;
        vector<double> offsets; 
        vector<double> z_offsets; 
        vector<Material> materials;
        vector<string> viss;
    };

    // Struct to support end-of-stave structures
    struct endOfStaveStruct{
        string name;
        double z_offset;
        double offset;
        vector<double> thicknesses;
        vector<double> lengths;
        vector<double> offsets; 
        vector<double> z_offsets; 
        vector<double> dxs; // Distance of end-of-stave structure to stave itself
        vector<double> xs; // Indicating whether end of stave struct should be place on pos x side (if x>0) or on neg x side (if x<0). Didn't work with using nsides 
        vector<Volume> volumes;
    };

    // --- Module information struct ---
    struct module_information{
        string name;

        vector<componentsStruct> components_vec;
        vector<endOfStaveStruct> endOfStaves;
        double sensor_z_offset;
        double sensor_offset;
        double sensor_thickness;
        vector<bool> sensor_sensitives;
        vector<double> sensor_xmin;
        vector<double> sensor_xmax;
        vector<double> sensor_ymin;
        vector<double> sensor_ymax;
        vector<string> sensor_names;
        double sensor_width;
        double sensor_length;
        vector<Volume> sensor_volumes;
        Material sensor_material;
    };
    list<module_information> module_information_list;

    // --- Collect module(s) information
    for(xml_coll_t mi(x_det,_U(module)); mi; ++mi, ++m_id)  {
        xml_comp_t x_mod   = mi;

        module_information m;
        m.name = x_mod.nameStr();

        // Components
        xml_coll_t c_components(x_mod,_U(components));
        for(c_components.reset(); c_components; ++c_components){
            componentsStruct components;
            components.name = xml_comp_t(c_components).nameStr();
            components.z_offset = xml_comp_t(c_components).z_offset(0);
            components.offset = xml_comp_t(c_components).offset(0);
            xml_coll_t c_component(c_components,_U(component));
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                components.thicknesses.push_back(component.thickness());
                components.widths.push_back(component.width());
                components.offsets.push_back(component.offset(0));
                components.z_offsets.push_back(component.z_offset(0));
                components.materials.push_back(theDetector.material(component.materialStr()));
                components.viss.push_back(component.visStr());
            }
            m.components_vec.push_back(components);
        }

        // End of stave structures
        xml_coll_t c_endOfStave(x_mod,_U(end_z));
        int iEndOfStave = 0;
        for(c_endOfStave.reset(); c_endOfStave; ++c_endOfStave,++iEndOfStave){
            endOfStaveStruct endOfStave;    
            endOfStave.z_offset = xml_comp_t(c_endOfStave).z_offset(0);
            endOfStave.offset = xml_comp_t(c_endOfStave).offset(0);
            endOfStave.name = xml_comp_t(c_endOfStave).nameStr();
            xml_coll_t c_component = xml_coll_t(c_endOfStave,_U(component));
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                endOfStave.thicknesses.push_back(component.thickness());
                endOfStave.lengths.push_back(component.length());
                endOfStave.offsets.push_back(component.offset(0));
                endOfStave.z_offsets.push_back(component.z_offset(0));
                endOfStave.dxs.push_back(component.dx());
                endOfStave.xs.push_back(component.x());

                Box ele_box = Box(component.width()/2., component.length()/2., component.thickness()/2.);
                Volume ele_vol = Volume( endOfStave.name + _toString(iEndOfStave, "_%d"), ele_box, theDetector.material(component.materialStr()));                    
                ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());

                endOfStave.volumes.push_back(ele_vol);
            }
            m.endOfStaves.push_back(endOfStave);
        }

        // Sensor
        xml_coll_t c_sensor(x_mod,_U(sensor));
        m.sensor_z_offset = xml_comp_t(c_sensor).z_offset(0);
        m.sensor_offset = xml_comp_t(c_sensor).offset(0);
        m.sensor_thickness = xml_comp_t(c_sensor).thickness();
        m.sensor_material = theDetector.material(xml_comp_t(c_sensor).materialStr());
        xml_coll_t c_component = xml_coll_t(c_sensor,_U(component));

        int iSensor = 0;
        for(c_component.reset(); c_component; ++c_component){
            xml_comp_t component = c_component;
            m.sensor_sensitives.push_back(component.isSensitive());
            m.sensor_xmin.push_back(component.xmin());
            m.sensor_xmax.push_back(component.xmax());
            m.sensor_ymin.push_back(component.ymin());
            m.sensor_ymax.push_back(component.ymax());
            m.sensor_names.push_back(component.nameStr("sensor"));

            Box ele_box = Box( abs(component.xmax()-component.xmin())/2., abs(component.ymax()-component.ymin())/2., m.sensor_thickness/2.);
            Volume ele_vol = Volume( m.sensor_names.back() + _toString(iSensor, "_%d"), ele_box, m.sensor_material);                    
            ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());

            if(m.sensor_sensitives.back())
                ele_vol.setSensitiveDetector(sens);
            m.sensor_volumes.push_back(ele_vol);
            iSensor++;
        }
        m.sensor_width  = *max_element(m.sensor_xmax.begin(), m.sensor_xmax.end()) - *min_element(m.sensor_xmin.begin(), m.sensor_xmin.end());
        m.sensor_length = *max_element(m.sensor_ymax.begin(), m.sensor_ymax.end()) - *min_element(m.sensor_ymin.begin(), m.sensor_ymin.end());
        printout(DEBUG, det_name, "Module: " + m.name + ", sensor width: " + to_string(m.sensor_width) + ", sensor length: " + to_string(m.sensor_length) );
        module_information_list.push_back(m);
    }
   
    vector<int> sides = {1};
    if(reflect){sides.push_back(-1);}

    printout(INFO, det_name, "Building of detector ...");
    for(auto & side : sides){
        string side_name = det_name + _toString(side,"_side%d");

        Assembly side_assembly(side_name);
        pv = envelope.placeVolume(side_assembly);
        pv.addPhysVolID("side", side);  

        for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
            xml_comp_t  x_layer(li);
            int layer_id        = x_layer.id();
            double rmin         = x_layer.rmin();
            double rmax         = x_layer.rmax();
            double dr           = x_layer.dr(0);
            double z            = x_layer.z();
            double layer_dz     = x_layer.dz(0);
            int nPetals         = x_layer.nphi();
            double phi0_layer   = x_layer.phi0(0);
            double reflect_rot  = x_layer.attr<double>(_Unicode(reflect_rot),0.0);
            
            string disk_name = side_name + _toString(layer_id,"_layer%d");

            // Either use assembly or volume, depending on whether motherVolThickness is given
            PlacedVolume whole_disk_volume_placed;
            Volume whole_disk_volume_v;
            Assembly whole_disk_volume_a;
            double disk_motherVolThickness = getAttrOrDefault(x_layer, _Unicode(motherVolThickness), double(0.0));
            if(disk_motherVolThickness>0.0){
                Tube whole_disk_tube = Tube(rmin, rmax, disk_motherVolThickness/2.);
                whole_disk_volume_v = Volume(disk_name, whole_disk_tube, theDetector.material("Air"));
                whole_disk_volume_v.setVisAttributes(theDetector, x_layer.visStr());
                whole_disk_volume_placed = side_assembly.placeVolume(whole_disk_volume_v, Position(0.0, 0.0, side*(z+disk_motherVolThickness/2.)));
            }
            else{
                //// Use just assembly to avoid overlap between disk volume and other volumes
                whole_disk_volume_a = Assembly(disk_name);
                whole_disk_volume_placed = side_assembly.placeVolume(whole_disk_volume_a, Position(0.0, 0.0, side*z));
            }


            whole_disk_volume_placed.addPhysVolID("layer", layer_id);  

            DetElement diskDE( sdet , disk_name, layer_id);
            diskDE.setPlacement( whole_disk_volume_placed );

            int iModule_tot = 0;
            for(int iPetal=0; iPetal<nPetals; iPetal++){
                double z_alternate_petal = (iPetal%2 == 0) ? 0.0 : layer_dz;

                string petal_name = disk_name + _toString(iPetal,"_petal%d");
                Assembly petal_assembly(petal_name);
                if(disk_motherVolThickness>0.)
                    pv = whole_disk_volume_v.placeVolume(petal_assembly);
                else
                    pv = whole_disk_volume_a.placeVolume(petal_assembly);

                int iStave = 0;
                for(xml_coll_t ri(x_layer,_U(stave)); ri; ++ri,++iStave)  {
                    xml_comp_t x_stave = ri;

                    int    nmodules     = x_stave.nmodules();
                    double r_stave      = x_stave.r();
                    double z_offset     = x_stave.z_offset(0);
                    double stave_dz     = x_stave.dz(0);
                    double step         = x_stave.step();   // Spacing of modules
                    string moduleStr    = x_stave.moduleStr();
                    double phi0_stave   = x_stave.phi0(0);
                    double stave_offset = x_stave.offset(0); // Offset of stave in r-phi
                    double phi     = 2*M_PI/nPetals*iPetal + phi0_layer + phi0_stave + (side == -1 ? reflect_rot : 0.0);

                    // Use the correct module
                    auto m = *find_if(module_information_list.cbegin(), module_information_list.cend(), [&moduleStr] (const module_information& module) {
                        return module.name == moduleStr;
                    });


         
                    // Place all components
                    RotationZYX rot( phi , 0, 0  );
                    double stave_length = nmodules*m.sensor_length + (nmodules-1)*step;
                    double r = rmin + m.sensor_width/2.0 + r_stave + ((iPetal%2 == 0) ? 0.0 : dr);

                    double x_pos = r*cos(phi) - stave_offset*sin(phi);
                    double y_pos = r*sin(phi) + stave_offset*cos(phi);
                    double z_pos = z_alternate_petal + z_offset; 
                    if(side == -1){z_pos = -z_pos;}
                    Position stave_pos = Position(x_pos, y_pos, z_pos);

                    string stave_name = petal_name + _toString(iStave,"_stave%d");

                    PlacedVolume whole_stave_volume_placed;
                    Volume whole_stave_volume_v;
                    Assembly whole_stave_volume_a;
                    double stave_motherVolThickness = getAttrOrDefault(x_stave, _Unicode(motherVolThickness), double(0.0));
                    double stave_motherVolWidth     = getAttrOrDefault(x_stave, _Unicode(motherVolWidth), double(0.0));
                    if(stave_motherVolThickness>0.0 && stave_motherVolWidth>0.0){
                        Box whole_stave_box = Box(stave_motherVolWidth/2., stave_length/2., stave_motherVolThickness/2.);
                        whole_stave_volume_v = Volume(stave_name, whole_stave_box, theDetector.material("Air"));
                        whole_stave_volume_v.setVisAttributes(theDetector, x_stave.visStr());
                        whole_stave_volume_placed = petal_assembly.placeVolume(whole_stave_volume_v, Transform3D(rot,stave_pos));
                    }
                    else{
                        //// Use just assembly to avoid overlap between stave volume and other volumes
                        whole_stave_volume_a = Assembly(stave_name);
                        whole_stave_volume_placed = petal_assembly.placeVolume(whole_stave_volume_a, Transform3D(rot,stave_pos));
                    }

                    // Place components
                    for(auto& component : m.components_vec){
                        Assembly component_assembly(stave_name + "_" + component.name);

                        if(stave_motherVolThickness>0.0 && stave_motherVolWidth>0.0)
                            pv = whole_stave_volume_v.placeVolume(component_assembly);
                        else
                            pv = whole_stave_volume_a.placeVolume(component_assembly);
                        
                        for(int i=0; i<int(component.thicknesses.size()); i++){
                            x_pos = component.offset + component.offsets[i];
                            y_pos = 0.0;
                            z_pos = component.z_offset + component.z_offsets[i] + component.thicknesses[i]/2.; 
                            if(side == -1){z_pos = -z_pos;}
                            Position pos(x_pos, y_pos, z_pos);

                            // Volumes for stave elements cannot be defined for all staves together as they can have different lengths
                            Box ele_box = Box(component.widths[i]/2., stave_length/2., component.thicknesses[i]/2.);
                            Volume ele_vol = Volume(component.name + _toString(i, "_%d"), ele_box, component.materials[i]);                    
                            ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.viss[i]);

                            pv = component_assembly.placeVolume(ele_vol, pos);
                        }
                        component_assembly->GetShape()->ComputeBBox();
                    }

                    // Place end of stave structures
                    int iEndOfStave = 0;
                    for(auto& endOfStave : m.endOfStaves){
                        Assembly endOfStave_assembly(stave_name + "_" + endOfStave.name + _toString(iEndOfStave,"_%d"));

                        if(stave_motherVolThickness>0.0 && stave_motherVolWidth>0.0)
                            pv = whole_stave_volume_v.placeVolume(endOfStave_assembly);
                        else
                            pv = whole_stave_volume_a.placeVolume(endOfStave_assembly);

                        for(int i=0; i<int(endOfStave.thicknesses.size()); i++){
                            x_pos = endOfStave.offset + endOfStave.offsets[i];
                            y_pos = endOfStave.xs[i]>0 ? stave_length/2.+endOfStave.lengths[i]/2.+endOfStave.dxs[i] : -(stave_length/2.+endOfStave.lengths[i]/2.+endOfStave.dxs[i]);
                            z_pos = endOfStave.z_offset + endOfStave.z_offsets[i] + endOfStave.thicknesses[i]/2.; 

                            if(side == -1){z_pos = -z_pos;}
                            Position pos(x_pos, y_pos, z_pos);

                            pv = endOfStave_assembly.placeVolume(endOfStave.volumes[i], pos);
                            iEndOfStave++;
                        }
                        endOfStave_assembly->GetShape()->ComputeBBox();
                    }

                    // Place sensor
                    for(int iModule=0; iModule<nmodules; iModule++){
                        double z_alternate_module = (iModule%2 == 0) ? 0.0 : stave_dz;
                        x_pos = m.sensor_offset;
                        y_pos = - stave_length/2. + m.sensor_length/2. + iModule*m.sensor_length + iModule*step;
                        z_pos = m.sensor_z_offset + z_alternate_module + m.sensor_thickness/2.; 
                        if(side == -1){z_pos = -z_pos;}
                        Position pos(x_pos, y_pos, z_pos);

                        string module_name = stave_name + _toString(iModule,"_module%d");
                        Assembly module_assembly(module_name);

                        if(stave_motherVolThickness>0.0 && stave_motherVolWidth>0.0)
                            pv = whole_stave_volume_v.placeVolume(module_assembly);
                        else
                            pv = whole_stave_volume_a.placeVolume(module_assembly);

                        pv.addPhysVolID("module", iModule_tot);
                        DetElement moduleDE(diskDE,module_name,x_det.id());
                        moduleDE.setPlacement(pv);

                        int iSensitive = 0;
                        for(int i=0; i<int(m.sensor_volumes.size()); i++){
                            x_pos = m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.;
                            y_pos = m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2.;
                            z_pos = 0;
                            Position pos_i(x_pos, y_pos, z_pos);

                            // Place active sensor parts
                            if(m.sensor_sensitives[i]) {
                                string sensor_name = module_name + _toString(iSensitive,"_sensor%d");
                                pv = module_assembly.placeVolume(m.sensor_volumes[i], pos+pos_i);

                                pv.addPhysVolID("sensor", iSensitive);
                                DetElement sensorDE(moduleDE,sensor_name,x_det.id());
                                sensorDE.setPlacement(pv);
                                iSensitive++;;
                            }
                            // Place passive sensor parts
                            else{
                                pv = module_assembly.placeVolume(m.sensor_volumes[i], pos+pos_i);
                            }
                        }
                        iModule_tot++;
                        module_assembly->GetShape()->ComputeBBox();
                    }
                    if(stave_motherVolThickness>0.0 && stave_motherVolWidth>0.0){}
                    else{
                        whole_stave_volume_a->GetShape()->ComputeBBox();
                    }
                }
                petal_assembly->GetShape()->ComputeBBox();
            }
            if(disk_motherVolThickness>0.0){}
            else{
                whole_disk_volume_a->GetShape()->ComputeBBox();
            }
        } 
        side_assembly->GetShape()->ComputeBBox();       
    }
    
    sdet.setAttributes(theDetector,envelope,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
    pv.addPhysVolID("system", x_det.id());

    printout(INFO, det_name, "Building of detector successfully completed.");
 
    return sdet;
}

DECLARE_DETELEMENT(VertexEndcap_detailed_o1_v02,create_detector)
// $Id: $
//====================================================================
//  Based on ZPlanarTracker module from F. Gaede

//  Tracking detector to describe the FCC-ee IDEA vertex detector barrel.
//  The vertex detector is assembled of stave structures which feature 
//  support and readout (flex) elements. Each stave features multiple 
//  individual modules, that consist of sensitive and insensitive
//  sensor elements.
//--------------------------------------------------------------------
//
//  Author     : Armin Ilg
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "XMLHandlerDB.h"
#include <exception>

using namespace std;

using dd4hep::Assembly;
using dd4hep::Box;
using dd4hep::Tube;
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

static Ref_t create_element(Detector& theDetector, xml_h e, SensitiveDetector sens)  {

    xml_det_t    x_det = e;
    int         m_id=0;
    std::string  det_name  = x_det.nameStr();

    DetElement sdet( det_name, x_det.id()  ) ;
    PlacedVolume pv;

    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

    envelope.setVisAttributes(theDetector, x_det.visStr());
    sens.setType("tracker");

    // Struct to support multiple readout or support layers (both using this struct)
    struct componentsStruct{
        string name;
        double r;
        double offset;
        double length;
        vector<double> thicknesses;
        vector<double> offsets; 
        vector<double> rs; 
        vector<Volume> volumes;
    };

    // Struct to support end-of-stave structures
    struct endOfStaveStruct{
        string name;
        double r;
        double offset;
        vector<double> thicknesses;
        vector<double> lengths;
        vector<double> dzs;
        vector<double> offsets; 
        vector<double> rs; 
        vector<Volume> volumes;
    };

    // --- Module information struct ---
    struct stave_information{
        string name;

        vector<componentsStruct> components_vec;
        vector<endOfStaveStruct> endOfStaves;
        double sensor_r;
        double sensor_offset;
        double sensor_thickness;
        Material sensor_material;
        vector<bool> sensor_sensitives;
        vector<double> sensor_xmin;
        vector<double> sensor_xmax;
        vector<double> sensor_ymin;
        vector<double> sensor_ymax;
        double sensor_width;
        double sensor_length;
        double stave_dr;
        vector<Volume> sensor_volumes;
    };
    list<stave_information> stave_information_list;

    // --- Collect stave(s) information
    for(xml_coll_t mi(x_det,_U(stave)); mi; ++mi, ++m_id)  {
        xml_comp_t x_stave   = mi;

        stave_information m;
        m.name = x_stave.nameStr();

        m.stave_dr = x_stave.dr(0); // Offset for every second module in r

        // Components
        xml_coll_t c_components(x_stave,_U(components));
        for(c_components.reset(); c_components; ++c_components){
            componentsStruct components;
            components.name = xml_comp_t(c_components).nameStr();
            components.r = xml_comp_t(c_components).r();
            components.length = xml_comp_t(c_components).length();
            components.offset = xml_comp_t(c_components).offset(0);
            xml_coll_t c_component(c_components,_U(component));
            int iComponent = 0;
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                components.thicknesses.push_back(component.thickness());
                components.offsets.push_back(component.offset(0));
                components.rs.push_back(component.r());

                Box ele_box = Box(component.thickness()/2., component.width()/2., components.length);
                Volume ele_vol = Volume(components.name + _toString(iComponent, "_%d"), ele_box, theDetector.material(component.materialStr()));               
                ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());

                components.volumes.push_back(ele_vol);
            }
            m.components_vec.push_back(components);
        }

        // End of stave structures
        xml_coll_t c_endOfStave(x_stave,_U(end_z));
        for(c_endOfStave.reset(); c_endOfStave; ++c_endOfStave){
            endOfStaveStruct endOfStave;    
            endOfStave.offset = xml_comp_t(c_endOfStave).offset(0);
            endOfStave.r =      xml_comp_t(c_endOfStave).r();
            endOfStave.name =   xml_comp_t(c_endOfStave).nameStr();

            xml_coll_t c_component = xml_coll_t(c_endOfStave,_U(component));
            int iEndOfStave = 0;
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                endOfStave.thicknesses.push_back(component.thickness());
                endOfStave.dzs.push_back(component.dz(0));
                endOfStave.offsets.push_back(component.offset(0));
                endOfStave.lengths.push_back(component.length());
                endOfStave.rs.push_back(component.r());

                Box ele_box = Box(component.thickness()/2., component.width()/2., component.length()/2.);
                Volume ele_vol = Volume(endOfStave.name + _toString(iEndOfStave, "_%d"), ele_box, theDetector.material(component.materialStr()));            
                ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());

                endOfStave.volumes.push_back(ele_vol);
                iEndOfStave++;
            }
            m.endOfStaves.push_back(endOfStave);
        }

        // Sensor
        xml_coll_t c_sensor(x_stave,_U(sensor));
        m.sensor_r = xml_comp_t(c_sensor).r();
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

            // Already create volumes for all sensor components as this is independent of number of sensors per layer
            Box ele_box = Box(m.sensor_thickness/2., abs(component.xmax()-component.xmin())/2., abs(component.ymax()-component.ymin())/2.);            
            Volume ele_vol = Volume("sensor" +  _toString(iSensor, "_%d"), ele_box, m.sensor_material);                    
            ele_vol.setAttributes(theDetector, x_det.regionStr(), x_det.limitsStr(), component.visStr());  

            if(m.sensor_sensitives.back())
                ele_vol.setSensitiveDetector(sens);
            m.sensor_volumes.push_back(ele_vol);
            iSensor++;
        }
        m.sensor_width  = *max_element(m.sensor_xmax.begin(), m.sensor_xmax.end()) - *min_element(m.sensor_xmin.begin(), m.sensor_xmin.end());
        m.sensor_length = *max_element(m.sensor_ymax.begin(), m.sensor_ymax.end()) - *min_element(m.sensor_ymin.begin(), m.sensor_ymin.end());
        cout << "Module: " << m.name << ", sensor width: " << to_string(m.sensor_width)  << ", sensor length: " << to_string(m.sensor_length) << endl;

        stave_information_list.push_back(m);
        cout << "Read stave information of stave " << m.name << endl;
    }

    int iModule_tot = 0;


    //=========  loop over layer elements in xml  ======================================

    cout << "Building vertex barrel detector " << det_name << "..." << endl;
    for(xml_coll_t c(e, _U(layer) ); c; ++c)  {

        xml_comp_t x_layer( c );

        // child elements: ladder, sensitive and periphery
        int layer_id = x_layer.id();
        int nLadders = x_layer.attr<double>(  _Unicode(nLadders) ) ;

        double dr           = x_layer.dr(0);     // Spacing in r for every second stave. 
        double layer_offset = x_layer.offset(0);

        double dphi = 2.*M_PI / double(nLadders);
        double phi0             = x_layer.phi0(0);
        string nameStr          = x_layer.nameStr();
        int    nmodules         = x_layer.nmodules();
        double step             = x_layer.step();   // Spacing of modules


        // Use the correct stave
        auto m = *find_if(stave_information_list.cbegin(), stave_information_list.cend(), [&nameStr] (const stave_information& stave) {
            return stave.name == nameStr;
        });    

        double stave_length = nmodules*m.sensor_length + (nmodules-1)*step;

        double motherVolThickness = getAttrOrDefault(x_layer, _Unicode(motherVolThickness), double(5.0));
        double motherVolLength = getAttrOrDefault(x_layer, _Unicode(motherVolLength), double(stave_length));


        // --- create an assembly and DetElement for the layer
        std::string layer_name = det_name+_toString(layer_id,"_layer%d");
        double motherVolRmin = getAttrOrDefault(x_layer, _Unicode(motherVolRmin), double(x_layer.r()));

        Tube whole_layer_tube = Tube(motherVolRmin, motherVolRmin+motherVolThickness, motherVolLength/2.);
        Volume whole_layer_volume = Volume(layer_name, whole_layer_tube, theDetector.material("Air"));
        whole_layer_volume.setVisAttributes(theDetector, x_det.visStr());
        PlacedVolume whole_layer_placed_volume = envelope.placeVolume(whole_layer_volume);

        whole_layer_placed_volume.addPhysVolID("system", x_det.id()).addPhysVolID("side",0).addPhysVolID("layer", layer_id);  

        DetElement layerDE( sdet , _toString(layer_id,"layer_%d"), layer_id );
        layerDE.setPlacement( whole_layer_placed_volume ) ;

    
        //--------- loop over ladders ---------------------------
        for(int iStave=0; iStave<nLadders; ++iStave) {
            double phi = phi0 + iStave * dphi  ;
            double layer_r = x_layer.r() + ((iStave%2 == 0) ? 0.0 : dr);  // Offset every second stave in r
            RotationZYX rot( phi , 0, 0  ) ;

            string stave_name = layer_name + _toString(iStave,"_stave%d");
            Assembly stave_assembly(stave_name);
            pv = whole_layer_volume.placeVolume(stave_assembly);

            double stave_length = nmodules*m.sensor_length + (nmodules-1)*step;
            
            // Place components
            for(auto& component : m.components_vec){
                Assembly component_assembly(stave_name + "_" + component.name);
                stave_assembly.placeVolume(component_assembly);
                for(int i=0; i<int(component.thicknesses.size()); i++){
                    double r_component = layer_r + component.r + component.rs[i] + component.thicknesses[i]/2.;
                    double r_offset_component = layer_offset + component.offset + component.offsets[i];
                    double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                    double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                    double z_pos = 0.0; 
                    Position pos(x_pos, y_pos, z_pos);

                    component_assembly.placeVolume(component.volumes[i], Transform3D(rot, pos) );
                }
            }

            // Place end of stave structures
            for(auto& endOfStave : m.endOfStaves){
                Assembly endOfStave_assembly(stave_name + "_endOfStave");
                stave_assembly.placeVolume(endOfStave_assembly);
                for(int i=0; i<int(endOfStave.thicknesses.size()); i++){
                    double r_component = layer_r + endOfStave.r + endOfStave.rs[i] + endOfStave.thicknesses[i]/2.;
                    double r_offset_component = layer_offset + endOfStave.offset + endOfStave.offsets[i];
                    double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                    double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                    double z_pos = stave_length/2.+endOfStave.lengths[i]/2.+endOfStave.dzs[i]; 

                    // Place it on both sides of the stave
                    vector<Position> positions = {Position(x_pos, y_pos, z_pos), Position(x_pos, y_pos, -z_pos)};
                    for(auto & pos : positions){                    
                        endOfStave_assembly.placeVolume(endOfStave.volumes[i], Transform3D(rot, pos) );
                    }
                }
            }

            // Place sensor
            for(int iModule=0; iModule<nmodules; iModule++){
                double r_component = layer_r + m.sensor_r + m.sensor_thickness/2. + (iModule%2 == 0 ? 0.0 : m.stave_dr);
                double r_offset_component = layer_offset + m.sensor_offset;
                double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                double z_pos = -(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + iModule*m.sensor_length + iModule*step;
                Position pos(x_pos, y_pos, z_pos);
                    
                string module_name = stave_name + _toString(iModule,"_module%d");
                Assembly module_assembly(module_name);
                pv = stave_assembly.placeVolume(module_assembly);
                pv.addPhysVolID("module", iModule_tot);
                DetElement moduleDE(layerDE,module_name,x_det.id());
                moduleDE.setPlacement(pv);

                // Place all sensor parts
                int iSensitive = 0;
                for(int i=0; i<int(m.sensor_volumes.size()); i++){
                    x_pos = 0.0;
                    y_pos = m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.;
                    z_pos = m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2.;
                    Position pos2(x_pos, y_pos, z_pos);

                    // Place active sensor parts
                    if(m.sensor_sensitives[i]) {
                        string sensor_name = module_name + _toString(iSensitive,"_sensor%d");
                        pv = module_assembly.placeVolume(m.sensor_volumes[i], Transform3D(rot, pos)*Translation3D(pos2));

                        pv.addPhysVolID("sensor", iSensitive);
                        DetElement sensorDE(moduleDE,sensor_name,x_det.id());
                        sensorDE.setPlacement(pv);
                        iSensitive++;
                    }
                    // Place passive sensor parts
                    else{
                        stave_assembly.placeVolume(m.sensor_volumes[i], Transform3D(rot, pos)*Translation3D(pos2));
                    }
                }
                iModule_tot++;
            }
            stave_assembly->GetShape()->ComputeBBox() ;
        }
    }

    sdet.setAttributes(theDetector,envelope,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

    pv.addPhysVolID("system", x_det.id()).addPhysVolID("side",0);

    return sdet;
}

DECLARE_DETELEMENT(VertexBarrel_detailed_o1_v01,create_element)

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
using dd4hep::Translation3D;
using dd4hep::Trapezoid;
using dd4hep::Volume;
using dd4hep::_toString;
//using dd4hep::rec::ZDiskPetalsData;
using dd4hep::Box;

static Ref_t create_detector(Detector& theDetector, xml_h e, SensitiveDetector sens)  {
    xml_det_t   x_det     = e;
    Material    vacuum    = theDetector.vacuum();
    string      det_name  = x_det.nameStr();
    bool        reflect   = x_det.reflect(false);
    DetElement  sdet        (det_name,x_det.id());
    int         m_id=0;
    map<string,Volume> modules;
    PlacedVolume pv;
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------
    
    
    envelope.setVisAttributes(theDetector.invisible()); sens.setType("tracker");

    // -------- reconstruction parameters  ----------------
 //   ZDiskPetalsData*  zDiskPetalsData = new ZDiskPetalsData ;
    std::map< std::string, double > moduleSensThickness;

    // Struct to support multiple readouts
    struct readoutStruct{
        double z_offset;
        double offset;
        vector<double> thicknesses;
        vector<double> widths;
        vector<double> offsets; 
        vector<double> z_offsets; 
        vector<Material> materials;
        vector<string> viss;
    };

    // Struct to support multiple multi-layer supports
    struct supportStruct{
        double z_offset;
        double offset;
        vector<double> thicknesses;
        vector<double> widths;
        vector<double> offsets; 
        vector<double> z_offsets; 
        vector<Material> materials;
        vector<string> viss;
    };

    // --- Module information struct ---
    struct module_information{
        string name;

        vector<readoutStruct> readouts;

        double sensor_z_offset;
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
        vector<string> sensor_viss;
        vector<Box> sensor_boxes;
        vector<PlacedVolume> sensitives;
        vector<Volume> sensitiveMotherVolumes;

        vector<Position> sensitive_pos;
        vector<Volume> sensitive_vol;

        vector<Position> passive_pos;
        vector<Volume> passive_vol;

        vector<supportStruct> supports;
    };
    list<module_information> module_information_list;

    // --- Collect module(s) information
    for(xml_coll_t mi(x_det,_U(module)); mi; ++mi, ++m_id)  {
        xml_comp_t x_mod   = mi;

        module_information m;
        m.name = x_mod.nameStr();

        // Readout
        xml_coll_t c_readout(x_mod,_U(readout));
        for(c_readout.reset(); c_readout; ++c_readout){
            readoutStruct readout;
            readout.z_offset = xml_comp_t(c_readout).z_offset();
            readout.offset = xml_comp_t(c_readout).offset();
            xml_coll_t c_component(c_readout,_U(component));
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                readout.thicknesses.push_back(component.thickness());
                readout.widths.push_back(component.width());
                readout.offsets.push_back(component.offset());
                readout.z_offsets.push_back(component.z_offset());
                readout.materials.push_back(theDetector.material(component.materialStr()));
                readout.viss.push_back(component.visStr());
            }
            m.readouts.push_back(readout);
        }

        // Support
        xml_coll_t c_support(x_mod,_U(support));
        for(c_support.reset(); c_support; ++c_support){
            supportStruct support;    
            support.z_offset = xml_comp_t(c_support).z_offset();
            support.offset = xml_comp_t(c_support).offset();
            xml_coll_t c_component = xml_coll_t(c_support,_U(component));
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                support.thicknesses.push_back(component.thickness());
                support.widths.push_back(component.width());
                support.offsets.push_back(component.offset());
                support.z_offsets.push_back(component.z_offset());
                support.materials.push_back(theDetector.material(component.materialStr()));
                support.viss.push_back(component.visStr());
            }
            m.supports.push_back(support);
        }

        // Sensor
        xml_coll_t c_sensor(x_mod,_U(sensor));
        m.sensor_z_offset = xml_comp_t(c_sensor).z_offset();
        m.sensor_offset = xml_comp_t(c_sensor).offset();
        m.sensor_thickness = xml_comp_t(c_sensor).thickness();
        m.sensor_material = theDetector.material(xml_comp_t(c_sensor).materialStr());
        xml_coll_t c_component = xml_coll_t(c_sensor,_U(component));

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
            m.sensor_boxes.push_back(ele_box);
            iComponent++;
        }
        m.sensor_width  = *max_element(m.sensor_xmax.begin(), m.sensor_xmax.end()) - *min_element(m.sensor_xmin.begin(), m.sensor_xmin.end());
        m.sensor_length = *max_element(m.sensor_ymax.begin(), m.sensor_ymax.end()) - *min_element(m.sensor_ymin.begin(), m.sensor_ymin.end());
        cout << "Module: " << m.name << ", sensor width: " << to_string(m.sensor_width)  << ", sensor length: " << to_string(m.sensor_length) << endl;

        int iSensitive = 0, iPassive = 0;
        for(int i=0; i<int(m.sensor_boxes.size()); i++){
            double x_pos = m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.;
            double y_pos = m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2.;
            double z_pos = 0;
            Position pos(x_pos, y_pos, z_pos);

            if(m.sensor_sensitives[i]) {
                Volume sensitive_vol(m.name + _toString(iSensitive, "sensor_%d"), m.sensor_boxes[i], m.sensor_material);
                sensitive_vol.setVisAttributes(theDetector.visAttributes(m.sensor_viss[i]));
                sensitive_vol.setSensitiveDetector(sens); 

                m.sensitive_vol.push_back(sensitive_vol);
                m.sensitive_pos.push_back(pos);
                iSensitive++;
            }
            else{
                Volume passive_vol = Volume(m.name + _toString(iPassive, "_passive%d"), m.sensor_boxes[i], m.sensor_material);                    
                passive_vol.setVisAttributes(theDetector.visAttributes(m.sensor_viss[i]));
                m.passive_vol.push_back(passive_vol);
                m.passive_pos.push_back(pos);
                iPassive++;
            }
        }
        moduleSensThickness[m.name] = m.sensor_thickness;

        module_information_list.push_back(m);
    }
   
    vector<int> sides = {1};
    if(reflect){sides.push_back(-1);}
 
    for(auto & side : sides){
        string side_name = det_name + _toString(side,"_side%d");
        Assembly side_assembly(side_name);
        DetElement sideDE(sdet, side_name, x_det.id());
        pv = envelope.placeVolume(side_assembly);
        pv.addPhysVolID("side", side);
        sideDE.setPlacement(pv);


        for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
            xml_comp_t  x_layer(li);
            int layer_id            = x_layer.id();
            double rmin         = x_layer.rmin();
            double dr           = x_layer.dr();
            double z            = x_layer.z();
            double layer_dz     = x_layer.dz();
            int nPetals         = x_layer.nphi();
            double phi0_layer   = x_layer.phi0();
            int mod_num = 0;

                       
            // -------- reconstruction parameters  ----------------
            //NOTE: Mostly Dummy information for event display/DDMarlinPandora
            //FIXME: have to see how to properly accommodate spirals and the fact that there's only
            //one ring
            //ZDiskPetalsData::LayerLayout thisLayer ;

            int numberOfRings=0; //check that only one ring is used
            
            string layer_name = side_name + _toString(layer_id,"_layer%d");
            Assembly layer_assembly(layer_name);
            DetElement layerDE( sideDE , layer_name, layer_id );
            pv = side_assembly.placeVolume(layer_assembly);
            pv.addPhysVolID("layer", layer_id );  
            layerDE.setPlacement( pv ) ;

            for(int iPetal=0; iPetal<nPetals; iPetal++){
                double z_alternate_petal = (iPetal%2 == 0) ? 0.0 : layer_dz;

                string petal_name = layer_name + _toString(iPetal,"_petal%d");
                // Assembly petal_assembly(petal_name);
                // DetElement petalDE( layerDE , petal_name, iPetal );
                // pv = layer_assembly.placeVolume(petal_assembly);
                // petalDE.setPlacement(pv);

                int iStave = 0;
                int nStaves = 0;
                for(xml_coll_t ri(x_layer,_U(stave)); ri; ++ri)  
                    nStaves+=1;
        
                for(xml_coll_t ri(x_layer,_U(stave)); ri; ++ri)  {
                    xml_comp_t x_stave = ri;

                    int    nmodules     = x_stave.nmodules();
                    double r_stave      = x_stave.r();
                    double z_offset     = x_stave.z_offset();
                    double stave_dz     = x_stave.dz(0);
                    double step         = x_stave.step();   // Spacing of modules
                    string moduleStr    = x_stave.moduleStr();
                    double phi0_stave   = x_stave.phi0();
                    double phi     = 2*M_PI/nPetals*iPetal + phi0_layer + phi0_stave;

                    // Use the correct module
                    auto m = *find_if(module_information_list.cbegin(), module_information_list.cend(), [&moduleStr] (const module_information& module) {
                        return module.name == moduleStr;
                    });
            
                    string stave_name = petal_name + _toString(iStave,"_stave%d");
                    // Assembly stave_assembly(stave_name);
                    // DetElement staveDE()
                    // pv = layer_assembly.placeVolume(stave_assembly);
                    
                    // Place all components
                    RotationZYX rot( phi , 0, 0  );

                    double stave_length = nmodules*m.sensor_length + (nmodules-1)*step;
                    double r = rmin + m.sensor_width/2.0 + r_stave + ((iPetal%2 == 0) ? 0.0 : dr);

                    
                    // Place support
                    int iSupport = 0;
                    for(auto& support : m.supports){
                        for(int i=0; i<int(support.thicknesses.size()); i++){
                            double r_component = r + support.offset + support.offsets[i];
                            double r_offset_component = 0.0;
                            double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                            double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                            double z_pos = z + z_alternate_petal + z_offset + support.z_offset + support.z_offsets[i] + support.thicknesses[i]/2.; 
                            if(side == -1){z_pos = -z_pos;}
                            Position pos(x_pos, y_pos, z_pos);

                            Box ele_box = Box( support.widths[i]/2., stave_length/2., support.thicknesses[i]/2.);
                            Volume ele_vol = Volume(_toString(int(support.thicknesses.size())*iSupport + i, "suport_%d"), ele_box, support.materials[i]);                    
                            ele_vol.setVisAttributes(theDetector.visAttributes(support.viss[i]));

                            pv = envelope.placeVolume(ele_vol, Transform3D(rot, pos) );
                        }
                        iSupport++;
                    }

                    // Place readout
                    int iReadout = 0;
                    for(auto& readout : m.readouts){
                        for(int i=0; i<int(readout.thicknesses.size()); i++){
                            double r_component = r + readout.offset + readout.offsets[i];
                            double r_offset_component = 0.0;
                            double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                            double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                            double z_pos = z + z_alternate_petal + z_offset + readout.z_offset + readout.z_offsets[i] + readout.thicknesses[i]/2.;
                            if(side == -1){z_pos = -z_pos;}
                            Position pos(x_pos, y_pos, z_pos);

                            Box ele_box = Box( readout.widths[i]/2., stave_length/2., readout.thicknesses[i]/2.);
                            Volume ele_vol = Volume( _toString(int(readout.thicknesses.size())*iReadout + i, "readout_%d"), ele_box, readout.materials[i]);                    
                            ele_vol.setVisAttributes(theDetector.visAttributes(readout.viss[i]));
                            pv = envelope.placeVolume(ele_vol, Transform3D(rot, pos) );
                        }
                        iReadout++;
                    }

                    // Place sensor
                    for(int iModule=0; iModule<nmodules; iModule++){
                        double z_alternate_module = (iModule%2 == 0) ? 0.0 : stave_dz;
                        double r_component = r + m.sensor_offset;
                        double r_offset_component = -(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + iModule*m.sensor_length + iModule*step;
                        double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                        double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                        double z_pos = z + z_alternate_petal + z_offset + m.sensor_z_offset + z_alternate_module + m.sensor_thickness/2.; 
                        if(side == -1){z_pos = -z_pos;}
                        Position pos(x_pos, y_pos, z_pos);

                        int iSensor=0;
                        // string module_name = stave_name + _toString(iModule,"_module%d");
                        string module_name = stave_name + _toString(mod_num,"_module%d");

                        // Place active sensor parts

                        for(int i=0; i<int(m.sensitive_vol.size()); i++)  {
                            string sensor_name = module_name + _toString(i,"_sensorMotherVolume%d");
                            
                            DetElement module(sdet,sensor_name,x_det.id());
                            pv = envelope.placeVolume(m.sensitive_vol[i], Transform3D(rot, pos)*Translation3D(m.sensitive_pos[i]));
                            m.sensitive_vol[i].setSensitiveDetector(sens);


                            pv.addPhysVolID("side", side ).addPhysVolID("layer", layer_id ).addPhysVolID("module",mod_num).addPhysVolID("sensor", iSensor);
                            // cout << "side: " << _toString(side) << "layer: " << _toString(layer_id) << "module: " << _toString(mod_num) << "sensor: " << _toString(iSensor) << endl; 
                            module.setPlacement(pv);
                            iSensor++;
                        }         
                        mod_num++;

                        // Place passive sensor parts
                        for(int i=0; i<int(m.passive_vol.size()); i++)  {
                            pv = envelope.placeVolume(m.passive_vol[i], Transform3D(rot, pos)*Translation3D(m.passive_pos[i]));
                        }
                    }
                    iStave++;
                }

            }
            layer_assembly->GetShape()->ComputeBBox() ;

        }
        side_assembly->GetShape()->ComputeBBox() ;
        
    }
    //attach data to detector
    // sdet.addExtension< ZDiskPetalsData >( zDiskPetalsData ) ;
    
    cout<<"Built vertex disks detector: " << std::endl;
    sdet.setAttributes(theDetector,envelope,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

    return sdet;
}

DECLARE_DETELEMENT(VertexDisks_detailed_o1_v01,create_detector)

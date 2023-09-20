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

#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include <exception>

#include <UTIL/BitField64.h>
#include <UTIL/BitSet32.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

using namespace std;

using dd4hep::Assembly;
using dd4hep::Box;
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
using dd4hep::Trapezoid;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::rec::Vector3D;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::VolPlane;
using dd4hep::rec::volSurfaceList;

static Ref_t create_element(Detector& theDetector, xml_h e, SensitiveDetector sens)  {

    xml_det_t    x_det = e;
    int         m_id=0;
    std::string  det_name  = x_det.nameStr();

    DetElement sdet( det_name, x_det.id()  ) ;
    PlacedVolume pv;

    // put the whole detector into an assembly
    //  - should be replaced by an envelope volume ...

    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;
    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

    envelope.setVisAttributes(theDetector, x_det.visStr());
    sens.setType("tracker");

    // Struct to support multiple readouts
    struct readoutStruct{
        double dr;
        double offset;
        vector<double> thicknesses;
        vector<double> widths;
        vector<double> offsets; 
        vector<double> drs; 
        vector<Material> materials;
        vector<string> viss;
    };

    // Struct to support multiple multi-layer supports
    struct supportStruct{
        double dr;
        double offset;
        vector<double> thicknesses;
        vector<double> widths;
        vector<double> offsets; 
        vector<double> drs; 
        vector<Material> materials;
        vector<string> viss;
    };

    // --- Module information struct ---
    struct stave_information{
        string name;

        vector<readoutStruct> readouts;

        double sensor_dr;
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
        vector<supportStruct> supports;
    };
    list<stave_information> stave_information_list;

    // --- Collect stave(s) information
    for(xml_coll_t mi(x_det,_U(stave)); mi; ++mi, ++m_id)  {
        xml_comp_t x_stave   = mi;

        stave_information m;
        m.name = x_stave.nameStr();

        // Readout
        xml_coll_t c_readout(x_stave,_U(readout));
        for(c_readout.reset(); c_readout; ++c_readout){
            readoutStruct readout;
            readout.dr = xml_comp_t(c_readout).dr();
            readout.offset = xml_comp_t(c_readout).offset();
            xml_coll_t c_component(c_readout,_U(component));
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                readout.thicknesses.push_back(component.thickness());
                readout.widths.push_back(component.width());
                readout.offsets.push_back(component.offset());
                readout.drs.push_back(component.dr());
                readout.materials.push_back(theDetector.material(component.materialStr()));
                readout.viss.push_back(component.visStr());
            }
            m.readouts.push_back(readout);
        }

        // Support
        xml_coll_t c_support(x_stave,_U(support));
        for(c_support.reset(); c_support; ++c_support){
            supportStruct support;    
            support.dr = xml_comp_t(c_support).dr();
            support.offset = xml_comp_t(c_support).offset();
            xml_coll_t c_component = xml_coll_t(c_support,_U(component));
            for(c_component.reset(); c_component; ++c_component){
                xml_comp_t component = c_component;
                support.thicknesses.push_back(component.thickness());
                support.widths.push_back(component.width());
                support.offsets.push_back(component.offset());
                support.drs.push_back(component.dr());
                support.materials.push_back(theDetector.material(component.materialStr()));
                support.viss.push_back(component.visStr());
            }
            m.supports.push_back(support);
        }

        // Sensor
        xml_coll_t c_sensor(x_stave,_U(sensor));
        m.sensor_dr = xml_comp_t(c_sensor).dr();
        m.sensor_offset = xml_comp_t(c_sensor).offset();
        m.sensor_thickness = xml_comp_t(c_sensor).thickness();
        m.sensor_material = theDetector.material(xml_comp_t(c_sensor).materialStr());
        xml_coll_t c_component = xml_coll_t(c_sensor,_U(component));

        for(c_component.reset(); c_component; ++c_component){
            xml_comp_t component = c_component;
            m.sensor_sensitives.push_back(component.isSensitive());
            m.sensor_xmin.push_back(component.xmin());
            m.sensor_xmax.push_back(component.xmax());
            m.sensor_ymin.push_back(component.ymin());
            m.sensor_ymax.push_back(component.ymax());
            m.sensor_viss.push_back(component.visStr());

            // Already create volumes for all sensor components as this is independent of number of sensors per layer
            Box ele_box = Box(m.sensor_thickness/2., abs(component.xmax()-component.xmin())/2., abs(component.ymax()-component.ymin())/2.);
            m.sensor_boxes.push_back(ele_box);
        }
        m.sensor_width  = *max_element(m.sensor_xmax.begin(), m.sensor_xmax.end()) - *min_element(m.sensor_xmin.begin(), m.sensor_xmin.end());
        m.sensor_length = *max_element(m.sensor_ymax.begin(), m.sensor_ymax.end()) - *min_element(m.sensor_ymin.begin(), m.sensor_ymin.end());
        cout << "Module: " << m.name << ", sensor width: " << to_string(m.sensor_width)  << ", sensor length: " << to_string(m.sensor_length) << endl;

        stave_information_list.push_back(m);
        cout << "Read stave information of stave " << m.name << endl;
    }



    //=========  loop over layer elements in xml  ======================================

    cout << "Building vertex barrel detector " << det_name << "..." << endl;
    for(xml_coll_t c(e, _U(layer) ); c; ++c)  {

        xml_comp_t x_layer( c );

        // child elements: ladder, sensitive and periphery
        int layer_id = x_layer.id();
        int nLadders = x_layer.attr<double>(  _Unicode(nLadders) ) ;
        double layer_r = x_layer.r();
        double layer_offset = x_layer.offset();

        double dphi = 2.*M_PI / double(nLadders);
        double phi0             = x_layer.phi0();
        string nameStr        = x_layer.nameStr();
        int    nmodules         = x_layer.nmodules();
        double step             = x_layer.step();   // Spacing of modules


        // Use the correct stave
        auto m = *find_if(stave_information_list.cbegin(), stave_information_list.cend(), [&nameStr] (const stave_information& stave) {
            return stave.name == nameStr;
        });    

        // --- create an assembly and DetElement for the layer

        std::string layer_name = det_name+_toString(layer_id,"_layer%d");
        Assembly layer_assembly(layer_name);
        DetElement layerDE( sdet , _toString(layer_id,"layer_%d"), layer_id );
        pv = envelope.placeVolume(  layer_assembly );
        pv.addPhysVolID("layer", layer_id);  
        layerDE.setPlacement( pv ) ;


        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        Vector3D u( 0. , 1. , 0. ) ;
        Vector3D v( 0. , 0. , 1. ) ;
        Vector3D n( 1. , 0. , 0. ) ;
        // Vector3D o( 0. , 0. , 0. ) ;

        //--------- loop over ladders ---------------------------
        for(int iStave=0; iStave<nLadders; ++iStave) {
            double phi = phi0 + iStave * dphi  ;

            RotationZYX rot( phi , 0, 0  ) ;

            string stave_name = layer_name + _toString(iStave,"_stave%d");
            Assembly stave_assembly(stave_name);
            DetElement staveDE(layerDE, stave_name, iStave);
            pv = layer_assembly.placeVolume(stave_assembly);
            pv.addPhysVolID("stave", iStave);
            staveDE.setPlacement(pv);

            string support_name = stave_name + "_support";
            Assembly support_assembly(support_name);
            // DetElement supportDE(staveDE, support_name, 0);
            pv = stave_assembly.placeVolume(support_assembly);
            // supportDE.setPlacement(pv);

            string readout_name = stave_name + "_readout";
            Assembly readout_assembly(readout_name);
            // DetElement readoutDE(staveDE, readout_name, 0);
            pv = stave_assembly.placeVolume(readout_assembly);
            // readoutDE.setPlacement(pv);

            // Place all components
            double stave_length = nmodules*m.sensor_length + (nmodules-1)*step;
            
            // Place support
            int iSupport = 0;
            for(auto& support : m.supports){
                for(int i=0; i<int(support.thicknesses.size()); i++){
                    double r_component = layer_r + support.dr + support.drs[i] + support.thicknesses[i]/2.;
                    double r_offset_component = layer_offset + support.offset + support.offsets[i];
                    double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                    double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                    double z_pos = 0.0; 
                    Position pos(x_pos, y_pos, z_pos);

                    Box ele_box = Box(support.thicknesses[i]/2., support.widths[i]/2., stave_length/2.);
                    Volume ele_vol = Volume(_toString(int(support.thicknesses.size())*iSupport + i, "support_%d"), ele_box, support.materials[i]);                    
                    ele_vol.setVisAttributes(theDetector.visAttributes(support.viss[i]));

                    pv = support_assembly.placeVolume(ele_vol, Transform3D(rot, pos) );
                }
                iSupport++;
            }

            // Place readout
            int iReadout = 0;
            for(auto& readout : m.readouts){
                for(int i=0; i<int(readout.thicknesses.size()); i++){
                    double r_component = layer_r + readout.dr + readout.drs[i] + readout.thicknesses[i]/2.;
                    double r_offset_component = layer_offset + readout.offset + readout.offsets[i];
                    double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                    double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                    double z_pos = 0.0; 
                    Position pos(x_pos, y_pos, z_pos);

                    Box ele_box = Box(readout.thicknesses[i]/2., readout.widths[i]/2., stave_length/2.);
                    Volume ele_vol = Volume( _toString(int(readout.thicknesses.size())*iReadout + i, "readout_%d"), ele_box, readout.materials[i]);                    
                    ele_vol.setVisAttributes(theDetector.visAttributes(readout.viss[i]));
                    pv = readout_assembly.placeVolume(ele_vol, Transform3D(rot, pos) );
                }
                iReadout++;
            }

            // Place sensor
            for(int iModule=0; iModule<nmodules; iModule++){
                double r_component = layer_r + m.sensor_dr + m.sensor_thickness/2.;
                double r_offset_component = layer_offset + m.sensor_offset;
                double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                double z_pos = -(nmodules-1)/2.*(m.sensor_length) - (nmodules-1)/2.*step + iModule*m.sensor_length + iModule*step;
                Position pos(x_pos, y_pos, z_pos);
                    
                string module_name = stave_name + _toString(iModule,"_module%d");
                Assembly module_assembly(module_name);
                DetElement moduleDE(staveDE, module_name, iModule);
                pv = stave_assembly.placeVolume(module_assembly);
                pv.addPhysVolID("module", iModule);
                moduleDE.setPlacement(pv);

                string sensors_name = module_name + "_sensors";
                // Assembly sensors_assembly(sensors_name);
                // DetElement sensorsDE(moduleDE, sensors_name, x_det.id());
                // pv = module_assembly.placeVolume(sensors_assembly);
                // sensorsDE.setPlacement(pv);

                // string passives_name = module_name + "_passives";
                // Assembly passives_assembly(passives_name);
                // DetElement passivesDE(moduleDE, passives_name, x_det.id());
                // pv = module_assembly.placeVolume(passives_assembly);
                // passivesDE.setPlacement(pv);

                // Place all sensor parts
                int iSensitive = 0, iPassive = 0;
                for(int i=0; i<int(m.sensor_boxes.size()); i++){
                    x_pos = 0.0;
                    y_pos = m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.;
                    z_pos = m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2.;
                    Position pos2(x_pos, y_pos, z_pos);

                    // Place active sensor parts
                    if(m.sensor_sensitives[i]) {
                        string sensor_name = sensors_name + _toString(iSensitive,"_sensor%d");
                        Volume sensitive_vol(sensor_name, m.sensor_boxes[i], m.sensor_material);
                        sensitive_vol.setVisAttributes(theDetector.visAttributes(m.sensor_viss[i]));
                        sensitive_vol.setSensitiveDetector(sens); 

                        double inner_thickness = -m.sensor_thickness/2.;
                        double outer_thickness = m.sensor_thickness/2.;
                        // VolPlane surf( sensitive_vol , SurfaceType(SurfaceType::Sensitive) , inner_thickness, outer_thickness, u,v,n);

                        DetElement sensorDE(sdet,sensor_name,x_det.id());
                        pv = module_assembly.placeVolume(sensitive_vol, Transform3D(rot, pos)*Translation3D(pos2));
                        sensitive_vol.setSensitiveDetector(sens);

                        pv.addPhysVolID("sensor", iSensitive);
                        sensorDE.setPlacement(pv);

                        // volSurfaceList( sensorDE )->push_back( surf ) ;
                        iSensitive++;
                    }
                    // Place passive sensor parts
                    else{
                        string passive_name = stave_name + _toString(iPassive,"_passive%d");
                        Volume passive_vol = Volume(passive_name, m.sensor_boxes[i], m.sensor_material);                    
                        passive_vol.setVisAttributes(theDetector.visAttributes(m.sensor_viss[i]));
                        pv = stave_assembly.placeVolume(passive_vol, Transform3D(rot, pos)*Translation3D(pos2));
                        iPassive++;
                    }
                }
            }
            stave_assembly->GetShape()->ComputeBBox() ;
        }
        layer_assembly->GetShape()->ComputeBBox() ;
    }

    sdet.setAttributes(theDetector,envelope,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

    // pv.addPhysVolID( "system", x_det.id() );
    // sdet.setPlacement(pv);

    return sdet;
}

DECLARE_DETELEMENT(VertexBarrel_detailed_o1_v01,create_element)

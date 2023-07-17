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
//using dd4hep::rec::ZPlanarData;
using dd4hep::rec::NeighbourSurfacesData;
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
  envelope.setVisAttributes(theDetector.invisible()); sens.setType("tracker");

  // for encoding
  std::string cellIDEncoding = sens.readout().idSpec().fieldDescription();
  UTIL::BitField64 encoder( cellIDEncoding );
  encoder.reset();
  encoder[lcio::LCTrackerCellID::subdet()] = x_det.id();
  encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::barrel;

  // ZPlanarData*  zPlanarData = new ZPlanarData ;
  NeighbourSurfacesData*  neighbourSurfacesData = new NeighbourSurfacesData() ;

  dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

  double minRadius = 1e99 ;
  double minZhalf = 1e99 ;

  std::map< std::string, double > moduleSensThickness;

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
      vector<PlacedVolume> sensitives;
      vector<Volume> sensitiveMotherVolumes;

      vector<Position> sensitive_pos;
      vector<Volume> sensitive_vol;

      vector<Position> passive_pos;
      vector<Volume> passive_vol;

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
          Box ele_box = Box(m.sensor_thickness/2., abs(component.xmax()-component.xmin())/2., abs(component.ymax()-component.ymin())/2.);
          m.sensor_boxes.push_back(ele_box);
          iComponent++;
      }
      m.sensor_width  = *max_element(m.sensor_xmax.begin(), m.sensor_xmax.end()) - *min_element(m.sensor_xmin.begin(), m.sensor_xmin.end());
      m.sensor_length = *max_element(m.sensor_ymax.begin(), m.sensor_ymax.end()) - *min_element(m.sensor_ymin.begin(), m.sensor_ymin.end());
      cout << "Module: " << m.name << ", sensor width: " << to_string(m.sensor_width)  << ", sensor length: " << to_string(m.sensor_length) << endl;

      int iSensitive = 0, iPassive = 0;
      for(int i=0; i<int(m.sensor_boxes.size()); i++){
          double x_pos = 0.0;
          double y_pos = m.sensor_xmin[i]+abs(m.sensor_xmax[i]-m.sensor_xmin[i])/2.;
          double z_pos = m.sensor_ymin[i]+abs(m.sensor_ymax[i]-m.sensor_ymin[i])/2.;
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
    double phi_tilt         = x_layer.phi_tilt();
    string nameStr        = x_layer.nameStr();
    int    nmodules         = x_layer.nmodules();
    double step             = x_layer.step();   // Spacing of modules

    std::string layer_name = det_name+_toString(layer_id,"_layer%d");
    int mod_num = 0;
    
    // Use the correct stave
    auto m = *find_if(stave_information_list.cbegin(), stave_information_list.cend(), [&nameStr] (const stave_information& stave) {
        return stave.name == nameStr;
    });    

    // --- create an assembly and DetElement for the layer

    Assembly layer_assembly( det_name + "_layer_" +_toString(layer_id,"_%d") );
    DetElement layerDE( sdet , _toString(layer_id,"layer_%d"), x_det.id() );
    pv = envelope.placeVolume(  layer_assembly );
    pv.addPhysVolID("layer", layer_id );  
    layerDE.setPlacement( pv ) ;

   
    // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
    Vector3D u( 0. , 1. , 0. ) ;
    Vector3D v( 0. , 0. , 1. ) ;
    Vector3D n( 1. , 0. , 0. ) ;

    SurfaceType type( SurfaceType::Sensitive ) ;
    sens.setType("tracker");

    //--------- loop over ladders ---------------------------
    for(int iStave=0; iStave<nLadders; ++iStave) {
        double phi = phi0 + iStave * dphi  ;

        RotationZYX rot( phi , 0, 0  ) ;
        // RotationZYX rot2( phi_tilt , 0, 0  ) ; // Tilt not working yet

        string stave_name = layer_name + _toString(iStave,"_stave%d");
        // Assembly stave_assembly(stave_name);
        // DetElement staveDE()
        // pv = layer_assembly.placeVolume(stave_assembly);
        
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
                double r_component = layer_r + readout.dr + readout.drs[i] + readout.thicknesses[i]/2.;
                double r_offset_component = layer_offset + readout.offset + readout.offsets[i];
                double x_pos = r_component*cos(phi) - r_offset_component*sin(phi);
                double y_pos = r_component*sin(phi) + r_offset_component*cos(phi);
                double z_pos = 0.0; 
                Position pos(x_pos, y_pos, z_pos);

                Box ele_box = Box(readout.thicknesses[i]/2., readout.widths[i]/2., stave_length/2.);
                Volume ele_vol = Volume( _toString(int(readout.thicknesses.size())*iReadout + i, "readout_%d"), ele_box, readout.materials[i]);                    
                ele_vol.setVisAttributes(theDetector.visAttributes(readout.viss[i]));
                pv = envelope.placeVolume(ele_vol, Transform3D(rot, pos) );
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
                
            int iSensor=0;
            string module_name = stave_name + _toString(mod_num,"_module%d");

            // Place active sensor parts

            for(int i=0; i<int(m.sensitive_vol.size()); i++)  {
                string sensor_name = module_name + _toString(i,"_sensorMotherVolume%d");
                
                DetElement module(sdet,sensor_name,x_det.id());
                pv = envelope.placeVolume(m.sensitive_vol[i], Transform3D(rot, pos)*Translation3D(m.sensitive_pos[i]));
                m.sensitive_vol[i].setSensitiveDetector(sens);


                pv.addPhysVolID("layer", layer_id ).addPhysVolID("module",mod_num).addPhysVolID("sensor", iSensor);
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

      // DetElement ladderDE( layerDE ,  laddername , x_det.id() );
      // ladderDE.setPlacement( pv ) ;

      // volSurfaceList( ladderDE )->push_back( surf ) ;


      ///////////////////
      encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::barrel;
      encoder[lcio::LCTrackerCellID::layer()] = layer_id;
      encoder[lcio::LCTrackerCellID::module()] = nLadders;
      encoder[lcio::LCTrackerCellID::sensor()] = 0; // there is no sensor defintion in VertexBarrel at the moment

      dd4hep::CellID cellID = encoder.lowWord(); // 32 bits

      //compute neighbours 

      int n_neighbours_module = 1; // 1 gives the adjacent modules (i do not think we would like to change this)

      int newmodule=0;

      for(int imodule=-n_neighbours_module; imodule<=n_neighbours_module; imodule++){ // neighbouring modules
		    
        if (imodule==0) continue; // cellID we started with
        newmodule = nLadders + imodule;
              
        //compute special case at the boundary  
        //general computation to allow (if necessary) more then adiacent neighbours (ie: +-2)
        if (newmodule < 0) newmodule = nLadders + newmodule;
        if (newmodule >= nLadders) newmodule = newmodule - nLadders;

        //encoding
        encoder[lcio::LCTrackerCellID::module()] = newmodule;
        encoder[lcio::LCTrackerCellID::sensor()] = 0;

        neighbourSurfacesData->sameLayer[cellID].push_back(encoder.lowWord());
 
      }
    }
    

    // is this needed ??
    layer_assembly->GetShape()->ComputeBBox() ;


    layer_assembly->GetShape()->ComputeBBox() ;
  }

#if 0  //-------- add an inscribing cylinder of air for tracking purposes -----------------
       //  this screw up the geometry and the material scan does not work anymore !!!!!?????

  double tube_thick =  1.0 * dd4hep::mm ;
  double inner_r    =  minRadius - 1.1 * tube_thick ;
  double outer_r    =  inner_r + tube_thick ;
  double z_half     =  minZhalf ; 
  
  Tube   tubeSolid (inner_r, outer_r, z_half ) ;
  Volume tube_vol( det_name+"_inner_cylinder_air", tubeSolid ,  theDetector.material("Air") ) ;
  
  envelope.placeVolume( tube_vol , Transform3D() ) ;
  
  Vector3D ocyl(  inner_r + 0.5*tube_thick , 0. , 0. ) ;
  
  VolCylinder cylSurf( tube_vol , SurfaceType( SurfaceType::Helper ) , 0.5*tube_thick  , 0.5*tube_thick , ocyl ) ;
  
  volSurfaceList( sdet )->push_back( cylSurf ) ;
  
#endif //----------------------------------------------------------------------------------



  // tracker.addExtension< ZPlanarData >( zPlanarData ) ;
  sdet.addExtension< NeighbourSurfacesData >( neighbourSurfacesData ) ;

  //set the vis Attribute (added by Thorben Quast)
  envelope.setVisAttributes(theDetector, x_det.visStr());
  
  pv.addPhysVolID( "system", x_det.id() ).addPhysVolID("side",0 )  ;
  
  sdet.setPlacement(pv);
       

  return sdet;
}

DECLARE_DETELEMENT(VertexBarrel_detailed_o1_v01,create_element)

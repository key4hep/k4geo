// $Id: $
//====================================================================
//  Based on ZPlanarTracker module from F. Gaede

//  Tracking detector to describe the FCC-ee IDEA vertex detector barrel.
//  The vertex detector is assembled of stave structures which feature 
//  support and readout (flex) elements. Each stave features multiple 
//  individual sensors, where each sensor features a sensitive and periphery
//  region.
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
using dd4hep::Trapezoid;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::rec::ZPlanarData;
using dd4hep::rec::NeighbourSurfacesData;
using dd4hep::rec::Vector3D;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::VolPlane;
using dd4hep::rec::volSurfaceList;

static Ref_t create_element(Detector& theDetector, xml_h e, SensitiveDetector sens)  {
  
  xml_det_t    x_det = e;
  std::string  det_name  = x_det.nameStr();

  // put the whole detector into an assembly
  //  - should be replaced by an envelope volume ...
  
  Assembly assembly( det_name+"_assembly" );
  DetElement tracker( det_name, x_det.id()  ) ;
  PlacedVolume pv;


  // for encoding
  std::string cellIDEncoding = sens.readout().idSpec().fieldDescription();
  UTIL::BitField64 encoder( cellIDEncoding );
  encoder.reset();
  encoder[lcio::LCTrackerCellID::subdet()] = x_det.id();
  encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::barrel;



  ZPlanarData*  zPlanarData = new ZPlanarData ;
  NeighbourSurfacesData*  neighbourSurfacesData = new NeighbourSurfacesData() ;

  dd4hep::xml::setDetectorTypeFlag( e, tracker ) ;

  double minRadius = 1e99 ;
  double minZhalf = 1e99 ;


  bool isStripDetector = false ;
  try {
    isStripDetector = x_det.attr<bool>( _Unicode(isStripDetector) ) ;

  } catch(const std::runtime_error &){}

  //=========  loop over layer elements in xml  ======================================

  cout << "Building vertex barrel detector " << det_name << "..." << endl;
  for(xml_coll_t c(e, _U(layer) ); c; ++c)  {
    
    xml_comp_t x_layer( c );
    
    // child elements: ladder, sensitive and periphery
    xml_comp_t x_sensitive( x_layer.child( _U(sensitive) ));
    xml_comp_t x_ladder(  x_layer.child( _U(ladder)  ));
    xml_comp_t x_periphery(  x_layer.child( _U(module_envelope)  ));
    xml_comp_t x_readout(  x_layer.child( _U(readout)  ));
    int layer_id = x_layer.id();
    int nLadders = x_layer.attr<double>(  _Unicode(nLadders) ) ;

    double dphi = 2.*M_PI / double(nLadders);

    std::string layername = det_name+_toString(layer_id,"_layer%d");
    

    // --- create an assembly and DetElement for the layer

    Assembly layer_assembly( det_name + "_layer_" +_toString(layer_id,"_%d") );
    DetElement layerDE( tracker , _toString(layer_id,"layer_%d"), x_det.id() );
    pv = assembly.placeVolume(  layer_assembly );
    pv.addPhysVolID("layer", layer_id );  
    layerDE.setPlacement( pv ) ;

    //--------------------------------

    // Support //
    double supp_offset    = x_ladder.offset();
    double supp_distance  = x_ladder.distance();
    double supp_zhalf     = x_ladder.length();
    double supp_width     = x_ladder.width();
    std::string supp_vis   = x_ladder.visStr();

    // Sensor //
    double sens_zhalf     = x_sensitive.length();
    double sens_offset    = x_sensitive.offset();
    double sens_distance  = x_sensitive.distance();
    double sens_thickness = x_sensitive.thickness();
    double sens_width     = x_sensitive.width();
    int sens_perLadder = x_sensitive.nModules();
    double sens_z_displace= x_sensitive.dz();
    double sens_length    = x_sensitive.z_length();

    std::string sens_vis  = x_sensitive.visStr() ;
    std::string sens_matS = x_sensitive.materialStr() ;

    // Periphery //
    double peri_width   = x_periphery.width();
    double peri_length  = x_periphery.length();    
    int peri_type       = x_periphery.type();   // 0: Don't use periphery, 1: put periphery in positive phi direction, -1: put periphery in negative phi direction
    std::string peri_vis= x_periphery.visStr();

    // Readout (e.g flex) //
    double readout_offset    = x_readout.offset();
    double readout_distance  = x_readout.distance();
    double readout_zhalf     = x_readout.length();
    std::string readout_vis  = x_readout.visStr();

    double phi0             = x_layer.phi0();
    double phi_tilt         = x_layer.phi_tilt();

    if( sens_distance < minRadius ) minRadius = sens_distance ;
    if( supp_distance < minRadius ) minRadius = supp_distance ;
    if( sens_zhalf < minZhalf ) minZhalf = sens_zhalf ;
    if( supp_zhalf < minZhalf ) minZhalf = supp_zhalf ;

    cout << "    Layer: " << to_string(layer_id) << ", zhalf: " << to_string(sens_zhalf) << ", minimum radius: " << to_string(minRadius) << endl;

    Material sens_mat     = theDetector.material( sens_matS ) ;

    //-------
    Box sens_box( sens_thickness/2., sens_width/2., sens_length/2. );
    Box peri_box_z( sens_thickness/2., sens_width/2., peri_length/2. );     // Periphery in z (making the sensor longer in z direction)
    Box peri_box_rphi( sens_thickness/2., peri_width/2., sens_length/2. );  // Periphery in r-phi (making the sensor longer in r-phi direction)
//    Box peri_box_rphi( sens_thickness/2., peri_width/2., sens_length/2. );  // Periphery in r-phi (making the sensor longer in r-phi direction)

    Volume sens_vol( layername+"_sens", sens_box, sens_mat );
    Volume peri_vol_z( layername+"_peri_z", peri_box_z, sens_mat );
    Volume peri_vol_rphi( layername+"_peri_rphi", peri_box_rphi, sens_mat );

    // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
    Vector3D u( 0. , 1. , 0. ) ;
    Vector3D v( 0. , 0. , 1. ) ;
    Vector3D n( 1. , 0. , 0. ) ;
    
    // Volume of stave-long support structures //
    double supp_total_thickness = 0.0;
    for(xml_coll_t mi(x_layer,_U(ladder)); mi; ++mi)  {
      xml_comp_t x_mod   = mi;
      xml_coll_t ci(x_mod,_U(support));
      for(ci.reset(); ci; ++ci)
        supp_total_thickness += xml_comp_t(ci).thickness();
    }

    // compute the inner and outer thicknesses that need to be assigned to the tracking surface
    // depending on wether the support is above or below the sensor
    double inner_thickness = ( sens_distance > supp_distance ?  ( sens_distance - supp_distance ) + sens_thickness/2  : sens_thickness/2 ) ;
    double outer_thickness = ( sens_distance > supp_distance ?    sens_thickness/2  :  ( supp_distance - sens_distance ) + supp_total_thickness - sens_thickness/2   ) ;

    SurfaceType type( SurfaceType::Sensitive ) ;

    if( isStripDetector )
      type.setProperty( SurfaceType::Measurement1D , true ) ;

    VolPlane surf( sens_vol , type , inner_thickness , outer_thickness , u,v,n ) ; //,o ) ;

    sens.setType("tracker");
    sens_vol.setSensitiveDetector(sens);

    sens_vol.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), sens_vis );
    peri_vol_z.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), peri_vis );
    peri_vol_rphi.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), peri_vis );

    //--------- loop over ladders ---------------------------
    for(int j=0; j<nLadders; ++j) {

      double phi = phi0 + j * dphi  ;


      RotationZYX rot( phi , 0, 0  ) ;
      RotationZYX rot2( phi_tilt , 0, 0  ) ; // Tilt not working yet
      std::vector<double> distances{supp_distance, sens_distance, readout_distance}; 
      double min_distance = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));

      std::string laddername = layername + _toString(j,"_ladder%d");
      // Assembly ladder_assembly(laddername);
      // pv = layer_assembly.placeVolume(ladder_assembly);

      // --- place support elements that are along whole stave -----
      int c_id, m_id = 0;
      double c_distance = supp_distance;
      double c_offset = readout_offset;
      double c_thick = 0.0;
      for(xml_coll_t mi(x_layer,_U(ladder)); mi; ++mi, ++m_id)  {
        xml_comp_t x_mod   = mi;
        xml_coll_t ci(x_mod,_U(support));
        for(ci.reset(), c_id=0; ci; ++ci, ++c_id){
          xml_comp_t comp     = ci;    // A component of the support
          c_thick      = comp.thickness();
          double c_width      = comp.width();
          c_offset     = readout_offset+comp.offset();
          Material c_mat      = theDetector.material(comp.materialStr());
          std::string c_name  = _toString(c_id,"_supp%d");
          Volume c_vol(layername+c_name, Box(c_thick/2.0,c_width/2.0, supp_zhalf), c_mat);
          c_vol.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), comp.visStr() );

          pv = layer_assembly.placeVolume(c_vol, Transform3D( rot, Position( ( c_distance + c_thick/2. ) * cos(phi)  - c_offset * sin( phi ) ,
                                                                             ( c_distance + c_thick/2. ) * sin(phi)  + c_offset * cos( phi ) ,
                                                                                                                                        0. ) ) 
                                                 * Transform3D( rot2, Position( sin(phi_tilt)*(c_distance+c_offset*sin(phi)-min_distance), (c_distance-c_offset*cos(phi)-min_distance)*(cos(phi_tilt)-1) , 0.) )) ;
          c_distance += c_thick;          
        }
      }

      // --- place readout elements that are along whole stave (e.g flex) -----
      m_id = 0;
      c_distance = readout_distance;
      for(xml_coll_t mi(x_layer,_U(readout)); mi; ++mi, ++m_id)  {
        xml_comp_t x_mod   = mi;
        xml_coll_t ci(x_mod,_U(component));
        for(ci.reset(), c_id=0; ci; ++ci, ++c_id){
          xml_comp_t comp     = ci;    // A component of the readout
          c_thick      = comp.thickness();
          double c_width      = comp.width();
          c_offset     = readout_offset+comp.offset();
          Material c_mat      = theDetector.material(comp.materialStr());
          std::string c_name  = _toString(c_id,"_readout%d");
          Volume c_vol(layername+c_name, Box(c_thick/2.0,c_width/2.0, readout_zhalf), c_mat);
          c_vol.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), comp.visStr() );
          pv = layer_assembly.placeVolume(c_vol, Transform3D( rot, Position( ( c_distance + c_thick/2. ) * cos(phi)  - c_offset * sin( phi ) ,
                                                                                                       ( c_distance + c_thick/2. ) * sin(phi)  + c_offset * cos( phi ) ,
                                                                                                                                                                                                                                               0. ) ) * Transform3D( rot2, Position( sin(phi_tilt)*(c_distance+c_offset*sin(phi)-min_distance), (c_distance-c_offset*cos(phi)-min_distance)*(cos(phi_tilt)-1), 0.) ) );
          c_distance += c_thick;          
        }
      }

      // --- place sensitive -----
      c_thick = sens_thickness ;
      c_offset = sens_offset ;
      c_distance = sens_distance ;      
      // --- place all modules ----

      for(int iMod=0; iMod<sens_perLadder; ++iMod) {
        double z_pos = - sens_zhalf + sens_length/2.0 + iMod*(sens_z_displace + sens_length); // Move by sens_length+sens_z_displace from module to module

        std::string modulename = laddername + _toString(iMod,"_module%d");

        // Module
        pv = layer_assembly.placeVolume( sens_vol,Transform3D(rot, Position((c_distance + c_thick/2.0)*cos(phi) - c_offset*sin(phi),
                                                                            (c_distance + c_thick/2.0)*sin(phi) + c_offset*cos(phi),
                                                                            z_pos))
                                                 *Transform3D( rot2, Position( sin(phi_tilt)*(c_distance-c_offset*sin(phi)-min_distance), (c_distance+c_offset*cos(phi)-min_distance)*(cos(phi_tilt)-1), 0.) ) );
  
        pv.addPhysVolID( "module" , j ).addPhysVolID("sensor", iMod ) ;

        DetElement moduleDE( layerDE , modulename, x_det.id() );
        moduleDE.setPlacement( pv ) ;
        volSurfaceList( moduleDE )->push_back( surf ) ;
      }

      for(int iMod=0; iMod<sens_perLadder; ++iMod) {
        double z_pos = - sens_zhalf + sens_length/2.0 + iMod*(sens_z_displace + sens_length); // Move by sens_length+sens_z_displace from module to module

        // Periphery
        if(peri_type != 0){
            if(peri_length>0){
                layer_assembly.placeVolume( peri_vol_z,Transform3D(rot, Position((c_distance + c_thick/2.0)*cos(phi) - c_offset*sin(phi),
                                                                            (c_distance + c_thick/2.0)*sin(phi) + c_offset*cos(phi),
                                                                            z_pos - (sens_length/2.0 + peri_length/2.0)*peri_type) )
                                                 *Transform3D( rot2, Position( sin(phi_tilt)*(c_distance-c_offset*sin(phi)-min_distance), (c_distance+c_offset*cos(phi)-min_distance)*(cos(phi_tilt)-1), 0.) ) );
            }
            if(peri_width>0){
                layer_assembly.placeVolume( peri_vol_rphi,Transform3D(rot, Position((c_distance + c_thick/2.0)*cos(phi) - c_offset*sin(phi) - (sens_width/2.0 + peri_width/2.0)*peri_type*sin(phi),
                                                                            (c_distance + c_thick/2.0)*sin(phi) + c_offset*cos(phi) + (sens_width/2.0 + peri_width/2.0)*peri_type*cos(phi),
                                                                            z_pos) )
                                                 *Transform3D( rot2, Position( sin(phi_tilt)*(c_distance - c_offset*sin(phi) - (sens_width/2.0 + peri_width/2.0)*peri_type*sin(phi) -min_distance), (c_distance + c_offset*cos(phi) + (sens_width/2.0 + peri_width/2.0)*peri_type*cos(phi) -min_distance)*(cos(phi_tilt)-1), 0.) ) );
            }
        }
      }
      

      // DetElement ladderDE( layerDE ,  laddername , x_det.id() );
      // ladderDE.setPlacement( pv ) ;

      // volSurfaceList( ladderDE )->push_back( surf ) ;


      ///////////////////

      //get cellID and fill map< cellID of surface, vector of cellID of neighbouring surfaces >

      //encoding

      encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::barrel;
      encoder[lcio::LCTrackerCellID::layer()] = layer_id;
      encoder[lcio::LCTrackerCellID::module()] = nLadders;
      encoder[lcio::LCTrackerCellID::sensor()] = 0; // there is no sensor defintion in VertexBarrel at the moment

      dd4hep::long64 cellID = encoder.lowWord(); // 32 bits

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
    

    //    tracker.setVisAttributes(theDetector, x_det.visStr(),laddervol);
    
    // is this needed ??
    layer_assembly->GetShape()->ComputeBBox() ;

    //-----------------------------------
    //  store the data in an extension to be used for reconstruction
    ZPlanarData::LayerLayout thisLayer ;

    thisLayer.sensorsPerLadder = sens_perLadder ;
    thisLayer.sensorsZDisplace = sens_z_displace ;
    thisLayer.lengthSensor     = sens_length ;

    thisLayer.distanceSupport  = supp_distance ;
    thisLayer.offsetSupport    = supp_offset ;
    thisLayer.thicknessSupport = supp_total_thickness ;
    thisLayer.zHalfSupport     = supp_zhalf ;
    thisLayer.widthSupport     = supp_width ;

    thisLayer.distanceSensitive  = sens_distance ;
    thisLayer.offsetSensitive    = sens_offset ;
    thisLayer.thicknessSensitive = sens_thickness ;
    thisLayer.zHalfSensitive     = sens_zhalf ;
    thisLayer.widthSensitive     = sens_width ;

    thisLayer.ladderNumber =  nLadders ;
    thisLayer.phi0         =  phi0 ;

    thisLayer.periWidth    = peri_width ;
    thisLayer.periLength   = peri_length ;
    thisLayer.periType     = peri_type ;

    zPlanarData->layers.push_back( thisLayer ) ;
    //-----------------------------------
  }

#if 0  //-------- add an inscribing cylinder of air for tracking purposes -----------------
       //  this screw up the geometry and the material scan does not work anymore !!!!!?????

  double tube_thick =  1.0 * dd4hep::mm ;
  double inner_r    =  minRadius - 1.1 * tube_thick ;
  double outer_r    =  inner_r + tube_thick ;
  double z_half     =  minZhalf ; 
  
  Tube   tubeSolid (inner_r, outer_r, z_half ) ;
  Volume tube_vol( det_name+"_inner_cylinder_air", tubeSolid ,  theDetector.material("Air") ) ;
  
  assembly.placeVolume( tube_vol , Transform3D() ) ;
  
  Vector3D ocyl(  inner_r + 0.5*tube_thick , 0. , 0. ) ;
  
  VolCylinder cylSurf( tube_vol , SurfaceType( SurfaceType::Helper ) , 0.5*tube_thick  , 0.5*tube_thick , ocyl ) ;
  
  volSurfaceList( tracker )->push_back( cylSurf ) ;
  
#endif //----------------------------------------------------------------------------------



  tracker.addExtension< ZPlanarData >( zPlanarData ) ;
  tracker.addExtension< NeighbourSurfacesData >( neighbourSurfacesData ) ;


  Volume mother =  theDetector.pickMotherVolume( tracker ) ;

  //set the vis Attribute (added by Thorben Quast)
  assembly.setVisAttributes(theDetector, x_det.visStr());

  pv = mother.placeVolume(assembly);
  
  pv.addPhysVolID( "system", x_det.id() ).addPhysVolID("side",0 )  ;
  
  tracker.setPlacement(pv);
       
  assembly->GetShape()->ComputeBBox() ;

  return tracker;
}

DECLARE_DETELEMENT(VertexBarrel_IDEA_o1_v01,create_element)

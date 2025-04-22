// $Id: $
//====================================================================
//  Simple tracking detector made from planar sensors that are parallel
//  to the z-axis. There are two materials per ladder: one sensitive
//  and one support.
//  Copied from ZPlanarTracker_geo.cpp with addition of Z segmentation.
//--------------------------------------------------------------------
//
//  Author     : N.Bartosik
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"

#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include <exception>

#include <UTIL/BitField64.h>
#include <UTIL/BitSet32.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

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
  std::string  name  = x_det.nameStr();

  // put the whole detector into an assembly
  //  - should be replaced by an envelope volume ...

  Assembly assembly( name+"_assembly" );

  DetElement tracker( name, x_det.id()  ) ;

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

  } catch(std::runtime_error& ){}

  //=========  loop over layer elements in xml  ======================================

  for(xml_coll_t c(e, _U(layer) ); c; ++c)  {

    xml_comp_t x_layer( c );

    // child elements: ladder and sensitive
    xml_comp_t x_sensitive( x_layer.child( _U(sensitive) ));
    xml_comp_t x_ladder(  x_layer.child( _U(ladder)  ));

    int layer_id = x_layer.id();
    int nLadders = x_layer.attr<double>(  _Unicode(nLadders) ) ;

    double dphi = 2.*M_PI / double(nLadders);

    std::string layername = name+_toString(layer_id,"_layer%d");


    // --- create an assembly and DetElement for the layer

    Assembly layer_assembly( "layer_assembly" +_toString(layer_id,"_%d") );

    DetElement layerDE( tracker , _toString(layer_id,"layer_%d"), x_det.id() );

    pv = assembly.placeVolume(  layer_assembly );

    pv.addPhysVolID("layer", layer_id );

    layerDE.setPlacement( pv ) ;

    //--------------------------------

    double supp_zhalf     = x_ladder.length();
    double supp_offset    = x_ladder.offset();
    double supp_distance  = x_ladder.distance();
    double supp_thickness = x_ladder.thickness();
    double supp_width     = x_ladder.width();

    std::string supp_vis  = x_ladder.visStr() ;
    std::string supp_matS = x_ladder.materialStr() ;

    double sens_zhalf     = x_sensitive.length();
    double sens_offset    = x_sensitive.offset();
    double sens_distance  = x_sensitive.distance();
    double sens_thickness = x_sensitive.thickness();
    double sens_width     = x_sensitive.width();
    int    sens_nmodules  = x_sensitive.nmodules();
    double sens_modlength = sens_zhalf*2.0 / sens_nmodules;

    std::string sens_vis  = x_sensitive.visStr() ;
    std::string sens_matS = x_sensitive.materialStr() ;

    double phi0           = x_layer.phi0() ;


    if( sens_distance < minRadius ) minRadius = sens_distance ;
    if( supp_distance < minRadius ) minRadius = supp_distance ;
    if( sens_zhalf < minZhalf ) minZhalf = sens_zhalf ;
    if( supp_zhalf < minZhalf ) minZhalf = supp_zhalf ;

    //-----------------------------------
    //  store the data in an extension to be used for reconstruction
    ZPlanarData::LayerLayout thisLayer ;

    thisLayer.sensorsPerLadder =  sens_nmodules ;
    thisLayer.lengthSensor     =  sens_modlength ;

    thisLayer.distanceSupport  = supp_distance ;
    thisLayer.offsetSupport    = supp_offset ;
    thisLayer.thicknessSupport = supp_thickness ;
    thisLayer.zHalfSupport     = supp_zhalf ;
    thisLayer.widthSupport     = supp_width ;

    thisLayer.distanceSensitive  = sens_distance ;
    thisLayer.offsetSensitive    = sens_offset ;
    thisLayer.thicknessSensitive = sens_thickness ;
    thisLayer.zHalfSensitive     = sens_zhalf ;
    thisLayer.widthSensitive     = sens_width ;

    thisLayer.ladderNumber =  nLadders ;
    thisLayer.phi0         =  phi0 ;

    zPlanarData->layers.push_back( thisLayer ) ;
    //-----------------------------------


    Material supp_mat     = theDetector.material( supp_matS ) ;
    Material sens_mat     = theDetector.material( sens_matS ) ;


    //-------
    Box supp_box( supp_thickness/2., supp_width/2., supp_zhalf );
    Box sens_box( sens_thickness/2., sens_width/2., sens_modlength/2.0 );

    Volume supp_vol( layername+"_supp", supp_box, supp_mat  );
    Volume sens_vol( layername+"_sens", sens_box, sens_mat );


    // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
    Vector3D u( 0. , 1. , 0. ) ;
    Vector3D v( 0. , 0. , 1. ) ;
    Vector3D n( 1. , 0. , 0. ) ;

    //    Vector3D o( 0. , 0. , 0. ) ;

    // compute the inner and outer thicknesses that need to be assigned to the tracking surface
    // depending on wether the support is above or below the sensor
    double inner_thickness = ( sens_distance > supp_distance ?  ( sens_distance - supp_distance ) + sens_thickness/2  : sens_thickness/2 ) ;
    double outer_thickness = ( sens_distance > supp_distance ?    sens_thickness/2  :  ( supp_distance - sens_distance ) + supp_thickness - sens_thickness/2   ) ;

    SurfaceType type( SurfaceType::Sensitive ) ;

    if( isStripDetector )
      type.setProperty( SurfaceType::Measurement1D , true ) ;

    VolPlane surf( sens_vol , type , inner_thickness , outer_thickness , u,v,n ) ; //,o ) ;

    //--------------------------------------------

    sens.setType("tracker");
    sens_vol.setSensitiveDetector(sens);

    sens_vol.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), sens_vis );
    supp_vol.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), supp_vis );


    //--------- loop over ladders ---------------------------

    for(int j=0; j<nLadders; ++j) {

      double phi = phi0 + j * dphi  ;
      RotationZYX rot( phi , 0, 0  ) ;


      // --- place support -----
      double lthick = supp_thickness ;
      double radius = supp_distance ;
      double offset = supp_offset ;

      pv = layer_assembly.placeVolume( supp_vol,Transform3D( rot, Position( ( radius + lthick/2. ) * cos(phi)  - offset * sin( phi ) ,
									    ( radius + lthick/2. ) * sin(phi)  + offset * cos( phi ) ,
									    0. ) ));


      // --- place sensitive -----
      lthick = sens_thickness ;
      radius = sens_distance ;
      offset = sens_offset ;

      //--------- loop over sensors ---------------------------
      for(int s=0; s<sens_nmodules; ++s) {

        pv = layer_assembly.placeVolume( sens_vol,Transform3D( rot, Position( ( radius + lthick/2. ) * cos(phi)  - offset * sin( phi ) ,
          ( radius + lthick/2. ) * sin(phi)  + offset * cos( phi ) ,
          -sens_zhalf + sens_modlength*(float(s)+0.5) ) ));



        //      pv.addPhysVolID("layer", layer_id ).addPhysVolID( "module" , j ).addPhysVolID("sensor", 0 )   ;
        pv.addPhysVolID( "module" , j ).addPhysVolID("sensor", s )   ;

        std::string sensorname = layername + _toString(j,"_ladder%d") + _toString(s,"_sensor%d");
        DetElement   sensorDE( layerDE ,  sensorname , x_det.id() );
        sensorDE.setPlacement( pv ) ;

        volSurfaceList( sensorDE )->push_back( surf ) ;

        ///////////////////

        //get cellID and fill map< cellID of surface, vector of cellID of neighbouring surfaces >

        //encoding

        encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::barrel;
        encoder[lcio::LCTrackerCellID::layer()] = layer_id;
        encoder[lcio::LCTrackerCellID::module()] = j;
        encoder[lcio::LCTrackerCellID::sensor()] = s;

        dd4hep::long64 cellID = encoder.lowWord(); // 32 bits

        //compute neighbours

        int n_neighbours_ladder = 1; // 1 gives the adjacent modules
        int n_neighbours_sensor = 1;
        int newladder=0;
        int newsensor=0;

        for(int iladder = -n_neighbours_ladder; iladder <= n_neighbours_ladder; iladder++){ // neighbouring ladders
          newladder = j + iladder;
          //compute special case at the boundary
          if (newladder < 0) newladder = nLadders + newladder;
          if (newladder >= nLadders) newladder = newladder - nLadders;

          for(int isensor = -n_neighbours_sensor; isensor <= n_neighbours_sensor; isensor++){ // neighbouring sensors
            if (iladder==0 && isensor==0) continue; // sensor we started with
            newsensor = s + isensor;
            // skip sensors outside boundaries
            if (newsensor < 0) continue;
            if (newsensor >= sens_nmodules) continue;


            //encoding
            encoder[lcio::LCTrackerCellID::module()] = newladder;
            encoder[lcio::LCTrackerCellID::sensor()] = newsensor;

            neighbourSurfacesData->sameLayer[cellID].push_back(encoder.lowWord());
          }

        }

        ///////////////////

      }

    }


    //    tracker.setVisAttributes(theDetector, x_det.visStr(),laddervol);

    // is this needed ??
    layer_assembly->GetShape()->ComputeBBox() ;

  }

#if 0  //-------- add an inscribing cylinder of air for tracking purposes -----------------
       //  this screw up the geometry and the material scan does not work anymore !!!!!?????

  double tube_thick =  1.0 * dd4hep::mm ;
  double inner_r    =  minRadius - 1.1 * tube_thick ;
  double outer_r    =  inner_r + tube_thick ;
  double z_half     =  minZhalf ;

  Tube   tubeSolid (inner_r, outer_r, z_half ) ;
  Volume tube_vol( name+"_inner_cylinder_air", tubeSolid ,  theDetector.material("Air") ) ;

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

DECLARE_DETELEMENT(ZSegmentedPlanarTracker,create_element)

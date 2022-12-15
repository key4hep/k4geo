// $Id: $
//====================================================================
//  Based on ZPlanarTracker module from F. Gaede

//  Add description!!!
//  Simple tracking detector made from planar sensors that are parallel
//  to the z-axis. There are two materials per ladder: one sensitive
//  and one support.
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

  } catch(const std::runtime_error &){}

  //=========  loop over layer elements in xml  ======================================

  for(xml_coll_t c(e, _U(layer) ); c; ++c)  {
    
    xml_comp_t x_layer( c );
    
    // child elements: ladder, sensitive and periphery
    xml_comp_t x_sensitive( x_layer.child( _U(sensitive) ));
    xml_comp_t x_ladder(  x_layer.child( _U(ladder)  ));
    xml_comp_t x_periphery(  x_layer.child( _U(module_envelope)  ));
    xml_comp_t x_flex(  x_layer.child( _U(module_component)  ));    
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

    // Support //
    double supp_zhalf     = x_ladder.length();
    double supp_offset    = x_ladder.offset();
    double supp_distance  = x_ladder.distance();
    double supp_thickness = x_ladder.thickness();
    double supp_width     = x_ladder.width();

    std::string supp_vis  = x_ladder.visStr() ;
    std::string supp_matS = x_ladder.materialStr() ;
    
    // Sensor //
    double sens_zhalf     = x_sensitive.length();
    double sens_offset    = x_sensitive.offset();
    double sens_distance  = x_sensitive.distance();
    double sens_thickness = x_sensitive.thickness();
    double sens_width     = x_sensitive.width();

    double sens_perLadder;
    double sens_z_displace;
    double sens_length;    

    std::string sens_vis  = x_sensitive.visStr() ;
    std::string sens_matS = x_sensitive.materialStr() ;

    try{
      sens_perLadder = x_sensitive.nModules();
      sens_z_displace= x_sensitive.dz();
      sens_length    = x_sensitive.z_length();
    }
    catch(...){
      sens_perLadder = 1;
      sens_z_displace= 0.0;
      sens_length    = 2.*sens_zhalf;
      std::cout << "nModules, dz or z_length not found in ZPlanarTracker instantiation, assuming nModules=1, dz=0.0 and z_length=2*length (one sensor per stave)" << std::endl;
    }

    // Periphery //
    double peri_width   = x_periphery.width();
    double peri_length  = x_periphery.length();    
    int peri_type       = x_periphery.type();   // 0: Don't use periphery, 1: put periphery in positive phi direction, -1: put periphery in negative phi direction
    std::string peri_vis= x_periphery.visStr();

    // Flex //
    int build_flex          = x_flex.type();    // 0: Don't build flex, 1: Build flex
    double flex_thickness   = x_flex.thickness();
    double flex_width       = x_flex.width();
    double flex_offset      = x_flex.offset();
    std::string flex_matS   = x_flex.materialStr();
    std::string flex_vis    = x_flex.visStr();
    std::cout << "type: " << build_flex << ", thickness: " << flex_thickness << ", width: "<< flex_width << ", offset: "<< flex_offset << ", mat: " << flex_matS << ", vis: " << flex_vis << std::endl;

    
    double phi0           = x_layer.phi0() ;


    if( sens_distance < minRadius ) minRadius = sens_distance ;
    if( supp_distance < minRadius ) minRadius = supp_distance ;
    if( sens_zhalf < minZhalf ) minZhalf = sens_zhalf ;
    if( supp_zhalf < minZhalf ) minZhalf = supp_zhalf ;

    //-----------------------------------
    //  store the data in an extension to be used for reconstruction 
    ZPlanarData::LayerLayout thisLayer ;

    thisLayer.sensorsPerLadder = sens_perLadder ;
    thisLayer.sensorsZDisplace = sens_z_displace ;
    thisLayer.lengthSensor     = sens_length ;
    
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

    thisLayer.periWidth    = peri_width ;
    thisLayer.periLength   = peri_length ;
    thisLayer.periType     = peri_type ;

    thisLayer.buildFlex         = build_flex ;
    thisLayer.flexThickness     = flex_thickness ;
    thisLayer.flexWidth         = flex_width ;
    thisLayer.flexOffset        = flex_offset ;
    
    zPlanarData->layers.push_back( thisLayer ) ;
    //-----------------------------------


    Material supp_mat     = theDetector.material( supp_matS ) ;
    Material sens_mat     = theDetector.material( sens_matS ) ;
    Material flex_mat     = theDetector.material( flex_matS ) ;

    //-------
    Box sens_box( sens_thickness/2., sens_width/2., sens_length/2. );
    Box supp_box( supp_thickness/2., supp_width/2., supp_zhalf );
    Box flex_box( flex_thickness/2., flex_width/2., supp_zhalf );
    Box peri_box_z( sens_thickness/2., sens_width/2., peri_length/2. );     // Periphery in z (making the sensor longer in z direction)
    Box peri_box_rphi( sens_thickness/2., peri_width/2., sens_length/2. );  // Periphery in r-phi (making the sensor longer in r-phi direction)
//    Box peri_box_rphi( sens_thickness/2., peri_width/2., sens_length/2. );  // Periphery in r-phi (making the sensor longer in r-phi direction)

    Volume supp_vol( layername+"_supp", supp_box, supp_mat  );
    Volume flex_vol( layername+"_flex", flex_box, flex_mat  );
    Volume sens_vol( layername+"_sens", sens_box, sens_mat );
    Volume peri_vol_z( layername+"_peri_z", peri_box_z, sens_mat );
    Volume peri_vol_rphi( layername+"_peri_rphi", peri_box_rphi, sens_mat );

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
    flex_vol.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), flex_vis );
    peri_vol_z.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), peri_vis );
    peri_vol_rphi.setAttributes( theDetector, x_det.regionStr(), x_det.limitsStr(), peri_vis );

    //--------- loop over ladders ---------------------------

    for(int j=0; j<nLadders; ++j) {

      double phi = phi0 + j * dphi  ;

      std::string laddername = layername + _toString(j,"_ladder%d");

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
      
      // --- place all modules ----
      for(int iMod=0; iMod<sens_perLadder; ++iMod) {
          double z_pos = - sens_zhalf + sens_length/2.0 + iMod*(sens_length +sens_z_displace); // Move by sens_length+sens_z_displace from module to module

//          for(int iSens_x; iSens_x<nSens_x; ++nSens_x){
//              for(int iSens_y; iSens_y<nSens_y; ++nSens_y){
                // Place sensitive part
                //
                // Place left periphery in x
                //
                // Place right periphery in x
                //
                // Place lower periphery in y
                //
                // Place higher periphery in y
                //
                // Place corner

//              }
//          }

          // Module
          pv = layer_assembly.placeVolume( sens_vol,Transform3D(rot, Position((radius + lthick/2.0)*cos(phi) - offset*sin(phi),
                                                                              (radius + lthick/2.0)*sin(phi) + offset*cos(phi),
                                                                              z_pos
                                                                              )) );

          // Periphery
          if(peri_type != 0){
              if(peri_length>0)
                  layer_assembly.placeVolume( peri_vol_z,Transform3D(rot, Position((radius + lthick/2.0)*cos(phi) - offset*sin(phi),
                                                                              (radius + lthick/2.0)*sin(phi) + offset*cos(phi),
                                                                              z_pos - (sens_length/2.0 + peri_length/2.0)*peri_type
                                                                              )) );  
              if(peri_width>0)
                  layer_assembly.placeVolume( peri_vol_rphi,Transform3D(rot, Position((radius + lthick/2.0)*cos(phi) - offset*sin(phi) - (sens_width/2.0 + peri_width/2.0)*peri_type*sin(phi),
                                                                              (radius + lthick/2.0)*sin(phi) + offset*cos(phi) + (sens_width/2.0 + peri_width/2.0)*peri_type*cos(phi),
                                                                              z_pos
                                                                              )) );
          }

          pv.addPhysVolID("layer", layer_id ).addPhysVolID( "module" , j ).addPhysVolID("sensor", iMod ) ;

          std::string modulename = laddername + _toString(iMod,"_module%d");
          DetElement   ladderDE( layerDE ,  modulename , x_det.id() );
          ladderDE.setPlacement( pv ) ;


          volSurfaceList( ladderDE )->push_back( surf ) ;
          std::cout << "Making sensor " << iMod << " of stave " << j << " with z position of " << std::to_string(z_pos) << std::endl;
      }

      // --- place flex ----
      if(build_flex == 1)
          layer_assembly.placeVolume( flex_vol,Transform3D( rot, Position( ( radius + lthick/2. + flex_thickness/2. ) * cos(phi)  - (offset + flex_offset) * sin( phi ) ,
                                        ( radius + lthick/2. + flex_thickness/2. ) * sin(phi)  + (offset + flex_offset) * cos( phi ) ,
                                        0. ) ));

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

      ///////////////////
      


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

DECLARE_DETELEMENT(VertexBarrel_o1_v01,create_element)

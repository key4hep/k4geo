// $Id: $
//====================================================================
//   A polyhedral endcap structure with surfaces along the sides
//   Used for tracking purposes (e.g. track state at the calorimeter face)
//  
//--------------------------------------------------------------------
//
//  Author     : F.Gaede
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"

#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::Transform3D;
using dd4hep::Volume;

using dd4hep::xml::_toString;

using dd4hep::rec::SurfaceType;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolPlane;
using dd4hep::rec::volSurfaceList;

static Ref_t create_element(Detector& theDetector, xml_h element, SensitiveDetector sens)  {
  
  xml_det_t    x_det = element;
  std::string  name  = x_det.nameStr();
  
  DetElement sdet( name, x_det.id()  ) ;
  
  PlacedVolume pv;
  

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;
  
  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
  
  //-----------------------------------------------------------------------------------
  
  // ----- read xml ----------------------

  xml_dim_t dim = x_det.dimensions();

  double inner_r    =  dim.rmin() ;
  double outer_r    =  dim.rmax() ;
  double z0          = dim.z0() ;
  double z1          = dim.z1() ;
  //  double phi0       =  dim.phi0() ;

  unsigned nsides   =  dim.nsides();
  
  double thick      =  z1 -  z0 ;
  double zpos       =  z0 + thick/2. ;


  Material mat      = envelope.material() ;      
  //--------------------------------------

  sens.setType("tracker");

  // base vectors for surfaces:
  Vector3D u(0,1,0) ;
  Vector3D v(1,0,0) ;
  Vector3D n(0,0,1) ;
  
  
  PolyhedraRegular phSolid ( nsides, inner_r, outer_r , 0.5*thick ) ;

  //==============================================================================

  DetElement fwdDE( sdet, name + std::string( "_fwd" )  , x_det.id() );

  Volume phVol( name + std::string("_vol") , phSolid ,  mat ) ;
    
  phVol.setSensitiveDetector(sens);
  
  //fixme: the drawing of endcap surfaces in a polyhedral shape does not work right now 
  //       -> set surface to be invisible for now
  VolPlane surf( phVol,SurfaceType(SurfaceType::Sensitive, 
                                   SurfaceType::Invisible), 
                 thick/4., thick/4., u,v,n ) ;
  
  volSurfaceList( fwdDE  )->push_back(  surf ) ;

  pv = envelope.placeVolume( phVol , Position( 0., 0., zpos ) ) ;

  pv.addPhysVolID("layer", 0 ).addPhysVolID( "side" , +1 )  ;
    
  fwdDE.setPlacement( pv ) ;

  //==============================================================================

  DetElement bwdDE( sdet, name + std::string( "_bwd" )  , x_det.id() );
  
  //  Volume phVol( name + std::string("_bwd") , phSolid ,  mat ) ;
  
  phVol.setSensitiveDetector(sens);
  
  
  volSurfaceList( bwdDE  )->push_back(  surf ) ;

  pv = envelope.placeVolume( phVol , Position( 0., 0., -zpos ) ) ;

  pv.addPhysVolID("layer", 0 ).addPhysVolID( "side" , -1 )  ;
    
  bwdDE.setPlacement( pv ) ;
  
  //--------------------------------------
  
  return sdet ;
}

DECLARE_DETELEMENT( PolyhedralEndcapSurfaces,create_element)

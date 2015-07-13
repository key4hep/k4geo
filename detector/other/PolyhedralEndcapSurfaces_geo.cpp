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


using namespace DD4hep;
using namespace DD4hep::Geometry;
using namespace DD4hep::DDRec ;
using namespace DDSurfaces ;

static Ref_t create_element(LCDD& lcdd, xml_h element, SensitiveDetector sens)  {
  
  xml_det_t    x_det = element;
  std::string  name  = x_det.nameStr();
  
  DetElement sdet( name, x_det.id()  ) ;
  
  PlacedVolume pv;
  

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = XML::createPlacedEnvelope( lcdd,  element , sdet ) ;
  
  if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;
  
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
  DDSurfaces::Vector3D u(1,0,0) ;
  DDSurfaces::Vector3D v(0,1,0) ;
  DDSurfaces::Vector3D n(0,0,1) ;
  
  
  PolyhedraRegular phSolid ( nsides, inner_r, outer_r , 0.5*thick ) ;

  //==============================================================================

  DetElement fwdDE( sdet, name + std::string( "_fwd" )  , x_det.id() );

  Volume phVol( name + std::string("_vol") , phSolid ,  mat ) ;
    
  phVol.setSensitiveDetector(sens);
  
  DDRec::VolPlane surf( phVol,DDSurfaces::SurfaceType(DDSurfaces::SurfaceType::Sensitive), thick/4., thick/4., u,v,n ) ;
  
  DDRec::volSurfaceList( fwdDE  )->push_back(  surf ) ;

  pv = envelope.placeVolume( phVol , Position( 0., 0., zpos ) ) ;

  pv.addPhysVolID("layer", 0 ).addPhysVolID( "side" , +1 )  ;
    
  fwdDE.setPlacement( pv ) ;

  //==============================================================================

  DetElement bwdDE( sdet, name + std::string( "_bwd" )  , x_det.id() );
  
  //  Volume phVol( name + std::string("_bwd") , phSolid ,  mat ) ;
  
  phVol.setSensitiveDetector(sens);
  
  // DDRec::VolPlane surf( phVol,DDSurfaces::SurfaceType(DDSurfaces::SurfaceType::Sensitive), thick/4., thick/4., u,v,n ) ;
  
  DDRec::volSurfaceList( bwdDE  )->push_back(  surf ) ;

  pv = envelope.placeVolume( phVol , Position( 0., 0., zpos ) ) ;

  pv.addPhysVolID("layer", 0 ).addPhysVolID( "side" , -1 )  ;
    
  bwdDE.setPlacement( pv ) ;
  
  //--------------------------------------
  
  return sdet ;
}

DECLARE_DETELEMENT( PolyhedralEndcapSurfaces,create_element)

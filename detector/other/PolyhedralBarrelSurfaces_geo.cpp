// $Id: $
//====================================================================
//   A polyhedral barrel structure with surfaces along the sides
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
  double z_half     =  dim.zhalf() ;
  double phi0       =  dim.phi0() ;
  unsigned nsides     =  dim.numsides();
  
  double alpha      =  M_PI / double( nsides ) ;

  double thick      =  outer_r - inner_r ;

  Material mat      = envelope.material() ;      
  //--------------------------------------

  double width_half = tan( alpha ) * inner_r ;


  sens.setType("tracker");

  // base vectors for surfaces:
  DDSurfaces::Vector3D u(0,1,0) ;
  DDSurfaces::Vector3D v(0,0,1) ;
  DDSurfaces::Vector3D n(1,0,0) ;
  
  for(unsigned i=0 ; i < nsides ; ++i ){

    DetElement moduleDE( sdet, name + _toString( i,"_module_%d")  , x_det.id() );
 
    double phi = phi0 + i * 2. * alpha  ;

    RotationZYX rot( phi , 0, 0  ) ;

    Box boxSolid ( thick/4. , width_half , z_half ) ;

    Volume boxVol( name + _toString( i,"_module_%d") , boxSolid ,  mat ) ;

    boxVol.setSensitiveDetector(sens);

    DDRec::VolPlane surf( boxVol,DDSurfaces::SurfaceType(DDSurfaces::SurfaceType::Sensitive), thick/4., thick/4., u,v,n ) ;
	      
    DDRec::volSurfaceList( moduleDE  )->push_back(  surf ) ;

    pv = envelope.placeVolume( boxVol , Transform3D( rot, Position( ( inner_r + thick/2. ) * cos( phi ) ,
								    ( inner_r + thick/2. ) * sin( phi ) , 
								    0  ) ) ) ;

    pv.addPhysVolID("layer", 0 ).addPhysVolID( "module" , i )  ;
    
    moduleDE.setPlacement( pv ) ;

  }



  //--------------------------------------
  
  return sdet ;
}

DECLARE_DETELEMENT( PolyhedralBarrelSurfaces,create_element)

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

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
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
  Vector3D u(0,1,0) ;
  Vector3D v(0,0,1) ;
  Vector3D n(1,0,0) ;
  
  for(unsigned i=0 ; i < nsides ; ++i ){

    DetElement moduleDE( sdet, name + _toString( i,"_module_%d")  , x_det.id() );
 
    double phi = phi0 + i * 2. * alpha  ;

    RotationZYX rot( phi , 0, 0  ) ;

    Box boxSolid ( thick/4. , width_half , z_half ) ;

    Volume boxVol( name + _toString( i,"_module_%d") , boxSolid ,  mat ) ;

    boxVol.setSensitiveDetector(sens);

    VolPlane surf( boxVol,SurfaceType(SurfaceType::Sensitive), thick/4., thick/4., u,v,n ) ;
	      
    volSurfaceList( moduleDE  )->push_back(  surf ) ;

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

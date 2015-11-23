#ifndef lcgeo_h
#define lcgeo_h

// file for some global configuration of lcgeo 
// so far only version number:

// define version macros for Lcgeo
#define LCGEO_MAJOR_VERSION 0
#define LCGEO_MINOR_VERSION 5

#define LCGEO_VERSION_GE( MAJV , MINV )  ( (  LCGEO_MAJOR_VERSION  > MAJV ) || ( (LCGEO_MAJOR_VERSION==MAJV) && ( LCGEO_MINOR_VERSION >= MINV ) ) )

#define LCGEO_VERSION_GT( MAJV , MINV )  ( (  LCGEO_MAJOR_VERSION  > MAJV ) || ( (LCGEO_MAJOR_VERSION==MAJV) && ( LCGEO_MINOR_VERSION >  MINV ) ) )

namespace lcgeo {

  /// return a string with the current lcgeo version in the form vXX-YY.
  inline std::string versionString(){
    std::string vs("vXX-YY") ;
    std::sprintf( &vs[0] , "v%2.2d-%2.2d", LCGEO_MAJOR_VERSION, LCGEO_MINOR_VERSION  ) ;
    return vs ;
  }
  
}

#endif

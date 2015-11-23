//====================================================================
//  DDSim - LC simulation based on DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//====================================================================
#ifndef Exceptions_h
#define Exceptions_h

#include <exception> 

namespace lcgeo {
  
  //define some exception to throw
  class GeometryException : public std::exception{
    
  protected:
    std::string message;
    GeometryException() {  /*no_op*/ ; } 
    
  public: 
    GeometryException( std::string text ){
      message = "GeometryException: " + text ;
    }
    virtual const char* what() const  throw() { return  message.c_str() ; } 
    virtual ~GeometryException() throw() { /*no_op*/; } 
  }; 
}
#endif

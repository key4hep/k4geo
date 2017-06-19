#include "XML/XMLDetector.h"
#include <string>

namespace {
  /** Wrapper class to replace the Database class used in Mokka to read the parameters.
   *  Assumes parameters are stored as attributes of the corresponding xml element.
   */
  struct XMLHandlerDB{
    xml_comp_t x_det ;
    /** C'tor initializes the handle */
  XMLHandlerDB(xml_comp_t det) : x_det(det) {}
  
    double fetchDouble( const char* _name){ return  x_det.attr<double>( dd4hep::xml::Strng_t(_name) )  ; }

    int    fetchInt( const char* _name){ return  x_det.attr<int>( dd4hep::xml::Strng_t(_name) ) ; }

    std::string fetchString( const char* _name){ return  x_det.attr<std::string>( dd4hep::xml::Strng_t(_name) ) ;}

    /** allow this to be used as a 'pointer' ( as was used for Mokka Database object)*/
    XMLHandlerDB* operator->() { return this ; }
  };

}

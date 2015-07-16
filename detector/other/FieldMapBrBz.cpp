
#include "FieldMapBrBz.h"


#include "DD4hep/Handle.inl"
#include <DD4hep/FieldTypes.h>


#include <DD4hep/DetFactoryHelper.h>
#include <XML/Utilities.h>


#include <TFile.h>


#include <string>
#include <stdexcept>

DD4HEP_INSTANTIATE_HANDLE(FieldMapBrBz);

FieldMapBrBz::FieldMapBrBz() {
  type = DD4hep::Geometry::CartesianField::MAGNETIC;
} //ctor


void FieldMapBrBz::fieldComponents(const double* , double* field) {
  field[0] = field[1] = 0.0;
  field[2] = 300 * dd4hep::tesla;
}


static DD4hep::Geometry::Ref_t create_FieldMap_rzBrBz(DD4hep::Geometry::LCDD& lcdd,
						      DD4hep::XML::Handle_t handle ) {
  DD4hep::XML::Component xmlParameter(handle);
  bool hasFilename = xmlParameter.hasAttr(_Unicode(filename));

  if (!hasFilename) {
    std::stringstream error;
    error << "FieldMap[ERROR]: For a FieldMap field at least the filename xml attribute MUST be set.";
    throw std::runtime_error(error.str());
  }
  std::string filename = xmlParameter.attr< std::string >(_Unicode(filename));

  DD4hep::Geometry::CartesianField obj;
  FieldMapBrBz* ptr = new FieldMapBrBz();
  ptr->_filename = filename;


  //Read the entries form the file in this place

  obj.assign(ptr, xmlParameter.nameStr(), xmlParameter.typeStr());
  return obj;
}
DECLARE_XMLELEMENT(FieldBrBz,create_FieldMap_rzBrBz)

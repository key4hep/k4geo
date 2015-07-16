#ifndef FieldMap_rzBrBz_h
#define FieldMap_rzBrBz_h 1

#include <DD4hep/FieldTypes.h>


#include <string>

class FieldMapBrBz: public DD4hep::Geometry::CartesianField::Object {
public:
  std::string _filename;

public:
  /// Initializing constructor
  FieldMapBrBz();
  /// Call to access the field components at a given location
  virtual void fieldComponents(const double* pos, double* field);
};


#endif // FieldMap_rzBrBz_h

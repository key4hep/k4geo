#ifndef FieldMap_XYZ_h
#define FieldMap_XYZ_h 1

#include <DD4hep/FieldTypes.h>

#include <string>
#include <vector>

class FieldMapXYZ: public DD4hep::Geometry::CartesianField::Object {
public:

  struct FieldValues_t {
    double Bx;
    double By;
    double Bz;
    FieldValues_t(double _Bx, double _By, double _Bz):
      Bx(_Bx), By(_By), Bz(_Bz) {}
  };

  int nX, nY, nZ;
  
  double xMin,xMax,xStep,xScale;
  double yMin,yMax,yStep,yScale;
  double zMin,zMax,zStep,zScale;
  
  double bScale;
  
  std::vector< FieldValues_t > fieldMap;

public:
  /// Initializing constructor
  FieldMapXYZ();
  
  /// Call to access the field components at a given location
  virtual void fieldComponents(const double* pos, double* field);
  
  /// Fiell the FieldMap from the the tree specified in the XML
  void fillFieldMapFromTree(const std::string& filename, const std::string& treename);
  
  int   GetGlobalIndex(double x, double y, double z);
  
};


#endif // FieldMap_XYZ_h

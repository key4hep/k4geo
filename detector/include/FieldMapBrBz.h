#ifndef FieldMap_rzBrBz_h
#define FieldMap_rzBrBz_h 1

#include <DD4hep/FieldTypes.h>

#include <string>
#include <vector>

class FieldMapBrBz: public DD4hep::Geometry::CartesianField::Object {
public:

  struct FieldValues_t {
    double Br;
    double Bz;
    FieldValues_t(double _Br, double _Bz):
      Br(_Br), Bz(_Bz) {}
  };

  int CoorsOrder;
  std::string StrCoorsOrder;
  int nRho, nZ;
  double rhoMin, rhoMax, rhoStep, rScale;
  double zMin,   zMax,   zStep,   zScale;
  double bScale;
  std::vector< FieldValues_t > fieldMap;

public:
  /// Initializing constructor
  FieldMapBrBz();
  /// Call to access the field components at a given location
  virtual void fieldComponents(const double* pos, double* field);
  /// Fill the FieldMap from the the tree specified in the XML
  void fillFieldMapFromTree(const std::string& filename, const std::string& treename, double coorUnits, double BfieldUnits);
  /// Get global index in the Field map
  int getGlobalIndex(const double r, const double z);
};


#endif // FieldMap_rzBrBz_h

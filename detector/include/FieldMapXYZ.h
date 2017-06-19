#ifndef FieldMap_XYZ_h
#define FieldMap_XYZ_h 1

#include <DD4hep/FieldTypes.h>

#include <string>
#include <vector>

class FieldMapXYZ: public dd4hep::CartesianField::Object {
public:

  struct FieldValues_t {
    double Bx;
    double By;
    double Bz;
    FieldValues_t(double _Bx, double _By, double _Bz):
      Bx(_Bx), By(_By), Bz(_Bz) {}
  };

  int coorsOrder;             // integer with the order with which variables are scanned in the fieldmap, 1(2) for RZ(ZR) order
  std::string  strCoorsOrder; // string  with the order with which variables are scanned in the fieldmap, RZ or ZR order
  std::string  ntupleName;    // tree name
  std::string  xVar;          // x  coordinate name in tree
  std::string  yVar;          // y  coordinate name in tree
  std::string  zVar;          // z  coordinate name in tree
  std::string  BxVar;         // Bx component  name in tree
  std::string  ByVar;         // By component  name in tree
  std::string  BzVar;         // Bz component  name in tree

  int nX, nY, nZ;                // bins in x, y and z coordinates in fieldmap
  int xOrdering;                 // x coordinate ordering, 1(-1) if from low-to-high (high-to-low)
  double xMin,xMax,xStep,xScale; // min, max, step-size and scale factor of x coordinate in fieldmap
  int yOrdering;                 // y coordinate ordering, 1(-1) if from low-to-high (high-to-low)
  double yMin,yMax,yStep,yScale; // min, max, step-size and scale factor of y coordinate in fieldmap
  int zOrdering;                 // z coordinate ordering, 1(-1) if from low-to-high (high-to-low)
  double zMin,zMax,zStep,zScale; // min, max, step-size and scale factor of z coordinate in fieldmap
  
  double bScale;                         //Bfield scale factor 
  std::vector< FieldValues_t > fieldMap; //List with the field map points

public:
  /// Initializing constructor
  FieldMapXYZ();
  
  /// Call to access the field components at a given location
  virtual void fieldComponents(const double* pos, double* field);
  /// Field the FieldMap from the the tree specified in the XML
  void fillFieldMapFromTree(const std::string& filename, double coorUnits, double BfieldUnits);
  /// Get global index in the Field map 
  int  getGlobalIndex(const int xBin, const int yBin, const int zBin);
  
};


#endif // FieldMap_XYZ_h

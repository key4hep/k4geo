#ifndef FieldMap_rzBrBz_h
#define FieldMap_rzBrBz_h 1

#include <DD4hep/FieldTypes.h>

#include <string>
#include <vector>

class FieldMapBrBz: public dd4hep::CartesianField::Object {
public:

  struct FieldValues_t {
    double Br;
    double Bz;
    FieldValues_t(double _Br, double _Bz):
      Br(_Br), Bz(_Bz) {}
  };

  int coorsOrder;             // integer with the order with which variables are scanned in the fieldmap, 1(2) for RZ(ZR) order
  std::string  strCoorsOrder; // string  with the order with which variables are scanned in the fieldmap, RZ or ZR order
  std::string  ntupleName;    // tree name
  std::string  rhoVar;        // rho  coordinate name in tree
  std::string  zVar;          // z    coordinate name in tree
  std::string  BrhoVar;       // Brho component  name in tree
  std::string  BzVar;         // Bz   component  name in tree

  int nRho, nZ;                            // bins in rho and z coordinates in fieldmap
  int rhoOrdering;                         // rho coordinate ordering, 1(-1) if from low-to-high (high-to-low)
  double rhoMin, rhoMax, rhoStep, rScale;  // min, max, step-size and scale factor of rho coordinate in fieldmap
  int zOrdering;                           // z   coordinate ordering, 1(-1) if from low-to-high (high-to-low)
  double zMin,   zMax,   zStep,   zScale;  // min, max, step-size and scale factor of z   coordinate in fieldmap

  double bScale;                           //Bfield scale factor
  std::vector< FieldValues_t > fieldMap;   //List with the field map points

public:
  /// Initializing constructor
  FieldMapBrBz();
  /// Call to access the field components at a given location
  virtual void fieldComponents(const double* pos, double* field);
  /// Field the FieldMap from the the tree specified in the XML
  void fillFieldMapFromTree(const std::string& filename, double coorUnits, double BfieldUnits);
  /// Get global index in the Field map
  int getGlobalIndex(const int rBin, const int zBin);
};


#endif // FieldMap_rzBrBz_h

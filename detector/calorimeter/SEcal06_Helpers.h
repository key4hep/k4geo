#ifndef __SECAL06_HELPERS_H__
#define __SECAL06_HELPERS_H__ 1

#include "DD4hep/DetFactoryHelper.h"

#include "DDRec/DetectorData.h"

#include "XML/Layering.h"
#include "XML/Utilities.h"

#include "DD4hep/Segmentations.h"

#include "DDSegmentation/MultiSegmentation.h"

#include <iostream>

#undef NDEBUG
#include <assert.h>

using std::cout;
using std::endl;

/*

we build general ECAL modules on a horizontal table. 

X-Y are on the table
Z is vertically upwards
slabs are aligned along X
the lower Z side will face the IP

2 shapes when looking down on table:

XYtype = 0:  (barrel and endcap)
     -----------------------------     ^ Y
     |                           |     |
     |                           |     |
     -----------------------------      ----> X

XYtype = 1:  (endcap)
     ----------------------
     |                      \
     |                        \
     |                          \
     |                           |   ^
     |                           |   | dY_kink
     -----------------------------   |

2 shapes when looking from side:

XZtype = 0:  (endcap)
   ---------------------------------
   |                               |
   |                               |
   ---------------------------------

XZtype = 1:  (barrel)
       --------------------------        ^ Z
      /                          \       |
     /                            \      |
    /                              \     |
    --------------------------------     ----> X

the slabs are prefectly aligned at -ve X side
"magic" wafers are at +ve X

local position: 0,0 defined as -veX,Y,Z corner

we start building from Z=0, the face nearest IP, moving in the +ve Z direction


D.Jeans update: 03/2015
dead space between slab end and module edge
other than 8-fold symmetry for barrel

 */

class SEcal06_Helpers {

 public:

  SEcal06_Helpers();

  ~SEcal06_Helpers() {
    delete _layering;
  }

  // SET THE PARAMETERS

  void setDet( xml_det_t* x_det ) {
    _x_det = x_det;
    _layering = new dd4hep::Layering(*x_det);
    _det_name  = x_det->nameStr();
  }

  void setAbsLayers( int nl1, double th1, int nl2, double th2, int nl3=0, double th3=0 );

  void setLayerConfig( std::vector < int > layerConfig ) {
    _layerConfig=layerConfig;
  }

  void setSegmentation(  dd4hep::Segmentation & seg ) {
    _geomseg = &seg;
  }

  void setSegmentation(  dd4hep::Segmentation * seg ) {
    _geomseg = seg;
  }

  void setNCells( int Ecal_cells_across_megatile, int Ecal_strips_across_megatile, int Ecal_strips_along_megatile) {
    _cells_across_megatile   =  Ecal_cells_across_megatile   ;
    _strips_across_megatile  =  Ecal_strips_across_megatile  ;
    _strips_along_megatile   =  Ecal_strips_along_megatile   ;
  }

  void setMagicMegatileStrategy ( int i ) {
    // magic megatile strategy
    //  0: no magic megatile at lend of slab
    //  1: magic megatile includes integer number of standard-sized cells
    //  2: last cell of magic megatile adjusted to fill up to end of slab (size between 1 and 2 times normal cell size)
    assert ( i>=0 && i<=2 );
    _magicMegatileStrategy = i;
  }

  void setPreshower( bool p ) {
    _preshower = p? 1 : 0; // do we have preshower layer? (if 1, first layer is sensitive; if 0, first layer is absorber)
  }

  void setCFthickness( double absWrap, double alvWall, double front, double back ) {
    _CF_absWrap = absWrap; // wrapping around absorber
    _CF_alvWall = alvWall; // alveolar wall
    _CF_front   = front;   // front (IP side) support plate
    _CF_back    = back;    // back support plate
  }

  void setModuleDimensions( int XYtype, // module shape in XY
			    int XZtype, // module shape in XZ
			    double dX_max, // maximum extent in X
			    double dY_kink = -999, // distance from lowerY edge to kink for XYtype=2
			    double angle = M_PI/4. // angle, if not rectangular. pi/4 for octagon
			    ) {
    if ( XYtype>0 ) assert ( XZtype==0 ); // endcap module has vertical edges
    if ( XZtype>0 ) assert ( XYtype==0 ); // barrel module should have rectangular shadow

    _module_XYtype  = XYtype;
    _module_XZtype  = XZtype;
    _module_dX_max  = dX_max;
    _module_dY_kink = dY_kink;
    _module_angle   = angle;

  }

  void setTowersUnits( std::vector <int> ntowers, // a vector of modules, containing the specified # of towers/alveolii
		       double towerWidth,
		       int unitsPerTower, // # megatiles/units per tower
		       double moduleDeadEdge, // dead region at edge of module
		       double towerDeadEdge, // dead region at edge of tower
		       double unitDeadEdge // width of dead region at edge of each unit (e.g. wafer's guard ring)
		       ) {
    _ntowers = ntowers;                 // number of towers (or alveoli) across Y
    _alveolus_total_dim_Y = towerWidth; // width of tower (including dead area)
    _unitsPerTower = unitsPerTower;     // how many units (megatiles/wafers) across one tower
    _moduleGap = moduleDeadEdge;        // distance from tower boundary to edge of unit
    _towerGap = towerDeadEdge;          // distance from tower boundary to edge of unit
    _unitDeadEdge = unitDeadEdge;       // width of dead area around edge of unit
  }

  void checkLayerConsistency();

  float getTotalThickness();

  void setTranslation( dd4hep::Position trans ) {_trans = trans;}

  // ---- this is the main workhorse
  void makeModule( dd4hep::Volume & mod_vol,  // the volume we'll fill
		   dd4hep::DetElement & stave_det, // the detector element
		   dd4hep::rec::LayeredCalorimeterData & caloData, // the reco data we'll fill
		   dd4hep::Detector & theDetector,
		   dd4hep::SensitiveDetector & sens
		   );

  void setPlugLength( float ll ) { _plugLength = ll; }


 private:

  void printSEcal06LayerInfo( dd4hep::rec::LayeredCalorimeterData::Layer & caloLayer);

  double getAbsThickness( unsigned int iAbsLay );

  int getNumberOfAbsorberLayers() {
    // total number of absorber layers
    return _nlayers1+_nlayers2+_nlayers2;
  }

  int getNumberOfStructuralAbsorberLayers() {
    // number of absorber layers in the module (not in the slabs)
    assert( _preshower==0 || _preshower==1 );
    return _preshower ? getNumberOfAbsorberLayers()/2 : (getNumberOfAbsorberLayers()-1)/2 ;
  }

  struct dimposXYStruct {
    double sizeX, sizeY;
    double posX, posY;
  };

  struct dxinfo {
    int normal_nX;
    
    double magic1_unitDX;
    int magic1_ncellsX;

    double magic2_unitDX;
  };

  xml_det_t* _x_det;
  dd4hep::Layering* _layering=NULL;
  std::string _det_name;

  std::vector <dimposXYStruct> getAbsPlateXYDimensions( double ztop=-999 );

  std::vector <dimposXYStruct> getSlabXYDimensions( double ztop=-999 );


  dd4hep::Position  getTranslatedPosition(double x, double y, double z) {
    return dd4hep::Position ( x, y, z ) + _trans;
  }

  dxinfo getNormalMagicUnitsInX( double dx_total, double dx_unit, double dx_cell, double dx_dead ,
				 int magicStrategy );

  void updateCaloLayers(double thickness,
			dd4hep::Material mat,
			bool isAbsorber,
			bool isSensitive,
			double cell_size_x=0, double cell_size_y=0,
			bool isFinal=false
			);

  const dd4hep::DDSegmentation::Segmentation* getSliceSegmentation( dd4hep::DDSegmentation::MultiSegmentation* multiSeg , int slice_number );

  dd4hep::Segmentation* _geomseg;

  dd4hep::Material _air_material;
  dd4hep::Material _carbon_fibre_material;
  dd4hep::Material _radiator_material;

  unsigned int _cells_across_megatile;
  unsigned int _strips_across_megatile;
  unsigned int _strips_along_megatile;   

  int _preshower;

  unsigned int _nlayers1;
  unsigned int _nlayers2;
  unsigned int _nlayers3;

  double _radiator_thickness1;
  double _radiator_thickness2;
  double _radiator_thickness3;

  std::vector < int > _layerConfig; // square or X or Y strips

  double _CF_absWrap;
  double _CF_alvWall;
  double _CF_front;
  double _CF_back;

  // overall module size
  int    _module_XYtype;
  int    _module_XZtype;
  double _module_dX_max;
  double _module_dY_total;
  double _module_dY_kink;
  double _module_angle;

  // internal details
  std::vector <int> _ntowers;
  double _towerGap;
  double _moduleGap;
  int    _unitsPerTower;
  double _unitDeadEdge;

  dd4hep::rec::LayeredCalorimeterData* _caloData;
  dd4hep::rec::LayeredCalorimeterData::Layer _caloLayer ; // this is the output info which is passed to reconstruction

  double _layer_thickness;
  double _layer_nRadiationLengths;
  double _layer_nInteractionLengths;

  double _totThick;

  double _alveolus_total_dim_Y;
  double _module_thickness;

  int _magicMegatileStrategy;

  dd4hep::Position _trans;

  std::vector <dimposXYStruct> _constantSlabXYDimensions;


  float _plugLength;


};

#endif

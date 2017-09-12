#include "SEcal05_Helpers.h"

SEcal05_Helpers::SEcal05_Helpers() {
  // constructor, initialise internal variables
  _x_det=NULL;
  _layering=NULL;
  _preshower=-1;

  _CF_absWrap=-1;
  _CF_alvWall=-1;
  _CF_front=-1;
  _CF_back=-1;

  _ntowers.clear();
  _towerGap=-999;
  _unitsPerTower=0;
  _unitDeadEdge=-999;

  _cells_across_megatile=0;
  _strips_across_megatile=0;
  _strips_along_megatile=0;

  _constantSlabXYDimensions.clear();

  _caloLayer.absorberThickness = -999;
  _caloLayer.sensitive_thickness = -999;
  _caloLayer.distance = 0;

  _totThick=0;

  _plugLength=0;

  _magicMegatileStrategy=-1;
}


void SEcal05_Helpers::setAbsLayers( int nl1, double th1, int nl2, double th2, int nl3, double th3 ) {
  // set the number and thicknesses of absorber layers
  _nlayers1=nl1;
  _nlayers2=nl2;
  _nlayers3=nl3;
  _radiator_thickness1=th1;
  _radiator_thickness2=th2;
  _radiator_thickness3=th3;
  return;
}

void SEcal05_Helpers::checkLayerConsistency() {
  // check requested number of absorber layers is cinsistent with preshower status 
  // we are constrained to have an even number of sensitive layers [ 2 sens. layers / slab ]
  assert( (_preshower==0 || _preshower==1) && "_preshower not set" );
  int n_total_abs_layers = _nlayers1 + _nlayers2 + _nlayers3;    // total number of Absorber layers
  // check that number of requested absober layers is consistent
  // if we want a preshower layer, total number of Absorber layers should be odd; otherwise even
  if ( _preshower==1 && n_total_abs_layers%2==0 ) {
    std::cout << "SEcal05_Helpers ERROR: inconsistent ECAL model !! if you request a preshower layer, the number of absorber layers = _nlayers1 + _nlayers2 + _nlayers3 must be odd" << std::endl;
    std::cout << " Ecal_PreshowerLayer = " << _preshower << " ; _nlayers1/2/3 = " << _nlayers1 << " " << _nlayers2 << " " << _nlayers3 << std::endl;
    assert(0);
  } else if ( _preshower==0 && n_total_abs_layers%2==1 ) {
    std::cout << "SEcal05_Helpers ERROR: inconsistent ECAL model !! if you request no preshower layer, the number of absorber layers = _nlayers1 + _nlayers2 + _nlayers3 must be even" << std::endl;
    std::cout << " Ecal_PreshowerLayer = " << _preshower << " ; _nlayers1/2/3 = " << _nlayers1 << " " << _nlayers2 << " " << _nlayers3 << std::endl;
    assert(0);
  }
  return;
}

float SEcal05_Helpers::getTotalThickness() {
  // calculate total thickness of ECAL
  checkLayerConsistency();
  assert( _x_det && _layering && "_x_det or _layering not set" );
  assert( (_preshower==0 || _preshower==1) && "_preshower not set" );
  assert( _CF_absWrap>=0. && _CF_alvWall>=0. && _CF_front>=0. && _CF_back>=0. && "CF thicknesses not set" );
  float totalThickness = _CF_front+_CF_back;// front and back supports

  // the absorber in the structure
  for (unsigned int i=0; i<_nlayers1+_nlayers2+_nlayers3; i++) {
    bool inStructure = _preshower ? i%2==1 : i%2==0 ;
    if ( inStructure ) {
      double thickness (_radiator_thickness1);
      if ( i>=_nlayers1 ) thickness = _radiator_thickness2;
      if ( i>=_nlayers1+_nlayers2 ) thickness = _radiator_thickness3;
      totalThickness += thickness + 2*_CF_absWrap; // the absorber and its wrapping
    }
  }
  // the slabs
  int l_num(0);
  for(xml_coll_t li(*_x_det,_U(layer)); li; ++li)  { // types of layers (i.e. thin/thick absorber) or "stack"
    xml_comp_t x_layer = li;
    // Loop over number of repeats for this layer.
    for (int j=0; j< x_layer.repeat(); j++)    {  // layers within this type (or "stack")
      float thisthick = _layering->layer(l_num)->thickness();
      totalThickness+=thisthick + 2*_CF_alvWall; // slab thickness, and the alveolar wall around it
      l_num++;
    }
  }
  return totalThickness;
}

void SEcal05_Helpers::printSEcal05LayerInfo( dd4hep::rec::LayeredCalorimeterData::Layer & caloLayer) {
  std::cout<<"SEcal05_Helpers === CALOLAYER printout: "                 << std::endl;
  std::cout<<"    caloLayer.distance: "                 << caloLayer.distance <<std::endl;
  std::cout<<"    caloLayer.inner_nRadiationLengths: "  << caloLayer.inner_nRadiationLengths <<std::endl;
  std::cout<<"    caloLayer.inner_nInteractionLengths: "<< caloLayer.inner_nInteractionLengths <<std::endl;
  std::cout<<"    caloLayer.inner_thickness: "          << caloLayer.inner_thickness <<std::endl;
  std::cout<<"    caloLayer.sensitive_thickness: "      << caloLayer.sensitive_thickness <<std::endl;
  std::cout<<"    caloLayer.cellSize0, 1: "             << caloLayer.cellSize0 << " " << caloLayer.cellSize1 << std::endl;
  std::cout<<"    caloLayer.outer_nRadiationLengths: "  << caloLayer.outer_nRadiationLengths <<std::endl;
  std::cout<<"    caloLayer.outer_nInteractionLengths: "<< caloLayer.outer_nInteractionLengths <<std::endl;
  std::cout<<"    caloLayer.outer_thickness: "          << caloLayer.outer_thickness <<std::endl;
  return;
}

double SEcal05_Helpers::getAbsThickness( unsigned int iAbsLay ) { //, int n1, int n2, int n3, double t1, double t2, double t3) {
  // get the thickness of a given absorber layer
  if ( iAbsLay < _nlayers1 ) return _radiator_thickness1;
  else if ( iAbsLay - _nlayers1 < _nlayers2 ) return _radiator_thickness2;
  else if ( iAbsLay - _nlayers1 - _nlayers2 < _nlayers3 ) return _radiator_thickness3;
  assert(0 && "impossible layer number"); // should never get here
  return -999.;
}



std::vector <SEcal05_Helpers::dimposXYStruct> SEcal05_Helpers::getAbsPlateXYDimensions( double ztop ) {
  // get the dimensions of absorber plates
  // for trapeziodal case, this depends on the Z position
  // CF wrap not taken into account, so this width constains both absorber and its wrapping

  if ( _module_XZtype==1 ) assert( ztop>=0 && "getAbsPlateXYDimensions: ztop not specified" ); // must specify Z position for case with non-uniform layers

  std::vector <dimposXYStruct> absorbersheets;

  if ( _module_XYtype == 1 ) { // each slab is different: calculate slab-by-slab (eg in endcap)
    absorbersheets = getSlabXYDimensions( ztop );
  } else {                    // single sheet over all slabs (eg in barrel)
    dimposXYStruct ss;
    if ( _module_XZtype==0 ) { // no taper
      ss.sizeX = _module_dX_max;
      ss.posX = ss.sizeX/2.;
    } else {   // size varies by laer

      //      ss.sizeX = _module_dX_max - 2.*ztop; // this assumes octagon

      ss.sizeX = _module_dX_max - 2.*ztop/tan(_module_angle); // for general shape

      ss.posX  = _module_dX_max/2.;   // keep it centered
    }
    ss.sizeY = _module_dY_total;
    ss.posY = _module_dY_total/2.;
    absorbersheets.push_back(ss);
  }
  return absorbersheets;
}

std::vector <SEcal05_Helpers::dimposXYStruct> SEcal05_Helpers::getSlabXYDimensions( double ztop ) {
  // calculate size and position of slab volumes
  //  size includes slab and dead area around it

  // ztop is Z of upper face of this slab (the face shortest in X)
  if ( _module_XZtype==1 ) assert( ztop>=0 && "getSlabXYDimensions: ztop not specified" ); // must specify Z position for case with non-uniform layers

  std::vector <dimposXYStruct> layerslabdims;
  if ( _module_XZtype == 0 && _constantSlabXYDimensions.size()>0 ) { // every layer is the same, and has already been calculated
    layerslabdims = _constantSlabXYDimensions;
  } else {
    dimposXYStruct ss;
    int totTow(0);
    for (size_t imod=0; imod<_ntowers.size(); imod++) { // loop over "modules" [nb, for barrel, we only make a single module; for endcaps, typically have 3 per quadrant]
      for (int itow = 0; itow<_ntowers[imod]; itow++) { // towers within the module
	// each slab has same width
        ss.sizeY = _alveolus_total_dim_Y; // this size include edge-of-tower dead space
        ss.posY = 
	  _moduleGap * ( 1 + 2*imod ) + // edge-of-module gaps from Y=0 to this slab
	  _alveolus_total_dim_Y * (totTow+0.5); // position of slab centre w.r.t (X,Y) = (0,0)

        if ( _module_XYtype==0 ) { // all slabs in a layer have same length

          if ( _module_XZtype == 0 ) { // all layers have same length
            ss.sizeX = _module_dX_max;
            ss.posX = ss.sizeX/2.;
          } else {                     // layers vary in length : barrel module

            // ss.sizeX = _module_dX_max - 2.*ztop; // assumes octagon

	    ss.sizeX = _module_dX_max - 2.*ztop/tan(_module_angle); // for general shape

            ss.posX = _module_dX_max/2.;  // slabs all centred
          }

	  ss.sizeX -= _plugLength; // DANIELHACK 
	  ss.posX += _plugLength/2.;  // DANIELHACK 



        } else if ( _module_XYtype==1 ) { // slabs within layer have different lengths : this is for endcap
          // placing in Y
          double upperY = ss.posY + 0.5*_alveolus_total_dim_Y; // Y position of upper edge

	  if ( upperY<=_module_dY_kink ) { // straight edge part under kink
	    ss.sizeX = _module_dX_max;
	  } else {            // sloping edge: take length at upperY, which is the shortest one
	    ss.sizeX = _module_dX_max - ( upperY - _module_dY_kink );
	  }

          ss.posX = ss.sizeX/2.; // this aligns the -X end of slab

	  ss.sizeX -=  _plugLength;  // DANIELHACK
	  ss.posX  +=  _plugLength/2.;  // DANIELHACK      

        } else {
	  cout << " SEcal05_Helpers ERROR _module_XYtype = " << _module_XYtype << "!!!" << endl;
	  assert(0);
	}
	layerslabdims.push_back(ss);
	totTow++;
      } // towers
    } // modules

    // if slab dimensions do not change layer-by-layer, memorise for next time
    if ( _module_XZtype == 0 ) _constantSlabXYDimensions = layerslabdims;
  }
  return layerslabdims;
}


void SEcal05_Helpers::updateCaloLayers(double thickness,
                                       dd4hep::Material mat,
                                       bool isAbsorber,
                                       bool isSensitive,
                                       double cell_size_x, double cell_size_y,
                                       bool isFinal
                                       ) {

  if ( isFinal ) { // add material before saving layer
    _layer_thickness           += (thickness);
    _layer_nRadiationLengths   += (thickness)/mat.radLength() ;
    _layer_nInteractionLengths += (thickness)/mat.intLength() ;
  }

  // should we finish off the current layer?
  if ( isFinal ||  // last slice of calo
       (isAbsorber && _caloLayer.sensitive_thickness > 0 )  // end of a caloLayer
       ) {

    // finalise caloLayer entry, add to caloData
    _caloLayer.outer_thickness           = _layer_thickness;
    _caloLayer.outer_nRadiationLengths   = _layer_nRadiationLengths;
    _caloLayer.outer_nInteractionLengths = _layer_nInteractionLengths;
    //_caloLayer.thickness                 = _caloLayer.inner_thickness + _caloLayer.outer_thickness ;

    // push it back
    _caloData->layers.push_back( _caloLayer ) ;
    //    printSEcal05LayerInfo( _caloLayer );

    // reset layer thicknesses
    _layer_thickness=0;
    _layer_nRadiationLengths=0;
    _layer_nInteractionLengths=0;

    // update the calolayer.distance ( DJeans 12 sep 2017)
    // here this is distance from ECAL start to the layer start; the Barrel and Endcap drivers add the distance from IP
    _caloLayer.distance += _caloLayer.inner_thickness+_caloLayer.outer_thickness;


  }

  if (!isFinal) {

    // first add half the material
    _layer_thickness           += (thickness/2.);
    _layer_nRadiationLengths   += (thickness/2.)/mat.radLength() ;
    _layer_nInteractionLengths += (thickness/2.)/mat.intLength() ;

    // then we store material before and after centre of this layer
    if ( isSensitive ) {
      _caloLayer.cellSize0 = cell_size_x;
      _caloLayer.cellSize1 = cell_size_y;
      _caloLayer.sensitive_thickness       = thickness ;
      _caloLayer.inner_nRadiationLengths   = _layer_nRadiationLengths ;
      _caloLayer.inner_nInteractionLengths = _layer_nInteractionLengths ;
      _caloLayer.inner_thickness           = _layer_thickness ;

      // reset layer thicknesses
      _layer_thickness=0;
      _layer_nRadiationLengths=0;
      _layer_nInteractionLengths=0;
    }

    // add in remaining half of layer
    _layer_thickness           += (thickness/2.);
    _layer_nRadiationLengths   += (thickness/2.)/mat.radLength() ;
    _layer_nInteractionLengths += (thickness/2.)/mat.intLength() ;
  }

  _totThick+=thickness;

  return;
}


SEcal05_Helpers::dxinfo SEcal05_Helpers::getNormalMagicUnitsInX( double dx_total, double dx_unit, double dx_cell, double dx_dead ,
								 int magicStrategy ) {

  dxinfo dxInf;
  dxInf.normal_nX=-1;
  dxInf.magic1_unitDX=-1;
  dxInf.magic1_ncellsX=-1;
  dxInf.magic2_unitDX=-1;

  int nNormalUnit = int( floor( dx_total / dx_unit ) );

  if ( magicStrategy==0 ) { // just integer number of standard units (wafers/megatiles)
    dxInf.normal_nX = nNormalUnit;
  } else {

    double extraSpace = dx_total - nNormalUnit*dx_unit;
    double extraSensitiveSpace = extraSpace - 2.*dx_dead;
    double extraNCells = extraSensitiveSpace / dx_cell;
 
    if ( magicStrategy==1 ) { // last magic megatile has integer number of standard-sized cells
      dxInf.normal_nX = nNormalUnit;
      if ( extraNCells > 1.0 ) { 
	int iext = int( floor( extraNCells ) );
	dxInf.magic1_unitDX = iext*dx_cell + 2.*dx_dead;
	dxInf.magic1_ncellsX = iext;
      }
    } else if ( magicStrategy==2 ) { 

      // one magic megatile with integer number of standard cells, 
      //  second with non-standard cell size to fill space exactly
      //   last magic cell size in range  shortestMacigCell<magicCellsSize/standardCellSize<shortestMagicCell

      // the shortest magic cell (in terms of the usual cell size)
      const double shortestMagicCell = 1.0; // this means that the last magic cell is at least as long as the standard cell

      extraSensitiveSpace = extraSpace - 4.*dx_dead; // because both the magic tiles may have a dead area
      extraNCells = extraSensitiveSpace / dx_cell;

      if ( extraNCells < shortestMagicCell ) { // too small to make magic unit. remove one normal unit
	nNormalUnit-=1;
	extraSpace = dx_total - nNormalUnit*dx_unit;
	extraSensitiveSpace = extraSpace - 4.*dx_dead;
	extraNCells = extraSensitiveSpace / dx_cell;
      }
      
      int imagic1 = int( floor( extraNCells ) ); // rounded down number of standard-sized cells

      double magic1size = imagic1>0 ? imagic1*dx_cell + 2*dx_dead : 0.;

      double magic2cellsize = extraSpace - magic1size - 2*dx_dead;

      if ( magic2cellsize/dx_cell < shortestMagicCell ) { // too small
	imagic1 -= 1; // remove last standard cell
	magic1size = imagic1>0 ? imagic1*dx_cell + 2*dx_dead : 0;
	magic2cellsize = extraSpace - magic1size - 2*dx_dead;
      }

      assert( magic2cellsize/dx_cell >= shortestMagicCell && magic2cellsize/dx_cell <= shortestMagicCell+1.0 && "problem in deciding magic unit size" );

      dxInf.normal_nX = nNormalUnit;
      dxInf.magic1_ncellsX = imagic1;
      dxInf.magic1_unitDX = imagic1>0 ? imagic1*dx_cell + 2*dx_dead : 0;

      dxInf.magic2_unitDX = magic2cellsize + 2*dx_dead;
      
    }
  }

  if ( magicStrategy==2 ) {
    // check it fills exactly
    if ( fabs( dxInf.normal_nX*dxInf.magic1_unitDX + dxInf.magic1_unitDX + dxInf.magic2_unitDX - dx_total ) < 0.01*dd4hep::mm ) {
      cout << " SEcal05_Helpers ERROR : " << endl;
      cout << dxInf.normal_nX*dxInf.magic1_unitDX << " " << dxInf.magic1_unitDX << " " << dxInf.magic2_unitDX << endl;
      cout << dxInf.normal_nX*dxInf.magic1_unitDX + dxInf.magic1_unitDX + dxInf.magic2_unitDX << " " << dx_total << endl;
      cout << dxInf.normal_nX*dxInf.magic1_unitDX + dxInf.magic1_unitDX + dxInf.magic2_unitDX - dx_total << endl;
      assert(0 && "magic unit does not fill space exactly");
    }

  }

  return dxInf;
}


void SEcal05_Helpers::makeModule( dd4hep::Volume & mod_vol,  // the volume we'll fill
				  dd4hep::DetElement & stave_det, // the detector element
				  dd4hep::rec::LayeredCalorimeterData & caloData, // the reco data we'll fill
				  dd4hep::Detector & theDetector,
				  dd4hep::SensitiveDetector & sens
				  ) {
  // make the module
  _caloData = &caloData;

  // calculate widths of module/tower/unit (in z-direction for barrel: across the slab/module
  assert ( _ntowers.size()>0 && _unitsPerTower>0 && "_ntowers or _unitsPerTower not set" );

  _module_dY_total=0;
  for (size_t i=0; i<_ntowers.size(); i++) {
    _module_dY_total = _ntowers[i]*_alveolus_total_dim_Y + 2.*_moduleGap;
  }

  double alveolus_active_dim_Y = _alveolus_total_dim_Y - 2.*_towerGap;
  double unit_dim_Y            = alveolus_active_dim_Y/_unitsPerTower;
  double unit_sensitive_dim_Y  = unit_dim_Y - 2*_unitDeadEdge; // width of the sensitive area in each sensor
  assert( unit_sensitive_dim_Y>0 && "negative-sized sensitive area..." ); // otherwise really weird!

  // get detector stuff
  assert ( _x_det && "_x_det not set");
  int det_id           = _x_det->id();
  xml_comp_t x_staves  = _x_det->staves();

  _carbon_fibre_material = theDetector.material("CarbonFiber");
  _radiator_material     = theDetector.material(x_staves.materialStr());
  _air_material          = theDetector.air();

  dd4hep::VisAttr _radiator_visatt      = theDetector.visAttributes( x_staves.visStr() );

  // get segmentation stuff
  assert (_geomseg && "segmentation not set");

  dd4hep::DDSegmentation::WaferGridXY* waferSeg = dynamic_cast< dd4hep::DDSegmentation::WaferGridXY* > ( _geomseg->segmentation() ) ;


  dd4hep::DDSegmentation::MegatileLayerGridXY* megatileSeg = dynamic_cast< dd4hep::DDSegmentation::MegatileLayerGridXY* > ( _geomseg->segmentation() ) ;
  assert( (waferSeg || megatileSeg) && "no segmentation found" );

  // set up the standard megatile size and offset
  if ( megatileSeg ) {
    megatileSeg->setMegaTileSizeXY( unit_sensitive_dim_Y, unit_sensitive_dim_Y );
    megatileSeg->setMegaTileOffsetXY( -unit_sensitive_dim_Y/2., -unit_sensitive_dim_Y/2. );
  }

  _module_thickness = getTotalThickness();

  double currentLayerBase_pos_Z  = 0; // this gets updated: it's the position of the bottom face of the next detector slice

  int layer_index = 0; // 1;// layer number - Daniel changes to c-type counting from 0...lets get rid of fortran habits!

  // to deal with initial layers
  bool isFrontFace = true;
  int myLayerNum = 0 ; // this is the sensitive layer number (one per sensitive)
  unsigned int absorber_index(0); // keep track of which absorber layer we're in

  // keep track of thickness within a layer
  // initialise
  _layer_thickness           = 0. ;
  _layer_nRadiationLengths   = 0. ;
  _layer_nInteractionLengths = 0. ;

  //--------------------------------
  // loop over the layers
  //---------------------------------

  for(xml_coll_t li(*_x_det,_U(layer)); li; ++li)  { // types of layers (i.e. thin/thick absorber) or "stack"
    xml_comp_t x_layer = li;
    //    cout << " ---- NEW LAYER TYPE " << layer_index << " repeat " <<  x_layer.repeat() << endl;

    // Loop over number of repeats for this layer type
    for (int j=0; j< x_layer.repeat(); j++)    {  // layers within this type (or "stack")
      std::string l_name = dd4hep::_toString(layer_index,"layer%d");

      // #########################
      // Build Structure Layer
      // #########################

      //----------------------------------------------------
      // position and thickness of absorber plates in the structure
      //----------------------------------------------------
      double radiator_dim_Z(0);
      double rad_pos_Z(0);
      double this_struct_CFthick_beforeAbs(0);
      double this_struct_CFthick_afterAbs(0);

      if ( isFrontFace ) { // the first part of the module depends on whether we have preshwer or not
        if ( _preshower==1 ) { // don't include W+CF wrapping; only one side of alveolus
          this_struct_CFthick_beforeAbs = 0;
          this_struct_CFthick_afterAbs = _CF_front + _CF_alvWall;
          radiator_dim_Z = 0;
        } else { // include W+CF wrapping; only one side of alveolus
          // cout << " -- no preshower, including absorber" << endl;
          this_struct_CFthick_beforeAbs = _CF_front + _CF_absWrap; // allow for initial CF front plate also in no preshower case - djeans 6 july 2017
          this_struct_CFthick_afterAbs = _CF_absWrap + _CF_alvWall;
          radiator_dim_Z = getAbsThickness( absorber_index++ );
          rad_pos_Z = radiator_dim_Z/2. + this_struct_CFthick_beforeAbs; // distance from top surface of structure to centre of radiator
        }
        isFrontFace=false;
      } else { // internal layer: include W+CF wrapping; both sides of alveolus
        this_struct_CFthick_beforeAbs = _CF_alvWall+_CF_absWrap;
        this_struct_CFthick_afterAbs = _CF_alvWall+_CF_absWrap;
        radiator_dim_Z = getAbsThickness( absorber_index++ );
        assert( radiator_dim_Z>0 && "no radiator!" );
        rad_pos_Z = this_struct_CFthick_beforeAbs + radiator_dim_Z/2.; // distance from top surface of structure to centre of radiator
      }


      // create the radiator volume
      if ( radiator_dim_Z>0 ) {       // only create the volume if we have radiator
        std::vector < dimposXYStruct > absorbersheets = getAbsPlateXYDimensions( currentLayerBase_pos_Z + this_struct_CFthick_beforeAbs + radiator_dim_Z ); // add CF before abs. added djeans 21 nov 2016
        // create and place the absorber sheets
        for ( size_t ipl=0; ipl<absorbersheets.size(); ipl++) {
          dimposXYStruct plSize = absorbersheets[ipl];
          dd4hep::Box    barrelStructureLayer_box( plSize.sizeX/2.,
                                                             plSize.sizeY/2. - _CF_absWrap,  // remove CF wrapping
                                                             radiator_dim_Z/2.);

          dd4hep::Volume barrelStructureLayer_vol( _det_name+"_"+l_name+"_"+dd4hep::_toString(int(ipl),"bs%02d"),
                                                             barrelStructureLayer_box, _radiator_material);

          barrelStructureLayer_vol.setVisAttributes( _radiator_visatt );

          // Position the layer.
          dd4hep::Position      bsl_pos = getTranslatedPosition(plSize.posX, plSize.posY, currentLayerBase_pos_Z + rad_pos_Z );

	  //          dd4hep::PlacedVolume  barrelStructureLayer_phv = 
	  mod_vol.placeVolume(barrelStructureLayer_vol, bsl_pos);

        } // loop over sheets
      }

      // update layer thickness until front face of absorber
      updateCaloLayers( this_struct_CFthick_beforeAbs, _carbon_fibre_material, false, false);
      updateCaloLayers( radiator_dim_Z, _radiator_material, true, false );
      updateCaloLayers( this_struct_CFthick_afterAbs, _carbon_fibre_material, false, false);

      // update position within module
      currentLayerBase_pos_Z += this_struct_CFthick_beforeAbs + radiator_dim_Z + this_struct_CFthick_afterAbs;

      // #########################
      // Build Slab
      // #########################

      //-------------------------------
      // first sum the various materials for the caloData/caloLayer
      //-------------------------------
      int myLayerNumTemp = myLayerNum;
      for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
        xml_comp_t  x_slice = si;
        double      s_thick = x_slice.thickness();
        dd4hep::Material slice_material  = theDetector.material( x_slice.materialStr() );
        if (x_slice.materialStr().compare(x_staves.materialStr()) == 0){
          // this is absorber material
          //  check it's consistent with what we expect from the detector parameters
          assert ( fabs( s_thick - getAbsThickness( absorber_index++ ) ) < 1e-5 && "inconsistent radiator thickness" );
          updateCaloLayers( s_thick, slice_material, true, false ); // absorber
        } else if ( x_slice.isSensitive() ) {
          double cell_size_x(0), cell_size_y(0);
          if ( waferSeg ) {
            cell_size_x = waferSeg->cellDimensions(0)[0];
            cell_size_y = waferSeg->cellDimensions(0)[1];
          } else if ( megatileSeg ) {
            // setup megatile
            int laytype = _layerConfig [ myLayerNumTemp%_layerConfig.size() ];
            if ( laytype==0 ) {
              megatileSeg->setMegaTileCellsXY( myLayerNumTemp, _cells_across_megatile , _cells_across_megatile );
            } else if  ( laytype==1 ) {
              megatileSeg->setMegaTileCellsXY( myLayerNumTemp, _strips_across_megatile , _strips_along_megatile ); // strips in one orientation
            } else if  ( laytype==2 ) {
              megatileSeg->setMegaTileCellsXY( myLayerNumTemp, _strips_along_megatile , _strips_across_megatile ); // and in the other
            } else {
              assert(0 && "unknown layer type");
            }
            cell_size_x = megatileSeg->cellDimensions(myLayerNumTemp, 0)[0]; // dummy wafer
            cell_size_y = megatileSeg->cellDimensions(myLayerNumTemp, 0)[1];
          }
          updateCaloLayers( s_thick, slice_material, false, true, cell_size_x, cell_size_y ); // sensitive
          myLayerNumTemp++;
        } else {
          updateCaloLayers( s_thick, slice_material, false, false );
        }
      }

      //------------------------------------
      // then actually construct the slabs
      //------------------------------------
      double slab_dim_Z = _layering->layer(layer_index)->thickness();
      double slab_pos_Z = currentLayerBase_pos_Z + slab_dim_Z/2.; // centre position of slab

      std::vector <dimposXYStruct> slabDims = getSlabXYDimensions( currentLayerBase_pos_Z + slab_dim_Z );

      for (size_t islab = 0; islab<slabDims.size(); islab++) {
        myLayerNumTemp = myLayerNum;
        double slab_dim_X  = slabDims[islab].sizeX;

        // make an air volume for the alveolus
        dd4hep::Box        l_box( slab_dim_X/2. ,
                                            slabDims[islab].sizeY/2. - _CF_alvWall,
                                            slab_dim_Z/2. );

        dd4hep::Volume     l_vol( _det_name+"_alveolus_"+l_name, l_box, _air_material);
	l_vol.setVisAttributes(theDetector.visAttributes( "GrayVis" ) );

        dd4hep::DetElement l_det( stave_det, l_name+dd4hep::_toString(int(islab),"tower%02d") , det_id );
        dd4hep::Position   l_pos = getTranslatedPosition(slabDims[islab].posX, slabDims[islab].posY, slab_pos_Z );
        dd4hep::PlacedVolume l_phv = mod_vol.placeVolume(l_vol,l_pos);
        l_phv.addPhysVolID("tower", int(islab) );
        l_det.setPlacement(l_phv);

        // then fill it with the slab sublayers
        int s_num(0);
        double s_pos_Z = -slab_dim_Z / 2.; // position of sub-layer with respect to centre of this slab

        for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  { // the sub-layers
          xml_comp_t  x_slice = si;
          std::string s_name  = dd4hep::_toString(s_num,"slice%d");
          double      s_thick = x_slice.thickness();

          dd4hep::Material slice_material  = theDetector.material( x_slice.materialStr() );

	  std::string vis_str = x_slice.visStr();

	  if ( !x_slice.isSensitive() ) { // not the sensitive slice: just a layer of stuff

	    dd4hep::Box      s_box( slab_dim_X/2. , slabDims[islab].sizeY/2. - _CF_alvWall, s_thick/2. );
	    dd4hep::Volume   s_vol(_det_name+"_"+l_name+"_"+s_name, s_box, slice_material);
	    s_vol.setVisAttributes(theDetector.visAttributes( vis_str ));

	    dd4hep::Position s_pos( 0, 0, s_pos_Z + s_thick/2. );
	    //	    dd4hep::PlacedVolume slice_phv = 
	    l_vol.placeVolume(s_vol, s_pos );

	  } else { // sensitive slice

            // Normal squared wafers - this is just the sensitive part
            // square piece of silicon, not including guard ring. guard ring material is not included

            dd4hep::Box WaferSiSolid( unit_sensitive_dim_Y/2., unit_sensitive_dim_Y/2., s_thick/2.);

	    // get the standard cell size in X for this layer
	    double cell_size_x = waferSeg ? waferSeg->cellDimensions(0)[0] : megatileSeg->cellDimensions(myLayerNumTemp, 0)[0];
	    double cell_size_y = waferSeg ? waferSeg->cellDimensions(0)[1] : megatileSeg->cellDimensions(myLayerNumTemp, 0)[1];

	    // work out how to make the magic units, if requested
	    dxinfo xseg = getNormalMagicUnitsInX( slab_dim_X,
						  unit_dim_Y,
						  cell_size_x,
						  _unitDeadEdge,
						  _magicMegatileStrategy );						  

            int n_wafers_x = xseg.normal_nX;
	    int ncellsy = int( unit_sensitive_dim_Y/cell_size_y + 0.5 ); // should be exactly integer. take nearest integer to account for possible rounding errors.

            int wafer_num = 0;

	    double wafer_pos_X = -slab_dim_X/2.; // keep track of wafer position in X

            for (int n_wafer_x = 0; n_wafer_x < n_wafers_x+2 ; n_wafer_x++) { // loop along slab, including magic megatile/wafer

              double megatile_size_x = unit_dim_Y;
	      int ncellsx(-1);

	      bool isMagic = n_wafer_x >= n_wafers_x;

	      if ( n_wafer_x==n_wafers_x ) { // first magic unit
		megatile_size_x = xseg.magic1_unitDX;
		ncellsx = xseg.magic1_ncellsX;
	      } else if ( n_wafer_x==n_wafers_x+1 ) { // second magic unit
		megatile_size_x = xseg.magic2_unitDX;
		ncellsx = 1;
	      }

	      // cout << "ismagic " << isMagic << " wafer " << n_wafer_x << " / " << n_wafers_x << " : sizex " << megatile_size_x << " nx " << ncellsx << endl;

	      if ( megatile_size_x<=0 ) continue;

	      // cout << " PASSED" << endl;

              double megatile_sensitive_size_x = megatile_size_x - 2*_unitDeadEdge;

              for (int n_wafer_Y = 0; n_wafer_Y < _unitsPerTower; n_wafer_Y++) { // loop across slab (usually 2 wafers, or 1 EBU) [ie along beam dir for barrel module]

                double wafer_pos_Y = -alveolus_active_dim_Y/2.0 + (n_wafer_Y+0.5)*unit_dim_Y;
                wafer_num++;
                std::string Wafer_name;
                if ( isMagic ) Wafer_name="magic";
                Wafer_name +=  dd4hep::_toString(wafer_num,"wafer%d");

                dd4hep::Box* box = isMagic ? new dd4hep::Box( megatile_sensitive_size_x/2,unit_sensitive_dim_Y/2,s_thick/2.) : &WaferSiSolid;
                dd4hep::Volume WaferSiLog(_det_name+"_"+l_name+"_"+s_name+"_"+Wafer_name,*box,slice_material);

		std::string wafer_vis_str = isMagic ? "YellowVis" : vis_str;
		
		WaferSiLog.setVisAttributes(theDetector.visAttributes( wafer_vis_str ));	  
                WaferSiLog.setSensitiveDetector(sens);

                dd4hep::Position w_pos(wafer_pos_X + megatile_size_x/2., wafer_pos_Y, s_pos_Z + s_thick/2. );
                dd4hep::PlacedVolume wafer_phv = l_vol.placeVolume(WaferSiLog, w_pos );
                wafer_phv.addPhysVolID("wafer", wafer_num);
		wafer_phv.addPhysVolID("layer", myLayerNumTemp );

                if ( isMagic ) {
                  if ( megatileSeg ) { // define the special megatile
                    megatileSeg->setSpecialMegaTile( myLayerNumTemp, wafer_num,
                                                     megatile_sensitive_size_x, unit_sensitive_dim_Y,
                                                     -megatile_size_x/2., -unit_sensitive_dim_Y/2., // the offset
                                                     ncellsx, ncellsy ); // the segmentation
                  }
                } else { // not magic
                  if ( waferSeg ) {  // Normal squared wafers, this waferOffsetX is 0.0 // if its an odd number of cells, need to do something?
                    waferSeg->setWaferOffsetX(myLayerNumTemp, wafer_num, 0.0);
                  }
                } // isMagic

              } // y-wafers

	      wafer_pos_X += megatile_size_x;

            } // x-wafers
            myLayerNumTemp++;
          } // end sensitive slice

          // Increment Z position of slice.
          s_pos_Z += s_thick;

          // Increment slice number.
          ++s_num;

        } // end slice within one slab
      } // slabs
      myLayerNum = myLayerNumTemp;
      currentLayerBase_pos_Z += slab_dim_Z;
      layer_index++; // this is the structural layer counter (one per slab, not per sensitive layer)
    } // layers in compact file
  } // layer types in compact file

    // add material after last slab. Just CF
  updateCaloLayers( _CF_alvWall + _CF_back, _carbon_fibre_material, false, false, -1, -1, true ); // the last layer

  return;
}

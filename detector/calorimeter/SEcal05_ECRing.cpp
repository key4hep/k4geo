//====================================================================
//  lcgeo - LC detector models in DD4hep
//--------------------------------------------------------------------
//  DD4hep Geometry driver for SiWEcalEndcaps
//  Ported from Mokka
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id$
//====================================================================

/* History:

// *******************************************************
// *                                                     *
// *                      Mokka                          *
// *   - the detailed geant4 simulation for ILC   -      *
// *                                                     *
// * For more information about Mokka, visit the         *
// *                                                     *
// *  Mokka.in2p3.fr  Mokka home page.                   *
// *                                                     *
// *******************************************************
//
// $Id$
// $Name: mokka-07-00 $
//
//
//
// SEcal04.cc
//
Shaojun Lu:  Ported from Mokka SEcal04 Endcaps part. Read the constants from XML
instead of the DB. Then build the Endcap in the same way with DD4hep
construct.
Inside SEcal04, some parameters, which used by Ecal Endcaps, come from
Ecal Barrel. They can be seen here again.
Start ECRing ...
*/


/*

  ==============
  SEcal05_ECRing
  ===============
  SEcal04_ECRing driver originally assumend 2 identical modules in +/- z
  this is not true, because lumical (and therefore hole in EC plug) is not centred, but offset in +z position

  SEcal05_ECRing driver corrects this (applies offset according to the crossing algle, to center it on the outgoing beampipe).
  it also deals with the presence or not of a preshower layer

  this is a hack of the 04 version, not a clean rewrite as was done for barrel and endcap, so it's not very clean & tidy....but I think it works

  DJeans jan 2017

  ========


  - fixed dd4hep::rec::LayeredCalorimeterData for case with no preshower
  - small fixes to layer thicknesses (now take from compact description) to get exactly same structure as main endcaps

  DJeans July 2017


  ======

  - fixes to LayeredCalorimeterData: "distance" defined to start of layer, some other issues
  - fix to use correct absorber thickness at transition between stacks
  - LayeredCalorimeterData material describtion: overall structure made in CF (previously air)

  DJeans Sep 2017

*/


#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace std;

using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationZYX;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Tube;
using dd4hep::Volume;
using dd4hep::_toString;

using dd4hep::rec::LayeredCalorimeterData;

//#define VERBOSE 1

// workaround for DD4hep v00-14 (and older)
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0
#endif

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {
  static double tolerance = 0e0;

  cout << "---------------------------------" << endl;
  cout << " creating Ecal ECRing ( SEcal05_ECRing ) " << endl;
  cout << "---------------------------------" << endl;

  xml_det_t     x_det     = element;
  string        det_name  = x_det.nameStr();
  Layering      layering (element);

  Material      air       = theDetector.air();
  Material      CF        = theDetector.material("CarbonFiber");

  int           det_id    = x_det.id();
  xml_comp_t    x_staves  = x_det.staves();
  DetElement    sdet      (det_name,det_id);

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , sdet ) ;

  dd4hep::xml::setDetectorTypeFlag( element, sdet ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  sens.setType("calorimeter");

  Material stave_material  = theDetector.material(x_staves.materialStr());

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();

  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0); //Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];

  //====================================================================
  //
  // Read all the constant from ILD_o1_v05.xml
  // Use them to build Ecal ECRing
  //
  //====================================================================

  //  read parametere from compact.xml file
  double Ecal_fiber_thickness_alveolus      = theDetector.constant<double>("Ecal_fiber_thickness_alveolus");
  double Ecal_fiber_thickness_structure     = theDetector.constant<double>("Ecal_fiber_thickness_structure");

  double Ecal_radiator_thickness1           = theDetector.constant<double>("Ecal_radiator_layers_set1_thickness");
  double Ecal_radiator_thickness2           = theDetector.constant<double>("Ecal_radiator_layers_set2_thickness");
  double Ecal_radiator_thickness3           = theDetector.constant<double>("Ecal_radiator_layers_set3_thickness");
  double Ecal_Barrel_halfZ                  = theDetector.constant<double>("Ecal_Barrel_halfZ");

  double Ecal_support_thickness             = theDetector.constant<double>("Ecal_support_thickness");
  double Ecal_front_face_thickness          = theDetector.constant<double>("Ecal_front_face_thickness");
  double Ecal_lateral_face_thickness        = theDetector.constant<double>("Ecal_lateral_face_thickness");

  double Ecal_EC_Ring_gap                   = theDetector.constant<double>("Ecal_EC_Ring_gap");

  double Ecal_cables_gap                    = theDetector.constant<double>("Ecal_cables_gap");
  double Lcal_outer_radius                  = theDetector.constant<double>("Lcal_outer_radius");
  double Ecal_Lcal_ring_gap                 = theDetector.constant<double>("Ecal_Lcal_ring_gap");
  double Ecal_endcap_center_box_size        = theDetector.constant<double>("Ecal_endcap_center_box_size");

  int    Ecal_nlayers1                      = theDetector.constant<int>("Ecal_nlayers1");
  int    Ecal_nlayers2                      = theDetector.constant<int>("Ecal_nlayers2");
  int    Ecal_nlayers3                      = theDetector.constant<int>("Ecal_nlayers3");

  double      EcalEndcapRing_inner_radius   = theDetector.constant<double>("EcalEndcapRing_inner_radius");
  double      EcalEndcapRing_outer_radius   = theDetector.constant<double>("EcalEndcapRing_outer_radius");
  double      EcalEndcapRing_min_z          = theDetector.constant<double>("EcalEndcapRing_min_z");
  double      EcalEndcapRing_max_z          = theDetector.constant<double>("EcalEndcapRing_max_z");


  double      crossing_angle = theDetector.constant<double>("ILC_Main_Crossing_Angle");
  bool        Ecal_Barrel_PreshowerLayer         = theDetector.constant<int>("Ecal_Barrel_Preshower") > 0;

  double      Ecal_barrel_thickness         = theDetector.constant<double>("Ecal_barrel_thickness");

  //========== fill data for reconstruction ============================
  dd4hep::rec::LayeredCalorimeterData* caloData = new dd4hep::rec::LayeredCalorimeterData ;
  caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = 0  ; // hard code cernter pipe hole
  caloData->outer_symmetry = 4  ; // outer box
  caloData->phi0 = 0 ;

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = EcalEndcapRing_inner_radius ;
  caloData->extent[1] = EcalEndcapRing_outer_radius ;
  caloData->extent[2] = EcalEndcapRing_min_z ;
  caloData->extent[3] = EcalEndcapRing_max_z ;


  //====================================================================
  //
  // general calculated parameters
  //
  //====================================================================

  int Number_of_Si_Layers_in_Barrel = 0;

#ifdef VERBOSE
  std::cout << " Ecal total number of Silicon layers = " << Number_of_Si_Layers_in_Barrel  << std::endl;
#endif

  int n_total_layers = Ecal_nlayers1 + Ecal_nlayers2 + Ecal_nlayers3;
  Number_of_Si_Layers_in_Barrel = n_total_layers; // take account of preshower on/off (djeans)
  if ( Ecal_Barrel_PreshowerLayer ) Number_of_Si_Layers_in_Barrel++;

  // calculate total thickness of ECAL
  // updated to take info from xml description - djeans: 12 July 2017
  double module_thickness = Ecal_support_thickness + Ecal_front_face_thickness;// front and back supports

  // the absorber in the structure
  for (int i=0; i<Ecal_nlayers1+Ecal_nlayers2+Ecal_nlayers3; i++) {
    bool inStructure = Ecal_Barrel_PreshowerLayer==1 ? i%2==1 : i%2==0 ;
    if ( inStructure ) {
      double thickness (Ecal_radiator_thickness1);
      if ( i>=Ecal_nlayers1 ) thickness = Ecal_radiator_thickness2;
      if ( i>=Ecal_nlayers1+Ecal_nlayers2 ) thickness = Ecal_radiator_thickness3;
      module_thickness += thickness + 2*Ecal_fiber_thickness_structure; // the absorber and its wrapping
    }
  }

  // the slabs
  int l_num2(0);
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  { // types of layers (i.e. thin/thick absorber) or "stack"
    xml_comp_t x_layer = li;
    // Loop over number of repeats for this layer.
    for (int j=0; j< x_layer.repeat(); j++)    {  // layers within this type (or "stack")
      float thisthick = layering.layer(l_num2)->thickness();
      module_thickness+=thisthick + 2*Ecal_fiber_thickness_alveolus; // slab thickness, and the alveolar wall around it
      l_num2++;
    }
  }


  if ( fabs( Ecal_barrel_thickness - module_thickness) > 0.1 ) {
    cout << "ERROR EC ecal thickness inconsistency: calculated, nominal = " << module_thickness << " " << Ecal_barrel_thickness << endl;
    assert(0);
  }

#ifdef VERBOSE
  std::cout << " module_thickness = " << module_thickness  << std::endl;
#endif


  double ECRingSiplateSize = Ecal_endcap_center_box_size
    - 2 * Ecal_EC_Ring_gap
    - 2 * Ecal_lateral_face_thickness;


  // central hole is not centred on detector axis, but on outgoing beampipe (as is the lumical)
  double hole_x_offset = tan(crossing_angle/2.) * (EcalEndcapRing_min_z + EcalEndcapRing_max_z)/2.;


  // ========= Create Ecal end cap ring   ====================================
  //  It will be the volume for palcing the Ecal endcaps alveolus(i.e. Layers).
  //  And the structure W plate.
  //  Itself will be placed into the world volume.
  // ==========================================================================

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~              EndcapRing                           ~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for(int module_num=0;module_num<2;module_num++) { // move module loop to here (modules not identical)

    DetElement    module_det(_toString(module_num,"module%d"),det_id);

    double thismodule_hole_dx = hole_x_offset;
    if ( module_num==0 ) thismodule_hole_dx*=-1;
    Position hole_position(thismodule_hole_dx, 0, 0);

    double Ecal_endcap_Tube_rmax = Lcal_outer_radius + Ecal_Lcal_ring_gap;

    Tube CenterECTub(0., Ecal_endcap_Tube_rmax, module_thickness);
    Box  CenterECBox((Ecal_endcap_center_box_size/2. - Ecal_EC_Ring_gap),
                     (Ecal_endcap_center_box_size/2. - Ecal_EC_Ring_gap),
                     module_thickness/2.);

    SubtractionSolid ECRingSolid( CenterECBox, CenterECTub, hole_position);

    //  Volume EnvLogECRing("ECRing",ECRingSolid,theDetector.material("g10"));
    Volume EnvLogECRing("ECRing",ECRingSolid,theDetector.material("CarbonFiber")); // DJeans 5-sep-2016

    //==============================================================
    // build the layer and place into the ECRing
    //==============================================================


    //-------------------------------------------------------
    // Radiator and towers placements inside the ECRing module
    //-------------------------------------------------------

    // We count the layers starting from IP and from 1,
    // so odd layers should be inside slabs and
    // even ones on the structure.
    //


    double distance_ip_to_front_face = Ecal_Barrel_halfZ + Ecal_cables_gap;
    double dist_from_front_face_tally(0);

    double z_floor =  - module_thickness/2;

    double l_pos_z = z_floor;

    dd4hep::rec::LayeredCalorimeterData::Layer caloLayer ;
    caloLayer.cellSize0 = cell_sizeX;
    caloLayer.cellSize1 = cell_sizeY;

    //-------------------- start loop over ECAL layers ----------------------
    // Loop over the sets of layer elements in the detector.
    double radiator_dim_y = Ecal_radiator_thickness1; //to be updated with slice radiator thickness

    double nRadiationLengths=0.;
    double nInteractionLengths=0.;
    double thickness_sum=0;

    int l_num = 1;
    bool isFirstSens = true;
    int myLayerNum = 0 ;
    bool isFrontFace=true;
    int absorber_index(0);

    for(xml_coll_t li(x_det,_U(layer)); li; ++li)  { // types of layer
      xml_comp_t x_layer = li;
      int repeat = x_layer.repeat();

      // Loop over number of repeats for this layer.
      for (int j=0; j<repeat; j++)    {

        // move structural layer to before the slabs - djeans

        // deal with preshower or lack thereof ------------------
        double radiator_dim_Z(0);
        //double rad_pos_Z(0);
        double this_struct_CFthick_beforeAbs(0);
        double this_struct_CFthick_afterAbs(0);
        double _CF_absWrap = Ecal_fiber_thickness_structure;
        double _CF_alvWall = Ecal_fiber_thickness_alveolus;
        if ( isFrontFace ) { // the first part of the module depends on whether we have preshwer or not
          if ( Ecal_Barrel_PreshowerLayer==1 ) { // don't include W+CF wrapping; only one side of alveolus
            this_struct_CFthick_beforeAbs = 0;
            this_struct_CFthick_afterAbs = Ecal_front_face_thickness + _CF_alvWall;
            radiator_dim_Z = 0;
          } else { // include W+CF wrapping; only one side of alveolus
            this_struct_CFthick_beforeAbs = Ecal_front_face_thickness + _CF_absWrap;
            this_struct_CFthick_afterAbs = _CF_absWrap + _CF_alvWall;
            radiator_dim_Z = Ecal_radiator_thickness1; // this line implicitly assumes Ecal_nlayers1>0...
            //rad_pos_Z = radiator_dim_Z/2. + this_struct_CFthick_beforeAbs; // distance from top surface of structure to centre of radiator
            absorber_index++;
          }

          thickness_sum       += radiator_dim_Z;
          nRadiationLengths   += radiator_dim_Z/stave_material.radLength();
          nInteractionLengths += radiator_dim_Z/stave_material.intLength();

          thickness_sum       += ( this_struct_CFthick_beforeAbs + this_struct_CFthick_afterAbs );
          nRadiationLengths   += ( this_struct_CFthick_beforeAbs + this_struct_CFthick_afterAbs ) / CF.radLength();
          nInteractionLengths += ( this_struct_CFthick_beforeAbs + this_struct_CFthick_afterAbs ) / CF.intLength();

          isFrontFace=false;
        } else { // internal layer: include W+CF wrapping; both sides of alveolus
          this_struct_CFthick_beforeAbs = _CF_alvWall+_CF_absWrap;
          this_struct_CFthick_afterAbs = _CF_alvWall+_CF_absWrap;

          if ( absorber_index<Ecal_nlayers1 ) radiator_dim_Z=Ecal_radiator_thickness1;
          else if ( absorber_index<Ecal_nlayers1+Ecal_nlayers2 ) radiator_dim_Z=Ecal_radiator_thickness2;
          else if ( absorber_index<Ecal_nlayers1+Ecal_nlayers2+Ecal_nlayers3 ) radiator_dim_Z=Ecal_radiator_thickness3;
          else {
            assert(0 && "ERROR cannot determine absorber thickness for layer ");
          }

          absorber_index++;

          assert( radiator_dim_Z>0 && "no radiator!" );
        }
        //-------------------------------

        radiator_dim_y = radiator_dim_Z; // ahh different axis conventions....

        l_pos_z += this_struct_CFthick_beforeAbs + radiator_dim_y/2;

        // #########################
        // BuildECRingStructureLayer
        // #########################

        Tube CenterECTubForSi(0.,
                              Lcal_outer_radius + Ecal_Lcal_ring_gap + 0.001, // + tolerance,
                              module_thickness,
                              0.,
                              2 * M_PI);


        string l_name = _toString(l_num,"layer%d");
        RotationZYX rot(0,0,0); // a null rotation

        if ( radiator_dim_y>0 ) {

          string bs_name="bs";
          Box        EndcapStructureLayer_box(ECRingSiplateSize/ 2. - tolerance, ECRingSiplateSize/ 2. - tolerance, radiator_dim_y/2.);
          SubtractionSolid  EndcapStructureLayerSolid( EndcapStructureLayer_box, CenterECTubForSi, hole_position);
          Volume     EndcapStructureLayer_vol(det_name+"_"+l_name+"_"+bs_name,EndcapStructureLayerSolid,theDetector.material(x_staves.materialStr()));
          EndcapStructureLayer_vol.setVisAttributes(theDetector.visAttributes( x_staves.visStr()));

          // this is where we place the tungsten into the envelope
          double bsl_pos_z = l_pos_z;
          Position   bsl_pos(0,0,bsl_pos_z);
          Transform3D bsl_tran3D(rot,bsl_pos);
          EnvLogECRing.placeVolume(EndcapStructureLayer_vol,bsl_tran3D);
        }

        // Increment to next layer Z position.
        l_pos_z += radiator_dim_y/2 + this_struct_CFthick_afterAbs;

        // now move to the alveolus and what's inside
        double l_thickness = layering.layer(l_num-1)->thickness();  // Layer's thickness. this is the internal thickness of the alveolus

        l_pos_z  += l_thickness/2.; // add in half the alveolus thickness

        int EC_Number_of_towers = 0;

        // We use the same method able to create the Barrel
        // Slabs, so we have to rotate it later when placing
        // into the EC modules.
        //
        // We use the same method able to create the Barrel
        // radiator plates between slabs, but with the good
        // dimensions to avoid to rotate it later when placing
        // into the EC modules.

        // While the towers have the same shape use the same
        // logical volumes and parameters.

        string tower_name = _toString(EC_Number_of_towers,"tower%d");

        // build the ECRing layer, a box with hole in the middle
        Box  ECRingSiBox( ECRingSiplateSize/ 2. - tolerance, ECRingSiplateSize/ 2. - tolerance, l_thickness/2.0-tolerance);

        SubtractionSolid ECRingSiSolid( ECRingSiBox, CenterECTubForSi, hole_position);

        Volume     l_vol(det_name+"_"+l_name+"_"+tower_name,ECRingSiSolid,air);
        DetElement layer(module_det, l_name+tower_name, det_id);

        l_vol.setVisAttributes(theDetector.visAttributes( x_layer.visStr() ));

        // Loop over the sublayers or slices for this layer.
        int s_num = 1;
        double s_pos_z = -(l_thickness / 2);

        //--------------------------------------------------------------------------------
        // BuildECRing keep the Alveolus structure, the same as the Barrel and Endcap
        //--------------------------------------------------------------------------------

        for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
          xml_comp_t x_slice = si;
          string     s_name  =  _toString(s_num,"slice%d");
          double     s_thick = x_slice.thickness();
          Material slice_material  = theDetector.material(x_slice.materialStr());

          double slab_dim_x = ECRingSiplateSize/2.-tolerance;
          double slab_dim_y = s_thick/2.-tolerance;
          double slab_dim_z = ECRingSiplateSize/2.-tolerance;

          Box        s_box(slab_dim_x,slab_dim_z,slab_dim_y);

          SubtractionSolid ECRingSiSliceSolid( s_box, CenterECTubForSi, hole_position);

          Volume     s_vol(det_name+"_"+l_name+"_"+s_name,ECRingSiSliceSolid,slice_material);
          DetElement slice(layer,s_name,det_id);

          s_vol.setVisAttributes(theDetector.visAttributes(x_slice.visStr()));

#ifdef VERBOSE
          std::cout<<"x_slice.materialStr(): "<< x_slice.materialStr() <<std::endl;
#endif
          if (x_slice.materialStr().compare(x_staves.materialStr()) == 0){
            radiator_dim_y = s_thick;

            absorber_index++;

#if DD4HEP_VERSION_GE( 0, 15 )
            caloLayer.outer_nRadiationLengths   = nRadiationLengths;
            caloLayer.outer_nInteractionLengths = nInteractionLengths;
            caloLayer.outer_thickness           = thickness_sum;
            caloLayer.distance                  = distance_ip_to_front_face + dist_from_front_face_tally;
#endif
            if (Ecal_Barrel_PreshowerLayer==0 || !isFirstSens ){ // add this layer to the caloData (unless its the first absorber layer of a no-preshower calo)
              if ( module_num==0 ) {
                caloData->layers.push_back( caloLayer ) ;
#ifdef VERBOSE
#if DD4HEP_VERSION_GE( 0, 15 )
                std::cout<<" caloLayer.distance: "<< caloLayer.distance <<std::endl;
                std::cout<<" caloLayer.inner_nRadiationLengths: "<< caloLayer.inner_nRadiationLengths <<std::endl;
                std::cout<<" caloLayer.inner_nInteractionLengths: "<< caloLayer.inner_nInteractionLengths <<std::endl;
                std::cout<<" caloLayer.inner_thickness: "<< caloLayer.inner_thickness <<std::endl;
                std::cout<<" caloLayer.sensitive_thickness: "<< caloLayer.sensitive_thickness <<std::endl;
                std::cout<<" caloLayer.cellSize0, 1: "             << caloLayer.cellSize0 << " " << caloLayer.cellSize1 << std::endl;
                std::cout<<" caloLayer.outer_nRadiationLengths: "<< caloLayer.outer_nRadiationLengths <<std::endl;
                std::cout<<" caloLayer.outer_nInteractionLengths: "<< caloLayer.outer_nInteractionLengths <<std::endl;
                std::cout<<" caloLayer.outer_thickness: "<< caloLayer.outer_thickness <<std::endl;
                std::cout<<" EcalECRing[1]==>caloLayer.inner_thickness + caloLayer.outer_thickness: "
                         << caloLayer.inner_thickness + caloLayer.outer_thickness <<std::endl;
#endif
#endif
                dist_from_front_face_tally += caloLayer.inner_thickness+caloLayer.outer_thickness;
              }
            }

            nRadiationLengths   = 0. ;
            nInteractionLengths = 0. ;
            thickness_sum       = 0. ;
            isFirstSens         = false;
          } // if radiator

          nRadiationLengths   += s_thick/(2.*slice_material.radLength());
          nInteractionLengths += s_thick/(2.*slice_material.intLength());
          thickness_sum       += s_thick/2;

          if ( x_slice.isSensitive() ) {
            s_vol.setSensitiveDetector(sens);

#if DD4HEP_VERSION_GE( 0, 15 )
            //Store "inner" quantities
            caloLayer.inner_nRadiationLengths   = nRadiationLengths;
            caloLayer.inner_nInteractionLengths = nInteractionLengths;
            caloLayer.inner_thickness           = thickness_sum;
            //Store sensitive slice thickness
            caloLayer.sensitive_thickness       = s_thick;
#ifdef VERBOSE
            std::cout<<" l_num: "<<l_num <<std::endl;
            std::cout<<" s_num: "<<s_num <<std::endl;
            std::cout<<" module_thickness: "<< module_thickness <<std::endl;
            std::cout<<" l_pos_z: "<< l_pos_z <<std::endl;
            std::cout<<" l_thickness: "<< l_thickness <<std::endl;
            std::cout<<" s_pos_z: "<< s_pos_z <<std::endl;
            std::cout<<" s_thick: "<< s_thick <<std::endl;
            std::cout<<" radiator_dim_y: "<< radiator_dim_y <<std::endl;
#endif
            caloLayer.absorberThickness = radiator_dim_y ;
#endif
            nRadiationLengths   = 0. ;
            nInteractionLengths = 0. ;
            thickness_sum       = 0. ;
          } // if sensitive

          nRadiationLengths   += s_thick/(2.*slice_material.radLength());
          nInteractionLengths += s_thick/(2.*slice_material.intLength());
          thickness_sum       += s_thick/2;

          slice.setAttributes(theDetector,s_vol,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());

          // Slice placement.
          PlacedVolume slice_phv = l_vol.placeVolume(s_vol,Position(0,0,s_pos_z+s_thick/2));

          if ( x_slice.isSensitive() ) {
            slice_phv.addPhysVolID("layer", myLayerNum++);
          }

          slice.setPlacement(slice_phv);
          // Increment Z position of slice.
          s_pos_z += s_thick;

          // Increment slice number.
          ++s_num;
        }

        // add material between the "slabs"

        // we need to use correct thickness here! especially at the interface between stacks.
        if      ( absorber_index < Ecal_nlayers1 ) radiator_dim_y=Ecal_radiator_thickness1;
        else if ( absorber_index < Ecal_nlayers1+Ecal_nlayers2 ) radiator_dim_y=Ecal_radiator_thickness2;
        else if ( absorber_index < Ecal_nlayers1+Ecal_nlayers2+Ecal_nlayers3 ) radiator_dim_y=Ecal_radiator_thickness3;
        else {
          radiator_dim_y=0;
        }
	
        // deal with final support plate
        float cf_thick = Ecal_fiber_thickness_alveolus;
        if ( absorber_index<Ecal_nlayers1+Ecal_nlayers2+Ecal_nlayers3 ) {
          cf_thick += Ecal_fiber_thickness_structure;
        } else {
          cf_thick += Ecal_support_thickness;
        }

        if ( module_num==0 ) {
#if DD4HEP_VERSION_GE( 0, 15 )
          caloLayer.distance = distance_ip_to_front_face + dist_from_front_face_tally;
          caloLayer.outer_nRadiationLengths   = nRadiationLengths + cf_thick/CF.radLength();
          caloLayer.outer_nInteractionLengths = nInteractionLengths + cf_thick/CF.intLength();
          caloLayer.outer_thickness           = thickness_sum + cf_thick;
#endif
          caloData->layers.push_back( caloLayer ) ;

#ifdef VERBOSE
#if DD4HEP_VERSION_GE( 0, 15 )
          std::cout<<" caloLayer.distance: "<< caloLayer.distance <<std::endl;
          std::cout<<" caloLayer.inner_nRadiationLengths: "<< caloLayer.inner_nRadiationLengths <<std::endl;
          std::cout<<" caloLayer.inner_nInteractionLengths: "<< caloLayer.inner_nInteractionLengths <<std::endl;
          std::cout<<" caloLayer.inner_thickness: "<< caloLayer.inner_thickness <<std::endl;
          std::cout<<" caloLayer.sensitive_thickness: "<< caloLayer.sensitive_thickness <<std::endl;
          std::cout<<" caloLayer.cellSize0, 1: "             << caloLayer.cellSize0 << " " << caloLayer.cellSize1 << std::endl;
          std::cout<<" caloLayer.outer_nRadiationLengths: "<< caloLayer.outer_nRadiationLengths <<std::endl;
          std::cout<<" caloLayer.outer_nInteractionLengths: "<< caloLayer.outer_nInteractionLengths <<std::endl;
          std::cout<<" caloLayer.outer_thickness: "<< caloLayer.outer_thickness <<std::endl;

          std::cout<<" EcalECRing[2]==>caloLayer.inner_thickness + caloLayer.outer_thickness: "
                   << caloLayer.inner_thickness + caloLayer.outer_thickness <<std::endl;
#endif
#endif
          dist_from_front_face_tally += caloLayer.inner_thickness+caloLayer.outer_thickness;
        }

        nRadiationLengths   = radiator_dim_y/stave_material.radLength() + cf_thick/CF.radLength();
        nInteractionLengths = radiator_dim_y/stave_material.intLength() + cf_thick/CF.intLength();
        thickness_sum       = radiator_dim_y + cf_thick;

        // this is where we place the slab (l_vol) into the envelope
        int i_stave = 1;
        int i_tower = 1;
        Position l_pos(0,0,l_pos_z);
        Transform3D tran3D(rot,l_pos);
        PlacedVolume layer_phv = EnvLogECRing.placeVolume(l_vol,tran3D);
        layer_phv.addPhysVolID("tower", i_tower);
        layer_phv.addPhysVolID("stave", i_stave);
        layer.setPlacement(layer_phv);

        l_pos_z +=   l_thickness/2.; // add in the other half of the alveolus

        ++l_num;

      }
    }


    // Set stave visualization.
    if (x_staves)   {
      EnvLogECRing.setVisAttributes(theDetector.visAttributes(x_staves.visStr()));
    }


    //====================================================================
    // Place Ecal Endcap module into the assembly envelope volume
    //====================================================================

    double EC_module_z_offset = (EcalEndcapRing_min_z + EcalEndcapRing_max_z)/2.;

    int module_id = ( module_num == 0 ) ? 0:6;
    double this_module_z_offset = ( module_id == 0 ) ? - EC_module_z_offset : EC_module_z_offset;
    double this_module_rotY = ( module_id == 0 ) ? M_PI:0;

    Position xyzVec(0,0,this_module_z_offset);
    RotationZYX rot(0,this_module_rotY,0);
    Rotation3D rot3D(rot);
    Transform3D tran3D(rot3D,xyzVec);

    PlacedVolume pv = envelope.placeVolume(EnvLogECRing,tran3D);
    pv.addPhysVolID("module",module_id); // z: -/+ 0/6

    DetElement sd = module_det;
    sd.setPlacement(pv);

  }

  sdet.addExtension< dd4hep::rec::LayeredCalorimeterData >( caloData ) ;

  return sdet;

}



DECLARE_DETELEMENT(SEcal05_ECRing, create_detector)


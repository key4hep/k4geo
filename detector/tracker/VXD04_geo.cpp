//====================================================================
//  lcgeo - LC detector models in DD4hep 
//--------------------------------------------------------------------
//  F.Gaede, DESY
//  $Id$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"

#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"
#include "XMLHandlerDB.h"

//#include "DDRec/DDGear.h"
//#define MOKKA_GEAR

#include <cmath>

using namespace std;
using namespace DD4hep;
using namespace dd4hep ;
using namespace DD4hep::Geometry;
using namespace DDRec;




/** Construction of VTX detector, ported from Mokka driver VXD04.cc
 *
 *  Mokka History:
 * - first implementation -- Damien Grandjean, April 2003
 * - fixed geometry overlap -- Adrian Vogel, 2005-12-12
 * - added optional GEAR output -- R. Lippe, DESY, 2006-09-04
 * -modification for double layer geometry -- Damien Grandjean, February 2008
 * -increased realism in the description of the ladders, the Be support and the cabling, added cooling tubes Y. Voutsinas, September 2011
 *
 *  @author: F.Gaede, DESY, Nov 2013
 *
 */
static Ref_t create_element(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
  
  xml_det_t    x_det = e;
  string       name  = x_det.nameStr();
  
  DetElement   vxd(  name, x_det.id()  ) ;
  
  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = XML::createPlacedEnvelope( lcdd,  e , vxd ) ;
  
  XML::setDetectorTypeFlag( e, vxd ) ;

  if( lcdd.buildType() == BUILD_ENVELOPE ) return vxd ;

  //-----------------------------------------------------------------------------------

  sens.setType("tracker");

  // --- create assembly and DetElement for support and service volumes 

  Assembly supp_assembly( name + "_support_assembly"  ) ;

  PlacedVolume pv = envelope.placeVolume( supp_assembly ) ;

  DetElement suppDE( vxd , name+"_support" , x_det.id() )  ;
  suppDE.setPlacement( pv ) ;
  //--------------------------------


  //######################################################################################################################################################################
  //  code ported from VXD04::construct() :
  //##################################
  //------------------------------------------
  //  See comments starting with '//**' for
  //     hints on porting issues
  //------------------------------------------
  
  // double sPhi = 0 * deg;
  // double dPhi = 360 * deg;
  //** DD4hep/TGeo seems to need rad (as opposed to the manual)
  double sPhi = 0 ;
  double dPhi = 2*M_PI;
  
  
  //****************************************
  // Layers
  //****************************************
  //
  // Common layer thickness parameters
  // db->exec("select * from layers_common_parameters;");
  XMLHandlerDB db(  x_det.child( _Unicode( layers_common_parameters ) ) );
			  
  double foam_spacer_thickness = db->fetchDouble("foam_spacer_thickness");
  double flex_cable_thickness = db->fetchDouble("flex_cable_thickness");
  double metal_traces_thickness = db->fetchDouble("metal_traces_thickness");
  double electronics_structure_thickness = db->fetchDouble("electronics_structure_thickness");
  double active_silicon_thickness =   db->fetchDouble("active_silicon_thickness");
  //unused:  double support_structure_radial_thickness =    db->fetchDouble("support_structure_radial_thickness");
  double end_electronics_half_z=    db->fetchDouble("end_electronics_half_z");
  //unused: double strip_final_beampipe_radious =    db->fetchDouble("strip_final_beampipe_radious");
  int    side_band_electronics_option=    db->fetchInt("side_band_electronics_option");
  std::string  flex_cable_material=    db->fetchString("flex_cable_material");
  std::string  metal_traces_material=    db->fetchString("metal_traces_material");
  std::string  foam_spacer_material=    db->fetchString("foam_spacer_material");
  std::string  cool_pipe_material=    db->fetchString("cool_pipe_material");
  int  end_ladd_electronics_option=    db->fetchInt("end_ladd_electronics_option"); 
  double  side_band_electronics_width=    db->fetchDouble("side_band_electronics_width");
  double  side_band_electronics_thickness=    db->fetchDouble("side_band_electronics_thickness");
  int  active_side_band_electronics_option=  db->fetchInt("active_side_band_electronics_option");
  double  layer_gap= db->fetchDouble("layer_gap");
  double  cool_pipe_outer_radius = db->fetchDouble("cool_pipe_outer_radius");
  double  cool_pipe_inner_radius =    db->fetchDouble("cool_pipe_inner_radius");
  double  external_kapton_thickness = db->fetchDouble("external_kapton_thickness");
  double  external_metal_thickness =  db->fetchDouble("external_metal_thickness");

  //Cryostat parameters
  // db->exec("SELECT * FROM cryostat;");
  db = XMLHandlerDB(  x_det.child( _Unicode( cryostat ) ) ) ;

  double rAlu   = db->fetchDouble("alu_skin_inner_radious") ;
  double drAlu  = db->fetchDouble("alu_skin_tickness") ;
  const double rSty   = db->fetchDouble("foam_inner_radious") ;
  const double drSty  = db->fetchDouble("foam_tickness") ;
  const double dzSty  = db->fetchDouble("foam_half_z") ;
  const double cryostat_apperture  = db->fetchDouble("cryostat_apperture") ;
  const double cryostat_apperture_radius  = db->fetchDouble("cryostat_apperture_radius") ;
  double rInner = db->fetchDouble("endplate_inner_radious") ;
  // double rAlu   = db->fetchDouble("alu_skin_inner_radious") * mm;
  // double drAlu  = db->fetchDouble("alu_skin_tickness") * mm;
  // const double rSty   = db->fetchDouble("foam_inner_radious") * mm;
  // const double drSty  = db->fetchDouble("foam_tickness") * mm;
  // const double dzSty  = db->fetchDouble("foam_half_z") * mm;
  // const double cryostat_apperture  = db->fetchDouble("cryostat_apperture") * mm;
  // const double cryostat_apperture_radius  = db->fetchDouble("cryostat_apperture_radius") * mm;
  // double rInner = db->fetchDouble("endplate_inner_radious") * mm;

  bool useCryo  = bool(db->fetchInt("cryostat_option"));

  // support shell parameters
  // db->exec("select * from support_shell;");
  db = XMLHandlerDB(  x_det.child( _Unicode( support_shell ) ) ) ;
  
  double shell_inner_radious = db->fetchDouble("inner_radious");
  double shell_half_z = db->fetchDouble("half_z");
  double shell_thickess = db->fetchDouble("thickess");
  double support_endplate_inner_radious = db->fetchDouble("endplate_inner_radious");
  double support_endplate_inner_radious_L1 = db->fetchDouble("endplate_inner_radius_L1");
  double support_endplate_outer_radious_L1 = db->fetchDouble("endplate_outer_radius_L1");
  //unused:  double offset_ladder_block = db->fetchDouble("offset_ladder_block");
  double beryllium_ladder_block_length = db->fetchDouble("beryllium_ladder_block_length");
  double beryllium_ladder_block_thickness = db->fetchDouble("beryllium_ladder_block_thickness");
  double beryllium_ladder_block_length2=0.;
  double shell_endplate_thickness = db->fetchDouble("shell_endplate_thickness");
  double forward_shell_half_z = db->fetchDouble("forward_shell_half_z");
  

  // ### fixme: SD ##############
  //   // The VXD Sensitive detector
  //   // Threshold is 20% of a MIP. For Si we have
  //   // 340 KeV/mm as MIP.
  //   theVXDSD =
  //     new TRKSiSD00("VXD",
  // 		  active_silicon_thickness * mm
  // 		  * 340 * keV
  // 		  * 0.2);
  //   RegisterSensitiveDetector(theVXDSD);
  

  //**fg: the encoder is no longer needed - replaced by physVolID() calls
  //   // setup the encoder 
  //   UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  //   encoder.reset() ;  // reset to 0
  //   encoder[ILDCellID0::subdet] = ILDDetID::VXD ;
  //   encoder[ILDCellID0::side] = 0 ;
  //   encoder[ILDCellID0::layer]  = 0 ;
  //   encoder[ILDCellID0::module] = 0 ;
  //   encoder[ILDCellID0::sensor] = 0 ;
  //   int cellID0 = encoder.lowWord() ;

  Material activeMaterial =  lcdd.material("G4_Si"); //silicon_2.33gccm"); 
  

  DDRec::ZPlanarData*  zPlanarData = new ZPlanarData ;

// #ifdef MOKKA_GEAR
//   // some variables for storing information for MOKKA_GEAR
//   // during the loop
//   std::vector<helpLayer> gearHelpLadders ;
//   std::vector<helpLayer> gearHelpSensitives ;
//   std::vector<int> gearHelpNumberLadders ;
//   std::vector<double> gearHelpPhi0 ;
//   double gearHelpGap = 0. ;
//   int gearHelpCount = 0 ;
//   int gearHelpType = 0 ;
// #endif
  

  // db->exec("select * from layer;");
  //   do{
  //**fg: get parameters for first layer - needed below
  double ladder_0_length = 0 ;

  for(xml_coll_t c( x_det ,_U(layer)); c; ++c)  {
    
    xml_comp_t  x_layer( c );
    db = XMLHandlerDB( x_layer )  ;

    int LayerId = db->fetchInt("id");
    double layer_radius = db->fetchDouble("layer_radius");
    double ladder_length  = db->fetchDouble("ladder_length");
    if( LayerId == 0 ) {
      ladder_0_length = ladder_length ; 
    }
    double ladder_width = db->fetchDouble("ladder_width");
    double support_width = db->fetchDouble("support_width");
    double ladder_gap = db->fetchDouble("ladder_gap");
    //unused:    double strip_line_final_z = db->fetchDouble("strip_line_final_z");
    double initial_kapton_striplines_thickness = db->fetchDouble("initial_kapton_striplines_thickness");
    double final_kapton_striplines_thickness = db->fetchDouble("final_kapton_striplines_thickness");
    double initial_metal_striplines_thickness = db->fetchDouble("initial_metal_striplines_thickness");
    double final_metal_striplines_thickness = db->fetchDouble("final_metal_striplines_thickness");

#ifdef LCIO_MODE
    ladder_gapVec.push_back(ladder_gap);
    StripLineFinalZ_Vec.push_back(strip_line_final_z);
#endif
    double nb_ladder = db->fetchDouble("nb_ladder");
    
    double phirot = 0.;
    double phirot2 = 0.;
    
    std::cout << " ############## layer : " << LayerId << " number of ladders : " << nb_ladder << std::endl ; 


    Assembly layer_assembly( _toString( LayerId , "layer_assembly_%d"  ) ) ;
    PlacedVolume pv = envelope.placeVolume( layer_assembly ) ;


    //replacing support ladder with flex cable (kapton+metal) & adding a foam spacer
    // ****************************************************************************************
    // **********************   flex  cable *****************************************
    // ****************************************************************************************
    
    Material flexCableMaterial =  lcdd.material( flex_cable_material ); 
    
    Box FlexCableSolid( ladder_width+(side_band_electronics_option*side_band_electronics_width/2.),
			ladder_length+(end_ladd_electronics_option*(2*end_electronics_half_z)) + beryllium_ladder_block_length*2.,
			flex_cable_thickness/2.);

    //** ----- Original Mokke/Geant4 code: -----
    // Box *FlexCableSolid
    //  = new Box("FlexCable",
    // 		ladder_width+(side_band_electronics_option*side_band_electronics_width/2.),
    // 		ladder_length+(end_ladd_electronics_option*(2*end_electronics_half_z)) + beryllium_ladder_block_length*2.,
    // 		flex_cable_thickness/2.);
    
    //    VisAttributes* flex_cableVisAtt = new VisAttributes(Colour(1.,0.,0.,1.0));   //red
    
    //**fg: we need distinct names for every instance of a Volume - so we append _(layer#) to the volume name
    Volume FlexCableLogical( _toString(LayerId,"FlexCable_%02d"), FlexCableSolid, flexCableMaterial ) ;
    
    // ----- Original Mokke/Geant4 code: -----
    // LogicalVolume *FlexCableLogical=
    //   new LogicalVolume(FlexCableSolid,
    //  			flexCableMaterial,
    //  			"FlexCable",
    //  			0,
    //  			0,
    //  			0);
    
    vxd.setVisAttributes(lcdd,  "RedVis" , FlexCableLogical);
    //** ----- Original Mokke/Geant4 code: -----
    //    FlexCableLogical->SetVisAttributes(flex_cableVisAtt);
    
    // ****************************************************************************************
    // **********************   metal traces  *****************************************
    // ****************************************************************************************
    
    Material metalTracesMaterial = lcdd.material( metal_traces_material); 
    
    Box MetalTracesSolid( ladder_width+(side_band_electronics_option*side_band_electronics_width/2.),
			  ladder_length+(end_ladd_electronics_option*(2*end_electronics_half_z)) + beryllium_ladder_block_length*2.,
			  metal_traces_thickness/2.);
    
    Volume MetalTracesLogical( _toString(LayerId,"MetalTraces_%02d") , MetalTracesSolid,metalTracesMaterial ) ;
    
    vxd.setVisAttributes(lcdd,  "GrayVis" , MetalTracesLogical) ;
    
    // ****************************************************************************************
    // **********************   foam spacer n support  *****************************************
    // ****************************************************************************************
    
    Material foamSpacerMaterial = lcdd.material( foam_spacer_material);
    
    Box FoamSpacerSolid( support_width+(side_band_electronics_option*side_band_electronics_width/2.),
			 ladder_length+(end_ladd_electronics_option*(2*end_electronics_half_z)) + beryllium_ladder_block_length*2. ,
			 foam_spacer_thickness/2.);
    
    Volume FoamSpacerLogical( _toString(LayerId,"FoamSpacer_%02d"), FoamSpacerSolid, foamSpacerMaterial) ;
      
    vxd.setVisAttributes(lcdd, "YellowVis", FoamSpacerLogical ) ;

    //here we place the physical volumes of both the flex cable (kapton & metal traces) and the foam spacer
    
    phirot = (2*pi)/nb_ladder;
    
    double ladder_clothest_approch = beryllium_ladder_block_thickness*2 +0.1;

    // calculate optimal offset, such that there is 0.1mm space between to the edge and the surface of two adjacent ladders.
    // in the case of ladders overlapped per superlayer
    /*
      double offset_phi=(1-cos(phirot))/sin(phirot)*layer_radius  
      -((ladder_width+(side_band_electronics_option*side_band_electronics_width/2.))
      +(ladder_clothest_approch+cos(phirot)*2*(foam_spacer_thickness+active_silicon_thickness+flex_cable_thickness+metal_traces_thickness))/sin(phirot));
    */
    // in the case of ladders overlapped per layer
    double offset_phi=(1-cos(phirot))/sin(phirot)*layer_radius  
      -((ladder_width+(side_band_electronics_option*side_band_electronics_width/2.))
	+(ladder_clothest_approch+cos(phirot)*2*(active_silicon_thickness+flex_cable_thickness+metal_traces_thickness))/sin(phirot));
    
    if (LayerId==0||LayerId==2||LayerId==4)  {  //------------------------------------------------------------------------
       
      for (double ladder_loop=0;ladder_loop<nb_ladder;ladder_loop++) {
	
	double phirot2 = ladder_loop*phirot;

	// RotationMatrix *rot = new RotationMatrix();
	// rot->rotateX(pi*0.5);
	// rot->rotateY(phirot2);
	RotationZYX rot( 0, phirot2 , (pi*0.5) ) ;
	
	pv = supp_assembly.placeVolume( FlexCableLogical,  
				   Transform3D( rot, Position(( layer_radius + metal_traces_thickness + (flex_cable_thickness/2.))*sin(phirot2)+offset_phi*cos(phirot2),
							      -(layer_radius + metal_traces_thickness + (flex_cable_thickness/2.))*cos(phirot2)+offset_phi*sin(phirot2),
							      0.))  );
	// Phys=
	//   new PVPlacement(rot,
	// 			ThreeVector((layer_radius + metal_traces_thickness + (flex_cable_thickness/2.))*sin(phirot2)+offset_phi*cos(phirot2),
	// 				      -(layer_radius + metal_traces_thickness + (flex_cable_thickness/2.))*cos(phirot2)+offset_phi*sin(phirot2),
	// 				      0.),
	// 			FlexCableLogical,
	// 			"FlexCable",
	// 			worldLog,
	// 			false,
	// 			0);
	       
	pv = supp_assembly.placeVolume( FoamSpacerLogical,  
				   Transform3D(  rot, Position((layer_radius + flex_cable_thickness + metal_traces_thickness + foam_spacer_thickness/2.)*sin(phirot2)+offset_phi*cos(phirot2),
							       -(layer_radius + flex_cable_thickness + metal_traces_thickness +  foam_spacer_thickness/2.)*cos(phirot2)+offset_phi*sin(phirot2),
							       0.))  );

	pv = supp_assembly.placeVolume( MetalTracesLogical,  Transform3D( rot,Position((layer_radius + (metal_traces_thickness/2))*sin(phirot2)+offset_phi*cos(phirot2),
										  -(layer_radius + (metal_traces_thickness/2.))*cos(phirot2)+offset_phi*sin(phirot2),
 										  0.))  );
      }
      
    } else if (LayerId==1||LayerId==3||LayerId==5) { //------------------------------------------------------------------------
      
      for (double ladder_loop=0;ladder_loop<nb_ladder;ladder_loop++) {
	
	double phirot2 = ladder_loop*phirot;
	
	RotationZYX rot( 0, phirot2 , (pi*0.5) ) ;
	
	pv = supp_assembly.placeVolume( FlexCableLogical,
				   Transform3D( rot, Position((layer_radius-(metal_traces_thickness + flex_cable_thickness/2.)+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
							      -(layer_radius-(metal_traces_thickness + flex_cable_thickness/2.)+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
							      0.))  ) ;

	pv = supp_assembly.placeVolume( FoamSpacerLogical, 
				   Transform3D( rot, Position((layer_radius + layer_gap - flex_cable_thickness -  metal_traces_thickness - foam_spacer_thickness/2.)*sin(phirot2)+offset_phi*cos(phirot2),
							      -(layer_radius + layer_gap - flex_cable_thickness - metal_traces_thickness - foam_spacer_thickness/2.)*cos(phirot2)+offset_phi*sin(phirot2),
							      0.))  );
	
	pv = supp_assembly.placeVolume( MetalTracesLogical,  Transform3D( rot,Position((layer_radius-(metal_traces_thickness/2)+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
										  -(layer_radius-(metal_traces_thickness/2.)+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
										  0.))  );
      }
    }

// #ifdef MOKKA_GEAR
    
//     //Definition of the VXDSupport composite material. It is going to be used only during the reconstruction stage, for tracking purposes. It consists by three layers: metal traces, flex cable and the foam spacer support with user defined materials and thicknesses. Here we define the element and calculate its effective radiation length, atomic number and atomic mass. For the simulation, the more realistic 3 different layers structure is being used.   

//     double MetalDensity = metalTracesMaterial->GetMaterial()->GetDensity()/(g/mm3);
//     double KaptonDensity = flexCableMaterial->GetMaterial()->GetDensity()/(g/mm3);
//     double FoamDensity = foamSpacerMaterial->GetMaterial()->GetDensity()/(g/mm3);

//     double VXDSupportThickness = metal_traces_thickness + flex_cable_thickness + foam_spacer_thickness;

//     //calculations of thickness fractions of each layer of the support
//     double metalTF = metal_traces_thickness / VXDSupportThickness;
//     double foamTF = foam_spacer_thickness / VXDSupportThickness;
//     double flexTF = flex_cable_thickness / VXDSupportThickness;

//     double elemVol = 1/(mm2);

//     double VXDSupportMass = foam_spacer_thickness*(elemVol)*FoamDensity + flex_cable_thickness*(elemVol)*KaptonDensity + metal_traces_thickness*(elemVol)*MetalDensity;

//     double VXDSupportDensity = VXDSupportMass/1/(mm3) ;

//     double foamFM = 100. * ((foam_spacer_thickness*(elemVol)*FoamDensity) / VXDSupportMass) ;
//     double kaptonFM = 100. * ((flex_cable_thickness*(elemVol)*KaptonDensity) / VXDSupportMass) ;
//     double metalFM = 100. * ((metal_traces_thickness*(elemVol)*MetalDensity) / VXDSupportMass) ;

//     //Calculation of an effective radiation length for the support based on the mass fraction of each material

//     double VXDSupportRadLen = 1. / ((metalTF/metalTracesMaterial->GetMaterial()->GetRadLen()) + (flexTF/flexCableMaterial->GetMaterial()->GetRadLen()) + (foamTF/foamSpacerMaterial->GetMaterial()->GetRadLen()));

//     //Calculation of the effective atomic number of the VXD support. The Z effectives are obtained from the formula: Zeff = Sum(Wi*Zi) where Wi are the mass fractions of the elements that consist the material 

//     Material *carbon = CGAGeometryManager::GetMaterial("carbon");
//     Material *silicon = CGAGeometryManager::GetMaterial("silicon");
//     Material *hydrogen = CGAGeometryManager::GetMaterial("H2");
//     Material *nitro = CGAGeometryManager::GetMaterial("N2");
//     Material *oxygen = CGAGeometryManager::GetMaterial("oxygen");

//     double C_Z = carbon->GetZ();
//     double Si_Z = silicon->GetZ();
//     double C_A = carbon->GetA()/g;
//     double Si_A = silicon->GetA()/g;
//     double H_Z = hydrogen->GetZ();
//     double H_A = hydrogen->GetA()/g;
//     double N_Z = nitro->GetZ();
//     double N_A = nitro->GetA()/g;
//     double O_Z = oxygen->GetZ();
//     double O_A = oxygen->GetA()/g;


//     double foamZeff = C_Z*(C_A/(C_A+Si_A)) + Si_Z*(Si_A/(C_A+Si_A));


//     double metalZ = metalTracesMaterial->GetZ();
//     double metalA = metalTracesMaterial->GetA()/g;
      
//     //Calculation of kapton effective Z - weight fractions for each element taken from NIST dB

//     double flexZeff = H_Z*0.026362 + C_Z*0.691133 + N_Z*0.073270 + O_Z*0.209235;

//     double VXDSupportZeff = (metalFM/100.)*metalZ + (kaptonFM/100.)*flexZeff + (foamFM/100.)*foamZeff;


//     //Calculation of the effective atomic mass of the VXD support. The Z effectives are obtained from the formula: Aeff = Zeff / (Z/A)eff where (Z/A)eff = Sum Wi*Zi/Ai

//     double metalZA = metalZ/metalA;
//     double foamZAeff = (C_A/(C_A+Si_A))*(C_Z/C_A) + (Si_A/(C_A+Si_A))*(Si_Z/Si_A);
//     double flexZAeff = (H_Z/H_A)*0.026362 + (C_Z/C_A)*0.691133 + (N_Z/N_A)*0.073270 + (O_Z/O_A)*0.209235;

//     double VXDSupportZAeff = (metalFM/100.)*metalZA + (kaptonFM/100.)*flexZAeff + (foamFM/100.)*foamZAeff;

//     double VXDSupportAeff = VXDSupportZeff / VXDSupportZAeff;

//     //Calculation of the effective nuclear interaction length of the VXD support

//     double VXDSupportIntLength = 1. / ((metalTF/metalTracesMaterial->GetNuclearInterLength()) + (flexTF/flexCableMaterial->GetNuclearInterLength()) + (foamTF/foamSpacerMaterial->GetNuclearInterLength()));

//     //Here we call the SimpleMaterial class of gear. The density should be converted to kg/m3
//     VXDSupportDensity = 1000000*VXDSupportDensity;

//     VXDSupportMaterial = new gear::SimpleMaterialImpl("VXDSupportMaterial", VXDSupportAeff, VXDSupportZeff, VXDSupportDensity, VXDSupportRadLen, VXDSupportIntLength );

//     //_________________________________________________________________________________________________________
//     //

//     helpLayer thisLadder ;
//     if (LayerId==2||LayerId==4||LayerId==6) 
//       { 
// 	thisLadder.distance  = layer_radius + layer_gap * 0.5 ;
//       }
//     if (LayerId==1||LayerId==3||LayerId==5) 
//       { 
// 	thisLadder.distance  = layer_radius  ;
//       }      
//     //      thisLadder.distance  = layer_radius ;
//     thisLadder.offset    = offset_phi ;
//     thisLadder.thickness = VXDSupportThickness ;
//     thisLadder.length    = ladder_length ;
//     thisLadder.width     = (ladder_width*2.)+(side_band_electronics_option*side_band_electronics_width) ;
//     thisLadder.radLength = VXDSupportMaterial->GetMaterial()->getRadLength()/mm ;

 
//     // find out type
//     if( side_band_electronics_option == 0 &&  end_ladd_electronics_option == 1) gearHelpType = gear::ZPlanarParametersImpl::CCD  ;
//     if( side_band_electronics_option == 1 &&  end_ladd_electronics_option == 0 ) gearHelpType = gear::ZPlanarParametersImpl::CMOS ;
//     if( side_band_electronics_option == 1 &&  end_ladd_electronics_option == 1) gearHelpType = gear::ZPlanarParametersImpl::HYBRID ;

// #endif

    DDRec::ZPlanarData::LayerLayout thisLayer ;
    
    if (LayerId==1||LayerId==3||LayerId==5) { 
      
      thisLayer.distanceSupport  = layer_radius + layer_gap * 0.5 ;
      
    }else if (LayerId==0||LayerId==2||LayerId==4) {
      
      thisLayer.distanceSupport  = layer_radius  ;
    }      
    
    thisLayer.offsetSupport    = offset_phi ;
    thisLayer.thicknessSupport = metal_traces_thickness + flex_cable_thickness + foam_spacer_thickness ;
    thisLayer.zHalfSupport    = ladder_length ;
    thisLayer.widthSupport     = (ladder_width*2.)+(side_band_electronics_option*side_band_electronics_width) ;
    //     thisLayer.radLength = VXDSupportMaterial->GetMaterial()->getRadLength()/mm ;
    
    


    // ****************************************************************************************
    // **********************   Berylium annulus block *****************************************
    // ****************************************************************************************
    
    //only one block per superlayer

    if (LayerId==1) {
      
      Box BerylliumAnnulusBlockSolid( ladder_width, beryllium_ladder_block_length, beryllium_ladder_block_thickness);
      
      Volume BerylliumAnnulusBlockLogical( _toString(LayerId,"BerylliumAnnulusBlock_%02d"), BerylliumAnnulusBlockSolid, lcdd.material("G4_Be")) ; //"beryllium") ) ;
      
      vxd.setVisAttributes(lcdd,  "CyanVis" , BerylliumAnnulusBlockLogical) ;


      //====== create the meassurement surface for Be annulus block ===================
      Vector3D u( 1. , 0. , 0. ) ;
      Vector3D v( 0. , 1. , 0. ) ;
      Vector3D n( 0. , 0. , 1. ) ;
      
      VolPlane surfAnnBlock( BerylliumAnnulusBlockLogical , SurfaceType(SurfaceType::Helper) , beryllium_ladder_block_thickness/2. ,  beryllium_ladder_block_thickness/2. , u,v,n ) ; //,o ) ;
      //============================================================

      
      for (double AnnulusBlock_loop=0;AnnulusBlock_loop<nb_ladder;AnnulusBlock_loop++) {

	std::string annBlockNameP =  _toString( LayerId , "BerylliumAnnulusBlock_%02d_posZ") + _toString( (int)AnnulusBlock_loop, "_%02d" ) ;
	std::string annBlockNameN =  _toString( LayerId , "BerylliumAnnulusBlock_%02d_negZ") + _toString( (int)AnnulusBlock_loop, "_%02d" ) ;
	
	phirot2 = phirot*AnnulusBlock_loop;

	RotationZYX rot( 0, phirot2 ,  pi*0.5 ) ;
	
	double ZAnnulusBlock = ladder_length + end_electronics_half_z + (beryllium_ladder_block_length*2.);
	    
	pv = supp_assembly.placeVolume( BerylliumAnnulusBlockLogical,  Transform3D( rot, Position((layer_radius+beryllium_ladder_block_thickness+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
											     -(layer_radius+beryllium_ladder_block_thickness+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
											     ZAnnulusBlock))  ) ;
	DetElement  annBlockPosZ( vxd , annBlockNameP  , x_det.id() );
	annBlockPosZ.setPlacement( pv ) ;
	volSurfaceList( annBlockPosZ )->push_back( surfAnnBlock ) ;
	
	pv = supp_assembly.placeVolume( BerylliumAnnulusBlockLogical,  Transform3D( rot, Position((layer_radius+beryllium_ladder_block_thickness+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
											     -(layer_radius+beryllium_ladder_block_thickness+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
											     -ZAnnulusBlock))  );
	DetElement  annBlockNegZ( vxd , annBlockNameN  , x_det.id() );
	annBlockNegZ.setPlacement( pv ) ;
	volSurfaceList( annBlockNegZ )->push_back( surfAnnBlock ) ;
      }	

    } else if (LayerId==3||LayerId==5)  { 
      
      beryllium_ladder_block_length2 = beryllium_ladder_block_length + (shell_half_z - (end_electronics_half_z *3.* end_ladd_electronics_option)-ladder_length);
      
      for (double AnnulusBlock_loop=0;AnnulusBlock_loop<nb_ladder;AnnulusBlock_loop++) {
	
	Box BerylliumAnnulusBlockSolid( ladder_width, beryllium_ladder_block_length2/2., beryllium_ladder_block_thickness);
	
	//**fg: need to create unique string name per Volume object
	std::string volName = _toString(LayerId,"BerylliumAnnulusBlock_%02d") ;
	volName +=  _toString( int(AnnulusBlock_loop), "_%02d");

	Volume BerylliumAnnulusBlockLogical( volName , BerylliumAnnulusBlockSolid, lcdd.material("G4_Be")) ; //"beryllium") ) ;

	//====== create the meassurement surface for Be annulus block ===================
	Vector3D u( 1. , 0. , 0. ) ;
	Vector3D v( 0. , 1. , 0. ) ;
	Vector3D n( 0. , 0. , 1. ) ;
	
	VolPlane surfAnnBlock( BerylliumAnnulusBlockLogical , SurfaceType(SurfaceType::Helper) , beryllium_ladder_block_thickness/2. ,  beryllium_ladder_block_thickness/2. , u,v,n ) ; //,o ) ;
	//============================================================

	std::string annBlockNameP =  _toString( LayerId , "BerylliumAnnulusBlock_%02d_posZ") + _toString( (int)AnnulusBlock_loop, "_%02d" ) ;
	std::string annBlockNameN =  _toString( LayerId , "BerylliumAnnulusBlock_%02d_negZ") + _toString( (int)AnnulusBlock_loop, "_%02d" ) ;
	
	vxd.setVisAttributes(lcdd,  "CyanVis" , BerylliumAnnulusBlockLogical) ;
	
	phirot2 = phirot*AnnulusBlock_loop;
	
	RotationZYX rot( 0, phirot2 , (pi*0.5) ) ;
	
	double ZAnnulusBlock2=shell_half_z -(beryllium_ladder_block_length2/2.);// - (shell_thickess/2.); 
	
	pv = supp_assembly.placeVolume( BerylliumAnnulusBlockLogical,  Transform3D( rot, Position((layer_radius+beryllium_ladder_block_thickness+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
											     -(layer_radius+beryllium_ladder_block_thickness+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
											     ZAnnulusBlock2))  );
	DetElement  annBlockPosZ( vxd , annBlockNameP  , x_det.id() );
	annBlockPosZ.setPlacement( pv ) ;
	volSurfaceList( annBlockPosZ )->push_back( surfAnnBlock ) ;

	pv = supp_assembly.placeVolume( BerylliumAnnulusBlockLogical,  Transform3D( rot, Position((layer_radius+beryllium_ladder_block_thickness+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
											     -(layer_radius+beryllium_ladder_block_thickness+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
											     -ZAnnulusBlock2)) ) ;
	DetElement  annBlockNegZ( vxd , annBlockNameN  , x_det.id() );
	annBlockNegZ.setPlacement( pv ) ;
	volSurfaceList( annBlockNegZ )->push_back( surfAnnBlock ) ;
      }
    }
    
    
    //****************************************************************************************
    // *********************************  Electronics   **********************************
    // ******************************  (dead Si layer ends)   ********************************
    //****************************************************************************************
    
    // *********************************  Electronics at the end of the ladder  **********************************
    
    if(end_ladd_electronics_option==1){
      
      Box ElectronicsEndSolid( ladder_width, end_electronics_half_z, electronics_structure_thickness/2. );
      
      Volume ElectronicsEndLogical(_toString(LayerId,"ElectronicsEnd_%02d"),ElectronicsEndSolid, activeMaterial ); //("silicon_2.33gccm")
      
      vxd.setVisAttributes(lcdd,  "GreenVis" , ElectronicsEndLogical );
      
      double end_ladd_electronic_offset_phi = offset_phi +(side_band_electronics_option * side_band_electronics_width/2.);
      
    //====== create the meassurement surface ===================
    Vector3D u( 1. , 0. , 0. ) ;
    Vector3D v( 0. , 1. , 0. ) ;
    Vector3D n( 0. , 0. , 1. ) ;

    double end_ladd_elec_thick = metal_traces_thickness + flex_cable_thickness + foam_spacer_thickness + electronics_structure_thickness;

    VolPlane surfEndElec( ElectronicsEndLogical , SurfaceType(SurfaceType::Helper) , end_ladd_elec_thick/2. ,  end_ladd_elec_thick/2. , u,v,n ) ; //,o ) ;
    //============================================================

      if (LayerId==1||LayerId==3||LayerId==5) {       
	
	for (double elec_loop=0; elec_loop<nb_ladder;elec_loop++) {

	  std::string elecEndLadNameP =  _toString( LayerId , "ElectronicsEnd_%02d_posZ") + _toString( (int)elec_loop, "_%02d" ) ;
	  std::string elecEndLadNameN =  _toString( LayerId , "ElectronicsEnd_%02d_negZ") + _toString( (int)elec_loop, "_%02d" ) ;
	  
	  phirot2 = phirot*elec_loop;
	  RotationZYX rot( 0, phirot2 , (pi*0.5) ) ;    
	  
	  double Z = ladder_length +end_electronics_half_z + (ladder_gap/2.);
	  
	  pv = layer_assembly.placeVolume( ElectronicsEndLogical,  
				     Transform3D( rot, Position((layer_radius+(electronics_structure_thickness/2.)+layer_gap)*sin(phirot2)+ end_ladd_electronic_offset_phi*cos(phirot2),
								-(layer_radius+(electronics_structure_thickness/2.)+layer_gap)*cos(phirot2)+ end_ladd_electronic_offset_phi*sin(phirot2),
								Z))  );

	  DetElement  elecEndLadDEposZ( vxd ,  elecEndLadNameP , x_det.id() );
	  elecEndLadDEposZ.setPlacement( pv ) ;
	  volSurfaceList( elecEndLadDEposZ )->push_back( surfEndElec ) ;
	  
	  pv = layer_assembly.placeVolume( ElectronicsEndLogical, 
				     Transform3D( rot, Position((layer_radius+(electronics_structure_thickness/2.)+layer_gap)*sin(phirot2)+ end_ladd_electronic_offset_phi*cos(phirot2),
								-(layer_radius+(electronics_structure_thickness/2.)+layer_gap)*cos(phirot2)+ end_ladd_electronic_offset_phi*sin(phirot2),
								-Z))  );

	  DetElement  elecEndLadDEnegZ( vxd ,  elecEndLadNameN , x_det.id() );
	  elecEndLadDEnegZ.setPlacement( pv ) ;
	  volSurfaceList( elecEndLadDEnegZ )->push_back( surfEndElec ) ;
	}
	
      } else if (LayerId==0||LayerId==2||LayerId==4)  {       
	
	for (double elec_loop=0; elec_loop<nb_ladder;elec_loop++) {

	  std::string elecEndLadNameP =  _toString( LayerId , "ElectronicsEnd_%02d_posZ") + _toString( (int)elec_loop, "_%02d" ) ;
	  std::string elecEndLadNameN =  _toString( LayerId , "ElectronicsEnd_%02d_negZ") + _toString( (int)elec_loop, "_%02d" ) ;
	  
	  phirot2 = phirot*elec_loop;
	  RotationZYX rot( 0, phirot2 , (pi*0.5) ) ;    
	  
	  double Z = ladder_length +end_electronics_half_z + (ladder_gap/2.);
	  
	  pv = layer_assembly.placeVolume( ElectronicsEndLogical,  
				     Transform3D( rot, Position((layer_radius-(electronics_structure_thickness/2.))*sin(phirot2)+ end_ladd_electronic_offset_phi*cos(phirot2),
								-(layer_radius-(electronics_structure_thickness/2.))*cos(phirot2)+ end_ladd_electronic_offset_phi*sin(phirot2),
								Z))  );

	  DetElement  elecEndLadDEposZ( vxd ,  elecEndLadNameP , x_det.id() );
	  elecEndLadDEposZ.setPlacement( pv ) ;
	  volSurfaceList( elecEndLadDEposZ )->push_back( surfEndElec ) ;
	  
	  pv = layer_assembly.placeVolume( ElectronicsEndLogical,
				     Transform3D( rot, Position((layer_radius-(electronics_structure_thickness/2.))*sin(phirot2)+ end_ladd_electronic_offset_phi*cos(phirot2),
								-(layer_radius-(electronics_structure_thickness/2.))*cos(phirot2)+ end_ladd_electronic_offset_phi*sin(phirot2),
								-Z))  );

	  DetElement  elecEndLadDEnegZ( vxd ,  elecEndLadNameN , x_det.id() );
	  elecEndLadDEnegZ.setPlacement( pv ) ;
	  volSurfaceList( elecEndLadDEnegZ )->push_back( surfEndElec ) ;

	}
      }
    }
    // *********************************  Electronics a long  the ladder  **********************************
    
    if(side_band_electronics_option==1){
      
      Box ElectronicsBandSolid( side_band_electronics_width/2., ladder_length/2., side_band_electronics_thickness/2. );
      
      Volume ElectronicsBandLogical(_toString(LayerId,"ElectronicsBand_%02d"), ElectronicsBandSolid, activeMaterial ) ;
      
      vxd.setVisAttributes(lcdd,  "GreenVis" , ElectronicsBandLogical ) ;
      
      //fixme: turn off sensitive sidebands for now - not sure they cause problems in the volume manager ....
      // active_side_band_electronics_option = 0 ;
      if(active_side_band_electronics_option==1)
	ElectronicsBandLogical.setSensitiveDetector(sens);
      
      
      double side_band_electronic_offset_phi = offset_phi - (side_band_electronics_option * ladder_width);
      
      if (LayerId==1||LayerId==3||LayerId==5) {       
	
	for (double elec_loop=0; elec_loop<nb_ladder;elec_loop++) {   
	  
	  phirot2 = phirot*elec_loop;
	  RotationZYX rot( 0, phirot2 , (pi*0.5) ) ;   
	  
	  double Z = (ladder_length* (1-side_band_electronics_option/2.)) + ladder_gap/2.;
	  
	  // encoder[ILDCellID0::layer]  =  LayerId -1;
	  // encoder[ILDCellID0::module] = elec_loop ;
	  // cellID0 = encoder.lowWord() ;  
	  
	  pv = layer_assembly.placeVolume( ElectronicsBandLogical,
				     Transform3D( rot, Position((layer_radius+(side_band_electronics_thickness/2.)+layer_gap)*sin(phirot2)+side_band_electronic_offset_phi*cos(phirot2),
								-(layer_radius+(side_band_electronics_thickness/2.)+layer_gap)*cos(phirot2)+side_band_electronic_offset_phi*sin(phirot2),
								Z))  ) ;
	  
	  //**fg: choose sensor 1 for sensitive electronics side band
	  if(active_side_band_electronics_option==1)
	    pv.addPhysVolID("layer", LayerId ).addPhysVolID( "module" , int(elec_loop)  ).addPhysVolID("sensor", 1 ).addPhysVolID("side", 1 )   ;

	  pv = layer_assembly.placeVolume( ElectronicsBandLogical,
				     Transform3D( rot, Position((layer_radius+(side_band_electronics_thickness/2.)+layer_gap)*sin(phirot2)+side_band_electronic_offset_phi*cos(phirot2),
								-(layer_radius+(side_band_electronics_thickness/2.)+layer_gap)*cos(phirot2)+side_band_electronic_offset_phi*sin(phirot2),
								-Z))  );

	  if(active_side_band_electronics_option==1)
	    pv.addPhysVolID("layer", LayerId ).addPhysVolID( "module" , int(elec_loop)  ).addPhysVolID("sensor", 1 ).addPhysVolID("side", -1 )   ;
	  
	}

      } else if (LayerId==1||LayerId==3||LayerId==5) {       
	    
	for (double elec_loop=0; elec_loop<nb_ladder;elec_loop++) { 
	  
	  phirot2 = phirot*elec_loop;
	  RotationZYX rot( 0, phirot2 , (pi*0.5) ) ;   

	  double Z = (ladder_length* (1-side_band_electronics_option/2.)) + ladder_gap/2.;
	      
	  // encoder[ILDCellID0::layer]  =  LayerId -1;
	  // encoder[ILDCellID0::module] = elec_loop ;
	  // cellID0 = encoder.lowWord() ;  

	  pv = layer_assembly.placeVolume( ElectronicsBandLogical,
				     Transform3D( rot, Position((layer_radius-(side_band_electronics_thickness/2.))*sin(phirot2)+side_band_electronic_offset_phi*cos(phirot2),
								-(layer_radius-(side_band_electronics_thickness/2.))*cos(phirot2)+side_band_electronic_offset_phi*sin(phirot2),
								Z))  );

	  if(active_side_band_electronics_option==1)
	    pv.addPhysVolID("layer", LayerId ).addPhysVolID( "module" , int(elec_loop)  ).addPhysVolID("sensor", 1 ).addPhysVolID("side", 1 )   ;

	  pv = layer_assembly.placeVolume( ElectronicsBandLogical,
				     Transform3D( rot, Position((layer_radius-(side_band_electronics_thickness/2.))*sin(phirot2)+side_band_electronic_offset_phi*cos(phirot2),
								-(layer_radius-(side_band_electronics_thickness/2.))*cos(phirot2)+side_band_electronic_offset_phi*sin(phirot2),
								-Z))  );
	  if(active_side_band_electronics_option==1)
	    pv.addPhysVolID("layer", LayerId ).addPhysVolID( "module" , int(elec_loop)  ).addPhysVolID("sensor", 1 ).addPhysVolID("side", -1 )   ;
	}
      }      
    }
    
    //****************************************************************************************
    //*******************************  Strip lines (Kapton + metal)  *************************
    //************ here the strip lines are still simulate by conical geometry ***************
    //************ the thickness varies linearly with z **************************************
    //****************************************************************************************
    
    double strip_line_start_z=0.;
    
    if (LayerId==0||LayerId==1) {
      strip_line_start_z = ladder_length + ladder_gap/2. +( end_electronics_half_z * 2.)+ shell_thickess + beryllium_ladder_block_length*2 ; // to avoid overlaps
    } else {
      strip_line_start_z = shell_half_z + shell_endplate_thickness;//ladder_length + ladder_gap/2. - end_electronics_half_z * 2.+ shell_thickess  ; // to avoid overlaps
    }
    //double strip_line_half_z = (strip_line_final_z - strip_line_start_z) / 2.;

    double strip_line_half_z = (dzSty - strip_line_start_z) / 2.;

    // std::cout << " ############## dzSty : " << dzSty << " strip_line_start_z  : " << strip_line_start_z 
    // 	      << "strip_line_half_z "  << strip_line_half_z << std::endl ; 
    //    std::cout << " ############## sPhi : " << sPhi << "  endPhi : " << sPhi+dPhi << std::endl;

    assert (strip_line_half_z>0);

    if (LayerId==0||LayerId==2||LayerId==4) {
      
      //Here we define the solid and logical volumes of the kapton strip lines
      ConeSegment KaptonLinesSolid( strip_line_half_z, 
				    layer_radius, // inside radius at  -fDz
				    layer_radius + initial_kapton_striplines_thickness, // outside radius at -fDz
				    cryostat_apperture + LayerId*final_kapton_striplines_thickness, // inside radius at  +fDz
				    cryostat_apperture + LayerId*final_kapton_striplines_thickness + final_kapton_striplines_thickness,
				    sPhi, sPhi+dPhi);
      
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      //**fg, NB:  TGeoConeSeg( dz,rmin0,rmax0,rmin1,rmax1, phi0, phi1 )  -  G4Cons( rmin0,rmin1,rmax0,rmax1, dz, phi0, delta_phi ) !!!!!!!!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      Volume KaptonLinesLogical( _toString(LayerId,"KaptonLines_%02d"), KaptonLinesSolid, flexCableMaterial ) ;
      
      vxd.setVisAttributes(lcdd,  "WhiteVis" , KaptonLinesLogical );
      
      //Here we define the solid and logical volumes of the metal traces of the strip lines
      ConeSegment MetalLinesSolid( strip_line_half_z,
				   layer_radius + initial_kapton_striplines_thickness, // inside radius at  -fDz
				   layer_radius + initial_kapton_striplines_thickness + initial_metal_striplines_thickness, // outside radius at -fDz
				   cryostat_apperture + LayerId*final_kapton_striplines_thickness + final_kapton_striplines_thickness, // inside radius at  +fDz
				   cryostat_apperture + LayerId*final_kapton_striplines_thickness + final_kapton_striplines_thickness + final_metal_striplines_thickness,
				   sPhi,  sPhi+dPhi);
      
      Volume MetalLinesLogical( _toString(LayerId,"MetalLines_%02d"), MetalLinesSolid, metalTracesMaterial );
      
      vxd.setVisAttributes(lcdd,  "GrayVis" , MetalLinesLogical) ;
      
      //here we place both the kapton and copper part of the strip lines
      double Z = strip_line_start_z + strip_line_half_z;
      
      pv = supp_assembly.placeVolume( KaptonLinesLogical,  Transform3D( RotationZYX() , Position(0., 0., Z) ) );
      pv = supp_assembly.placeVolume( MetalLinesLogical,   Transform3D( RotationZYX() , Position(0., 0., Z) ) );
      
      RotationZYX rot( 0, 0 , pi ) ;   // the same but other side
      
      pv = supp_assembly.placeVolume( KaptonLinesLogical,  Transform3D( rot , Position(0., 0., -Z) ) );
      pv = supp_assembly.placeVolume( MetalLinesLogical,   Transform3D( rot , Position(0., 0., -Z) ) );
    }

    //****************************************************************************************
    //*** Here we place the cabling going outside the VXD ************************************
    //****************************************************************************************
    
    double external_cable_length = (drAlu + drSty)/2.;
    double ExternalCablesZ = dzSty + drSty/2. + drAlu/2. ;

    //** G4Tubs( const G4String& pName, G4double pRMin, G4double pRMax, G4double pDz, G4double pSPhi,G4double pDPhi );
    //** TGeoTubeSeg(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2)

    //kapton part
    Tube ExternalKaptonCablesSolid(cryostat_apperture, cryostat_apperture + 3*external_kapton_thickness/2.,  external_cable_length, sPhi, sPhi+dPhi);
    //The reason for the factor three is that the thickness refer to the thickness of each single cable, and we have three cables in total, one per layer

    Volume ExternalKaptonCablesLogical(_toString(LayerId,"ExternalKaptonCables_%02d"), ExternalKaptonCablesSolid, lcdd.material("G4_KAPTON") );
 
    vxd.setVisAttributes(lcdd,  "WhiteVis" , ExternalKaptonCablesLogical );

    //metal part
    Tube ExternalMetalCablesSolid(cryostat_apperture - 3*external_metal_thickness/2., cryostat_apperture, external_cable_length,  sPhi, sPhi+dPhi);

    Volume ExternalMetalCablesLogical( _toString(LayerId,"ExternalMetalCables_%02d"), ExternalMetalCablesSolid ,  lcdd.material("G4_Cu") ) ;
 
    vxd.setVisAttributes(lcdd,  "GrayVis" , ExternalMetalCablesLogical );


    pv = supp_assembly.placeVolume( ExternalKaptonCablesLogical,  Transform3D( RotationZYX() , Position(0., 0.,  ExternalCablesZ) ) );
    pv = supp_assembly.placeVolume( ExternalKaptonCablesLogical,  Transform3D( RotationZYX() , Position(0., 0., -ExternalCablesZ) ) );

    pv = supp_assembly.placeVolume( ExternalMetalCablesLogical,   Transform3D( RotationZYX() , Position(0., 0.,  ExternalCablesZ) ) );
    pv = supp_assembly.placeVolume( ExternalMetalCablesLogical,   Transform3D( RotationZYX() , Position(0., 0., -ExternalCablesZ) ) );

    //****************************************************************************************
    //*******************************  Cooling Pipes (Titanium )  ********************************
    //****************************************************************************************      
    
    //endplate cooling pipes
    
    double ZEndPlateCoolPipes = shell_half_z + shell_endplate_thickness;

    double ZEndPlateCoolPipesL1  = ladder_length +  ((end_electronics_half_z*end_ladd_electronics_option) * 2) + shell_thickess + (beryllium_ladder_block_length*2) ;

    //**TGeoTorus(double r, double rmin, double rmax, double phi, double delta_phi)    
    //**G4Torus(const G4String &pName, G4double pRmin, G4double pRmax, G4double pRtor, G4double pSPhi, G4double pDPhi )

    Torus CoolPipeSolid(layer_radius + layer_gap + cool_pipe_outer_radius/2., cool_pipe_inner_radius,
			cool_pipe_outer_radius, sPhi, dPhi);
    
    Volume CoolPipeLogical(_toString(LayerId,"CoolPipe_%02d"), CoolPipeSolid, lcdd.material("G4_Ti")) ; //titanium") ) ;
    
    vxd.setVisAttributes(lcdd,  "MagentaVis" , CoolPipeLogical );

    //one cooling pipe for each double layer
 
    if (LayerId==3 || LayerId==5) {

      pv = supp_assembly.placeVolume( CoolPipeLogical, Transform3D( RotationZYX() , Position(0., 0.,   ZEndPlateCoolPipes+cool_pipe_outer_radius  )) ) ; 
      pv = supp_assembly.placeVolume( CoolPipeLogical, Transform3D( RotationZYX() , Position(0., 0., -(ZEndPlateCoolPipes+cool_pipe_outer_radius) )) ) ; 
      
    } else if (LayerId==1)  { 
      
      pv = supp_assembly.placeVolume( CoolPipeLogical, Transform3D( RotationZYX() , Position(0., 0.,   ZEndPlateCoolPipesL1 + cool_pipe_outer_radius + shell_thickess) )) ;
      pv = supp_assembly.placeVolume( CoolPipeLogical, Transform3D( RotationZYX() , Position(0., 0., -(ZEndPlateCoolPipesL1 + cool_pipe_outer_radius + shell_thickess)) ));
    }

    //***************************************************************************************************************
    // *** cooling pipe connecting the pipes at the central be support endplate to the layer 1 support endplate  ****
    //***************************************************************************************************************

    if (LayerId==1){

      double thetaTube = atan((support_endplate_inner_radious - (layer_radius + layer_gap + 2*cool_pipe_outer_radius)) / (shell_half_z - ZEndPlateCoolPipesL1)) ;

      std::cout << "############## thetaTube = " << thetaTube << std::endl ;

      double CoolPipeLength = (shell_half_z - shell_thickess/2.) - ZEndPlateCoolPipesL1;
      
      Tube CoolPipeTubeSolid( cool_pipe_inner_radius, cool_pipe_outer_radius, CoolPipeLength/2., sPhi, sPhi+dPhi );
      
      Volume CoolPipeTubeLogical( _toString(LayerId,"CoolPipeTube_%02d"), CoolPipeTubeSolid, lcdd.material("G4_Ti")) ; //titanium") ) ; 
      
      vxd.setVisAttributes(lcdd,  "MagentaVis" , CoolPipeTubeLogical );
      
      //**fg: reversed sign compared to original VXD04.cc as this created obvious overlaps in TGeo (also in Mokka ???)
      RotationZYX rm( 0., 0., -thetaTube);
      RotationZYX rm2(0., 0.,  thetaTube);
      
      pv = supp_assembly.placeVolume( CoolPipeTubeLogical, Transform3D( rm,  Position( 0., 
										  (layer_radius + layer_gap + support_endplate_inner_radious)/2.,
										  (ZEndPlateCoolPipesL1 + 3*cool_pipe_outer_radius + CoolPipeLength/2.)) )) ;
      
      pv = supp_assembly.placeVolume( CoolPipeTubeLogical, Transform3D( rm2, Position( 0., 
										  -(layer_radius + layer_gap + support_endplate_inner_radious)/2.,
										  (ZEndPlateCoolPipesL1 + 3*cool_pipe_outer_radius + CoolPipeLength/2.)) )) ;
      
      pv = supp_assembly.placeVolume( CoolPipeTubeLogical, Transform3D( rm2, Position( 0., 
										  (layer_radius + layer_gap + support_endplate_inner_radious)/2.,
										  -(ZEndPlateCoolPipesL1 + 3*cool_pipe_outer_radius + CoolPipeLength/2.)) )) ;
      
      pv = supp_assembly.placeVolume( CoolPipeTubeLogical, Transform3D( rm,  Position( 0., 
										  -(layer_radius + layer_gap + support_endplate_inner_radious)/2.,
										  -(ZEndPlateCoolPipesL1 + 3*cool_pipe_outer_radius + CoolPipeLength/2.)) )) ;
      
    }
      
    //****************************************************************************************
    // *******************************  Si Active layer  *************************************
    //****************************************************************************************

    Box SiActiveLayerSolid( ladder_width, ladder_length/2., active_silicon_thickness/2. );
      
    Volume SiActiveLayerLogical( _toString( LayerId , "SiActiveLayer_%02d"),SiActiveLayerSolid, activeMaterial ) ;

    vxd.setVisAttributes(lcdd,  "BlueVis" , SiActiveLayerLogical ) ;

    SiActiveLayerLogical.setSensitiveDetector(sens);
   

    //====== create the meassurement surface ===================
    Vector3D u( 1. , 0. , 0. ) ;
    Vector3D v( 0. , 1. , 0. ) ;
    Vector3D n( 0. , 0. , 1. ) ;
    //    Vector3D o( 0. , 0. , 0. ) ;

    double supp_thick = metal_traces_thickness + flex_cable_thickness + foam_spacer_thickness;

    //fg: set inner and outer thickness to sens/2+sup as  every other layer is flipped -> will just add a bit of air ...
    VolPlane surf( SiActiveLayerLogical , SurfaceType(SurfaceType::Sensitive) , active_silicon_thickness/2 + supp_thick  , active_silicon_thickness/2 + supp_thick , u,v,n ) ; //,o ) ;
    //============================================================


    double active_offset_phi = offset_phi +(side_band_electronics_option * side_band_electronics_width/2.); 
      
    for (double active_loop=0;active_loop<nb_ladder;active_loop++){
	
      phirot2 =  phirot*active_loop;
      RotationZYX rot( 0. , phirot2, (pi*0.5) ) ;
	
      double Z = ladder_length/2.+ ladder_gap;
	
      // encoder[ILDCellID0::layer]  =  LayerId -1;
      // encoder[ILDCellID0::module] = active_loop ;
      // cellID0 = encoder.lowWord() ;  

      std::string ladderNameP =  _toString( LayerId , "SiActiveLayer_%02d_posZ") + _toString( (int)active_loop, "_%02d" ) ;
      std::string ladderNameN =  _toString( LayerId , "SiActiveLayer_%02d_negZ") + _toString( (int)active_loop, "_%02d" ) ;


      if (LayerId==1 ) {

	std::cout << "  ***** place ladder "  << int(active_loop) << " for " << _toString( LayerId , "SiActiveLayer_%02d") << std::endl ;
      }

      if (LayerId==1||LayerId==3||LayerId==5) {
	
	pv = layer_assembly.placeVolume( SiActiveLayerLogical,  Transform3D( rot, Position((layer_radius+(active_silicon_thickness/2.)+layer_gap)*sin(phirot2)+active_offset_phi*cos(phirot2),
										     -(layer_radius+(active_silicon_thickness/2.)+layer_gap)*cos(phirot2)+active_offset_phi*sin(phirot2),
										     Z)) ) ;
	
	pv.addPhysVolID("layer", LayerId ).addPhysVolID( "module" , int(active_loop)  ).addPhysVolID("sensor", 0 ).addPhysVolID("side", 1 )   ;

	DetElement   ladderDEposZ( vxd ,  ladderNameP , x_det.id() );
	ladderDEposZ.setPlacement( pv ) ;
	volSurfaceList( ladderDEposZ )->push_back( surf ) ;


	pv = layer_assembly.placeVolume( SiActiveLayerLogical,  Transform3D( rot, Position((layer_radius+(active_silicon_thickness/2.)+layer_gap)*sin(phirot2)+active_offset_phi*cos(phirot2),
										     -(layer_radius+(active_silicon_thickness/2.)+layer_gap)*cos(phirot2)+active_offset_phi*sin(phirot2),
										     -Z)) );
	
	pv.addPhysVolID("layer", LayerId ).addPhysVolID( "module" , int(active_loop)  ).addPhysVolID("sensor", 0 ).addPhysVolID("side", -1 )   ;

	DetElement   ladderDEnegZ( vxd ,   ladderNameN , x_det.id() );
	ladderDEnegZ.setPlacement( pv ) ;
	volSurfaceList( ladderDEnegZ )->push_back( surf ) ;
	

      } else if (LayerId==0||LayerId==2||LayerId==4) { 

	pv = layer_assembly.placeVolume( SiActiveLayerLogical,  Transform3D( rot, Position((layer_radius-(active_silicon_thickness/2.))*sin(phirot2)+active_offset_phi*cos(phirot2),
										     -(layer_radius-(active_silicon_thickness/2.))*cos(phirot2)+active_offset_phi*sin(phirot2),
										     Z)) ) ;
	
	pv.addPhysVolID("layer", LayerId ).addPhysVolID( "module" , int(active_loop)  ).addPhysVolID("sensor", 0 ).addPhysVolID("side", 1 )   ;
	
	DetElement   ladderDEposZ( vxd ,  ladderNameP , x_det.id() );
	ladderDEposZ.setPlacement( pv ) ;
	volSurfaceList( ladderDEposZ )->push_back( surf ) ;

	
	pv = layer_assembly.placeVolume( SiActiveLayerLogical,  Transform3D( rot, Position((layer_radius-(active_silicon_thickness/2.))*sin(phirot2)+active_offset_phi*cos(phirot2),
	 									     -(layer_radius-(active_silicon_thickness/2.))*cos(phirot2)+active_offset_phi*sin(phirot2),
	 									     -Z)) );
	DetElement   ladderDEnegZ( vxd ,  ladderNameN , x_det.id() );
	ladderDEnegZ.setPlacement( pv ) ;
	volSurfaceList( ladderDEnegZ )->push_back( surf ) ;
	
	pv.addPhysVolID("layer", LayerId ).addPhysVolID( "module" , int(active_loop)  ).addPhysVolID("sensor", 0 ).addPhysVolID("side", -1 )   ;
	
      }		  

    }
      
    // #ifdef MOKKA_GEAR
    //       // sensitive layer
    //       helpLayer thisSens ;
    //       if (LayerId==2||LayerId==4||LayerId==6) 
    // 	{ 
    // 	  thisSens.distance  = layer_gap + layer_radius;
    // 	}
    //       if (LayerId==1||LayerId==3||LayerId==5) 
    // 	{ 
    // 	  thisSens.distance  = layer_radius  - active_silicon_thickness ;
    // 	}
    //       thisSens.offset    = active_offset_phi ;
    //       thisSens.thickness = active_silicon_thickness ;
    //       thisSens.length    = ladder_length ;
    //       if (active_side_band_electronics_option==1) {
    // 	thisSens.width     = ladder_width*2.+side_band_electronics_width ;
    //       }
    //       else  {
    // 	thisSens.width     = ladder_width*2.; 
    //       }
    //       thisSens.radLength = (SiActiveLayerLogical->GetMaterial())->GetRadLen()/mm ;
      
    //       // save information for gear
    //       gearHelpLadders.push_back( thisLadder );
    //       gearHelpSensitives.push_back( thisSens ) ;
    //       gearHelpNumberLadders.push_back( (int) nb_ladder ) ;
      
    //       // fg: here we start with the first ladder at -pi/2 (i.e. the negative y-axis)
    //       gearHelpPhi0.push_back( -pi/2. ) ;
      
    //       gearHelpGap = std::max( gearHelpGap , ladder_gap ) ;
    //       gearHelpCount ++ ;
    // #endif


    //--- fill DDRec::ZPlanarData

    if (LayerId==1||LayerId==3||LayerId==5)  { 

      thisLayer.distanceSensitive  = layer_gap + layer_radius;

    } else if (LayerId==0||LayerId==2||LayerId==4) { 

      thisLayer.distanceSensitive  = layer_radius  - active_silicon_thickness ;
    }
    //Adjusted by Thorben Quast to provide consistency with CLIC while drawing for CED
    //Please check as the prefactor is simply chosen for the geometry to match the hits
    thisLayer.offsetSensitive    = 0.5*active_offset_phi ;
    thisLayer.thicknessSensitive = active_silicon_thickness ;
    thisLayer.zHalfSensitive    = ladder_length ;

    if (active_side_band_electronics_option==1) {

      thisLayer.widthSensitive     = ladder_width*2.+side_band_electronics_width ;

    } else  {

      thisLayer.widthSensitive      = ladder_width*2.; 
    }

    thisLayer.ladderNumber = (int) nb_ladder  ;
    thisLayer.phi0 =  -pi/2.  ;
    
    zPlanarData->layers.push_back( thisLayer ) ;

    // gearHelpGap = std::max( gearHelpGap , ladder_gap ) ;
    // gearHelpCount ++ ;



  } // --- end loop over layers ----------------------------------------------------------------------------------------


  
  //****************************************
  // Outer support shell
  //****************************************
  
  // ************central tube************
  
  Tube SupportShellSolid( shell_inner_radious, shell_inner_radious+shell_thickess, shell_half_z, sPhi, sPhi+dPhi );

  Volume SupportShellLogical("SupportShell", SupportShellSolid,  lcdd.material("G4_Be")) ; //"beryllium") ) ;

  vxd.setVisAttributes(lcdd,  "CyanVis" , SupportShellLogical ) ;
  

  //====== create a material surface for the support shell ===================
  
  Vector3D osupshell( shell_inner_radious+shell_thickess/2. , 0. , 0.  ) ;
  
  VolCylinder surfSupShell( SupportShellLogical , SurfaceType(SurfaceType::Helper) , shell_thickess/2. , shell_thickess/2., osupshell ) ;
  
  volSurfaceList( suppDE )->push_back( surfSupShell ) ;
  //============================================================



  pv = supp_assembly.placeVolume( SupportShellLogical ) ;
  //  pv = supp_assembly.placeVolume( SupportShellLogical, Transform3D( RotationZYX(), Position() ) ) ;
  

  // ************support endplates************

  double support_endplate_half_z = shell_endplate_thickness/2;
  
  Tube EndPlateShellSolid( support_endplate_inner_radious,  shell_inner_radious+shell_thickess, support_endplate_half_z, sPhi, sPhi+dPhi ) ;

  Volume EndPlateShellLogical("EndPlateShell_outer", EndPlateShellSolid,  lcdd.material("G4_Be")) ; //"beryllium") ) ;

  vxd.setVisAttributes(lcdd,  "CyanVis" , EndPlateShellLogical ) ;
  
  double ZEndPlateShell = shell_half_z + shell_endplate_thickness/2.;// + (beryllium_ladder_block_length*2);


  DetElement endplateFwdDE( suppDE , name+"_endplate_fwd" , x_det.id() )  ;
  DetElement endplateBwdDE( suppDE , name+"_endplate_bwd" , x_det.id() )  ;

  pv = supp_assembly.placeVolume( EndPlateShellLogical, Transform3D( RotationZYX(), Position(0., 0.,  ZEndPlateShell ) ) ) ;
  endplateFwdDE.setPlacement( pv ) ;
  pv = supp_assembly.placeVolume( EndPlateShellLogical, Transform3D( RotationZYX(), Position(0., 0., -ZEndPlateShell ) ) ) ;
  endplateBwdDE.setPlacement( pv ) ;

  // --- add a helper surface for the outer part of the endplate shell ---------------------------------

  Vector3D up( 1. , 0. , 0. ) ;
  Vector3D vp( 0. , 1. , 0. ) ;
  Vector3D np( 0. , 0. , 1. ) ;

  // need to set the origin of this helper plane to be inside the material ( otherwise it would pick up the vacuum at the origin)
  Vector3D o_endplate( 0. ,   0.5 * (  support_endplate_inner_radious + shell_inner_radious+shell_thickess )    , 0. ) ;

  VolPlane surfEndplate( EndPlateShellLogical , SurfaceType(SurfaceType::Helper) , support_endplate_half_z , support_endplate_half_z, up,vp,np, o_endplate ) ;
  volSurfaceList( endplateFwdDE )->push_back( surfEndplate ) ;
  volSurfaceList( endplateBwdDE )->push_back( surfEndplate ) ;

  //-----------------------------------------------------------------------------------------------------
    


  
  // ************support endplates for the layer 1************
  
  double support_endplate_half_z_L1 = shell_thickess/2;
  double ladder_length = ladder_0_length ;
  
  Tube EndPlateShellSolidL1( support_endplate_inner_radious_L1, support_endplate_outer_radious_L1, support_endplate_half_z_L1, sPhi, sPhi+dPhi ) ;
  
  Volume EndPlateShellLogicalL1("EndPlateShell_inner", EndPlateShellSolidL1,  lcdd.material("G4_Be")) ; //"beryllium") ) ;
  
  vxd.setVisAttributes(lcdd,  "CyanVis" , EndPlateShellLogicalL1 ) ;
  
  double ZEndPlateShell2 = ladder_length +  ((end_electronics_half_z*end_ladd_electronics_option) * 2) + shell_thickess/2. + (beryllium_ladder_block_length*2) ;
  
  pv = supp_assembly.placeVolume( EndPlateShellLogicalL1, Transform3D( RotationZYX(), Position(0., 0.,   ZEndPlateShell2 ) ) ) ;
  pv = supp_assembly.placeVolume( EndPlateShellLogicalL1, Transform3D( RotationZYX(), Position(0., 0.,  -ZEndPlateShell2 ) ) ) ;
  
  //**** beryllium support shell cone ************************************************

  double support_cone_half_z = (shell_half_z - (ZEndPlateShell2 + shell_thickess/2.))/2.;
  
  ConeSegment SupportShellCone( support_cone_half_z, 
				support_endplate_outer_radious_L1, support_endplate_outer_radious_L1 + shell_thickess, 
				support_endplate_inner_radious,    support_endplate_inner_radious    + shell_thickess,
				sPhi, sPhi+dPhi);
  
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  //**fg, NB:  TGeoConeSeg( dz,rmin0,rmax0,rmin1,rmax1, phi0, phi1 )  -  G4Cons( rmin0,rmin1,rmax0,rmax1, dz, phi0, delta_phi ) !!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
  Volume SupportConeLogical("SupportCone", SupportShellCone,  lcdd.material("G4_Be")) ; //"beryllium") ) ;
  
  vxd.setVisAttributes(lcdd,  "CyanVis" , SupportConeLogical ) ;
  
  double ZCone = ladder_length +  ((end_electronics_half_z*end_ladd_electronics_option) * 2) + shell_thickess + (beryllium_ladder_block_length*2) + support_cone_half_z;

  // ================================== adding tracking surfaces for beryllium support cone ====================================
  const double shell_theta = atan2( -1.*support_endplate_outer_radious_L1 + support_endplate_inner_radious, 2.*support_cone_half_z);  

  Vector3D o_shellcon(0.5 * (support_endplate_outer_radious_L1 + shell_thickess + support_endplate_inner_radious) , 0., 0. ) ;
  Vector3D shellcone_angle (  1. , 0. , shell_theta, Vector3D::spherical ) ;

  VolCone suppShellCone( SupportConeLogical, SurfaceType( SurfaceType::Helper ) , 0.5*shell_thickess  , 0.5*shell_thickess , shellcone_angle, o_shellcon );
  // ==============================================================================================================================
  
  pv = supp_assembly.placeVolume( SupportConeLogical, Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,  ZCone ) ) ) ;
  DetElement suppShellConePosDE ( suppDE , "posShellCone" , x_det.id() )  ;
  suppShellConePosDE.setPlacement( pv );
  volSurfaceList( suppShellConePosDE )->push_back( suppShellCone );

  pv = supp_assembly.placeVolume( SupportConeLogical, Transform3D( RotationZYX( 0, 0, pi ), Position(0., 0., -ZCone ) ) ) ;
  DetElement suppShellConeNegDE ( suppDE , "negShellCone" , x_det.id() )  ;
  suppShellConeNegDE.setPlacement( pv );
  volSurfaceList( suppShellConeNegDE )->push_back( suppShellCone );

  //*** beryllium support forward part **************************************************************************
  
  double supportForZ = shell_half_z + shell_endplate_thickness + forward_shell_half_z;
  
  Tube SupportForSolid( support_endplate_inner_radious, support_endplate_inner_radious + shell_endplate_thickness, forward_shell_half_z,  sPhi, sPhi+dPhi ) ;

  Volume SupportForLogical("SupportFor", SupportForSolid,  lcdd.material("G4_Be")) ; //"beryllium") ) ;

  vxd.setVisAttributes(lcdd,  "CyanVis" , SupportForLogical ) ;

  pv = supp_assembly.placeVolume( SupportForLogical, Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,  supportForZ ) ) ) ;
  pv = supp_assembly.placeVolume( SupportForLogical, Transform3D( RotationZYX( 0, 0, pi ), Position(0., 0., -supportForZ ) ) ) ;

  
  //*** Cryostat ***************************************************************


  double aluEndcapZ = dzSty + drSty + drAlu / 2;
  double styEndcapZ = dzSty + drSty / 2;
  
  double aluHalfZ = dzSty + drSty;
  
  if (useCryo) {
    Material aluMaterial = lcdd.material( "G4_Al" ) ;
    //          VisAttributes *aluVisAttributes = new VisAttributes(Colour(0.5, 0.5, 0.5)); 
    
    Material styMaterial = lcdd.material("styropor");
    


    Tube aluBarrelSolid( rAlu, rAlu + drAlu, aluHalfZ, sPhi, sPhi+dPhi);
    Volume aluBarrelLog( "CryostatAluSkinBarrel", aluBarrelSolid, aluMaterial );
    vxd.setVisAttributes(lcdd,  "GrayVis" , aluBarrelLog ) ;
    pv = supp_assembly.placeVolume( aluBarrelLog ) ;
    
    Tube styBarrelSolid(  rSty, rSty + drSty, dzSty, sPhi, sPhi+dPhi);
    Volume styBarrelLog( "CryostatFoamBarrel", styBarrelSolid, styMaterial );
    vxd.setVisAttributes(lcdd,  "LightGrayVis", styBarrelLog ) ;
    pv = supp_assembly.placeVolume( styBarrelLog ) ;

    //====== create a material surface for the cryostat barrel ===================

    double rc =  ( rAlu + drAlu /2.) ;
    Vector3D oc( rc , 0. , 0.  ) ;

    double outer_thick = drAlu/2. ;
    double inner_thick = drAlu/2. + drSty ;

    VolCylinder surfC( aluBarrelLog , SurfaceType(SurfaceType::Helper) , inner_thick , outer_thick, oc ) ;

    volSurfaceList( suppDE )->push_back( surfC ) ;
    //============================================================


    //Aluminium + styropor endplates for the cryostat
    //Create an apperture at the cryostat endcap for the cabling and the cooling pipes
    

    // in order to assign material surfaces to the endcap of the cryostat
    // we need four DetElements as the endcap consists of four distinct volumes ...
    DetElement suppFwdInDE ( suppDE , name+"_support_fwd_inner" , x_det.id() )  ;
    DetElement suppBwdInDE ( suppDE , name+"_support_bwd_inner" , x_det.id() )  ;
    DetElement suppFwdOutDE( suppDE , name+"_support_fwd_outer" , x_det.id() )  ;
    DetElement suppBwdOutDE( suppDE , name+"_support_bwd_outer" , x_det.id() )  ;

    Tube aluEndcapSolidInner(  rInner, cryostat_apperture - cryostat_apperture_radius, drAlu / 2, sPhi, sPhi+dPhi);
    Volume aluEndcapInnerLog( "CryostatAluSkinEndPlateInner", aluEndcapSolidInner, aluMaterial );
    vxd.setVisAttributes(lcdd,  "GrayVis" , aluEndcapInnerLog ) ;
    pv = supp_assembly.placeVolume( aluEndcapInnerLog , Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,   aluEndcapZ ) ) ) ;
    suppFwdInDE.setPlacement( pv ) ;
    pv = supp_assembly.placeVolume( aluEndcapInnerLog , Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,  -aluEndcapZ ) ) ) ;
    suppBwdInDE.setPlacement( pv ) ;

    Tube aluEndcapSolidOuter(   cryostat_apperture + cryostat_apperture_radius, rAlu + drAlu, drAlu / 2, sPhi, sPhi+dPhi);
    Volume aluEndcapOuterLog( "CryostatAluSkinEndPlateOuter", aluEndcapSolidOuter, aluMaterial );
    vxd.setVisAttributes(lcdd,  "GrayVis" , aluEndcapOuterLog ) ;
    pv = supp_assembly.placeVolume( aluEndcapOuterLog , Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,   aluEndcapZ ) ) ) ;
    suppFwdOutDE.setPlacement( pv ) ;
    pv = supp_assembly.placeVolume( aluEndcapOuterLog , Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,  -aluEndcapZ ) ) ) ;
    suppBwdOutDE.setPlacement( pv ) ;
    
    Tube styEndcapSolidInner(  rInner, cryostat_apperture - cryostat_apperture_radius, drSty / 2, sPhi, sPhi+dPhi);
    Volume styEndcapInnerLog( "CryostatFoamEndPlateInner", styEndcapSolidInner, styMaterial );
    vxd.setVisAttributes(lcdd,  "WhiteVis" , styEndcapInnerLog ) ;
    pv = supp_assembly.placeVolume( styEndcapInnerLog , Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,   styEndcapZ ) ) ) ;
    pv = supp_assembly.placeVolume( styEndcapInnerLog , Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,  -styEndcapZ ) ) ) ;

    Tube styEndcapSolidOuter(  cryostat_apperture + cryostat_apperture_radius, rSty + drSty, drSty / 2, sPhi, sPhi+dPhi);
    Volume styEndcapOuterLog( "CryostatFoamEndPlateOuter", styEndcapSolidOuter, styMaterial );
    vxd.setVisAttributes(lcdd,  "WhiteVis" , styEndcapOuterLog ) ;
    pv = supp_assembly.placeVolume( styEndcapOuterLog , Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,   styEndcapZ ) ) ) ;
    pv = supp_assembly.placeVolume( styEndcapOuterLog , Transform3D( RotationZYX( 0, 0, 0  ), Position(0., 0.,  -styEndcapZ ) ) ) ;

    //====== create a material surface for the cryostat endcap ===================
    // Vector3D up( 1. , 0. , 0. ) ;
    // Vector3D vp( 0. , 1. , 0. ) ;
    // Vector3D np( 0. , 0. , 1. ) ;

    // need to set the origin of this helper plane to be inside the material ( otherwise it would pick up the vacuum at the origin)
    double mid_r = 0.5 * ( cryostat_apperture - cryostat_apperture_radius +  rInner  ) ;
    Vector3D op_i( 0. , mid_r , 0. ) ;

    VolPlane surfPi( aluEndcapInnerLog , SurfaceType(SurfaceType::Helper) , inner_thick , outer_thick, up,vp,np, op_i ) ;
    volSurfaceList( suppFwdInDE )->push_back( surfPi ) ;
    volSurfaceList( suppBwdInDE )->push_back( surfPi ) ;
    
    mid_r = 0.5 * ( cryostat_apperture + cryostat_apperture_radius + rSty + drSty ) ;
    Vector3D op_o( 0. , mid_r , 0. ) ;

    VolPlane surfPo( aluEndcapOuterLog , SurfaceType(SurfaceType::Helper) , inner_thick , outer_thick, up,vp,np, op_o) ;
    volSurfaceList( suppFwdOutDE )->push_back( surfPo ) ;
    volSurfaceList( suppBwdOutDE )->push_back( surfPo ) ;

    //============================================================

  }
  
  //     cout<<"rAlu "<< rAlu<<endl;
  //     cout<<"drAlu "<< drAlu<<endl;
  //     cout<<"aluHalfZ "<< aluHalfZ<<endl;
  //     cout<<"drSty "<< drSty<<endl;
  //     cout<<"rInner "<<rInner <<endl;
  //     cout<<"+aluEndcapZ "<<+aluEndcapZ <<endl;
  //     cout<<"shell_inner_radious "<<shell_inner_radious <<endl; 
  //     cout<<"foam inner radius "<<rSty <<endl; 
  //     cout << "database name =" << dbName << endl;
  
  // #ifdef MOKKA_GEAR
  //     // -------write data to gear
  
  //     // get gear manager
  //     MokkaGear* gearMgr = MokkaGear::getMgr() ;
  
  //     // construct VXDParameters
  //     gear::ZPlanarParametersImpl* vxdParams = 
  //       new gear::ZPlanarParametersImpl(gearHelpType ,                                        // vxd type 
  // 				      shell_inner_radious ,                                 // inner radius
  // 				      shell_inner_radious+shell_thickess ,                  // outer radius
  // 				      shell_half_z ,                                        // half length
  // 				      gearHelpGap ,                                         // shell gap
  // 				      (SupportShellLogical->GetMaterial())->GetRadLen()/mm ) ; // shell rad length
  
  //     // add all layers
  //     for( int i = 0 ; i < gearHelpCount ; i++ ) {
  //       vxdParams->addLayer( gearHelpNumberLadders[i] , gearHelpPhi0[i] ,
  // 			   gearHelpLadders[i].distance , gearHelpLadders[i].offset,gearHelpLadders[i].thickness ,
  // 			   gearHelpLadders[i].length , gearHelpLadders[i].width , gearHelpLadders[i].radLength ,
  // 			   gearHelpSensitives[i].distance, gearHelpSensitives[i].offset , gearHelpSensitives[i]. thickness , 
  // 			   gearHelpSensitives[i].length , gearHelpSensitives[i].width , gearHelpSensitives[i].radLength ) ;
  //       gearMgr->setVXDParameters( vxdParams ) ;
  //     }
  
  // #endif
  
  zPlanarData->rInnerShell = shell_inner_radious  ;
  zPlanarData->rOuterShell = shell_inner_radious+shell_thickess ;
  zPlanarData->zHalfShell  = shell_half_z ;
  zPlanarData->gapShell    = 0. ;
  //######################################################################################################################################################################
  
  vxd.addExtension< ZPlanarData >( zPlanarData ) ;

  
  //--------------------------------------
  
  
  vxd.setVisAttributes( lcdd, x_det.visStr(), envelope );

  return vxd;
}
DECLARE_DETELEMENT(VXD04,create_element)

//**************************************************************************
// \file    Grainita_ECAL_Barrel_v02.cpp
// \brief:  Grainita ECAL full detector geometry, with pointing structure
// \author: Fangyi Guo (fangyi.guo@cern.ch) start from DD4HEP tutorial in 
//          DRD Calo (https://github.com/DRDCalo/DD4hepTutorials)
// \date:   April 2026
//**************************************************************************

// Includers from DD4hep
#include "DDRec/Vector3D.h"
#include <DD4hep/DetFactoryHelper.h>

using namespace dd4hep;

// Build geometry
//
static Ref_t create_detector(Detector &description, xml_h e,
                             SensitiveDetector sens) {
  std::cout << "--> Grainita_ECAL_Barrel_v02::create_detector() start" << std::endl;

  // Get info from the xml file
  //
  sens.setType("calorimeter");
  xml_det_t x_det = e;
  std::string det_name = x_det.nameStr();
  std::cout << "--> Going to create " << det_name << ", with ID: " << x_det.id() << std::endl;

  //xml_dim_t x_dim = x_det.dimensions();
  //const double CaloX = x_dim.x();
  //const double CaloY = x_dim.y();
  //const double CaloZ = x_dim.z();

  //Define readout
  //Readout readout = sens.readout();
  //Segmentation seg = readout.segmentation();


  //Define material
  Material  air(description.material("Air"));
  Material  MatCarbonfiber(description.material("CarbonFiber"));
  Material  MatCrystal(description.material( x_det.attr<std::string>(_U(material)) ));
  //Material  MatWLSfiber(description.material("WLSFiber"));
  //Material  MatPMMAcladding(description.material("PMMA"));


  // Detector global size
  auto innerR = description.constant<double>("InnerRadius");
  auto crystal_thick = description.constant<double>("CrystalThickness");
  auto back_space = description.constant<double>("BackSpace");
  auto halfZ = description.constant<double>("HalfZ");
  auto nSec_phi = description.constant<int>("NumberOfPhiSectors");
  auto nModule_z = description.constant<int>("NumberOfZModules");

  auto Cframe_thick = description.constant<double>("FrameThickness");
  auto Cseg_thick = description.constant<double>("SegThickness");
  auto module_tilt = description.constant<double>("ModuleTiltAngle"); //Unit in degree
  double module_tilt_rad = module_tilt*M_PI/180.;

  //auto fiber_r = description.constant<double>("FiberRadius");
  //auto cladding_thick = description.constant<double>("CladdingThickness");
  //auto fiber_pitch = description.constant<double>("FiberPitchSize"); 


  double dphi_sec = 2*M_PI/nSec_phi;
  double outerR = (innerR + crystal_thick + back_space + 2*Cframe_thick)/cos(dphi_sec/2.); 
  std::cout<<"  -- Dimension: inner R "<<innerR<<", outer R "<<outerR<<", inner half Z "<<halfZ<<std::endl;
  std::cout<<"     Z segmentation: "<<nModule_z<<", Phi segmentation: "<<nSec_phi<<std::endl;

  //Calculate the shift from phi direction tilt. 
  // Note: The tilt works on the outer vertex.
  // The section of sector in phi-plane is defined with 4 points (in R-theta coordinate):
  //    Inner: A(Rin, -dphi_sec/2.),  D(Rin, dphi_sec/2.)
  //    Original outer withouth tilt: (Rout, -dphi_sec/2.), (Rout, dphi_sec/2.)
  //    Outer with tilt: B(Rout, -dphi_sec/2. + tilt), C(Rout, dphi_sec/2. + tilt)
  // BUT this gives non-parallel inner and outer surface. We force the outer surface be parallel with inner surface
  // in the plane of x=Rout*cos(dphi_sec/2.). 
  //    (in x-y coordinate)
  //    B' (Rout*cos(dphi_sec/2.), y_B'), C'(Rout*cos(dphi_sec/2.), y_C')
  // where B' and C' are from extending AB and DC to x=Rout*cos(dphi_sec/2.) plane. 
  // To get y_B' and y_C', we calculated the function of AB and DC. 
  //
  //        \
  //     \   \ /|C'
  //      \  // |
  //       \//  |
  //      D |   |
  //        |   |
  //       A --- B'

  //Slope and intercept of line AB
  double tmp_slope1 = -(outerR*sin(dphi_sec/2.-module_tilt_rad) - innerR*tan(dphi_sec/2.)) / 
                      (outerR*cos(dphi_sec/2.-module_tilt_rad) - innerR);
  double tmp_intercept1 = -innerR*tan(dphi_sec/2) - innerR*tmp_slope1;

  //Slope and intercept of line DC
  double tmp_slope2 = (outerR*sin(dphi_sec/2+module_tilt_rad) - innerR*tan(dphi_sec/2)) /
                      (outerR*cos(dphi_sec/2+module_tilt_rad) - innerR);
  double tmp_intercept2 = innerR*tan(dphi_sec/2) - innerR*tmp_slope2;

  //Calculate the shift regarding to the non-tilt point (Rout, +-dphi_sec/2.).  
  double shift_phi_neg = fabs( fabs(tmp_slope1*outerR*cos(dphi_sec/2.)+tmp_intercept1 ) - outerR*sin(dphi_sec/2.) );
  double shift_phi_pos = fabs( fabs(tmp_slope2*outerR*cos(dphi_sec/2.)+tmp_intercept2 ) - outerR*sin(dphi_sec/2.) );

  std::cout<<"-- Consider a tilt angle in phi to avoid pointing cracks. Tilt angle (in deg): "<<module_tilt<<", in rad: "<<module_tilt_rad<<std::endl;
  std::cout<<"Input vars: Rin "<<innerR<<", Rout "<<outerR<<", phi "<<dphi_sec/2.<<", delta "<<module_tilt_rad<<std::endl;
  std::cout<<"  Calculate the new boundary lines: "<<std::endl;
  std::cout<<"    in -y side: k = "<<tmp_slope1<<", b = "<<tmp_intercept1<<std::endl;
  std::cout<<"    in +y side: k = "<<tmp_slope2<<", b = "<<tmp_intercept2<<std::endl;
  std::cout<<"  Calculated +y direction shift: "<<shift_phi_pos<<std::endl;
  std::cout<<"  Calculated -y direction shift: "<<shift_phi_neg<<std::endl;
  double halfZ_out = outerR*cos(dphi_sec/2.)/tan(atan(innerR / halfZ) - module_tilt_rad);
  double tilt_rad = atan( (outerR*sin(dphi_sec/2.)+shift_phi_pos)/(outerR*cos(dphi_sec/2.)) );
std::cout<<"  halfZ out = "<<outerR*cos(dphi_sec/2.)<<" / tan("<<atan(innerR / halfZ) - module_tilt_rad<<") "<<std::endl;
  outerR = sqrt( pow(outerR*cos(dphi_sec/2.), 2) + pow((outerR*sin(dphi_sec/2.) + shift_phi_pos) , 2) );
  std::cout<<"  Outer radius and halfZ considering this: "<<outerR<<", "<<halfZ_out<<std::endl;

  // Create the geometry
  DetElement ECAL(det_name, x_det.id());
  Volume worldVol = description.pickMotherVolume(ECAL);

  TGeoTube* EcalBarrel = new TGeoTube(innerR, outerR, halfZ_out);
  Volume EcalBarrelVol("EcalBarrel", EcalBarrel, air);
  //EcalBarrelVol.setVisAttributes(description, "CaloVis");

  // ====== Carbon fiber frame as supporting, at inner and outer of barrel. ====== //
  PolyhedraRegular innerFrame(nSec_phi, -dphi_sec/2., innerR, innerR + Cframe_thick, 2*halfZ);
  Volume innerFrame_vol("InnerCarbonFrame", innerFrame, MatCarbonfiber);

  PolyhedraRegular outerFrame(nSec_phi, tilt_rad, 
                              outerR*cos(dphi_sec/2.)-Cframe_thick, 
                              outerR*cos(dphi_sec/2.), 
                              2*halfZ_out);
  Volume outerFrame_vol("outerCarbonFrame", outerFrame, MatCarbonfiber);

  innerFrame_vol.setVisAttributes(description, "CarbonFrameVis");
  outerFrame_vol.setVisAttributes(description, "CarbonFrameVis");
  EcalBarrelVol.placeVolume(innerFrame_vol);
  EcalBarrelVol.placeVolume(outerFrame_vol);

  
  // ======== Define a sector (backspace is included in sector)======== //
  double innerR_sector = innerR + Cframe_thick;
  double outerR_sector = (innerR_sector + crystal_thick + back_space)/cos(dphi_sec/2.);;
  halfZ_out = outerR_sector*cos(dphi_sec/2.)/tan(atan(innerR_sector / halfZ) - module_tilt_rad);
std::cout<<"  halfZ out = "<<outerR_sector*cos(dphi_sec/2.)<<" / tan("<<atan(innerR_sector / halfZ) - module_tilt_rad<<") = "<<halfZ_out<<std::endl;

  // Sector inner and outer width
  double inner_width_neg = fabs(tmp_slope1*innerR_sector + tmp_intercept1);
  double inner_width_pos = fabs(tmp_slope2*innerR_sector + tmp_intercept2);
  double outer_width_neg = fabs(tmp_slope1*outerR_sector*cos(dphi_sec/2.) + tmp_intercept1);
  double outer_width_pos = fabs(tmp_slope2*outerR_sector*cos(dphi_sec/2.) + tmp_intercept2);

  std::cout<<std::endl;
  std::cout<<"  Calculate the new boundary lines in Sector: "<<std::endl;
  std::cout<<"    in -y side: k = "<<tmp_slope1<<", b = "<<tmp_intercept1<<std::endl;
  std::cout<<"    in +y side: k = "<<tmp_slope2<<", b = "<<tmp_intercept2<<std::endl;
  std::cout<<"  Calculated +y direction shift: "<<shift_phi_pos<<std::endl;
  std::cout<<"  Calculated -y direction shift: "<<shift_phi_neg<<std::endl;


  double vertices_sector[16] = {
    halfZ,  inner_width_pos,
    halfZ,  -inner_width_neg,
    -halfZ, -inner_width_neg,
    -halfZ, inner_width_pos,

    halfZ_out,  outer_width_pos,
    halfZ_out,  -outer_width_neg,
    -halfZ_out, -outer_width_neg,
    -halfZ_out, outer_width_pos,
  };

  EightPointSolid Sector("Sector", (crystal_thick+back_space)/2., vertices_sector );
  Volume Sector_vol("sector_vol", Sector, air);
  Sector_vol.setVisAttributes(description, "CaloVis");


  /// ------ Define module along Z axis ------  ///
  double theta_min = atan(innerR_sector/halfZ);
  double theta_max = M_PI - theta_min; 
  double theta_module = (theta_max-theta_min)/nModule_z; 
  for(int iz=0; iz<nModule_z; iz++){

    double theta_min_module = theta_min + iz*theta_module; 
    double theta_max_module = theta_min_module + theta_module; 
    double tilt_min = theta_min_module<M_PI/2. ? -module_tilt*M_PI/180. : module_tilt*M_PI/180.; 
    double tilt_max = theta_max_module<M_PI/2. ? -module_tilt*M_PI/180. : module_tilt*M_PI/180.; 
    double tilt_min_local = atan( (outerR_sector*cos(dphi_sec/2.)-innerR_sector)/ 
                                    (outerR_sector*cos(dphi_sec/2.)/tan(theta_min_module+tilt_min) - 
                                    innerR_sector/tan(theta_min_module)) );
    double tilt_max_local = atan( (outerR_sector*cos(dphi_sec/2.)-innerR_sector)/
                                    (outerR_sector*cos(dphi_sec/2.)/tan(theta_max_module+tilt_max) -
                                    innerR_sector/tan(theta_max_module)) );

    std::cout<<"Module #"<<iz<<" Geometry parameters: "<<std::endl;
    std::cout<<"Sector outerR - innerR = " << (outerR_sector*cos(dphi_sec/2.)-innerR_sector);
    std::cout<<", width (min) = "<<fabs(outerR_sector*cos(dphi_sec/2.)/tan(theta_min_module+tilt_min) -
                                    innerR_sector/tan(theta_min_module));
    std::cout<<", width (max) = "<<fabs(outerR_sector*cos(dphi_sec/2.)/tan(theta_max_module+tilt_max) -
                                    innerR_sector/tan(theta_max_module))<<std::endl;

    // Module height: keeping the total length = crystal thickness, so h = thick*sin(theta)
    // Exception in the most central module: h = thick. 
    double height_module = (theta_max_module<M_PI/2.) ? crystal_thick*fabs(sin(tilt_max_local)) : crystal_thick*fabs(sin(tilt_min_local));
    if( (theta_max_module-M_PI/2.)*(theta_min_module-M_PI/2.)<0 ) height_module = crystal_thick;

    //Module width: 
    double outer_width_module_neg = fabs(tmp_slope1*(innerR_sector+height_module) + tmp_intercept1);
    double outer_width_module_pos = fabs(tmp_slope2*(innerR_sector+height_module) + tmp_intercept2);


    std::cout<<"  Theta range: "<<theta_min_module<<", "<<theta_max_module<<", tilt angle "<<tilt_min<<", "<<tilt_max<<", local tilt angle "<<tilt_min_local<<", "<<tilt_max_local<<std::endl;
    std::cout<<"  height: "<<height_module<<", inner width "<<inner_width_neg<<", "<<inner_width_pos;
    std::cout<<", outer width "<<outer_width_module_neg<<", "<<outer_width_module_pos<<std::endl;
    
    double vertices_frame[16] = {
      innerR_sector/tan(theta_min_module), inner_width_pos,
      innerR_sector/tan(theta_min_module), -inner_width_neg,
      innerR_sector/tan(theta_max_module), -inner_width_neg,
      innerR_sector/tan(theta_max_module), inner_width_pos,

      innerR_sector/tan(theta_min_module)+height_module/tan(tilt_min_local), outer_width_module_pos,
      innerR_sector/tan(theta_min_module)+height_module/tan(tilt_min_local), -outer_width_module_neg,
      innerR_sector/tan(theta_max_module)+height_module/tan(tilt_max_local), -outer_width_module_neg,
      innerR_sector/tan(theta_max_module)+height_module/tan(tilt_max_local), outer_width_module_pos,
    };

    std::cout<<"  Inner surface vertices: "<<std::endl;
    for(int i=0; i<8; i+=2) std::cout<<"    ["<<vertices_frame[i]<<", "<<vertices_frame[i+1]<<"], "<<std::endl;
    std::cout<<"  Outer surface vertices: "<<std::endl;
    for(int i=8; i<16; i+=2) std::cout<<"    ["<<vertices_frame[i]<<", "<<vertices_frame[i+1]<<"], "<<std::endl;


    double vertices_crystal[16] = {
      innerR_sector/tan(theta_min_module)-Cseg_thick, inner_width_pos-Cseg_thick,
      innerR_sector/tan(theta_min_module)-Cseg_thick, -inner_width_neg+Cseg_thick,
      innerR_sector/tan(theta_max_module)+Cseg_thick, -inner_width_neg+Cseg_thick,
      innerR_sector/tan(theta_max_module)+Cseg_thick, inner_width_pos-Cseg_thick,

      innerR_sector/tan(theta_min_module)+height_module/tan(tilt_min_local)-Cseg_thick, outer_width_module_pos-Cseg_thick,
      innerR_sector/tan(theta_min_module)+height_module/tan(tilt_min_local)-Cseg_thick, -outer_width_module_neg+Cseg_thick,
      innerR_sector/tan(theta_max_module)+height_module/tan(tilt_max_local)+Cseg_thick, -outer_width_module_neg+Cseg_thick,
      innerR_sector/tan(theta_max_module)+height_module/tan(tilt_max_local)+Cseg_thick, outer_width_module_pos-Cseg_thick,
    };

    EightPointSolid ModuleFrame(Form("Module_%d", iz), height_module/2., vertices_frame );
    Volume ModuleFrame_vol(Form("module_vol%d", iz), ModuleFrame, MatCarbonfiber);
    ModuleFrame_vol.setVisAttributes(description, "CarbonFrameVis");

    EightPointSolid ModuleCrystal(Form("ModuleCrystal_%d", iz), height_module/2., vertices_crystal );
    Volume ModuleCrystal_vol(Form("moduleCrystal_vol%d", iz), ModuleCrystal, MatCrystal);
    ModuleCrystal_vol.setVisAttributes(description, "SensitiveVis");
    ModuleCrystal_vol.setSensitiveDetector(sens);

    PlacedVolume pv = ModuleFrame_vol.placeVolume(ModuleCrystal_vol);
    pv.addPhysVolID("stave", iz);

    Transform3D tr(Translation3D(0, 0, -(crystal_thick+back_space-height_module)/2.));
    Sector_vol.placeVolume(ModuleFrame_vol, tr);
  }



  //========  Loop to place sectors  ========  //
  for(int isec_phi=0; isec_phi<nSec_phi; isec_phi++){
    double phi = isec_phi * dphi_sec;
    Transform3D tr(RotationY(M_PI/2.0), Translation3D(innerR_sector+(crystal_thick+back_space)/2., 0., 0.));
    PlacedVolume pv = EcalBarrelVol.placeVolume(Sector_vol, RotationZ(phi) * tr);
    pv.addPhysVolID("sector", isec_phi);
  }
  
  // Finalize geometry
  PlacedVolume EcalBarrel_plv = worldVol.placeVolume(EcalBarrelVol);
  EcalBarrel_plv.addPhysVolID("system", x_det.id());
  ECAL.setPlacement(EcalBarrel_plv);

  std::cout << "--> Grainita_ECAL_Barrel_v02::create_detector() end" << std::endl;
  return ECAL;
}

DECLARE_DETELEMENT(Grainita_ECAL_Barrel_v02, create_detector)

//**************************************************************************

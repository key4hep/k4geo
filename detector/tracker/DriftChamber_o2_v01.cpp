//----------------------------------
//         DCH detector o2
//----------------------------------

/*!
 *  \brief     Detector constructor of DCH o2
 *  \details   This code creates full geometry of IDEA DCH subdetector
 *  \author    Alvaro Tolosa-Delgado alvaro.tolosa.delgado@cern.ch
 *  \author    Brieuc Francois       brieuc.francois@cern.ch
 *  \version   2
 *  \date      2024
 *  \pre       DD4hep compiled with Geant4+Qt
 */


#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Shapes.h"
#include "DD4hep/Printout.h"
#include "TMath.h"

#include <map>

namespace DCH_o2 {

// use alias of types to show more clear what the variable is
// if everything is double, the code is not readable
/// type for layer number
using DCH_layer = int;
/// tpye for lengths
using MyLength_t = double;
/// tpye for angles
using MyAngle_t = double;


// hardcode values in the code at compilation time
constexpr double PI = TMath::Pi();
constexpr double TWOPI = TMath::TwoPi();

// ancillary class to store DCH excel table
// to be filled before building the geometry!
// each DCH_info corresponds to one entry in the table
// global variables are static members of the class
class DCH_info
{
public:
    // global variables of DCH
    // inline before static allow initialization at this point
    inline static MyLength_t dch_Lhalf = {0};
    inline static MyLength_t dch_rin_z0 = {0};
    inline static MyLength_t dch_rin_z0_guard = {0};
    inline static MyLength_t dch_rout_z0 = {0};
    inline static MyLength_t dch_rout_z0_guard = {0};
    inline static int dch_ncell0 = {0};
    inline static int dch_ncell_increment = {0};
    inline static int dch_ncell_per_sector = {0};
    inline static int dch_nlayersPerSuperlayer = {0};
    inline static int dch_nsuperlayers = {0};
    inline static int dch_nlayers = {0};
    inline static MyAngle_t  dch_twist_angle = {0};

    // each member corresopnds to a column in the excel table
    /// layer number
    DCH_layer layer = {0};
    /// 2x number of wires in that layer
    int nwires = {0};
    /// cell parameter
    double height_z0 = {0.};
    /// cell parameter
    double width_z0 = {0.};

    /// radius (cylindral coord) of sensitive wire
    MyLength_t radius_sw_z0 = {0.};

    /// radius (cylindral coord) of 'down' field wires
    MyLength_t radius_fdw_z0 = {0.};

    /// radius (cylindral coord) of 'up' field wires
    MyLength_t radius_fuw_z0 = {0.};

    MyLength_t Radius_zLhalf(MyLength_t r_z0) const {
        return r_z0/cos(dch_twist_angle/2/dd4hep::rad);
    }


    // some quantityies are derived from previous ones
    // see excel and paper for equations
    ///  stereo angle is positive for odd layer number
    bool IsStereoPositive() const {
        return (1 == layer%2);
    }
    /// calculate sign based on IsStereoPositive
    inline int StereoSign() const {
        return (IsStereoPositive()*2 -1 );
    }

    /// separation between wires
    MyLength_t Pitch_z0(MyLength_t r_z0) const {
        return TWOPI*r_z0/nwires;
    };

    /// tan(stereoangle) = R(z=0)   / (L/2) * tan( twist_angle/2)
    inline MyAngle_t stereoangle_z0(MyLength_t r_z0) const {
        return atan( r_z0/dch_Lhalf*tan(dch_twist_angle/2));
    }

    /// tan(stereoangle) = R(z=L/2) / (L/2) * sin( twist_angle/2)
    inline MyAngle_t stereoangle_zLhalf(MyLength_t r_zLhalf) const {
        return atan( r_zLhalf/dch_Lhalf*sin(dch_twist_angle/2));
    }

    /// WireLength = 2*dch_Lhalf/cos(atan(Pitch_z0(r_z0)/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0))
    inline MyLength_t WireLength(MyLength_t r_z0) const {
        return  2*dch_Lhalf/cos(atan(Pitch_z0(r_z0)/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0)) ;
    };

    /// map to store excel table
    inline static std::map<DCH_layer, DCH_info> database;
    static bool IsDatabaseEmpty() {
        return (0 == database.size() );
    }

    static void Fill_DCH_info_database(dd4hep::Detector &desc);
    static void Show_DCH_info_database();

};

/// Function to build ARC endcaps
static dd4hep::Ref_t create_DCH_o2_v01(dd4hep::Detector &desc, dd4hep::xml::Handle_t handle, dd4hep::SensitiveDetector sens)
{
    dd4hep::xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    dd4hep::DetElement det(detName, detID);
    sens.setType("tracker");

    DCH_info::Fill_DCH_info_database(desc);
    if( DCH_info::IsDatabaseEmpty() )
        throw std::runtime_error("Empty database");
    DCH_info::Show_DCH_info_database();

    auto gasvol_material = desc.material("GasHe_90Isob_10");

    MyLength_t dch_SWire_thickness = desc.constantAsDouble("dch_SWire_thickness");
    MyLength_t dch_FSideWire_thickness = desc.constantAsDouble("dch_FSideWire_thickness");
    MyLength_t dch_FCentralWire_thickness = desc.constantAsDouble("dch_FCentralWire_thickness");

    /* Geometry tree:
     * Wall (tube) -> Gas (tube) -> Layer_1 (hyp) -> cell_1 (twisted tube)
     *                                            -> cell_... (twisted tube)
     *
     *                          -> Layer_... (hyp) -> cell_1 (twisted tube)
     *                                             -> cell_... (twisted tube)
     *
     * Layers represent a segmentation in radius
     * Sectors represent a segmentation in phi
     */

    MyLength_t safety_r_interspace   = 1    * dd4hep::nm;
    MyLength_t safety_z_interspace   = 1    * dd4hep::nm;
    MyLength_t safety_phi_interspace = 1e-6 * dd4hep::rad;


    MyLength_t vessel_thickness = desc.constantAsDouble("dch_vessel_thickness");
    dd4hep::Tube vessel_s(  DCH_info::dch_rin_z0 - vessel_thickness,
                            DCH_info::dch_rout_z0+ vessel_thickness,
                            DCH_info::dch_Lhalf  + vessel_thickness);
    dd4hep::Volume vessel_v (detName+"_vessel", vessel_s,  desc.material("Silicon") );
    vessel_v.setVisAttributes( desc.visAttributes("dch_no_vis") );


    dd4hep::Tube gas_s( DCH_info::dch_rin_z0,
                        DCH_info::dch_rout_z0,
                        DCH_info::dch_Lhalf + 2*safety_z_interspace );
    dd4hep::Volume gas_v (detName+"_gas", gas_s, gasvol_material );
    gas_v.setVisAttributes( desc.visAttributes("dch_no_vis") );
    vessel_v.placeVolume(gas_v);
//     dd4hep::Assembly vessel_v("vessel_v");

    for(const auto& [ilayer, l]  : DCH_info::database )
    {

//             if(ilayer>1)break;

//         MyLength_t rmin = l.Radius_zLhalf(l.radius_fdw_z0+safety_r_interspace);
//         MyLength_t rmax = l.Radius_zLhalf(l.radius_fuw_z0-safety_r_interspace);
//         MyLength_t dz = DCH_info::dch_Lhalf;
//         dd4hep::TwistedTube layer_s( DCH_info::dch_twist_angle,
//                                      rmin,
//                                      rmax,
//                                      dz,
//                                      1,
//                                      PI*dd4hep::rad
//                                    );
//         dd4hep::Hyperboloid layer_s(0.,0., 1*dd4hep::m, 30*dd4hep::deg, 1*dd4hep::m);

        // // // // // // // // // // // // // // // // // // // // /
        // // // // // INITIALIZATION OF THE LAYER // // // // // //
        // // // // // // // // // // // // // // // // // // // //
        // Hyperboloid parameters:
        /// inner radius at z=0
        MyLength_t rin   = l.radius_fdw_z0+safety_r_interspace;
        /// inner stereoangle, calculated from rin(z=0)
        MyAngle_t  stin  = l.stereoangle_z0(rin);
        /// outer radius at z=0
        MyLength_t rout  = l.radius_fuw_z0-safety_r_interspace;
        /// outer stereoangle, calculated from rout(z=0)
        MyAngle_t  stout = l.stereoangle_z0(rout);
        /// half-length
        MyLength_t dz    = DCH_info::dch_Lhalf + safety_z_interspace;

        dd4hep::Hyperboloid layer_s(rin, stin, rout, stout, dz);


        std::string layer_name = detName+"_layer"+std::to_string(ilayer);
        dd4hep::Volume layer_v ( layer_name , layer_s, gasvol_material );
        //layer_v.setVisAttributes( desc.visAttributes( Form("dch_layer_vis%d", ilayer%22) ) );
        layer_v.setVisAttributes( desc.visAttributes( "dch_no_vis") );
        auto layer_pv = gas_v.placeVolume(layer_v);
        layer_pv.addPhysVolID("layer", ilayer);

        dd4hep::DetElement layer_DE(det,layer_name+"DE", ilayer);
        layer_DE.setPlacement(layer_pv);
        // Assign the system ID to our mother volume


        // // // // // // // // // // // // // // // // // // // //
        // // // // // SEGMENTATION OF THE LAYER  // // // // // //
        // // // // // INTO CELLS (TWISTED TUBES) // // // // // //
        // // // // // // // // // // // // // // // // // // // //

        // ncells in this layer = 2x number of wires
        int ncells = l.nwires/2;
        MyAngle_t phi_step = (TWOPI/ncells)*dd4hep::rad;

        // unitary cell (Twisted tube) is repeated for each layer l.nwires/2 times
        // Twisted tube parameters
        MyAngle_t cell_twistangle    = l.StereoSign() * DCH_info::dch_twist_angle;
        MyLength_t cell_rin_z0       = l.radius_fdw_z0 + 2*safety_r_interspace;
        MyLength_t cell_rout_z0      = l.radius_fuw_z0 - 2*safety_r_interspace;
        MyLength_t cell_rin_zLhalf   = l.Radius_zLhalf(cell_rin_z0);
        MyLength_t cell_rout_zLhalf  = l.Radius_zLhalf(cell_rout_z0);
        MyLength_t cell_dz           = DCH_info::dch_Lhalf;
        MyAngle_t cell_phi_width     = phi_step - safety_phi_interspace;
        dd4hep::TwistedTube cell_s( cell_twistangle, cell_rin_zLhalf, cell_rout_zLhalf, cell_dz, 1, cell_phi_width);

        // initialize cell volume
        std::string cell_name = detName+"_layer"+std::to_string(ilayer)+"_cell";
        dd4hep::Volume cell_v (cell_name, cell_s, gasvol_material );
        cell_v.setSensitiveDetector(sens);
        cell_v.setVisAttributes( desc.visAttributes( "dch_gas_vis" /* Form("dch_layer_vis%d", ilayer%22)*/ ) );

        // // // // // // // // // // // // // // // // // // // //
        // // // // // // POSITIONING OF WIRES // // // // // // //
        // // // // // // // // // // // // // // // // // // // //
        {
            // // // // // // // // // // // // // // // // // // // //
            // // // // // // POSITIONING OF SENSE WIRES // // // // //
            // // // // // // // // // // // // // // // // // // // //
            // average radius to position sense wire
            MyLength_t cell_rave_z0 = 0.5*(cell_rin_z0+cell_rout_z0);
            MyLength_t cell_swire_radius = dch_SWire_thickness/2;
            MyLength_t swlength = 0.5*l.WireLength(cell_rave_z0) - cell_swire_radius*cos(l.stereoangle_z0(cell_rave_z0)) - safety_z_interspace;
            dd4hep::Tube swire_s(0., dch_SWire_thickness, swlength);
            dd4hep::Volume swire_v(cell_name+"_swire", swire_s, desc.material("W"));
            // Change sign of stereo angle to place properly the wire inside the twisted tube
            dd4hep::RotationX stereoTr( (-1.)*l.StereoSign()*l.stereoangle_z0(cell_rave_z0) );
            dd4hep::Transform3D swireTr ( stereoTr * dd4hep::Translation3D(cell_rave_z0,0.,0.) );
            cell_v.placeVolume(swire_v,swireTr);

            // // // // // // // // // // // // // // // // // // // //
            // // // // // // POSITIONING OF FIELD WIRES // // // // //
            // // // // // // // // // // // // // // // // // // // //
            //
            //  The following sketch represents the crossection of a DCH cell, where
            //      O symbol = Field wires, the number in parenthesis is used as ID
            //      X symbol = sense wire
            //
            //   ^ radius
            //
            //   O(1)---O(4)---O(6)    radius_z0 = l.radius_fuw_z0 == (++l).radius_fdw_z0
            //
            //   O(2)   X      O(7)    radius_z0 = average(l.radius_fuw_z0, l.radius_fdw_z0)
            //
            //   O(3)---O(5)---O(8)    radius_z0 = l.radius_fdw_z0 == (--l).radius_fuw_z0
            //
            //   --> phi axis
            //
            //  In the previous sketch, the wires are shared among several cells.
            //  Since we are using an actual shape to contain each cell,
            //  it is not feasible.
            //
            //  As a workaround, we introduce an offset in phi and radially to the cell center,
            //  in such a manner that the wires are fully contained in one cell.
            //  The following code implements the following sketch:
            //
            //   O(1)---O(4)---    radius_z0 = l.radius_fuw_z0 - wire_thickness/2
            //
            //   O(2)   X          radius_z0 = average(l.radius_fuw_z0, l.radius_fdw_z0)
            //
            //   O(3)---O(5)---    radius_z0 = l.radius_fdw_z0 + wire_thickness/2
            //
            //  phi_offset(n) = atan(  wire_thickness/2 / radius_z0 )
            //
            //  notice that the field wires are offcentered with respect to the sense wire
            //  by about 20um/1cm ~ 0.1 mrad, which is not expected to have any impact

            /// encapsulate the calculation of the phi offset into a function
            /// since it will be different for each field wire
            /// it includes the safety phi distance
            auto fwire_phi_offset = [&](MyLength_t radial_distance, MyLength_t wire_radius)->MyAngle_t
            {
                return atan(wire_radius/radial_distance)*dd4hep::rad + safety_phi_interspace;
            };

            // // // // // // // // // // // // // // // // // // // //
            // // // // // // POSITIONING OF F WIRE 2 // // // // // //
            // // // // // // REQUIRES OFFSET ON PHI  // // // // // //
            // // // // // // // // // // // // // // // // // // // //
            {
                MyLength_t fwire_radius = dch_FCentralWire_thickness/2;
                MyLength_t fwire_r_z0   = cell_rave_z0;
                MyAngle_t  fwire_stereo =  (-1.)*l.StereoSign()*l.stereoangle_z0(fwire_r_z0);
                MyAngle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                MyLength_t fwire_length = 0.5*l.WireLength(fwire_r_z0) - fwire_radius*cos(l.stereoangle_z0(fwire_r_z0)) - safety_z_interspace;;

                dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                dd4hep::Volume fwire_v(cell_name+"_f2wire", fwire_s, desc.material("W"));
                // Change sign of stereo angle to place properly the wire inside the twisted tube
                dd4hep::RotationX fwireStereoTr( fwire_stereo );
                dd4hep::RotationZ fwirePhoTr( fwire_phi );
                dd4hep::Transform3D fwireTr ( fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
                cell_v.placeVolume(fwire_v,fwireTr);
            }

            // // // // // // // // // // // // // // // // // // // //
            // // // // // // POSITIONING OF F WIRE 1    // // // // //
            // // // // // // REQUIRES OFFSET ON PHI & R // // // // //
            // // // // // // // // // // // // // // // // // // // //
            {
                MyLength_t fwire_radius = dch_FSideWire_thickness/2;
                // decrease radial distance, move it closer to the sense wire
                MyLength_t fwire_r_z0   = cell_rout_z0 - fwire_radius;
                MyAngle_t  fwire_stereo =  (-1.)*l.StereoSign()*l.stereoangle_z0(fwire_r_z0);
                MyAngle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                MyLength_t fwire_length = 0.5*l.WireLength(fwire_r_z0) - fwire_radius*cos(l.stereoangle_z0(fwire_r_z0)) - safety_z_interspace;;

                dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                dd4hep::Volume fwire_v(cell_name+"_f1wire", fwire_s, desc.material("W"));
                // Change sign of stereo angle to place properly the wire inside the twisted tube
                dd4hep::RotationX fwireStereoTr( fwire_stereo );
                dd4hep::RotationZ fwirePhoTr( fwire_phi );
                dd4hep::Transform3D fwireTr ( fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
                cell_v.placeVolume(fwire_v,fwireTr);
            }
            // // // // // // // // // // // // // // // // // // // //
            // // // // // // POSITIONING OF F WIRE 3    // // // // //
            // // // // // // REQUIRES OFFSET ON PHI & R // // // // //
            // // // // // // // // // // // // // // // // // // // //
            {
                MyLength_t fwire_radius = dch_FSideWire_thickness/2;
                // increase radial distance, move it closer to the sense wire
                MyLength_t fwire_r_z0   = cell_rin_z0 + fwire_radius;
                MyAngle_t  fwire_stereo =  (-1.)*l.StereoSign()*l.stereoangle_z0(fwire_r_z0);
                MyAngle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                MyLength_t fwire_length = 0.5*l.WireLength(fwire_r_z0) - fwire_radius*cos(l.stereoangle_z0(fwire_r_z0)) - safety_z_interspace;;

                dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                dd4hep::Volume fwire_v(cell_name+"_f3wire", fwire_s, desc.material("W"));
                // Change sign of stereo angle to place properly the wire inside the twisted tube
                dd4hep::RotationX fwireStereoTr( fwire_stereo );
                dd4hep::RotationZ fwirePhoTr( fwire_phi );
                dd4hep::Transform3D fwireTr ( fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
                cell_v.placeVolume(fwire_v,fwireTr);
            }

        }/// end building wires
        for(int nphi = 0; nphi < 3 /*ncells*/; ++nphi)
        {
            // TODO: check if staggering is just + 0.5*cell_phi_width*(ilayer%2);
            // phi positioning, adding offset for odd ilayers
            MyAngle_t cell_phi_angle = phi_step * nphi + 0.5*cell_phi_width*(ilayer%2);
            // conversion of RotationZ into Transform3D using constructor)
            dd4hep::Transform3D cellTr { dd4hep::RotationZ(cell_phi_angle) };
            auto cell_pv = layer_v.placeVolume(cell_v, cellTr);
            cell_pv.addPhysVolID("nphi", nphi);
            cell_pv.addPhysVolID("stereosign", l.StereoSign() );

            dd4hep::DetElement cell_DE(layer_DE,cell_name+std::to_string(nphi)+"DE", nphi);
            cell_DE.setPlacement(cell_pv);
        }


    }




    // Place our mother volume in the world
    dd4hep::Volume wVol = desc.pickMotherVolume(det);
    dd4hep::PlacedVolume vessel_pv = wVol.placeVolume(vessel_v);
    // Associate the wall to the detector element.
    det.setPlacement(vessel_pv);
    // Assign the system ID to our mother volume
    vessel_pv.addPhysVolID("system", detID);



    // Box wall_sh( DCH_info::dch_rin_z0  - wall_thickness ,
    //               DCH_info::dch_rout_z0 + wall_thickness,
    //               DCH_info::dch_Lhalf   + wall_thickness
    //              );
    // Volume wall_v( "DCH_wall", wall_sh, desc.material("Air"));
    // wall_v.setSensitiveDetector(sens);
    // wall_v.setVisAttributes(desc.visAttributes("cooling_vis"));
    /*
     *        Tube gas_sh( DCH_info::dch_rin_z0  ,
     *                     DCH_info::dch_rout_z0 ,
     *                     DCH_info::dch_Lhalf
     *                     );
     *        Volume gas_v( "DCH_wall", gas_sh, desc.material("Air"));*/

    // // Place our mother volume in the world
    // Volume wVol = desc.pickMotherVolume(det);
    // PlacedVolume wall_pv = wVol.placeVolume(wall_sh);
    //
    // // Assign the system ID to our mother volume
    // wall_pv.addPhysVolID("system", detID);
    //
    // // Associate the silicon Placed Volume to the detector element.
    // det.setPlacement(wall_sh);

    // // Create the mother Detector element to be returned at the end
    // double twist_angle = 45*dd4hep::degree;
    // double rmin = 1*dd4hep::cm;
    // double rmax = 5*dd4hep::cm;
    // double dz = 5*dd4hep::cm;
    // double dphi = 90*dd4hep::deg;
    // int nsegments = 1;
    // TwistedTube myshape( twist_angle,  rmin,  rmax, -dz, dz, dphi);
    // // Define volume (shape+material)
    // Volume siVol(detName +"_sensor", myshape, desc.material("Silicon"));
    // siVol.setVisAttributes(desc.visAttributes("sensor_vis"));
    // siVol.setSensitiveDetector(sens);
    //
    // double wireR = 0.5*rmin + 0.5*rmax; //0.5*(rmax+rmin);
    // double wirePhi = 0.3*dphi;
    // double stereoangle = atan( wireR/dz*sin(twist_angle/2/rad) );
    //
    // double wireThickness = 0.1*rmin;
    // Tube ancshape( 0., wireThickness, dz/cos(stereoangle/rad) - wireThickness*tan(stereoangle/rad) );
    // Volume ancVol(detName +"ancshape", ancshape, desc.material("Silicon"));
    // ancVol.setVisAttributes(desc.visAttributes("cooling_vis"));
    //
    // auto ancshapeTrp = RotationZ(0.3*dphi)*Translation3D( wireR,0,0) * RotationX(-stereoangle/rad);
    // auto ancshapeTr0 = RotationZ(  0     )*Translation3D( wireR,0,0) * RotationX(-stereoangle/rad);
    // auto ancshapeTrm = RotationZ(-0.3*dphi)*Translation3D( wireR,0,0) * RotationX(-stereoangle/rad);
    //
    // // Place our mother volume in the world
    // Volume wVol = desc.pickMotherVolume(det);
    //
    // PlacedVolume siPV = wVol.placeVolume(siVol);
    // wVol.placeVolume(ancVol, ancshapeTrp);
    // wVol.placeVolume(ancVol, ancshapeTr0);
    // wVol.placeVolume(ancVol, ancshapeTrm);
    //
    // // Assign the system ID to our mother volume
    // siPV.addPhysVolID("system", detID);
    //
    // // Associate the silicon Placed Volume to the detector element.
    // det.setPlacement(siPV);

    return det;

}

void DCH_info::Fill_DCH_info_database(dd4hep::Detector & desc)
{
    // do not fill twice the database
    if( not IsDatabaseEmpty() ) return;

    // DCH outer geometry dimensions
    DCH_info::dch_rin_z0 = desc.constantAsDouble("dch_inner_cyl_R_z0");
    DCH_info::dch_rout_z0  = desc.constantAsDouble("dch_outer_cyl_R_z0");
    DCH_info::dch_Lhalf = desc.constantAsDouble("dch_Lhalf");

    // guard wires position, fix position
    DCH_info::dch_rin_z0_guard  = desc.constantAsDouble("dch_in_guard_z0");
    DCH_info::dch_rout_z0_guard = desc.constantAsDouble("dch_out_guard_zL2");

    /// number of cells of first layer
    DCH_info::dch_ncell0 = desc.constantAsLong("dch_ncell");
    /// increment of cell number by each superlayer
    DCH_info::dch_ncell_increment = desc.constantAsLong("dch_ncell_increment");
    /// layer cells are grouped into sectors
    DCH_info::dch_ncell_per_sector = desc.constantAsLong("dch_ncell_per_sector");
    // if dch_ncell_per_sector is not divisor of dch_ncell0 and dch_ncell_increment
    // trow an error
    if( 0 != (dch_ncell0 % dch_ncell_per_sector) || 0 != (dch_ncell_increment % dch_ncell_per_sector) )
        throw std::runtime_error("dch_ncell_per_sector is not divisor of dch_ncell0 or dch_ncell_increment");

    // wires/cells are grouped in layers (and superlayers)
    DCH_info::dch_nsuperlayers = desc.constantAsLong("dch_nsuperlayers");
    DCH_info::dch_nlayersPerSuperlayer = desc.constantAsLong("dch_nlayersPerSuperlayer");
    /// nlayers = nsuperlayers * nlayersPerSuperlayer
    /// default: 112 = 14 * 8
    DCH_info::dch_nlayers = DCH_info::dch_nsuperlayers * DCH_info::dch_nlayersPerSuperlayer;

    MyAngle_t dch_alpha = desc.constantAsDouble("dch_alpha");
    DCH_info::dch_twist_angle = 2*dch_alpha;

    // retrieve some initial values for the first layer
    double dch_first_width = desc.constantAsDouble("dch_first_width");
    MyLength_t dch_first_sense_r = desc.constantAsDouble("dch_first_sense_r");

    // intialize layer 1 from input parameters
    {
        DCH_info layer1_info;
        layer1_info.layer         = 1;
        layer1_info.nwires        = 2*DCH_info::dch_ncell0;
        layer1_info.height_z0     = dch_first_width;
        layer1_info.radius_sw_z0  = dch_first_sense_r;
        layer1_info.radius_fdw_z0 = dch_first_sense_r - 0.5*dch_first_width;
        layer1_info.radius_fuw_z0 = dch_first_sense_r + 0.5*dch_first_width;
        layer1_info.width_z0      = TWOPI*dch_first_sense_r/DCH_info::dch_ncell0;

        DCH_info::database.emplace(1, layer1_info);
    }

    // some parameters of the following layer are calculated based on the previous ones
    // the rest are left as methods of DCH_info class
    // loop over all layers, skipping the first one
    for(int ilayer = 2; ilayer<= DCH_info::dch_nlayers; ++ilayer)
    {
        // initialize here an object that will contain one row from the excel table
        DCH_info layer_info;

        // the loop counter actually corresponds to the layer number
        layer_info.layer = ilayer;
        // calculate number of superlayer
        // WARNING: division of integers on purpose!
        int nsuperlayer_minus_1 = (ilayer-1)/dch_nlayersPerSuperlayer;
        layer_info.nwires = 2*(DCH_info::dch_ncell0 + DCH_info::dch_ncell_increment*nsuperlayer_minus_1 );

        // the previous layer info is needed to calculate parameters of current layer
        auto previousLayer = DCH_info::database.at(ilayer-1);

        //calculate height_z0, radius_sw_z0
        {
            double h  = previousLayer.height_z0;
            double ru = previousLayer.radius_fuw_z0;
            double rd = previousLayer.radius_fdw_z0;

            if(0 == nsuperlayer_minus_1)
                layer_info.height_z0 = h*ru/rd;
            else
                layer_info.height_z0 = TWOPI*ru/(0.5*layer_info.nwires - PI);

            layer_info.radius_sw_z0 = 0.5*layer_info.height_z0 + ru;
        }

        //calculate radius_fdw_z0, radius_fuw_z0, width_z0
        layer_info.radius_fdw_z0 = previousLayer.radius_fuw_z0;
        layer_info.radius_fuw_z0 = previousLayer.radius_fuw_z0 + layer_info.height_z0;
        layer_info.width_z0 = TWOPI*layer_info.radius_sw_z0/(0.5*layer_info.nwires);

        // according to excel, width_z0 == height_z0
        if(fabs(layer_info.width_z0 - layer_info.height_z0)>1e-4)
            throw std::runtime_error("fabs(l.width_z0 - l.height_z0)>1e-4");



        DCH_info::database.emplace(ilayer, layer_info);
    }

    std::cout << "\t+ Total size of DCH database = " << DCH_info::database.size() << std::endl;
    return;
}

void DCH_info::Show_DCH_info_database()
{
    if( IsDatabaseEmpty() )
        throw std::runtime_error("Can not show empty database!");

    std::cout << "\n";
    std::cout << "Global parameters of DCH:\n";
    std::cout << "\tHalf length/mm = " << dch_Lhalf/dd4hep::mm << '\n';
    std::cout << "\tRadius in/mm  = " << dch_rin_z0/dd4hep::mm << '\n';
    std::cout << "\tRadius out/mm = " << dch_rout_z0/dd4hep::mm<< '\n';
    std::cout << "\tRadius guard in/mm  = " << dch_rin_z0_guard/dd4hep::mm << '\n';
    std::cout << "\tRadius guard out/mm = " << dch_rout_z0_guard/dd4hep::mm << '\n';
    std::cout << "\n";
    std::cout << "\tTwis angle (2*alpha) / deg = " << dch_twist_angle/dd4hep::deg << '\n';
    std::cout << "\n";
    std::cout << "\tN superlayers = " << dch_nsuperlayers << '\n';
    std::cout << "\tN layers per superlayer = " << dch_nlayersPerSuperlayer << '\n';
    std::cout << "\tN layers = " << dch_nlayers << '\n';
    std::cout << "\n";
    std::cout << "\tN cells layer1 = " << dch_ncell0 << '\n';
    std::cout << "\tN cells increment per superlayer = " << dch_ncell_increment << '\n';
    std::cout << "\tN cells per sector = " << dch_ncell_per_sector << '\n';
    std::cout << "\n";
    std::cout << "Layer parameters of DCH:\n";
    std::cout
            << "\t" << "layer"
            << "\t" << "nwires"
            << "\t" << "height_z0/mm"
            << "\t" << "width_z0/mm"
            << "\t" << "radius_fdw_z0/mm"
            << "\t" << "radius_sw_z0/mm"
            << "\t" << "radius_fuw_z0/mm"
            << "\t" << "stereoangle_z0/deg"
            << "\t" << "Pitch_z0/mm"
            << "\t" << "radius_sw_zLhalf/mm"
            << "\t" << "WireLength/mm"
            << "\n" << std::endl;

    for(const auto& [nlayer, l]  : DCH_info::database )
    {
        std::cout
                << "\t" << l.layer
                << "\t" << l.nwires
                << "\t" << l.height_z0/dd4hep::mm
                << "\t" << l.width_z0/dd4hep::mm
                << "\t" << l.radius_fdw_z0/dd4hep::mm
                << "\t" << l.radius_sw_z0/dd4hep::mm
                << "\t" << l.radius_fuw_z0/dd4hep::mm
                << "\t" << l.StereoSign()*l.stereoangle_z0(l.radius_sw_z0)/dd4hep::deg
                << "\t" << l.Pitch_z0(l.radius_sw_z0)/dd4hep::mm
                << "\t" << l.Radius_zLhalf(l.radius_sw_z0)/dd4hep::mm
                << "\t" << l.WireLength(l.radius_sw_z0)/dd4hep::mm
                << "\n" << std::endl;
    }
    return;
}



}; // end DCH_o2 namespace
DECLARE_DETELEMENT(DriftChamber_o2_v01_T, DCH_o2::create_DCH_o2_v01)

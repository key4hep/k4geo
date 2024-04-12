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

namespace kk{
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
/* Data structure to store DCH excel table
 * to be filled before building the geometry!
 *
 * Global variables are members of this class
 * Helper class DCH_info contains one entry in the table
 * which corresponds to one layer
 *
 * @author A. Tolosa Delgado, CERN
 * @date April 2024
 * @version Drift chamber v2
 *
 */
class DCH_info
{
public:
    // global variables of DCH
    MyLength_t dch_Lhalf = {0};
    MyLength_t dch_rin = {0};
    MyLength_t dch_rin_z0_guard = {0};
    MyLength_t dch_rout = {0};
    MyLength_t dch_rout_z0_guard = {0};
    int dch_ncell0 = {0};
    int dch_ncell_increment = {0};
    int dch_ncell_per_sector = {0};
    int dch_nlayersPerSuperlayer = {0};
    int dch_nsuperlayers = {0};
    int dch_nlayers = {0};
    MyAngle_t  dch_twist_angle = {0};

    /// Get number of cells in a given layer
    ///   ncells = 2x number of wires
    inline int Get_ncells(int ilayer){return database.at(ilayer).nwires/2;}

    /// Get phi width for the twisted tube and the step (phi distance between cells)
    inline MyAngle_t Get_phi_width(int ilayer){return (TWOPI/Get_ncells(ilayer))*dd4hep::rad;}

    ///
    // TODO: check if staggering is just + 0.25*cell_phi_width*(ilayer%2);
    // phi positioning, adding offset for odd ilayers
    MyAngle_t Get_cell_phi_angle(int ilayer, int nphi){ return (Get_phi_width(ilayer) * (nphi + 0.25*(ilayer%2)));}

    /// calculate superlayer for a given ilayer.
    /// WARNING: division of integers on purpose!
    int Get_nsuperlayer_minus_1(int ilayer){ return int((ilayer-1)/dch_nlayersPerSuperlayer);}

    /// Calculate radius at z=L/2 given at z=0
    MyLength_t Radius_zLhalf(MyLength_t r_z0) const {
        return r_z0/cos(dch_twist_angle/2/dd4hep::rad);
    }

    /// tan(stereoangle) = R(z=0)   / (L/2) * tan( twist_angle/2)
    inline MyAngle_t stereoangle_z0(MyLength_t r_z0) const {
        return atan( r_z0/dch_Lhalf*tan(dch_twist_angle/2));
    }

    /// tan(stereoangle) = R(z=L/2) / (L/2) * sin( twist_angle/2)
    inline MyAngle_t stereoangle_zLhalf(MyLength_t r_zLhalf) const {
        return atan( r_zLhalf/dch_Lhalf*sin(dch_twist_angle/2));
    }

    /// WireLength = 2*dch_Lhalf/cos(atan(Pitch_z0(r_z0)/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0))
    inline MyLength_t WireLength(int nlayer, MyLength_t r_z0) const {
        auto Pitch_z0 = database.at(nlayer).Pitch_z0(r_z0);
        return  2*dch_Lhalf/cos(atan(Pitch_z0/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0)) ;
    };

    /// Internal helper struct for defining the layer layout
    struct DCH_info_layer
    {
        // each member corresponds to a column in the excel table
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

        /// separation between wires (along the circle)
        MyLength_t Pitch_z0(MyLength_t r_z0) const {
            return TWOPI*r_z0/nwires;
        };

    };
    /// map to store excel table
    std::map<DCH_layer, DCH_info_layer> database;
    bool IsDatabaseEmpty() { return (0 == database.size() ); }

    void Fill_DCH_info_database(dd4hep::Detector &desc);
    void Show_DCH_info_database();

};
} // end namespace kk

namespace DCH_v2 {

using namespace kk;

/// Function to build DCH
static dd4hep::Ref_t create_DCH_o1_v02(dd4hep::Detector &desc, dd4hep::xml::Handle_t handle, dd4hep::SensitiveDetector sens)
{
    dd4hep::xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    dd4hep::DetElement det(detName, detID);
    sens.setType("tracker");

    kk::DCH_info DCH_i;
    DCH_i.Fill_DCH_info_database(desc);
    if( DCH_i.IsDatabaseEmpty() )
        throw std::runtime_error("Empty database");

    bool printExcelTable = detElem.attr<bool>(_Unicode(printExcelTable));
    if(printExcelTable)
        DCH_i.Show_DCH_info_database();

    auto gasElem    = detElem.child("gas");
    auto gasvolMat  = desc.material(gasElem.attr<std::string>(_Unicode(material)));
    auto gasvolVis  = desc.visAttributes(gasElem.attr<std::string>(_Unicode(vis)));

    auto vesselElem = detElem.child("vessel");
    auto vesselSkinMat  = desc.material(vesselElem.attr<std::string>(_Unicode(material)));
    auto vesselSkinVis  = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(vis)));


    auto wiresElem = detElem.child("wires");
    auto wiresMat  = desc.material(wiresElem.attr<std::string>(_Unicode(material)));
    auto wiresVis  = desc.visAttributes(wiresElem.attr<std::string>(_Unicode(vis)));
    bool buildSenseWires  = wiresElem.attr<bool>(_Unicode(buildSenseWires));
    bool buildFieldWires  = wiresElem.attr<bool>(_Unicode(buildFieldWires));

    MyLength_t dch_SWire_thickness        = wiresElem.attr<double>(_Unicode(SWire_thickness)) ;
    MyLength_t dch_FSideWire_thickness    = wiresElem.attr<double>(_Unicode(FSideWire_thickness)) ;
    MyLength_t dch_FCentralWire_thickness = wiresElem.attr<double>(_Unicode(FCentralWire_thickness)) ;

    /* Geometry tree:
     * Wall (tube) -> Gas (tube) -> Layer_1 (hyp) -> cell_1 (twisted tube)
     *                                            -> cell_... (twisted tube)
     *
     *                          -> Layer_... (hyp) -> cell_1 (twisted tube)
     *                                             -> cell_... (twisted tube)
     *
     * Layers represent a segmentation in radius
     * Sectors represent a segmentation in phi
     * Each cell corresponds to a Detector Element
     */

    MyLength_t safety_r_interspace   = 1    * dd4hep::nm;
    MyLength_t safety_z_interspace   = 1    * dd4hep::nm;
    MyLength_t safety_phi_interspace = 1e-6 * dd4hep::rad;


    MyLength_t vessel_thickness = desc.constantAsDouble("DCH_vessel_thickness");
    dd4hep::Tube vessel_s(  DCH_i.dch_rin - vessel_thickness,
                            DCH_i.dch_rout+ vessel_thickness,
                            DCH_i.dch_Lhalf  + vessel_thickness);
    dd4hep::Volume vessel_v (detName+"_vessel", vessel_s,  vesselSkinMat );
    vessel_v.setVisAttributes( vesselSkinVis );


    dd4hep::Tube gas_s( DCH_i.dch_rin,
                        DCH_i.dch_rout,
                        DCH_i.dch_Lhalf + 2*safety_z_interspace );
    dd4hep::Volume gas_v (detName+"_gas", gas_s, gasvolMat );
    gas_v.setVisAttributes( gasvolVis );
    vessel_v.placeVolume(gas_v);

    for(const auto& [ilayer, l]  : DCH_i.database )
    {

//             if(ilayer>1)break;

        // // // // // // // // // // // // // // // // // // // // /
        // // // // // INITIALIZATION OF THE LAYER // // // // // //
        // // // // // // // // // // // // // // // // // // // //
        // Hyperboloid parameters:
        /// inner radius at z=0
        MyLength_t rin   = l.radius_fdw_z0+safety_r_interspace;
        /// inner stereoangle, calculated from rin(z=0)
        MyAngle_t  stin  = DCH_i.stereoangle_z0(rin);
        /// outer radius at z=0
        MyLength_t rout  = l.radius_fuw_z0-safety_r_interspace;
        /// outer stereoangle, calculated from rout(z=0)
        MyAngle_t  stout = DCH_i.stereoangle_z0(rout);
        /// half-length
        MyLength_t dz    = DCH_i.dch_Lhalf + safety_z_interspace;

        dd4hep::Hyperboloid layer_s(rin, stin, rout, stout, dz);


        std::string layer_name = detName+"_layer"+std::to_string(ilayer);
        dd4hep::Volume layer_v ( layer_name , layer_s, gasvolMat );
        layer_v.setVisAttributes( desc.visAttributes( Form("dch_layer_vis%d", ilayer%22) ) );
        auto layer_pv = gas_v.placeVolume(layer_v);
        layer_pv.addPhysVolID("layer", ilayer);
        // add superlayer bitfield
        int nsuperlayer_minus_1 = DCH_i.Get_nsuperlayer_minus_1(ilayer);
        layer_pv.addPhysVolID("superlayer", nsuperlayer_minus_1 );

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
        MyAngle_t cell_twistangle    = l.StereoSign() * DCH_i.dch_twist_angle;
        MyLength_t cell_rin_z0       = l.radius_fdw_z0 + 2*safety_r_interspace;
        MyLength_t cell_rout_z0      = l.radius_fuw_z0 - 2*safety_r_interspace;
        MyLength_t cell_rin_zLhalf   = DCH_i.Radius_zLhalf(cell_rin_z0);
        MyLength_t cell_rout_zLhalf  = DCH_i.Radius_zLhalf(cell_rout_z0);
        MyLength_t cell_dz           = DCH_i.dch_Lhalf;
        MyAngle_t cell_phi_width     = phi_step - safety_phi_interspace;
        dd4hep::TwistedTube cell_s( cell_twistangle, cell_rin_zLhalf, cell_rout_zLhalf, cell_dz, 1, cell_phi_width);

        // initialize cell volume
        std::string cell_name = detName+"_layer"+std::to_string(ilayer)+"_cell";
        dd4hep::Volume cell_v (cell_name, cell_s, gasvolMat );
        cell_v.setSensitiveDetector(sens);
        cell_v.setVisAttributes( desc.visAttributes( "dch_no_vis_nodaughters" ) );

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
            MyLength_t swlength = 0.5*DCH_i.WireLength(ilayer,cell_rave_z0)
                                - cell_swire_radius*cos(DCH_i.stereoangle_z0(cell_rave_z0))
                                - safety_z_interspace;
            if(buildSenseWires)
            {
                dd4hep::Tube swire_s(0., dch_SWire_thickness, swlength);
                dd4hep::Volume swire_v(cell_name+"_swire", swire_s, wiresMat);
                swire_v.setVisAttributes( wiresVis );
                // Change sign of stereo angle to place properly the wire inside the twisted tube
                dd4hep::RotationX stereoTr( (-1.)*l.StereoSign()*DCH_i.stereoangle_z0(cell_rave_z0) );
                dd4hep::Transform3D swireTr ( stereoTr * dd4hep::Translation3D(cell_rave_z0,0.,0.) );
                cell_v.placeVolume(swire_v,swireTr);
            }
            if(buildFieldWires)
            {
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
                // // // // // // REQUIRES OFFSET OF PHI  // // // // // //
                // // // // // // // // // // // // // // // // // // // //
                {
                    MyLength_t fwire_radius = dch_FCentralWire_thickness/2;
                    MyLength_t fwire_r_z0   = cell_rave_z0;
                    MyAngle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i.stereoangle_z0(fwire_r_z0);
                    MyAngle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                    MyLength_t fwire_length = 0.5*DCH_i.WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i.stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f2wire", fwire_s, wiresMat );
                    fwire_v.setVisAttributes( wiresVis );
                    // Change sign of stereo angle to place properly the wire inside the twisted tube
                    dd4hep::RotationX fwireStereoTr( fwire_stereo );
                    dd4hep::RotationZ fwirePhoTr( fwire_phi );
                    dd4hep::Transform3D fwireTr ( fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
                    cell_v.placeVolume(fwire_v,fwireTr);
                }

                // // // // // // // // // // // // // // // // // // // //
                // // // // // // POSITIONING OF F WIRE 1    // // // // //
                // // // // // // REQUIRES OFFSET OF PHI & R // // // // //
                // // // // // // // // // // // // // // // // // // // //
                {
                    MyLength_t fwire_radius = dch_FSideWire_thickness/2;
                    // decrease radial distance, move it closer to the sense wire
                    MyLength_t fwire_r_z0   = cell_rout_z0 - fwire_radius;
                    MyAngle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i.stereoangle_z0(fwire_r_z0);
                    MyAngle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                    MyLength_t fwire_length = 0.5*DCH_i.WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i.stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f1wire", fwire_s, wiresMat );
                    fwire_v.setVisAttributes( wiresVis );
                    // Change sign of stereo angle to place properly the wire inside the twisted tube
                    dd4hep::RotationX fwireStereoTr( fwire_stereo );
                    dd4hep::RotationZ fwirePhoTr( fwire_phi );
                    dd4hep::Transform3D fwireTr ( fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
                    cell_v.placeVolume(fwire_v,fwireTr);
                }
                // // // // // // // // // // // // // // // // // // // //
                // // // // // // POSITIONING OF F WIRE 3    // // // // //
                // // // // // // REQUIRES OFFSET OF PHI & R // // // // //
                // // // // // // // // // // // // // // // // // // // //
                {
                    MyLength_t fwire_radius = dch_FSideWire_thickness/2;
                    // increase radial distance, move it closer to the sense wire
                    MyLength_t fwire_r_z0   = cell_rin_z0 + fwire_radius;
                    MyAngle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i.stereoangle_z0(fwire_r_z0);
                    MyAngle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                    MyLength_t fwire_length = 0.5*DCH_i.WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i.stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f3wire", fwire_s, wiresMat );
                    fwire_v.setVisAttributes( wiresVis );
                    // Change sign of stereo angle to place properly the wire inside the twisted tube
                    dd4hep::RotationX fwireStereoTr( fwire_stereo );
                    dd4hep::RotationZ fwirePhoTr( fwire_phi );
                    dd4hep::Transform3D fwireTr ( fwirePhoTr * fwireStereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
                    cell_v.placeVolume(fwire_v,fwireTr);
                }
                // // // // // // // // // // // // // // // // // // // //
                // // // // // // POSITIONING OF F WIRE 5    // // // // //
                // // // // // // REQUIRES OFFSET OF R       // // // // //
                // // // // // // // // // // // // // // // // // // // //
                {
                    MyLength_t fwire_radius = dch_FSideWire_thickness/2;
                    // increase radial distance, move it closer to the sense wire
                    MyLength_t fwire_r_z0   = cell_rin_z0 + fwire_radius;
                    MyAngle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i.stereoangle_z0(fwire_r_z0);
                    MyLength_t fwire_length = 0.5*DCH_i.WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i.stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f5wire", fwire_s, wiresMat );
                    fwire_v.setVisAttributes( wiresVis );
                    // Change sign of stereo angle to place properly the wire inside the twisted tube
                    dd4hep::RotationX fwireStereoTr( fwire_stereo );
                    dd4hep::Transform3D fwireTr ( fwireStereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
                    cell_v.placeVolume(fwire_v,fwireTr);
                }
                // // // // // // // // // // // // // // // // // // // //
                // // // // // // POSITIONING OF F WIRE 4    // // // // //
                // // // // // // REQUIRES OFFSET OF R       // // // // //
                // // // // // // // // // // // // // // // // // // // //
                {
                    MyLength_t fwire_radius = dch_FSideWire_thickness/2;
                    // increase radial distance, move it closer to the sense wire
                    MyLength_t fwire_r_z0   = cell_rout_z0 - fwire_radius;
                    MyAngle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i.stereoangle_z0(fwire_r_z0);
                    MyLength_t fwire_length = 0.5*DCH_i.WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i.stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f4wire", fwire_s, wiresMat );
                    fwire_v.setVisAttributes( wiresVis );
                    // Change sign of stereo angle to place properly the wire inside the twisted tube
                    dd4hep::RotationX fwireStereoTr( fwire_stereo );
                    dd4hep::Transform3D fwireTr ( fwireStereoTr * dd4hep::Translation3D(fwire_r_z0,0.,0.) );
                    cell_v.placeVolume(fwire_v,fwireTr);
                }
            }// end building field wires
        }/// end building wires
        for(int nphi = 0; nphi < ncells; ++nphi)
        {
            // TODO: check if staggering is just + 0.25*cell_phi_width*(ilayer%2);
            // phi positioning, adding offset for odd ilayers
            MyAngle_t cell_phi_angle = phi_step * nphi + 0.25*cell_phi_width*(ilayer%2);
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

    return det;

}

}; // end DCH_o2 namespace

namespace kk{
void DCH_info::Fill_DCH_info_database(dd4hep::Detector & desc)
{
    // do not fill twice the database
    if( not this->IsDatabaseEmpty() ) return;

    // DCH outer geometry dimensions
    this->dch_rin   = desc.constantAsDouble("DCH_inner_cyl_R");
    this->dch_rout  = desc.constantAsDouble("DCH_outer_cyl_R");
    this->dch_Lhalf = desc.constantAsDouble("DCH_Lhalf");

    // guard wires position, fix position
    this->dch_rin_z0_guard  = desc.constantAsDouble("DCH_in_guard_z0");
    this->dch_rout_z0_guard = desc.constantAsDouble("DCH_out_guard_zL2");

    /// number of cells of first layer
    this->dch_ncell0 = desc.constantAsLong("DCH_ncell");
    /// increment of cell number by each superlayer
    this->dch_ncell_increment = desc.constantAsLong("DCH_ncell_increment");
    /// layer cells are grouped into sectors
    this->dch_ncell_per_sector = desc.constantAsLong("DCH_ncell_per_sector");
    // if dch_ncell_per_sector is not divisor of dch_ncell0 and dch_ncell_increment
    // trow an error
    if( 0 != (dch_ncell0 % dch_ncell_per_sector) || 0 != (dch_ncell_increment % dch_ncell_per_sector) )
        throw std::runtime_error("dch_ncell_per_sector is not divisor of dch_ncell0 or dch_ncell_increment");

    // wires/cells are grouped in layers (and superlayers)
    this->dch_nsuperlayers = desc.constantAsLong("DCH_nsuperlayers");
    this->dch_nlayersPerSuperlayer = desc.constantAsLong("DCH_nlayersPerSuperlayer");
    /// nlayers = nsuperlayers * nlayersPerSuperlayer
    /// default: 112 = 14 * 8
    this->dch_nlayers = this->dch_nsuperlayers * this->dch_nlayersPerSuperlayer;

    MyAngle_t dch_alpha = desc.constantAsDouble("DCH_alpha");
    this->dch_twist_angle = 2*dch_alpha;

    // retrieve some initial values for the first layer
    double dch_first_width = desc.constantAsDouble("DCH_first_width");
    MyLength_t dch_first_sense_r = desc.constantAsDouble("DCH_first_sense_r");

    // intialize layer 1 from input parameters
    {
        DCH_info_layer layer1_info;
        layer1_info.layer         = 1;
        layer1_info.nwires        = 2*this->dch_ncell0;
        layer1_info.height_z0     = dch_first_width;
        layer1_info.radius_sw_z0  = dch_first_sense_r;
        layer1_info.radius_fdw_z0 = dch_first_sense_r - 0.5*dch_first_width;
        layer1_info.radius_fuw_z0 = dch_first_sense_r + 0.5*dch_first_width;
        layer1_info.width_z0      = TWOPI*dch_first_sense_r/this->dch_ncell0;

        this->database.emplace(layer1_info.layer, layer1_info);
    }

    // some parameters of the following layer are calculated based on the previous ones
    // the rest are left as methods of DCH_info or DCH_info_layer class
    // loop over all layers, skipping the first one
    for(int ilayer = 2; ilayer<= this->dch_nlayers; ++ilayer)
    {
        // initialize here an object that will contain one row from the excel table
        DCH_info_layer layer_info;

        // the loop counter actually corresponds to the layer number
        layer_info.layer = ilayer;
        // calculate number of superlayer
        // WARNING: division of integers on purpose!
        int nsuperlayer_minus_1 = (ilayer-1)/dch_nlayersPerSuperlayer;
        layer_info.nwires = 2*(this->dch_ncell0 + this->dch_ncell_increment*nsuperlayer_minus_1 );

        // the previous layer info is needed to calculate parameters of current layer
        auto previousLayer = this->database.at(ilayer-1);

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



        this->database.emplace(ilayer, layer_info);
    }

    std::cout << "\t+ Total size of DCH database = " << DCH_info::database.size() << std::endl;
    return;
}

void DCH_info::Show_DCH_info_database()
{
    if( this->IsDatabaseEmpty() )
        throw std::runtime_error("Can not show empty database!");

    std::cout << "\n";
    std::cout << "Global parameters of DCH:\n";
    std::cout << "\tHalf length/mm = " << dch_Lhalf/dd4hep::mm << '\n';
    std::cout << "\tRadius in/mm  = " << dch_rin/dd4hep::mm << '\n';
    std::cout << "\tRadius out/mm = " << dch_rout/dd4hep::mm<< '\n';
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
                << "\t" << l.StereoSign()*this->stereoangle_z0(l.radius_sw_z0)/dd4hep::deg
                << "\t" << l.Pitch_z0(l.radius_sw_z0)/dd4hep::mm
                << "\t" << this->Radius_zLhalf(l.radius_sw_z0)/dd4hep::mm
                << "\t" << this->WireLength(l.layer,l.radius_sw_z0)/dd4hep::mm
                << "\n" << std::endl;
    }
    return;
}

} // close namespace kk

DECLARE_DETELEMENT(DriftChamber_o1_v02_T, DCH_v2::create_DCH_o1_v02)

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
#include "DD4hep/Detector.h"

#include "DCH_info.h"

namespace DCH_v2 {

using DCH_length_t = dd4hep::rec::DCH_info_struct::DCH_length_t;
using DCH_angle_t  = dd4hep::rec::DCH_info_struct::DCH_angle_t;
using DCH_layer    = dd4hep::rec::DCH_info_struct::DCH_layer;

/// Function to build DCH
static dd4hep::Ref_t create_DCH_o1_v02(dd4hep::Detector &desc, dd4hep::xml::Handle_t handle, dd4hep::SensitiveDetector sens)
{
    dd4hep::xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    dd4hep::DetElement det(detName, detID);
    sens.setType("tracker");

    dd4hep::rec::DCH_info * DCH_i = new dd4hep::rec::DCH_info();
    {
        // DCH outer geometry dimensions
        DCH_i->Set_rin  ( desc.constantAsDouble("DCH_inner_cyl_R") );
        DCH_i->Set_rout ( desc.constantAsDouble("DCH_outer_cyl_R") );
        DCH_i->Set_lhalf( desc.constantAsDouble("DCH_Lhalf")       );

        // guard wires position, fix position
        DCH_i->Set_guard_rin_at_z0 ( desc.constantAsDouble("DCH_guard_inner_r_at_z0" ) );
        DCH_i->Set_guard_rout_at_zL2( desc.constantAsDouble("DCH_guard_outer_r_at_zL2") );

        DCH_angle_t dch_alpha = desc.constantAsDouble("DCH_alpha");
        DCH_i->Set_twist_angle( 2*dch_alpha );

        DCH_i->Set_nsuperlayers(desc.constantAsLong("DCH_nsuperlayers"));
        DCH_i->Set_nlayersPerSuperlayer(desc.constantAsLong("DCH_nlayersPerSuperlayer"));

        DCH_i->Set_ncell0(desc.constantAsLong("DCH_ncell"));
        DCH_i->Set_ncell_increment(desc.constantAsLong("DCH_ncell_increment"));
        DCH_i->Set_ncell_per_sector(desc.constantAsLong("DCH_ncell_per_sector"));

        DCH_i->Set_first_width  ( desc.constantAsDouble("DCH_first_width")   );
        DCH_i->Set_first_sense_r( desc.constantAsDouble("DCH_first_sense_r") );

    }
    DCH_i->BuildLayerDatabase();
    if( DCH_i->IsDatabaseEmpty() )
        throw std::runtime_error("Empty database");

    bool printExcelTable = detElem.attr<bool>(_Unicode(printExcelTable));
    if(printExcelTable)
        DCH_i->Show_DCH_info_database(std::cout);

    auto gasElem    = detElem.child("gas");
    auto gasvolMat  = desc.material(gasElem.attr<std::string>(_Unicode(material)));
    auto gasvolVis  = desc.visAttributes(gasElem.attr<std::string>(_Unicode(vis)));

    auto vesselElem = detElem.child("vessel");
    auto vesselSkinMat  = desc.material(vesselElem.attr<std::string>(_Unicode(material)));
    auto vesselSkinVis  = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(vis)));


    auto wiresElem = detElem.child("wires");
    auto wiresVis  = desc.visAttributes(wiresElem.attr<std::string>(_Unicode(vis)));
    bool buildSenseWires  = wiresElem.attr<bool>(_Unicode(buildSenseWires));
    bool buildFieldWires  = wiresElem.attr<bool>(_Unicode(buildFieldWires));

    DCH_length_t dch_SWire_thickness        = wiresElem.attr<double>(_Unicode(SWire_thickness)) ;
    DCH_length_t dch_FSideWire_thickness    = wiresElem.attr<double>(_Unicode(FSideWire_thickness)) ;
    DCH_length_t dch_FCentralWire_thickness = wiresElem.attr<double>(_Unicode(FCentralWire_thickness)) ;

    auto dch_SWire_material        = desc.material(wiresElem.attr<std::string>(_Unicode(SWire_material)) ) ;
    auto dch_FSideWire_material    = desc.material(wiresElem.attr<std::string>(_Unicode(FSideWire_material)) ) ;
    auto dch_FCentralWire_material = desc.material(wiresElem.attr<std::string>(_Unicode(FCentralWire_material)) ) ;

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

    DCH_length_t safety_r_interspace   = 1    * dd4hep::nm;
    DCH_length_t safety_z_interspace   = 1    * dd4hep::nm;
    DCH_length_t safety_phi_interspace = 1e-6 * dd4hep::rad;


    DCH_length_t vessel_thickness = desc.constantAsDouble("DCH_vessel_thickness");
    dd4hep::Tube vessel_s(  DCH_i->rin - vessel_thickness,
                            DCH_i->rout+ vessel_thickness,
                            DCH_i->Lhalf  + vessel_thickness);
    dd4hep::Volume vessel_v (detName+"_vessel", vessel_s,  vesselSkinMat );
    vessel_v.setVisAttributes( vesselSkinVis );
    vessel_v.setRegion  ( desc, detElem.regionStr() );
    vessel_v.setLimitSet( desc, detElem.limitsStr() );


    dd4hep::Tube gas_s( DCH_i->rin,
                        DCH_i->rout,
                        DCH_i->Lhalf + 2*safety_z_interspace );
    dd4hep::Volume gas_v (detName+"_gas", gas_s, gasvolMat );
    gas_v.setVisAttributes( gasvolVis );
    vessel_v.placeVolume(gas_v);

    for(const auto& [ilayer, l]  : DCH_i->database )
    {

        // // // // // // // // // // // // // // // // // // // // /
        // // // // // INITIALIZATION OF THE LAYER // // // // // //
        // // // // // // // // // // // // // // // // // // // //
        // Hyperboloid parameters:
        /// inner radius at z=0
        DCH_length_t rin   = l.radius_fdw_z0+safety_r_interspace;
        /// inner stereoangle, calculated from rin(z=0)
        DCH_angle_t  stin  = DCH_i->stereoangle_z0(rin);
        /// outer radius at z=0
        DCH_length_t rout  = l.radius_fuw_z0-safety_r_interspace;
        /// outer stereoangle, calculated from rout(z=0)
        DCH_angle_t  stout = DCH_i->stereoangle_z0(rout);
        /// half-length
        DCH_length_t dz    = DCH_i->Lhalf + safety_z_interspace;

        dd4hep::Hyperboloid layer_s(rin, stin, rout, stout, dz);


        std::string layer_name = detName+"_layer"+std::to_string(ilayer);
        dd4hep::Volume layer_v ( layer_name , layer_s, gasvolMat );
        layer_v.setVisAttributes( desc.visAttributes( Form("dch_layer_vis%d", ilayer%22) ) );
        auto layer_pv = gas_v.placeVolume(layer_v);
        layer_pv.addPhysVolID("layer", ilayer);
        // add superlayer bitfield
        int nsuperlayer_minus_1 = DCH_i->Get_nsuperlayer_minus_1(ilayer);
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
        DCH_angle_t phi_step = (TMath::TwoPi()/ncells)*dd4hep::rad;

        // unitary cell (Twisted tube) is repeated for each layer l.nwires/2 times
        // Twisted tube parameters
        DCH_angle_t cell_twistangle    = l.StereoSign() * DCH_i->twist_angle;
        DCH_length_t cell_rin_z0       = l.radius_fdw_z0 + 2*safety_r_interspace;
        DCH_length_t cell_rout_z0      = l.radius_fuw_z0 - 2*safety_r_interspace;
        DCH_length_t cell_rin_zLhalf   = DCH_i->Radius_zLhalf(cell_rin_z0);
        DCH_length_t cell_rout_zLhalf  = DCH_i->Radius_zLhalf(cell_rout_z0);
        DCH_length_t cell_dz           = DCH_i->Lhalf;
        DCH_angle_t cell_phi_width     = phi_step - safety_phi_interspace;
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
            DCH_length_t cell_rave_z0 = 0.5*(cell_rin_z0+cell_rout_z0);
            DCH_length_t cell_swire_radius = dch_SWire_thickness/2;
            DCH_length_t swlength = 0.5*DCH_i->WireLength(ilayer,cell_rave_z0)
                                - cell_swire_radius*cos(DCH_i->stereoangle_z0(cell_rave_z0))
                                - safety_z_interspace;
            if(buildSenseWires)
            {
                dd4hep::Tube swire_s(0., dch_SWire_thickness, swlength);
                dd4hep::Volume swire_v(cell_name+"_swire", swire_s, dch_SWire_material);
                swire_v.setVisAttributes( wiresVis );
                // Change sign of stereo angle to place properly the wire inside the twisted tube
                dd4hep::RotationX stereoTr( (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(cell_rave_z0) );
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
                auto fwire_phi_offset = [&](DCH_length_t radial_distance, DCH_length_t wire_radius)->DCH_angle_t
                {
                    return atan(wire_radius/radial_distance)*dd4hep::rad + safety_phi_interspace;
                };

                // // // // // // // // // // // // // // // // // // // //
                // // // // // // POSITIONING OF F WIRE 2 // // // // // //
                // // // // // // REQUIRES OFFSET OF PHI  // // // // // //
                // // // // // // // // // // // // // // // // // // // //
                {
                    DCH_length_t fwire_radius = dch_FCentralWire_thickness/2;
                    DCH_length_t fwire_r_z0   = cell_rave_z0;
                    DCH_angle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
                    DCH_angle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                    DCH_length_t fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f2wire", fwire_s, dch_FCentralWire_material );
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
                    DCH_length_t fwire_radius = dch_FSideWire_thickness/2;
                    // decrease radial distance, move it closer to the sense wire
                    DCH_length_t fwire_r_z0   = cell_rout_z0 - fwire_radius;
                    DCH_angle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
                    DCH_angle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                    DCH_length_t fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f1wire", fwire_s, dch_FSideWire_material );
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
                    DCH_length_t fwire_radius = dch_FSideWire_thickness/2;
                    // increase radial distance, move it closer to the sense wire
                    DCH_length_t fwire_r_z0   = cell_rin_z0 + fwire_radius;
                    DCH_angle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
                    DCH_angle_t  fwire_phi    = -cell_phi_width/2 + fwire_phi_offset( fwire_r_z0, fwire_radius);
                    DCH_length_t fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f3wire", fwire_s, dch_FSideWire_material );
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
                    DCH_length_t fwire_radius = dch_FSideWire_thickness/2;
                    // increase radial distance, move it closer to the sense wire
                    DCH_length_t fwire_r_z0   = cell_rin_z0 + fwire_radius;
                    DCH_angle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
                    DCH_length_t fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f5wire", fwire_s, dch_FSideWire_material );
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
                    DCH_length_t fwire_radius = dch_FSideWire_thickness/2;
                    // increase radial distance, move it closer to the sense wire
                    DCH_length_t fwire_r_z0   = cell_rout_z0 - fwire_radius;
                    DCH_angle_t  fwire_stereo =  (-1.)*l.StereoSign()*DCH_i->stereoangle_z0(fwire_r_z0);
                    DCH_length_t fwire_length = 0.5*DCH_i->WireLength(ilayer, fwire_r_z0)
                                            - fwire_radius*cos(DCH_i->stereoangle_z0(fwire_r_z0))
                                            - safety_z_interspace;

                    dd4hep::Tube fwire_s(0., fwire_radius, fwire_length);
                    dd4hep::Volume fwire_v(cell_name+"_f4wire", fwire_s, dch_FSideWire_material );
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
            DCH_angle_t cell_phi_angle = phi_step * nphi + 0.25*cell_phi_width*(ilayer%2);
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
    det.addExtension<dd4hep::rec::DCH_info>(DCH_i);

    return det;

}

}; // end DCH_v2 namespace


DECLARE_DETELEMENT(DriftChamber_o1_v02_T, DCH_v2::create_DCH_o1_v02)

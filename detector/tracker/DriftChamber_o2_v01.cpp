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
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "TMath.h"

#include <map>

using namespace dd4hep;

namespace DCH_o2{

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
    struct DCH_info{

        // global variables of DCH
        // inline before static allow initialization at this point
        inline static MyLength_t dch_Lhalf = {0};
        inline static MyLength_t dch_rin_z0 = {0};
        inline static MyLength_t dch_rin_z0_guard = {0};
        inline static MyLength_t dch_rout_z0 = {0};
        inline static MyLength_t dch_rout_z0_guard = {0};
        inline static int dch_ncell0 = {0};
        inline static int dch_ncell_increment = {0};
        inline static int dch_nlayersPerSuperlayer = {0};
        inline static int dch_nsuperlayers = {0};
        inline static int dch_nlayers = {0};
        inline static MyAngle_t  dch_twist_angle = {0};

        // each member corresopnds to a column in the excel table
        /// layer number
        DCH_layer layer = {0};
        /// number of wires in that layer
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
        bool IsStereoPositive() const { return (1 == layer%2);}
        /// calculate sign based on IsStereoPositive
        inline int StereoSign() const {return (IsStereoPositive()*2 -1 );}

        /// separation between wires
        MyLength_t Pitch_z0(MyLength_t r_z0) const {return TWOPI*r_z0/nwires;};

        /// tan(stereoangle) = R(z=0)   / (L/2) * tan( twist_angle/2)
        inline MyAngle_t stereoangle_z0(MyLength_t r_z0) const {return atan( r_z0/dch_Lhalf*tan(dch_twist_angle/2));}

        /// tan(stereoangle) = R(z=L/2) / (L/2) * sin( twist_angle/2)
        inline MyAngle_t stereoangle_zLhalf(MyLength_t r_zLhalf) const {return atan( r_zLhalf/dch_Lhalf*sin(dch_twist_angle/2));}

        /// WireLength = 2*dch_Lhalf/cos(atan(Pitch_z0(r_z0)/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0))
        inline MyLength_t WireLength(MyLength_t r_z0) const {return  2*dch_Lhalf/cos(atan(Pitch_z0(r_z0)/(2*dch_Lhalf)))/cos(stereoangle_z0(r_z0)) ;};

        /// map to store excel table
        inline static std::map<DCH_layer, DCH_info> database;
        static bool IsDatabaseEmpty() {return (0 == database.size() );}

        static void Fill_DCH_info_database(Detector &desc);
        static void Show_DCH_info_database();

    };

    /// Function to build ARC endcaps
    static Ref_t create_DCH_o2_v01(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
    {
        xml::DetElement detElem = handle;
        std::string detName = detElem.nameStr();
        int detID = detElem.id();
        DetElement det(detName, detID);
        sens.setType("tracker");


        DCH_info::Fill_DCH_info_database(desc);

        // Create the mother Detector element to be returned at the end
        double twist_angle = 45*dd4hep::degree;
        double rmin = 1*dd4hep::cm;
        double rmax = 5*dd4hep::cm;
        double dz = 5*dd4hep::cm;
        double dphi = 90*dd4hep::deg;
        int nsegments = 1;
        TwistedTube myshape( twist_angle,  rmin,  rmax, -dz, dz, dphi);
        // Define volume (shape+material)
        Volume siVol(detName +"_sensor", myshape, desc.material("Silicon"));
        siVol.setVisAttributes(desc.visAttributes("sensor_vis"));
        siVol.setSensitiveDetector(sens);

        double wireR = 0.5*rmin + 0.5*rmax; //0.5*(rmax+rmin);
        double wirePhi = 0.3*dphi;
        double stereoangle = atan( wireR/dz*sin(twist_angle/2/rad) );

        double wireThickness = 0.1*rmin;
        Tube ancshape( 0., wireThickness, dz/cos(stereoangle/rad) - wireThickness*tan(stereoangle/rad) );
        Volume ancVol(detName +"ancshape", ancshape, desc.material("Silicon"));
        ancVol.setVisAttributes(desc.visAttributes("cooling_vis"));

        auto ancshapeTrp = RotationZ(0.3*dphi)*Translation3D( wireR,0,0) * RotationX(-stereoangle/rad);
        auto ancshapeTr0 = RotationZ(  0     )*Translation3D( wireR,0,0) * RotationX(-stereoangle/rad);
        auto ancshapeTrm = RotationZ(-0.3*dphi)*Translation3D( wireR,0,0) * RotationX(-stereoangle/rad);

        // Place our mother volume in the world
        Volume wVol = desc.pickMotherVolume(det);

        PlacedVolume siPV = wVol.placeVolume(siVol);
        wVol.placeVolume(ancVol, ancshapeTrp);
        wVol.placeVolume(ancVol, ancshapeTr0);
        wVol.placeVolume(ancVol, ancshapeTrm);

        // Assign the system ID to our mother volume
        siPV.addPhysVolID("system", detID);

        // Associate the silicon Placed Volume to the detector element.
        det.setPlacement(siPV);

        return det;

    }

    void DCH_info::Fill_DCH_info_database(Detector & desc)
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
        DCH_info::dch_ncell0 = desc.constantAsDouble("dch_ncell");
        /// increment of cell number by each superlayer
        DCH_info::dch_ncell_increment = desc.constantAsDouble("dch_ncell_increment");

        // wires/cells are grouped in layers (and superlayers)
        DCH_info::dch_nsuperlayers = desc.constantAsLong("dch_nsuperlayers");
        DCH_info::dch_nlayersPerSuperlayer = desc.constantAsLong("dch_nlayersPerSuperlayer");
        /// nlayers = nsuperlayers * nlayersPerSuperlayer
        /// default: 112 = 14 * 8
        DCH_info::dch_nlayers = DCH_info::dch_nsuperlayers * DCH_info::dch_nlayersPerSuperlayer;

        MyAngle_t dch_alpha = desc.constantAsDouble("dch_alpha");
        MyAngle_t dch_twist_angle = 2*dch_alpha;
        DCH_info::dch_twist_angle = dch_twist_angle;

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
        DCH_info::Show_DCH_info_database();
        return;
    }

    void DCH_info::Show_DCH_info_database()
    {
        if( IsDatabaseEmpty() )
            throw std::runtime_error("Can not show empty database!");
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
            << "\t" << "WireLength/mm"
            << "\n" << std::endl;

        for(const auto& [nlayer, l]  : DCH_info::database )
        {
            std::cout
                << "\t" << l.layer
                << "\t" << l.nwires
                << "\t" << l.height_z0/mm
                << "\t" << l.width_z0/mm
                << "\t" << l.radius_fdw_z0/mm
                << "\t" << l.radius_sw_z0/mm
                << "\t" << l.radius_fuw_z0/mm
                << "\t" << l.StereoSign()*l.stereoangle_z0(l.radius_sw_z0)/deg
                << "\t" << l.Pitch_z0(l.radius_sw_z0)/mm
                << "\t" << l.WireLength(l.radius_sw_z0)/mm
                << "\n" << std::endl;
        }
        return;
    }



}; // end DCH_o2 namespace
DECLARE_DETELEMENT(DriftChamber_o2_v01, DCH_o2::create_DCH_o2_v01)

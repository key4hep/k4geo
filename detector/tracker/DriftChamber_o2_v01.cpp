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

    // use alias of types to show more clear what the variable cotains
    // if everything is double, the code is not readable

    /// type for layer number
    using DCH_layer = int;
    /// tpye for lengths
    using MyLength_t = double;
    using MyAngle_t = double;
    /// tpye for angles

    // hardcode values in the code
    constexpr double PI = TMath::Pi();
    constexpr double TWOPI = TMath::TwoPi();

    // ancillary class to store DCH excel table
    // to be filled before building the geometry!
    struct DCH_info{
        inline static MyLength_t dch_Lhalf = {0};
        inline static MyAngle_t  dch_twist_angle = {0};

        DCH_layer layer = {0};
        int nwires = {0};
        double height_z0 = {0.};
        double width_z0 = {0.};

        MyLength_t radius_sw_z0 = {0.};
        MyLength_t radius_fdw_z0 = {0.};
        MyLength_t radius_fuw_z0 = {0.};


        double pitch_z0(MyLength_t r_z0){return TWOPI*r_z0/nwires;};

        bool IsStereoPositive = {true};
        int  StereoSign(){return (IsStereoPositive*2 -1 );}

        // tan(stereoangle) = R(z=0)   / (L/2) * tan( twist_angle/2)
        MyAngle_t stereoangle_z0(MyLength_t r_z0){return atan( r_z0/dch_Lhalf*tan(dch_twist_angle/2));}

        // tan(stereoangle) = R(z=L/2) / (L/2) * sin( twist_angle/2)
        MyAngle_t stereoangle_zLhalf(MyLength_t r_zLhalf){return atan( r_zLhalf/dch_Lhalf*sin(dch_twist_angle/2));}

        /// map to store excel table
        inline static std::map<DCH_layer, DCH_info> database;

    };

    // legacy, initialization of static before C++17
    // since C++17, initialization with inline keyword
    // MyLength_t DCH_info::dch_Lhalf = 0.0;
    // MyAngle_t DCH_info::dch_twist_angle = 0.0;
    // std::map<DCH_layer, DCH_info> DCH_info::database = std::map<DCH_layer, DCH_info>();



    void Fill_DCH_info_database(Detector &desc);


    /// Function to build ARC endcaps
    static Ref_t create_DCH_o2_v01(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
    {
        xml::DetElement detElem = handle;
        std::string detName = detElem.nameStr();
        int detID = detElem.id();
        DetElement det(detName, detID);
        sens.setType("tracker");


        Fill_DCH_info_database(desc);

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

    void Fill_DCH_info_database(Detector & desc)
    {
        // do not fill twice the database
        if( 0 != DCH_info::database.size() ) return;

        double dch_inner_cyl_R_z0 = desc.constantAsDouble("dch_inner_cyl_R_z0");
        double dch_in_guard_z0 = desc.constantAsDouble("dch_in_guard_z0");
        double dch_outer_cyl_R_z0 = desc.constantAsDouble("dch_outer_cyl_R_z0");
        double dch_out_guard_zL2 = desc.constantAsDouble("dch_out_guard_zL2");
        double dch_Lhalf = desc.constantAsDouble("dch_Lhalf");
        double dch_first_sense_r = desc.constantAsDouble("dch_first_sense_r");
        double dch_ncell = desc.constantAsDouble("dch_ncell");
        double dch_ncell_increment = desc.constantAsDouble("dch_ncell_increment");
        double dch_first_width = desc.constantAsDouble("dch_first_width");
        double dch_alpha = desc.constantAsDouble("dch_alpha");
        double dch_twist_angle = 2*dch_alpha;
        int nlayers = desc.constantAsLong("dch_nlayers");
        int dch_nsuperlayers = desc.constantAsLong("dch_nsuperlayers");
        int dch_nlayersPerSuperlayer = desc.constantAsLong("dch_nlayersPerSuperlayer");

        DCH_info::dch_Lhalf = dch_Lhalf;
        DCH_info::dch_twist_angle = dch_twist_angle;

        // intialize layer 1 from input parameters
        {
            DCH_info layer1_info;
            layer1_info.layer         = 1;
            layer1_info.nwires        = 2*dch_ncell;
            layer1_info.height_z0     = dch_first_width;
            layer1_info.radius_sw_z0  = dch_first_sense_r;
            layer1_info.radius_fdw_z0 = dch_first_sense_r - 0.5*dch_first_width;
            layer1_info.radius_fuw_z0 = dch_first_sense_r + 0.5*dch_first_width;
            layer1_info.width_z0      = TWOPI*dch_first_sense_r/dch_ncell;

            DCH_info::database.emplace(1, layer1_info);
        }

        for(int ilayer = 2; ilayer<= nlayers; ++ilayer)
        {
            DCH_info layer_info;
            layer_info.layer = ilayer;
            /// WARNING: division of integers on purpose!
            int nsuperlayer_minus_1 = (ilayer-1)/dch_nlayersPerSuperlayer;
            layer_info.nwires = 2*(dch_ncell + dch_ncell_increment*nsuperlayer_minus_1 );

            auto previousLayer = DCH_info::database.at(ilayer-1);
            //calculate height_z0
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

            layer_info.radius_fdw_z0 = previousLayer.radius_fuw_z0;
            layer_info.radius_fuw_z0 = previousLayer.radius_fuw_z0 + layer_info.height_z0;
            layer_info.width_z0 = TWOPI*layer_info.radius_sw_z0/(0.5*layer_info.nwires);

            // check if width_z0 == height_z0
            if(fabs(layer_info.width_z0 - layer_info.height_z0)>1e-4)
                throw std::runtime_error("fabs(l.width_z0 - l.height_z0)>1e-4");

            // stereo angle is positive for odd layer number
            layer_info.IsStereoPositive = (1 == ilayer%2);

            DCH_info::database.emplace(ilayer, layer_info);

            auto & l = layer_info;
            std::cout << "\t" << l.layer
            << "\t" << l.nwires
            << "\t" << l.height_z0/mm
            << "\t" << l.width_z0/mm
            << "\t" << l.radius_fdw_z0/mm
            << "\t" << l.radius_sw_z0/mm
            << "\t" << l.radius_fuw_z0/mm
            << "\t" << l.StereoSign()*l.stereoangle_z0(l.radius_sw_z0)/deg
            << std::endl<< std::endl;





        }

        std::cout << "++Total size of database = " << DCH_info::database.size() << std::endl;



    }


}; // end DCH_o2 namespace
DECLARE_DETELEMENT(DriftChamber_o2_v01, DCH_o2::create_DCH_o2_v01)
